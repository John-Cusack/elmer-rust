//! The `elmer-solver` backends of a partitioned run.
//!
//! * [`ParAssembly`]: `elmer-solver`'s default assembly, plus
//!   `EnforceDirichletConditions`' parallel branch: Dirichlet rows are
//!   communicated to every rank holding a copy (`CommunicateDirichletBCs`)
//!   and scaled by the diagonal summed over ranks.
//! * [`linear_registry`]: the HUTI methods as distributed solves
//!   (`elmer_dist::dist_solve`: split matrix, block-local preconditioner,
//!   parallel diagonal scaling).
//! * [`prepare_slot`]: the per-solver keywords that route `elmer-solver`'s
//!   serial steps to the parallel ones.

use elmer_core::{
    normalize_key, CsrMatrix, ElemId, IterativeSolver, LinalgError, LinearOperator, LocalSystem, MeshView, ParamMap,
    Params, Preconditioner, SolveStats, SolveStatus, Value,
};
use elmer_dist::{dist_solve, DistError, DistSolveOptions};
use elmer_solver::{AssemblyBackend, ElmerAssembly, LinearSolverRegistry, Perm, SolverError};

use crate::context::RankContext;

/// Where [`prepare_slot`] keeps the section's own `Linear System Scaling`
/// once it has switched `elmer-solver`'s serial scaling off.
pub const SCALING_KEY: &str = "mpi run: linear system scaling";

/// The HUTI methods `elmer_dist::dist_solve` runs (IterSolve's names).
pub const PARALLEL_METHODS: [&str; 6] = ["bicgstab", "bicgstab2", "tfqmr", "cgs", "cg", "gmres"];

fn dist_err(e: DistError) -> SolverError {
    match e {
        DistError::Linalg(l) => SolverError::Linalg(l),
        other => SolverError::fatal("ParallelIter", other.to_string()),
    }
}

fn no_context() -> SolverError {
    SolverError::NotWired("parallel backend used outside a partitioned run".into())
}

/// Switches a solver section to the parallel path, given its variable's
/// dofs and this rank's `n_ghost` (owned nodes are `n_ghost..`):
///
/// * `Linear System Scaling` is moved to [`SCALING_KEY`] and set false, so
///   `elmer-solver` does not scale with the unsummed local diagonal; the
///   distributed solve scales with the summed one, as `ScaleLinearSystem`
///   does in parallel.
/// * `Nonlinear System Consistent Norm` counts each shared dof once, on its
///   owner (`ComputeNorm`, SolverBasics.F90:3216). `Norm Permutation =
///   n_ghost + 1` restricts `elmer-solver`'s parallel norm to the owned
///   nodes, which [`crate::ParMesh`] numbers last; for one dof per node that
///   is that norm. Several dofs, or a `Norm Permutation` of the case's own,
///   are refused with the consistent norm.
pub fn prepare_slot(params: &mut ParamMap, dofs: usize, n_ghost: usize) -> Result<(), SolverError> {
    let scaling = params.logical("linear system scaling").unwrap_or(true);
    params.insert(SCALING_KEY, Value::Logical(scaling));
    params.insert("Linear System Scaling", Value::Logical(false));
    if params.logical("nonlinear system consistent norm") == Some(true) {
        if params.get("norm permutation").is_some() {
            return Err(SolverError::Unsupported("Norm Permutation with Nonlinear System Consistent Norm in parallel".into()));
        }
        if dofs != 1 {
            return Err(SolverError::Unsupported(format!(
                "Nonlinear System Consistent Norm in parallel for a {dofs}-dof variable"
            )));
        }
        if n_ghost > 0 {
            params.insert("Norm Permutation", Value::Integer(vec![n_ghost as i32 + 1]));
        }
    }
    Ok(())
}

/// The section's `Linear System Scaling` as it was before [`prepare_slot`].
fn scaling_of(params: &ParamMap) -> bool {
    params.logical(SCALING_KEY).or_else(|| params.logical("linear system scaling")).unwrap_or(true)
}

/// [`ElmerAssembly`] with the parallel Dirichlet conditions. Every matrix it
/// creates is recorded in the [`RankContext`] with its dof layout.
pub struct ParAssembly;

impl AssemblyBackend for ParAssembly {
    fn build_perm(&self, mesh: &dyn MeshView, active: &[ElemId], optimize_bandwidth: bool) -> Perm {
        ElmerAssembly.build_perm(mesh, active, optimize_bandwidth)
    }

    fn build_perm_for(&self, mesh: &dyn MeshView, active: &[ElemId], params: &ParamMap) -> Result<Perm, SolverError> {
        ElmerAssembly.build_perm_for(mesh, active, params)
    }

    fn create_matrix(&self, mesh: &dyn MeshView, active: &[ElemId], perm: &Perm, dofs: usize) -> CsrMatrix<f64> {
        let a = ElmerAssembly.create_matrix(mesh, active, perm, dofs);
        if let Some(ctx) = RankContext::current() {
            ctx.register(&perm.perm, dofs, &a);
        }
        a
    }

    fn add_local(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], slots: &[u32], local: &LocalSystem) {
        ElmerAssembly.add_local(a, rhs, slots, local);
    }

    fn apply_dirichlet(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], fixed: &[(u32, f64)], params: &ParamMap) {
        if let Err(e) = self.try_apply_dirichlet(a, rhs, fixed, params) {
            panic!("apply_dirichlet: {e}");
        }
    }

    /// `EnforceDirichletConditions` (BoundaryConditionUtils.F90:3269) in
    /// parallel. Collective.
    fn try_apply_dirichlet(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        fixed: &[(u32, f64)],
        params: &ParamMap,
    ) -> Result<(), SolverError> {
        let ctx = RankContext::current().ok_or_else(no_context)?;
        let layout = ctx.layout_of(a).map_err(dist_err)?;
        let ex = &*ctx.ex;
        let n = a.nrows;

        // CommunicateDirichletBCs: a shared dof constrained on any rank is
        // constrained on all its holders, with the value they agree on.
        let mut count = vec![0.0; n];
        let mut value = vec![0.0; n];
        for &(row, v) in fixed {
            let r = row as usize;
            if r >= n {
                return Err(SolverError::Linalg(LinalgError::Dim(format!("Dirichlet row {r} >= {n}"))));
            }
            // A repeated row keeps its last value, as in the serial path.
            count[r] = 1.0;
            value[r] = v;
        }
        layout.sum_shared(ex, &mut count);
        layout.sum_shared(ex, &mut value);
        let all: Vec<(u32, f64)> =
            (0..n).filter(|&r| count[r] > 0.0).map(|r| (r as u32, value[r] / count[r])).collect();
        let mut total = [all.iter().filter(|&&(r, _)| layout.owned_pos()[r as usize] != elmer_core::NO_DOF).count() as f64];
        ex.sum_f64(&mut total);
        if total[0] == 0.0 {
            return Ok(());
        }

        let scaled = params.logical("linear system dirichlet scaling").unwrap_or_else(|| scaling_of(params));
        let method = params.string("linear system scaling method").map(normalize_key).unwrap_or_else(|| "diagonal".into());
        let diag_scaled = scaled
            && match method.as_str() {
                "diagonal" => true,
                "none" => false,
                m => return Err(SolverError::Unsupported(format!("Linear System Scaling Method = {m} in parallel"))),
            };
        // ScaleLinearSystem(ApplyScaling=.FALSE.) before any row changes.
        let s = if diag_scaled { Some(summed_diagonal_scaling(ex, &layout, a)) } else { None };

        let mode = elmer_assembly::DirichletMode {
            symmetric: params.logical("linear system symmetric").unwrap_or(false),
            scaling: elmer_assembly::DirichletScaling::None,
        };
        elmer_assembly::apply_dirichlet(a, rhs, &all, mode)?;
        if let Some(d) = s {
            // The rows are now e_k with rhs v_k; scale them by 1/d_k^2.
            for &(row, _) in &all {
                let k = row as usize;
                let mut dk = d[k];
                if dk.abs() <= f64::MIN_POSITIVE {
                    dk = 1.0;
                }
                let s = 1.0 / (dk * dk);
                for v in &mut a.values[a.row_ptr[k]..a.row_ptr[k + 1]] {
                    *v *= s;
                }
                rhs[k] *= s;
            }
        }
        Ok(())
    }
}

/// `ScaleLinearSystem` diagonal factors `1/sqrt|a_ii|` from the diagonal
/// summed over the ranks holding the dof (`ParallelSumVector`), with the
/// absolute row sum, summed likewise, standing in for a zero diagonal
/// (MatrixScaling.F90:122-200). Collective.
fn summed_diagonal_scaling(ex: &dyn elmer_dist::Exchange, layout: &elmer_dist::DofLayout, a: &CsrMatrix<f64>) -> Vec<f64> {
    let tiny = f64::MIN_POSITIVE;
    let diag_of = |i: usize| {
        (a.row_ptr[i]..a.row_ptr[i + 1]).find(|&k| a.col_idx[k] as usize == i).map_or(0.0, |k| a.values[k])
    };
    let n = a.nrows;
    let mut d: Vec<f64> = (0..n).map(diag_of).collect();
    layout.sum_shared(ex, &mut d);
    let mut zero = [if d.iter().any(|v| v.abs() <= tiny) { 1.0 } else { 0.0 }];
    ex.max_f64(&mut zero);
    if zero[0] > 0.0 {
        for (i, di) in d.iter_mut().enumerate() {
            *di = if di.abs() <= tiny {
                a.values[a.row_ptr[i]..a.row_ptr[i + 1]].iter().fold(0.0, |s, v| s + v.abs())
            } else {
                diag_of(i)
            };
        }
        layout.sum_shared(ex, &mut d);
    }
    d.iter().map(|&di| if di.abs() > tiny { 1.0 / di.abs().sqrt() } else { 1.0 }).collect()
}

/// The linear solvers of a partitioned run: [`PARALLEL_METHODS`] as
/// distributed solves. Direct methods and the IterativeMethods.F90 family
/// are not registered, so selecting them is an error (the run refuses them
/// up front).
pub fn linear_registry() -> LinearSolverRegistry {
    let mut r = LinearSolverRegistry::new();
    for m in PARALLEL_METHODS {
        r.register_iterative_with_matrix(m, dist_iterative);
        r.register_iterative(m, needs_matrix);
    }
    // The preconditioner is built per rank inside the distributed solve;
    // these entries only keep the registry from warning about the name.
    r.register_precond("diagonal", unused_precond);
    r.register_precond("ilut", unused_precond);
    r.register_precond("ilu", unused_precond);
    for k in 0..=9 {
        r.register_precond(&format!("ilu{k}"), unused_precond);
    }
    r
}

fn needs_matrix(p: &ParamMap) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError> {
    Err(SolverError::NotWired(format!(
        "parallel iterative method {} (built from the assembled matrix)",
        p.string("linear system iterative method").unwrap_or("")
    )))
}

fn unused_precond(_: &ParamMap, _: &CsrMatrix<f64>) -> Result<Box<dyn Preconditioner<f64> + Send>, SolverError> {
    Ok(Box::new(elmer_core::IdentityPrecond))
}

fn dist_iterative(p: &ParamMap, a: &CsrMatrix<f64>) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError> {
    let mut opts = DistSolveOptions::from_params(p)?;
    opts.scaling = scaling_of(p)
        && p.string("linear system scaling method").is_none_or(|m| normalize_key(m) == "diagonal");
    Ok(Box::new(DistIterative { opts, a: a.clone() }))
}

/// One distributed solve of the matrix it was built from.
struct DistIterative {
    opts: DistSolveOptions,
    a: CsrMatrix<f64>,
}

impl IterativeSolver<f64> for DistIterative {
    fn solve(
        &mut self,
        _a: &dyn LinearOperator<f64>,
        _m: &dyn Preconditioner<f64>,
        b: &[f64],
        x: &mut [f64],
        _comm: &dyn elmer_core::Communicator,
    ) -> Result<SolveStats, LinalgError> {
        let to_linalg = |e: DistError| match e {
            DistError::Linalg(l) => l,
            other => LinalgError::Unsupported(other.to_string()),
        };
        let ctx = RankContext::current()
            .ok_or_else(|| LinalgError::Unsupported("parallel solve outside a partitioned run".into()))?;
        let layout = ctx.layout_of(&self.a).map_err(to_linalg)?;
        let stats = dist_solve(&*ctx.ex, &layout, &self.a, b, x, &self.opts, None).map_err(to_linalg)?;
        // IterSolve.F90:1652-1659.
        if self.opts.huti.abort_not_converged && matches!(stats.status, SolveStatus::MaxIterations | SolveStatus::Diverged) {
            return Err(LinalgError::Breakdown {
                method: "IterSolve",
                detail: format!("failed convergence tolerances after {} iterations", stats.iterations),
            });
        }
        Ok(stats)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prepare_slot_moves_scaling_and_restricts_the_consistent_norm() {
        let mut p = ParamMap::new();
        p.insert("Nonlinear System Consistent Norm", Value::Logical(true));
        prepare_slot(&mut p, 1, 7).unwrap();
        assert_eq!(p.logical("linear system scaling"), Some(false));
        assert!(scaling_of(&p));
        assert_eq!(p.integer("norm permutation"), Some(8));

        let mut q = ParamMap::new();
        q.insert("Linear System Scaling", Value::Logical(false));
        prepare_slot(&mut q, 3, 7).unwrap();
        assert!(!scaling_of(&q));
        assert_eq!(q.get("norm permutation"), None);
    }

    #[test]
    fn consistent_norm_of_a_vector_is_refused() {
        let mut p = ParamMap::new();
        p.insert("Nonlinear System Consistent Norm", Value::Logical(true));
        assert!(matches!(prepare_slot(&mut p, 2, 7), Err(SolverError::Unsupported(_))));
    }
}
