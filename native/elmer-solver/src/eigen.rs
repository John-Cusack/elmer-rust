//! `Eigen Analysis` of one solver (SolveCore.F90:1355-1495, `SolveEigenSystem`
//! and EigenSolve.F90:140 `ArpackEigenSolve`) with `elmer-eigen` as the
//! eigensolver.
//!
//! The solver assembles its stiffness into the system matrix and its mass
//! with [`crate::SolverCtx::default_update_mass`] (`Matrix % MassValues`, on
//! the same pattern). `DefaultDirichletBCs` zeroes the Dirichlet rows of the
//! mass, and with `Linear System Symmetric` the columns too
//! ([`dirichlet_mass`]: `CRS_ZeroRow`, `CRS_SetSymmDirichlet`). The solve
//! then:
//!
//! 1. scales `K` and `M` with the same diagonal (`ScaleLinearSystem(Solver,
//!    A)`, `Linear System Scaling` default on);
//! 2. solves `K x = lambda M x` in shift-invert mode, the inner solves with
//!    the section's linear solver (a direct one is factored once, as
//!    `Linear System Refactorize = False` after the first solve in
//!    `ArpackEigenSolve`);
//! 3. scales the eigenvectors back (`BackScaleLinearSystem(EigenScaling)`),
//!    sets the variable's norm from its (unchanged) values (`ComputeNorm`)
//!    and normalises the modes when asked (`ScaleEigenVectors`, after the
//!    back scaling);
//! 4. returns without `ComputeChange`, as `SolveSystem` does for an eigen
//!    system.
//!
//! The modes are kept in [`crate::Variable::eigen_values`] and
//! [`crate::Variable::eigen_vectors`].
//!
//! Not reproduced: ARPACK's own iteration. `elmer-eigen` converges to the
//! same eigenpairs within `Eigen System Convergence Tolerance`, but its
//! vectors are not ARPACK's to round-off and their signs are arbitrary.
//! Refused: complex, damped, lumped-mass, stability and parallel eigen
//! systems, block and multigrid inner solvers, and scaling methods other
//! than the diagonal one.

use std::cell::RefCell;

use elmer_core::{CsrMatrix, DirectSolver, LinalgError, ParamMap, Params};
use elmer_eigen::{EigenConfig, Mass, Normalize};

use crate::backend::LinearSolverChoice;
use crate::error::SolverError;
use crate::solver::{solve_linear, BlockEnv, SolverCtx};

const CALLER: &str = "EigenSolve";

/// `CRS_ZeroRow` / `CRS_SetSymmDirichlet` on `Matrix % MassValues`: the
/// mass rows of the Dirichlet dofs are zero, and with `Linear System
/// Symmetric` their columns as well.
pub(crate) fn dirichlet_mass(m: &mut CsrMatrix<f64>, fixed: &[(u32, f64)], params: &ParamMap) -> Result<(), SolverError> {
    if fixed.is_empty() {
        return Ok(());
    }
    let symmetric = crate::assembly::dirichlet_mode(params)?.symmetric;
    let mut on = vec![false; m.nrows];
    for &(r, _) in fixed {
        if let Some(o) = on.get_mut(r as usize) {
            *o = true;
        }
    }
    for r in 0..m.nrows {
        let row = m.row_ptr[r]..m.row_ptr[r + 1];
        if on[r] {
            m.values[row].fill(0.0);
        } else if symmetric {
            for k in row {
                if on[m.col_idx[k] as usize] {
                    m.values[k] = 0.0;
                }
            }
        }
    }
    Ok(())
}

/// Whether this visit's solve is an eigen solve (`EigenAnalysis` of
/// `SolveSystem`: `NOFEigenValues > 0` and `Eigen Analysis` on now).
pub(crate) fn active(ctx: &SolverCtx<'_>) -> bool {
    let p = ctx.params();
    ctx.nof_eigen() > 0 && p.logical("eigen analysis") == Some(true) && p.logical("skip eigen analysis") != Some(true)
}

fn unsupported(what: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("Eigen Analysis: {}", what.into()))
}

/// The inner solves of an iterative `Linear System Solver`: each `OP` call
/// solves the shifted system with the section's Krylov method, without the
/// outer scaling (`IterSolver` inside `ArpackEigenSolve` does not scale).
struct IterativeInner<'a> {
    linear: &'a crate::LinearSolverRegistry,
    params: ParamMap,
    dofs: usize,
    simulation: &'a ParamMap,
    a: RefCell<Option<CsrMatrix<f64>>>,
    log: &'a dyn elmer_core::Logger,
    comm: &'a dyn elmer_core::Communicator,
}

impl DirectSolver<f64> for IterativeInner<'_> {
    fn factor(&mut self, a: &CsrMatrix<f64>) -> Result<(), LinalgError> {
        *self.a.borrow_mut() = Some(a.clone());
        Ok(())
    }

    fn solve_in_place(&self, b: &mut [f64]) -> Result<(), LinalgError> {
        let mut a = self.a.borrow_mut();
        let a = a.as_mut().ok_or_else(|| LinalgError::Unsupported("inner solve before factor".into()))?;
        let bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        let mut x = vec![0.0; b.len()];
        let env = BlockEnv { dofs: self.dofs, ndeg: self.dofs, simulation: self.simulation, schur: None, multigrid: None, alias: None };
        solve_linear(self.linear, &env, &self.params, a, b, &mut x, bnorm, false, self.log, self.comm)
            .map_err(|e| LinalgError::Unsupported(format!("inner linear solve: {e}")))?;
        b.copy_from_slice(&x);
        Ok(())
    }
}

/// The eigen solve of `SolveSystem`; returns the variable's norm.
pub(crate) fn solve(ctx: &mut SolverCtx<'_>) -> Result<f64, SolverError> {
    let (mesh, log, comm, ix) = (ctx.mesh, ctx.log, ctx.comm, ctx.slot);
    let p = ctx.params().clone();
    if comm.size() > 1 {
        return Err(unsupported("parallel eigen systems (ParallelArpackEigenSolve)"));
    }
    if p.logical("linear system complex") == Some(true) && p.logical("linear system skip complex") != Some(true) {
        return Err(unsupported("complex systems (ArpackEigenSolveComplex)"));
    }
    if p.logical("lumped mass matrix") == Some(true) {
        return Err(unsupported("a lumped mass matrix"));
    }
    if p.logical("linear system row equilibration") == Some(true)
        || p.string("linear system scaling method").is_some_and(|m| elmer_core::normalize_key(m) != "diagonal")
    {
        return Err(unsupported("scaling other than the diagonal one"));
    }
    let kind = p.string("linear system solver").map(elmer_core::normalize_key).unwrap_or_default();
    if kind != "direct" && kind != "iterative" {
        return Err(unsupported(format!("Linear System Solver = {kind}")));
    }
    let symmetric = p.logical("linear system symmetric").unwrap_or(false);
    let mut cfg = EigenConfig::from_params(&p, symmetric).map_err(|e| match e {
        LinalgError::Unsupported(m) => unsupported(m),
        e => SolverError::from(e),
    })?;
    // ScaleEigenVectors runs after the back scaling, on the original mass.
    let normalize = std::mem::replace(&mut cfg.normalize, Normalize::None);

    let sim = &mut *ctx.sim;
    let slot = &mut sim.slots[ix];
    let id = slot.variable.ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    let sys = slot.system.as_mut().ok_or_else(|| SolverError::fatal(CALLER, "no assembled system"))?;
    let mut mass = sys.mass.take().ok_or_else(|| SolverError::fatal(CALLER, "no mass matrix (DefaultUpdateMass)"))?;
    let k = &mut sys.matrix;
    let n = k.nrows;

    // ScaleLinearSystem(Solver, A): K and MassValues with one diagonal.
    let scale = p.logical("linear system skip scaling") != Some(true) && p.logical("linear system scaling").unwrap_or(true);
    let scaling = scale.then(|| {
        let opts = elmer_sparse::ScaleOptions { method: elmer_sparse::ScalingMethod::Diagonal, apply: true, rhs_scaling: false };
        elmer_sparse::scale_system(k, None, None, &opts)
    });
    if let Some(rs) = &scaling {
        let mut v = std::mem::take(&mut mass.values);
        rs.scale_values(&mass, &mut v);
        mass.values = v;
    }

    let result = {
        let dofs = sim.vars.by_id(id).dofs;
        let direct = match kind.as_str() {
            "direct" => match sim.backends.linear.select_for(&p, k, log)? {
                LinearSolverChoice::Direct(d) => Some(d),
                LinearSolverChoice::Iterative { .. } => None,
            },
            _ => None,
        };
        let r = match direct {
            Some(mut d) => elmer_eigen::eigs(k, Mass::Consistent(&mass), &mut *d, &cfg),
            None => {
                let mut q = p.clone();
                q.insert("Linear System Scaling", elmer_core::Value::Logical(false));
                let mut inner = IterativeInner {
                    linear: &sim.backends.linear,
                    params: q,
                    dofs,
                    simulation: &sim.model.simulation,
                    a: RefCell::new(None),
                    log,
                    comm,
                };
                elmer_eigen::eigs(k, Mass::Consistent(&mass), &mut inner, &cfg)
            }
        };
        r
    };

    // BackScaleLinearSystem(EigenScaling = .TRUE.), then K and M unscaled.
    let mut result = match result {
        Ok(r) => r,
        Err(e) => {
            if let Some(rs) = &scaling {
                rs.unscale(k, None, None);
                let mut v = std::mem::take(&mut mass.values);
                rs.unscale_values(&mass, &mut v);
                mass.values = v;
            }
            sys.mass = Some(mass);
            return Err(SolverError::fatal(CALLER, e.to_string()));
        }
    };
    if let Some(rs) = &scaling {
        let d = rs.diag();
        for v in &mut result.vectors {
            for (x, s) in v.iter_mut().zip(d) {
                *x *= *s;
            }
        }
        rs.unscale(k, None, None);
        let mut v = std::mem::take(&mut mass.values);
        rs.unscale_values(&mass, &mut v);
        mass.values = v;
    }
    elmer_eigen::log_summary(&result, cfg.tol, log);
    if normalize != Normalize::None {
        elmer_eigen::scale_eigen_vectors(Mass::Consistent(&mass), &mut result.vectors, normalize == Normalize::Unity);
    }
    sys.mass = Some(mass);

    // Norm = ComputeNorm(Solver, n, x) of the unchanged solution vector.
    let var = sim.vars.by_id_mut(id);
    let spec = crate::norm::NormSpec::from_params(&p, var.dofs);
    let (vals, m) = crate::norm::select_norm_entries(&spec, &var.values, n.min(var.values.len()), var.perm.as_deref(), mesh.n_nodes(), log);
    let norm = crate::norm::compute_norm(&spec, &vals, m, comm);
    var.conv.norm = norm;
    var.eigen_values = result.values;
    var.eigen_vectors = result.vectors;
    log.info(
        "SolveLinearSystem",
        &format!("Repointing {} eigenvalue components for: {}", var.dofs, var.name),
        3,
    );
    Ok(norm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::Value;

    fn full3() -> CsrMatrix<f64> {
        CsrMatrix {
            nrows: 3,
            ncols: 3,
            row_ptr: vec![0, 3, 6, 9],
            col_idx: vec![0, 1, 2, 0, 1, 2, 0, 1, 2],
            values: vec![1.0; 9],
            diag: None,
        }
    }

    /// Dirichlet dof 1: its mass row is zero; with `Linear System
    /// Symmetric` its column too (`CRS_SetSymmDirichlet`).
    #[test]
    fn dirichlet_rows_and_symmetric_columns_of_the_mass_are_zero() {
        let mut m = full3();
        dirichlet_mass(&mut m, &[(1, 0.5)], &ParamMap::new()).unwrap();
        assert_eq!(m.values, [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let mut m = full3();
        let mut p = ParamMap::new();
        p.insert("Linear System Symmetric", Value::Logical(true));
        dirichlet_mass(&mut m, &[(1, 0.5)], &p).unwrap();
        assert_eq!(m.values, [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]);
    }
}
