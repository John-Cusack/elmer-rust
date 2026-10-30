//! [`ElmerAssembly`]: the [`AssemblyBackend`] built on `elmer-assembly`
//! (Elmer `CreateMatrix`, `UpdateGlobalEquations`,
//! `EnforceDirichletConditions`).

use elmer_assembly::{CsrAssembler, DirichletMode, DirichletScaling, PermOptions};
use elmer_core::{normalize_key, CsrMatrix, ElemId, GlobalAssembler, LocalSystem, MeshView, ParamMap, Params};

use crate::backend::{AssemblyBackend, Perm};
use crate::error::SolverError;

/// The default assembly backend.
#[derive(Copy, Clone, Debug, Default)]
pub struct ElmerAssembly;

/// `active(e)` over bulk elements from the slot's element list.
fn active_mask(mesh: &dyn MeshView, active: &[ElemId]) -> Vec<bool> {
    let mut m = vec![false; mesh.n_bulk() + mesh.n_boundary()];
    for e in active {
        m[e.0 as usize] = true;
    }
    m
}

/// `Optimize Bandwidth` (default true, MainUtils.F90:1724) and `Optimize
/// Bandwidth Use Always` (default false).
/// A `No Matrix` solver keeps its initial permutation: `CreateMatrix`
/// returns before `OptimizeBandwidth` (ElementUtils.F90:2229).
pub fn perm_options(params: &ParamMap) -> PermOptions {
    PermOptions {
        optimize_bandwidth: params.logical("optimize bandwidth").unwrap_or(true)
            && params.logical("no matrix") != Some(true),
        use_optimized_always: params.logical("optimize bandwidth use always").unwrap_or(false),
    }
}

/// The Dirichlet matrix options of a solver section: `Linear System
/// Symmetric` (default false); `Linear System Dirichlet Scaling`, falling
/// back to `Linear System Scaling` (default true); `Linear System Scaling
/// Method` (default `diagonal`). Other scaling methods are refused.
pub fn dirichlet_mode(params: &ParamMap) -> Result<DirichletMode, SolverError> {
    let scaled = params
        .logical("linear system dirichlet scaling")
        .or_else(|| params.logical("linear system scaling"))
        .unwrap_or(true);
    let method = params
        .string("linear system scaling method")
        .map(normalize_key)
        .unwrap_or_else(|| "diagonal".into());
    let scaling = if !scaled {
        DirichletScaling::None
    } else {
        match method.as_str() {
            "diagonal" => DirichletScaling::Diagonal,
            "none" => DirichletScaling::None,
            m => return Err(SolverError::Unsupported(format!("Linear System Scaling Method = {m}"))),
        }
    };
    Ok(DirichletMode {
        symmetric: params.logical("linear system symmetric").unwrap_or(false),
        scaling,
    })
}

fn to_perm(p: elmer_assembly::Perm) -> Perm {
    Perm {
        perm: p.perm,
        n_slots: p.n_dofs,
    }
}

impl AssemblyBackend for ElmerAssembly {
    fn build_perm(&self, mesh: &dyn MeshView, active: &[ElemId], optimize_bandwidth: bool) -> Perm {
        let m = active_mask(mesh, active);
        to_perm(elmer_assembly::build_perm(mesh, &|e: ElemId| m[e.0 as usize], optimize_bandwidth))
    }

    fn build_perm_for(&self, mesh: &dyn MeshView, active: &[ElemId], params: &ParamMap) -> Result<Perm, SolverError> {
        let m = active_mask(mesh, active);
        let (p, _) = elmer_assembly::build_perm_with(mesh, &|e: ElemId| m[e.0 as usize], perm_options(params));
        Ok(to_perm(p))
    }

    fn build_perm_coupled(
        &self,
        mesh: &dyn MeshView,
        active: &[ElemId],
        params: &ParamMap,
        extra: &[(u32, Vec<u32>)],
    ) -> Result<Option<Perm>, SolverError> {
        let m = active_mask(mesh, active);
        let (p, _) = elmer_assembly::build_perm_coupled(mesh, &|e: ElemId| m[e.0 as usize], perm_options(params), extra);
        Ok(Some(to_perm(p)))
    }

    fn create_matrix(&self, mesh: &dyn MeshView, active: &[ElemId], perm: &Perm, dofs: usize) -> CsrMatrix<f64> {
        let m = active_mask(mesh, active);
        let p = elmer_assembly::Perm::from_perm(perm.perm.clone());
        elmer_assembly::create_matrix(mesh, &|e: ElemId| m[e.0 as usize], &p, dofs)
    }

    fn add_local(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], slots: &[u32], local: &LocalSystem) {
        CsrAssembler::new(a, rhs).add_local(slots, local);
    }

    /// Legacy entry point; panics where [`try_apply_dirichlet`](AssemblyBackend::try_apply_dirichlet)
    /// would return an error. The framework never calls it.
    fn apply_dirichlet(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], fixed: &[(u32, f64)], params: &ParamMap) {
        if let Err(e) = self.try_apply_dirichlet(a, rhs, fixed, params) {
            panic!("apply_dirichlet: {e}");
        }
    }

    fn try_apply_dirichlet(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        fixed: &[(u32, f64)],
        params: &ParamMap,
    ) -> Result<(), SolverError> {
        let mode = dirichlet_mode(params)?;
        elmer_assembly::apply_dirichlet(a, rhs, fixed, mode)?;
        Ok(())
    }
}
