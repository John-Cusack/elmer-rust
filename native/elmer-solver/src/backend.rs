//! Pluggable backends the framework dispatches to: matrix structure and
//! assembly (implemented by `elmer-assembly`), time discretisation
//! (`elmer-timeint`) and linear solvers (`elmer-krylov`,
//! `elmer-itermethods`, `elmer-direct`). They are trait objects and
//! registries so this crate builds before those crates merge.

use std::collections::BTreeMap;

use elmer_core::{
    normalize_key, CsrMatrix, DirectSolver, ElemId, IterativeSolver, LocalSystem, MeshView, ParamMap,
    Preconditioner,
};

use crate::error::SolverError;

/// A dof permutation for one solver: `perm[node]` is the slot or
/// [`elmer_core::NO_DOF`]; `n_slots` slots are used.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Perm {
    pub perm: Vec<u32>,
    pub n_slots: usize,
}

/// Matrix structure and assembly (Elmer `CreateMatrix`, `DefaultUpdate*`,
/// `DefaultDirichletBCs`). [`crate::assembly::ElmerAssembly`] implements it
/// with `elmer-assembly` and is installed by default.
///
/// The framework calls [`build_perm_for`](Self::build_perm_for) and
/// [`try_apply_dirichlet`](Self::try_apply_dirichlet). `build_perm` and
/// `apply_dirichlet` are the first API's signatures, kept for one merge
/// cycle (swarm plan 4.4) so existing implementations still compile; the
/// provided methods fall back to them.
pub trait AssemblyBackend: Send + Sync {
    /// The permutation for the nodes of the `active` bulk elements, with
    /// Elmer's `Optimize Bandwidth` decision applied when `optimize_bandwidth`.
    fn build_perm(&self, mesh: &dyn MeshView, active: &[ElemId], optimize_bandwidth: bool) -> Perm;

    /// [`build_perm`](Self::build_perm) from the solver section: `Optimize
    /// Bandwidth` (default true) and, where supported, `Optimize Bandwidth
    /// Use Always`.
    fn build_perm_for(&self, mesh: &dyn MeshView, active: &[ElemId], params: &ParamMap) -> Result<Perm, SolverError> {
        use elmer_core::Params;
        if params.logical("optimize bandwidth use always") == Some(true) {
            return Err(SolverError::Unsupported(
                "Optimize Bandwidth Use Always with this assembly backend".into(),
            ));
        }
        Ok(self.build_perm(mesh, active, params.logical("optimize bandwidth").unwrap_or(true)))
    }

    /// The CSR pattern (zero values, diagonal index set) for `dofs` dofs per slot.
    fn create_matrix(&self, mesh: &dyn MeshView, active: &[ElemId], perm: &Perm, dofs: usize) -> CsrMatrix<f64>;

    /// Adds `local.stiff` and `local.force`. `slots` holds one entry per
    /// element node: its perm slot, or [`elmer_core::NO_DOF`] to skip it
    /// (Elmer's `Perm(Indexes)`); local row `dofs * i + c` is component `c`
    /// of node `i`.
    fn add_local(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], slots: &[u32], local: &LocalSystem);

    /// Fixes each `(row, value)` with Elmer's Dirichlet semantics, read from
    /// the solver section (`Linear System Symmetric`, Dirichlet scaling).
    fn apply_dirichlet(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], fixed: &[(u32, f64)], params: &ParamMap);

    /// Fallible [`apply_dirichlet`](Self::apply_dirichlet).
    fn try_apply_dirichlet(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        fixed: &[(u32, f64)],
        params: &ParamMap,
    ) -> Result<(), SolverError> {
        self.apply_dirichlet(a, rhs, fixed, params);
        Ok(())
    }
}

/// The time-history context a local time discretisation needs.
pub struct TimeTerms<'a> {
    pub dt: f64,
    /// The solver section (`Timestepping Method`, `BDF Order`, ...).
    pub params: &'a ParamMap,
    /// The `Simulation` section (method and order fallbacks).
    pub simulation: &'a ParamMap,
    /// Committed previous values at the element's dofs, newest first
    /// (`PrevValues(L, 1..)`), all the variable keeps.
    pub prev: &'a [Vec<f64>],
    /// The current iterate at the element's dofs (`Values(L)`).
    pub cur: &'a [f64],
    /// `Solver % DoneTime`: timesteps initialised so far (start-up order).
    pub done_time: usize,
    /// `Solver % Order` (BDF order, 1 otherwise).
    pub order: usize,
    /// Previous step sizes, newest first (`Timestep size` PrevValues).
    pub dts: &'a [f64],
}

/// Local time discretisation (`Default1stOrderTime`, `Add1stOrderTime`).
/// `elmer-timeint` provides the implementation.
pub trait TimeBackend: Send + Sync {
    /// Folds `local.mass` (and `local.damp` for 2nd order) into
    /// `local.stiff`/`local.force` for one step.
    fn add_time_terms(&self, local: &mut LocalSystem, terms: &TimeTerms<'_>) -> Result<(), SolverError>;
}

pub type IterativeFactory = fn(&ParamMap) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError>;
pub type DirectFactory = fn(&ParamMap) -> Result<Box<dyn DirectSolver<f64> + Send>, SolverError>;
pub type PrecondFactory =
    fn(&ParamMap, &CsrMatrix<f64>) -> Result<Box<dyn Preconditioner<f64> + Send>, SolverError>;

/// Linear-system dispatch by keyword. Keys are normalised method names:
/// `Linear System Iterative Method` (`bicgstab`, `cg`, `gmres`, ...),
/// `Linear System Direct Method` (`umfpack`, `mumps`, `banded`, ...) and
/// `Linear System Preconditioning` (`ilu0`, `diagonal`, `none`, ...).
#[derive(Default)]
pub struct LinearSolverRegistry {
    iterative: BTreeMap<String, IterativeFactory>,
    direct: BTreeMap<String, DirectFactory>,
    precond: BTreeMap<String, PrecondFactory>,
}

/// A linear solver selected from the solver section.
pub enum LinearSolverChoice {
    Iterative {
        solver: Box<dyn IterativeSolver<f64> + Send>,
        /// Factory for the preconditioner, built after assembly.
        precond: Option<PrecondFactory>,
    },
    Direct(Box<dyn DirectSolver<f64> + Send>),
}

impl LinearSolverRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn register_iterative(&mut self, method: &str, f: IterativeFactory) {
        self.iterative.insert(normalize_key(method), f);
    }

    pub fn register_direct(&mut self, method: &str, f: DirectFactory) {
        self.direct.insert(normalize_key(method), f);
    }

    pub fn register_precond(&mut self, name: &str, f: PrecondFactory) {
        self.precond.insert(normalize_key(name), f);
    }

    /// Picks the solver the section asks for. `Linear System Solver`
    /// defaults to `iterative` and the method to `bicgstab` (IterSolve.F90:724);
    /// the preconditioner defaults to `none` (IterSolve.F90:1053), which maps to `None`.
    pub fn select(&self, params: &ParamMap) -> Result<LinearSolverChoice, SolverError> {
        use elmer_core::Params;
        let kind = params
            .string("linear system solver")
            .map(normalize_key)
            .unwrap_or_else(|| "iterative".into());
        match kind.as_str() {
            "direct" => {
                let m = params
                    .string("linear system direct method")
                    .map(normalize_key)
                    .unwrap_or_else(|| "umfpack".into());
                let f = self
                    .direct
                    .get(&m)
                    .ok_or_else(|| SolverError::NotWired(format!("direct method {m}")))?;
                Ok(LinearSolverChoice::Direct(f(params)?))
            }
            "iterative" => {
                let m = params
                    .string("linear system iterative method")
                    .map(normalize_key)
                    .unwrap_or_else(|| "bicgstab".into());
                let f = self
                    .iterative
                    .get(&m)
                    .ok_or_else(|| SolverError::NotWired(format!("iterative method {m}")))?;
                let p = params
                    .string("linear system preconditioning")
                    .map(normalize_key)
                    .unwrap_or_else(|| "none".into());
                let precond = if p == "none" {
                    None
                } else {
                    Some(
                        *self
                            .precond
                            .get(&p)
                            .ok_or_else(|| SolverError::NotWired(format!("preconditioner {p}")))?,
                    )
                };
                Ok(LinearSolverChoice::Iterative { solver: f(params)?, precond })
            }
            other => Err(SolverError::Unsupported(format!("linear system solver {other}"))),
        }
    }
}

/// The installed backends of a [`crate::Simulation`]. The default has
/// [`crate::assembly::ElmerAssembly`], [`crate::time::ElmerTime`], an
/// empty linear solver registry and no MATC engine.
pub struct Backends {
    pub assembly: Option<Box<dyn AssemblyBackend>>,
    pub time: Option<Box<dyn TimeBackend>>,
    pub linear: LinearSolverRegistry,
    /// MATC for `MATC` keyword values (the MATC binding, S1).
    pub matc: Option<Box<dyn elmer_expr::MatcEngine>>,
}

impl Default for Backends {
    fn default() -> Self {
        Self {
            assembly: Some(Box::new(crate::assembly::ElmerAssembly)),
            time: Some(Box::new(crate::time::ElmerTime)),
            linear: LinearSolverRegistry::default(),
            matc: None,
        }
    }
}
