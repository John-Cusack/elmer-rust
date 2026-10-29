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
/// `DefaultDirichletBCs`). `elmer-assembly` provides the implementation.
pub trait AssemblyBackend: Send + Sync {
    /// The permutation for the nodes of `active` (bulk) elements, with
    /// Elmer's `Optimize Bandwidth` decision applied when `optimize_bandwidth`.
    fn build_perm(&self, mesh: &dyn MeshView, active: &[ElemId], optimize_bandwidth: bool) -> Perm;

    /// The CSR pattern (zero values) for `dofs` dofs per slot.
    fn create_matrix(&self, mesh: &dyn MeshView, active: &[ElemId], perm: &Perm, dofs: usize) -> CsrMatrix<f64>;

    /// Adds `local.stiff` and `local.force` at the global `dof_index` rows/cols.
    fn add_local(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], dof_index: &[u32], local: &LocalSystem);

    /// Fixes each `(dof, value)` with Elmer's Dirichlet semantics (keeping
    /// symmetry or eliminating rows per `params`).
    fn apply_dirichlet(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], fixed: &[(u32, f64)], params: &ParamMap);
}

/// The time-history context a local time discretisation needs.
pub struct TimeTerms<'a> {
    pub dt: f64,
    /// The solver section (`Timestepping Method`, `BDF Order`, ...).
    pub params: &'a ParamMap,
    /// Previous values at the element's dofs, newest first.
    pub prev: &'a [Vec<f64>],
    /// 1-based time step number (BDF start-up order).
    pub timestep: usize,
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

/// The installed backends of a [`crate::Simulation`].
#[derive(Default)]
pub struct Backends {
    pub assembly: Option<Box<dyn AssemblyBackend>>,
    pub time: Option<Box<dyn TimeBackend>>,
    pub linear: LinearSolverRegistry,
}
