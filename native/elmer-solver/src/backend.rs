//! Pluggable backends the framework dispatches to: matrix structure and
//! assembly (implemented by `elmer-assembly`), time discretisation
//! (`elmer-timeint`) and linear solvers (`elmer-krylov`,
//! `elmer-itermethods`, `elmer-direct`). They are trait objects and
//! registries so this crate builds before those crates merge.

use std::collections::BTreeMap;
use std::sync::Arc;

use elmer_core::{
    normalize_key, CsrMatrix, DirectSolver, ElemId, ElementBasis, IterativeSolver, LocalSystem, MeshView, ParamMap,
    Preconditioner, QuadPoint, Quadrature,
};

use elmer_expr::UdfRegistry;

use crate::error::SolverError;
use crate::model::Model;

/// A dof permutation for one solver: `perm[node]` is the slot or
/// [`elmer_core::NO_DOF`]; `n_slots` slots are used.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Perm {
    pub perm: Vec<u32>,
    pub n_slots: usize,
}

/// Non-nodal dofs of one solver (Elmer `Element = "p:..."`): edge, face
/// and bubble dofs next to the nodal ones. `elmer-pdofs` provides it
/// through a [`DofsBackend`]; a solver without one is nodal.
///
/// Slots below the node count are the nodes (node `n` is slot `n` of
/// [`perm`](Self::perm)), so nodal lookups through the variable's Perm
/// keep working.
pub trait SolverDofs: Send + Sync {
    /// Perm over every dof slot (Elmer's `Ndeg` entries, [`elmer_core::NO_DOF`]
    /// where unused) and the number of slots used.
    fn perm(&self) -> Perm;
    /// Zero-valued pattern for `dofs` components per slot (`CreateMatrix`).
    fn create_matrix(&self, dofs: usize) -> CsrMatrix<f64>;
    /// Perm slots of element `e`'s dofs, `Perm(Indexes(1:nd))` of
    /// `GetElementDOFs`; [`elmer_core::NO_DOF`] where absent.
    fn element_slots(&self, e: ElemId, out: &mut Vec<u32>);
    /// Slots beyond the nodes of boundary element `e` that
    /// `DefaultDirichletBCs` fixes to zero when a boundary condition gives
    /// the field a value that is not a procedure (the edge dofs of the
    /// parent's edge); empty when there are none.
    fn boundary_extra_slots(&self, e: ElemId, out: &mut Vec<u32>);
    /// `ElementInfo` of element `e`.
    fn element_basis(&self, e: ElemId) -> Result<Box<dyn ElementBasis + Send + Sync>, String>;
    /// `GaussPoints(Element)` of element `e`.
    fn element_rule(&self, e: ElemId) -> Result<Vec<QuadPoint>, String>;
    /// `GaussPoints(Element, np)`: the `np`-point rule of element `e` on its
    /// p reference element (`Number of Integration Points` through
    /// `GaussPointsAdapt`). The default refuses.
    fn element_rule_np(&self, e: ElemId, np: usize) -> Result<Vec<QuadPoint>, String> {
        let _ = e;
        Err(format!("GaussPoints with {np} points on non-nodal elements is not provided"))
    }
    /// The node order the basis of `e` expects (Elmer reorders p boundary
    /// elements to their parent edge, `AssignLocalNumber`).
    fn element_nodes(&self, e: ElemId) -> Vec<u32>;
    /// `GetElementNOFBDOFs`: bubbles a kernel condenses locally (0 when the
    /// bubbles are in the global system).
    fn local_bubbles(&self, e: ElemId) -> usize;
}

/// Builds a solver's [`SolverDofs`].
pub trait DofsBackend: Send + Sync {
    /// `None` when the solver section is nodal; the framework then numbers
    /// nodes with the [`AssemblyBackend`] as before. `active` are the
    /// solver's bulk elements, `dofs` its `Variable DOFs`.
    fn for_solver(
        &self,
        mesh: &dyn MeshView,
        model: &Model,
        solver: &ParamMap,
        active: &[ElemId],
        dofs: usize,
    ) -> Option<Result<Box<dyn SolverDofs>, SolverError>>;

    /// [`for_solver`](Self::for_solver) for the solver with 1-based id
    /// `solver_id` (Elmer's `Solver % SolverId`), which element definitions
    /// outside the solver section name: `Element{id}` in an Equation,
    /// `Solver id: Element` in a Body. The framework calls this one; the
    /// default ignores the id.
    fn for_solver_id(
        &self,
        mesh: &dyn MeshView,
        model: &Model,
        solver_id: usize,
        solver: &ParamMap,
        active: &[ElemId],
        dofs: usize,
    ) -> Option<Result<Box<dyn SolverDofs>, SolverError>> {
        let _ = solver_id;
        self.for_solver(mesh, model, solver, active, dofs)
    }
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

    /// [`build_perm_for`](Self::build_perm_for) with extra list-matrix rows
    /// `(node, coupled nodes)`, as `MakeListMatrix` adds radiation
    /// couplings before `OptimizeBandwidth`. `None`: the backend keeps the
    /// permutation it has.
    fn build_perm_coupled(
        &self,
        _mesh: &dyn MeshView,
        _active: &[ElemId],
        _params: &ParamMap,
        _extra: &[(u32, Vec<u32>)],
    ) -> Result<Option<Perm>, SolverError> {
        Ok(None)
    }

    /// The CSR pattern (zero values, diagonal index set) for `dofs` dofs per slot.
    fn create_matrix(&self, mesh: &dyn MeshView, active: &[ElemId], perm: &Perm, dofs: usize) -> CsrMatrix<f64>;

    /// Adds `local.stiff` and `local.force`. `slots` holds one entry per
    /// element node: its perm slot, or [`elmer_core::NO_DOF`] to skip it
    /// (Elmer's `Perm(Indexes)`); local row `dofs * i + c` is component `c`
    /// of node `i`, with `dofs = local.n / slots.len()`.
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

/// Local time discretisation (`Default1stOrderTime`, `Add1stOrderTime`;
/// `Default2ndOrderTime`, `Add2ndOrderTime`). `elmer-timeint` provides the
/// implementation.
pub trait TimeBackend: Send + Sync {
    /// Folds `local.mass` into `local.stiff`/`local.force` for one step of a
    /// first-order equation.
    fn add_time_terms(&self, local: &mut LocalSystem, terms: &TimeTerms<'_>) -> Result<(), SolverError>;

    /// Folds `local.mass` and `local.damp` into `local.stiff`/`local.force`
    /// for one step of a second-order equation. `terms.prev` holds the seven
    /// `PrevValues` columns of a `Time Derivative Order = 2` variable.
    fn add_2nd_order_time_terms(&self, _local: &mut LocalSystem, _terms: &TimeTerms<'_>) -> Result<(), SolverError> {
        Err(SolverError::NotWired("second-order time discretisation".into()))
    }
}

pub type IterativeFactory = fn(&ParamMap) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError>;
/// An iterative method built from the assembled matrix (Jacobi, Richardson,
/// SGS read its entries).
pub type IterativeMatrixFactory =
    fn(&ParamMap, &CsrMatrix<f64>) -> Result<Box<dyn IterativeSolver<f64> + Send>, SolverError>;
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
    iterative_matrix: BTreeMap<String, IterativeMatrixFactory>,
    /// Serves iterative method names not registered explicitly.
    iterative_default: Option<IterativeFactory>,
    direct: BTreeMap<String, DirectFactory>,
    precond: BTreeMap<String, PrecondFactory>,
}

/// Preconditioners IterSolve knows (IterSolve.F90:1074-1093) that are not
/// ported: refused instead of IterSolve's unknown-name fallback to none.
const KNOWN_UNPORTED_PRECONDITIONERS: &[&str] = &["bilu", "multigrid", "vanka", "auxiliary space solver", "slave", "circuit"];

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

    /// Registers a method built from the assembled matrix. It is used by
    /// [`select_for`](Self::select_for), which has the matrix, and takes
    /// precedence over a [`register_iterative`](Self::register_iterative)
    /// entry of the same name there.
    pub fn register_iterative_with_matrix(&mut self, method: &str, f: IterativeMatrixFactory) {
        self.iterative_matrix.insert(normalize_key(method), f);
    }

    /// The factory for iterative method names with no explicit entry (it
    /// reads the name itself, as IterSolve's `CASE DEFAULT` does).
    pub fn set_iterative_default(&mut self, f: IterativeFactory) {
        self.iterative_default = Some(f);
    }

    pub fn register_direct(&mut self, method: &str, f: DirectFactory) {
        self.direct.insert(normalize_key(method), f);
    }

    pub fn register_precond(&mut self, name: &str, f: PrecondFactory) {
        self.precond.insert(normalize_key(name), f);
    }

    /// Picks the solver the section asks for. `Linear System Solver` is
    /// required (SolveCore.F90:1575); the direct method defaults to `banded`
    /// (DirectSolve.F90) and the iterative one to `bicgstab`
    /// (IterSolve.F90:724). `Linear System Preconditioning` defaults to
    /// `none` (IterSolve.F90:1053); an unknown name is warned about and
    /// treated as `none`, as IterSolve.F90:1095-1097 does.
    pub fn select(&self, params: &ParamMap, log: &dyn elmer_core::Logger) -> Result<LinearSolverChoice, SolverError> {
        self.select_impl(params, None, log)
    }

    /// [`select`](Self::select) once the system is assembled, so methods
    /// registered with [`register_iterative_with_matrix`](Self::register_iterative_with_matrix)
    /// can be built.
    pub fn select_for(
        &self,
        params: &ParamMap,
        matrix: &CsrMatrix<f64>,
        log: &dyn elmer_core::Logger,
    ) -> Result<LinearSolverChoice, SolverError> {
        self.select_impl(params, Some(matrix), log)
    }

    fn select_impl(
        &self,
        params: &ParamMap,
        matrix: Option<&CsrMatrix<f64>>,
        log: &dyn elmer_core::Logger,
    ) -> Result<LinearSolverChoice, SolverError> {
        use elmer_core::Params;
        let kind = params.string("linear system solver").map(normalize_key).ok_or_else(|| {
            SolverError::fatal("SolveLinearSystem", "Give \"Linear System Solver\", e.g. \"iterative\" or \"direct\"")
        })?;
        match kind.as_str() {
            "direct" => {
                let m = params
                    .string("linear system direct method")
                    .map(normalize_key)
                    .unwrap_or_else(|| "banded".into());
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
                let with_matrix = match (matrix, self.iterative_matrix.get(&m)) {
                    (Some(a), Some(f)) => Some(f(params, a)?),
                    _ => None,
                };
                let f = self
                    .iterative
                    .get(&m)
                    .or(self.iterative_default.as_ref())
                    .ok_or_else(|| SolverError::NotWired(format!("iterative method {m}")));
                let p = params
                    .string("linear system preconditioning")
                    .map(normalize_key)
                    .unwrap_or_else(|| "none".into());
                let precond = match (p.as_str(), self.precond.get(&p)) {
                    ("none", _) => None,
                    (_, Some(f)) => Some(*f),
                    (_, None) if KNOWN_UNPORTED_PRECONDITIONERS.iter().any(|k| p.starts_with(k)) => {
                        return Err(SolverError::NotWired(format!("preconditioner {p}")))
                    }
                    (_, None) => {
                        log.warn("IterSolve", &format!("Unknown preconditioner type: {p}, feature disabled."));
                        None
                    }
                };
                let solver = match with_matrix {
                    Some(s) => s,
                    None => f?(params)?,
                };
                Ok(LinearSolverChoice::Iterative { solver, precond })
            }
            other => Err(SolverError::Unsupported(format!("Linear System Solver = {other}"))),
        }
    }
}

/// The installed backends of a [`crate::Simulation`]. The default has
/// [`crate::assembly::ElmerAssembly`], [`crate::time::ElmerTime`],
/// [`crate::linsys::elmer_registry`], no MATC engine and an empty
/// function registry.
pub struct Backends {
    pub assembly: Option<Box<dyn AssemblyBackend>>,
    pub time: Option<Box<dyn TimeBackend>>,
    pub linear: LinearSolverRegistry,
    /// MATC for `MATC` keyword values (the MATC binding, S1). Share the
    /// front end's evaluator with [`crate::Simulation::set_matc`].
    pub matc: Option<Box<dyn elmer_expr::MatcEngine>>,
    /// `Procedure` functions of keyword values (`Real Procedure "lib"
    /// "entry"`), looked up before the process-wide registry. Fixtures that
    /// read the model take [`crate::Simulation::live_model`].
    pub udfs: Arc<UdfRegistry>,
    /// Basis functions and Gauss rules for the integrals the framework
    /// itself computes (`Coupled Solver` integral constraints). `elmer-fe`
    /// provides them.
    pub fe: Option<FeBackend>,
    /// Non-nodal (p-element) dof numbering, consulted for every solver;
    /// `None` numbers every solver by nodes.
    pub dofs: Option<Box<dyn DofsBackend>>,
}

/// `ElementInfo` and `GaussPoints` for the framework's own integrals.
pub struct FeBackend {
    pub basis: Box<dyn ElementBasis + Send + Sync>,
    pub quadrature: Box<dyn Quadrature + Send + Sync>,
}

impl Default for Backends {
    fn default() -> Self {
        Self {
            assembly: Some(Box::new(crate::assembly::ElmerAssembly)),
            time: Some(Box::new(crate::time::ElmerTime)),
            linear: crate::linsys::elmer_registry(),
            matc: None,
            udfs: Arc::new(UdfRegistry::new()),
            fe: None,
            dofs: None,
        }
    }
}
