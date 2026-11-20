//! Pluggable backends the framework dispatches to: matrix structure and
//! assembly (implemented by `elmer-assembly`), time discretisation
//! (`elmer-timeint`) and linear solvers (`elmer-krylov`,
//! `elmer-itermethods`, `elmer-direct`). They are trait objects and
//! registries so this crate builds before those crates merge.

use std::collections::BTreeMap;
use std::sync::Arc;

use elmer_core::{
    normalize_key, CsrMatrix, DirectSolver, ElemId, ElementBasis, IterativeSolver, LocalSystem, Logger, MeshView, ParamMap,
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
    /// [`Self::element_basis`] evaluated as the scalar `ElementInfo` does,
    /// for a solver whose Fortran calls `ElementInfo` (not
    /// `ElementInfoVec`) and needs its rounding (elmer-fe-hp `scalar.rs`).
    /// The default is [`Self::element_basis`].
    fn element_basis_scalar(&self, e: ElemId) -> Result<Box<dyn ElementBasis + Send + Sync>, String> {
        self.element_basis(e)
    }
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
    /// A discontinuous Galerkin solver's `Element % DGIndexes` (the slots
    /// are `perm[dg index]`); `None` for other layouts.
    fn dg_layout(&self) -> Option<&crate::DgIndexes> {
        None
    }
    /// `GetElementNOFBDOFs`: bubbles a kernel condenses locally (0 when the
    /// bubbles are in the global system).
    fn local_bubbles(&self, e: ElemId) -> usize;
    /// `Element % BDOFs`: the element's bubble dofs, in the global system
    /// or not.
    fn element_bdofs(&self, _e: ElemId) -> usize {
        0
    }
    /// `Element % EdgeIndexes` of element `e` (0-based global edges) with
    /// each edge's `NodeIndexes`; `None` when the mesh has no edges.
    fn element_edges(&self, _e: ElemId) -> Option<Vec<(u32, [u32; 2])>> {
        None
    }
    /// `Mesh % NumberOfEdges` (0 without edges).
    fn n_edges(&self) -> usize {
        0
    }
    /// `Find_Face(Mesh, Parent, Element)` of boundary element `e`, parent
    /// `Left` else `Right`: the mesh face's `NodeIndexes` and its 0-based
    /// global `EdgeIndexes`. `None` when the mesh has no faces.
    fn boundary_face(&self, e: ElemId) -> Option<(Vec<u32>, Vec<u32>)> {
        let _ = e;
        None
    }
    /// The parts of Elmer's `Ndeg` beyond the nodes that this layout gives
    /// the mesh: `MaxBDOFs * NumberOfBulkElements`, and `MaxEDOFs *
    /// NumberOfEdges + MaxFDOFs * NumberOfFaces`.
    fn ndeg_parts(&self) -> (usize, usize) {
        (0, 0)
    }
    /// `GaussPoints(Face)` of the mesh edge between `nodes` when the
    /// layout gives it p definitions; `None` for the nodal rule.
    fn edge_rule(&self, nodes: [u32; 2]) -> Option<Result<Vec<QuadPoint>, String>> {
        let _ = nodes;
        None
    }
    /// `isActivePElement(Element, Solver)`. The default is `false`.
    fn is_p_element(&self, e: ElemId) -> bool {
        let _ = e;
        false
    }
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

/// A matrix row set as a whole: `ZeroRow(A, row)`, then each `(col,
/// value)`, and the right-hand side `rhs`.
#[derive(Clone, Debug, PartialEq)]
pub struct DirectRow {
    pub row: u32,
    pub entries: Vec<(u32, f64)>,
    pub rhs: f64,
}

/// One nodal periodic projector: rows on the boundary of BC `bc`
/// (1-based), `inv_perm[i]` the mesh node of row `i`, `cols` mesh nodes
/// (0-based, [`elmer_core::NO_DOF`] for none) in CSR `rows`.
#[derive(Clone, Debug, Default)]
pub struct NodalProjectorView {
    pub bc: usize,
    pub inv_perm: Vec<u32>,
    pub rows: Vec<usize>,
    pub cols: Vec<u32>,
    pub values: Vec<f64>,
}

/// Sets each [`DirectRow`] (`ZeroRow` + `SetMatrixElement`); a column off
/// the row's pattern is an error.
pub fn set_direct_rows(a: &mut CsrMatrix<f64>, rhs: &mut [f64], rows: &[DirectRow]) -> Result<(), SolverError> {
    for r in rows {
        let row = r.row as usize;
        let span = a.row_ptr[row]..a.row_ptr[row + 1];
        a.values[span.clone()].iter_mut().for_each(|v| *v = 0.0);
        for &(col, v) in &r.entries {
            let at = span
                .clone()
                .find(|&k| a.col_idx[k] == col)
                .ok_or_else(|| SolverError::fatal("SetMatrixElement", format!("no entry ({row}, {col}) in the matrix")))?;
            a.values[at] = v;
        }
        rhs[row] = 0.0;
        rhs[row] += r.rhs;
    }
    Ok(())
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
    /// After a `RUN` block reload: the boundary-condition sections as they
    /// now read, for a backend that keeps its own copy (the periodic
    /// assembly's `Periodic BC Offset`). Nothing by default.
    fn reload_boundary_conditions(&mut self, _bcs: &[ParamMap]) {}

    /// The permutation for the nodes of the `active` bulk elements, with
    /// Elmer's `Optimize Bandwidth` decision applied when `optimize_bandwidth`.
    fn build_perm(&self, mesh: &dyn MeshView, active: &[ElemId], optimize_bandwidth: bool) -> Perm;

    /// [`build_perm`](Self::build_perm) from the solver section: `Optimize
    /// Bandwidth` (default true) and, where supported, `Optimize Bandwidth
    /// Use Always`.
    /// [`Self::build_perm_for`] with Elmer's permutation size `Ndeg` when it
    /// exceeds the node count (p dofs of other solvers on the mesh). The
    /// default ignores it.
    fn build_perm_ndeg(&self, mesh: &dyn MeshView, active: &[ElemId], params: &ParamMap, _ndeg: usize) -> Result<Perm, SolverError> {
        self.build_perm_for(mesh, active, params)
    }

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

    /// A solver extended the pattern of a matrix made by
    /// [`create_matrix`](Self::create_matrix) over `perm`
    /// (`UpdateMatrixTopologyWithFactors`, which ends in
    /// `ParallelInitMatrix`). Nothing by default.
    fn pattern_changed(&self, _perm: &Perm, _dofs: usize, _a: &CsrMatrix<f64>) {}

    /// Adds `local.stiff` and `local.force`. `slots` holds one entry per
    /// element node: its perm slot, or [`elmer_core::NO_DOF`] to skip it
    /// (Elmer's `Perm(Indexes)`); local row `dofs * i + c` is component `c`
    /// of node `i`, with `dofs = local.n / slots.len()`.
    fn add_local(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], slots: &[u32], local: &LocalSystem);

    /// Fixes each `(row, value)` with Elmer's Dirichlet semantics, read from
    /// the solver section (`Linear System Symmetric`, Dirichlet scaling).
    fn apply_dirichlet(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], fixed: &[(u32, f64)], params: &ParamMap);

    /// [`try_apply_dirichlet`](Self::try_apply_dirichlet) with whole rows
    /// set first (`rows`: the normal-tangential rows `SetElementValues`
    /// writes directly, BoundaryConditionUtils.F90:1782-1802), after any
    /// pass of the backend's own (the periodic dependencies) and before the
    /// Dirichlet values.
    fn try_apply_dirichlet_rows(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        rows: &[DirectRow],
        fixed: &[(u32, f64)],
        params: &ParamMap,
    ) -> Result<(), SolverError> {
        set_direct_rows(a, rhs, rows)?;
        self.try_apply_dirichlet(a, rhs, fixed, params)
    }

    /// [`try_apply_dirichlet_rows`](Self::try_apply_dirichlet_rows) with the
    /// solver's `(A % BulkValues, A % BulkRHS)` copy (`Calculate Loads`),
    /// which a backend's own pass updates as it does `a` and `rhs`
    /// (`SetPeriodicBoundariesPass1` moves the slave rows onto the masters
    /// in both). The bulk copy is left alone by default.
    fn try_apply_dirichlet_rows_bulk(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        bulk: (Option<&mut [f64]>, Option<&mut [f64]>),
        rows: &[DirectRow],
        fixed: &[(u32, f64)],
        params: &ParamMap,
    ) -> Result<(), SolverError> {
        let _ = bulk;
        self.try_apply_dirichlet_rows(a, rhs, rows, fixed, params)
    }

    /// The nodal periodic projectors (`Model % BCs(i) % PMatrix` of
    /// `PROJECTOR_TYPE_NODAL`) the backend applies, for the normal averaging
    /// of normal-tangential systems. None by default.
    fn nodal_periodic_projectors(&self) -> Vec<NodalProjectorView> {
        Vec::new()
    }

    /// The nodal weights of `Periodic BC Coefficient <name>`
    /// (`CalculateNodalWeights`, by permutation slot) for the jump rows of
    /// the periodic passes. Ignored by default.
    fn set_periodic_jump_weights(&self, _name: &str, _weights: &[f64]) {}

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
#[derive(Default, Clone)]
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
    /// The `Conforming BC` node map (`Mesh % PeriodicPerm`), which a solver
    /// with `Apply Conforming BCs` numbers its dofs by; `elmer-mortar`
    /// builds it.
    pub conforming: Option<Arc<crate::conforming::ConformingMap>>,
    /// Mortar contact for solvers with `Apply Contact BCs`
    /// (`elmer-mortar`'s `DetermineContact`).
    pub contact: Option<Arc<dyn ContactBackend>>,
    /// Sums over the mesh nodes shared between partitions; `None` in a
    /// serial run.
    pub nodes: Option<Box<dyn NodeInterface>>,
    /// The eigen solve of a partitioned run (`ParallelArpackEigenSolve`);
    /// `None` refuses parallel eigen systems.
    pub eigen: Option<Arc<dyn ParEigenBackend>>,
    /// Point evaluation on the other partitions (the parallel part of
    /// `InterpolateMeshToMesh`); `None` in a serial run.
    pub points: Option<Box<dyn PointExchange>>,
}

/// The parallel part of `InterpolateMeshToMesh` (InterpolateMeshToMesh.F90):
/// the points a rank did not find in its own piece of the source mesh go
/// to every other rank whose source bounding box (padded by a tenth of its
/// largest extent) contains them; each of those evaluates them on its own
/// piece, and the values found come back.
pub trait PointExchange: Send + Sync {
    /// Collective. `bbox` is this rank's source mesh box (`[min x, min y,
    /// min z, max x, max y, max z]`, padded), `None` for an empty piece;
    /// `points` are this rank's unfound points; `eval` evaluates points
    /// another rank sent (`None` where not found here). Returns, for each
    /// of `points`, the value of the highest rank that found it.
    #[allow(clippy::type_complexity)]
    fn exchange(
        &self,
        bbox: Option<[f64; 6]>,
        points: &[[f64; 3]],
        eval: &mut dyn FnMut(&[[f64; 3]]) -> Vec<Option<f64>>,
    ) -> Result<Vec<Option<f64>>, SolverError>;
}

/// `ScaleLinearSystem`, `ParallelArpackEigenSolve` (ParallelEigenSolve.F90:82)
/// and `BackScaleLinearSystem(EigenScaling)` of a partitioned run, for the
/// real consistent-mass eigen system `k x = lambda mass x` whose rank-local
/// pieces are `k` and `mass` (both restored on return). The vectors come
/// back over all local dofs. Collective.
pub trait ParEigenBackend: Send + Sync {
    #[allow(clippy::too_many_arguments)]
    fn solve(
        &self,
        params: &ParamMap,
        k: &mut CsrMatrix<f64>,
        mass: &mut CsrMatrix<f64>,
        cfg: &elmer_eigen::EigenConfig,
        symmetric: bool,
        arpack: &mut elmer_eigen::arpack::saupd::Arpack,
        log: &dyn Logger,
    ) -> Result<elmer_eigen::EigenResult, SolverError>;
}

/// Nodal exchange between the partitions of a parallel run.
pub trait NodeInterface: Send + Sync {
    /// Adds the other partitions' values of each shared mesh node to this
    /// rank's: its own value first, then the other holders' in ascending
    /// rank order (`SendInterface` / `RecvInterface` of
    /// `CalculateBodyAverage`, MeshBasics.F90). `v` has one value per local
    /// mesh node. Collective.
    fn sum_shared(&self, v: &mut [f64]) -> Result<(), SolverError>;
    /// Every copy of a shared node takes the least of all holders' values
    /// (`ExchangeSourceVecInt` with `op = 1`). Collective.
    fn min_shared(&self, v: &mut [f64]) -> Result<(), SolverError>;
    /// Every copy of a shared node takes the greatest of all holders'
    /// values (`op = 2`). Collective.
    fn max_shared(&self, v: &mut [f64]) -> Result<(), SolverError>;
}

/// The three places the library calls into mortar contact during a
/// solver's nonlinear iteration (DefUtils.F90). The backend keeps its own
/// per-solver state.
pub trait ContactBackend: Send + Sync {
    /// `DefaultInitialize` (DefUtils.F90:3924-3928), before the n-t setup:
    /// `Displace Mesh At Init` moves the mesh by the solution, so that the
    /// n-t normals and the projectors are those of the displaced mesh.
    fn begin_initialize(&self, _ctx: &mut crate::SolverCtx<'_>) -> Result<(), SolverError> {
        Ok(())
    }
    /// `DefaultInitialize` (DefUtils.F90:3924-3952): the `Mortar BC
    /// Nonlinear` projectors, generated on the mesh displaced by the
    /// current solution.
    fn initialize(&self, ctx: &mut crate::SolverCtx<'_>) -> Result<(), SolverError>;
    /// `DefaultFinishBoundaryAssembly` (DefUtils.F90:7933-7937):
    /// `DetermineContact` on the assembled system, before the Dirichlet
    /// conditions.
    fn finish_boundary(&self, ctx: &mut crate::SolverCtx<'_>) -> Result<(), SolverError>;
    /// `DefaultSolve` (DefUtils.F90:4543): `GenerateConstraintMatrix` from
    /// the contact projectors, set with
    /// [`crate::SolverCtx::set_generated_constraints`].
    fn before_solve(&self, ctx: &mut crate::SolverCtx<'_>) -> Result<(), SolverError>;
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
            conforming: None,
            contact: None,
            nodes: None,
            eigen: None,
            points: None,
        }
    }
}
