//! [`PeriodicAssembly`]: `elmer-solver`'s default assembly backend with the
//! periodic projectors wired in, so a `Simulation` runs `Periodic BC` cases
//! without changes to the framework.
//!
//! - [`build_perm_for`](AssemblyBackend::build_perm_for) and
//!   [`create_matrix`](AssemblyBackend::create_matrix) add the projector
//!   couplings to the node graph before the bandwidth decision and the CSR
//!   pattern (`MakeListMatrix` with `DoProjectors`, the default for every
//!   solver).
//! - [`try_apply_dirichlet`](AssemblyBackend::try_apply_dirichlet) applies
//!   [`set_periodic_boundaries`] to every component of the solver's variable,
//!   then the Dirichlet values, as `SetDirichletBoundaries` sets the periodic
//!   rows before `EnforceDirichletConditions`.
//!
//! Bulk (`Calculate Loads`) copies are not updated, and a block system whose
//! size does not match its variable's permutation is refused when periodic.

use std::collections::HashMap;
use std::sync::Mutex;

use elmer_assembly::{
    bandwidth, bandwidth_reordered, initial_permutation, list_matrix, rcm_order, ListMatrix,
    Perm as AsmPerm,
};
use elmer_core::{CsrMatrix, ElemId, LocalSystem, MeshView, ParamMap, Params, NO_DOF};
use elmer_solver::assembly::perm_options;
use elmer_solver::sim::parse_variable_spec;
use elmer_solver::{AssemblyBackend, ElmerAssembly, Perm, SolverError};

use crate::constrain::{
    add_projector_topology, periodic_active, set_periodic_boundaries, PeriodicTarget,
};
use crate::projector::Projector;
use crate::MortarError;

impl From<MortarError> for SolverError {
    fn from(e: MortarError) -> Self {
        match e {
            MortarError::Unsupported(s) => SolverError::Unsupported(s),
            e => SolverError::fatal("SetPeriodicBoundaries", e.to_string()),
        }
    }
}

/// The default assembly plus the boundary projectors of one mesh.
pub struct PeriodicAssembly {
    inner: ElmerAssembly,
    bcs: Vec<ParamMap>,
    projectors: Vec<Option<Projector>>,
    n_nodes: usize,
    /// Solver `Variable` keyword -> its node permutation, recorded by
    /// `build_perm_for` for the Dirichlet step, which only sees the solver
    /// section.
    perms: Mutex<HashMap<String, Vec<u32>>>,
}

impl PeriodicAssembly {
    /// `bcs[i]` and `projectors[i]` belong to boundary condition `i + 1`
    /// (the output of [`crate::generate_periodic_projectors`]).
    pub fn new(bcs: Vec<ParamMap>, projectors: Vec<Option<Projector>>, n_nodes: usize) -> Self {
        PeriodicAssembly {
            inner: ElmerAssembly,
            bcs,
            projectors,
            n_nodes,
            perms: Mutex::new(HashMap::new()),
        }
    }

    fn bc_views(&self) -> Vec<&dyn Params> {
        self.bcs.iter().map(|b| b as &dyn Params).collect()
    }

    fn graph(
        &self,
        mesh: &dyn MeshView,
        active: &dyn Fn(ElemId) -> bool,
        reorder: &AsmPerm,
    ) -> ListMatrix {
        let mut list = list_matrix(mesh, active, reorder);
        add_projector_topology(&mut list, &reorder.perm, &self.bc_views(), &self.projectors);
        list
    }
}

/// `active(e)` from the slot's element list, as `ElmerAssembly` builds it.
fn active_mask(mesh: &dyn MeshView, active: &[ElemId]) -> Vec<bool> {
    let mut m = vec![false; mesh.n_bulk() + mesh.n_boundary()];
    for e in active {
        m[e.0 as usize] = true;
    }
    m
}

impl AssemblyBackend for PeriodicAssembly {
    fn build_perm(&self, mesh: &dyn MeshView, active: &[ElemId], optimize_bandwidth: bool) -> Perm {
        let mut params = ParamMap::new();
        params.insert(
            "optimize bandwidth",
            elmer_core::Value::Logical(optimize_bandwidth),
        );
        self.build_perm_for(mesh, active, &params)
            .expect("build_perm without Optimize Bandwidth Use Always cannot fail")
    }

    /// `CreateMatrix`'s permutation (`OptimizeBandwidth` over the graph
    /// including the projector couplings).
    fn build_perm_for(
        &self,
        mesh: &dyn MeshView,
        active: &[ElemId],
        params: &ParamMap,
    ) -> Result<Perm, SolverError> {
        let m = active_mask(mesh, active);
        let act = |e: ElemId| m[e.0 as usize];
        let init = initial_permutation(mesh, &act);
        let list = self.graph(mesh, &act, &init);
        let opts = perm_options(params);
        let mut perm = init.perm;
        let n_slots = init.n_dofs;
        if opts.optimize_bandwidth && n_slots > 0 {
            let new_index = rcm_order(&list);
            let initial = bandwidth(&list) + 1;
            let optimized = bandwidth_reordered(&list, &new_index) + 1;
            if !(initial < optimized && !opts.use_optimized_always) {
                for p in perm.iter_mut().filter(|p| **p != NO_DOF) {
                    *p = new_index[*p as usize];
                }
            }
        }
        if let Some(v) = params.string("variable") {
            self.perms
                .lock()
                .unwrap()
                .insert(v.to_string(), perm.clone());
        }
        Ok(Perm { perm, n_slots })
    }

    fn create_matrix(
        &self,
        mesh: &dyn MeshView,
        active: &[ElemId],
        perm: &Perm,
        dofs: usize,
    ) -> CsrMatrix<f64> {
        let m = active_mask(mesh, active);
        let p = AsmPerm::from_perm(perm.perm.clone());
        let list = self.graph(mesh, &|e: ElemId| m[e.0 as usize], &p);
        elmer_assembly::initialize_matrix(&list, None, dofs)
    }

    fn add_local(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        slots: &[u32],
        local: &LocalSystem,
    ) {
        self.inner.add_local(a, rhs, slots, local);
    }

    fn apply_dirichlet(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        fixed: &[(u32, f64)],
        params: &ParamMap,
    ) {
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
        if let Some(var) = params.string("variable") {
            let spec = parse_variable_spec(var, params.integer("variable dofs"))?;
            let names = if spec.components.is_empty() {
                vec![spec.name.clone()]
            } else {
                spec.components
            };
            let bcs = self.bc_views();
            let periodic = names
                .iter()
                .any(|n| bcs.iter().any(|bc| periodic_active(*bc, n)));
            if periodic {
                let perms = self.perms.lock().unwrap();
                let perm = perms.get(var).ok_or_else(|| {
                    SolverError::fatal(
                        "SetPeriodicBoundaries",
                        format!("no permutation recorded for {var}"),
                    )
                })?;
                let n_slots = perm.iter().filter(|&&p| p != NO_DOF).count();
                if n_slots * names.len() != a.nrows {
                    return Err(SolverError::Unsupported(format!(
                        "periodic BCs on a {}-row system for {var} ({n_slots} nodes x {} dofs)",
                        a.nrows,
                        names.len()
                    )));
                }
                for (dof, name) in names.iter().enumerate() {
                    let t = PeriodicTarget {
                        name,
                        dof,
                        ndofs: names.len(),
                        perm,
                    };
                    set_periodic_boundaries(
                        a,
                        rhs,
                        (None, None),
                        &bcs,
                        &self.projectors,
                        &t,
                        self.n_nodes,
                    )?;
                }
            }
        }
        self.inner.try_apply_dirichlet(a, rhs, fixed, params)
    }
}
