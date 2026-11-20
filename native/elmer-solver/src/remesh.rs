//! Moving a running simulation onto a new mesh, as adaptive mesh
//! refinement does (Adaptive.F90 `RefineMesh`): every variable is
//! interpolated from the old mesh (`VariableGet( NewMesh % Variables, Name,
//! ThisOnly = .FALSE. )`, with the mesh projector of
//! `InterpolateMeshToMesh`), then every solver is moved to the new mesh
//! (`UpdateSolverMesh`) keeping its state.
//!
//! The interface:
//!
//! - a solver calls [`SolverCtx::request_mesh_change`] with the new mesh;
//! - [`Simulation::run`] finishes the steady-state iteration on the old
//!   mesh, saves it (`SaveCurrent` of the mesh the solution lives on), then
//!   calls [`Simulation::swap_mesh`] and goes on with the new mesh;
//! - [`SolverCtx::set_adaptive_finished`] marks the mesh converged
//!   (`Mesh % AdaptiveFinished`), which overrides the steady-state
//!   convergence of the solver (SolveEquations, MainUtils.F90:3502);
//! - the runner gives a [`MeshWriterFactory`] that makes the result writer
//!   of each new mesh (each mesh writes its own `<mesh name>/case.result`).
//!
//! Refused for now: transient runs, coupled and block systems, solvers on
//! other than nodal variables, p and DG layouts, conforming assembly, n-t
//! systems and more than one partition (and the runner refuses periodic
//! assembly).

use std::sync::Arc;

use elmer_core::{Communicator, FieldKind, Logger, MeshView, Params, NO_DOF};
use elmer_interp::field::FieldData;
use elmer_interp::transfer::{build_projector, Extras};
use elmer_interp::{InterpOptions, Locator};

use crate::error::SolverError;
use crate::sim::{SaveCurrent, Simulation};
use crate::variable::{Variable, VariableStore};

/// A new mesh a solver asks the run to move to after the current save.
#[derive(Clone)]
pub struct MeshChange {
    pub mesh: Arc<dyn MeshView + Send + Sync>,
    /// `NewMesh % Name`, the directory of its results.
    pub name: String,
}

/// How the variables go to the new mesh. The default interpolates every
/// nodal field (RefineMesh's `VariableGet( NewMesh % Variables, Name,
/// .FALSE. )`).
#[derive(Clone, Default, Debug)]
pub struct MeshTransfer {
    /// Fields left zero on the new mesh (`VariableGet(..., ThisOnly)`):
    /// `PrevNorm` takes the old norm.
    pub zero: Vec<String>,
    /// Fields not carried (their solver makes them on the new mesh).
    pub skip: Vec<String>,
    /// `Adaptive Interpolate = False`: the other single-dof fields are new
    /// zero fields (`VariableAddVector`, output on), `PrevNorm` the old norm.
    pub no_scalar_interp: bool,
}

/// Wraps a new mesh before the run moves to it (the runner's movable
/// coordinates, so that `Displace Mesh` and the mesh mappers move the
/// current mesh).
pub type MeshWrapper = Box<dyn FnMut(Arc<dyn MeshView + Send + Sync>) -> Arc<dyn MeshView + Send + Sync>>;

/// `mesh` with the coordinates `coords` (the old mesh as RefineMesh saw it).
pub(crate) struct WithCoords<'a> {
    pub mesh: &'a dyn MeshView,
    pub coords: &'a [[f64; 3]],
}

impl MeshView for WithCoords<'_> {
    fn dim(&self) -> usize {
        self.mesh.dim()
    }
    fn n_nodes(&self) -> usize {
        self.mesh.n_nodes()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.mesh.n_bulk()
    }
    fn n_boundary(&self) -> usize {
        self.mesh.n_boundary()
    }
    fn element(&self, e: elmer_core::ElemId) -> elmer_core::ElementRef<'_> {
        self.mesh.element(e)
    }
    fn element_edges(&self, e: elmer_core::ElemId) -> Option<&[u32]> {
        self.mesh.element_edges(e)
    }
    fn element_faces(&self, e: elmer_core::ElemId) -> Option<&[u32]> {
        self.mesh.element_faces(e)
    }
}

/// Makes the result writer of a mesh from its name.
pub type MeshWriterFactory =
    Box<dyn FnMut(&str, &Simulation) -> Result<Option<Box<dyn SaveCurrent>>, SolverError>>;

/// Reads an Elmer mesh directory the way the runner reads its first mesh
/// (`LoadMesh2` of `External_ReMesh`).
pub type MeshLoader = Box<dyn FnMut(&std::path::Path) -> Result<Box<dyn MeshView>, SolverError>>;

/// The run's adaptive-mesh state.
#[derive(Default)]
pub struct AdaptiveMeshes {
    pub(crate) pending: Option<MeshChange>,
    /// The old mesh's coordinates when the change was asked for: RefineMesh
    /// interpolates there, before the solver displaces its mesh.
    pub(crate) pending_coords: Option<Vec<[f64; 3]>>,
    /// The variable of the solver that asked for the change.
    pub(crate) pending_var: Option<String>,
    /// `Mesh % AdaptiveFinished` of the current mesh.
    pub(crate) finished: bool,
    /// `Mesh % AdaptiveDepth` of the current mesh.
    pub(crate) depth: usize,
    pub(crate) writer: Option<MeshWriterFactory>,
    pub(crate) wrapper: Option<MeshWrapper>,
    pub(crate) transfer: MeshTransfer,
    pub(crate) loader: Option<MeshLoader>,
    /// Name of the current mesh, empty for the run's first mesh.
    pub(crate) name: String,
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("adaptive mesh change: {what}"))
}

impl Simulation {
    /// Installs the mesh reader of [`SolverCtx::load_mesh`](crate::SolverCtx::load_mesh).
    pub fn set_mesh_loader(&mut self, f: MeshLoader) {
        self.adaptive.loader = Some(f);
    }

    /// Installs the factory of the per-mesh result writers.
    pub fn set_mesh_writer_factory(&mut self, f: MeshWriterFactory) {
        self.adaptive.writer = Some(f);
    }

    /// Installs the wrapper of new meshes ([`MeshWrapper`]).
    pub fn set_mesh_wrapper(&mut self, f: MeshWrapper) {
        self.adaptive.wrapper = Some(f);
    }

    /// The mesh the run moves to: `mesh` through the installed wrapper.
    pub(crate) fn wrap_mesh(&mut self, mesh: Arc<dyn MeshView + Send + Sync>) -> Arc<dyn MeshView + Send + Sync> {
        match self.adaptive.wrapper.as_mut() {
            Some(w) => w(mesh),
            None => mesh,
        }
    }

    /// `Mesh % AdaptiveDepth` of the current mesh (0 for the first).
    pub fn adaptive_depth(&self) -> usize {
        self.adaptive.depth
    }

    /// `Mesh % AdaptiveFinished` of the current mesh.
    pub fn adaptive_finished(&self) -> bool {
        self.adaptive.finished
    }

    /// The pending mesh change, if a solver asked for one.
    pub fn has_pending_mesh(&self) -> bool {
        self.adaptive.pending.is_some()
    }

    /// Checks that the run can move to another mesh at all.
    fn check_swappable(&self, comm: &dyn Communicator) -> Result<(), SolverError> {
        if !self.meshes.is_empty() {
            return Err(unsupported("solver-specific meshes"));
        }
        if comm.size() > 1 {
            return Err(unsupported("more than one partition"));
        }
        if self.is_transient() || self.is_scanning() {
            return Err(unsupported("a transient or scanning run"));
        }
        // Periodic assembly is refused by the runner, which installs it.
        if self.backends.conforming.is_some() {
            return Err(unsupported("conforming assembly"));
        }
        for slot in &self.slots {
            if crate::coupled::is_coupled_system(&slot.params)
                || slot.params.logical("block solver") == Some(true)
            {
                return Err(unsupported("coupled or block systems"));
            }
            if slot.nt.is_some() || !slot.children.is_empty() || slot.mg_meshes.is_some() {
                return Err(unsupported(
                    "n-t systems, child systems or multigrid meshes",
                ));
            }
        }
        for (_, v) in self.vars.iter() {
            match v.kind {
                FieldKind::Nodal | FieldKind::Global => {}
                k => return Err(unsupported(&format!("the {k:?} variable {}", v.name))),
            }
        }
        Ok(())
    }

    /// Moves the run from `old` to `new`: the slots are set up again on the
    /// new mesh with their solver instances and run state kept, and every
    /// nodal variable is interpolated (`InterpolateMeshToMesh` with the mesh
    /// projector). Variables only the old mesh had get the identity
    /// permutation, as RefineMesh's `VariableAdd` without a `Perm` gives
    /// them; their norm (and the solvers') is that of the interpolated
    /// values, as RefineMesh sets it.
    pub fn swap_mesh(
        &mut self,
        old: &dyn MeshView,
        new: &dyn MeshView,
        log: &dyn Logger,
        comm: &dyn Communicator,
    ) -> Result<(), SolverError> {
        self.check_swappable(comm)?;
        let opts = InterpOptions::from_params(&self.model.simulation, None);
        let loc = Locator::new(old, opts)
            .map_err(|e| SolverError::fatal("InterpolateMeshToMesh", e.to_string()))?;
        let fe = self
            .backends
            .fe
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
        let extras = Extras {
            body_mask: None,
            new_mask_perm: None,
            old_dg: None,
            new_dg: None,
        };
        let projector = build_projector(&loc, new, fe.basis.as_ref(), &extras, log)
            .map_err(|e| SolverError::fatal("InterpolateMeshToMesh", e.to_string()))?;

        // The new store keeps the globals; the slots rebuild their fields.
        let old_vars = std::mem::replace(&mut self.vars, VariableStore::new());
        // A solver's own global variable (e.g. an output section's
        // `vtu_internal_dummy`) comes back with its slot; its value is
        // carried after the rebuild.
        let slot_owned: Vec<String> = self
            .slots
            .iter()
            .filter_map(|s| s.variable)
            .map(|id| old_vars.by_id(id).name.clone())
            .collect();
        for (id, v) in old_vars.iter() {
            if v.kind == FieldKind::Global && !slot_owned.contains(&v.name) {
                let comps = old_vars.component_names(id);
                self.vars.add_with_components(v.clone(), &comps)?;
            }
        }
        self.dg_indexes = None;
        self.clear_node_owned();
        self.p_ndeg_parts = None;
        self.p_element_nodes = None;
        let old_slots = std::mem::take(&mut self.slots);
        for mut o in old_slots {
            let Some(mut s) = self.build_slot(o.id - 1, new, None, false)? else {
                continue;
            };
            s.instance = o.instance.take();
            // `Solver % Values` lives on across RefineMesh, with what the
            // solver added to it (FlowSolve's `Nonlinear System Norm Dofs`).
            s.params = std::mem::take(&mut o.params);
            s.times_visited = o.times_visited;
            s.done_time = o.done_time;
            s.dt = o.dt;
            s.state = o.state;
            s.null_norm = o.null_norm;
            s.order = o.order;
            s.started = o.started;
            s.finished = o.finished;
            s.lin_converged = o.lin_converged;
            self.slots.push(s);
        }

        for name in &slot_owned {
            let Some(old) = old_vars.get(name) else { continue };
            if old.kind != FieldKind::Global {
                continue;
            }
            match self.vars.id(name) {
                Some(nid) => {
                    let nv = self.vars.by_id_mut(nid);
                    if nv.values.len() == old.values.len() {
                        nv.values.copy_from_slice(&old.values);
                    }
                }
                None => {
                    let comps = old_vars.component_names(old_vars.id(name).expect("present"));
                    self.vars.add_with_components(old.clone(), &comps)?;
                }
            }
        }
        // Every old nodal field onto the new mesh, component by component.
        let n_new = new.n_nodes();
        let transfer = std::mem::take(&mut self.adaptive.transfer);
        let named = |list: &[String], n: &str| list.iter().any(|x| x.eq_ignore_ascii_case(n));
        for (id, ov) in old_vars.iter() {
            if ov.kind != FieldKind::Nodal || ov.values.len() == ov.dofs {
                continue;
            }
            if named(&transfer.skip, &ov.name) || ov.name.starts_with("coordinate ") {
                continue;
            }
            let zero = named(&transfer.zero, &ov.name);
            if zero || (transfer.no_scalar_interp && ov.dofs == 1) {
                // A zero field on the new mesh, `PrevNorm` the old norm.
                if self.vars.id(&ov.name).is_none() {
                    let comps = old_vars.component_names(id);
                    let mut v = Variable::new(&ov.name, FieldKind::Nodal, ov.dofs, None, vec![0.0; n_new * ov.dofs]);
                    v.output = if zero { ov.output } else { true };
                    v.secondary = ov.secondary;
                    self.vars.add_with_components(v, &comps)?;
                }
                let nid = self.vars.id(&ov.name).expect("added");
                let nv = self.vars.by_id_mut(nid);
                nv.values.iter_mut().for_each(|x| *x = 0.0);
                nv.conv.prev_norm = ov.conv.norm;
                continue;
            }
            if self.vars.id(&ov.name).is_none() {
                let comps = old_vars.component_names(id);
                let perm = match &ov.perm {
                    None => None,
                    Some(p) if p.iter().all(|&k| k != NO_DOF) => None,
                    Some(_) => return Err(unsupported(&format!(
                        "{}: a field on part of the mesh that no solver of the new mesh defines",
                        ov.name
                    ))),
                };
                let mut v = Variable::new(
                    &ov.name,
                    FieldKind::Nodal,
                    ov.dofs,
                    perm,
                    vec![0.0; n_new * ov.dofs],
                );
                v.output = ov.output;
                v.secondary = ov.secondary;
                self.vars.add_with_components(v, &comps)?;
            }
            let nid = self.vars.id(&ov.name).expect("added");
            let nv = self.vars.by_id_mut(nid);
            let dofs = ov.dofs;
            if nv.dofs != dofs {
                return Err(unsupported(&format!(
                    "{}: different dofs on the new mesh",
                    ov.name
                )));
            }
            // The refining solver's vector variable (`Var % DOFs > 1`, RefineMesh's
            // "Interpolate vectors" path) on a mesh an earlier RefineMesh made
            // keeps the identity Perm it was interpolated with, while the
            // solver (UpdateSolverMesh) writes it in its own,
            // bandwidth-optimised order: the next RefineMesh takes the value
            // at equation j for node j. Scalars keep their Perm.
            let refining = self.adaptive.depth > 0
                && ov.dofs > 1
                && self.adaptive.pending_var.as_deref().is_some_and(|n| n.eq_ignore_ascii_case(&ov.name));
            let old_perm: Option<Vec<u32>> = match ov.perm.as_deref() {
                // An explicit identity: the projection still places the
                // values through the new Perm.
                Some(p) if refining && p.iter().all(|&k| k != NO_DOF) => Some((0..p.len() as u32).collect()),
                p => p.map(<[u32]>::to_vec),
            };
            let new_perm: Option<Vec<u32>> = nv.perm.as_deref().map(|p| p.to_vec());
            let nslots_old = ov.values.len() / dofs;
            let nslots_new = nv.values.len() / dofs;
            for k in 0..dofs {
                let ou: Vec<f64> = (0..nslots_old).map(|s| ov.values[dofs * s + k]).collect();
                let old_f = FieldData {
                    name: ov.name.clone(),
                    location: elmer_interp::Location::Nodes,
                    values: ou,
                    perm: old_perm.clone(),
                    prev_values: Vec::new(),
                };
                let mut new_f = vec![FieldData {
                    name: ov.name.clone(),
                    location: elmer_interp::Location::Nodes,
                    values: (0..nslots_new).map(|s| nv.values[dofs * s + k]).collect(),
                    perm: new_perm.clone(),
                    prev_values: Vec::new(),
                }];
                // `Var % Values = 0` before the projection.
                new_f[0].values.iter_mut().for_each(|x| *x = 0.0);
                projector.apply(std::slice::from_ref(&old_f), &mut new_f);
                for (s, x) in new_f[0].values.iter().enumerate() {
                    nv.values[dofs * s + k] = *x;
                }
            }
            // RefineMesh: `NewVar % Norm = SQRT( SUM( Values**2 ) / k )` over
            // the slots with a dof.
            // RefineMesh interpolates the refining solver's variable with the
            // identity permutation, then `UpdateSolverMesh` gives it the
            // solver's (bandwidth-optimised) permutation without moving the
            // values: until the next solve the value at equation `j` is the
            // one of node `j`.
            if self.adaptive.pending_var.as_deref().is_some_and(|n| n.eq_ignore_ascii_case(&ov.name)) {
                if let Some(p) = new_perm.as_deref().filter(|p| p.iter().enumerate().any(|(i, &k)| k as usize != i)) {
                    if p.iter().any(|&k| k == NO_DOF) || p.len() != n_new {
                        return Err(unsupported(&format!(
                            "{}: the refining solver on part of the new mesh with a reordered permutation",
                            ov.name
                        )));
                    }
                    let by_node: Vec<f64> = (0..n_new)
                        .flat_map(|n| (0..dofs).map(move |d| (n, d)))
                        .map(|(n, d)| nv.values[dofs * p[n] as usize + d])
                        .collect();
                    nv.values[..by_node.len()].copy_from_slice(&by_node);
                }
            }
            let k = match &nv.perm {
                Some(p) => p.iter().filter(|&&x| x != NO_DOF).count(),
                None => nv.values.len(),
            };
            let mut s = 0.0;
            for x in &nv.values {
                s += x.powi(2);
            }
            nv.conv.norm = (s / k.max(1) as f64).sqrt();
        }
        self.adaptive.pending_var = None;
        for ix in 0..self.slots.len() {
            if let Some(id) = self.slots[ix].variable {
                let v = self.vars.by_id(id);
                if v.kind != FieldKind::Nodal && v.kind != FieldKind::Global {
                    return Err(unsupported("a solver whose variable is not nodal"));
                }
            }
        }
        Ok(())
    }

    /// The writer of the mesh named `name`, from the installed factory.
    pub(crate) fn mesh_writer(&mut self, name: &str) -> Result<(), SolverError> {
        if let Some(mut f) = self.adaptive.writer.take() {
            let w = f(name, self);
            self.adaptive.writer = Some(f);
            if let Some(w) = w? {
                self.set_save_current(w);
            }
        }
        Ok(())
    }
}
