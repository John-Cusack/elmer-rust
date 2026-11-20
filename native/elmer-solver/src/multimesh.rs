//! Solver-specific meshes (`Solver i :: Mesh = "<dir>"`): solvers that run
//! on meshes of their own, side by side with the run's first mesh, as
//! Elmer's `Model % Meshes` list.
//!
//! Each mesh has its own variable list (`Mesh % Variables`), so the
//! [`VariableStore`] of a slot's mesh is the one in [`Simulation::vars`]
//! while that slot runs, and the others are parked; the global variables
//! (`Time`, `Timestep`, ...) are shared, as AddTimeEtc points every mesh's
//! globals at the same storage, so they move with each switch.
//!
//! A field read on a mesh that does not have it (`globalT = Equals localT`
//! on the mesh of `globalT`) is `VariableGet( Mesh % Variables, Name,
//! ThisOnly = .FALSE. )`: a copy of the field on the reading mesh, with the
//! identity permutation and the values interpolated from the mesh the field
//! lives on (`InterpolateMeshToMesh`; nodes not found keep zero). The copy is
//! made when the case is set up (the Dirichlet initialisation reads it
//! first) and interpolated again before a solver of its mesh runs once the
//! field's own solver has changed it (`Var % Valid`). The copies the
//! reading mesh needs are those its solvers' Dirichlet conditions depend
//! on; a cross-mesh dependency elsewhere is not found, as before.
//!
//! The same mesh projector and interpolation as the mesh change of
//! [`crate::remesh`] (`elmer-interp`).

use std::collections::HashMap;
use std::sync::Arc;

use elmer_core::{FieldKind, Logger, MeshView, Params, Value};
use elmer_interp::field::FieldData;
use elmer_interp::transfer::{build_projector, Extras, Projector};
use elmer_interp::{InterpOptions, Locator};

use crate::error::SolverError;
use crate::sim::{SaveCurrent, Simulation};
use crate::variable::{VarHandle, Variable, VariableStore};

/// A mesh of the solvers that name it (`Mesh % Name`).
pub struct SolverMesh {
    pub name: String,
    pub mesh: Arc<dyn MeshView>,
}

/// An interpolated copy of a field on a mesh that reads it.
struct Copy {
    /// The reading mesh (0 is the run's first mesh).
    mesh: usize,
    name: String,
    /// The mesh the field lives on.
    source: usize,
    stale: bool,
}

/// The run's solver-specific meshes and their parked variable lists.
#[derive(Default)]
pub(crate) struct Meshes {
    /// Mesh `m` is `list[m - 1]`.
    list: Vec<SolverMesh>,
    /// The variable list of each mesh but the active one, by mesh index.
    parked: Vec<Option<VariableStore>>,
    pub(crate) active: usize,
    copies: Vec<Copy>,
    /// The `Output File` writer of each solver mesh, by `m - 1`.
    writers: Vec<Option<Box<dyn SaveCurrent>>>,
    projectors: HashMap<(usize, usize), Projector>,
}

/// A slot's mesh: the run's first mesh, borrowed, or a solver mesh.
pub(crate) enum MeshRef<'a> {
    First(&'a dyn MeshView),
    Solver(Arc<dyn MeshView>),
}

impl MeshRef<'_> {
    pub(crate) fn get(&self) -> &dyn MeshView {
        match self {
            MeshRef::First(m) => *m,
            MeshRef::Solver(m) => &**m,
        }
    }
}

impl Meshes {
    pub(crate) fn is_empty(&self) -> bool {
        self.list.is_empty()
    }
}

impl Simulation {
    /// Adds the mesh named `name` (a solver's `Mesh` keyword) for solvers
    /// on it; call before [`setup`](Simulation::setup). Returns its index.
    pub fn add_solver_mesh(&mut self, name: &str, mesh: Arc<dyn MeshView>) -> usize {
        let mm = &mut self.meshes;
        if mm.parked.is_empty() {
            mm.parked.push(None);
        }
        mm.list.push(SolverMesh { name: name.to_string(), mesh });
        mm.parked.push(Some(VariableStore::new()));
        mm.writers.push(None);
        mm.list.len()
    }

    /// The solver meshes, in the order they were added.
    pub fn solver_meshes(&self) -> &[SolverMesh] {
        &self.meshes.list
    }

    /// Installs the `Output File` writer of solver mesh `m` (1-based).
    pub fn set_solver_mesh_writer(&mut self, m: usize, writer: Box<dyn SaveCurrent>) {
        if let Some(w) = self.meshes.writers.get_mut(m.wrapping_sub(1)) {
            *w = Some(writer);
        }
    }

    /// The mesh index of a solver section: 0 without `Mesh`, else the
    /// solver mesh of that name.
    pub(crate) fn mesh_index_of(&self, params: &elmer_core::ParamMap) -> Result<usize, SolverError> {
        let Some(name) = params.string("mesh") else { return Ok(0) };
        let name = name.trim();
        self.meshes
            .list
            .iter()
            .position(|m| m.name == name)
            .map(|i| i + 1)
            .ok_or_else(|| SolverError::Unsupported(format!("Solver Mesh \"{name}\" was not loaded")))
    }

    /// Mesh `m`, with `first` the run's first mesh.
    pub(crate) fn mesh_ref<'a>(&self, m: usize, first: &'a dyn MeshView) -> MeshRef<'a> {
        if m == 0 {
            MeshRef::First(first)
        } else {
            MeshRef::Solver(Arc::clone(&self.meshes.list[m - 1].mesh))
        }
    }

    /// Makes mesh `m`'s variable list the active one ([`Simulation::vars`]),
    /// carrying the global variables over.
    pub(crate) fn enter_mesh(&mut self, m: usize) {
        let a = self.meshes.active;
        if a == m || self.meshes.list.is_empty() {
            return;
        }
        let globals: Vec<(String, f64)> = self
            .vars
            .iter()
            .filter(|(_, v)| v.kind == FieldKind::Global)
            .filter_map(|(_, v)| v.values.first().map(|&x| (v.name.clone(), x)))
            .collect();
        let next = self.meshes.parked[m].take().expect("a parked variable list");
        self.meshes.parked[a] = Some(std::mem::replace(&mut self.vars, next));
        for (name, x) in globals {
            self.vars.set_scalar(&name, x);
        }
        self.meshes.active = m;
    }

    /// The copies mesh `m` needs: the variables its solvers' Dirichlet
    /// conditions depend on that live on another mesh only. Made with zero
    /// values on the run's first visit to mesh `m`'s Dirichlet values.
    pub(crate) fn make_mesh_copies(&mut self, m: usize, mesh: &dyn MeshView) -> Result<(), SolverError> {
        let own: Vec<String> = self
            .slots
            .iter()
            .filter(|s| s.mesh_ix == m)
            .filter_map(|s| s.params.string("variable").map(|v| crate::sim::parse_variable_spec(v, s.params.integer("variable dofs"))))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flat_map(|spec| {
                let mut v = vec![spec.name.clone()];
                v.extend(spec.components);
                v
            })
            .collect();
        let mut wanted: Vec<String> = Vec::new();
        // The fields an output solver of this mesh names (`Scalar Field n`,
        // `Vector Field n` of ResultOutputSolver, `Variable n` of SaveData).
        for slot in self.slots.iter().filter(|s| s.mesh_ix == m) {
            for (k, v) in &slot.params.0 {
                let mut w = k.rsplitn(2, ' ');
                let (Some(num), Some(head)) = (w.next(), w.next()) else { continue };
                if !num.bytes().all(|c| c.is_ascii_digit()) || !["scalar field", "vector field", "variable"].contains(&head) {
                    continue;
                }
                if let Value::Str(name) | Value::File(name) = v {
                    let name = elmer_core::normalize_key(name.trim_matches('"'));
                    if !wanted.contains(&name) {
                        wanted.push(name);
                    }
                }
            }
        }
        for bc in &self.model.boundary_conditions {
            for (k, v) in &bc.0 {
                if !own.iter().any(|o| o == k) {
                    continue;
                }
                if let Value::Dependent(d) = v {
                    for dep in &d.variables {
                        let dep = elmer_core::normalize_key(dep);
                        if !wanted.contains(&dep) {
                            wanted.push(dep);
                        }
                    }
                }
            }
        }
        for name in wanted {
            if self.vars.resolve(&name).is_some() {
                continue;
            }
            // The mesh the field lives on: a solver's own variable there.
            let source = self.slots.iter().find_map(|s| {
                let spec = crate::sim::parse_variable_spec(s.params.string("variable")?, s.params.integer("variable dofs")).ok()?;
                (spec.name == name || spec.components.contains(&name)).then_some(s.mesh_ix)
            });
            let Some(source) = source.filter(|&s| s != m) else { continue };
            let dofs = {
                let store = self.meshes.parked[source].as_ref().expect("parked");
                match store.resolve(&name) {
                    Some(VarHandle::Primary(id)) => store.by_id(id).dofs,
                    Some(VarHandle::Component { .. }) => 1,
                    None => continue,
                }
            };
            if dofs != 1 {
                return Err(SolverError::Unsupported(format!(
                    "solver meshes: an interpolated copy of {name} with {dofs} components"
                )));
            }
            let n = mesh.n_nodes();
            let perm: Arc<[u32]> = (0..n as u32).collect::<Vec<_>>().into();
            self.vars.add(Variable::new(&name, FieldKind::Nodal, 1, Some(perm), vec![0.0; n]))?;
            self.meshes.copies.push(Copy { mesh: m, name, source, stale: true });
        }
        Ok(())
    }

    /// Interpolates the stale copies on the active mesh (`mesh`) from the
    /// meshes their fields live on; `first` is the run's first mesh.
    pub(crate) fn refresh_mesh_copies(&mut self, mesh: &dyn MeshView, first: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        let m = self.meshes.active;
        for i in 0..self.meshes.copies.len() {
            let c = &self.meshes.copies[i];
            if c.mesh != m || !c.stale {
                continue;
            }
            let (name, source) = (c.name.clone(), c.source);
            if !self.meshes.projectors.contains_key(&(source, m)) {
                let held = (source > 0).then(|| Arc::clone(&self.meshes.list[source - 1].mesh));
                let src: &dyn MeshView = match &held {
                    Some(m) => &**m,
                    None => first,
                };
                let opts = InterpOptions::from_params(&self.model.simulation, None);
                let loc = Locator::new(src, opts).map_err(|e| SolverError::fatal("InterpolateMeshToMesh", e.to_string()))?;
                let fe = self.backends.fe.as_ref().ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
                let extras = Extras { body_mask: None, new_mask_perm: None, old_dg: None, new_dg: None };
                let p = build_projector(&loc, mesh, fe.basis.as_ref(), &extras, log)
                    .map_err(|e| SolverError::fatal("InterpolateMeshToMesh", e.to_string()))?;
                self.meshes.projectors.insert((source, m), p);
            }
            let old = {
                let store = self.meshes.parked[source].as_ref().expect("parked");
                let (id, offset, stride) = match store.resolve(&name) {
                    Some(VarHandle::Primary(id)) => (id, 0, 1),
                    Some(VarHandle::Component { parent, component }) => (parent, component, store.by_id(parent).dofs),
                    None => continue,
                };
                let v = store.by_id(id);
                let slots = v.values.len() / stride;
                FieldData {
                    name: name.clone(),
                    location: elmer_interp::Location::Nodes,
                    values: (0..slots).map(|s| v.values[stride * s + offset]).collect(),
                    perm: v.perm.as_deref().map(|p| p.to_vec()),
                    prev_values: Vec::new(),
                }
            };
            let id = self.vars.id(&name).expect("the copy");
            let nv = self.vars.by_id_mut(id);
            let mut new = vec![FieldData {
                name: name.clone(),
                location: elmer_interp::Location::Nodes,
                values: vec![0.0; nv.values.len()],
                perm: nv.perm.as_deref().map(|p| p.to_vec()),
                prev_values: Vec::new(),
            }];
            self.meshes.projectors[&(source, m)].apply(std::slice::from_ref(&old), &mut new);
            nv.values.copy_from_slice(&new[0].values);
            self.meshes.copies[i].stale = false;
        }
        Ok(())
    }

    /// After slot `ix` ran: its field's copies on other meshes are stale
    /// (`InvalidateVariable`).
    pub(crate) fn invalidate_mesh_copies(&mut self, ix: usize) {
        let m = self.slots[ix].mesh_ix;
        let Some(id) = self.slots[ix].variable else { return };
        let name = self.vars.by_id(id).name.clone();
        let comps = self.vars.component_names(id);
        for c in &mut self.meshes.copies {
            if c.source == m && (c.name == name || comps.contains(&c.name)) {
                c.stale = true;
            }
        }
    }

    /// The equation (`Var % Solver`'s) a field of the active mesh belongs
    /// to, for the `Output File` header: a slot of this mesh, or for a copy
    /// the slot of the field's own mesh.
    pub fn equation_of(&self, id: crate::variable::VarId) -> Option<String> {
        let m = self.meshes.active;
        let name = &self.vars.by_id(id).name;
        let slot = match self.meshes.copies.iter().find(|c| c.mesh == m && &c.name == name) {
            Some(c) => self.slots.iter().find(|s| {
                s.mesh_ix == c.source
                    && s.params.string("variable").is_some_and(|v| {
                        crate::sim::parse_variable_spec(v, s.params.integer("variable dofs")).is_ok_and(|sp| sp.name == c.name)
                    })
            }),
            None => self.slots.iter().find(|s| s.mesh_ix == m && s.variable == Some(id)),
        };
        slot.and_then(|s| self.model.solvers.get(s.id - 1)).and_then(|p| p.string("equation")).map(elmer_core::normalize_key)
    }

    /// The variable list of slot `ix`'s mesh.
    pub fn slot_vars(&self, ix: usize) -> &VariableStore {
        let m = self.slots[ix].mesh_ix;
        if m == self.meshes.active {
            &self.vars
        } else {
            self.meshes.parked[m].as_ref().expect("a parked variable list")
        }
    }

    /// The norm of slot `ix`'s variable, from its mesh's variable list.
    pub fn slot_norm(&self, ix: usize) -> Option<f64> {
        let slot = &self.slots[ix];
        let id = slot.variable?;
        let store = if slot.mesh_ix == self.meshes.active { &self.vars } else { self.meshes.parked[slot.mesh_ix].as_ref()? };
        Some(store.by_id(id).conv.norm)
    }

    /// `SaveCurrent` of every mesh: the first mesh's writer, then each
    /// solver mesh's with its variable list.
    pub(crate) fn save_meshes(&mut self, first: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        if let Some(mut writer) = self.save_current.take() {
            self.enter_mesh(0);
            let saved = writer.save(self, first, log);
            self.save_current = Some(writer);
            saved?;
        }
        for m in 1..=self.meshes.list.len() {
            let Some(mut writer) = self.meshes.writers[m - 1].take() else { continue };
            self.enter_mesh(m);
            let mesh = Arc::clone(&self.meshes.list[m - 1].mesh);
            let saved = writer.save(self, &*mesh, log);
            self.meshes.writers[m - 1] = Some(writer);
            self.enter_mesh(0);
            saved?;
        }
        Ok(())
    }
}
