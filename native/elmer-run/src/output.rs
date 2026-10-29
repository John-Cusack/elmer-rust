//! VTU output with `ResultOutputSolver` semantics (`VtuOutputSolver.F90`)
//! as a scheduled solver: every `ResultOutputSolver` section with VTU
//! format, including the `Post File = "*.vtu"` auto-solver
//! (ElmerSolver.F90 `AddVtuOutputSolverHack`, which `elmer-sif` usually has
//! already appended), becomes a [`VtuSolver`] in its slot, so
//! `elmer-solver` visits it as its `Exec Solver` says (always, after
//! saving, each output step of a transient run). Each visit writes one
//! piece with `elmer-vtu`, numbered as `VtuOutputSolver` numbers them.

use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use elmer_core::{ElemId, ElementRef, FieldKind, FieldRef, FieldSource, MeshView, ParamMap, Params};
use elmer_solver::{
    solver::procedure_of,
    variable::{Snapshot, VarHandle, VarId, VariableStore},
    Model, Simulation, Solver, SolverCtx, SolverError,
};
use elmer_vtu::{post_file_solver, FileIndexer, GeometryIds, VtuNaming, VtuOptions};

/// One VTU output section, checked and with its directory resolved.
#[derive(Clone, Debug)]
pub struct VtuOutput {
    /// Index of the section in `Model::solvers` (solver id - 1).
    pub index: usize,
    opts: VtuOptions,
    dir: PathBuf,
}

pub fn is_result_output(params: &ParamMap) -> bool {
    procedure_of(params).is_some_and(|(lib, entry)| {
        lib.eq_ignore_ascii_case("resultoutputsolve") && entry.eq_ignore_ascii_case("resultoutputsolver")
    })
}

/// Finds the model's VTU outputs: a Simulation `Post File` is first turned
/// into an appended `ResultOutputSolver` section (and removed, as Fortran's
/// `ListRemove` does); then every `ResultOutputSolver` section is checked.
/// `mesh_dir` is the default output directory (`SolverOutputDirectory(...,
/// UseMeshDir = .TRUE.)`), `case_dir` the base of a relative `Output
/// Directory`. Refuses what the composition cannot write correctly, before
/// any side effect.
pub fn prepare_vtu_output(model: &mut Model, case_dir: &Path, mesh_dir: &Path) -> Result<Vec<VtuOutput>, SolverError> {
    if let Some(post) = model.simulation.string("post file").map(str::to_string) {
        match post_file_solver(&model.simulation, &model.solvers) {
            Some(params) => model.solvers.push(params),
            None if post.contains(".vtu") => {} // A VTU ResultOutputSolver already exists.
            None => return Err(SolverError::Unsupported(format!("Post File = {post} (only .vtu is composed)"))),
        }
        model.simulation.0.remove("post file");
    }
    let mut outputs = Vec::new();
    for (index, params) in model.solvers.iter().enumerate() {
        if !is_result_output(params) {
            continue;
        }
        let id = index + 1;
        if params.string("exec solver").is_some_and(|e| e.eq_ignore_ascii_case("never")) {
            continue;
        }
        let vtu = params.logical("vtu format") == Some(true)
            || params.string("output format").is_some_and(|f| f.eq_ignore_ascii_case("vtu"));
        if !vtu {
            return Err(SolverError::Unsupported(format!("Solver {id}: ResultOutputSolver other than VTU")));
        }
        let mut opts = VtuOptions::from_params(params);
        opts.apply_simulation(&model.simulation);
        if opts.time_collection || opts.part_collection {
            opts.unsupported.push("vtu time/part collection (.pvd)".into());
        }
        if !opts.unsupported.is_empty() {
            return Err(SolverError::Unsupported(format!("Solver {id}: VTU output: {}", opts.unsupported.join(", "))));
        }
        if opts.output_file_name.contains('/') {
            return Err(SolverError::Unsupported(format!("Solver {id}: a directory in the VTU output file name")));
        }
        if opts.save_geometry_ids {
            opts = opts.with_geometry(geometry(model));
        }
        let dir = match params
            .string("output directory")
            .or_else(|| model.simulation.string("output directory"))
        {
            Some(d) => case_dir.join(d),
            None => mesh_dir.to_path_buf(),
        };
        outputs.push(VtuOutput { index, opts, dir });
    }
    Ok(outputs)
}

/// Body and BC facts for `GeometryIds`. The writer sees [`MappedTags`], so
/// BC `i` has tag `i` (`Model % BCs(i) % Tag = i`, ModelDescription.F90:658).
fn geometry(model: &Model) -> GeometryIds {
    let ids = |secs: &[ParamMap]| secs.iter().map(|s| s.integer("geometry id")).collect();
    let names = |secs: &[ParamMap]| secs.iter().map(|s| s.string("name").map(str::to_string)).collect();
    GeometryIds {
        n_bodies: model.bodies.len(),
        bc_tags: (1..=model.boundary_conditions.len() as i32).collect(),
        body_geometry_id: ids(&model.bodies),
        bc_geometry_id: ids(&model.boundary_conditions),
        body_names: names(&model.bodies),
        bc_names: names(&model.boundary_conditions),
    }
}

/// The mesh as Fortran's writers see it after `MapBodiesAndBCs`
/// (MeshLoad.F90:330-455): a boundary element's tag is the BC that targets
/// it, or 0 when none does; a bulk element's tag is its body. The solver
/// keeps the file tags and maps them on lookup, so only the writer needs
/// this view.
pub(crate) struct MappedTags<'a> {
    mesh: &'a dyn MeshView,
    model: &'a Model,
}

impl<'a> MappedTags<'a> {
    pub(crate) fn new(mesh: &'a dyn MeshView, model: &'a Model) -> Self {
        MappedTags { mesh, model }
    }
}

impl MeshView for MappedTags<'_> {
    fn dim(&self) -> usize {
        self.mesh.dim()
    }
    fn n_nodes(&self) -> usize {
        self.mesh.n_nodes()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.mesh.coord(node)
    }
    fn n_bulk(&self) -> usize {
        self.mesh.n_bulk()
    }
    fn n_boundary(&self) -> usize {
        self.mesh.n_boundary()
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let mut el = self.mesh.element(e);
        el.tag = if (e.0 as usize) < self.mesh.n_bulk() {
            self.model.body_of_tag(el.tag).map_or(el.tag, |b| b as i32)
        } else {
            self.model.bc_of_tag(el.tag).map_or(0, |b| b as i32)
        };
        el
    }
    fn element_edges(&self, e: ElemId) -> Option<&[u32]> {
        self.mesh.element_edges(e)
    }
    fn element_faces(&self, e: ElemId) -> Option<&[u32]> {
        self.mesh.element_faces(e)
    }
}

impl VtuOutput {
    /// `CreateListForSaving` (Lists.F90:10031) when the section names no
    /// `Scalar Field 1` / `Vector Field 1`: every output variable of the
    /// mesh's variable list ([`ElmerVariables`]) in order. Multi-dof
    /// variables are saved through their components; components that make
    /// up a `dim`- or 3-vector (`velocity 1..3`) become one vector field.
    /// `Displacement` and `Mesh Update` are vectors, and the `Mesh Update`
    /// complement of `Displacement` is refused (`elmer-vtu` has no `Vector
    /// Field n Complement`). Call after `Simulation::setup`, before the run;
    /// `dim` is the mesh dimension (`Mesh % MeshDim`).
    pub fn select_fields(&mut self, sim: &Simulation, dim: usize) -> Result<(), SolverError> {
        if !self.opts.scalar_fields.is_empty() || !self.opts.vector_fields.is_empty() {
            return Ok(());
        }
        let enforce_vectors = sim.model.solvers[self.index].logical("enforce vectors").unwrap_or(true);
        let (scalars, vectors) = fields_for_saving(&sim.vars, enforce_vectors, dim)?;
        self.opts.scalar_fields = scalars;
        self.opts.vector_fields = vectors;
        Ok(())
    }

    /// The solver to install in this output's slot; it appends each file it
    /// writes to `written`.
    pub fn into_solver(self, written: Arc<Mutex<Vec<PathBuf>>>) -> VtuSolver {
        VtuSolver { out: self, indexer: FileIndexer::new(), written }
    }
}

/// `VtuOutputSolver` for one section: one `.vtu` piece per visit.
pub struct VtuSolver {
    out: VtuOutput,
    indexer: FileIndexer,
    written: Arc<Mutex<Vec<PathBuf>>>,
}

impl Solver for VtuSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let opts = &self.out.opts;
        let index = self.indexer.advance(opts);
        let naming = VtuNaming { no_fileindex: opts.no_fileindex, ..VtuNaming::default() };
        let name = naming.file_name(&opts.output_file_name, ".vtu", index, None, false);
        std::fs::create_dir_all(&self.out.dir)?;
        let path = self.out.dir.join(name);
        let mapped = MappedTags { mesh: ctx.mesh, model: &ctx.sim.model };
        elmer_vtu::write_vtu(&path, &mapped, &ElmerVariables::new(&ctx.sim.vars), opts)?;
        self.written.lock().unwrap().push(path);
        Ok(())
    }
}

/// [`VtuOutput::select_fields`]'s lists: `(Scalar Field n, Vector Field n)`.
fn fields_for_saving(
    store: &VariableStore,
    enforce_vectors: bool,
    dim: usize,
) -> Result<(Vec<String>, Vec<String>), SolverError> {
    let exists = |name: &str| store.resolve(name).is_some();
    let variables = ElmerVariables::new(store);
    let (mut scalars, mut vectors) = (Vec::new(), Vec::new());
    let fields = variables.fields();
    let displacement_v = fields.iter().any(|f| f.name == "displacement");
    for f in &fields {
        if f.kind == FieldKind::Global || f.values.len() == f.dofs {
            continue;
        }
        let (mut name, mut is_vector) = (f.name.to_string(), f.dofs > 1);
        let set = match f.name {
            "coordinate 1" | "coordinate 2" | "coordinate 3" => false,
            "mesh update 1" | "mesh update 2" | "mesh update 3" => false,
            "mesh update" => !exists("displacement"),
            "displacement" => {
                if exists("mesh update") {
                    return Err(SolverError::Unsupported(
                        "VTU output of Displacement with its Mesh Update complement".into(),
                    ));
                }
                true
            }
            _ if f.dofs > 1 => false,
            _ => {
                let mut set = true;
                // A trailing component index: `velocity 2`.
                let index = f.name.rsplit_once(' ').filter(|(_, t)| t.bytes().all(|b| b.is_ascii_digit()));
                let comp = index.and_then(|(_, t)| t.parse::<usize>().ok()).unwrap_or(0);
                let base = index.map_or(f.name, |(b, _)| b);
                let parent = index.and_then(|(b, _)| store.resolve(b)).map(|_| store.get(base));
                if let Some(parent) = parent {
                    let dofs = parent.map_or(1, |v| v.dofs);
                    is_vector = dofs == dim || dofs == 3;
                    set = comp == 1 || !is_vector;
                } else if enforce_vectors && comp == 1 {
                    if exists(&format!("{base} {dim}")) {
                        is_vector = !exists(&format!("{base} 4"));
                    }
                } else if enforce_vectors && (comp == 2 || comp == 3) && exists(&format!("{base} 1")) {
                    set = exists(&format!("{base} 4"));
                    if !set && comp == 2 && dim == 3 {
                        set = !exists(&format!("{base} {dim}"));
                    }
                }
                if set && is_vector {
                    name = base.to_string();
                }
                set && !(name == "displacement" && displacement_v)
            }
        };
        if set {
            if is_vector {
                vectors.push(name);
            } else {
                scalars.push(name);
            }
        }
    }
    Ok((scalars, vectors))
}


/// Elmer's `Mesh % Variables` list as the writer sees it: each output
/// variable followed by its components, which Elmer keeps as variables of
/// their own (`VariableAddVector`: `pot` and `temperature`, `potential` for
/// `POT[Temperature:1 Potential:1]`; `displacement 1..n` otherwise).
/// `elmer-solver` keeps components as aliases, so their values are copied
/// out of the parent here.
struct ElmerVariables<'a> {
    store: &'a VariableStore,
    /// `(name, parent, values)` of each component.
    components: Vec<(String, VarId, Vec<f64>)>,
}

impl<'a> ElmerVariables<'a> {
    fn new(store: &'a VariableStore) -> Self {
        let components = store
            .names_in_order()
            .into_iter()
            .filter_map(|(name, h)| match h {
                VarHandle::Component { parent, component } => {
                    let v = store.by_id(parent);
                    let values = v.values.iter().skip(component).step_by(v.dofs).copied().collect();
                    Some((name, parent, values))
                }
                VarHandle::Primary(_) => None,
            })
            .collect();
        ElmerVariables { store, components }
    }
}

impl FieldSource for ElmerVariables<'_> {
    fn fields(&self) -> Vec<FieldRef<'_>> {
        let mut out = Vec::new();
        for (id, v) in self.store.iter() {
            if !v.output {
                continue;
            }
            out.push(v.field_ref());
            out.extend(self.components.iter().filter(|(_, p, _)| *p == id).map(|(name, _, values)| FieldRef {
                name,
                kind: v.kind,
                dofs: 1,
                perm: v.perm.as_deref(),
                values,
            }));
        }
        out
    }

    fn time(&self) -> f64 {
        Snapshot { store: self.store }.time()
    }

    fn step(&self) -> usize {
        Snapshot { store: self.store }.step()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_solver::variable::Variable;

    fn nodal(name: &str, dofs: usize, values: Vec<f64>) -> Variable {
        Variable::new(name, FieldKind::Nodal, dofs, None, values)
    }

    /// `POT[Temperature:1 Potential:1]` is saved through its two named
    /// components; a 2-vector through `velocity 1` as one vector field.
    #[test]
    fn components_are_fields_and_vectors_are_merged() {
        let mut store = VariableStore::new();
        store
            .add_with_components(nodal("pot", 2, vec![1.0, 10.0, 2.0, 20.0]), &["temperature".into(), "potential".into()])
            .unwrap();
        store.add(nodal("velocity", 2, vec![0.0; 4])).unwrap();
        store.add(nodal("coordinate 1", 1, vec![0.0, 1.0])).unwrap();
        store.add(Variable::new("time", FieldKind::Global, 1, None, vec![0.0])).unwrap();

        let (scalars, vectors) = fields_for_saving(&store, true, 2).unwrap();
        assert_eq!(scalars, ["temperature", "potential"]);
        assert_eq!(vectors, ["velocity"]);

        let variables = ElmerVariables::new(&store);
        let fields = variables.fields();
        let potential = fields.iter().find(|f| f.name == "potential").unwrap();
        assert_eq!(potential.values, [10.0, 20.0]);
    }

    #[test]
    fn displacement_with_mesh_update_complement_is_refused() {
        let mut store = VariableStore::new();
        store.add(nodal("displacement", 2, vec![0.0; 4])).unwrap();
        store.add(nodal("mesh update", 2, vec![0.0; 4])).unwrap();
        assert!(matches!(fields_for_saving(&store, true, 2), Err(SolverError::Unsupported(_))));
    }
}
