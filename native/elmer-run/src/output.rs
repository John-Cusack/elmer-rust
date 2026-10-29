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

use elmer_core::{ElemId, ElementRef, FieldKind, FieldSource, MeshView, ParamMap, Params};
use elmer_solver::{
    solver::procedure_of, variable::Snapshot, Model, Simulation, Solver, SolverCtx, SolverError,
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
    /// `Scalar Field 1` / `Vector Field 1`: every output variable in
    /// registration order. Only the plain scalar branch is ported; vector and
    /// component variables (`dofs > 1`, names ending in an index,
    /// `displacement`/`mesh update` pairing) are refused. Call after
    /// `Simulation::setup`, before the run.
    pub fn select_fields(&mut self, sim: &Simulation) -> Result<(), SolverError> {
        if !self.opts.scalar_fields.is_empty() || !self.opts.vector_fields.is_empty() {
            return Ok(());
        }
        let snapshot = Snapshot { store: &sim.vars };
        for f in snapshot.fields() {
            if f.kind == FieldKind::Global || f.values.len() == f.dofs {
                continue;
            }
            let indexed = f.name.rsplit_once(' ').is_some_and(|(_, t)| t.bytes().all(|b| b.is_ascii_digit()));
            if f.name.starts_with("coordinate ") && indexed {
                continue;
            }
            if f.dofs > 1 || indexed || matches!(f.name, "displacement" | "mesh update") {
                return Err(SolverError::Unsupported(format!(
                    "VTU output of vector or component variable {} without Scalar/Vector Field keywords",
                    f.name
                )));
            }
            self.opts.scalar_fields.push(f.name.to_string());
        }
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
        elmer_vtu::write_vtu(&path, &mapped, &Snapshot { store: &ctx.sim.vars }, opts)?;
        self.written.lock().unwrap().push(path);
        Ok(())
    }
}
