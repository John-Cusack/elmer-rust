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
    /// `Reference Values` (first column): the visit's norm is then
    /// `AscBinCompareNorm` of the written piece against them.
    reference: Option<[f64; 6]>,
    /// `Eigen Analysis` / `Constraint Modes Analysis` of the section: one
    /// piece per eigen or constraint mode ([`VtuSolver::write_mode_files`]).
    mode_files: bool,
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
/// Directory`. `parallel` is whether the run has several ranks. Refuses
/// what the composition cannot write correctly, before any side effect.
pub fn prepare_vtu_output(
    model: &mut Model,
    case_dir: &Path,
    mesh_dir: &Path,
    parallel: bool,
) -> Result<Vec<VtuOutput>, SolverError> {
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
        // `Skip Halo Elements` / `Save Halo Elements Only` are read only in a
        // parallel run; serially both are off (SaveUtils.F90:373-379).
        if !parallel {
            opts.unsupported.retain(|u| u != "skip halo elements" && u != "save halo elements only");
        }
        if opts.time_collection || opts.part_collection {
            opts.unsupported.push("vtu time/part collection (.pvd)".into());
        }
        // ResultOutputSolver.F90:349: the solver's norm becomes
        // `AscBinCompareNorm(Reference Values)`, a pseudonorm of the written
        // data ([`asc_bin_sums`]), reproduced for one binary serial piece.
        let reference = match params.get("reference values") {
            None => None,
            Some(v) => match reference_values(v) {
                Some(r) if !opts.ascii && !parallel => Some(r),
                Some(_) => {
                    opts.unsupported.push("Reference Values with ASCII or parallel output".into());
                    None
                }
                None => {
                    opts.unsupported.push("Reference Values with fewer than 6 values".into());
                    None
                }
            },
        };
        // VtuOutputSolver.F90:561: each mode to its own file. `Ignore
        // Eigenmodes` turns it off (and the mode fields with it).
        let ignore_modes = params.logical("ignore eigenmodes") == Some(true);
        let mode_files = !ignore_modes
            && (params.logical("eigen analysis") == Some(true) || params.logical("constraint modes analysis") == Some(true));
        if mode_files || ignore_modes {
            opts.unsupported.retain(|u| u != "eigen analysis" && u != "constraint modes analysis");
        }
        if mode_files {
            if reference.is_some() || parallel {
                opts.unsupported.push("eigen or constraint mode files with Reference Values or in parallel".into());
            }
            if params.string("eigen vector component").is_some_and(|c| !c.eq_ignore_ascii_case("re")) {
                opts.unsupported.push("Eigen Vector Component other than re".into());
            }
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
        outputs.push(VtuOutput { index, opts, dir, reference, mode_files });
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
        VtuSolver { out: self, indexer: FileIndexer::new(), written, part: 0, parts: 1 }
    }
}

/// `VtuOutputSolver` for one section: one `.vtu` piece per visit.
pub struct VtuSolver {
    out: VtuOutput,
    indexer: FileIndexer,
    written: Arc<Mutex<Vec<PathBuf>>>,
    /// This rank and the number of ranks writing pieces (`Part`, `ParEnv % PEs`).
    part: usize,
    parts: usize,
}

impl VtuSolver {
    /// The pieces of rank `part` of `parts` ranks that each hold a whole
    /// mesh (`Single Mesh`): `<base>_<np>np<part+1>.vtu`, and rank 0 also
    /// writes the `.pvtu` listing every piece (`WritePvtuFile` with all
    /// partitions active). One rank keeps the serial names.
    pub fn in_job(mut self, part: usize, parts: usize) -> Self {
        self.part = part;
        self.parts = parts;
        self
    }
}

/// `MaxModes` of VtuOutputSolver (F90:506-528): the section's `Eigen System
/// Values`, else the largest `Eigen System Values` / `Harmonic System
/// Values` (or `Scanning Loops` with `Save Scanning Modes`) of any solver.
fn max_modes(model: &Model, p: &ParamMap) -> usize {
    if let Some(n) = p.integer("eigen system values") {
        return n.max(0) as usize;
    }
    let mut m = 0;
    for s in &model.solvers {
        m = m.max(s.integer("eigen system values").unwrap_or(0));
        m = m.max(s.integer("harmonic system values").unwrap_or(0));
        if s.logical("save scanning modes") == Some(true) {
            m = m.max(s.integer("scanning loops").unwrap_or(0));
        }
    }
    m.max(0) as usize
}

/// `MaxModes2`: the largest `NumberOfConstraintModes` of a solver variable.
fn max_constraint_modes(sim: &Simulation) -> usize {
    sim.slots
        .iter()
        .filter_map(|s| s.variable)
        .map(|id| sim.vars.by_id(id).constraint_modes.len())
        .max()
        .unwrap_or(0)
}

/// The primary variable a `Scalar Field` / `Vector Field` name saves.
fn field_variable(store: &VariableStore, name: &str) -> Option<VarId> {
    let of = |h: VarHandle| match h {
        VarHandle::Primary(id) => id,
        VarHandle::Component { parent, .. } => parent,
    };
    store.resolve(name).or_else(|| store.resolve(&format!("{name} 1"))).map(of)
}

impl VtuSolver {
    /// `Eigen Analysis` / `Constraint Modes Analysis` output
    /// (VtuOutputSolver.F90:561, :1028-1062): `FileIndex` 1 to `MaxModes`
    /// holds eigenmode `FileIndex` of each saved variable that has
    /// eigenvectors, the following `MaxModes2` files constraint mode
    /// `FileIndex - MaxModes` of each one that has constraint modes, under
    /// the variable's own name. A variable without such a mode is not in
    /// that file; one with too few is left out with Elmer's warning.
    fn write_mode_files(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        const CALLER: &str = "VtuOutputSolver";
        ctx.info(CALLER, "Saving each mode to different file", 4);
        let params = &ctx.sim.model.solvers[self.out.index];
        let mm = max_modes(&ctx.sim.model, params);
        let mm2 = max_constraint_modes(ctx.sim);
        let opts = &self.out.opts;
        let naming = VtuNaming { pes: 1, no_fileindex: opts.no_fileindex, parallel_base: opts.partition_numbering, ..VtuNaming::default() };
        std::fs::create_dir_all(&self.out.dir)?;
        let store = &ctx.sim.vars;
        let mut warnings = Vec::new();
        for file_index in 1..=mm + mm2 {
            let mut over: std::collections::HashMap<VarId, Vec<f64>> = std::collections::HashMap::new();
            let mut keep = |names: &[String]| -> Vec<String> {
                names
                    .iter()
                    .filter(|name| {
                        let Some(id) = field_variable(store, name) else { return false };
                        if over.contains_key(&id) {
                            return true;
                        }
                        let v = store.by_id(id);
                        let (modes, ind, what): (Vec<Vec<f64>>, usize, &str) = if file_index <= mm {
                            if v.eigen_vectors.is_empty() {
                                return false;
                            }
                            (v.eigen_vectors.iter().map(|m| m.iter().map(|z| z.re).collect()).collect(), file_index, "eigenmodes")
                        } else {
                            if v.constraint_modes.is_empty() {
                                return false;
                            }
                            (v.constraint_modes.clone(), file_index - mm, "constraint modes")
                        };
                        if ind > modes.len() {
                            warnings.push(format!("Too few {what} ({ind},{}) in {name}", modes.len()));
                            return false;
                        }
                        over.insert(id, modes[ind - 1].clone());
                        true
                    })
                    .cloned()
                    .collect()
            };
            let mut fopts = opts.clone();
            fopts.scalar_fields = keep(&opts.scalar_fields);
            fopts.vector_fields = keep(&opts.vector_fields);
            let variables = ModeVariables::new(store, &over);
            let path = self.out.dir.join(naming.file_name(&opts.output_file_name, ".vtu", file_index as i32, None, false));
            let mapped = MappedTags { mesh: ctx.mesh, model: &ctx.sim.model };
            elmer_vtu::write_vtu(&path, &mapped, &variables, &fopts)?;
            self.written.lock().unwrap().push(path);
        }
        for w in warnings {
            ctx.warn(CALLER, &w);
        }
        Ok(())
    }
}

impl Solver for VtuSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if self.out.mode_files {
            return self.write_mode_files(ctx);
        }
        // Without `Eigen Analysis`, Elmer adds `<name> EigenMode<i>` /
        // `ConstraintMode<i>` arrays to the one piece; not written here.
        if !ctx.sim.model.solvers[self.out.index].logical("ignore eigenmodes").unwrap_or(false) {
            let mm = max_modes(&ctx.sim.model, &ctx.sim.model.solvers[self.out.index]);
            let has_modes = ctx.sim.vars.iter().any(|(_, v)| {
                v.output && ((mm > 0 && !v.eigen_vectors.is_empty()) || !v.constraint_modes.is_empty())
            });
            if has_modes {
                return Err(SolverError::Unsupported(
                    "VTU output: eigen or constraint modes as extra fields of one piece (no Eigen Analysis in the output section)".into(),
                ));
            }
        }
        let opts = &self.out.opts;
        let index = self.indexer.advance(opts);
        let naming = VtuNaming {
            pes: self.parts,
            no_fileindex: opts.no_fileindex,
            parallel_base: opts.partition_numbering,
            ..VtuNaming::default()
        };
        std::fs::create_dir_all(&self.out.dir)?;
        let variables = ElmerVariables::new(&ctx.sim.vars);
        let part = if self.parts > 1 {
            // VtuOutputSolver.F90:622-627 then :633: the .pvtu first.
            if self.part == 0 {
                let pvtu = self.out.dir.join(naming.file_name(&opts.output_file_name, ".pvtu", index, None, false));
                let pieces: Vec<String> =
                    (1..=self.parts).map(|p| naming.file_name(&opts.output_file_name, ".vtu", index, Some(p), true)).collect();
                elmer_vtu::write_pvtu(&pvtu, &variables, opts, &pieces)?;
                self.written.lock().unwrap().push(pvtu);
            }
            Some(self.part + 1)
        } else {
            None
        };
        let name = naming.file_name(&opts.output_file_name, ".vtu", index, part, false);
        let path = self.out.dir.join(name);
        let mapped = MappedTags { mesh: ctx.mesh, model: &ctx.sim.model };
        let Some(reference) = self.out.reference else {
            elmer_vtu::write_vtu(&path, &mapped, &variables, opts)?;
            self.written.lock().unwrap().push(path);
            return Ok(());
        };
        // `AscBinInitNorm` at each visit, then the sums of what this visit
        // writes (no file when nothing is selected: all sums zero).
        let mut sums = [0.0; 6];
        if let Some(bytes) = elmer_vtu::render_vtu(&mapped, &variables, opts, None)? {
            std::fs::write(&path, &bytes)?;
            self.written.lock().unwrap().push(path);
            let params = &ctx.sim.model.solvers[self.out.index];
            let buffer = match params.integer("binary output buffer size") {
                Some(b) if b > 0 => b as usize,
                _ => ctx.mesh.n_nodes().max(1000),
            };
            sums = asc_bin_sums(&bytes, buffer)
                .ok_or_else(|| SolverError::fatal("ResultOutputSolver", "written VTU piece is not binary appended"))?;
        }
        // `Solver % Variable % Norm = Nrm`: the Post File hack's
        // `vtu_internal_dummy`, or the null variable of a section without one.
        let nrm = asc_bin_compare_norm(&reference, &sums);
        let id = ctx.solver_id();
        if let Some(slot) = ctx.sim.slots.iter_mut().find(|s| s.id == id) {
            match slot.variable {
                Some(v) => ctx.sim.vars.by_id_mut(v).conv.norm = nrm,
                None => slot.null_norm = nrm,
            }
        }
        ctx.info("ResultOutputSolver", &format!("Calculate Pseudonorm:{}", elmer_solver::fortran_fmt::es(nrm, 15, 6)), 4);
        Ok(())
    }
}

/// The first column of `Reference Values` (`RefResults(:,1)`); `None` with
/// fewer than the six entries `AscBinCompareNorm` reads.
fn reference_values(v: &elmer_core::Value) -> Option<[f64; 6]> {
    let col: Vec<f64> = match v {
        elmer_core::Value::Real { cols, data, .. } => data.iter().step_by((*cols).max(1)).copied().collect(),
        elmer_core::Value::Integer(v) => v.iter().map(|&i| f64::from(i)).collect(),
        _ => return None,
    };
    col.get(..6)?.try_into().ok()
}

/// The `AscBinOutputUtils` checksums (GeneralUtils.F90:2810-3051) of one
/// binary VTU piece written by `elmer-vtu`, whose bytes are Fortran's:
/// `[scount, icount, rcount, ssum, isum, rsum]`. `buffer` is `BufferSize`
/// (`Binary Output Buffer Size`, else `MAX(1000, Mesh % NumberOfNodes)`).
///
/// * Strings: one `AscBinStrWrite` per XML line, except a `DataArray` tag,
///   which is two calls, three when it has both a `Name` and
///   `NumberOfComponents`; then `</AppendedData>`, `</VTKFile>` and the
///   empty `' '`. `ssum` is their length, every XML byte but the `_`.
/// * Data arrays go through a buffer of `buffer` values flushed when full
///   and at the end of each array: every integer flush adds one to `icount`
///   and its sum to `isum`; every real flush adds its count to `rcount` and
///   *sets* `rsum` to its sum (`Rsum = SUM(...)`), so `rsum` is the sum of
///   the last real flush.
///
/// `None` when `bytes` is not a binary appended piece.
pub fn asc_bin_sums(bytes: &[u8], buffer: usize) -> Option<[f64; 6]> {
    const HEAD: &[u8] = b"<AppendedData encoding=\"raw\">\n";
    let head = bytes.windows(HEAD.len()).position(|w| w == HEAD)?;
    let xml = std::str::from_utf8(&bytes[..head + HEAD.len()]).ok()?;
    let data = bytes.get(head + HEAD.len() + 1..)?;
    let attr = |tag: &str, key: &str| -> Option<String> {
        let at = tag.find(&format!(" {key}=\""))? + key.len() + 3;
        Some(tag[at..].split('"').next()?.to_string())
    };
    let buffer = buffer.max(1);
    let (mut scount, mut icount, mut rcount) = (0usize, 0usize, 0usize);
    let (mut isum, mut rsum) = (0i64, 0.0f64);
    for line in xml.split_inclusive('\n') {
        if !line.trim_start().starts_with("<DataArray") {
            scount += 1;
            continue;
        }
        scount += if line.contains(" Name=\"") && line.contains(" NumberOfComponents=\"") { 3 } else { 2 };
        let offset: usize = attr(line, "offset")?.parse().ok()?;
        let len = u32::from_le_bytes(data.get(offset..offset + 4)?.try_into().ok()?) as usize;
        let raw = data.get(offset + 4..offset + 4 + len)?;
        match attr(line, "type")?.as_str() {
            "Int32" => {
                for chunk in raw.chunks(4 * buffer) {
                    icount += 1;
                    isum += chunk.chunks_exact(4).map(|b| i64::from(i32::from_le_bytes(b.try_into().unwrap()))).sum::<i64>();
                }
            }
            "Float64" => {
                for chunk in raw.chunks(8 * buffer) {
                    rcount += chunk.len() / 8;
                    rsum = chunk.chunks_exact(8).fold(0.0, |s, b| s + f64::from_le_bytes(b.try_into().unwrap()));
                }
            }
            "Float32" => {
                for chunk in raw.chunks(4 * buffer) {
                    rcount += chunk.len() / 4;
                    rsum = f64::from(chunk.chunks_exact(4).fold(0.0f32, |s, b| s + f32::from_le_bytes(b.try_into().unwrap())));
                }
            }
            _ => return None,
        }
    }
    const TAIL: &str = "\n</AppendedData>\n</VTKFile>\n";
    scount += 3;
    let ssum = xml.len() + TAIL.len();
    Some([scount as f64, icount as f64, rcount as f64, ssum as f64, isum as f64, rsum])
}

/// `AscBinCompareNorm` (GeneralUtils.F90:3053): the mean over the six
/// checksums of `max(c, 1/c)` with `c` this over the reference (`1 + |this|`
/// where the reference is zero); one when they agree.
pub fn asc_bin_compare_norm(reference: &[f64; 6], this: &[f64; 6]) -> f64 {
    let mut norm = 0.0;
    for (&r, &t) in reference.iter().zip(this) {
        let c = if r.abs() > f64::EPSILON {
            let c = t / r;
            if t.abs() > f64::EPSILON {
                c.max(1.0 / c)
            } else {
                c
            }
        } else {
            1.0 + t.abs()
        };
        norm += c;
    }
    norm / 6.0
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

/// [`ElmerVariables`] for one mode file: only the variables in `over`, with
/// the mode's values in place of their own (their components too).
struct ModeVariables<'a> {
    store: &'a VariableStore,
    over: &'a std::collections::HashMap<VarId, Vec<f64>>,
    components: Vec<(String, VarId, Vec<f64>)>,
}

impl<'a> ModeVariables<'a> {
    fn new(store: &'a VariableStore, over: &'a std::collections::HashMap<VarId, Vec<f64>>) -> Self {
        let components = store
            .names_in_order()
            .into_iter()
            .filter_map(|(name, h)| match h {
                VarHandle::Component { parent, component } => {
                    let vals = over.get(&parent)?;
                    let dofs = store.by_id(parent).dofs;
                    Some((name, parent, vals.iter().skip(component).step_by(dofs).copied().collect()))
                }
                VarHandle::Primary(_) => None,
            })
            .collect();
        ModeVariables { store, over, components }
    }
}

impl FieldSource for ModeVariables<'_> {
    fn fields(&self) -> Vec<FieldRef<'_>> {
        let mut out = Vec::new();
        for (id, v) in self.store.iter() {
            let Some(vals) = self.over.get(&id) else { continue };
            if !v.output {
                continue;
            }
            out.push(FieldRef { values: vals, ..v.field_ref() });
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

    /// A mode file shows only the variables that have the mode, with the
    /// mode's values under the variable's and its components' own names.
    #[test]
    fn mode_files_show_the_mode_values_under_the_field_names() {
        let mut store = VariableStore::new();
        let pot = store
            .add_with_components(nodal("pot", 2, vec![1.0, 10.0, 2.0, 20.0]), &["temperature".into(), "potential".into()])
            .unwrap();
        store.add(nodal("other", 1, vec![7.0, 8.0])).unwrap();
        assert_eq!(field_variable(&store, "potential"), Some(pot));
        assert_eq!(field_variable(&store, "pot"), Some(pot));
        let over = std::collections::HashMap::from([(pot, vec![5.0, 50.0, 6.0, 60.0])]);
        let variables = ModeVariables::new(&store, &over);
        let fields = variables.fields();
        let names: Vec<&str> = fields.iter().map(|f| f.name).collect();
        assert_eq!(names, ["pot", "temperature", "potential"]);
        assert_eq!(fields[0].values, [5.0, 50.0, 6.0, 60.0]);
        assert_eq!(fields[1].values, [5.0, 6.0]);
        assert_eq!(fields[2].values, [50.0, 60.0]);
    }

    /// Buffers of 2: reals flush [1, 2] then [4] (`rsum` is the last flush,
    /// 4), integers [1, 2] then [3] (two flushes, sum 6).
    #[test]
    fn asc_bin_sums_follow_the_buffered_writes() {
        let mut xml = String::new();
        xml.push_str("<?xml version=\"1.0\"?>\n");
        xml.push_str("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        xml.push_str("        <DataArray type=\"Float64\" Name=\"t\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\"/>\n");
        xml.push_str("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"28\"/>\n");
        xml.push_str("<AppendedData encoding=\"raw\">\n");
        let mut bytes = xml.clone().into_bytes();
        bytes.push(b'_');
        bytes.extend_from_slice(&24u32.to_le_bytes());
        for x in [1.0f64, 2.0, 4.0] {
            bytes.extend_from_slice(&x.to_le_bytes());
        }
        bytes.extend_from_slice(&12u32.to_le_bytes());
        for i in [1i32, 2, 3] {
            bytes.extend_from_slice(&i.to_le_bytes());
        }
        bytes.extend_from_slice(b"\n</AppendedData>\n</VTKFile>\n");
        let sums = asc_bin_sums(&bytes, 2).unwrap();
        assert_eq!(sums, [11.0, 2.0, 3.0, (xml.len() + 28) as f64, 6.0, 4.0]);
        assert_eq!(asc_bin_compare_norm(&sums, &sums), 1.0);
        let mut off = sums;
        off[1] = 4.0;
        assert_eq!(asc_bin_compare_norm(&sums, &off), (5.0 + 2.0) / 6.0);
        assert_eq!(asc_bin_sums(b"<VTKFile/>", 2), None);
    }

    #[test]
    fn displacement_with_mesh_update_complement_is_refused() {
        let mut store = VariableStore::new();
        store.add(nodal("displacement", 2, vec![0.0; 4])).unwrap();
        store.add(nodal("mesh update", 2, vec![0.0; 4])).unwrap();
        assert!(matches!(fields_for_saving(&store, true, 2), Err(SolverError::Unsupported(_))));
    }
}
