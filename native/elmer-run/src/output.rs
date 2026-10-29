//! The `Post File = "*.vtu"` auto-solver (ElmerSolver.F90 `AddVtuOutputSolverHack`)
//! written with `elmer-vtu` once the run finishes.
//!
//! Only the steady case is composed: the output solver runs `after saving`,
//! which for a steady run is once, after the last steady iteration. Transient
//! runs with a `Post File` are refused until the output solver is a scheduled
//! `Solver` slot.

use std::path::{Path, PathBuf};

use elmer_core::{FieldKind, FieldSource, MeshView, Params, Value};
use elmer_solver::{variable::Snapshot, Model, Simulation, SolverError};
use elmer_vtu::{post_file_solver, FileIndexer, GeometryIds, VtuNaming, VtuOptions};

/// A configured VTU output.
#[derive(Clone, Debug)]
pub struct VtuOutput {
    opts: VtuOptions,
    dir: PathBuf,
}

/// Turns the Simulation `Post File` into a VTU output, removing the keyword
/// as Fortran's `ListRemove` does. `mesh_dir` is the default output
/// directory (`SolverOutputDirectory(..., UseMeshDir = .TRUE.)`), `case_dir`
/// the base of a relative `Output Directory`. Refuses what the composition
/// cannot write correctly.
pub fn take_post_file(model: &mut Model, case_dir: &Path, mesh_dir: &Path) -> Result<Option<VtuOutput>, SolverError> {
    let Some(post) = model.simulation.string("post file").map(str::to_string) else {
        return Ok(None);
    };
    let Some(params) = post_file_solver(&model.simulation, &model.solvers) else {
        return Err(SolverError::Unsupported(format!("Post File = {post} (only .vtu is composed)")));
    };
    model.simulation.0.remove("post file");
    let transient = model
        .simulation
        .string("simulation type")
        .is_some_and(|s| s.eq_ignore_ascii_case("transient"));
    if transient {
        return Err(SolverError::Unsupported("Post File in a transient simulation".into()));
    }

    let mut opts = VtuOptions::from_params(&params);
    opts.apply_simulation(&model.simulation);
    if !opts.unsupported.is_empty() {
        return Err(SolverError::Unsupported(format!("VTU output: {}", opts.unsupported.join(", "))));
    }
    if opts.output_file_name.contains('/') {
        return Err(SolverError::Unsupported("a directory in the Post File name".into()));
    }
    if opts.save_geometry_ids {
        opts = opts.with_geometry(geometry(model)?);
    }
    let dir = match params
        .string("output directory")
        .or_else(|| model.simulation.string("output directory"))
    {
        Some(d) => case_dir.join(d),
        None => mesh_dir.to_path_buf(),
    };
    Ok(Some(VtuOutput { opts, dir }))
}

/// Body and BC facts for `GeometryIds`. A BC's tag is the mesh boundary it
/// targets (Fortran maps `Constraint` to the BC section at load time), so a
/// BC that targets several boundaries cannot be expressed by `elmer-vtu`'s
/// one-tag-per-BC `GeometryIds` and is refused.
fn geometry(model: &Model) -> Result<GeometryIds, SolverError> {
    let mut bc_tags = Vec::new();
    for (i, bc) in model.boundary_conditions.iter().enumerate() {
        match bc.get("target boundaries") {
            Some(Value::Integer(t)) if t.len() == 1 => bc_tags.push(t[0]),
            None => bc_tags.push(0),
            Some(_) => {
                return Err(SolverError::Unsupported(format!(
                    "Save Geometry Ids with several Target Boundaries in Boundary Condition {}",
                    i + 1
                )))
            }
        }
    }
    let ids = |secs: &[elmer_core::ParamMap]| secs.iter().map(|s| s.integer("geometry id")).collect();
    let names = |secs: &[elmer_core::ParamMap]| secs.iter().map(|s| s.string("name").map(str::to_string)).collect();
    Ok(GeometryIds {
        n_bodies: model.bodies.len(),
        bc_tags,
        body_geometry_id: ids(&model.bodies),
        bc_geometry_id: ids(&model.boundary_conditions),
        body_names: names(&model.bodies),
        bc_names: names(&model.boundary_conditions),
    })
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

    /// Writes the one steady-state piece and returns its path.
    pub fn write(&self, sim: &Simulation, mesh: &dyn MeshView) -> Result<PathBuf, SolverError> {
        let index = FileIndexer::new().advance(&self.opts);
        let naming = VtuNaming { no_fileindex: self.opts.no_fileindex, ..VtuNaming::default() };
        let name = naming.file_name(&self.opts.output_file_name, ".vtu", index, None, false);
        std::fs::create_dir_all(&self.dir)?;
        let path = self.dir.join(name);
        elmer_vtu::write_vtu(&path, mesh, &Snapshot { store: &sim.vars }, &self.opts)?;
        Ok(path)
    }
}
