//! Declared component names (`Variable = Name[Comp A:1 Comp B:1]`): names
//! with spaces parse, and Dirichlet conditions are read by those names as
//! `ComponentName` gives them to `DefaultDirichletBCs`. The Helmholtz case
//! of `tests/data/complex/direct.sif` with `Pressure[Pressure Re:1 Pressure
//! Im:1]` and its boundary values under those names solves exactly as with
//! `-dofs 2 Pressure` and `Pressure 1` / `Pressure 2`. HelmholtzSolve is the
//! test-local stand-in.

mod support;

use std::path::Path;

use elmer_core::{SerialComm, Value};
use elmer_solver::{Model, Simulation, SolverRegistry};
use support::poisson::helmholtz;
use support::{sif_model, Capture};

fn run(model: Model) -> Simulation {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join("1d"), 1);
    let mut reg = SolverRegistry::new();
    reg.register("HelmholtzSolve", "HelmholtzSolver", helmholtz);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    sim
}

#[test]
fn dirichlet_values_by_declared_component_names() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/complex/direct.sif");
    let want = run(sif_model(&path));
    let want = &want.vars.get("pressure").unwrap().values;

    let mut model = sif_model(&path);
    model.solvers[0].insert("Variable", Value::Str("pressure[pressure re:1 pressure im:1]".into()));
    for bc in &mut model.boundary_conditions {
        for (from, to) in [("pressure 1", "pressure re"), ("pressure 2", "pressure im")] {
            if let Some(v) = bc.0.remove(from) {
                bc.insert(to, v);
            }
        }
    }
    let sim = run(model);
    let id = sim.vars.id("pressure").unwrap();
    assert_eq!(sim.vars.component_names(id), vec!["pressure re", "pressure im"]);
    assert_eq!(&sim.vars.get("pressure").unwrap().values, want);
    // The old names are not components any more.
    assert!(sim.vars.resolve("pressure 1").is_none());
}
