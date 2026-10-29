//! What `FlowSolver` refuses instead of approximating, and FlowSolve's
//! keyword defaults that change the run.

mod support;

use elmer_core::{ParamMap, Value};
use elmer_solver::SolverError;
use support::cases::{data, flow_solver, ints, logical, model, nodal, real, run, sec, st, steady};
use support::mesh::TestMesh;

fn stokes(extra: &[(&str, Value)], material: ParamMap, bc: Vec<(&str, Value)>) -> Result<elmer_solver::Simulation, SolverError> {
    let mesh = TestMesh::read(&data("stokes_channel"), 2);
    let mut v = vec![("Flow Model", st("Stokes")), ("Stabilize", logical(true))];
    v.extend(extra.iter().cloned());
    let mut walls = vec![("Target Boundaries", ints(&[1, 3])), ("Velocity 1", real(0.0)), ("Velocity 2", real(0.0))];
    walls.extend(bc);
    run(model(steady(), flow_solver(2, &v), material, None, vec![sec(&walls)]), &mesh)
}

fn water() -> ParamMap {
    sec(&[("Density", real(1.0)), ("Viscosity", real(1.0))])
}

fn assert_unsupported(r: Result<elmer_solver::Simulation, SolverError>, what: &str) {
    match r {
        Err(SolverError::Unsupported(m)) => assert!(m.to_lowercase().contains(what), "{m}"),
        Err(e) => panic!("expected Unsupported({what}), got {e}"),
        Ok(_) => panic!("expected Unsupported({what}), got a run"),
    }
}

#[test]
fn refuses_unported_options() {
    let it = ("Nonlinear System Max Iterations", ints(&[1]));
    assert_unsupported(
        stokes(&[it.clone(), ("Stabilization Method", st("Bubbles"))], water(), vec![]),
        "bubbles",
    );
    let mut gas = water();
    gas.insert("Compressibility Model", st("Perfect Gas"));
    assert_unsupported(stokes(&[it.clone()], gas, vec![]), "compressibility");
    assert_unsupported(
        stokes(&[it.clone()], water(), vec![("Normal-Tangential Velocity", logical(true))]),
        "normal-tangential",
    );
    assert_unsupported(
        stokes(&[it.clone()], water(), vec![("Surface Tension Coefficient", real(1.0))]),
        "surface tension",
    );
    let mut nn = water();
    nn.insert("Viscosity Model", st("power law"));
    assert_unsupported(stokes(&[it], nn, vec![]), "viscosity model");
}

/// `Nonlinear System Max Iterations` is read without a default, so FlowSolve
/// without it runs no iteration at all: the field keeps its start value
/// (1e-6 for `Equation = Navier-Stokes`, boundary values where Dirichlet
/// conditions initialize it).
#[test]
fn absent_max_iterations_runs_nothing() {
    let sim = stokes(&[], water(), vec![]).unwrap();
    let mesh = TestMesh::read(&data("stokes_channel"), 2);
    for (k, v) in nodal(&sim, mesh.coords.len()).iter().enumerate() {
        let on_wall = mesh.coords[k][1] == 0.0 || mesh.coords[k][1] == 1.0;
        let u = if on_wall { 0.0 } else { 1.0e-6 };
        assert_eq!(v, &vec![u, u, 1.0e-6], "node {}", k + 1);
    }
}
