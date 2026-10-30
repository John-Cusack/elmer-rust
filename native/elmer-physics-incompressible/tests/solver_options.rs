//! What the solver refuses instead of approximating, and the keyword
//! defaults of `IncompressibleNSSolver_Init0`/`_init`.

mod support;

use elmer_core::{ParamMap, Params, Value};
use elmer_physics_incompressible::solver_defaults;
use elmer_solver::SolverError;
use support::cases::{data, ints, logical, model, ns_solver, real, run, sec, st, steady};
use support::mesh::TestMesh;

fn step_mesh() -> TestMesh {
    TestMesh::read(&data("stokes_tri"), 2)
}

fn walls() -> Vec<ParamMap> {
    vec![sec(&[
        ("Target Boundaries", ints(&[1, 2, 3])),
        ("Velocity 1", real(0.0)),
        ("Velocity 2", real(0.0)),
    ])]
}

fn refused(
    solver_extra: &[(&str, Value)],
    material: &[(&str, Value)],
    bcs: Vec<ParamMap>,
    needle: &str,
) {
    let mesh = step_mesh();
    let mut m = vec![("Density", real(1.0)), ("Viscosity", real(1.0))];
    m.extend(material.iter().cloned());
    let solver = ns_solver(solver_extra);
    match run(model(steady("Cartesian 2D"), solver, sec(&m), None, bcs), &mesh) {
        Err(SolverError::Unsupported(msg)) => assert!(msg.contains(needle), "{msg}"),
        Err(e) => panic!("expected a refusal naming {needle:?}, got {e}"),
        Ok(_) => panic!("expected a refusal naming {needle:?}"),
    }
}

#[test]
fn default_mini_element_is_refused() {
    refused(&[("Nonlinear System Max Iterations", ints(&[1]))], &[], walls(), "MINI");
}

#[test]
fn bubble_element_with_stabilization_is_refused() {
    refused(
        &[("Pressure Stabilization", logical(true)), ("Element", st("p:1 -tri b:1"))],
        &[],
        walls(),
        "Element",
    );
}

#[test]
fn viscosity_model_is_refused() {
    refused(
        &[("Pressure Stabilization", logical(true))],
        &[("Viscosity Model", st("glen"))],
        walls(),
        "Viscosity Model",
    );
}

#[test]
fn compressibility_is_refused() {
    refused(
        &[("Pressure Stabilization", logical(true))],
        &[("Compressibility Model", st("Thermal"))],
        walls(),
        "Compressibility",
    );
}

#[test]
fn old_pressure_traction_keyword_is_refused() {
    let mut bcs = walls();
    bcs.push(sec(&[("Target Boundaries", ints(&[2])), ("Pressure 1", real(1.0))]));
    refused(&[("Pressure Stabilization", logical(true))], &[], bcs, "Surface Traction 1");
}

#[test]
fn friction_boundary_is_refused() {
    let mut bcs = walls();
    bcs.push(sec(&[("Target Boundaries", ints(&[3])), ("Friction Coefficient", real(1.0))]));
    refused(&[("Pressure Stabilization", logical(true))], &[], bcs, "Friction");
}

#[test]
fn unknown_element_size_is_fatal() {
    let mesh = step_mesh();
    let solver = ns_solver(&[
        ("Pressure Stabilization", logical(true)),
        ("Pressure Stabilization Element Size", st("diagonal")),
    ]);
    let m = sec(&[("Density", real(1.0)), ("Viscosity", real(1.0))]);
    let e = run(model(steady("Cartesian 2D"), solver, m, None, walls()), &mesh).err().expect("must fail");
    assert!(e.to_string().contains("Unknown \"Pressure Stabilization Element Size\""), "{e}");
}

#[test]
fn defaults_follow_init0_and_init() {
    let mut p = sec(&[("Pressure Stabilization", logical(true)), ("Flow Model", st("Stokes"))]);
    solver_defaults(&mut p, 2, false);
    assert_eq!(p.string("element"), Some("n:1"));
    assert_eq!(p.string("variable"), Some("Flow Solution[Velocity:2 Pressure:1]"));
    assert_eq!(p.integer("nonlinear system norm dofs"), Some(2));
    assert_eq!(p.logical("relative pressure relaxation"), Some(true));
    assert_eq!(p.logical("stokes flow"), Some(true));
    assert_eq!(p.logical("bubbles in global system"), Some(false));

    // ListAddNew: what the section states wins.
    let mut q = sec(&[("Stabilize", logical(true)), ("Element", st("n:1")), ("GradP Discretization", logical(true))]);
    solver_defaults(&mut q, 3, true);
    assert_eq!(q.string("variable"), Some("Flow Solution[Velocity:3 Pressure:1]"));
    assert_eq!(q.logical("gradp discretization"), Some(true));

    // No stabilization: _Init0's MINI element, and _Init's bubble rules.
    let mut r = ParamMap::new();
    solver_defaults(&mut r, 2, false);
    assert_eq!(r.string("element"), Some("p:1 -tri b:1 -tetra b:1 -quad b:3 -brick b:4 -prism b:4 -pyramid b:4"));
    assert_eq!(r.string("element integration points"), Some("-tetra 24 -prism 85"));
}
