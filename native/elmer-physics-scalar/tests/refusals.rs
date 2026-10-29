//! What the solvers do not implement is refused when the solver is set up,
//! instead of being dropped.

mod support;

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_solver::{Model, SolverError};
use support::cases::*;
use support::meshio::read_elmer_mesh;
use support::simmesh::SimMesh;

fn model(
    procedure: &str,
    solver: &[(&str, Value)],
    material: &[(&str, Value)],
    bc: &[(&str, Value)],
    dofs: i32,
) -> Model {
    let mut s = vec![
        ("Equation", st("Test")),
        ("Variable", st("Potential")),
        ("Variable DOFs", ints(&[dofs])),
        ("Procedure", st(procedure)),
        ("Linear System Solver", st("Direct")),
    ];
    s.extend(solver.iter().cloned());
    let mut b = vec![("Target Boundaries", ints(&[1]))];
    b.extend(bc.iter().cloned());
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian 2D")),
            ("Simulation Type", st("Steady State")),
            ("Angular Frequency", real(100.0)),
        ]),
        solvers: vec![sec(&s)],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[("Equation", ints(&[1])), ("Material", ints(&[1]))])],
        materials: vec![sec(material)],
        boundary_conditions: vec![sec(&b)],
        ..Model::default()
    }
}

fn setup(m: Model) -> Result<(), SolverError> {
    let mesh = SimMesh::new(read_elmer_mesh(&data("helmholtz_fem"), 2));
    let (mut sim, reg) = simulation(m);
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm)
}

fn refused(m: Model, what: &str) {
    match setup(m) {
        Err(SolverError::Unsupported(msg)) => assert!(msg.contains(what), "{msg}"),
        other => panic!("expected Unsupported({what}), got {other:?}"),
    }
}

const ELSTAT: &str = "\"StatElecSolve\" \"StatElecSolver\"";
const CURRENT: &str = "\"StatCurrentSolve\" \"StatCurrentSolver\"";
const HELMHOLTZ: &str = "\"HelmholtzSolve\" \"HelmholtzSolver\"";

#[test]
fn statelec_refusals() {
    let mat = [("Relative Permittivity", real(1.0))];
    let t = Value::Logical(true);
    refused(
        model(
            ELSTAT,
            &[],
            &[("Relative Permittivity Im", real(0.1))],
            &[],
            1,
        ),
        "Relative Permittivity Im",
    );
    refused(
        model(ELSTAT, &[], &mat, &[("Electric Infinity BC", t.clone())], 1),
        "Electric Infinity BC",
    );
    refused(
        model(ELSTAT, &[("Legacy Assembly", t.clone())], &mat, &[], 1),
        "legacy assembly",
    );
    refused(
        model(
            ELSTAT,
            &[("Calculate Capacitance Matrix", t.clone())],
            &mat,
            &[],
            1,
        ),
        "capacitance",
    );
    refused(
        model(
            ELSTAT,
            &[("Calculate Elemental Fields", t.clone())],
            &mat,
            &[],
            1,
        ),
        "elemental",
    );
    refused(
        model(ELSTAT, &[], &mat, &[("Capacitance Body", ints(&[1]))], 1),
        "capacitance body",
    );
    refused(model(ELSTAT, &[], &mat, &[], 2), "Variable DOFs");
}

#[test]
fn statcurrent_refusals() {
    let mat = [("Electric Conductivity", real(1.0))];
    let t = Value::Logical(true);
    refused(
        model(CURRENT, &[("Harmonic Mode", t.clone())], &mat, &[], 1),
        "harmonic mode",
    );
    refused(
        model(CURRENT, &[], &mat, &[("Farfield Potential", real(0.0))], 1),
        "Farfield Potential",
    );
    refused(
        model(
            CURRENT,
            &[("Calculate Conductivity Matrix", t.clone())],
            &mat,
            &[],
            1,
        ),
        "conductivity matrix",
    );
}

#[test]
fn helmholtz_refusals() {
    let mat = [("Density", real(1.0)), ("Sound Speed", real(340.0))];
    let t = Value::Logical(true);
    refused(
        model(HELMHOLTZ, &[("Bubbles", t.clone())], &mat, &[], 2),
        "Bubbles",
    );
    refused(
        model(
            HELMHOLTZ,
            &[],
            &[("Density", real(1.0)), ("Convection Velocity 1", real(1.0))],
            &[],
            2,
        ),
        "Convection Velocity",
    );
    refused(
        model(HELMHOLTZ, &[], &mat, &[("Flow Interface", t.clone())], 2),
        "flow interface",
    );
    let one_dof = setup(model(HELMHOLTZ, &[], &mat, &[], 1));
    assert!(
        matches!(one_dof, Err(SolverError::Fatal { .. })),
        "{one_dof:?}"
    );
}

#[test]
fn registered_factories_build_once_a_provider_is_installed() {
    install();
    let mut reg = elmer_solver::SolverRegistry::new();
    elmer_physics_scalar::register(&mut reg);
    for (lib, entry) in [
        ("StatElecSolve", "StatElecSolver"),
        ("StatCurrentSolve", "StatCurrentSolver"),
        ("HelmholtzSolve", "HelmholtzSolver"),
    ] {
        assert!(reg.get(lib, entry).unwrap()(&ParamMap::new()).is_ok());
    }
}
