//! Models mirroring the `tests/data/<case>/*.sif` cases, built in code.
#![allow(dead_code)]

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_physics_elastic::{install_fe, register, FeServices};
use elmer_solver::{Model, Simulation, SolverRegistry};

use super::sim::{dense_factory, read_mesh, sec, DenseAssembly, Mesh};
use super::{Gauss, P1Basis};

pub fn try_run(model: Model, mesh: &Mesh) -> Result<Simulation, elmer_solver::SolverError> {
    let _ = install_fe(FeServices { basis: Box::new(P1Basis), quad: Box::new(Gauss) });
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(model);
    sim.backends.assembly = Some(Box::new(DenseAssembly));
    sim.backends.linear.register_direct("umfpack", dense_factory);
    sim.setup(mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

pub fn run(model: Model, mesh: &Mesh) -> Simulation {
    try_run(model, mesh).unwrap()
}

pub fn solver_params(extra: &[(&str, Value)]) -> ParamMap {
    let mut s = vec![
        ("Equation", Value::Str("Stress Solver".into())),
        ("Procedure", Value::Str("\"StressSolve\" \"StressSolver\"".into())),
        ("Variable", Value::Str("Displacement".into())),
        ("Linear System Solver", Value::Str("Direct".into())),
        ("Linear System Direct Method", Value::Str("umfpack".into())),
        ("Calculate Stresses", Value::Logical(true)),
        ("Calculate Strains", Value::Logical(true)),
        ("Calculate Principal", Value::Logical(true)),
    ];
    s.extend(extra.iter().map(|(k, v)| (*k, v.clone())));
    sec(&s)
}

pub fn simulation(coord: &str) -> ParamMap {
    sec(&[
        ("Coordinate System", Value::Str(coord.into())),
        ("Simulation Type", Value::Str("Steady State".into())),
        ("Steady State Max Iterations", Value::Integer(vec![1])),
    ])
}

pub fn bc(tag: i32, entries: &[(&str, Value)]) -> ParamMap {
    let mut e = vec![("Target Boundaries", Value::Integer(vec![tag]))];
    e.extend(entries.iter().map(|(k, v)| (*k, v.clone())));
    sec(&e)
}

pub fn body(force: bool) -> ParamMap {
    let mut e = vec![("Equation", Value::Integer(vec![1])), ("Material", Value::Integer(vec![1]))];
    if force {
        e.push(("Body Force", Value::Integer(vec![1])));
    }
    sec(&e)
}

pub fn a2d() -> (Model, Mesh) {
    let m = Model {
        simulation: simulation("Cartesian 2D"),
        solvers: vec![solver_params(&[("Variable DOFs", Value::Integer(vec![2]))])],
        equations: vec![sec(&[("Active Solvers", Value::Integer(vec![1])), ("Plane Stress", Value::Logical(false))])],
        bodies: vec![body(true)],
        body_forces: vec![sec(&[
            ("Stress Bodyforce 1", Value::real(1.0e3)),
            ("Stress Bodyforce 2", Value::real(-2.0e3)),
        ])],
        materials: vec![sec(&[
            ("Density", Value::real(1.0)),
            ("Youngs Modulus", Value::real(70e9)),
            ("Poisson Ratio", Value::real(0.3)),
        ])],
        boundary_conditions: vec![
            bc(1, &[("Displacement 1", Value::real(0.0)), ("Displacement 2", Value::real(0.0))]),
            bc(2, &[("Surface Traction 2", Value::real(-5.0e7)), ("Normal Surface Traction", Value::real(1.0e6))]),
        ],
        ..Model::default()
    };
    (m, read_mesh("a2d", 2))
}

pub fn e2d() -> (Model, Mesh) {
    let m = Model {
        simulation: simulation("Axi Symmetric"),
        solvers: vec![solver_params(&[("Variable DOFs", Value::Integer(vec![2]))])],
        equations: vec![sec(&[("Active Solvers", Value::Integer(vec![1])), ("Plane Stress", Value::Logical(false))])],
        bodies: vec![body(true)],
        body_forces: vec![sec(&[
            ("Stress Bodyforce 1", Value::real(1.0e3)),
            ("Stress Bodyforce 2", Value::real(-2.0e3)),
            ("Stress Pressure", Value::real(5.0e2)),
        ])],
        materials: vec![sec(&[
            ("Density", Value::real(1.0)),
            ("Youngs Modulus", Value::real(70e9)),
            ("Poisson Ratio", Value::real(0.3)),
        ])],
        boundary_conditions: vec![
            bc(1, &[("Displacement 1", Value::real(0.0)), ("Displacement 2", Value::real(0.0))]),
            bc(2, &[("Surface Traction 2", Value::real(-5.0e7)), ("Normal Surface Traction", Value::real(1.0e6))]),
            bc(3, &[("Spring", Value::real(2.0e9)), ("Normal Surface Traction", Value::real(3.0e5))]),
            bc(
                4,
                &[
                    ("Spring 1", Value::real(1.0e9)),
                    ("Spring 2", Value::real(4.0e9)),
                    ("Spring 12", Value::real(5.0e8)),
                    ("Surface Traction 1", Value::real(2.0e5)),
                ],
            ),
        ],
        ..Model::default()
    };
    (m, read_mesh("e2d", 2))
}

pub fn c3d() -> (Model, Mesh) {
    #[rustfmt::skip]
    let c = vec![
        1346.154, 576.923, 576.923, 0.0, 0.0, 0.0,
        576.923, 1346.154, 576.923, 0.0, 0.0, 0.0,
        576.923, 576.923, 1346.154, 200.0, 0.0, 0.0,
        0.0, 0.0, 200.0, 384.615, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 384.615, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 384.615,
    ];
    let m = Model {
        simulation: simulation("Cartesian"),
        solvers: vec![solver_params(&[("Variable DOFs", Value::Integer(vec![3]))])],
        equations: vec![sec(&[("Active Solvers", Value::Integer(vec![1]))])],
        bodies: vec![body(false)],
        materials: vec![sec(&[
            ("Density", Value::real(1.0)),
            ("Youngs Modulus", Value::Real { rows: 6, cols: 6, data: c }),
        ])],
        boundary_conditions: vec![
            bc(5, &[("Displacement 1", Value::real(0.0)), ("Displacement 2", Value::real(0.0)), ("Displacement 3", Value::real(0.0))]),
            bc(6, &[("Displacement 1", Value::real(1.0e-3)), ("Displacement 2", Value::real(2.0e-3)), ("Displacement 3", Value::real(0.5e-3))]),
        ],
        ..Model::default()
    };
    (m, read_mesh("c3d", 3))
}

