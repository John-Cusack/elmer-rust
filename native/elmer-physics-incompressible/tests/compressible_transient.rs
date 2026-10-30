//! `tests/data/thermal_transient_quad/case.sif` (and its perfect gas twin)
//! against the pristine
//! ElmerSolver: ModelPDE (`elmer-physics-advdiff`, a diffusion equation for
//! the variable `Temperature`) heats ThermalCompress_vec's cavity from the
//! walls each BDF1 step and IncompressibleNS follows with `Compressibility
//! Model = Thermal` on the MINI element, so the density changes in time:
//! the d(rho)/dt load of the continuity equation, with the bubbles
//! condensed with their history.

mod support;

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_pdofs::PDofsBackend;
use elmer_physics_incompressible::{register, solver_defaults};
use elmer_solver::{Model, Simulation, SolverRegistry};
use support::cases::{data, ints, ns_solver, real, sec, st, table};
use support::mesh::TestMesh;

fn lines(path: &std::path::Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .filter(|l: &Vec<String>| !l.is_empty())
        .collect()
}

fn bc(tags: &[i32], entries: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![("Target Boundaries", ints(tags))];
    v.extend(entries.iter().cloned());
    sec(&v)
}

#[test]
fn transient_thermal_compressibility_matches_elmer() {
    let model = [("Compressibility Model", st("Thermal")), ("Heat Expansion Coefficient", real(0.25))];
    transient_cavity("thermal_transient_quad", &model, 0.005155763633228127);
}

/// The same with `Compressibility Model = Perfect Gas`: the pressure mass
/// `-(1/(p + p0)) (dp/dt, q)` as well as the temperature's d(rho)/dt load.
#[test]
fn transient_perfect_gas_compressibility_matches_elmer() {
    let model = [
        ("Heat Capacity", real(1000.0)),
        ("Compressibility Model", st("Perfect Gas")),
        ("Specific Heat Ratio", real(5.0 / 3.0)),
        ("Reference Pressure", real(1.0e5)),
    ];
    transient_cavity("perfect_gas_transient_quad", &model, 6.323421412565369e-05);
}

/// The transient cavity against the pristine `case`, the material's
/// compressibility keywords `model`.
fn transient_cavity(case: &str, model: &[(&str, Value)], ref_norm: f64) {
    let dir = data(case);
    let mesh = TestMesh::read(&dir, 2);
    let heat = sec(&[
        ("Equation", st("Temperature")),
        ("Procedure", st("\"ModelPDE\" \"AdvDiffSolver\"")),
        ("Variable", st("Temperature")),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Linear System Solver", st("Direct")),
        ("Linear System Direct Method", st("UMFPack")),
    ]);
    let mut flow = ns_solver(&[
        ("Equation", st("Navier-Stokes")),
        ("Nonlinear System Max Iterations", ints(&[10])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
        ("Nonlinear System Newton After Iterations", ints(&[20])),
        ("Nonlinear System Newton After Tolerance", real(1.0e-3)),
    ]);
    solver_defaults(&mut flow, 2, true);
    let wall = |t: i32, temp: Option<Value>| {
        let mut v = vec![("Velocity 1", real(0.0)), ("Velocity 2", real(0.0))];
        v.extend(temp.map(|x| ("Temperature", x)));
        bc(&[t], &v)
    };
    let mut first = wall(1, None);
    first.insert("Pressure", real(0.0));
    let bcs = vec![
        first,
        wall(2, Some(real(293.5))),
        wall(3, Some(table(1, &[(0.0, 293.5), (1.0, 292.5)]))),
        wall(4, Some(real(292.5))),
    ];
    let mut material = vec![
        ("Density", real(1.0)),
        ("Viscosity", real(1.0)),
        ("Diffusion Coefficient", real(1.0)),
        ("Time Derivative Coefficient", real(1.0)),
        ("Reference Temperature", real(293.0)),
    ];
    material.extend(model.iter().cloned());
    let model = Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian 2D")),
            ("Simulation Type", st("Transient")),
            ("Timestepping Method", st("BDF")),
            ("BDF Order", ints(&[1])),
            ("Timestep Sizes", real(0.1)),
            ("Timestep Intervals", ints(&[3])),
            ("Steady State Max Iterations", ints(&[1])),
        ]),
        solvers: vec![heat, flow],
        equations: vec![sec(&[("Active Solvers", ints(&[1, 2]))])],
        bodies: vec![sec(&[
            ("Target Bodies", ints(&[1])),
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
            ("Body Force", ints(&[1])),
            ("Initial Condition", ints(&[1])),
        ])],
        materials: vec![sec(&material)],
        body_forces: vec![sec(&[("Flow BodyForce 2", real(-10.0))])],
        initial_conditions: vec![sec(&[("Temperature", real(293.0))])],
        boundary_conditions: bcs,
        ..Model::default()
    };
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    elmer_physics_advdiff::register_model_pde(&mut reg);
    let mut sim = Simulation::new(model);
    sim.backends.dofs = Some(Box::new(PDofsBackend::new()));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();

    let u = sim.vars.get("flow solution").unwrap();
    let perm = u.perm.as_deref().unwrap();
    let rows = lines(&dir.join("perm.dat"));
    assert_eq!(rows.len(), 441);
    let want: Vec<Vec<f64>> = lines(&dir.join("solution.dat"))
        .iter()
        .map(|l| l[1..].iter().map(|x| x.parse().unwrap()).collect())
        .collect();
    let mut scale = [0.0f64; 3];
    let mut err = [0.0f64; 3];
    for (l, w) in rows.iter().zip(&want) {
        let (slot, row): (usize, u32) = (l[0].parse().unwrap(), l[1].parse().unwrap());
        assert_eq!(perm[slot - 1] + 1, row, "Perm({slot})");
        let s = perm[slot - 1] as usize * 3;
        for c in 0..3 {
            scale[c] = scale[c].max(w[c].abs());
            err[c] = err[c].max((u.values[s + c] - w[c]).abs());
        }
    }
    let norm = u.conv.norm;
    eprintln!("{case}: norm {norm:.17e}, max error / scale {:?}", [0, 1, 2].map(|c| err[c] / scale[c]));
    assert!((norm - ref_norm).abs() / ref_norm < 1e-10, "{case}: norm {norm:e}");
    for c in 0..3 {
        assert!(err[c] <= 1e-10 * scale[c], "{case}: component {c} max error {} (scale {})", err[c], scale[c]);
    }
}
