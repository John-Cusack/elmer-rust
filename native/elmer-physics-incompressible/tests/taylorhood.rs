//! The Taylor-Hood pair (`Element = "p:2"`, no `Pressure Stabilization`)
//! against the pristine ElmerSolver: `StokesStabilized/quad_taylorhood.sif`
//! `tri_taylorhood.sif` and `tet_taylorhood.sif` with a direct solve
//! (`tests/data/README.md`), the p dofs from `elmer-pdofs`. Serendipity
//! quads and p:2 triangles and tetrahedra have no bubbles: the velocity is quadratic, the edge pressures are pinned to
//! zero. Every slot's row, velocity and pressure of the `case.result`, and
//! the norm the run printed (velocity only).

mod support;

use elmer_core::{NullLogger, ParamMap, Params, SerialComm, Value};
use elmer_pdofs::PDofsBackend;
use elmer_physics_incompressible::{register, solver_defaults};
use elmer_solver::{Simulation, SolverRegistry};
use support::cases::{data, ints, model, ns_solver, profile, real, sec, st, steady, table};
use support::mesh::TestMesh;

/// The whitespace-split lines of a fixture, `#` comments skipped.
fn lines(path: &std::path::Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .filter(|l: &Vec<String>| !l.is_empty())
        .collect()
}

/// `Element = p:2`, the Taylor-Hood pair.
fn p2() -> [(&'static str, Value); 1] {
    [("Element", st("p:2"))]
}

fn bc(tags: &[i32], entries: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![("Target Boundaries", ints(tags))];
    v.extend(entries.iter().cloned());
    sec(&v)
}

/// The backward facing step of `quad_taylorhood.sif` / `tri_taylorhood.sif`.
fn step(mesh: &TestMesh) -> (Vec<ParamMap>, f64) {
    let bcs = vec![
        bc(&[1], &[("Velocity 1", profile(mesh, 1, 2, |y| (2.0 - y) * (y - 1.0))), ("Velocity 2", real(0.0))]),
        bc(&[2], &[("Velocity 2", real(0.0))]),
        bc(&[3], &[("Velocity 1", real(0.0)), ("Velocity 2", real(0.0))]),
    ];
    (bcs, 1.0e3)
}

/// The step with loads on the outflow: `External Pressure` 100 at y = 0 to
/// 300 at y = 2 along the outward normal and `Surface Traction 1 = 10`.
fn step_loads(mesh: &TestMesh) -> (Vec<ParamMap>, f64) {
    let (mut bcs, viscosity) = step(mesh);
    bcs[1].insert("External Pressure", table(2, &[(0.0, 100.0), (2.0, 300.0)]));
    bcs[1].insert("Surface Traction 1", real(10.0));
    (bcs, viscosity)
}

/// The step of `tests/data/ns_tri` (Navier-Stokes): viscosity 0.1.
fn ns_step(mesh: &TestMesh) -> (Vec<ParamMap>, f64) {
    (step(mesh).0, 0.1)
}

/// The driven cavity of `tet_taylorhood.sif`.
fn cavity(mesh: &TestMesh) -> (Vec<ParamMap>, f64) {
    let bcs = vec![
        bc(
            &[5],
            &[
                ("Velocity 1", profile(mesh, 5, 2, |y| y * (1.0 - y))),
                ("Velocity 2", real(0.0)),
                ("Velocity 3", real(0.0)),
            ],
        ),
        bc(&[6], &[("Velocity 1", real(0.0)), ("Velocity 2", real(0.0)), ("Velocity 3", real(0.0))]),
    ];
    (bcs, 1.0)
}

type Setup = fn(&TestMesh) -> (Vec<ParamMap>, f64);

fn run(
    case: &str,
    dim: usize,
    setup: Setup,
    simulation: ParamMap,
    extra: &[(&str, Value)],
) -> Result<Simulation, elmer_solver::SolverError> {
    let dir = data(case);
    let mesh = TestMesh::read(&dir, dim);
    let mut keys = vec![
        ("Stokes Flow", Value::Logical(true)),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-3)),
    ];
    keys.extend(extra.iter().cloned());
    let mut solver = ns_solver(&keys);
    // _Init0 (the MINI element without an Element) before the dofs are made.
    let transient = simulation.string("simulation type").is_some_and(|t| t.eq_ignore_ascii_case("transient"));
    solver_defaults(&mut solver, dim, transient);
    let (bcs, viscosity) = setup(&mesh);
    let material = sec(&[("Density", real(1.0)), ("Viscosity", real(viscosity))]);
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(model(simulation, solver, material, None, bcs));
    sim.backends.dofs = Some(Box::new(PDofsBackend::new()));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(&mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

fn run_case(case: &str, dim: usize, setup: Setup, extra: &[(&str, Value)], n_slots: usize, ref_norm: f64) {
    let csys = if dim == 3 { "Cartesian 3D" } else { "Cartesian 2D" };
    run_case_with(case, dim, setup, steady(csys), extra, n_slots, ref_norm, dim + 1);
}

#[allow(clippy::too_many_arguments)]
fn run_case_with(
    case: &str,
    dim: usize,
    setup: Setup,
    simulation: ParamMap,
    extra: &[(&str, Value)],
    n_slots: usize,
    ref_norm: f64,
    compared: usize,
) {
    let dir = data(case);
    let sim = run(case, dim, setup, simulation, extra).unwrap();
    let dofs = dim + 1;

    let u = sim.vars.get("flow solution").unwrap();
    assert_eq!(u.values.len(), dofs * n_slots, "{case}: dofs");
    let norm = u.conv.norm;

    // perm.dat: `slot row`; solution.dat: `slot u1 .. p`, x(Perm(slot)).
    let perm = u.perm.as_deref().unwrap();
    let rows = lines(&dir.join("perm.dat"));
    assert_eq!(rows.len(), n_slots, "{case}: Perm");
    for l in &rows {
        let (slot, row): (usize, u32) = (l[0].parse().unwrap(), l[1].parse().unwrap());
        assert_eq!(perm[slot - 1] + 1, row, "{case}: Perm({slot})");
    }
    let want: Vec<Vec<f64>> = lines(&dir.join("solution.dat"))
        .iter()
        .map(|l| l[1..].iter().map(|x| x.parse().unwrap()).collect())
        .collect();
    assert_eq!(want.len(), n_slots);
    let mut scale = vec![0.0f64; dofs];
    let mut err = vec![0.0f64; dofs];
    for (slot, w) in want.iter().enumerate() {
        assert_eq!(w.len(), dofs);
        let s = perm[slot] as usize * dofs;
        for c in 0..dofs {
            scale[c] = scale[c].max(w[c].abs());
            err[c] = err[c].max((u.values[s + c] - w[c]).abs());
        }
    }
    let rel: Vec<f64> = (0..dofs).map(|c| err[c] / scale[c]).collect();
    eprintln!("{case}: norm {norm:.17e}, max error / scale per component {rel:?}");
    assert!((norm - ref_norm).abs() / ref_norm < 1e-10, "{case}: norm {norm:.17e} vs {ref_norm:.17e}");
    for c in 0..compared {
        assert!(err[c] <= 1e-10 * scale[c], "{case}: component {c} max error {} (scale {})", err[c], scale[c]);
    }
}

/// 531 serendipity quads: 595 nodes + 1125 edges.
#[test]
fn taylorhood_quad_matches_elmer() {
    run_case("taylorhood_quad", 2, step, &p2(), 1720, 0.04399609386373855);
}

/// 1062 triangles: 595 nodes + 1656 edges.
#[test]
fn taylorhood_tri_matches_elmer() {
    run_case("taylorhood_tri", 2, step, &p2(), 2251, 0.03851002030822979);
}

/// 386 p:2 tetrahedra (143 nodes + 660 edges) on the 5-point rule of
/// `Element Integration Points = "-tetra 5"`.
#[test]
fn taylorhood_tet_matches_elmer() {
    let extra = [("Element", st("p:2")), ("Element Integration Points", st("-tetra 5"))];
    run_case("taylorhood_tet", 3, cavity, &extra, 803, 0.025134495063016765);
}

/// Full p:2 quads (`Serendipity P Elements = False`) have a bubble each,
/// which a Stokes run condenses (`LCondensate`) without pinning the edge
/// pressures: 595 nodes + 1125 edges in the system. The quadratic pressure
/// then has spurious modes the direct solvers resolve differently, so the
/// velocity (and its norm) is compared, not the pressure.
#[test]
fn taylorhood_quad_full_matches_elmer() {
    let mut simulation = steady("Cartesian 2D");
    simulation.insert("Serendipity P Elements", Value::Logical(false));
    run_case_with("taylorhood_quad_full", 2, step, simulation, &p2(), 1720, 0.04394492021142661, 2);
}

/// `Relative Integration Order` (`GaussPointsAdapt`) is not ported: refused.
#[test]
fn relative_integration_order_is_refused() {
    let extra = [("Element", st("p:2")), ("Relative Integration Order", ints(&[1]))];
    let err = run("taylorhood_quad", 2, step, steady("Cartesian 2D"), &extra).err().expect("refused");
    assert!(err.to_string().contains("Relative Integration Order"), "{err}");
}

/// `StokesStabilized/quad.sif` and `tri.sif` without `Pressure
/// Stabilization`: `_Init0`'s MINI element (`p:1 -tri b:1 -quad b:3`), the
/// velocity bubbles condensed element by element (`LCondensate`), so only
/// the 595 nodes are in the system.
#[test]
fn mini_quad_matches_elmer() {
    run_case("mini_quad", 2, step, &[], 595, 0.07554915049888652);
}

#[test]
fn mini_tri_matches_elmer() {
    run_case("mini_tri", 2, step, &[], 595, 0.07587897980190471);
}

/// `tests/data/ns_tri` without `Pressure Stabilization`: Navier-Stokes on
/// the MINI element, Picard then Newton after 3 iterations, the bubbles
/// condensed in every iteration (the convection taken from the nodal
/// velocity, the gradient from all functions with the bubbles at zero).
#[test]
fn mini_navier_stokes_tri_matches_elmer() {
    let extra = [
        ("Stokes Flow", Value::Logical(false)),
        ("Nonlinear System Max Iterations", ints(&[30])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
        ("Nonlinear System Newton After Iterations", ints(&[3])),
        ("Nonlinear System Newton After Tolerance", real(1.0e-3)),
    ];
    run_case("mini_ns_tri", 2, ns_step, &extra, 595, 0.07592011845809431);
}

/// `tests/data/transient_tri` without `Pressure Stabilization`: the MINI
/// element over 4 steps of 0.1, Newton after 2 iterations. `LCondensate`
/// forms the time derivative of the bubble-augmented system itself (BDF1,
/// though the run asks for BDF2) with each element's bubble history.
#[test]
fn mini_transient_tri_matches_elmer() {
    let extra = transient_ns();
    run_case_with("mini_transient_tri", 2, ns_step, transient_bdf2(), &extra, 595, 0.07440496513539256, 3);
}

/// `transient_tri`'s time stepping: 4 BDF2 steps of 0.1.
fn transient_bdf2() -> ParamMap {
    let mut simulation = steady("Cartesian 2D");
    for (k, v) in [
        ("Simulation Type", st("Transient")),
        ("Timestepping Method", st("BDF")),
        ("BDF Order", ints(&[2])),
        ("Timestep Sizes", real(0.1)),
        ("Timestep Intervals", ints(&[4])),
    ] {
        simulation.insert(k, v);
    }
    simulation
}

/// `transient_tri`'s Navier-Stokes iteration: 3, Newton after 2.
fn transient_ns() -> Vec<(&'static str, Value)> {
    vec![
        ("Stokes Flow", Value::Logical(false)),
        ("Nonlinear System Max Iterations", ints(&[3])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
        ("Nonlinear System Newton After Iterations", ints(&[2])),
    ]
}

/// `transient_tri` with `Element = p:2` and no stabilization: Taylor-Hood
/// Navier-Stokes over 4 BDF2 steps, `Default1stOrderTime` on the mass
/// matrix with the edge pressures pinned (595 nodes + 1656 edges).
#[test]
fn taylorhood_transient_tri_matches_elmer() {
    let mut extra = transient_ns();
    extra.push(("Element", st("p:2")));
    run_case_with("taylorhood_transient_tri", 2, ns_step, transient_bdf2(), &extra, 2251, 0.037828652451776607, 3);
}

/// `tests/data/thermal_quad/case.sif`: ThermalCompress_vec's cavity with
/// `Compressibility Model = Thermal` and the temperature an exported
/// variable held at its initial condition `293.5 - x` (a two-point table,
/// exact on the unit square): rho(T) at the integration points and the
/// `-(1/rho) d(rho)/dx_i u_i` continuity term, on the MINI element.
#[test]
fn thermal_compressibility_matches_elmer() {
    let material = sec(&[
        ("Density", real(1.0)),
        ("Viscosity", real(1.0)),
        ("Reference Temperature", real(293.0)),
        ("Compressibility Model", st("Thermal")),
        ("Heat Expansion Coefficient", real(0.25)),
    ]);
    compressible_cavity("thermal_quad", material, 0.00445569022547017);
}

/// `tests/data/perfect_gas_quad/case.sif`: the same cavity with
/// `Compressibility Model = Perfect Gas` (`Heat Capacity` 1000, `Specific
/// Heat Ratio` 5/3, `Reference Pressure` 1e5): rho from the previous
/// iterate's pressure and the temperature, through the Picard iterations.
#[test]
fn perfect_gas_compressibility_matches_elmer() {
    let material = sec(&[
        ("Density", real(1.0)),
        ("Viscosity", real(1.0)),
        ("Heat Capacity", real(1000.0)),
        ("Reference Temperature", real(293.0)),
        ("Compressibility Model", st("Perfect Gas")),
        ("Specific Heat Ratio", real(5.0 / 3.0)),
        ("Reference Pressure", real(1.0e5)),
    ]);
    compressible_cavity("perfect_gas_quad", material, 5.183198317697339e-05);
}

/// ThermalCompress_vec's cavity on the MINI element with `material`, the
/// temperature an exported variable at its initial condition `293.5 - x`,
/// against the pristine `case`.
fn compressible_cavity(case: &str, material: ParamMap, ref_norm: f64) {
    let dir = data(case);
    let mesh = TestMesh::read(&dir, 2);
    let mut solver = ns_solver(&[
        ("Equation", st("Navier-Stokes")),
        ("Exported Variable 1", st("Temperature")),
        ("Nonlinear System Max Iterations", ints(&[10])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
        ("Nonlinear System Newton After Iterations", ints(&[20])),
        ("Nonlinear System Newton After Tolerance", real(1.0e-3)),
    ]);
    solver_defaults(&mut solver, 2, false);
    let wall = |t: i32| bc(&[t], &[("Velocity 1", real(0.0)), ("Velocity 2", real(0.0))]);
    let mut first = wall(1);
    first.insert("Pressure", real(0.0));
    let bcs = vec![first, wall(2), wall(3), wall(4)];
    let body_force = sec(&[("Flow BodyForce 2", real(-10.0))]);
    let mut m = model(steady("Cartesian 2D"), solver, material, Some(body_force), bcs);
    m.bodies[0].insert("Initial Condition", ints(&[1]));
    m.initial_conditions = vec![sec(&[("Temperature", table(1, &[(0.0, 293.5), (1.0, 292.5)]))])];
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(m);
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

/// `fem/tests/IncompressibleNSAxi` (`tests/data/axi_pipe`): axisymmetric
/// Stokes with Pressure Stabilization on `Element = "p:1"`, `Specific Load`
/// `Flow Bodyforce 2 = -1/r`. Its MATC conditions are linear along each
/// boundary (two-point tables, exact there) and the load a table through
/// the nodal radii (the values it is interpolated from).
#[test]
fn axisymmetric_stokes_matches_elmer() {
    axi_pipe("axi_pipe", false, 1.3600542597525465);
}

/// `tests/data/axi_pipe_loads`: the outer wall (boundary 2, r = 1.5)
/// loaded instead of held: `External Pressure` 1 to 3 along z, `Surface
/// Traction 2 = 0.5`, `Slip Coefficient 2 = 4`, integrated over r ds.
#[test]
fn axisymmetric_boundary_loads_match_elmer() {
    axi_pipe("axi_pipe_loads", true, 1.3245395366921995);
}

fn axi_pipe(case: &str, loads: bool, ref_norm: f64) {
    let dir = data(case);
    let mesh = TestMesh::read(&dir, 2);
    let mut solver = ns_solver(&[
        ("Equation", st("Navier-Stokes")),
        ("Stokes Flow", Value::Logical(true)),
        ("Specific Load", Value::Logical(true)),
        ("Element", st("p:1")),
        ("Pressure Stabilization", Value::Logical(true)),
        ("Nonlinear System Max Iterations", ints(&[1])),
    ]);
    solver_defaults(&mut solver, 2, false);
    let material = sec(&[("Density", real(1.0)), ("Viscosity", real(1.0))]);
    let v1 = || table(1, &[(0.5, 0.5), (1.5, 1.5)]);
    let bcs = vec![
        bc(&[1], &[("Velocity 1", v1()), ("Velocity 2", table(1, &[(0.5, 0.5), (1.5, 1.5)]))]),
        if loads {
            bc(
                &[2],
                &[
                    ("External Pressure", table(2, &[(0.0, 1.0), (2.0, 3.0)])),
                    ("Surface Traction 2", real(0.5)),
                    ("Slip Coefficient 2", real(4.0)),
                ],
            )
        } else {
            bc(&[2], &[("Velocity 1", v1()), ("Velocity 2", table(2, &[(0.0, 1.5), (2.0, -2.5)]))])
        },
        bc(
            &[3],
            &[("Velocity 1", v1()), ("Velocity 2", table(1, &[(0.5, -3.5), (1.5, -2.5)])), ("Pressure", real(0.0))],
        ),
        bc(&[4], &[("Velocity 1", v1()), ("Velocity 2", table(2, &[(0.0, 0.5), (2.0, -3.5)]))]),
    ];
    let mut radii: Vec<f64> = mesh.coords.iter().map(|c| c[0]).collect();
    radii.sort_by(f64::total_cmp);
    radii.dedup();
    let load = table(1, &radii.iter().map(|&r| (r, -1.0 / r)).collect::<Vec<_>>());
    let body_force = sec(&[("Flow Bodyforce 1", real(0.0)), ("Flow Bodyforce 2", load)]);
    let m = model(steady("Axi Symmetric"), solver, material, Some(body_force), bcs);
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(m);
    sim.backends.dofs = Some(Box::new(PDofsBackend::new()));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();

    let u = sim.vars.get("flow solution").unwrap();
    let perm = u.perm.as_deref().unwrap();
    let rows = lines(&dir.join("perm.dat"));
    assert_eq!(rows.len(), 187);
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
        assert!(err[c] <= 1e-10 * scale[c], "component {c} max error {} (scale {})", err[c], scale[c]);
    }
}

/// `tests/data/taylorhood_quad_loads`: Taylor-Hood with the outflow loads,
/// integrated over the p:2 boundary edges (their edge functions take the
/// load as well).
#[test]
fn taylorhood_boundary_loads_match_elmer() {
    run_case("taylorhood_quad_loads", 2, step_loads, &p2(), 1720, 0.04399059383796517);
}

/// `tests/data/ns_tri` with a non-newtonian `Viscosity Model` (Picard only,
/// the Newton keywords removed; `tests/data/visc_*/case.sif`): the
/// viscosity at the points from the previous iterate's shear rate, under
/// the Brezzi-Pitkaranta stabilization (its tau reads it too).
fn viscosity_case(case: &str, model: &[(&str, Value)], ref_norm: f64) {
    let dir = data(case);
    let mesh = TestMesh::read(&dir, 2);
    let mut solver = ns_solver(&[
        ("Pressure Stabilization", Value::Logical(true)),
        ("Pressure Stabilization Element Size", st("longest edge")),
        ("Nonlinear System Max Iterations", ints(&[60])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-10)),
    ]);
    solver_defaults(&mut solver, 2, false);
    let mut material = vec![("Density", real(1.0)), ("Viscosity", real(0.1))];
    material.extend(model.iter().cloned());
    let m = model_of(steady("Cartesian 2D"), solver, sec(&material), step(&mesh).0);
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(m);
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();

    let u = sim.vars.get("flow solution").unwrap();
    let perm = u.perm.as_deref().unwrap();
    let rows = lines(&dir.join("perm.dat"));
    assert_eq!(rows.len(), 595);
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

fn model_of(simulation: ParamMap, solver: ParamMap, material: ParamMap, bcs: Vec<ParamMap>) -> elmer_solver::Model {
    model(simulation, solver, material, None, bcs)
}

#[test]
fn power_law_viscosity_matches_elmer() {
    let m = [
        ("Viscosity Model", st("power law")),
        ("Viscosity Exponent", real(0.7)),
        ("Critical Shear Rate", real(0.05)),
    ];
    viscosity_case("visc_powerlaw", &m, 0.07532484369893568);
}

#[test]
fn carreau_viscosity_matches_elmer() {
    let m = [
        ("Viscosity Model", st("carreau")),
        ("Viscosity Difference", real(0.2)),
        ("Viscosity Exponent", real(0.5)),
        ("Viscosity Transition", real(1.5)),
    ];
    viscosity_case("visc_carreau", &m, 0.07577798970745593);
}

#[test]
fn cross_viscosity_matches_elmer() {
    let m = [
        ("Viscosity Model", st("cross")),
        ("Viscosity Difference", real(0.2)),
        ("Viscosity Exponent", real(1.3)),
        ("Viscosity Transition", real(2.0)),
    ];
    viscosity_case("visc_cross", &m, 0.07552704554913707);
}

/// Not converged after the 60 Picard iterations (as in the pristine run,
/// change 1.7e-6): the 60th iterate is compared.
#[test]
fn powell_eyring_viscosity_matches_elmer() {
    let m = [
        ("Viscosity Model", st("powell eyring")),
        ("Viscosity Difference", real(0.05)),
        ("Viscosity Transition", real(1.5)),
    ];
    viscosity_case("visc_powell", &m, 0.07420013848890271);
}
