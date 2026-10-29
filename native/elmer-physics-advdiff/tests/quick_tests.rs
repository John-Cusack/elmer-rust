//! The quick tests that use the two solvers, run end to end through
//! `elmer-solver`'s `Simulation` and compared with the full-precision
//! solutions of the pristine ElmerSolver (see `data/README.md`):
//! `adv_diff1`, `adv_diff2`, `AdvDiffSoret`, `AdvDiffSoretAbsMass`, and a
//! linear-element variant of `ModelPDE` (the quick test itself uses p:2).

mod support;

use elmer_core::{ParamMap, Value};
use elmer_solver::Model;
use support::backends::SimMesh;
use support::cases::{data, ints, max_diff, nodal, read_solution, real, rms, run, sec, st};
use support::meshio::read_elmer_mesh;

fn iterative_ilu1() -> Vec<(&'static str, Value)> {
    vec![
        ("Linear System Solver", st("Iterative")),
        ("Linear System Iterative Method", st("BiCGStab")),
        ("Linear System Max Iterations", ints(&[200])),
        ("Linear System Convergence Tolerance", real(1.0e-8)),
        ("Linear System Preconditioning", st("ILU1")),
        ("Linear System Residual Output", ints(&[0])),
    ]
}

fn advdiff_solver(variable: &str, extra: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![
        ("Equation", st("Advection Diffusion")),
        (
            "Procedure",
            st("\"AdvectionDiffusion\" \"AdvectionDiffusionSolver\""),
        ),
        ("Variable", st(variable)),
        ("Variable DOFs", ints(&[1])),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Steady State Convergence Tolerance", real(1.0e-4)),
    ];
    v.extend(extra.iter().cloned());
    sec(&v)
}

fn transient(dt: f64, steps: i32) -> ParamMap {
    sec(&[
        ("Coordinate System", st("Cartesian 2D")),
        ("Simulation Type", st("transient")),
        ("Timestepping Method", st("BDF")),
        ("BDF Order", ints(&[2])),
        ("Timestep Sizes", real(dt)),
        ("Timestep Intervals", ints(&[steps])),
        ("Output Intervals", ints(&[0])),
        ("Steady State Max Iterations", ints(&[2])),
    ])
}

fn body() -> ParamMap {
    sec(&[
        ("Target Bodies", ints(&[1])),
        ("Equation", ints(&[1])),
        ("Material", ints(&[1])),
    ])
}

/// Tolerance on pointwise differences from the pristine solutions: Elmer's own
/// linear solves stop at 1e-8, so the references carry that much noise.
const TOL: f64 = 5.0e-7;

#[test]
fn adv_diff1_matches_elmer() {
    // 408 elements, no convection (so no bubbles), BDF2, ILU1 BiCGStab.
    let mesh = SimMesh::new(read_elmer_mesh(&data("adv_diff1"), 2));
    let n = mesh.mesh.coords.len();
    let mut v = iterative_ilu1();
    v.push(("Bubbles", Value::Logical(false)));
    let model = Model {
        simulation: transient(2.0, 10),
        solvers: vec![advdiff_solver("Poison", &v)],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![body()],
        materials: vec![sec(&[
            ("Density", real(1.3)),
            ("Poison Diffusivity", real(1.0e-3)),
        ])],
        boundary_conditions: vec![sec(&[
            ("Target Boundaries", ints(&[1])),
            ("Poison", real(1.0)),
        ])],
        ..Model::default()
    };
    let sim = run(model, &mesh).unwrap();
    let got = nodal(&sim, "poison", n);
    let want = read_solution(&data("adv_diff1").join("solution.dat"));
    assert!(max_diff(&got, &want) < TOL, "{:e}", max_diff(&got, &want));
    assert!(
        (rms(&got) - 0.42696309972986873).abs() < 1e-7,
        "{}",
        rms(&got)
    );
}

/// `adv_diff2`'s model: 404 elements, absolute mass, constant convection and
/// a mass-transfer boundary. `bubbles` is the Elmer default; `steady` drops
/// the BDF2 time stepping.
fn adv_diff2_model(method: Method, steady: bool) -> Model {
    let mut v = iterative_ilu1();
    match method {
        Method::Bubbles => {}
        Method::Galerkin => v.push(("Bubbles", Value::Logical(false))),
        Method::Supg => v.push(("Stabilize", Value::Logical(true))),
    }
    let simulation = if steady {
        sec(&[
            ("Coordinate System", st("Cartesian 2D")),
            ("Simulation Type", st("steady state")),
            ("Output Intervals", ints(&[0])),
            ("Steady State Max Iterations", ints(&[2])),
        ])
    } else {
        transient(0.1, 5)
    };
    Model {
        simulation,
        solvers: vec![advdiff_solver("Poison", &v)],
        equations: vec![sec(&[
            ("Active Solvers", ints(&[1])),
            ("Concentration Units", st("Absolute Mass")),
            ("Convection", st("Constant")),
        ])],
        bodies: vec![body()],
        materials: vec![sec(&[
            ("Density", real(1.3)),
            ("Poison Diffusivity", real(1.0e-3)),
            ("Convection Velocity 1", real(1.0)),
            ("Convection Velocity 2", real(0.0)),
            ("Convection Velocity 3", real(0.0)),
        ])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", ints(&[1])), ("Poison", real(1.0))]),
            sec(&[
                ("Target Boundaries", ints(&[2])),
                ("Mass Transfer Coefficient", real(0.1)),
                ("External Concentration", real(0.0)),
            ]),
        ],
        ..Model::default()
    }
}

#[derive(Copy, Clone)]
enum Method {
    Bubbles,
    Galerkin,
    Supg,
}

fn check_adv_diff2(method: Method, steady: bool, mesh_dir: &str, file: &str, norm: f64) {
    let mesh = SimMesh::new(read_elmer_mesh(&data(mesh_dir), 2));
    let n = mesh.mesh.coords.len();
    let sim = run(adv_diff2_model(method, steady), &mesh).unwrap();
    let got = nodal(&sim, "poison", n);
    let want = read_solution(&data(mesh_dir).join(file));
    assert!(max_diff(&got, &want) < TOL, "{:e}", max_diff(&got, &want));
    assert!((rms(&got) - norm).abs() < 1e-7, "{}", rms(&got));
}

#[test]
fn adv_diff2_matches_elmer() {
    // The default residual-free bubbles, condensed after BDF2 on the nodal
    // block, with a mass-transfer boundary.
    check_adv_diff2(
        Method::Bubbles,
        false,
        "adv_diff2",
        "solution.dat",
        0.8766881677902325,
    );
}

#[test]
fn adv_diff2_galerkin_convection_matches_elmer() {
    check_adv_diff2(
        Method::Galerkin,
        false,
        "adv_diff2",
        "nb_solution.dat",
        0.8691287700487503,
    );
}

#[test]
fn adv_diff2_steady_bubbles_match_elmer() {
    check_adv_diff2(
        Method::Bubbles,
        true,
        "adv_diff2",
        "st_solution.dat",
        0.9980513902691243,
    );
}

#[test]
fn adv_diff2_supg_on_triangles_matches_elmer() {
    // `Stabilize = True` on 303 elements, with `hK` and `mK` from elmer-fe's
    // StabParam; second derivatives of linear triangles vanish.
    check_adv_diff2(
        Method::Supg,
        false,
        "adv_diff2_tri",
        "solution.dat",
        0.8844887712095778,
    );
}

fn soret_model(absolute: bool, with_soret: bool) -> Model {
    let mut eq = vec![
        ("Active Solvers", ints(&[1, 2])),
        // Per-variable convection keeps HeatSolve convection-free; the
        // velocity is zero in the case either way.
        ("Conc Convection", st("Constant")),
    ];
    if absolute {
        eq.push(("Concentration Units", st("Absolute Mass")));
    }
    let mut mat = vec![
        ("Density", real(1000.0)),
        ("Heat Conductivity", real(0.6)),
        ("Heat Capacity", real(4000.0)),
        ("Conc Diffusivity", real(1.0e-9)),
        ("Temperature Field Variable", st("Temperature")),
        ("Reference Temperature", real(300.0)),
        ("Convection Velocity 1", real(0.0)),
        ("Convection Velocity 2", real(0.0)),
    ];
    if with_soret {
        mat.push(("Conc Soret Diffusivity", real(5.0e-12)));
    }
    let mut heat = vec![
        ("Equation", st("Heat Equation")),
        ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
        ("Variable", st("Temperature")),
        ("Variable DOFs", ints(&[1])),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Steady State Convergence Tolerance", real(1.0e-10)),
    ];
    heat.push(("Linear System Solver", st("Direct")));
    let mut ad = vec![("Bubbles", Value::Logical(false))];
    ad.push(("Linear System Solver", st("Direct")));
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian 2D")),
            ("Simulation Type", st("steady state")),
            ("Steady State Max Iterations", ints(&[2])),
            ("Output Intervals", ints(&[0])),
        ]),
        solvers: vec![sec(&heat), advdiff_solver("Conc", &ad)],
        equations: vec![sec(&eq)],
        bodies: vec![body()],
        materials: vec![sec(&mat)],
        boundary_conditions: vec![
            sec(&[
                ("Target Boundaries", ints(&[1])),
                ("Temperature", real(300.0)),
                ("Conc", real(1.0)),
            ]),
            sec(&[
                ("Target Boundaries", ints(&[3])),
                ("Temperature", real(400.0)),
            ]),
        ],
        ..Model::default()
    }
}

fn check_soret(absolute: bool) {
    let mesh = SimMesh::new(read_elmer_mesh(&data("adv_diff_soret"), 2));
    let n = mesh.mesh.coords.len();
    let sim = run(soret_model(absolute, true), &mesh).unwrap();
    let t = nodal(&sim, "temperature", n);
    let c = nodal(&sim, "conc", n);
    let t_want = read_solution(&data("adv_diff_soret").join("temp.dat"));
    let file = if absolute { "conc_abs.dat" } else { "conc.dat" };
    let c_want = read_solution(&data("adv_diff_soret").join(file));
    assert!(
        max_diff(&t, &t_want) < 1e-9 * 400.0,
        "{:e}",
        max_diff(&t, &t_want)
    );
    assert!(max_diff(&c, &c_want) < 1e-10, "{:e}", max_diff(&c, &c_want));
    assert!((rms(&c) - 0.7638989898318287).abs() < 1e-10, "{}", rms(&c));

    // The Soret flux is what moves the solution off the plain diffusion one.
    let plain = run(soret_model(absolute, false), &mesh).unwrap();
    let c0 = nodal(&plain, "conc", n);
    assert!(max_diff(&c, &c0) > 1e-4, "{:e}", max_diff(&c, &c0));
}

#[test]
fn adv_diff_soret_relative_mass_matches_elmer() {
    check_soret(false);
}

#[test]
fn adv_diff_soret_absolute_mass_matches_elmer() {
    check_soret(true);
}

fn model_pde_model(transient_case: bool) -> Model {
    let mut sim = vec![
        ("Coordinate System", st("Cartesian")),
        ("Steady State Max Iterations", ints(&[1])),
        ("Output Intervals", ints(&[0])),
    ];
    let mut mat = vec![
        ("Convection Velocity 1", real(1.0)),
        ("Convection Velocity 2", real(0.0)),
        ("diffusion coefficient", real(1.0)),
        ("convection coefficient", real(9.0)),
    ];
    if transient_case {
        sim.extend([
            ("Simulation Type", st("transient")),
            ("Timestepping Method", st("BDF")),
            ("BDF Order", ints(&[2])),
            ("Timestep Sizes", real(0.05)),
            ("Timestep Intervals", ints(&[6])),
        ]);
        mat.push(("time derivative coefficient", real(2.0)));
        mat.push(("reaction coefficient", real(0.5)));
    } else {
        sim.push(("Simulation Type", st("steady")));
        mat.push(("time derivative coefficient", real(0.0)));
    }
    Model {
        simulation: sec(&sim),
        solvers: vec![sec(&[
            ("Equation", st("ModelPDE")),
            ("Variable", st("Field")),
            ("Procedure", st("\"ModelPDE\" \"AdvDiffSolver\"")),
            ("Linear System Solver", st("Direct")),
            ("Steady State Convergence Tolerance", real(1e-9)),
        ])],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[
            ("Target Bodies", ints(&[1])),
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
            ("Body Force", ints(&[1])),
        ])],
        materials: vec![sec(&mat)],
        body_forces: vec![sec(&[("Field Source", real(1.0))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", ints(&[1])), ("Field", real(0.0))]),
            sec(&[
                ("Target Boundaries", ints(&[2])),
                ("Robin Coefficient", real(10.0)),
                ("External Field", real(5.0)),
            ]),
            sec(&[
                ("Target Boundaries", ints(&[3, 4])),
                ("Field Flux", real(10.0)),
            ]),
        ],
        ..Model::default()
    }
}

#[test]
fn model_pde_steady_matches_elmer() {
    // The quick test's case with linear elements instead of `Element = p:2`.
    let mesh = SimMesh::new(read_elmer_mesh(&data("model_pde"), 2));
    let n = mesh.mesh.coords.len();
    let sim = run(model_pde_model(false), &mesh).unwrap();
    let got = nodal(&sim, "field", n);
    let want = read_solution(&data("model_pde").join("lin_solution.dat"));
    assert!(max_diff(&got, &want) < 1e-10, "{:e}", max_diff(&got, &want));
    assert!((rms(&got) - 1.7972256169743914).abs() < 1e-11);
}

#[test]
fn model_pde_number_of_integration_points_on_nodal_elements() {
    // GaussPointsAdapt with Number of Integration Points = 4 on bilinear
    // quads is the type's default rule: the linear case's solution.
    let mesh = SimMesh::new(read_elmer_mesh(&data("model_pde"), 2));
    let n = mesh.mesh.coords.len();
    let mut model = model_pde_model(false);
    model.solvers[0].insert("Number of Integration Points", ints(&[4]));
    let sim = run(model, &mesh).unwrap();
    let got = nodal(&sim, "field", n);
    let want = read_solution(&data("model_pde").join("lin_solution.dat"));
    assert!(max_diff(&got, &want) < 1e-10, "{:e}", max_diff(&got, &want));
}

#[test]
fn model_pde_transient_with_reaction_matches_elmer() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("model_pde"), 2));
    let n = mesh.mesh.coords.len();
    let sim = run(model_pde_model(true), &mesh).unwrap();
    let got = nodal(&sim, "field", n);
    let want = read_solution(&data("model_pde").join("evol_solution.dat"));
    assert!(max_diff(&got, &want) < 1e-10, "{:e}", max_diff(&got, &want));
}

/// The quick tests `ModelPDE` (`Element = "p:2"`, full elements) and
/// `SD_ModelPDE` (the same on serendipity elements), with the p dofs from
/// `elmer-pdofs`: Perm, every slot's value and the norm against the pristine
/// runs (`p2.sif`, `sd_p2.sif`).
fn model_pde_p(serendipity: bool) -> (Vec<u32>, Vec<f64>, f64) {
    let mesh = SimMesh::new(read_elmer_mesh(&data("model_pde"), 2));
    let mut model = model_pde_model(false);
    model.solvers[0].insert("Element", st("p:2"));
    if !serendipity {
        model
            .simulation
            .insert("Serendipity P Elements", Value::Logical(false));
    }
    let (mut sim, reg) = support::cases::simulation(model);
    sim.backends.dofs = Some(Box::new(elmer_pdofs::PDofsBackend::new()));
    sim.setup(
        &mesh,
        &reg,
        &elmer_core::NullLogger,
        &elmer_core::SerialComm,
    )
    .unwrap();
    sim.run(&mesh, &elmer_core::NullLogger, &elmer_core::SerialComm)
        .unwrap();
    let v = sim.vars.get("field").unwrap();
    let perm = v.perm.as_deref().unwrap().to_vec();
    (perm, v.values.clone(), v.conv.norm)
}

fn check_model_pde_p(serendipity: bool, case: &str, n_dofs: usize, norm: f64) {
    let (perm, values, got_norm) = model_pde_p(serendipity);
    assert_eq!(values.len(), n_dofs, "{case}");
    let want_perm = read_solution(&data("model_pde").join(format!("{case}_perm.dat")));
    let ours: Vec<f64> = perm.iter().map(|&p| f64::from(p + 1)).collect();
    assert_eq!(ours, want_perm, "{case}: Perm");
    // x(Perm(i)) slot by slot, from the run's case.result.
    let want = read_solution(&data("model_pde").join(format!("{case}_solution.dat")));
    let got: Vec<f64> = perm.iter().map(|&p| values[p as usize]).collect();
    assert!(
        max_diff(&got, &want) < 1e-10,
        "{case}: {:e}",
        max_diff(&got, &want)
    );
    assert!((got_norm - norm).abs() < 1e-10, "{case}: norm {got_norm}");
}

#[test]
fn model_pde_p2_matches_elmer() {
    // 441 nodes + 840 edges + 400 global bubbles; Reference Norm 9.23486039E-01.
    check_model_pde_p(false, "p2", 1681, 0.9234860388553119);
}

#[test]
fn sd_model_pde_p2_matches_elmer() {
    // Serendipity p:2: no bubbles, the economic 8-point quad rule; Reference
    // Norm 1.0578593.
    check_model_pde_p(true, "sd_p2", 1281, 1.057859263392748);
}

/// `(slot, row)` pairs (1-based) of a `case.result` Perm block: the used
/// slots only, in slot order.
fn read_perm_pairs(path: &std::path::Path) -> Vec<(usize, u32)> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .map(|l| {
            let t: Vec<&str> = l.split_whitespace().collect();
            (t[0].parse().unwrap(), t[1].parse().unwrap())
        })
        .collect()
}

/// The quick tests `QuadraturesBrick`, `QuadraturesWedge` and their
/// serendipity `SD_` variants: ModelPDE with `Element = "p:2"` on a 2x2x2
/// mesh of `case` (`quadratures_brick` / `quadratures_wedge`), a flux BC on
/// 3D boundary faces, run as their RUN blocks do with the default rule and
/// the listed `Number of Integration Points`. Each run gives the pristine
/// Perm (`<prefix>perm.dat`, used slots), values (`<prefix>solution.dat`,
/// in that slot order) and norm.
fn check_quadratures(case: &str, serendipity: bool, nps: &[Option<i32>], n_dofs: usize, norm: f64) {
    let prefix = if serendipity { "sd_" } else { "" };
    let mesh = SimMesh::new(read_elmer_mesh(&data(case), 3));
    let want_perm = read_perm_pairs(&data(case).join(format!("{prefix}perm.dat")));
    let want = read_solution(&data(case).join(format!("{prefix}solution.dat")));
    for &np in nps {
        let mut solver = vec![
            ("Equation", st("ModelPDE")),
            ("Variable", st("Field")),
            ("Procedure", st("\"ModelPDE\" \"AdvDiffSolver\"")),
            ("Element", st("p:2")),
            ("Linear System Solver", st("Direct")),
            ("Steady State Convergence Tolerance", real(1e-9)),
        ];
        if let Some(n) = np {
            solver.push(("Number of Integration Points", ints(&[n])));
        }
        let mut simulation = vec![
            ("Coordinate System", st("Cartesian")),
            ("Simulation Type", st("steady")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[0])),
        ];
        if !serendipity {
            simulation.push(("Serendipity P Elements", Value::Logical(false)));
        }
        let model = Model {
            simulation: sec(&simulation),
            solvers: vec![sec(&solver)],
            equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
            bodies: vec![sec(&[
                ("Target Bodies", ints(&[1])),
                ("Equation", ints(&[1])),
                ("Material", ints(&[1])),
                ("Body Force", ints(&[1])),
            ])],
            materials: vec![sec(&[("diffusion coefficient", real(1.0))])],
            body_forces: vec![sec(&[("Field Source", real(0.0))])],
            boundary_conditions: vec![
                sec(&[("Target Boundaries", ints(&[1]))]),
                sec(&[("Target Boundaries", ints(&[2])), ("Field", real(0.0))]),
                sec(&[("Target Boundaries", ints(&[3])), ("Field Flux", real(1.0))]),
            ],
            ..Model::default()
        };
        let (mut sim, reg) = support::cases::simulation(model);
        sim.backends.dofs = Some(Box::new(elmer_pdofs::PDofsBackend::new()));
        sim.setup(
            &mesh,
            &reg,
            &elmer_core::NullLogger,
            &elmer_core::SerialComm,
        )
        .unwrap();
        sim.run(&mesh, &elmer_core::NullLogger, &elmer_core::SerialComm)
            .unwrap();
        let v = sim.vars.get("field").unwrap();
        assert_eq!(v.values.len(), n_dofs, "{np:?}");
        let perm = v.perm.as_deref().unwrap();
        let ours: Vec<(usize, u32)> = perm
            .iter()
            .enumerate()
            .filter(|(_, &p)| p != elmer_core::NO_DOF)
            .map(|(i, &p)| (i + 1, p + 1))
            .collect();
        assert_eq!(ours, want_perm, "{case} {np:?}: Perm");
        let got: Vec<f64> = want_perm
            .iter()
            .map(|&(slot, _)| v.values[perm[slot - 1] as usize])
            .collect();
        assert!(
            max_diff(&got, &want) < 1e-12,
            "{np:?}: {:e}",
            max_diff(&got, &want)
        );
        assert!(
            (v.conv.norm - norm).abs() < 1e-12,
            "{np:?}: norm {}",
            v.conv.norm
        );
    }
}

#[test]
fn quadratures_brick_matches_elmer() {
    // 27 nodes + 54 edges + 36 faces + 8 bubbles; the solution is linear.
    check_quadratures(
        "quadratures_brick",
        false,
        &[None, Some(27), Some(64)],
        125,
        0.3,
    );
}

#[test]
fn sd_quadratures_brick_matches_elmer() {
    // Serendipity p:2: 27 nodes + 54 edges; Reference Norm 0.37267800.
    check_quadratures(
        "quadratures_brick",
        true,
        &[None, Some(27), Some(64), Some(8)],
        81,
        0.37267799624996323,
    );
}

/// The point counts of the wedge cases' RUN blocks after the first.
const WEDGE_NPS: [i32; 15] = [64, 44, 28, 27, 24, 21, 18, 16, 15, 14, 12, 11, 10, 8, 44];

#[test]
fn quadratures_wedge_matches_elmer() {
    // 16 wedges: 27 nodes + 66 edges + 32 quad faces (the 24 triangle
    // faces carry no p:2 dofs but keep their Perm slots: 149 slots).
    let mut nps = vec![None];
    nps.extend(WEDGE_NPS.iter().map(|&n| Some(n)));
    check_quadratures("quadratures_wedge", false, &nps, 125, 0.29999999999999927);
}

#[test]
fn sd_quadratures_wedge_matches_elmer() {
    // Serendipity p:2: 27 nodes + 66 edges; Reference Norm 3.47804172E-01.
    let mut nps = vec![None];
    nps.extend(WEDGE_NPS[..14].iter().map(|&n| Some(n)));
    nps.extend([7, 6, 5, 44].iter().map(|&n| Some(n)));
    check_quadratures("quadratures_wedge", true, &nps, 93, 0.3478041718201255);
}

fn setup_error(model: Model, mesh: &SimMesh) -> String {
    let (mut sim, reg) = support::cases::simulation(model);
    sim.setup(mesh, &reg, &elmer_core::NullLogger, &elmer_core::SerialComm)
        .unwrap_err()
        .to_string()
}

#[test]
fn unsupported_advdiff_cases_are_refused_at_setup() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("adv_diff2"), 2));

    let mut m = adv_diff2_model(Method::Bubbles, true);
    m.equations[0].insert("Concentration Units", st("Mass To Max Solubility"));
    assert!(setup_error(m, &mesh).contains("Concentration Units"));

    let mut m = adv_diff2_model(Method::Bubbles, true);
    m.materials[0].insert("Compressibility Model", st("Perfect Gas 1"));
    assert!(setup_error(m, &mesh).contains("Compressibility Model"));

    let mut m = adv_diff2_model(Method::Bubbles, true);
    m.simulation
        .insert("Coordinate System", st("Axi Symmetric"));
    assert!(setup_error(m, &mesh).contains("Coordinate System"));

    let mut m = adv_diff2_model(Method::Bubbles, true);
    m.boundary_conditions[1].insert("Poison Solubility Change Boundary", Value::Logical(true));
    assert!(setup_error(m, &mesh).contains("Solubility Change Boundary"));

    // Elmer's own fatal errors.
    let mut m = adv_diff2_model(Method::Supg, true);
    m.solvers[0].insert("Bubbles", Value::Logical(true));
    assert!(setup_error(m, &mesh).contains("Choose either Bubbles or Stabilize"));

    let mut m = adv_diff2_model(Method::Bubbles, true);
    m.solvers[0].insert("Stabilization Method", st("streamline"));
    assert!(setup_error(m, &mesh).contains("Unknown stabilization method"));
}

#[test]
fn supg_on_quadrilaterals_is_refused_until_second_derivatives_exist() {
    // SU(p) needs ddBasisddx, which elmer-fe does not provide for 404.
    let mesh = SimMesh::new(read_elmer_mesh(&data("adv_diff2"), 2));
    let err = run(adv_diff2_model(Method::Supg, true), &mesh)
        .err()
        .expect("must be refused")
        .to_string();
    assert!(err.contains("second derivatives"), "{err}");
}

#[test]
fn unsupported_model_pde_cases_are_refused_at_setup() {
    let mesh = SimMesh::new(read_elmer_mesh(&data("model_pde"), 2));

    let mut m = model_pde_model(false);
    m.solvers[0].insert("Element", st("p:2"));
    assert!(setup_error(m, &mesh).contains("Element"));

    let mut m = model_pde_model(false);
    m.solvers[0].insert("CutFEM", Value::Logical(true));
    assert!(setup_error(m, &mesh).contains("CutFEM"));

    let mut m = model_pde_model(false);
    m.solvers[0].insert("Relative Integration Order", ints(&[1]));
    assert!(setup_error(m, &mesh).contains("Relative Integration Order"));
}
