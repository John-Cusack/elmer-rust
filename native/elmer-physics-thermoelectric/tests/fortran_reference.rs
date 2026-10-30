//! `ThermoElectricSolver` against the pristine ElmerSolver on the
//! `fem/tests/ThermoElectric*` cases in `tests/data` (see
//! `tests/data/README.md`): every nodal temperature and potential of the last
//! time level of its `.result`, and the norm Elmer prints.

mod support;

use elmer_solver::Simulation;
use support::cases::{array, data, ints, logical, material, model, nodal, norm, read_solution, real, run, sec, st, steady, te_solver};
use support::mesh::TestMesh;

fn compare(case: &str, sim: &Simulation, mesh: &TestMesh, ref_norm: Option<f64>, norm_tol: f64, tol: f64) {
    let got = nodal(sim, mesh.coords.len());
    let want = read_solution(&data(case).join("solution.dat"));
    assert_eq!(got.len(), want.len());
    let mut err = [0.0f64; 2];
    let mut scale = [0.0f64; 2];
    for (g, w) in got.iter().zip(&want) {
        for c in 0..2 {
            err[c] = err[c].max((g[c] - w[c]).abs());
            scale[c] = scale[c].max(w[c].abs());
        }
    }
    eprintln!(
        "{case}: max relative error [{:.1e}, {:.1e}]",
        err[0] / scale[0].max(1e-300),
        err[1] / scale[1].max(1e-300)
    );
    // Elmer's norm counts the bubble dofs it keeps in the global system, so
    // the bubble cases compare the nodal values only.
    if let Some(ref_norm) = ref_norm {
        let rel_norm = (norm(sim) - ref_norm).abs() / ref_norm;
        eprintln!("{case}: norm {:.17e} (ref {ref_norm:.17e}, rel {rel_norm:.1e})", norm(sim));
        assert!(rel_norm < norm_tol, "{case}: norm {} vs {ref_norm}", norm(sim));
    }
    for c in 0..2 {
        assert!(err[c] <= tol * scale[c].max(1.0e-300), "{case}: component {c} max error {} (scale {})", err[c], scale[c]);
    }
}

fn body_force() -> elmer_core::ParamMap {
    sec(&[("Current Source", real(100.0))])
}

fn dirichlet(tags: &[i32]) -> elmer_core::ParamMap {
    sec(&[("Target Boundaries", ints(tags)), ("Potential", real(0.0)), ("Temperature", real(300.0))])
}

/// The Robin and current-density conditions of the `AnisoBC`, `Bubble` and
/// `Stabilize` cases.
fn aniso_bcs() -> Vec<elmer_core::ParamMap> {
    vec![
        dirichlet(&[1]),
        sec(&[
            ("Target Boundaries", ints(&[2, 4])),
            ("Electric Transfer Coefficient", array(2, 2, &[2.0, 0.5, 0.5, 1.0])),
            ("External Potential", real(5.0)),
        ]),
        sec(&[
            ("Target Boundaries", ints(&[3])),
            ("Electric Flux Vector 1", real(30.0)),
            ("Electric Flux Vector 2", real(50.0)),
        ]),
    ]
}

/// `ThermoElectric/ptest.sif`: a Seebeck material with a current source on a
/// square held at 300 K and 0 V; Newton from the second iteration.
#[test]
fn plain_matches_elmer() {
    let mesh = TestMesh::read(&data("plain"), 2);
    let solver = te_solver(&[("Element", st("p:1"))]);
    let sim = run(model(steady(), solver, material(), body_force(), vec![dirichlet(&[1])]), &mesh).unwrap();
    compare("plain", &sim, &mesh, Some(221.28739379604136), 1e-11, 1e-10);
}

/// `ThermoElectricAnisoBC`: the anisotropic Robin condition on the sides and
/// the current density vector on the top.
#[test]
fn aniso_bc_matches_elmer() {
    let mesh = TestMesh::read(&data("aniso_bc"), 2);
    let sim = run(model(steady(), te_solver(&[]), material(), body_force(), aniso_bcs()), &mesh).unwrap();
    compare("aniso_bc", &sim, &mesh, Some(1140.1563699081962), 1e-11, 1e-10);
}

/// `ThermoElectricStabilize`: SUPG stabilization of the Peltier drift.
#[test]
fn stabilize_matches_elmer() {
    let mesh = TestMesh::read(&data("stabilize"), 2);
    let solver = te_solver(&[("Stabilize", logical(true))]);
    let sim = run(model(steady(), solver, material(), body_force(), aniso_bcs()), &mesh).unwrap();
    compare("stabilize", &sim, &mesh, Some(1140.4826593463742), 1e-11, 1e-10);
}


/// The nonlinear loop settings of the bubble cases: converged to round-off
/// so that the global (Elmer) and condensed (here) bubbles agree.
fn tight() -> [(&'static str, elmer_core::Value); 2] {
    [("Nonlinear System Convergence Tolerance", real(1.0e-12)), ("Nonlinear System Max Iterations", ints(&[100]))]
}

/// `ThermoElectricBubble`: `Element = "p:1 b:4"`, four interior bubble modes
/// per quadrilateral, which Elmer keeps in the global system and this crate
/// condenses element by element.
#[test]
fn bubble_matches_elmer() {
    let mesh = TestMesh::read(&data("bubble"), 2);
    let mut extra = vec![("Element", st("p:1 b:4"))];
    extra.extend(tight());
    let sim = run(model(steady(), te_solver(&extra), material(), body_force(), aniso_bcs()), &mesh).unwrap();
    compare("bubble", &sim, &mesh, None, 0.0, 1e-8);
}

/// `ThermoElectricBubbleTransient`: the same element over five backward-Euler
/// steps, which carry the bubble values from step to step.
#[test]
fn bubble_transient_matches_elmer() {
    let mesh = TestMesh::read(&data("bubble_transient"), 2);
    let mut extra = vec![("Element", st("p:1 b:4"))];
    extra.extend(tight());
    let sim_sec = support::cases::simulation(&[
        ("Simulation Type", st("Transient")),
        ("Timestepping Method", st("BDF")),
        ("BDF Order", ints(&[1])),
        ("Timestep Sizes", real(0.1)),
        ("Timestep Intervals", ints(&[5])),
    ]);
    let sim = run(model(sim_sec, te_solver(&extra), material(), body_force(), aniso_bcs()), &mesh).unwrap();
    compare("bubble_transient", &sim, &mesh, None, 0.0, 1e-8);
}

/// The whitespace-split numbers of a fixture, `#` comments skipped.
fn rows(path: &std::path::Path) -> Vec<Vec<f64>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|t| t.parse().unwrap()).collect())
        .collect()
}

/// A bubble case with the p dofs from `elmer-pdofs`: the bubbles in the
/// global system, as in Elmer, so every slot of the `case.result`
/// (`perm.dat`, `slots.dat`: 121 nodes + 100 x 4 bubbles) and the norm
/// Elmer prints, which counts the bubbles, are compared.
fn global_bubble_case(case: &str, simulation: elmer_core::ParamMap, ref_norm: f64) {
    use elmer_core::{NullLogger, SerialComm};
    let mesh = TestMesh::read(&data(case), 2);
    let mut extra = vec![("Element", st("p:1 b:4"))];
    extra.extend(tight());
    let mut reg = elmer_solver::SolverRegistry::new();
    elmer_physics_thermoelectric::register(&mut reg);
    let mut sim = Simulation::new(model(simulation, te_solver(&extra), material(), body_force(), aniso_bcs()));
    sim.backends.dofs = Some(Box::new(elmer_pdofs::PDofsBackend::new()));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();

    let v = sim.vars.get("pot").unwrap();
    assert_eq!(v.values.len(), 2 * 521, "{case}");
    let perm = v.perm.as_deref().unwrap();
    let want_perm = rows(&data(case).join("perm.dat"));
    let want = rows(&data(case).join("slots.dat"));
    assert_eq!(want_perm.len(), 521, "{case}");
    let mut scale = [0.0f64; 2];
    for w in &want {
        for c in 0..2 {
            scale[c] = scale[c].max(w[1 + c].abs());
        }
    }
    let mut err = [0.0f64; 2];
    for (p, w) in want_perm.iter().zip(&want) {
        let slot = p[0] as usize;
        assert_eq!(f64::from(perm[slot - 1] + 1), p[1], "{case}: Perm of slot {slot}");
        let row = perm[slot - 1] as usize;
        for c in 0..2 {
            err[c] = err[c].max((v.values[2 * row + c] - w[1 + c]).abs());
        }
    }
    let norm = norm(&sim);
    let rel = [err[0] / scale[0], err[1] / scale[1]];
    eprintln!("{case} (global): norm {norm:.17e}, max relative error [{:.1e}, {:.1e}]", rel[0], rel[1]);
    assert!((norm - ref_norm).abs() / ref_norm < 1e-11, "{case}: norm {norm} vs {ref_norm}");
    for c in 0..2 {
        assert!(err[c] <= 1e-10 * scale[c], "{case}: component {c} max error {} (scale {})", err[c], scale[c]);
    }
}

#[test]
fn bubble_in_global_system_matches_elmer() {
    global_bubble_case("bubble", steady(), 549.0240106292581);
}

/// The transient case: the bubbles are global dofs, so `Default1stOrderTime`
/// integrates them as Elmer does; the last of the five steps is compared.
#[test]
fn bubble_transient_in_global_system_matches_elmer() {
    let sim_sec = support::cases::simulation(&[
        ("Simulation Type", st("Transient")),
        ("Timestepping Method", st("BDF")),
        ("BDF Order", ints(&[1])),
        ("Timestep Sizes", real(0.1)),
        ("Timestep Intervals", ints(&[5])),
    ]);
    global_bubble_case("bubble_transient", sim_sec, 343.1731238726303);
}
