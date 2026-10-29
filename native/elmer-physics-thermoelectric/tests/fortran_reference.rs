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
