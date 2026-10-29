//! Steady `WhitneyAVSolver` against the pristine ElmerSolver reference
//! norms of the quick tests copied under `tests/data` (see the README there).

use std::path::PathBuf;
use std::sync::Mutex;

use elmer_physics_magneto::{solve_steady, Case};

/// The `.sif` `INCLUDE` paths resolve against the working directory, as
/// ElmerSolver runs in the case directory; tests change it one at a time.
static CWD: Mutex<()> = Mutex::new(());

fn load(case: &str, sif: &str) -> Case {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(case);
    let _guard = CWD.lock().unwrap_or_else(|e| e.into_inner());
    std::env::set_current_dir(&dir).unwrap();
    Case::load(&dir.join(sif)).unwrap()
}

/// `CompareToReferenceSolution`'s default relative tolerance.
const REL_TOL: f64 = 1e-5;

#[test]
fn mgdyn_airgap_solver1_norm() {
    let case = load("mgdyn_airgap", "airgap.sif");
    let sol = solve_steady(&case, 1).unwrap();
    // The pristine run prints Norm = 7.08538973E-05; the test's stored
    // reference is 7.08538965E-05.
    let reference = 7.08538965e-05;
    let rel = (sol.norm - reference).abs() / reference;
    eprintln!(
        "mgdyn_airgap: {} nodal + {} edge DOFs, norm {:.9e}, rel err {rel:.3e}",
        sol.n_nodal, sol.n_edges, sol.norm
    );
    assert!(rel < REL_TOL, "norm {} vs {reference}", sol.norm);
}

#[test]
fn unported_options_are_refused() {
    use elmer_core::Value;
    use elmer_physics_magneto::MagnetoError;
    for (key, value) in [
        ("Use Piola Transform", Value::Logical(true)),
        ("Use Lagrange Gauge", Value::Logical(true)),
        ("Linear System Solver", Value::Str("Direct".into())),
        ("Fix Input Current Density", Value::Logical(true)),
    ] {
        let mut case = load("mgdyn_airgap", "airgap.sif");
        case.solvers[0].as_mut().unwrap().insert(key, value);
        match solve_steady(&case, 1) {
            Err(MagnetoError::Unsupported(_)) => {}
            Err(e) => panic!("{key}: unexpected error {e}"),
            Ok(_) => panic!("{key}: accepted"),
        }
    }
}

#[test]
fn mgdyn_airgap_harmonic_solver1_converged_norm() {
    // The case's 1e-6 linear tolerance stops BiCGStab(l) inside a stretch
    // where the residual history reacts to rounding: a 1e-15 relative
    // perturbation of the matrix changes it from iteration 15 on. The norm
    // at that tolerance (pristine 2.42423001E-03, here 2.4242639e-3) then
    // depends on the last bits of the assembly. Converged to 1e-10 both
    // runs take 45 iterations with the same history up to iteration 14,
    // and the pristine run prints Norm = 2.42424125E-03.
    let mut case = load("mgdyn_airgap_harmonic", "airgap.sif");
    case.solvers[0]
        .as_mut()
        .unwrap()
        .insert("Linear System Convergence Tolerance", elmer_core::Value::real(1e-10));
    let sol = elmer_physics_magneto::solve_harmonic(&case, 1).unwrap();
    let reference = 2.42424125e-03;
    let rel = (sol.norm - reference).abs() / reference;
    eprintln!(
        "mgdyn_airgap_harmonic: {} nodal + {} edge complex DOFs, omega {}, {} iterations to {:.3e}, norm {:.9e}, rel err {rel:.3e}",
        sol.n_nodal, sol.n_edges, sol.omega, sol.stats.iterations, sol.stats.residual, sol.norm
    );
    assert_eq!(sol.stats.iterations, 45);
    assert!(rel < REL_TOL, "norm {} vs {reference}", sol.norm);
}

#[test]
fn unported_harmonic_options_are_refused() {
    use elmer_core::Value;
    use elmer_physics_magneto::MagnetoError;
    for (key, value) in [
        ("Linear System Solver", Value::Str("Direct".into())),
        ("Use Piola Transform", Value::Logical(true)),
        ("Variable", Value::Str("P".into())),
    ] {
        let mut case = load("mgdyn_airgap_harmonic", "airgap.sif");
        case.solvers[0].as_mut().unwrap().insert(key, value);
        match elmer_physics_magneto::solve_harmonic(&case, 1) {
            Err(MagnetoError::Unsupported(_) | MagnetoError::Input(_)) => {}
            Err(e) => panic!("{key}: unexpected error {e}"),
            Ok(_) => panic!("{key}: accepted"),
        }
    }
    let mut case = load("mgdyn_airgap_harmonic", "airgap.sif");
    case.materials[2].as_mut().unwrap().insert("Electric Conductivity im", Value::real(1.0));
    assert!(matches!(
        elmer_physics_magneto::solve_harmonic(&case, 1),
        Err(MagnetoError::Unsupported(_))
    ));
}
