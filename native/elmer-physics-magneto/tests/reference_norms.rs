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
