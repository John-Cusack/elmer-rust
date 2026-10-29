//! `Apply Explicit Control`, checked against the pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) run of `fem/tests/HeatControlExplicit` from
//! its `case.sif` on its ElmerGrid mesh (`tests/data/expcontrol/`): the
//! temperature at the control node, kept in `cpar` before each assembly and
//! after each visit, drives the test's `AlternatingSource` (here a test-local
//! function), which flips sign whenever `cpar` leaves [-0.07, 0.07]. The
//! flips happen at the pristine run's times with its `cpar` values, Solver
//! 1's reference norm passes, and the final `|cpar|` is Solver 2's
//! (SaveScalars, left out) reference. HeatSolve is the test-local stand-in.

mod support;

use std::path::Path;
use std::sync::{Arc, Mutex};

use elmer_core::SerialComm;
use elmer_expr::Udf;
use elmer_solver::{Simulation, SolverRegistry};
use support::poisson::heat_nl;
use support::{sif_model, Capture};

/// The sign and the (time, cpar) of each flip.
static STATE: Mutex<(f64, Vec<(f64, f64)>)> = Mutex::new((1.0, Vec::new()));

#[test]
fn source_follows_the_control_value() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/expcontrol");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let mut model = sif_model(&dir.join("case.sif"));
    model.solvers.truncate(1);
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    // AlternatingSource.F90: f = sgn, sgn flipping at |tcenter| > 0.07.
    sim.backends.udfs.register_fortran(
        "AlternatingSource",
        "AlternatingSource",
        Udf::Node(Arc::new(|_node, x: &[f64]| {
            let (tcenter, time) = (x[0], x[1]);
            let mut st = STATE.lock().unwrap();
            if st.0 == 1.0 && tcenter > 0.07 {
                st.0 = -1.0;
                st.1.push((time, tcenter));
            }
            if st.0 == -1.0 && tcenter < -0.07 {
                st.0 = 1.0;
                st.1.push((time, tcenter));
            }
            st.0
        })),
    );
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));

    // Pristine: ` Neg source:  0.19...  7.062_197_848_033_637e-2` and so on.
    let want = [
        (0.19, 7.062_197_848_033_637e-2),
        (0.40, -7.017_077_757_793_04e-2),
        (0.61, 7.018_097_013_536_626e-2),
        (0.82, -7.018_073_987_155_621e-2),
    ];
    let got = STATE.lock().unwrap().1.clone();
    assert_eq!(got.len(), want.len(), "{got:?}");
    for ((t, c), (wt, wc)) in got.iter().zip(want) {
        assert!((t - wt).abs() < 1e-9, "{t} vs {wt}");
        assert!((c - wc).abs() < 1e-6 * wc.abs(), "{c} vs {wc}");
    }
    let cpar = sim.vars.get("cpar").unwrap().values[0];
    assert!((cpar.abs() - 6.86000303E-02).abs() < 1e-6 * 6.86000303E-02, "{cpar}");
}
