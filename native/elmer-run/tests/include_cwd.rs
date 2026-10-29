//! heateq_spline_table from its pristine `.sif`, whose `include cond.dat`
//! resolves against the working directory as in Elmer (ElmerSolver runs in
//! the case directory). This binary holds one test, so changing the
//! process's working directory affects nothing else.

mod support;

use elmer_core::SerialComm;
use elmer_run::{run_case, Case, Frontend};
use support::case::{scratch_case, sif_frontend};
use support::logger;

#[test]
fn heateq_spline_table() {
    let dir = scratch_case("heateq_spline_table", "sif-spline");
    std::env::set_current_dir(&dir).unwrap();
    let frontend = sif_frontend();
    let loaded = frontend.load(&dir.join("TempDist.sif")).unwrap();
    let (log, buf) = logger();
    let case = Case {
        model: loaded.model,
        reloads: loaded.reloads,
        mesh: loaded.mesh.as_ref(),
        case_dir: loaded.case_dir,
        mesh_dir: loaded.mesh_dir,
    };
    let out = run_case(case, &log, &SerialComm).unwrap_or_else(|e| panic!("{e}\n{}", buf.text()));
    let Some(&(_, norm)) = out.norms.iter().find(|(id, _)| *id == 1) else {
        panic!("{:?}", out.norms)
    };
    let reference = 1.5131;
    let rel = (norm - reference).abs() / reference;
    println!("heateq_spline_table: native norm {norm:.17e}, reference {reference:e}, rel {rel:e}");
    // `Reference Norm Tolerance` default.
    assert!(rel <= 1.0e-5, "norm {norm}, rel {rel:e}");
    assert_eq!(out.passed, Some(true));
}
