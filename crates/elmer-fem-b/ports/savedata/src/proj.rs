//! `saveprojection_` + `saveprojection_init_` — port of SaveProjection.F90.
//!
//! Faithfulness notes: the entry is projector-matrix traversal driven by
//! Elmer calls, so each variable's projection (GenerateProjectors output,
//! MortarBC/BC PMatrix walk, ToSlave/ToMaster/Normalize arms, WHERE
//! normalization) runs in one Fortran wrapper off captured handles; Rust
//! drives the variable scan, accumulates the combined L2 norm exactly as
//! the source (`Nrm = SUM(target**2)` per variable, single SQRT at the
//! end), and applies the size-1 norm write. Formats (`ES12.3`) are
//! transcribed verbatim.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn sdp_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sdp_init();
    fn sdp_info_start();
    fn sdp_gen_projectors();
    fn sdp_count_vars() -> c_int;
    fn sdp_project_var(idx: c_int) -> c_double;
    fn sdp_set_norm(nrm: c_double);
}

/// `saveprojection_init_` — Variable setup.
#[allow(non_snake_case)]
pub unsafe extern "C" fn saveprojection_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        sdp_begin(model, solver);
        sdp_init();
    }
}

/// `saveprojection_` — the projector.
pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sdp_begin(model, solver) };
    unsafe { sdp_info_start() };
    unsafe { sdp_gen_projectors() };

    let no_var = unsafe { sdp_count_vars() };
    if no_var == 0 {
        return;
    }
    let mut nrm: f64 = 0.0;
    for i in 1..=no_var {
        nrm += unsafe { sdp_project_var(i) };
    }
    let nrm = nrm.sqrt();
    unsafe { sdp_set_norm(nrm) };
}
