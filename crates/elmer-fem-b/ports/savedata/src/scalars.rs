// SaveScalars pilot: Rust driver calling Fortran accessor.
use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn sds_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sds_savescalars();
}

pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sds_begin(model, solver) };
    unsafe { sds_savescalars() };
}

#[allow(non_snake_case)]
pub unsafe extern "C" fn savescalars_init_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    // SaveScalars_init is a no-op in source; nothing to do.
}
