// SaveLine pilot: Rust driver calling Fortran phase accessors.
use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn sdl_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sdl_init();
    fn sdl_setup(novar: *mut c_int, noresults: *mut c_int, maxboundary: *mut c_int);
    fn sdl_alloc(nresults: *mut c_int);
    fn sdl_saveline();
}

pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sdl_begin(model, solver) };
    unsafe { sdl_init() };
    let mut novar: c_int = 0;
    let mut noresults: c_int = 0;
    let mut maxboundary: c_int = 0;
    unsafe { sdl_setup(&mut novar, &mut noresults, &mut maxboundary) };
    if novar == 0 || noresults == 0 {
        return;
    }
    let mut nresults: c_int = 0;
    unsafe { sdl_alloc(&mut nresults) };
    if nresults < 0 {
        return;
    }
    unsafe { sdl_saveline() };
}

#[allow(non_snake_case)]
pub unsafe extern "C" fn saveline_init_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    // SaveLine_init is a no-op in source (just returns); nothing to do.
}
