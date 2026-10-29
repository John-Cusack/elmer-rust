//! `savedependence_` + `savedependence_init_` — port of SaveDependence.F90.
//!
//! Faithfulness notes: the entry is file I/O driven by Elmer/MATC calls,
//! so the wrappers own the file unit, the grid, and the norm accumulator
//! (same lifetimes as the source locals); Rust drives the parallel guard,
//! the row loop, and the finish. Formats (`I6,ES15.6`, `2ES15.6`, `ES15.6`)
//! are transcribed verbatim. Norm accumulates only the `Show Norm Index`
//! expression, `SQRT(sum/n)`, then sets the solver variable — exactly as
//! the source.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn sdd_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sdd_parallel_skip() -> c_int;
    fn sdd_info_start();
    fn sdd_init();
    fn sdd_setup(n: *mut c_int, nopar: *mut c_int, normind: *mut c_int, takeder: *mut c_int);
    fn sdd_write_row(i: c_int, n: c_int, normind: c_int);
    fn sdd_finish(n: c_int, normind: c_int);
}

/// `savedependence_init_` — Show Norm Index variable setup.
#[allow(non_snake_case)]
pub unsafe extern "C" fn savedependence_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        sdd_begin(model, solver);
        sdd_init();
    }
}

/// `savedependence_` — the table writer.
pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sdd_begin(model, solver) };
    if unsafe { sdd_parallel_skip() } != 0 {
        return;
    }
    unsafe { sdd_info_start() };

    let (mut n, mut nopar, mut normind, mut _takeder) = (0, 0, 0, 0);
    unsafe { sdd_setup(&mut n, &mut nopar, &mut normind, &mut _takeder) };
    if nopar == 0 {
        return;
    }
    for i in 1..=n {
        unsafe { sdd_write_row(i, n, normind) };
    }
    unsafe { sdd_finish(n, normind) };
}
