// Port of `fem/src/umf4_f77wrapper.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/umf4_f77wrapper.c. True object symbols: `umf4def`,
// `umf4sym`, `umf4num`, `umf4solr`, `umf4sol`, `umf4fnum`, `umf4fsym` plus
// the `SuiteSparse_long` variants `umf4_l_def`, `umf4_l_sym`,
// `umf4_l_num`, `umf4_l_solr`, `umf4_l_sol`, `umf4_l_fnum`, `umf4_l_fsym`
// (`-DUSE_ISO_C_BINDINGS` plain names; `HAVE_UMFPACK` is defined).
//
// Thin pass-throughs to the UMFPACK `di` (int) / `dl` (long) entry
// points, declared extern and resolved at the final ElmerSolver link
// exactly like the C object. Statuses are discarded (`(void)` casts)
// and `Control[7]` (`UMFPACK_IRSTEP`) is zeroed in the `_sol` variants,
// matching the original statement by statement.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int, c_long};

extern "C" {
    fn umfpack_di_defaults(control: *mut c_double);
    fn umfpack_di_symbolic(
        n_row: c_int, n_col: c_int,
        ap: *const c_int, ai: *const c_int, ax: *const c_double,
        symbolic: *mut *mut c_void, control: *const c_double, info: *mut c_double,
    ) -> c_int;
    fn umfpack_di_numeric(
        ap: *const c_int, ai: *const c_int, ax: *const c_double,
        symbolic: *mut c_void, numeric: *mut *mut c_void,
        control: *const c_double, info: *mut c_double,
    ) -> c_int;
    fn umfpack_di_solve(
        sys: c_int, ap: *const c_int, ai: *const c_int, ax: *const c_double,
        x: *mut c_double, b: *const c_double, numeric: *mut c_void,
        control: *const c_double, info: *mut c_double,
    ) -> c_int;
    fn umfpack_di_free_numeric(numeric: *mut *mut c_void);
    fn umfpack_di_free_symbolic(symbolic: *mut *mut c_void);
    fn umfpack_dl_defaults(control: *mut c_double);
    fn umfpack_dl_symbolic(
        n_row: c_long, n_col: c_long,
        ap: *const c_long, ai: *const c_long, ax: *const c_double,
        symbolic: *mut *mut c_void, control: *const c_double, info: *mut c_double,
    ) -> c_long;
    fn umfpack_dl_numeric(
        ap: *const c_long, ai: *const c_long, ax: *const c_double,
        symbolic: *mut c_void, numeric: *mut *mut c_void,
        control: *const c_double, info: *mut c_double,
    ) -> c_long;
    fn umfpack_dl_solve(
        sys: c_long, ap: *const c_long, ai: *const c_long, ax: *const c_double,
        x: *mut c_double, b: *const c_double, numeric: *mut c_void,
        control: *const c_double, info: *mut c_double,
    ) -> c_long;
    fn umfpack_dl_free_numeric(numeric: *mut *mut c_void);
    fn umfpack_dl_free_symbolic(symbolic: *mut *mut c_void);
}

const UMFPACK_IRSTEP: usize = 7;

#[export_name = "umf4def"]
pub unsafe extern "C" fn port_umf4def(control: *mut c_double) {
    umfpack_di_defaults(control);
}

#[export_name = "umf4sym"]
pub unsafe extern "C" fn port_umf4sym(
    m: *const c_int, n: *const c_int,
    ap: *const c_int, ai: *const c_int, ax: *const c_double,
    symbolic: *mut *mut c_void, control: *const c_double, info: *mut c_double,
) {
    umfpack_di_symbolic(*m, *n, ap, ai, ax, symbolic, control, info);
}

#[export_name = "umf4num"]
pub unsafe extern "C" fn port_umf4num(
    ap: *const c_int, ai: *const c_int, ax: *const c_double,
    symbolic: *const *mut c_void, numeric: *mut *mut c_void,
    control: *const c_double, info: *mut c_double,
) {
    umfpack_di_numeric(ap, ai, ax, *symbolic, numeric, control, info);
}

#[export_name = "umf4solr"]
pub unsafe extern "C" fn port_umf4solr(
    sys: *const c_int, ap: *const c_int, ai: *const c_int, ax: *const c_double,
    x: *mut c_double, b: *const c_double, numeric: *mut *mut c_void,
    control: *const c_double, info: *mut c_double,
) {
    umfpack_di_solve(*sys, ap, ai, ax, x, b, *numeric, control, info);
}

#[export_name = "umf4sol"]
pub unsafe extern "C" fn port_umf4sol(
    sys: *const c_int, x: *mut c_double, b: *const c_double,
    numeric: *mut *mut c_void, control: *mut c_double, info: *mut c_double,
) {
    *control.add(UMFPACK_IRSTEP) = 0.0;
    umfpack_di_solve(
        *sys,
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        x,
        b,
        *numeric,
        control,
        info,
    );
}

#[export_name = "umf4fnum"]
pub unsafe extern "C" fn port_umf4fnum(numeric: *mut *mut c_void) {
    umfpack_di_free_numeric(numeric);
}

#[export_name = "umf4fsym"]
pub unsafe extern "C" fn port_umf4fsym(symbolic: *mut *mut c_void) {
    umfpack_di_free_symbolic(symbolic);
}

#[export_name = "umf4_l_def"]
pub unsafe extern "C" fn port_umf4_l_def(control: *mut c_double) {
    umfpack_dl_defaults(control);
}

#[export_name = "umf4_l_sym"]
pub unsafe extern "C" fn port_umf4_l_sym(
    m: *const c_long, n: *const c_long,
    ap: *const c_long, ai: *const c_long, ax: *const c_double,
    symbolic: *mut *mut c_void, control: *const c_double, info: *mut c_double,
) {
    umfpack_dl_symbolic(*m, *n, ap, ai, ax, symbolic, control, info);
}

#[export_name = "umf4_l_num"]
pub unsafe extern "C" fn port_umf4_l_num(
    ap: *const c_long, ai: *const c_long, ax: *const c_double,
    symbolic: *const *mut c_void, numeric: *mut *mut c_void,
    control: *const c_double, info: *mut c_double,
) {
    umfpack_dl_numeric(ap, ai, ax, *symbolic, numeric, control, info);
}

#[export_name = "umf4_l_solr"]
pub unsafe extern "C" fn port_umf4_l_solr(
    sys: *const c_long, ap: *const c_long, ai: *const c_long, ax: *const c_double,
    x: *mut c_double, b: *const c_double, numeric: *mut *mut c_void,
    control: *const c_double, info: *mut c_double,
) {
    umfpack_dl_solve(*sys, ap, ai, ax, x, b, *numeric, control, info);
}

#[export_name = "umf4_l_sol"]
pub unsafe extern "C" fn port_umf4_l_sol(
    sys: *const c_long, x: *mut c_double, b: *const c_double,
    numeric: *mut *mut c_void, control: *mut c_double, info: *mut c_double,
) {
    *control.add(UMFPACK_IRSTEP) = 0.0;
    umfpack_dl_solve(
        *sys,
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        x,
        b,
        *numeric,
        control,
        info,
    );
}

#[export_name = "umf4_l_fnum"]
pub unsafe extern "C" fn port_umf4_l_fnum(numeric: *mut *mut c_void) {
    umfpack_dl_free_numeric(numeric);
}

#[export_name = "umf4_l_fsym"]
pub unsafe extern "C" fn port_umf4_l_fsym(symbolic: *mut *mut c_void) {
    umfpack_dl_free_symbolic(symbolic);
}
