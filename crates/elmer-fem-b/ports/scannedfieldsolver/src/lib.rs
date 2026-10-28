//! Port of `fem/src/modules/ScannedFieldSolver.F90` (Elmer @ 9f6af2f85,
//! the campaign fork) for lane T021-B, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ScannedFieldSolver.F90` — two F77-global
//! entries, `scannedfieldsolver_init_` (solver-parameter setup: counts the
//! pre-existing `Exported Variable` entries, counts `Field Name` entries,
//! adds the `-nooutput SF_dummy` variable, `Number of Fields`, `Scan Max`,
//! and per-field `Scan i name` / `-nooutput name Dummy` exports) and
//! `scannedfieldsolver_` (copies each named field into its `scan i name`
//! variable, accumulates the `name Dummy` sum, and writes the sum back on
//! the last scan). The emitted dynamic symbols are exactly
//! `scannedfieldsolver_` and `scannedfieldsolver_init_` (proven by `nm -D`
//! on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/scannedfieldsolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive whole-array
//! updates (`ScanField % Values = Field % Values`, the sum accumulation,
//! the final write-back), `i2s`-composed names, and `Info`/`Fatal` text
//! stay in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ScannedFieldSolver.F90 statement order):
//! - `SAVE` semantics replicated: the `First` flag plus `FieldVars`,
//!   `FieldNames`, `NofFields`, and `Mesh` live as wrapper-side SAVE state
//!   with identical lifetimes (allocated once on the first call, reused
//!   across scan-loop and timestep re-entry, exactly as the original's
//!   internal-procedure SAVE).
//! - The `First` block's `N = Mesh % MaxElementDOFs` sizing touch, the
//!   `STAT=istat` allocation check (`Fatal('ScannedFieldSolver', ...)`),
//!   and the per-field `VariableGet` + `Fatal(... ' not associated!')`
//!   are transcribed into the wrapper.
//! - `Model`, `dt`, `TransientSimulation` are unused by the bodies, exactly
//!   as in the original (received and ignored; `Model` only supplies the
//!   mesh handle on the first call).
//! - The original's `SumFieldVar` association check tests `ScanFieldVar`
//!   rather than `SumFieldVar` (with the `' not associated'` message
//!   lacking the `!`); replicated exactly (bug-compatible).
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn sfa_init_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sfa_init_count_exported() -> c_int;
    fn sfa_init_count_fields() -> c_int;
    fn sfa_init_add_variable();
    fn sfa_init_add_nfields(n: c_int);
    fn sfa_init_scan_solver() -> c_int;
    fn sfa_init_scanmax(solverint: c_int) -> c_int;
    fn sfa_init_field_exports(fieldint: c_int, scanmax: c_int, i: *mut c_int);
    fn sfa_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sfa_first_block();
    fn sfa_info_header();
    fn sfa_scan_index() -> c_int;
    fn sfa_scanmax() -> c_int;
    fn sfa_nfields() -> c_int;
    fn sfa_field_select(fi: c_int);
    fn sfa_field_scan(fi: c_int, scanint: c_int);
    fn sfa_field_sum(fi: c_int);
    fn sfa_field_combine(scanint: c_int, scanmax: c_int);
}

/// Entry-sharing guard only; all SAVE state lives in the Fortran wrappers.
static STATE: Mutex<()> = Mutex::new(());

/// `scannedfieldsolver_init_` — the solver init entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient flag
/// are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "scannedfieldsolver_init_"]
pub unsafe extern "C" fn scanned_field_solver_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe { sfa_init_begin(model, solver) };

    // Count pre-existing "Exported Variable " entries (first free index).
    let mut i = unsafe { sfa_init_count_exported() };
    // Count "Field Name " entries (NofFields).
    let nof_fields = unsafe { sfa_init_count_fields() };

    unsafe {
        sfa_init_add_variable();
        sfa_init_add_nfields(nof_fields);
    }

    let scan_solver = unsafe { sfa_init_scan_solver() };
    let scan_max = unsafe { sfa_init_scanmax(scan_solver) };

    for field in 1..=nof_fields {
        unsafe { sfa_init_field_exports(field, scan_max, &mut i) };
    }
}

/// `scannedfieldsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient flag
/// are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "scannedfieldsolver_"]
pub unsafe extern "C" fn scanned_field_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe { sfa_begin(model, solver) };

    // IF (First) block: mesh sizing touch + one-time allocation.
    unsafe { sfa_first_block() };

    unsafe { sfa_info_header() };

    let scan_int = unsafe { sfa_scan_index() };
    let scan_max = unsafe { sfa_scanmax() };
    let nof_fields = unsafe { sfa_nfields() };

    for field in 1..=nof_fields {
        unsafe {
            sfa_field_select(field);
            sfa_field_scan(field, scan_int);
            sfa_field_sum(field);
            sfa_field_combine(scan_int, scan_max);
        }
    }
}
