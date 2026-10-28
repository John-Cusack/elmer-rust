//! Port of `fem/src/modules/MeshChecksum.F90` (Elmer @ fork) for lane
//! T021-B, at the solver-plugin boundary (whole `.so` → Rust cdylib with
//! the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/MeshChecksum.F90` — two F77-global
//! entries `meshchecksum_init_` + `meshchecksum_` (proven by `nm -D` on
//! the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/meshchecksum_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset.
//! Order-sensitive reductions (`SUM`) and list-directed `PRINT` stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors MeshChecksum.F90 statement for statement):
//! - `MeshChecksum_init` only calls `ListAddNewString` for the Variable
//!   spec; the body only reads mesh/element components, prints, and —
//!   when `Reference Values` is present — prints the reference, forms
//!   the pseudo-norm loop, assigns `Solver % Variable % Values`, and
//!   adds `Reference Norm`. All of that lives in the wrappers; Rust only
//!   sequences the `Found` branch.
//! - The original keeps no SAVE state (all locals), so the port keeps
//!   none either; the `Mutex` is only a `static`-sharing guard for the
//!   wrapper's saved solver handle, never contended across threads.
//! - `Model`, `dt`, `TransientSimulation` are unused by both entries,
//!   exactly as in the original (received and ignored).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn mca_begin(h_solver: *const c_void);
    fn mca_init_var();
    fn mca_info();
    fn mca_checksums(cs: *mut c_double);
    fn mca_print_this(cs: *const c_double);
    fn mca_get_ref(ref_vals: *mut c_double) -> c_int;
    fn mca_print_ref(ref_vals: *const c_double);
    fn mca_pseudonorm(cs: *const c_double, ref_vals: *const c_double) -> c_double;
    fn mca_print_pseudo(pseudo: c_double);
    fn mca_set_pseudo(pseudo: c_double);
}

static GUARD: Mutex<()> = Mutex::new(());

/// `meshchecksum_init_` — the solver init entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "meshchecksum_init_"]
pub unsafe extern "C" fn meshchecksum_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    // CALL ListAddNewString( Solver % Values,'Variable',
    //     '-nooutput -global meshcheck_var')
    unsafe { mca_begin(solver) };
    unsafe { mca_init_var() };
}

/// `meshchecksum_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "meshchecksum_"]
pub unsafe extern "C" fn meshchecksum_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe { mca_begin(solver) };

    // CALL Info(Caller,'Checking for mesh consistency')
    unsafe { mca_info() };

    // CheckSum accumulation over bulk + boundary elements and nodes.
    let mut cs = [0.0 as c_double; 8];
    unsafe { mca_checksums(cs.as_mut_ptr()) };

    // PRINT *,'Checksums for file output:'
    // PRINT *,'ThisResults:',NINT(CheckSum(1:7)),CheckSum(8)
    unsafe { mca_print_this(cs.as_ptr()) };

    // WrkArray => ListGetConstRealArray( Params,'Reference Values',Found )
    let mut ref_vals = [0.0 as c_double; 8];
    let found = unsafe { mca_get_ref(ref_vals.as_mut_ptr()) };
    if found != 0 {
        // PRINT *,'RefResults:',NINT(RefSum(1:7)),RefSum(8)
        unsafe { mca_print_ref(ref_vals.as_ptr()) };
        // PseudoNorm loop over the 8 checksum slots.
        let pseudo = unsafe { mca_pseudonorm(cs.as_ptr(), ref_vals.as_ptr()) };
        // PRINT *,'PseudoNorm:',PseudoNorm
        unsafe { mca_print_pseudo(pseudo) };
        // Solver % Variable % Values = PseudoNorm (+ Reference Norm entry).
        unsafe { mca_set_pseudo(pseudo) };
    }
}
