//! Port of `fem/src/modules/ElectricForce.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/ElectricForce.F90` — one F77-global
//! entry `statelecforce_` with one `CONTAINS` internal. The emitted
//! dynamic symbol is exactly `statelecforce_` (proven by `nm -D` on the
//! reference `.so`); the internal emits only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/electricforce_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. Order-
//! sensitive reductions (`SUM`, `MATMUL`) and all formatted/list-directed
//! I/O (`WRITE` to `Message`, `Info`) stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ElectricForce.F90 statement for statement):
//! - The `Potential Field Name` lookup replicates the original's `stat`
//!   logic exactly, including the quirk that a custom (found) name skips
//!   the `VariableGet` lookup; the default path resolves `'Potential'`
//!   and fatals identically when absent.
//! - The `FirstTime` local the original sets-but-never-reads is not
//!   replicated (no observable state).
//! - Entry-level flags (`CalculateNodal/Field/Moment`, displaced mode)
//!   and the entry locals (nodes, permittivity, accumulators) live as
//!   wrapper-side SAVE state with identical per-invocation lifetimes
//!   (allocated at entry, consumed by the sweep, deallocated in the
//!   teardown, as in the original).
//! - The `ElementCode == 101` skip and the active-boundary `CYCLE` are
//!   driven from Rust on wrapper-returned flags, matching the original's
//!   skip order (code check before the boundary test).
//! - The nodal-branch `InvalidateVariable` call names
//!   `'Electric Force Density'` (not `'Electric Nodal Force'`), exactly
//!   as in the original.
//! - `Model` is used by the body (current element, BCs, bodies,
//!   materials, constants, simulation results); `dt` and the transient
//!   flag are unused, exactly as in the original (received and ignored).
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn efe_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn efe_pot_setup();
    fn efe_nodal_setup();
    fn efe_field_setup();
    fn efe_alloc();
    fn efe_displaced_setup();
    fn efe_const_setup();
    fn efe_zero();
    fn efe_counts(nbulk: *mut c_int, nbnd: *mut c_int);
    fn efe_elem_begin(t: c_int, n: *mut c_int, is101: *mut c_int);
    fn efe_elem_active() -> c_int;
    fn efe_elem_compute();
    fn efe_normalize();
    fn efe_report();
    fn efe_finish();
}

// No first-mesh SAVE arrays: the original allocates and frees its work
// arrays on every invocation, so the port keeps no cross-call state. The
// mutex only guards `static` sharing, never contended across threads.
static GUARD: Mutex<()> = Mutex::new(());

/// `statelecforce_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "statelecforce_"]
pub unsafe extern "C" fn statelecforce_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe {
        efe_begin(model, solver);
        efe_pot_setup();
        efe_nodal_setup();
        efe_field_setup();
        efe_alloc();
        efe_displaced_setup();
        efe_const_setup();
        efe_zero();

        // Boundary sweep: bulk+1 ..= bulk+boundary, as in the original.
        let (mut nbulk, mut nbnd) = (0, 0);
        efe_counts(&mut nbulk, &mut nbnd);
        for t in (nbulk + 1)..=(nbulk + nbnd) {
            let (mut n, mut is101) = (0, 0);
            efe_elem_begin(t, &mut n, &mut is101);
            // Point-element skip precedes the boundary test, as in the
            // original.
            if is101 != 0 {
                continue;
            }
            if efe_elem_active() == 0 {
                continue;
            }
            efe_elem_compute();
        }

        efe_normalize();
        efe_report();
        efe_finish();
    }
}
