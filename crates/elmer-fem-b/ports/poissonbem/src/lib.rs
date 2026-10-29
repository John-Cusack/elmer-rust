//! Port of `fem/src/modules/PoissonBEM.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` →
//! Rust cdylib with the same dynamic surface).
//!
//! Unit: `fortran:fem/src/modules/PoissonBEM.F90` — one F77-global entry
//! `poissonbemsolver_` plus the co-entries `matvec_` and `precond_`
//! (assumed-size array dummies, so they cross the C boundary as plain
//! pointers) and the `globmat` module descriptor (`Matrix`), which the
//! accessor mirrors to keep the `B __globmat_MOD_matrix` symbol. The
//! emitted dynamic surface is exactly `T matvec_`, `T poissonbemsolver_`,
//! `T precond_`, `B __globmat_MOD_matrix` (proven by `nm -D` on the
//! reference `.so`); the internal procedures emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/poissonbem_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset.
//! Order-sensitive reductions (`SUM`), formatted I/O (`PRINT`), the
//! vector-subscripted BC loads, and the HUTI iterative solve (including
//! its `AddrFunc` callback wiring) stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence
//! matches. The `matvec_`/`precond_` entries only forward `(n, u, v)`
//! to those wrappers; `n` is read from the HUTI `ipar(3)` slot exactly
//! as the original's `HUTI_NDIM` macro does.
//!
//! Faithfulness notes (mirrors PoissonBEM.F90 statement for statement):
//! - `SAVE` semantics replicated: all first-call-sized state
//!   (`PotentialPerm`, `BoundaryPerm`, `ElementNodes`, `VolumeSource`,
//!   `VolumeForce`, `PotentialKnown`, `Matrix`, `Pot`, `Flx`, `Load`,
//!   `Diagonal`, `BoundaryNodes`, `AllocationsDone`) lives as
//!   wrapper-side SAVE state with identical lifetimes; the
//!   per-call pointers (`Potential`, `ForceVector`) are re-associated in
//!   `pbb_begin` on every entry, exactly like the original's locals.
//! - Allocation failure calls the same `Fatal('PoissonBEMSolver', ...)`
//!   via the `STAT=istat` check (no Rust-side allocation exists to fail).
//! - `dt` and `TransientSimulation` are unused by the body, exactly as
//!   in the original (received and ignored); the dead `Norm` read is
//!   skipped (pure component read, no call); the unused `EquationName`
//!   query is kept as a touch call.
//! - `Skip Volume Computing` is the one config branch evaluated in Rust
//!   (the switch result crosses); the direct/iterative dispatch stays in
//!   the `SolveFull` transcription, driven by the same ValueList string.
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads. The
//!   guard is deliberately never held across `pbb_solve`: the HUTI
//!   iterative path calls back into `matvec_`/`precond_` on this same
//!   thread, and those entries lock the guard themselves (the mutex is
//!   not reentrant, so holding it across the solve would deadlock).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn pbb_begin(h_model: *const c_void, h_solver: *const c_void);
    fn pbb_allocate();
    fn pbb_zero();
    fn pbb_bc_pass();
    fn pbb_volume_pass();
    fn pbb_assemble_pass();
    fn pbb_diagonal_fix();
    fn pbb_print_assembly();
    fn pbb_solve();
    fn pbb_extract();
    fn pbb_skip_volume() -> c_int;
    fn pbb_volume_compute();
    fn pbb_finalize();
    fn pbb_matvec(n: c_int, u: *mut c_double, v: *mut c_double);
    fn pbb_precond(n: c_int, u: *mut c_double, v: *mut c_double);
}

/// Sharing guard for the wrapper-side SAVE state (never contended: the
/// solver entries run serially; never held across `pbb_solve`, which
/// re-enters `matvec_`/`precond_` on this thread).
static STATE: Mutex<()> = Mutex::new(());

/// `poissonbemsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "poissonbemsolver_"]
pub unsafe extern "C" fn poissonbem_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    {
        let _guard = STATE.lock().unwrap();
        // Lines 105-114: cartesian check (Fatal inside), pointer
        // re-association, EquationName touch.
        pbb_begin(model, solver);
        // Lines 116-170: first-call permutation + SAVE allocation.
        pbb_allocate();
        // Lines 175-181: at = CPUTime(); zero fills.
        pbb_zero();
        // Lines 190-221: BC scan (Potential/variable-name/Flux branches).
        pbb_bc_pass();
        // Lines 228-257: volume source + scatter into ForceVector.
        pbb_volume_pass();
        // Lines 262-276: boundary assembly.
        pbb_assemble_pass();
        // Lines 278-284: PotentialKnown diagonal fix.
        pbb_diagonal_fix();
        // Lines 287-288: assembly timing print.
        pbb_print_assembly();
    }
    // Lines 292-297: st = CPUTime(); direct/iterative solve (unguarded:
    // the iterative path re-enters matvec_/precond_ on this thread).
    pbb_solve();
    {
        let _guard = STATE.lock().unwrap();
        // Lines 301-327: Pot/Flx extraction + Flux variable scatter.
        pbb_extract();
        // Line 329: volume-computing switch.
        if pbb_skip_volume() == 0 {
            // Lines 331-345: interior potential from boundary data.
            pbb_volume_compute();
        }
        // Lines 350-356: norm, InvalidateVariable, solve timing print.
        pbb_finalize();
    }
}

/// `matvec_` — HUTI matrix-vector callback (`v = Matrix*u`). `n` is the
/// HUTI dimension slot (`ipar(3)`), read exactly as the original's
/// `HUTI_NDIM` macro does; the loop itself runs in the wrapper.
///
/// # Safety
///
/// `u`/`v` must point at `n` readable/writable doubles, `ipar` at a HUTI
/// parameter block with the dimension at index 2.
#[export_name = "matvec_"]
pub unsafe extern "C" fn matvec_entry(
    u: *mut c_double,
    v: *mut c_double,
    ipar: *mut c_int,
) {
    let _guard = STATE.lock().unwrap();
    let n = *ipar.add(2);
    pbb_matvec(n, u, v);
}

/// `precond_` — HUTI diagonal preconditioner (`u = v/diag(Matrix)`).
/// Same `ipar(3)` dimension convention as `matvec_`.
///
/// # Safety
///
/// Same contract as [`matvec_entry`].
#[export_name = "precond_"]
pub unsafe extern "C" fn precond_entry(
    u: *mut c_double,
    v: *mut c_double,
    ipar: *mut c_int,
) {
    let _guard = STATE.lock().unwrap();
    let n = *ipar.add(2);
    pbb_precond(n, u, v);
}
