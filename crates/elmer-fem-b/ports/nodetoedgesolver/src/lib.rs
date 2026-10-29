//! Port of `fem/src/modules/NodeToEdgeSolver.F90` (Elmer @ 9f6af2f85,
//! the lane fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/NodeToEdgeSolver.F90` — four F77-global
//! entries `extrudedrestart_init0_`, `extrudedrestart_`,
//! `nodetoedgefield_init0_`, `nodetoedgefield_`. The emitted dynamic
//! symbols are exactly those four (proven by `nm -D` on the reference
//! `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/nodetoedgesolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MAXVAL`, `COUNT`), formatted `WRITE`s, `PRINT`s, whole-array
//! updates, and the `GOTO 1` imaginary-part redo stay in Fortran wrappers
//! (or in the same Rust loop order), so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors NodeToEdgeSolver.F90 statement for
//! statement):
//! - The original has no `SAVE` state except the `LocalMatrix`
//!   `TYPE(Nodes_t), SAVE :: Nodes`, mirrored as wrapper-side SAVE state
//!   with identical lifetime (refilled by `GetElementNodes` on every
//!   element, as in the original). All other per-call allocatables
//!   (`Stiff`/`Force`/`Anodal`, the `ExtrudedRestart` handles) are
//!   per-invocation scratch, sized from `MaxElementDOFs` exactly like the
//!   original's `ALLOCATE` (plus a grow-only guard that only triggers
//!   where the original would overflow silently; the leading sections
//!   used by the wrappers are unaffected, so valid inputs never observe
//!   it).
//! - The `STAT=istat` allocation check is unchecked in the original, so
//!   the port's `try_reserve` failure calls the same-caller `Fatal`
//!   instead of aborting blindly.
//! - `ExtrudedRestart`'s `VariableAddVector` adoption (which takes over
//!   the freshly allocated `Perm`/`Values`) is mirrored by wrapper-side
//!   SAVE handles that are re-nulled on every variable, exactly matching
//!   the original's automatic (de)allocation lifetimes.
//! - The `PeriodicFlip` `WHERE` block (including its `dofs > 1`
//!   uninitialized-index shape) is transcribed verbatim into the
//!   wrapper; the proof stays on the defined scalar path.
//! - The edge-target copy loop (with its permutation-mismatch `Fatal`s)
//!   stays in the wrapper; with no `Edge Variable` present both arms take
//!   the `Warn` path, as in the reference runs.
//! - `Model`, `dt`, `TransientSimulation` are unused by the
//!   `NodeToEdgeField` body (received and ignored, `Model` is only used
//!   by `ExtrudedRestart` for the `Solvers` handle), exactly as in the
//!   original.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn nte_xinit0();
    fn nte_xbegin(h_model: *const c_void, h_solver: *const c_void);
    fn nte_xstart(novar: *mut c_int);
    fn nte_xmap(i: c_int);
    fn nte_xfinish();
    fn nte_ninit0();
    fn nte_nbegin(h_solver: *const c_void);
    fn nte_nsetup(dofs: *mut c_int, dim: *mut c_int, iscomplex: *mut c_int);
    fn nte_nedge();
    fn nte_nstyle();
    fn nte_nmaxdofs() -> c_int;
    fn nte_npass(second: c_int);
    fn nte_ninit();
    fn nte_nactive() -> c_int;
    fn nte_nelem(t: c_int, n: *mut c_int, nd: *mut c_int);
    fn nte_nanodal(ld: c_int, n: c_int, anodal: *mut c_double);
    fn nte_naflip(ld: c_int, n: c_int, anodal: *mut c_double);
    fn nte_nip(ld: c_int, stiff: *mut c_double, force: *mut c_double) -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn nte_npoint(
        t: c_int,
        ld: c_int,
        n: c_int,
        nd: c_int,
        anodal: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn nte_nupdate(ld: c_int, stiff: *mut c_double, force: *mut c_double);
    fn nte_nbulk();
    fn nte_nsolve() -> c_double;
    fn nte_ntail();
    fn nte_ndone();
    fn nte_nfatal_alloc();
}

/// Grow the scratch arrays when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its comment assumes the first mesh is big enough); growth only
/// triggers where the original would overflow silently, and the leading
/// sections used by the wrappers are unaffected, so valid inputs never
/// observe it.
fn ensure_state(stiff: &mut Vec<c_double>, force: &mut Vec<c_double>, anodal: &mut Vec<c_double>, need: usize) {
    if need * need <= stiff.len() {
        return;
    }
    if stiff.try_reserve(need * need - stiff.len()).is_err()
        || force.try_reserve(need - force.len()).is_err()
        || anodal.try_reserve(3 * need - anodal.len()).is_err()
    {
        unsafe { nte_nfatal_alloc() };
        std::process::abort();
    }
    stiff.resize(need * need, 0.0);
    force.resize(need, 0.0);
    anodal.resize(3 * need, 0.0);
}

/// Sharing guard for the wrapper-side SAVE handles (never contended:
/// the entries run serially).
static STATE: Mutex<()> = Mutex::new(());

/// `extrudedrestart_init0_` — flags the solver mesh a local copy with no
/// output. Argument order and by-reference travel match the Fortran dummy
/// list `(Model, Solver, dt, TransientSimulation)`; only the solver
/// context (via `GetSolverParams`) is used, as in the original.
///
/// # Safety
///
/// `model`, `solver`, `dt`, `transient` must be valid reads if non-null
/// (they are received and ignored).
#[export_name = "extrudedrestart_init0_"]
pub unsafe extern "C" fn extrudedrestart_init0_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe { nte_xinit0() };
}

/// `extrudedrestart_` — maps result variables between the solver meshes.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "extrudedrestart_"]
pub unsafe extern "C" fn extrudedrestart_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe { nte_xbegin(model, solver) };
    let mut novar: c_int = 0;
    unsafe { nte_xstart(&mut novar) };
    for i in 1..=novar {
        unsafe { nte_xmap(i) };
    }
    unsafe { nte_xfinish() };
}

/// `nodetoedgefield_init0_` — sets the edge-element defaults for the
/// solver. Argument order and by-reference travel match the Fortran dummy
/// list `(Model, Solver, dt, TransientSimulation)`; only the solver
/// context (via `GetSolverParams`) is used, as in the original.
///
/// # Safety
///
/// `model`, `solver`, `dt`, `transient` must be valid reads if non-null
/// (they are received and ignored).
#[export_name = "nodetoedgefield_init0_"]
pub unsafe extern "C" fn nodetoedgefield_init0_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe { nte_ninit0() };
}

/// `nodetoedgefield_` — projects a nodal field to an Hcurl field.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `model`, `dt`, `transient` must be valid reads if non-null (they are
/// ignored).
#[export_name = "nodetoedgefield_"]
pub unsafe extern "C" fn nodetoedgefield_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe { nte_nbegin(solver) };
    let (mut dofs, mut dim, mut iscomplex) = (0, 0, 0);
    unsafe { nte_nsetup(&mut dofs, &mut dim, &mut iscomplex) };
    let _ = (dofs, dim);
    unsafe { nte_nedge() };
    unsafe { nte_nstyle() };

    // First-mesh permanent sizing (bug-compatible: sized once).
    let mut ld = unsafe { nte_nmaxdofs() } as usize;
    let mut stiff = Vec::<c_double>::new();
    let mut force = Vec::<c_double>::new();
    let mut anodal = Vec::<c_double>::new();
    if stiff.try_reserve(ld * ld).is_err()
        || force.try_reserve(ld).is_err()
        || anodal.try_reserve(3 * ld).is_err()
    {
        unsafe { nte_nfatal_alloc() };
        std::process::abort();
    }
    stiff.resize(ld * ld, 0.0);
    force.resize(ld, 0.0);
    anodal.resize(3 * ld, 0.0);

    // The imaginary-part redo (`GOTO 1`): a second pass with the
    // imaginary offset reusing the constant bulk matrix.
    let passes: c_int = if iscomplex != 0 { 2 } else { 1 };
    for pass in 0..passes {
        unsafe { nte_npass(pass) };
        // 1 CALL DefaultInitialize(Solver, ReadySystemMatrix)
        unsafe { nte_ninit() };

        let active = unsafe { nte_nactive() };
        for t in 1..=active {
            let (mut n, mut nd) = (0, 0);
            unsafe { nte_nelem(t, &mut n, &mut nd) };
            let need = (n as usize).max(nd as usize);
            if need > ld {
                ensure_state(&mut stiff, &mut force, &mut anodal, need);
                ld = need;
            }
            let ldc = ld as c_int;
            unsafe { nte_nanodal(ldc, n, anodal.as_mut_ptr()) };
            unsafe { nte_naflip(ldc, n, anodal.as_mut_ptr()) };
            // Get element local matrix and rhs vector:
            //----------------------------------------
            let nip = unsafe { nte_nip(ldc, stiff.as_mut_ptr(), force.as_mut_ptr()) };
            for tt in 1..=nip {
                unsafe {
                    nte_npoint(
                        tt,
                        ldc,
                        n,
                        nd,
                        anodal.as_mut_ptr(),
                        stiff.as_mut_ptr(),
                        force.as_mut_ptr(),
                    )
                };
            }
            // Update global matrix and rhs vector from local matrix & vector:
            //---------------------------------------------------------------
            unsafe { nte_nupdate(ldc, stiff.as_mut_ptr(), force.as_mut_ptr()) };
        }

        unsafe { nte_nbulk() };
        // And finally, solve (Norm is local and unused, as in the original).
        let _norm: c_double = unsafe { nte_nsolve() };
        unsafe { nte_ntail() };
    }

    unsafe { nte_ndone() };
}
