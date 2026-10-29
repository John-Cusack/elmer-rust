//! Port of `fem/src/modules/NormalSolver.F90` (Elmer @ 9f6af2f85, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/NormalSolver.F90` — two F77-global
//! entries `normalsolver_` and `normalsolver_init_` with two `CONTAINS`
//! internals (`BulkAssembly`, `EnforceSymmetry`). The emitted dynamic
//! symbols are exactly `normalsolver_` + `normalsolver_init_` (proven by
//! `nm -D` on the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/normalsolver_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Order-sensitive reductions (`SUM`,
//! `DOT_PRODUCT`-style double loops), formatted `WRITE`s to `Message`,
//! and whole-array/strided updates stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors NormalSolver.F90 statement for statement):
//! - `SAVE` semantics replicated: entry state (solver/mesh/params/result
//!   handles, `ForceVector`+`SaveRHS`, `Values`, `at0/at1`, the dead
//!   `Visited` flag) and the internals' SAVE `Nodes` live as wrapper-side
//!   SAVE state with identical lifetimes (consumed before the next fill,
//!   as in the original); `BulkAssembly`'s per-call `IntegStuff` is
//!   re-filled per element, matching the original's per-element lifetime.
//! - Allocation has no `STAT=` check in the original, so plain
//!   `ALLOCATE` is mirrored (no `Fatal`-on-alloc path exists to match).
//! - The component loop (`RHS => ForceVector(:,i)`, the `SetD` rename +
//!   `Values` restore + `DefaultDirichletBCs`, `Values = 0`, the
//!   discarded `DefaultSolve` norm, the `(i::DOFs)` copy) is driven from
//!   Rust with one wrapper per statement group; all pointer traffic and
//!   `CHARACTER` manipulation stays in Fortran.
//! - Normalization, the vector-projection `BLOCK` (including its `GotIt`
//!   overwrite quirk), and `EnforceSymmetry` run wholly in wrappers so
//!   their reductions, `WRITE`s, and whole-array updates are bit-exact.
//! - Defined-case `DOFs` handling: the original threads a `Dofs` local
//!   that is only assigned when the result variable pre-exists (the
//!   `DefaultVariableAdd` branch leaves it unread, which is undefined and
//!   crashes the reference in practice). The wrappers read
//!   `NrmSol % DOFs` at each use site instead — the same value in every
//!   defined case, deterministic where the original is undefined. Proof
//!   configs always take the pre-existing-variable branch (as does the
//!   upstream `normals` test, via `Exported Variable 1`).
//! - `dt`/`TransientSimulation` cross by value into the projection
//!   wrapper (the only reader, exactly as in the original); `Model` is
//!   received and ignored, exactly as in the original.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn noa_begin(h_solver: *const c_void, dim: *mut c_int);
    fn noa_info_head();
    fn noa_assoc_matrix() -> c_int;
    fn noa_active_perm() -> c_int;
    fn noa_result_variable();
    fn noa_constant_bulk(cb: *mut c_int, inuse: *mut c_int);
    fn noa_default_initialize();
    fn noa_tick0();
    fn noa_forcevec_alloc();
    fn noa_forcevec_free();
    fn noa_finish_bulk(cb: c_int);
    fn noa_finish_assembly();
    fn noa_info_assembly();
    fn noa_info_solution();
    fn noa_setd_begin() -> c_int;
    fn noa_setd_free(setd: c_int);
    fn noa_comp_rhs(i: c_int);
    fn noa_comp_dirichlet(i: c_int, setd: c_int);
    fn noa_comp_solve() -> c_double;
    fn noa_comp_copy(i: c_int);
    fn noa_rhs_restore();
    fn noa_normalize();
    fn noa_bulk_alloc();
    fn noa_bulk_free();
    fn noa_bulk_nelems() -> c_int;
    fn noa_bulk_elem(elem: c_int);
    fn noa_bulk_zero();
    fn noa_bulk_ipbegin() -> c_int;
    fn noa_bulk_ipoint(t: c_int);
    fn noa_bulk_tomatrix();
    fn noa_bulk_toforce();
    fn noa_project(dt: c_double, transient: c_int);
    fn noa_want_symmetry() -> c_int;
    fn noa_enforce_symmetry();
    fn noa_init_begin(h_solver: *const c_void);
    fn noa_init_variable();
    fn noa_init_dim() -> c_int;
    fn noa_init_result(dim: c_int);
    fn noa_init_timeorder();
}

static STATE: Mutex<()> = Mutex::new(());

/// `normalsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` is unused by the
/// body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they reach only the
/// projection block, as in the original).
#[export_name = "normalsolver_"]
pub unsafe extern "C" fn normalsolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    let mut dim: c_int = 0;
    unsafe { noa_begin(solver, &mut dim) };

    unsafe { noa_info_head() };

    // IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { noa_assoc_matrix() } == 0 {
        return;
    }
    // IF (COUNT(Solver % Variable % Perm > 0) <= 0) RETURN
    if unsafe { noa_active_perm() } <= 0 {
        return;
    }

    unsafe { noa_result_variable() };

    let (mut cb, mut inuse) = (0, 0);
    unsafe { noa_constant_bulk(&mut cb, &mut inuse) };
    unsafe { noa_default_initialize() };

    unsafe { noa_tick0() };
    unsafe { noa_forcevec_alloc() };

    // BulkAssembly.
    unsafe {
        noa_bulk_alloc();
    }
    let nelems = unsafe { noa_bulk_nelems() };
    for elem in 1..=nelems {
        unsafe { noa_bulk_elem(elem) };
        unsafe { noa_bulk_zero() };
        let nip = unsafe { noa_bulk_ipbegin() };
        for t in 1..=nip {
            unsafe { noa_bulk_ipoint(t) };
        }
        if inuse == 0 {
            unsafe { noa_bulk_tomatrix() };
        }
        unsafe { noa_bulk_toforce() };
    }
    unsafe { noa_bulk_free() };

    if cb != 0 {
        unsafe { noa_finish_bulk(1) };
    } else {
        unsafe { noa_finish_bulk(0) };
    }
    unsafe { noa_finish_assembly() };

    unsafe { noa_info_assembly() };

    let setd = unsafe { noa_setd_begin() };
    for i in 1..=dim {
        unsafe {
            noa_comp_rhs(i);
            noa_comp_dirichlet(i, setd);
        }
        // And finally, solve (Norm is local and unused, as in the original).
        let _norm: c_double = unsafe { noa_comp_solve() };
        unsafe { noa_comp_copy(i) };
    }
    unsafe {
        noa_rhs_restore();
        noa_setd_free(setd);
    }

    // Normalize the computed normals so that the length is always one.
    unsafe { noa_normalize() };

    unsafe { noa_forcevec_free() };

    // Optionally use the solver to project a vector field to the normal
    // direction (whole BLOCK lives in the wrapper, GotIt quirk included).
    let dtv: c_double = if dt.is_null() { 0.0 } else { unsafe { *dt } };
    let tr: c_int = if transient.is_null() {
        0
    } else {
        unsafe { *transient }
    };
    unsafe { noa_project(dtv, tr) };

    if unsafe { noa_want_symmetry() } != 0 {
        unsafe { noa_enforce_symmetry() };
    }

    unsafe { noa_info_solution() };
}

/// `normalsolver_init_` — the solver init entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; only the solver params
/// list is touched, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "normalsolver_init_"]
pub unsafe extern "C" fn normalsolver_init_entry(
    _model: *const c_void,
    solver: *const c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe { noa_init_begin(solver) };
    unsafe { noa_init_variable() };
    let dim = unsafe { noa_init_dim() };
    unsafe { noa_init_result(dim) };
    unsafe { noa_init_timeorder() };
}
