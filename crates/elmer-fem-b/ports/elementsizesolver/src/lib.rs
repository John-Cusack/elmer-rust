//! Port of `fem/src/modules/ElementSizeSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/ElementSizeSolver.F90` — one F77-global
//! entry (`elementsizesolver_`) with one `CONTAINS` internal. The emitted
//! dynamic symbol is exactly `elementsizesolver_` (proven by `nm -D` on
//! the reference `.so`); the internal emits only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/elementsize_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (`**` powers, `MIN`/`MAX` reductions, the
//! elemental matrix/vector updates) and the `ES12.4` reports stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors ElementSizeSolver.F90 statement for
//! statement):
//! - `STIFF`/`FORCE` are allocated at `Mesh % MaxElementNodes` sizing on
//!   every invocation and deallocated at the end; the `SAVE`d
//!   `AllocationsDone` flag the original sets-but-never-reads is not
//!   replicated (no observable state survives the deallocate).
//! - The internal procedure's `SAVE Nodes` lives as wrapper-side SAVE
//!   state (consumed before the next fill, as in the original).
//! - `Power = 1.0 / ListGetCReal(...)` defaults to `1.0` when the keyword
//!   is absent; IEEE division, exact.
//! - Allocation failure calls the same `Fatal('PoissonSolve', ...)` the
//!   original calls (bug-compatible program name) via `try_reserve`.
//! - `Norm = DefaultSolve()` is kept and discarded; `Model`, `dt`, and
//!   the transient flag are unused by the body, exactly as in the
//!   original.
//! - No Rust threads: the solver entry runs serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn ees_begin(h_solver: *mut c_void);
    fn ees_info_start();
    fn ees_max_elem_nodes() -> c_int;
    fn ees_params_begin();
    fn ees_read_exponent(power: *mut c_double, gotit: *mut c_int);
    fn ees_read_noweight(val: *mut c_int, gotit: *mut c_int);
    fn ees_sweep_begin(power: c_double, noweigh: c_int);
    fn ees_nof_active() -> c_int;
    fn ees_default_initialize();
    fn ees_begin_active(t: c_int, n: *mut c_int);
    fn ees_lm_begin();
    fn ees_lm_ip_n() -> c_int;
    fn ees_lm_ip(t: c_int, n: c_int, stiff: *mut c_double, force: *mut c_double);
    fn ees_default_update(nmax: c_int, g: *mut c_double, f: *mut c_double);
    fn ees_default_finish_bulk();
    fn ees_default_finish();
    fn ees_default_dirichlet();
    fn ees_default_solve() -> c_double;
    fn ees_report();
    fn ees_fatal_alloc();
}

/// `elementsizesolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "elementsizesolver_"]
pub unsafe extern "C" fn elementsizesolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ees_begin(solver);
        ees_info_start();
        let nmax = ees_max_elem_nodes() as usize;
        let mut stiff = Vec::<c_double>::new();
        let mut force = Vec::<c_double>::new();
        if stiff.try_reserve(nmax * nmax).is_err() || force.try_reserve(nmax).is_err() {
            ees_fatal_alloc();
        }
        stiff.resize(nmax * nmax, 0.0);
        force.resize(nmax, 0.0);
        ees_params_begin();
        let mut raw: c_double = 0.0;
        let mut gotit: c_int = 0;
        ees_read_exponent(&mut raw, &mut gotit);
        let power: c_double = if gotit != 0 { 1.0 / raw } else { 1.0 };
        let mut noweight: c_int = 0;
        let mut _found: c_int = 0;
        ees_read_noweight(&mut noweight, &mut _found);
        ees_sweep_begin(power, noweight);
        ees_default_initialize();
        let active = ees_nof_active();
        let mut n: c_int = 0;
        for t in 1..=active {
            ees_begin_active(t, &mut n);
            ees_lm_begin();
            let nn = n as usize;
            stiff[..nn * nn].fill(0.0);
            force[..nn].fill(0.0);
            let nip = ees_lm_ip_n();
            for tt in 1..=nip {
                ees_lm_ip(tt, n, stiff.as_mut_ptr(), force.as_mut_ptr());
            }
            ees_default_update(nmax as c_int, stiff.as_mut_ptr(), force.as_mut_ptr());
        }
        ees_default_finish_bulk();
        ees_default_finish();
        ees_default_dirichlet();
        let _norm: c_double = ees_default_solve();
        ees_report();
    }
}
