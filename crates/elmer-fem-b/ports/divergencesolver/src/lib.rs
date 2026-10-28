//! Port of `fem/src/modules/DivergenceSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/DivergenceSolver.F90` — two F77-global
//! entries (`divergencesolver_`, `divergencesolver_init_`) with one
//! `CONTAINS` internal. The emitted dynamic symbols are exactly those two
//! (proven by `nm -D` on the reference `.so`); the internal emits only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/divergence_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the Galerkin assembly, `SQRT`, `MIN`/`MAX`
//! reductions, the `RealTime`-based reports) stays in Fortran wrappers,
//! so bit-identity holds by construction wherever the call sequence
//! matches.
//!
//! Faithfulness notes (mirrors DivergenceSolver.F90 statement for
//! statement):
//! - The variable/matrix/DOFs guard maps to the same `Warn` + early
//!   return / `Fatal` calls with identical messages; the `RETURN` paths
//!   return from the entry without assembling or solving.
//! - The `Divergence Coefficient` keyword read feeds only a commented-out
//!   line in the original; the pure read is documented as skipped (no
//!   side effects either way).
//! - The `Visited` SAVE flag the original sets-but-never-reads is not
//!   replicated (no observable state).
//! - `SAVE Coeff, Nodes`: `Nodes` lives as wrapper-side SAVE state
//!   (consumed before the next fill, as in the original); `Coeff` backs
//!   only the commented-out line and is documented as skipped.
//! - Work arrays are allocated at
//!   `MAX(MaxElementDOFs, MaxElementNodes)` sizing once per invocation
//!   and dropped at return, matching the original's per-call
//!   allocate/deallocate pair. Allocation failure calls the same
//!   `Fatal('DivergenceSolver', ...)` via `try_reserve`.
//! - The `Level=5` timing reports print only when the `.sif` allows
//!   output level 5 (verified content-identical in the proof configs via
//!   `Max Output Level = 4`, which keeps the deterministic `Level=4`
//!   `Result Norm` line).
//! - `Norm = DefaultSolve()` is kept and reported, exactly as in the
//!   original. `Model`, `dt`, and the transient flag are unused by the
//!   entry body itself (received and ignored), exactly as in the
//!   original.
//! - No Rust threads: the solver entry runs serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn dvs_begin(h_solver: *mut c_void);
    fn dvs_info_start();
    fn dvs_var_check() -> c_int;
    fn dvs_warn_zero();
    fn dvs_warn_nomatrix();
    fn dvs_fatal_dofs();
    fn dvs_fatal_novar();
    fn dvs_fatal_alloc();
    fn dvs_sweep_begin();
    fn dvs_alloc_n() -> c_int;
    fn dvs_nof_active_elems() -> c_int;
    fn dvs_elem(
        elem: c_int,
        nmax: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        vx: *mut c_double,
        vy: *mut c_double,
        vz: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
    );
    fn dvs_default_finish_bulk();
    fn dvs_default_finish();
    fn dvs_default_dirichlet();
    fn dvs_tail() -> c_double;
    fn dvs_init_begin(h_solver: *mut c_void);
    fn dvs_init_body();
}

/// `divergencesolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `DivergenceSolver_Init` (default `Variable`, time
/// derivative order, linear-system defaults).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "divergencesolver_init_"]
pub unsafe extern "C" fn divergencesolver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        dvs_init_begin(solver);
        dvs_init_body();
    }
}

/// `divergencesolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`divergencesolver_init_entry`].
#[export_name = "divergencesolver_"]
pub unsafe extern "C" fn divergencesolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        dvs_begin(solver);
        dvs_info_start();
        match dvs_var_check() {
            1 => {
                dvs_warn_zero();
                return;
            }
            2 => {
                dvs_warn_nomatrix();
                return;
            }
            3 => dvs_fatal_dofs(),
            4 => dvs_fatal_novar(),
            _ => {}
        }
        dvs_sweep_begin();
        let nmax = dvs_alloc_n() as usize;
        let mut stiff = Vec::<c_double>::new();
        let mut force = Vec::<c_double>::new();
        let mut vx = Vec::<c_double>::new();
        let mut vy = Vec::<c_double>::new();
        let mut vz = Vec::<c_double>::new();
        let mut basis = Vec::<c_double>::new();
        let mut dbasis = Vec::<c_double>::new();
        if stiff.try_reserve(nmax * nmax).is_err()
            || force.try_reserve(nmax).is_err()
            || vx.try_reserve(nmax).is_err()
            || vy.try_reserve(nmax).is_err()
            || vz.try_reserve(nmax).is_err()
            || basis.try_reserve(nmax).is_err()
            || dbasis.try_reserve(nmax * 3).is_err()
        {
            dvs_fatal_alloc();
        }
        stiff.resize(nmax * nmax, 0.0);
        force.resize(nmax, 0.0);
        vx.resize(nmax, 0.0);
        vy.resize(nmax, 0.0);
        vz.resize(nmax, 0.0);
        basis.resize(nmax, 0.0);
        dbasis.resize(nmax * 3, 0.0);
        let nelem = dvs_nof_active_elems();
        for elem in 1..=nelem {
            dvs_elem(
                elem,
                nmax as c_int,
                stiff.as_mut_ptr(),
                force.as_mut_ptr(),
                vx.as_mut_ptr(),
                vy.as_mut_ptr(),
                vz.as_mut_ptr(),
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
            );
        }
        dvs_default_finish_bulk();
        dvs_default_finish();
        dvs_default_dirichlet();
        let _norm: c_double = dvs_tail();
    }
}
