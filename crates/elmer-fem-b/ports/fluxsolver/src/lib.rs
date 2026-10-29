//! Port of `fem/src/modules/FluxSolver.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T022-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/FluxSolver.F90` — two F77-global entries
//! (`fluxsolver_`, `fluxsolver_init_`) with three `CONTAINS` internals. The
//! emitted dynamic symbols are exactly those two (proven by `nm -D` on the
//! reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/fluxsolver_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the Galerkin assembly, `MATMUL`/`SUM` reductions,
//! `SQRT` magnitude post-processing) and the `RealTime`-based reports stay
//! in Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors FluxSolver.F90 statement for statement):
//! - The matrix/variable guard maps to the same silent early returns
//!   (no message, exactly as in the original).
//! - The `Dofs == 0` case maps to the same `Warn` + return.
//! - `SAVE Conductivity, Nodes` live as wrapper-side SAVE state with
//!   identical lifetimes (refilled per element before use, as in the
//!   original); `LocalJumps`' SAVE nodes likewise. `ForceVector` is
//!   allocated per invocation and deallocated after the solve loop,
//!   matching the original's allocate/deallocate pair; allocation failure
//!   calls the same `Fatal('FluxSolver', ...)` via `try_reserve`
//!   (the `STAT=istat` check), never aborting.
//! - `AutoNorm`/`Skip Compute Nonlinear Change` handling, the
//!   `Flux Variable` → `Target Variable` → `Target Expression` fallback
//!   chain (with `UsePot`), the `Flux Coefficient` default, the
//!   `Constant Bulk Matrix` + `BulkValues` reuse path, the per-component
//!   `DefaultSolve` loop with `ListRemove` on the last component, the
//!   `Enforce Positive Magnitude` clamp, and the `Flux_abs`/`Grad_abs`
//!   post-processing all run in the wrappers in the original order.
//! - The `Level=5` timing reports print only when the `.sif` allows output
//!   level 5 (verified content-identical in the proof configs via
//!   `Max Output Level = 4`, which keeps the deterministic `Level=4`
//!   `Result Norm` line).
//! - `Model`, `dt`, and the transient flag are unused by the entry bodies,
//!   exactly as in the original (received and ignored).
//! - No Rust threads: the solver entry runs serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn fls_begin(h_solver: *mut c_void);
    fn fls_info_start();
    fn fls_guard() -> c_int;
    fn fls_sweep_begin() -> c_int;
    fn fls_warn_nofield();
    fn fls_fatal_alloc();
    fn fls_alloc_n() -> c_int;
    fn fls_nof_active() -> c_int;
    fn fls_fields_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn fls_elem(
        elem: c_int,
        nmax: c_int,
        ndofs: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        locpot: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
    );
    fn fls_faces(nmax: c_int, ndofs: c_int, stiff: *mut c_double, force: *mut c_double);
    fn fls_finish();
    fn fls_solve() -> c_double;
    fn fls_init_begin(h_solver: *mut c_void);
    fn fls_init_body();
}

/// `fluxsolver_init_` — the `_Init` entry Elmer calls once at solver setup.
/// Mirrors `FluxSolver_Init` (requested-field defaults, `-nooutput`
/// variable, exported flux/gradient variables, time derivative order,
/// linear-system defaults).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "fluxsolver_init_"]
pub unsafe extern "C" fn fluxsolver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        fls_init_begin(solver);
        fls_init_body();
    }
}

/// `fluxsolver_` — the solver entry. Argument order and by-reference travel
/// match the Fortran dummy list `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`fluxsolver_init_entry`].
#[export_name = "fluxsolver_"]
pub unsafe extern "C" fn fluxsolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        fls_begin(solver);
        fls_info_start();
        if fls_guard() != 0 {
            return;
        }
        let dofs = fls_sweep_begin();
        if dofs == 0 {
            fls_warn_nofield();
            return;
        }
        let nmax = fls_alloc_n() as usize;
        let dofs_usize = dofs as usize;
        let mut stiff = Vec::<c_double>::new();
        let mut force = Vec::<c_double>::new();
        let mut locpot = Vec::<c_double>::new();
        let mut basis = Vec::<c_double>::new();
        let mut dbasis = Vec::<c_double>::new();
        if stiff.try_reserve(nmax * nmax).is_err()
            || force.try_reserve(dofs_usize * nmax).is_err()
            || locpot.try_reserve(nmax).is_err()
            || basis.try_reserve(nmax).is_err()
            || dbasis.try_reserve(nmax * 3).is_err()
        {
            fls_fatal_alloc();
        }
        stiff.resize(nmax * nmax, 0.0);
        force.resize(dofs_usize * nmax, 0.0);
        locpot.resize(nmax, 0.0);
        basis.resize(nmax, 0.0);
        dbasis.resize(nmax * 3, 0.0);
        if fls_fields_begin() != 0 {
            fls_fatal_alloc();
        }
        let nelem = fls_nof_active();
        for elem in 1..=nelem {
            fls_elem(
                elem,
                nmax as c_int,
                dofs,
                stiff.as_mut_ptr(),
                force.as_mut_ptr(),
                locpot.as_mut_ptr(),
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
            );
        }
        fls_faces(
            nmax as c_int,
            dofs,
            stiff.as_mut_ptr(),
            force.as_mut_ptr(),
        );
        fls_finish();
        let _norm: c_double = fls_solve();
    }
}
