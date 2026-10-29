//! Port of `fem/src/modules/EliminatePeriodic.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the linear-solve hook boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/EliminatePeriodic.F90` — one recursive
//! INTEGER FUNCTION `EliminatePeriodic` (8 dummies, int result). The
//! emitted dynamic symbol is exactly `eliminateperiodic_` (proven by
//! `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, linsolve-hook gate shape): the hook is reached
//! through the loader's C trampoline (`DoLinSolveProcs` in `Load.c`),
//! which calls the raw address with eight plain pointers and takes an
//! int result — so the entry is a plain `extern "C"` function. The
//! derived-type dummies (`Model`, `Solver`, the `Matrix_t` pointer)
//! cross as opaque addresses. Every component access and every Elmer
//! library call goes through the generated-style Fortran `BIND(C)`
//! wrappers in `accessors/eliminateperiodic_acc.F90`, compiled against
//! the built `.mod` files of the target configure — the port itself
//! contains no layout assumption and no pinned offset.
//!
//! Faithfulness notes (mirrors EliminatePeriodic.F90 statement for
//! statement):
//! - `SAVE` design: the original's `ColN` is a non-SAVE allocatable
//!   local filled on the first call and read on every later call; the
//!   wrapper keeps it as the SAVE module array `acc_coln` with identical
//!   contents (the author's evident intent), so re-entry is
//!   deterministic. The `A % EMatrix` / `C % EMatrix` links persist
//!   through the matrix exactly as in the original, and the re-entry
//!   refresh (`C_Trans` / `RedA`) follows the original's pointer
//!   chain verbatim — including that `RedA` starts disassociated on
//!   re-entry and is rebuilt by the multiply.
//! - The no-periodic first call (original 184-187) and the
//!   reduced-system solve (original 322) both call the built-tree
//!   `EliminateDirichlet`, exactly like the reference `.so` (which
//!   carries `NEEDED EliminateDirichlet.so`); both arms of the proof
//!   therefore delegate to the identical routine.
//! - Structure surgery (projection-matrix build, CAC^T reduction), the
//!   eigen setup/map-back, the Norm `SQRT(SUM(x**2)/n)` reduction, and
//!   all formatted `WRITE`/`PRINT`/`Info` traffic stay in Fortran, so
//!   bit-identity holds by construction wherever the call sequence
//!   matches. The reduced-system Perm remap (original 303-313) runs here
//!   through scalar accessors with the original's evaluation order
//!   reproduced exactly (`DOFs*(j-1)+1`, `(ColN-1)/DOFs+1` truncating
//!   division on positive values, exact `> 0` comparisons).
//! - The `Level=5` timing reports print only when the `.sif` allows
//!   output level 5 (the proof configs use `Max Output Level = 4`, which
//!   keeps the deterministic lines); the `CPUtime()` reads are still
//!   issued inside the wrappers, exactly as in the original.
//! - `Model` is unused by the body except as the delegation handle,
//!   exactly as in the original (received and passed through).
//! - No Rust threads: the hook runs serially inside the solve, so there
//!   is no `static` state on the Rust side at all.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn elp_begin(
        h_model: *const c_void,
        h_solver: *const c_void,
        h_a: *const c_void,
        h_b: *const c_void,
        h_x: *const c_void,
        n: c_int,
    );
    fn elp_need_build() -> c_int;
    fn elp_delegated() -> c_int;
    fn elp_build(dofs: c_int, n: c_int, b: *mut c_double, x: *mut c_double, norm: *mut c_double)
        -> c_int;
    fn elp_relink();
    fn elp_eig_setup() -> c_int;
    fn elp_reduce();
    fn elp_pre_solve();
    fn elp_perm_alloc();
    fn elp_vperm_size() -> c_int;
    fn elp_vperm_get(i: c_int) -> c_int;
    fn elp_perm_set(i: c_int, v: c_int);
    fn elp_coln_get(d: c_int) -> c_int;
    fn elp_solve(dofs: c_int, norm: *mut c_double);
    fn elp_mapback(n: c_int, norm: *mut c_double);
    fn elp_finish() -> c_int;
}

/// The reduced-system Perm remap (original lines 303-313): for each
/// solver-dof row with a positive variable perm, map it through the
/// projection column map. `acc_perm` is zero-filled by `elp_perm_alloc`;
/// entries with no kept column stay zero, exactly as in the original.
///
/// # Safety
///
/// The wrapper-side `acc_perm` storage must be set up (`elp_perm_alloc`);
/// reads go through the entry's solver/col maps.
unsafe fn perm_remap(dofs: c_int) {
    unsafe {
        let n = elp_vperm_size();
        let mut i: c_int = 1;
        while i <= n {
            let j = elp_vperm_get(i);
            if j > 0 {
                let dof = dofs * (j - 1) + 1;
                if elp_coln_get(dof) > 0 {
                    elp_perm_set(i, (elp_coln_get(dof) - 1) / dofs + 1);
                }
            }
            i += 1;
        }
    }
}

/// `eliminateperiodic_` — the `Before Linsolve` hook. Argument order
/// and by-reference travel match the Fortran dummy list
/// `(Model, Solver, A, b, x, n, DOFs, Norm)` as passed by the loader's C
/// trampoline; the INTEGER FUNCTION result returns the hook status (1 =
/// system already solved, skip the normal solve).
///
/// # Safety
///
/// `solver`/`a` must be the Elmer `Solver_t`/`Matrix_t` addresses the
/// loader passes; `b`/`x` must point at `n` readable/writable
/// `c_double`s; `n`/`dofs` must be valid reads; `norm` must be a valid
/// write.
#[export_name = "eliminateperiodic_"]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn eliminateperiodic_entry(
    model: *mut c_void,
    solver: *mut c_void,
    a: *mut c_void,
    b: *mut c_double,
    x: *mut c_double,
    n: *const c_int,
    dofs: *const c_int,
    norm: *mut c_double,
) -> c_int {
    unsafe {
        // `A` is a `TYPE(Matrix_t), POINTER` dummy: gfortran passes the
        // address of the pointer variable, so dereference once to reach
        // the matrix (verified empirically against gfortran 13.3.0).
        // All other dummies are non-pointer: plain addresses.
        let matrix = *(a as *const *const c_void);
        elp_begin(
            model,
            solver,
            matrix,
            b as *const c_void,
            x as *const c_void,
            *n,
        );
        if elp_need_build() != 0 {
            let rv = elp_build(*dofs, *n, b, x, norm);
            if elp_delegated() != 0 {
                return rv;
            }
        } else {
            elp_relink();
        }
        elp_eig_setup();
        elp_reduce();
        elp_pre_solve();
        elp_perm_alloc();
        perm_remap(*dofs);
        elp_solve(*dofs, norm);
        elp_mapback(*n, norm);
        elp_finish()
    }
}
