//! Port of `fem/src/modules/EliminateDirichlet.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the linear-solve hook boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/EliminateDirichlet.F90` — one recursive
//! INTEGER FUNCTION `EliminateDirichlet` (8 dummies, int result). The
//! emitted dynamic symbol is exactly `eliminatedirichlet_` (proven by
//! `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, linsolve-hook gate shape): the hook is reached
//! through the loader's C trampoline (`DoLinSolveProcs` in `Load.c`),
//! which calls the raw address with eight plain pointers and takes an
//! int result — so the entry is a plain `extern "C"` function. The
//! derived-type dummies (`Model`, `Solver`, the `Matrix_t` pointer) cross
//! as opaque addresses. Every component access and every Elmer library
//! call goes through the generated-style Fortran `BIND(C)` wrappers in
//! `accessors/eliminatedirichlet_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset.
//!
//! Faithfulness notes (mirrors EliminateDirichlet.F90 statement for
//! statement):
//! - `SAVE Bmatrix, AllocationsDone`: `Bmatrix` is re-pointed at
//!   `A % EMatrix` on every entry and `A % EMatrix` is freed + nullified
//!   before return, so no cross-call state survives; `AllocationsDone`
//!   is set-but-never-read in the original and is not replicated (no
//!   observable state). The wrapper module holds the per-call handles
//!   (`Solver`, `A`, `b`, `x`, `BMatrix`, `Permutation`, `F`, `U`) with
//!   identical lifetimes (consumed before the next fill, as in the
//!   original).
//! - The commented-out early return (`Dirichlet <= 0`) is dead code in
//!   the original and stays absent here.
//! - Structure surgery (reduced-matrix build), the complex fill, the
//!   eigen map-back, and all Elmer library calls stay in Fortran
//!   wrappers, so bit-identity holds by construction wherever the call
//!   sequence matches. The real fill and the real map-back run here
//!   through scalar accessors with the original's evaluation order
//!   reproduced exactly (`(Values(l)*b(t))/Values(Diag(t))`, no
//!   reassociation; exact `==`/`/= 0` comparisons).
//! - The `Level=5` timing reports print only when the `.sif` allows
//!   output level 5 (the proof configs use `Max Output Level = 4`, which
//!   keeps the deterministic lines); the `CPUtime()` reads are still
//!   issued inside the wrappers, exactly as in the original.
//! - `Model` is unused by the body, exactly as in the original
//!   (received and ignored).
//! - No Rust threads: the hook runs serially inside the solve, so there
//!   is no `static` state on the Rust side at all.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn eld_begin(
        h_solver: *const c_void,
        h_a: *const c_void,
        h_b: *const c_void,
        h_x: *const c_void,
        n: c_int,
    );
    fn eld_mass() -> c_int;
    fn eld_damp() -> c_int;
    fn eld_ensure_md(mass: c_int, damp: c_int);
    fn eld_need_build() -> c_int;
    fn eld_build(mass: c_int, damp: c_int);
    fn eld_setup_fu();
    fn eld_is_complex() -> c_int;
    fn eld_a_nrows() -> c_int;
    fn eld_a_row(i: c_int) -> c_int;
    fn eld_a_row1(i: c_int) -> c_int;
    fn eld_a_col(l: c_int) -> c_int;
    fn eld_a_val(l: c_int) -> c_double;
    fn eld_a_mass(l: c_int) -> c_double;
    fn eld_a_damp(l: c_int) -> c_double;
    fn eld_a_diagval(t: c_int) -> c_double;
    fn eld_a_is_diag(i: c_int, l: c_int) -> c_int;
    fn eld_perm(i: c_int) -> c_int;
    fn eld_b_get(t: c_int) -> c_double;
    fn eld_x_get(i: c_int) -> c_double;
    fn eld_x_set(i: c_int, v: c_double);
    fn eld_f_get(j: c_int) -> c_double;
    fn eld_f_set(j: c_int, v: c_double);
    fn eld_u_get(j: c_int) -> c_double;
    fn eld_u_set(j: c_int, v: c_double);
    fn eld_bm_val_set(k: c_int, v: c_double);
    fn eld_bm_mass_set(k: c_int, v: c_double);
    fn eld_bm_damp_set(k: c_int, v: c_double);
    fn eld_fill_complex(n: c_int, b: *mut c_double, x: *mut c_double, mass: c_int, damp: c_int);
    fn eld_pre_solve();
    fn eld_n_eigen() -> c_int;
    fn eld_mapback_eigen(n: c_int, b: *mut c_double);
    fn eld_solve(dofs: c_int, nrm: *mut c_double);
    fn eld_finish() -> c_int;
}

/// The real-system fill (original lines 290-319): scatter `B`/`X` into
/// the reduced right-hand side and initial guess, folding Dirichlet
/// columns into `F` and copying the kept block into `Bmatrix % Values`.
/// `k` is the running reduced-matrix write cursor (starts at 1).
///
/// # Safety
///
/// The wrapper-side `F`, `U`, and `Bmatrix` storage must be set up
/// (`eld_setup_fu`); reads/writes go through the entry's `b`/`x`.
unsafe fn fill_real(mass: bool, damp: bool) {
    unsafe {
        let na = eld_a_nrows();
        let mut j: c_int = 0;
        let mut k: c_int = 1;
        let mut i: c_int = 1;
        while i <= na {
            if eld_perm(i) != 0 {
                j += 1;
                eld_f_set(j, eld_b_get(i));
                eld_u_set(j, eld_x_get(i));
                let r0 = eld_a_row(i);
                let r1 = eld_a_row1(i);
                let mut l = r0;
                while l < r1 {
                    let t = eld_a_col(l);
                    if eld_perm(t) == 0 {
                        let dec = eld_a_val(l) * eld_b_get(t) / eld_a_diagval(t);
                        eld_f_set(j, eld_f_get(j) - dec);
                    } else {
                        let mut add = eld_a_val(l) != 0.0;
                        if mass {
                            add |= eld_a_mass(l) != 0.0;
                        }
                        if damp {
                            add |= eld_a_damp(l) != 0.0;
                        }
                        add |= eld_a_is_diag(i, l) != 0;
                        if add {
                            eld_bm_val_set(k, eld_a_val(l));
                            if mass {
                                eld_bm_mass_set(k, eld_a_mass(l));
                            }
                            if damp {
                                eld_bm_damp_set(k, eld_a_damp(l));
                            }
                            k += 1;
                        }
                    }
                    l += 1;
                }
            }
            i += 1;
        }
    }
}

/// The real-system map-back (original lines 416-425): Dirichlet rows
/// take `B(i)/Values(Diag(i))`, kept rows take the reduced solution.
/// # Safety
///
/// Same contract as [`fill_real`]; runs before `eld_finish`.
unsafe fn mapback_real() {
    unsafe {
        let na = eld_a_nrows();
        let mut j: c_int = 0;
        let mut i: c_int = 1;
        while i <= na {
            if eld_perm(i) == 0 {
                eld_x_set(i, eld_b_get(i) / eld_a_diagval(i));
            } else {
                j += 1;
                eld_x_set(i, eld_u_get(j));
            }
            i += 1;
        }
    }
}

/// `eliminatedirichlet_` — the `Before Linsolve` hook. Argument order
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
#[export_name = "eliminatedirichlet_"]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn eliminatedirichlet_entry(
    _model: *mut c_void,
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
        eld_begin(solver, matrix, b as *const c_void, x as *const c_void, *n);
        let mass = eld_mass() != 0;
        let damp = eld_damp() != 0;
        eld_ensure_md(mass as c_int, damp as c_int);
        if eld_need_build() != 0 {
            eld_build(mass as c_int, damp as c_int);
        }
        eld_setup_fu();
        if eld_is_complex() != 0 {
            eld_fill_complex(*n, b, x, mass as c_int, damp as c_int);
        } else {
            fill_real(mass, damp);
        }
        eld_pre_solve();
        eld_solve(*dofs, norm);
        if eld_n_eigen() > 0 {
            eld_mapback_eigen(*n, b);
        } else {
            mapback_real();
        }
        eld_finish()
    }
}
