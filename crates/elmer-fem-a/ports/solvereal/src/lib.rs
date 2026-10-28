//! Port of `fem/src/solve_real.F90` (Elmer @ 9f6af2f85) for lane-femcore-b,
//! through the gfortran ABI shim (ADR-027, T017 shape A).
//!
//! Unit: `fortran:fem/src/solve_real.F90` — two file-scope external
//! subprograms. The emitted gfortran symbols are `solvelapack_` and
//! `solvelapacksym_`; every dummy travels by reference (a pointer to the
//! caller's argument).
//!
//! The unit is a LAPACK dense-solver wrapper: the LU factorization and
//! the triangular solves stay in LAPACK (`dgetrf_`/`dgetrs_`,
//! `dpotrf_`/`dpotrs_`), so the numerics are bit-identical by
//! construction. What the port owns is the wrapper itself — the argument
//! marshalling, the `N <= 0` early return and the `info` diagnostic
//! prints.
//!
//! Faithfulness notes:
//! - The Fortran declares `A(N*N)` and `X(N)` and passes `LDA = N` to
//!   LAPACK. The port passes the same `LDA = N`, so every shape behaves
//!   identically, not just the square systems the campaign's callers use.
//! - A nonzero `info` only prints: the Fortran has no `STOP` here, so the
//!   port prints the same list-directed line and continues into the next
//!   LAPACK call, exactly like the original.
//! - `PRINT*` list-directed output is reproduced byte for byte: a leading
//!   blank, the literal verbatim, then the integer right-justified in a
//!   field of 12 columns (gfortran's list-directed integer width).

use std::os::raw::{c_char, c_double, c_int};

extern "C" {
    fn dgetrf_(
        m: *const c_int,
        n: *const c_int,
        a: *mut c_double,
        lda: *const c_int,
        ipiv: *mut c_int,
        info: *mut c_int,
    );
    fn dgetrs_(
        trans: *const c_char,
        n: *const c_int,
        nrhs: *const c_int,
        a: *mut c_double,
        lda: *const c_int,
        ipiv: *mut c_int,
        b: *mut c_double,
        ldb: *const c_int,
        info: *mut c_int,
        trans_len: usize,
    );
    fn dpotrf_(
        uplo: *const c_char,
        n: *const c_int,
        a: *mut c_double,
        lda: *const c_int,
        info: *mut c_int,
        uplo_len: usize,
    );
    fn dpotrs_(
        uplo: *const c_char,
        n: *const c_int,
        nrhs: *const c_int,
        a: *mut c_double,
        lda: *const c_int,
        b: *mut c_double,
        ldb: *const c_int,
        info: *mut c_int,
        uplo_len: usize,
    );
}

/// `NO TRANS` / lower triangle — LAPACK's one-character arguments, passed
/// by reference with gfortran's hidden length argument.
static NO_TRANS: [c_char; 1] = [b'N' as c_char];
static LOWER: [c_char; 1] = [b'L' as c_char];

/// The Fortran's list-directed diagnostic: a leading blank, the literal
/// verbatim, then the integer right-justified in a field of 12 columns.
fn lapack_note(message: &str, info: c_int) {
    println!(" {message}{info:>12}");
}

/// `SOLVELAPACK`: general LU factorization with partial pivoting
/// (`DGETRF`) and the triangular solves (`DGETRS`).
/// # Safety
///
/// `n` must point to a valid `c_int`. `a` must point to a writable
/// `A(N*N)` array and `x` to a writable `X(N)` array, both filled the way
/// the Fortran caller fills them (column-major). The pointers are handed
/// to LAPACK unchanged.
#[export_name = "solvelapack_"]
pub unsafe extern "C" fn port_solve_lapack(
    n: *const c_int,
    a: *mut c_double,
    x: *mut c_double,
) {
    let n = unsafe { *n };
    if n <= 0 {
        return;
    }
    let mut ipiv = vec![0 as c_int; n as usize];
    let mut info: c_int = 0;
    unsafe {
        dgetrf_(&n, &n, a, &n, ipiv.as_mut_ptr(), &mut info);
    }
    if info != 0 {
        lapack_note("DGETRF: ", info);
    }
    info = 0;
    let nrhs: c_int = 1;
    unsafe {
        dgetrs_(
            NO_TRANS.as_ptr(),
            &n,
            &nrhs,
            a,
            &n,
            ipiv.as_mut_ptr(),
            x,
            &n,
            &mut info,
            1,
        );
    }
    if info != 0 {
        lapack_note("DGETRS: ", info);
    }
}

/// `SOLVELAPACKSYM`: symmetric Cholesky factorization (`DPOTRF`, lower
/// triangle) and the triangular solves (`DPOTRS`).
/// # Safety
///
/// `n` must point to a valid `c_int`. `a` must point to a writable
/// `A(N*N)` array and `x` to a writable `X(N)` array, both filled the way
/// the Fortran caller fills them (column-major). The pointers are handed
/// to LAPACK unchanged.
#[export_name = "solvelapacksym_"]
pub unsafe extern "C" fn port_solve_lapack_sym(
    n: *const c_int,
    a: *mut c_double,
    x: *mut c_double,
) {
    let n = unsafe { *n };
    if n <= 0 {
        return;
    }
    let mut info: c_int = 0;
    unsafe {
        dpotrf_(LOWER.as_ptr(), &n, a, &n, &mut info, 1);
    }
    if info != 0 {
        lapack_note("DPOTRF: ", info);
    }
    info = 0;
    let nrhs: c_int = 1;
    unsafe {
        dpotrs_(
            LOWER.as_ptr(),
            &n,
            &nrhs,
            a,
            &n,
            x,
            &n,
            &mut info,
            1,
        );
    }
    if info != 0 {
        lapack_note("DPOTRS: ", info);
    }
}
