//! Port of `fem/src/solve_cmplx.F90` (Elmer @ c0b4bc6, the campaign fork)
//! for lane T020, through the gfortran ABI shim (ADR-027, T017 shape A,
//! T020 `KIND=` legacy-parse fix).
//!
//! Unit: `fortran:fem/src/solve_cmplx.F90` — one file-scope external
//! subprogram. The emitted gfortran symbol is `solvelapack_cmplx_`; every
//! dummy travels by reference (a pointer to the caller's argument).
//!
//! The unit is the complex twin of the ported `solve_real` wrapper: the
//! LU factorization and the triangular solves stay in LAPACK (`zgetrf_` /
//! `zgetrs_`), so the numerics are bit-identical by construction. What
//! the port owns is the wrapper itself — the argument marshalling, the
//! `N <= 0` early return and the `info` diagnostic prints.
//!
//! Faithfulness notes:
//! - The Fortran declares `COMPLEX(KIND=dp) A(N*N)` and `X(N)` and passes
//!   `LDA = N` to LAPACK. The port passes the same `LDA = N`, so every
//!   shape behaves identically, not just the square systems the
//!   campaign's callers use.
//! - A nonzero `info` only prints: the Fortran has no `STOP` here, so the
//!   port prints the same list-directed line and continues into the next
//!   LAPACK call, exactly like the original.
//! - `PRINT*` list-directed output is reproduced byte for byte: a leading
//!   blank, the literal verbatim, then the integer right-justified in a
//!   field of 12 columns (gfortran's list-directed integer width — the
//!   same convention the `solvereal` port pins).
//! - The Fortran's `COMPLEX(KIND=dp)` dummies are `complex*16`: pairs of
//!   `f64` in memory. The port hands the caller's pointers to LAPACK
//!   unchanged and never interprets the lanes itself.

use std::os::raw::{c_char, c_double, c_int};

extern "C" {
    fn zgetrf_(
        m: *const c_int,
        n: *const c_int,
        a: *mut c_double,
        lda: *const c_int,
        ipiv: *mut c_int,
        info: *mut c_int,
    );
    fn zgetrs_(
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
}

/// `NO TRANS` — LAPACK's one-character argument, passed by reference with
/// gfortran's hidden length argument.
static NO_TRANS: [c_char; 1] = [b'N' as c_char];

/// The Fortran's list-directed diagnostic: a leading blank, the literal
/// verbatim, then the integer right-justified in a field of 12 columns.
fn lapack_note(message: &str, info: c_int) {
    println!(" {message}{info:>12}");
}

/// `SOLVELAPACK_CMPLX`: general LU factorization with partial pivoting
/// (`ZGETRF`) and the triangular solves (`ZGETRS`) on `complex*16` data.
/// # Safety
///
/// `n` must point to a valid `c_int`. `a` must point to a writable
/// `A(N*N)` complex array (2·N·N `f64` lanes) and `x` to a writable `X(N)`
/// complex array, both filled the way the Fortran caller fills them
/// (column-major). The pointers are handed to LAPACK unchanged.
#[export_name = "solvelapack_cmplx_"]
pub unsafe extern "C" fn port_solve_lapack_cmplx(
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
        zgetrf_(&n, &n, a, &n, ipiv.as_mut_ptr(), &mut info);
    }
    if info != 0 {
        lapack_note("ZGETRF: ", info);
    }
    info = 0;
    let nrhs: c_int = 1;
    unsafe {
        zgetrs_(
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
        lapack_note("ZGETRS: ", info);
    }
}
