//! Port of `fem/src/SolveSBand.F90` (Elmer @ 9f6af2f85) for lane-fem-a,
//! through the gfortran ABI shim (ADR-027, T017 shape A).
//!
//! Unit: `fortran:fem/src/SolveSBand.F90` — two file-scope external
//! subprograms. The emitted gfortran symbols are `solvesbandlapack_` and
//! `solvecomplexsbandlapack_`; every dummy travels by reference (a
//! pointer to the caller's argument).
//!
//! The unit is a LAPACK symmetric band-solver wrapper: the factorization
//! and the triangular solves stay in LAPACK (`dpbtrf_`/`dpbtrs_`,
//! `zpbtrf_`/`zpbtrs_`), so the numerics are bit-identical by
//! construction. What the port owns is the wrapper itself — the argument
//! marshalling, the `N <= 0` early return, the `info` check and the
//! failure path.
//!
//! Faithfulness notes:
//! - The Fortran declares `A(Band,N)` and `X(N,M)` and passes `LDB = N`
//!   to LAPACK; the port passes the same values so every shape behaves
//!   identically.
//! - `PRINT*` list-directed output is reproduced byte for byte: a
//!   leading blank, the literal verbatim, then the integer right-justified
//!   in a field of 12 columns (gfortran's list-directed integer width).
//! - `STOP 1` writes `STOP 1` to stderr and exits with status 1.

use std::io::Write;
use std::os::raw::{c_char, c_double, c_int};

/// gfortran's `COMPLEX(KIND=dp)`: two 8-byte reals, real part first.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct C64 {
    re: c_double,
    im: c_double,
}

extern "C" {
    fn dpbtrf_(
        uplo: *const c_char,
        n: *const c_int,
        kd: *const c_int,
        ab: *mut c_double,
        ldab: *const c_int,
        info: *mut c_int,
        uplo_len: usize,
    );
    fn dpbtrs_(
        uplo: *const c_char,
        n: *const c_int,
        kd: *const c_int,
        nrhs: *const c_int,
        ab: *mut c_double,
        ldab: *const c_int,
        b: *mut c_double,
        ldb: *const c_int,
        info: *mut c_int,
        uplo_len: usize,
    );
    fn zpbtrf_(
        uplo: *const c_char,
        n: *const c_int,
        kd: *const c_int,
        ab: *mut C64,
        ldab: *const c_int,
        info: *mut c_int,
        uplo_len: usize,
    );
    fn zpbtrs_(
        uplo: *const c_char,
        n: *const c_int,
        kd: *const c_int,
        nrhs: *const c_int,
        ab: *mut C64,
        ldab: *const c_int,
        b: *mut C64,
        ldb: *const c_int,
        info: *mut c_int,
        uplo_len: usize,
    );
}

/// `LOWER` — LAPACK's `UPLO` argument, one character, passed by reference
/// with gfortran's hidden length argument.
static LOWER: [c_char; 1] = [b'L' as c_char];

/// The failure path of the Fortran: the list-directed `PRINT*` line on
/// stdout, then `STOP 1` on stderr with exit status 1.
fn sband_failure(message: &str, info: c_int) -> ! {
    println!(" {message}{info:>12}");
    let _ = std::io::stdout().flush();
    eprintln!("STOP 1");
    std::process::exit(1);
}

/// `SOLVESBANDLAPACK`: symmetric band Cholesky factorization (`DPBTRF`)
/// and the band triangular solves (`DPBTRS`).
#[export_name = "solvesbandlapack_"]
/// # Safety
///
/// `n`, `m`, `subband` and `band` must point to valid `c_int`s. `a` must
/// point to a writable `A(Band,N)` array and `x` to a writable `X(N,M)`
/// array, both filled the way the Fortran caller fills them (LAPACK
/// symmetric band storage). The pointers are handed to LAPACK unchanged.
pub unsafe extern "C" fn port_solve_sband_lapack(
    n: *const c_int,
    m: *const c_int,
    a: *mut c_double,
    x: *mut c_double,
    subband: *const c_int,
    band: *const c_int,
) {
    let (n, m, subband, band) = unsafe { (*n, *m, *subband, *band) };
    if n <= 0 {
        return;
    }
    let mut info: c_int = 0;
    unsafe {
        dpbtrf_(
            LOWER.as_ptr(),
            &n,
            &subband,
            a,
            &band,
            &mut info,
            1,
        );
    }
    if info != 0 {
        sband_failure(
            "ERROR: SolveSymmetricBand: singular matrix. LAPACK DPBTRF iinfo: ",
            info,
        );
    }
    info = 0;
    unsafe {
        dpbtrs_(
            LOWER.as_ptr(),
            &n,
            &subband,
            &m,
            a,
            &band,
            x,
            &n,
            &mut info,
            1,
        );
    }
    if info != 0 {
        sband_failure(
            "ERROR: SolveSymmetricBand: singular matrix. LAPACK DPBTRS iinfo: ",
            info,
        );
    }
}

/// `SOLVECOMPLEXSBANDLAPACK`: the complex counterpart
/// (`ZPBTRF`/`ZPBTRS`).
#[export_name = "solvecomplexsbandlapack_"]
/// # Safety
///
/// `n`, `m`, `subband` and `band` must point to valid `c_int`s. `a` must
/// point to a writable complex `A(Band,N)` array and `x` to a writable
/// complex `X(N,M)` array, both filled the way the Fortran caller fills
/// them (LAPACK symmetric band storage). The pointers are handed to
/// LAPACK unchanged.
pub unsafe extern "C" fn port_solve_complex_sband_lapack(
    n: *const c_int,
    m: *const c_int,
    a: *mut C64,
    x: *mut C64,
    subband: *const c_int,
    band: *const c_int,
) {
    let (n, m, subband, band) = unsafe { (*n, *m, *subband, *band) };
    if n <= 0 {
        return;
    }
    let mut info: c_int = 0;
    unsafe {
        zpbtrf_(
            LOWER.as_ptr(),
            &n,
            &subband,
            a,
            &band,
            &mut info,
            1,
        );
    }
    if info != 0 {
        sband_failure(
            "ERROR: SolveSymmetricBand: singular matrix. LAPACK ZPBTRF iinfo: ",
            info,
        );
    }
    info = 0;
    unsafe {
        zpbtrs_(
            LOWER.as_ptr(),
            &n,
            &subband,
            &m,
            a,
            &band,
            x,
            &n,
            &mut info,
            1,
        );
    }
    if info != 0 {
        sband_failure(
            "ERROR: SolveSymmetricBand: singular matrix. LAPACK ZPBTRS iinfo: ",
            info,
        );
    }
}
