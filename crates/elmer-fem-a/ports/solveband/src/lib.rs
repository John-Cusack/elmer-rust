//! Port of `fem/src/SolveBand.F90` (Elmer @ 9f6af2f85) for lane-fem-a,
//! through the gfortran ABI shim (ADR-027, T017 shape A).
//!
//! Unit: `fortran:fem/src/SolveBand.F90` — two file-scope external
//! subprograms. The emitted gfortran symbols are `solvebandlapack_` and
//! `solvecomplexbandlapack_`; every dummy travels by reference (a pointer
//! to the caller's argument).
//!
//! The unit is a LAPACK band-solver wrapper: the factorization and the
//! triangular solves stay in LAPACK (`dgbtrf_`/`dgbtrs_`,
//! `zgbtrf_`/`zgbtrs_`), so the numerics are bit-identical by
//! construction. What the port owns is the wrapper itself — the argument
//! marshalling, the `N <= 0` early return, the `info` check and the
//! failure path.
//!
//! Faithfulness notes:
//! - The Fortran declares `A(Band,N)` and `X(M,N)` and passes `LDB = N`
//!   to LAPACK. With `NRHS = 1` (the only shape the campaign's callers
//!   use) the leading dimension does not affect the access pattern; the
//!   port passes the same `LDB = N` so every other shape behaves
//!   identically too.
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
    fn dgbtrf_(
        m: *const c_int,
        n: *const c_int,
        kl: *const c_int,
        ku: *const c_int,
        ab: *mut c_double,
        ldab: *const c_int,
        ipiv: *mut c_int,
        info: *mut c_int,
    );
    fn dgbtrs_(
        trans: *const c_char,
        n: *const c_int,
        kl: *const c_int,
        ku: *const c_int,
        nrhs: *const c_int,
        ab: *mut c_double,
        ldab: *const c_int,
        ipiv: *mut c_int,
        b: *mut c_double,
        ldb: *const c_int,
        info: *mut c_int,
        trans_len: usize,
    );
    fn zgbtrf_(
        m: *const c_int,
        n: *const c_int,
        kl: *const c_int,
        ku: *const c_int,
        ab: *mut C64,
        ldab: *const c_int,
        ipiv: *mut c_int,
        info: *mut c_int,
    );
    fn zgbtrs_(
        trans: *const c_char,
        n: *const c_int,
        kl: *const c_int,
        ku: *const c_int,
        nrhs: *const c_int,
        ab: *mut C64,
        ldab: *const c_int,
        ipiv: *mut c_int,
        b: *mut C64,
        ldb: *const c_int,
        info: *mut c_int,
        trans_len: usize,
    );
}

/// `NO TRANS` — LAPACK's `TRANS` argument, one character, passed by
/// reference with gfortran's hidden length argument.
static NO_TRANS: [c_char; 1] = [b'N' as c_char];

/// The failure path of the Fortran: the list-directed `PRINT*` line on
/// stdout, then `STOP 1` on stderr with exit status 1.
fn band_failure(message: &str, info: c_int) -> ! {
    println!(" {message}{info:>12}");
    let _ = std::io::stdout().flush();
    eprintln!("STOP 1");
    std::process::exit(1);
}

/// `SOLVEBANDLAPACK`: general band LU factorization with partial pivoting
/// (`DGBTRF`) and the band triangular solves (`DGBTRS`).
/// # Safety
///
/// `n`, `m`, `subband` and `band` must point to valid `c_int`s. `a` must
/// point to a writable `A(Band,N)` array and `x` to a writable `X(M,N)`
/// array, both filled the way the Fortran caller fills them (LAPACK band
/// storage). The pointers are handed to LAPACK unchanged.
#[export_name = "solvebandlapack_"]
pub unsafe extern "C" fn port_solve_band_lapack(
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
    let mut ipiv = vec![0 as c_int; n as usize];
    let mut info: c_int = 0;
    unsafe {
        dgbtrf_(
            &n,
            &n,
            &subband,
            &subband,
            a,
            &band,
            ipiv.as_mut_ptr(),
            &mut info,
        );
    }
    if info != 0 {
        band_failure("ERROR: SolveBand: singular matrix. LAPACK DGBTRF info: ", info);
    }
    info = 0;
    unsafe {
        dgbtrs_(
            NO_TRANS.as_ptr(),
            &n,
            &subband,
            &subband,
            &m,
            a,
            &band,
            ipiv.as_mut_ptr(),
            x,
            &n,
            &mut info,
            1,
        );
    }
    if info != 0 {
        band_failure("ERROR: SolveBand: singular matrix. LAPACK DGBTRS info: ", info);
    }
}

/// `SOLVECOMPLEXBANDLAPACK`: the complex counterpart (`ZGBTRF`/`ZGBTRS`).
/// # Safety
///
/// `n`, `m`, `subband` and `band` must point to valid `c_int`s. `a` must
/// point to a writable complex `A(Band,N)` array and `x` to a writable
/// complex `X(M,N)` array, both filled the way the Fortran caller fills
/// them (LAPACK band storage). The pointers are handed to LAPACK
/// unchanged.
#[export_name = "solvecomplexbandlapack_"]
pub unsafe extern "C" fn port_solve_complex_band_lapack(
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
    let mut ipiv = vec![0 as c_int; n as usize];
    let mut info: c_int = 0;
    unsafe {
        zgbtrf_(
            &n,
            &n,
            &subband,
            &subband,
            a,
            &band,
            ipiv.as_mut_ptr(),
            &mut info,
        );
    }
    if info != 0 {
        band_failure(
            "ERROR: SolveBand: singular matrix. LAPACK ZGBTRF info: ",
            info,
        );
    }
    info = 0;
    unsafe {
        zgbtrs_(
            NO_TRANS.as_ptr(),
            &n,
            &subband,
            &subband,
            &m,
            a,
            &band,
            ipiv.as_mut_ptr(),
            x,
            &n,
            &mut info,
            1,
        );
    }
    if info != 0 {
        band_failure(
            "ERROR: SolveBand: singular matrix. LAPACK ZGBTRS info: ",
            info,
        );
    }
}
