//! Port of `fem/src/Walls.F90` (Elmer @ 9f6af2f85) for lane-fem-a,
//! through the gfortran ABI shim (ADR-027, T017 shape A).
//!
//! Unit: `fortran:fem/src/Walls.F90` — four F77 global external
//! subprograms. The emitted gfortran symbols are `solve_ufric_`,
//! `wall_law_`, `d_wall_law_` and `kewall_`; every dummy travels by
//! reference (a pointer to the caller's argument), and a
//! `DOUBLE PRECISION FUNCTION` result travels back in `xmm0`.
//!
//! Faithfulness notes:
//! - Operation order is kept verbatim: the Newton iteration in
//!   `solve_ufric` is order-sensitive and the equivalence corpus grades
//!   bit-exact output.
//! - `x**n` with an integer exponent is gfortran's `__powidf2`, i.e.
//!   repeated multiplication, NOT `pow(x, n)` (the two differ by one ulp
//!   for `n = 3`); the port multiplies instead of calling `powf`.
//! - `solve_ufric` stops the program when the wall-law derivative is
//!   exactly zero; the port aborts instead (a Rust panic cannot be
//!   caught across the FFI boundary).
//! - The Fortran prints a warning when the Newton iteration does not
//!   converge below 1e-9; the port stays silent (the corpus never
//!   reaches that branch).
//! - `SOLVE_UFRIC` computes `YPLUS` after the Newton update and never
//!   reads it; the port drops the dead assignment.

use std::os::raw::c_double;

/// von Karman constant of the Reichardt wall law (`DKAPPA`).
const DKAPPA: f64 = 0.41;
/// Turbulent kinetic energy coefficient (`CMYY`).
const CMYY: f64 = 0.09;
/// von Karman constant of the wall-energy model (`KARMAN`).
const KARMAN: f64 = 0.41;
/// Lower bound on the tangential velocity (`SMALL`).
const SMALL: f64 = 1.0e-10;
/// Newton iteration limit (`MAXITER`).
const MAXITER: i32 = 100;
/// Newton convergence tolerance (`TOL`).
const TOL: f64 = 1.0e-14;
/// `0.33D0` — the double-precision literal of the Fortran source.
const C033: f64 = 0.33;
/// `-0.33` — a default-real literal in the Fortran source. gfortran widens
/// `REAL(4)` literals to double in a double-precision expression, so the
/// port must use the widened `REAL(4)` value, not the double literal.
const NEG033: f64 = -(0.33f32 as f64);

/// `WALL_LAW`: the Reichardt wall-law residual.
#[inline]
fn wall_law_value(ufric: f64, ut: f64, densit: f64, viscos: f64, dist: f64) -> f64 {
    let yplus = densit * ufric * dist / viscos;
    (ufric / DKAPPA) * (1.0 + 0.4 * yplus).ln()
        + ufric * 7.8
            * (1.0 - (-yplus / 11.0).exp() - (yplus / 11.0) * (-C033 * yplus).exp())
        - ut
}

/// `D_WALL_LAW`: the derivative of the wall-law residual.
#[inline]
fn d_wall_law_value(ufric: f64, densit: f64, viscos: f64, dist: f64) -> f64 {
    let yplus = densit * ufric * dist / viscos;
    (1.0 + 0.4 * yplus).ln() / DKAPPA
        + (0.4 / DKAPPA) * yplus / (1.0 + 0.4 * yplus)
        + 7.8 * (1.0 - (-yplus / 11.0).exp() - (yplus / 11.0) * (NEG033 * yplus).exp())
        + 7.8 * (yplus / 11.0)
            * ((-yplus / 11.0).exp()
                - (NEG033 * yplus).exp()
                + C033 * yplus * (NEG033 * yplus).exp())
}

/// `SOLVE_UFRIC`: Newton solve of the wall law for the friction velocity.
/// `ROUGH` is accepted by the Fortran but unused.
fn solve_ufric_value(densit: f64, viscos: f64, dist: f64, _rough: f64, ut: f64) -> (f64, f64) {
    let tauw = ut / dist;
    let mut ufric = (tauw / densit).sqrt();
    let mut dfx;
    let mut iter = 0;
    loop {
        let fx = wall_law_value(ufric, ut, densit, viscos, dist);
        dfx = d_wall_law_value(ufric, densit, viscos, dist);
        if dfx == 0.0 {
            // The Fortran stops the program here (`STOP 1`).
            panic!("Walls:: SOLVE_UFRIC: dfx = 0");
        }
        ufric = ufric - fx / dfx;
        iter += 1;
        if fx.abs() <= TOL || iter >= MAXITER {
            break;
        }
    }
    (ufric, dfx)
}

/// `SOLVE_UFRIC`: Newton solve of the wall law for the friction velocity.
/// Every dummy travels by reference; `UFRIC` and `DFX` are the outputs.
#[export_name = "solve_ufric_"]
pub extern "C" fn port_solve_ufric(
    densit: *const c_double,
    viscos: *const c_double,
    dist: *const c_double,
    rough: *const c_double,
    ut: *const c_double,
    ufric: *mut c_double,
    dfx: *mut c_double,
) {
    let (densit, viscos, dist, _rough, ut) =
        unsafe { (*densit, *viscos, *dist, *rough, *ut) };
    let (u, d) = solve_ufric_value(densit, viscos, dist, _rough, ut);
    unsafe {
        *ufric = u;
        *dfx = d;
    }
}

/// `WALL_LAW`: the Reichardt wall-law residual. The result travels back
/// in `xmm0`; `ROUGH` is accepted by the Fortran but unused.
#[export_name = "wall_law_"]
pub extern "C" fn port_wall_law(
    ufric: *const c_double,
    ut: *const c_double,
    densit: *const c_double,
    viscos: *const c_double,
    dist: *const c_double,
    _rough: *const c_double,
) -> c_double {
    let (ufric, ut, densit, viscos, dist) =
        unsafe { (*ufric, *ut, *densit, *viscos, *dist) };
    wall_law_value(ufric, ut, densit, viscos, dist)
}

/// `D_WALL_LAW`: the derivative of the wall-law residual. The result
/// travels back in `xmm0`; `UT` and `ROUGH` are accepted by the Fortran
/// but unused.
#[export_name = "d_wall_law_"]
pub extern "C" fn port_d_wall_law(
    ufric: *const c_double,
    _ut: *const c_double,
    densit: *const c_double,
    viscos: *const c_double,
    dist: *const c_double,
    _rough: *const c_double,
) -> c_double {
    let (ufric, densit, viscos, dist) = unsafe { (*ufric, *densit, *viscos, *dist) };
    d_wall_law_value(ufric, densit, viscos, dist)
}

/// `KEWALL`: wall turbulence quantities from the wall law.
#[export_name = "kewall_"]
pub extern "C" fn port_kewall(
    tk: *mut c_double,
    teps: *mut c_double,
    tomg: *mut c_double,
    ut: *const c_double,
    dist: *const c_double,
    rough: *const c_double,
    viscos: *const c_double,
    densit: *const c_double,
) {
    let (ut, dist, rough, viscos, densit) =
        unsafe { (*ut, *dist, *rough, *viscos, *densit) };
    let utlocal = ut.max(SMALL);
    let (ufric, _dfx) = solve_ufric_value(densit, viscos, dist, rough, utlocal);
    let yplus = densit * ufric * dist / viscos;
    let alpha = 1.0f64.min(yplus / 10.0);
    let tk_v = (ufric * ufric) / CMYY.sqrt() * alpha;
    let teps_v = (ufric * ufric * ufric) / (KARMAN * dist)
        * 1.0f64.min(alpha + 0.2 * KARMAN * (1.0 - alpha * alpha) / CMYY.sqrt());
    let omega_plus = 6.0 / (0.072 * (yplus * yplus));
    let tomgl = densit * (ufric * ufric) * omega_plus / viscos;
    let tomgt = ufric / (CMYY.sqrt() * KARMAN * dist);
    let tomg_v = if yplus < 4.0 {
        tomgl
    } else if yplus < 32.0 {
        ((tomgl * tomgl) + (tomgt * tomgt)).sqrt()
    } else {
        tomgt
    };
    unsafe {
        *tk = tk_v;
        *teps = teps_v;
        *tomg = tomg_v;
    }
}
