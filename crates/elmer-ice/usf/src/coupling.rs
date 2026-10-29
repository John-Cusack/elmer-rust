//! `horizontalvelo_`, `overburdenpressure_` — port of
//! USF_CouplingGlaDS_SSA.F90 (GlaDS/SSA coupling: horizontal velocity and
//! overburden pressure).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - No SAVE state in either entry (all locals set before every read).
//! - `ub = SQRT(velo(1)**2+velo(2)**2)` always reads indices 1..2 even
//!   when `DIM == 1` (`velo` is zero-filled first) — verbatim.
//! - `print *, "Ice pressure value :", IcePress` on negative pressure
//!   goes through the Fortran accessor (list-directed formatting +
//!   stdout, byte-identical).

use crate::support::*;

pub unsafe fn horizontal_velo(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _x = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let dim = unsafe { usf_coord_dim() } as usize;
    let h_flow = unsafe { variable_get(h_vars, "ssavelocity", false, false) };
    if h_flow.is_null() {
        unsafe { fatal("HorizontalVelo", "Need SSA Solver, variable >SSAVelocity< not associated !!!") };
    }
    let fp = unsafe { var_perm_at(h_flow, node) };
    let mut velo = [0.0 as c_double; 2];
    // 2D-sized scratch, filled `DIM` deep like the source (2D-only entries).
    for (i, v) in velo.iter_mut().enumerate().take(dim) {
        *v = unsafe { var_value_at_idx(h_flow, dim as c_int * (fp - 1) + 1 + i as c_int) };
    }
    // N.B. integer `**2` is repeated multiplication, not `pow` — the
    // reference evaluates `velo(1)*velo(1) + velo(2)*velo(2)` exactly.
    (velo[0] * velo[0] + velo[1] * velo[1]).sqrt()
}

pub unsafe fn overburden_pressure(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _x = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let _dim = unsafe { usf_coord_dim() };
    let mut h_c: *const c_void = std::ptr::null();
    unsafe { usf_get_constants(&mut h_c) };
    let mut got = false;
    let gravity = unsafe { const_real(h_c, "Gravity Norm", &mut got) };
    if !got {
        unsafe { fatal("OverburdenPressure", "Keyword >Gravity Norm< not found in Constant section") };
    }
    let rhoi = unsafe { const_real(h_c, "Ice Density", &mut got) };
    if !got {
        unsafe { fatal("OverburdenPressure", "Keyword >Ice Density< not found in Constant section") };
    }
    let h_zb = unsafe { variable_get(h_vars, "zb", false, false) };
    if h_zb.is_null() {
        unsafe { fatal("OverburdenPressure", "Variable >Zb< not associated !!!") };
    }
    let h_zs = unsafe { variable_get(h_vars, "zs", false, false) };
    if h_zs.is_null() {
        unsafe { fatal("OverburdenPressure", "Variable >Zs< not associated !!!") };
    }
    let zb = unsafe { var_value_at(h_zb, node) };
    let zs = unsafe { var_value_at(h_zs, node) };
    let ice_press = rhoi * gravity * (zs - zb);
    if ice_press < 0.0 {
        unsafe { print_ice_pressure(ice_press) };
    }
    ice_press
}
