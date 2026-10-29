//! `getfrictionheat_`, `getfrictionloads_` — port of
//! USF_GetFrictionHeating.F90 (basal friction heat and friction loads).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated: `getFrictionHeat` keeps `FirstTime`, `DIM`, `Ind`
//!   (the `WRITE(FunctionName,...)` is a dead store — `FunctionName` is
//!   never read — documented, not replicated); `getFrictionLoads` keeps
//!   `FirstTime`, `DIM`, `UseMask` plus a `last` slot for the RESULT.
//! - `frictionLoad` with an active mask (`UseMask`, mask value > 0) is
//!   never assigned: the non-SAVE RESULT local keeps its static storage,
//!   replicated via the state slot (init 0.0).
//! - Fixed-length `CHARACTER(128)` names are replicated blank-padded
//!   (`pad128`, truncate-then-pad), because `(A,A)` message WRITEs embed
//!   the padding; `TRIM(...)` call sites pass the trimmed form, exactly
//!   like the source's `LEN_TRIM`-based lookup.
//! - All order-sensitive reductions (`SUM(a*b)`, `SUM(a**2)`, `**` with a
//!   real exponent) evaluate in the Fortran accessors (`usf_dot`,
//!   `usf_sum_pow2`, `usf_pow_dp`); `SQRT` is same-libm.
//! - The `(A,I3)` WRITE with no integer item is replicated verbatim
//!   through `usf_write_msg_equation_default` (missing-item behavior
//!   included).
//! - `Sig` rows are copied to a contiguous temp before `usf_dot` (copy
//!   is exact; the strided slice never crosses).

use crate::support::*;

fn pad128(s: &str) -> String {
    let b = s.as_bytes();
    let m = b.len().min(128);
    let mut out = String::from_utf8_lossy(&b[..m]).into_owned();
    while out.len() < 128 {
        out.push(' ');
    }
    out
}

struct HeatState {
    first: bool,
    dim: c_int,
    ind: [[c_int; 3]; 3],
}

static HEAT_STATE: Mutex<HeatState> = Mutex::new(HeatState {
    first: true,
    dim: 0,
    ind: [[0; 3]; 3],
});

pub unsafe fn friction_heat(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _dummy = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    {
        let mut st = HEAT_STATE.lock().unwrap();
        if st.first {
            st.dim = unsafe { usf_coord_dim() };
            st.first = false;
            st.ind[0][0] = 1;
            st.ind[1][1] = 2;
            st.ind[2][2] = 3;
            st.ind[0][1] = 4;
            st.ind[1][0] = 4;
            st.ind[2][1] = 5;
            st.ind[1][2] = 5;
            st.ind[0][2] = 6;
            st.ind[2][0] = 6;
        }
    }
    let st = HEAT_STATE.lock().unwrap();
    let (dim, ind) = (st.dim as usize, st.ind);
    drop(st);
    let h_flow = unsafe { variable_get(h_vars, "Flow Solution", true, false) };
    let h_stress = unsafe { variable_get(h_vars, "Stress", true, false) };
    let h_normal = unsafe { variable_get(h_vars, "Normal Vector", true, false) };
    let fp = unsafe { var_perm_at(h_flow, node) };
    let spp = unsafe { var_perm_at(h_stress, node) };
    let np = unsafe { var_perm_at(h_normal, node) };
    let mut normal = [0.0 as c_double; 3];
    let mut velo = [0.0 as c_double; 3];
    for i in 0..dim {
        // SAFETY: same 1-based slot arithmetic the source uses.
        normal[i] = unsafe { var_value_at_idx(h_normal, dim as c_int * (np - 1) + 1 + i as c_int) };
        velo[i] = unsafe {
            var_value_at_idx(h_flow, (dim as c_int + 1) * (fp - 1) + 1 + i as c_int)
        };
    }
    let un = unsafe { dot(velo.as_ptr(), normal.as_ptr(), dim) };
    let ut = unsafe { sum_pow2(velo.as_ptr(), dim) - usf_pow_dp(un, 2.0) }.sqrt();
    let mut sig = [0.0 as c_double; 9];
    for i in 0..dim {
        for j in 0..dim {
            let slot = 2 * dim as c_int * (spp - 1) + ind[i][j];
            sig[j * 3 + i] = unsafe { var_value_at_idx(h_stress, slot) };
        }
    }
    let mut sn = [0.0 as c_double; 3];
    for i in 0..dim {
        let mut row = [0.0 as c_double; 3];
        for j in 0..dim {
            row[j] = sig[j * 3 + i];
        }
        sn[i] = unsafe { dot(row.as_ptr(), normal.as_ptr(), dim) };
    }
    let snn = unsafe { dot(sn.as_ptr(), normal.as_ptr(), dim) };
    let snt = unsafe { sum_pow2(sn.as_ptr(), dim) - usf_pow_dp(snn, 2.0) }.sqrt();
    ut * snt
}

struct LoadState {
    first: bool,
    dim: c_int,
    use_mask: bool,
    last: c_double,
}

static LOAD_STATE: Mutex<LoadState> = Mutex::new(LoadState {
    first: true,
    dim: 0,
    use_mask: false,
    last: 0.0,
});

const FUNCTION_NAME: &str = "USF_GetFrictionHeating(getFrictionLoads)";

pub unsafe fn friction_loads(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _dummy = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let mut warned = false;
    let (flow_name, loads_name, mask_name) = {
        let mut st = LOAD_STATE.lock().unwrap();
        if st.first {
            st.dim = unsafe { usf_coord_dim() };
        }
        // Boundary-element parent resolution (outbody/Right/Left).
        let mut h_elem: *const c_void = std::ptr::null();
        unsafe { usf_current_element(&mut h_elem) };
        let h_parent = unsafe { boundary_parent(h_elem) };
        let mut got = false;
        let h_eq = unsafe { equation_of(h_parent, &mut got) };
        let first = st.first;
        let (mut fs, mut fl, mn);
        if h_eq.is_null() || !got {
            if first {
                unsafe { warn(FUNCTION_NAME, &write_msg_equation_default()) };
                warned = true;
            }
            fs = pad128("Flow Solution");
            fl = pad128("Flow Solution Loads");
            st.use_mask = false;
            mn = pad128("");
        } else {
            let mut g = false;
            fs = pad128(&unsafe { get_string(h_eq, "Flow Solution Name", &mut g) });
            if !g {
                // `WRITE(FlowSolutionName,'(A)') 'Flow Solution'`.
                let d = pad128("Flow Solution");
                if first {
                    let msg = format!("Using default name for flow solution: {}", d);
                    unsafe { warn(FUNCTION_NAME, &msg) };
                    warned = true;
                }
                fs = d;
            }
            // N.B. the source reassigns unconditionally here.
            let mut g2 = false;
            fl = pad128(&unsafe { get_string(h_eq, "Flow Loads Name", &mut g2) });
            if !g2 {
                let d = pad128(&format!("{} Loads", fs.trim_end()));
                if first {
                    let msg = format!("Using default name for flow solution loads: {}", d);
                    unsafe { warn(FUNCTION_NAME, &msg) };
                    warned = true;
                }
                fl = d;
            }
            let mut um = false;
            mn = pad128(&unsafe { get_string(h_eq, "Friction Load Mask", &mut um) });
            st.use_mask = um;
            if um {
                let msg = format!(
                    ">Friction Load Mask< found and set to {}",
                    mn.trim_end()
                );
                unsafe { info(FUNCTION_NAME, &msg, 20) };
            }
        }
        if warned && first {
            unsafe { warn(FUNCTION_NAME, "All Warnings will be further omitted") };
        }
        st.first = false;
        (fs, fl, mn)
    };
    let st = LOAD_STATE.lock().unwrap();
    let (dim, use_mask) = (st.dim as usize, st.use_mask);
    drop(st);
    let h_flow = unsafe { variable_get(h_vars, flow_name.trim_end(), true, false) };
    let h_fload = unsafe { variable_get(h_vars, loads_name.trim_end(), true, false) };
    let h_normal = unsafe { variable_get(h_vars, "Normal Vector", true, false) };
    if use_mask {
        let h_mask = unsafe { variable_get(h_vars, mask_name.trim_end(), true, false) };
        if unsafe { var_value_at(h_mask, node) } <= 0.0 {
            LOAD_STATE.lock().unwrap().last = 0.0;
            return 0.0;
        }
        // Mask value > 0: the source assigns nothing — stale RESULT.
        return LOAD_STATE.lock().unwrap().last;
    }
    let fp = unsafe { var_perm_at(h_flow, node) };
    let lp = unsafe { var_perm_at(h_fload, node) };
    let np = unsafe { var_perm_at(h_normal, node) };
    let mut normal = [0.0 as c_double; 3];
    let mut velo = [0.0 as c_double; 3];
    let mut fload = [0.0 as c_double; 3];
    for i in 0..dim {
        normal[i] = unsafe { var_value_at_idx(h_normal, dim as c_int * (np - 1) + 1 + i as c_int) };
        velo[i] = unsafe {
            var_value_at_idx(h_flow, (dim as c_int + 1) * (fp - 1) + 1 + i as c_int)
        };
        fload[i] = unsafe {
            var_value_at_idx(h_fload, (dim as c_int + 1) * (lp - 1) + 1 + i as c_int)
        };
    }
    let normalvelocity = unsafe { dot(velo.as_ptr(), normal.as_ptr(), dim) };
    let mut tang = [0.0 as c_double; 3];
    for i in 0..dim {
        tang[i] = velo[i] - normalvelocity * normal[i];
    }
    // MAX((-1.0 * SUM(tang*flowLoad)), 0.0) — NaN yields 0.0 in both arms.
    let val = -unsafe { dot(tang.as_ptr(), fload.as_ptr(), dim) };
    let load = if val > 0.0 { val } else { 0.0 };
    LOAD_STATE.lock().unwrap().last = load;
    load
}
