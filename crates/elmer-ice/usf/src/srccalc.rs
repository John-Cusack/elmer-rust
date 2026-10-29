//! `sourcecalc_` — port of USF_SourceCalcCalving.F90 (calving-front
//! melt source from internal/surface melt variables).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated: `FirstTime`, `HPSolver` (0 = never matched, passed
//!   straight through — the reference faults identically on `Solvers(0)`),
//!   `Calving`.
//! - `GetLogical(Model % Simulation, ...)` reads the simulation list;
//!   the hydraulic-potential solver search compares `==` case-sensitively
//!   exactly like the source.
//! - `Found`/`SMLogical`/`IMLogical` are plain locals reset every call;
//!   `IMVarName`/`SMVarName` are set before every read, like the source.
//! - NaN handling uses `is_nan` (IEEE, identical to `ISNAN`).

use crate::support::*;

struct SrcState {
    first: bool,
    hp_solver: c_int,
    calving: bool,
}

static STATE: Mutex<SrcState> = Mutex::new(SrcState {
    first: true,
    hp_solver: 0,
    calving: false,
});

pub unsafe fn source_calc(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _x = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_bf: *const c_void = std::ptr::null();
    unsafe { usf_get_bodyforce_of(h_elem, &mut h_bf) };
    {
        let mut st = STATE.lock().unwrap();
        if st.first {
            let h_sim = unsafe { simulation_values() };
            let mut found = false;
            let calv = unsafe { get_logical(h_sim, "Calving", &mut found) };
            st.calving = if found { calv } else { false };
            let nsolvers = unsafe { usf_number_of_solvers() };
            let mut i = 1;
            while i <= nsolvers {
                // SAFETY: same call sequence the reference makes; a
                // null solver Variable faults identically in both arms.
                if unsafe { solver_variable_name(i) } == "hydraulic potential" {
                    st.hp_solver = i;
                    break;
                }
                i += 1;
            }
            unsafe { calculate_nodal_weights(st.hp_solver, "Weights") };
            st.first = false;
        }
    }
    let st = STATE.lock().unwrap();
    let (hp, calving) = (st.hp_solver, st.calving);
    drop(st);
    let mut found = false;
    let im_logical = unsafe { get_logical(h_bf, "Internal Melt", &mut found) };
    if !found {
        unsafe { fatal("SourceCalc", "You've not specified if you want internal melt, silly") };
    }
    let (h_imvar, h_weights) = if im_logical {
        let mut got = false;
        let im_name = unsafe { get_string(h_bf, "Internal Melt Variable Name", &mut got) };
        let h_vars = if calving {
            unsafe { solver_mesh_variables(hp) }
        } else {
            let mut h: *const c_void = std::ptr::null();
            unsafe { usf_model_variables(&mut h) };
            h
        };
        let hi = unsafe { variable_get(h_vars, &im_name, true, true) };
        let hw = unsafe { variable_get(h_vars, "Weights", true, true) };
        (hi, hw)
    } else {
        (std::ptr::null(), std::ptr::null())
    };
    let sm_name = unsafe { get_string(h_bf, "Surface Melt Variable Name", &mut found) };
    let sm_logical = unsafe { get_logical(h_bf, "Surface Melt", &mut found) };
    if !found {
        unsafe { fatal("SourceCalc", "You've not defined surface melt, you muppet") };
    }
    let h_smvar = if sm_logical {
        let h_vars = if calving {
            unsafe { solver_mesh_variables(hp) }
        } else {
            let mut h: *const c_void = std::ptr::null();
            unsafe { usf_model_variables(&mut h) };
            h
        };
        unsafe { variable_get(h_vars, &sm_name, true, true) }
    } else {
        std::ptr::null()
    };
    let internal_melt = if im_logical {
        let p = unsafe { var_perm_at(h_imvar, node) };
        if p > 0 {
            let val = unsafe { var_value_at_idx(h_imvar, p) };
            // Cases: val>0 -> 0; val==0 -> 0 (source first arm);
            // val<0 -> compute (third arm); NaN -> 0 (second arm, since
            // NaN>=0.0 is false). Collapsed without changing any value.
            if val >= 0.0 || val.is_nan() {
                0.0
            } else {
                let w = unsafe { var_value_at_idx(h_weights, var_perm_at(h_weights, node)) };
                (val.abs() / w) / 333.55
            }
        } else {
            0.0
        }
    } else {
        0.0
    };
    let surface_melt = if sm_logical {
        // SAFETY: `SMVar % Perm(NodeNumber)` may be disassociated only
        // when the source would fault too; pass-through.
        let p = unsafe { var_perm_at(h_smvar, node) };
        let val = unsafe { var_value_at_idx(h_smvar, p) };
        if calving {
            // NaN>0.0 is false, so the source's first two arms both yield
            // 0.0; collapsed without changing any value.
            if val > 0.0 { val } else { 0.0 }
        } else if val > 0.0 {
            val
        } else {
            0.0
        }
    } else {
        0.0
    };
    internal_melt + surface_melt
}
