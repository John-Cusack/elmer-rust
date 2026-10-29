//! `lateralfriction_x_`, `lateralfriction_y_`, `lateralfriction_` — port
//! of USF_LateralFriction.F90 (2D lateral drag `-g + K*|u|^(m-1)*u`).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated: `g(2)`, `mm`, `FlowSolverName`, `auxReal`
//!   (sized `MaxElementNodes`, reused per call). The `material =>`
//!   pointer is dead in the source (never read); documented, not
//!   replicated.
//! - `g(1)/g(2)` and `mm` each overwrite the shared `GotIt`: only the
//!   SECOND call's `Found` gates the Fatal, verbatim.
//! - `NVelo` accumulates `Velo(i)**2.0` with a real exponent (pow per
//!   element, then `SQRT`), verbatim order.
//! - `LateralFriction_x/y` forward to the 4-arg form with axis 1/2 —
//!   same-module calls, like the source's internal references.
//! - The misaddressed `CALL Info('ShapeFactorGravity', ...)` on the
//!   missing-velocity path is replicated verbatim (wrong caller and all).

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

struct LateralState {
    first: bool,
    g: [c_double; 2],
    mm: c_double,
    flow_name: String,
    aux: Vec<c_double>,
}

static LATERAL_STATE: Mutex<LateralState> = Mutex::new(LateralState {
    first: true,
    g: [0.0; 2],
    mm: 0.0,
    flow_name: String::new(),
    aux: Vec::new(),
});

pub unsafe fn lateral_friction_x(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let x = unsafe { *v };
    lateral_friction(m, node, x, 1)
}

pub unsafe fn lateral_friction_y(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let x = unsafe { *v };
    lateral_friction(m, node, x, 2)
}

pub(crate) fn lateral_friction(m: *mut c_void, node: c_int, x: c_double, axis: c_int) -> c_double {
    let _ = x;
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let mut h_bf: *const c_void = std::ptr::null();
    unsafe { usf_get_bodyforce(&mut h_bf) };
    {
        let mut st = LATERAL_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            let maxnodes = unsafe { usf_max_element_nodes() } as usize;
            st.aux = vec![0.0; maxnodes];
            let mut got = false;
            st.g[0] = unsafe { const_real(h_bf, "Lateral Friction Gravity 1", &mut got) };
            st.g[1] = unsafe { const_real(h_bf, "Lateral Friction Gravity 2", &mut got) };
            if !got {
                unsafe { fatal("LateralFriction", "Gravity vector must be specified") };
            }
            st.mm = unsafe { const_real(h_bf, "Lateral Friction Exponent", &mut got) };
            if !got {
                unsafe { fatal("LateralFriction", "Lateral Friction Exponent must be defined") };
            }
            let mut g2 = false;
            let fsn = unsafe { get_string(h_bf, "Flow Solver Name", &mut g2) };
            st.flow_name = pad128(if g2 { &fsn } else { "Flow Solution" });
        }
    }
    let st = LATERAL_STATE.lock().unwrap();
    let (mm, g) = (st.mm, st.g);
    let flow_name = st.flow_name.trim_end().to_string();
    drop(st);
    let h_flow = unsafe { variable_get(h_vars, &flow_name, false, false) };
    if h_flow.is_null() {
        // Wrong-caller INFO is live in the source; the null handle then
        // faults on the very next dereference — pass-through, never guarded.
        unsafe { info("ShapeFactorGravity", "No variable for velocity associated.", 4) };
    }
    let fp = unsafe { var_perm_at(h_flow, node) };
    // 2D-only: stride 3, indices 1..2 — verbatim.
    let mut velo = [0.0 as c_double; 2];
    let mut nvelo = 0.0 as c_double;
    for (i, vv) in velo.iter_mut().enumerate() {
        // SAFETY: same 1-based slot arithmetic the source uses.
        *vv = unsafe { var_value_at_idx(h_flow, 3 * (fp - 1) + 1 + i as c_int) };
        nvelo += unsafe { usf_pow_dp(*vv, 2.0) };
    }
    nvelo = nvelo.sqrt();
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let n = unsafe { element_type_nodes(h_elem) } as usize;
    let mut st = LATERAL_STATE.lock().unwrap();
    let mut idx = vec![0 as c_int; n];
    unsafe { usf_element_nodeindexes(h_elem, idx.as_mut_ptr(), n as c_int) };
    let mut buf = vec![0.0 as c_double; n];
    let mut found = false;
    unsafe { get_real_into(h_bf, "Lateral Friction Coefficient", &mut found, &mut buf) };
    st.aux[..n].copy_from_slice(&buf);
    let mut pos = n; // EXIT-scan; `n+1` on no-match like the source.
    for (k, &id) in idx.iter().enumerate() {
        if node == id {
            pos = k;
            break;
        }
    }
    let mut kspring = st.aux[pos];
    drop(st);
    if nvelo > 1.0e-6 && (mm - 1.0).abs() > 1.0e-6 {
        kspring *= unsafe { usf_pow_dp(nvelo, mm - 1.0) };
    }
    let a = (axis - 1) as usize;
    g[a] - kspring * velo[a]
}
