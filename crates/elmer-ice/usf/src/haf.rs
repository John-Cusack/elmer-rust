//! `calculate_haf_` — port of USF_Haf.F90 (volume above flotation).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated: `FirstTime`, `UseGroundedMask`, `thickName`,
//!   `GroundedMaskName`; `rhoRatio`/`Thickness` are source-LOCALS read
//!   without a SAVE, which gfortran keeps in static storage
//!   (`-fno-automatic` default) — replicated as statics so the
//!   first-call value persists exactly like the reference.
//! - `GroundedMaskName` starts as zero-filled static storage in the
//!   reference (BSS); the port starts it empty — both arms only use it
//!   when `GotIt` was true, and both `Fatal` identically otherwise.
//! - Name buffers truncate to `MAX_NAME_LEN` (128) exactly like the
//!   fixed-length source declarations.

use crate::support::*;

const MAX_NAME_LEN: usize = 128;

struct HafState {
    first: bool,
    use_mask: bool,
    thick_name: String,
    mask_name: String,
    rho_ratio: c_double,
}

static STATE: Mutex<HafState> = Mutex::new(HafState {
    first: true,
    use_mask: false,
    // Zero-filled static storage, exactly like the reference's BSS SAVE
    // (Fortran CHARACTER(128) zero bytes); only used when the key lookup
    // missed, in which case both arms Fatal identically.
    thick_name: String::new(),
    mask_name: String::new(),
    rho_ratio: 0.0,
});

/// Truncate to `MAX_NAME_LEN` bytes, like the fixed-length source
/// declaration (keys are ASCII in every test).
fn truncate128(s: &str) -> String {
    let b = s.as_bytes();
    let m = b.len().min(MAX_NAME_LEN);
    String::from_utf8_lossy(&b[..m]).into_owned()
}

pub unsafe fn calculate_haf(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node_height = unsafe { *v };
    let node = unsafe { *n };
    unsafe { usf_begin(m) };
    {
        let mut st = STATE.lock().unwrap();
        if st.first {
            st.first = false;
            let mut h_elem: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_elem) };
            let mut h_bc: *const c_void = std::ptr::null();
            unsafe { usf_get_bc(h_elem, &mut h_bc) };
            if h_bc.is_null() {
                unsafe { fatal("USF_Haf", "No BC Found") };
            }
            let mut got = false;
            let gm = unsafe { get_string(h_bc, "Haf GroundedMask", &mut got) };
            if got {
                st.use_mask = true;
                st.mask_name = truncate128(&gm);
            } else {
                st.use_mask = false;
            }
            let th = unsafe { get_string(h_bc, "Haf Thickness", &mut got) };
            if !got {
                unsafe { fatal("USF_Haf", "Can not find >Haf Thickness< in BC") };
            }
            st.thick_name = truncate128(&th);
            let mut got2 = false;
            st.rho_ratio = unsafe { const_real(h_bc, "Haf rhoRatio", &mut got2) };
            if !got2 {
                unsafe { fatal("USF_Haf", "Can not find >Haf rhoRatio< in BC") };
            }
        }
    }
    let st = STATE.lock().unwrap();
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    // UnFoundFatal=.TRUE. in the source: missing variables Fatal inside.
    let h_mask = unsafe { variable_get(h_vars, &st.mask_name, true, false) };
    let h_thick = unsafe { variable_get(h_vars, &st.thick_name, true, false) };
    let thickness = unsafe { var_value_at(h_thick, node) };
    let mask_v = if st.use_mask {
        unsafe { var_value_at(h_mask, node) }
    } else {
        0.0
    };
    let rho_ratio = st.rho_ratio;
    drop(st);
    let mut haf = if node_height > 0.0 {
        thickness
    } else {
        thickness + node_height / rho_ratio
    };
    // `IF (UseGroundedMask)` re-read under lock; value fixed after init.
    if STATE.lock().unwrap().use_mask && mask_v < -0.5 {
        haf = 0.0;
    }
    haf
}
