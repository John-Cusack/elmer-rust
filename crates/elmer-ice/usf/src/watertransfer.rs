//! `epltoids_` — port of USF_WaterTransfer.F90 (EPL/IDS transfer).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated: all 14 allocatables (sized `3N`/`N` with
//!   `N = Model % MaxElementNodes` at first call), `FirstTime`, `Found`,
//!   `AllocationsDone` — including the first-call `DEALLOCATE` guarded by
//!   `AllocationsDone` (dead on first entry, kept verbatim) and the
//!   `WRITE(Message)+INFO(Level=4)` allocation note.
//! - `g` is zeroed over the WHOLE `(3,M)` store each associated call and
//!   `Gravity` keeps its previous values when no BodyForce is associated
//!   (SAVE, no recompute) — replicated via the state struct.
//!   divide-then-sum over the full store, including the zeroed tail —

use crate::support::*;

struct WaterState {
    first: bool,
    alloc_done: bool,
    m: usize,
    g: Vec<c_double>,
    density: Vec<c_double>,
    gravity: Vec<c_double>,
    epl_comp: Vec<c_double>,
    epl_porosity: Vec<c_double>,
    epl_thick: Vec<c_double>,
    epl_storing: Vec<c_double>,
    ids_comp: Vec<c_double>,
    ids_porosity: Vec<c_double>,
    ids_thick: Vec<c_double>,
    ids_storing: Vec<c_double>,
    ids_trans: Vec<c_double>,
    upper: Vec<c_double>,
    leak: Vec<c_double>,
    // `Transfer` is a non-SAVE RESULT local: gfortran static storage
    // retains the previous call's value when no assignment runs.
    transfer: c_double,
}

static STATE: Mutex<WaterState> = Mutex::new(WaterState {
    first: true,
    alloc_done: false,
    transfer: 0.0,
    m: 0,
    g: Vec::new(),
    density: Vec::new(),
    gravity: Vec::new(),
    epl_comp: Vec::new(),
    epl_porosity: Vec::new(),
    epl_thick: Vec::new(),
    epl_storing: Vec::new(),
    ids_comp: Vec::new(),
    ids_porosity: Vec::new(),
    ids_thick: Vec::new(),
    ids_storing: Vec::new(),
    ids_trans: Vec::new(),
    upper: Vec::new(),
    leak: Vec::new(),
});

pub unsafe fn epl_to_ids(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _x = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_mask = unsafe { variable_get(h_vars, "Open EPL", true, false) };
    if unsafe { var_value_at(h_mask, node) } >= 0.0 {
        STATE.lock().unwrap().transfer = 0.0;
        return 0.0;
    }
    let h_epl = unsafe { variable_get(h_vars, "EPLHead", true, false) };
    let h_ids = unsafe { variable_get(h_vars, "IDSHead", true, false) };
    let epl_at = unsafe { var_value_at(h_epl, node) };
    let ids_at = unsafe { var_value_at(h_ids, node) };
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_mat: *const c_void = std::ptr::null();
    unsafe { usf_get_material(h_elem, &mut h_mat) };
    let nn = unsafe { usf_get_element_nof_nodes(h_elem) } as usize;
    let mut h_c: *const c_void = std::ptr::null();
    unsafe { usf_get_constants(&mut h_c) };
    let mut h_bf: *const c_void = std::ptr::null();
    unsafe { usf_get_bodyforce(&mut h_bf) };
    // Element indexes (1-based Fortran) for the ListGetReal calls.
    let mut idx = vec![0 as c_int; nn];
    unsafe { usf_element_nodeindexes(h_elem, idx.as_mut_ptr(), nn as c_int) };
    // Material / body-force preconditions, verbatim messages.
    let mut st = STATE.lock().unwrap();
    if st.first {
        st.first = false;
        if st.alloc_done {
            st.g.clear();
            st.density.clear();
            st.gravity.clear();
            st.epl_comp.clear();
            st.epl_porosity.clear();
            st.epl_thick.clear();
            st.epl_storing.clear();
            st.ids_comp.clear();
            st.ids_porosity.clear();
            st.ids_thick.clear();
            st.ids_storing.clear();
            st.ids_trans.clear();
            st.upper.clear();
            st.leak.clear();
        }
        let msize = unsafe { usf_max_element_nodes() } as usize;
        st.m = msize;
        st.g = vec![0.0; 3 * msize];
        st.density = vec![0.0; msize];
        st.gravity = vec![0.0; msize];
        st.epl_comp = vec![0.0; msize];
        st.epl_porosity = vec![0.0; msize];
        st.epl_thick = vec![0.0; msize];
        st.epl_storing = vec![0.0; msize];
        st.ids_comp = vec![0.0; msize];
        st.ids_porosity = vec![0.0; msize];
        st.ids_thick = vec![0.0; msize];
        st.ids_storing = vec![0.0; msize];
        st.ids_trans = vec![0.0; msize];
        st.upper = vec![0.0; msize];
        st.leak = vec![0.0; msize];
        unsafe { info("WaterTransfer", "Memory allocation done", 4) };
        st.alloc_done = true;
    }
    if h_mat.is_null() {
        unsafe { fatal("WaterTransfer", "No Material found for boundary element no. ") };
    }
    let f = {
        let mut ff: c_int = 0;
        let id = unsafe { usf_get_material_id(h_elem, &mut ff) };
        let _ = id;
        ff != 0
    };
    if !f {
        unsafe { fatal("WaterTransfer", "No Material ID found for boundary element no. ") };
    }
    let mut found = false;
    let watcomp = {
        let w = unsafe { const_real(h_c, "Water Compressibility", &mut found) };
        if !found { 5.04e-4 } else { w }
    };
    let mut buf = vec![0.0 as c_double; nn];
    unsafe { list_get_real_into(h_mat, "Water Density", &idx, true, &mut found, &mut buf) };
    st.density[..nn].copy_from_slice(&buf);
    if !h_bf.is_null() {
        let bf_id = unsafe { usf_get_bodyforce_id() };
        let mut h_bfv: *const c_void = std::ptr::null();
        unsafe { usf_bodyforce_values(bf_id, &mut h_bfv) };
        let mut f1 = false;
        let b1 = unsafe { const_real(h_bfv, "Flow BodyForce 1", &mut f1) };
        let mut f2 = false;
        let b2 = unsafe { const_real(h_bfv, "Flow BodyForce 2", &mut f2) };
        let mut f3 = false;
        let b3 = unsafe { const_real(h_bfv, "Flow BodyForce 3", &mut f3) };
        for a in st.g.iter_mut() {
            *a = 0.0;
        }
        let m = st.m;
        for j in 0..nn {
            st.g[j] = b1;
            st.g[m + j] = b2;
            st.g[2 * m + j] = b3;
        }
        // Gravity(1:N) = SQRT(SUM(g**2.0/N)) — scalar broadcast.
        let s = unsafe { usf_sum_g2_over_n(st.g.as_ptr(), (3 * st.m) as c_int, nn as c_int) };
        let gv = s.sqrt();
        st.gravity[..nn].fill(gv);
    }
    unsafe { list_get_real_into(h_mat, "EPL Compressibility", &idx, true, &mut found, &mut buf) };
    st.epl_comp[..nn].copy_from_slice(&buf);
    unsafe { list_get_real_into(h_mat, "EPL Porosity", &idx, true, &mut found, &mut buf) };
    st.epl_porosity[..nn].copy_from_slice(&buf);
    unsafe { list_get_real_into(h_mat, "EPL Thickness", &idx, true, &mut found, &mut buf) };
    st.epl_thick[..nn].copy_from_slice(&buf);
    unsafe { list_get_real_into(h_mat, "IDS Compressibility", &idx, true, &mut found, &mut buf) };
    st.ids_comp[..nn].copy_from_slice(&buf);
    unsafe { list_get_real_into(h_mat, "IDS Porosity", &idx, true, &mut found, &mut buf) };
    st.ids_porosity[..nn].copy_from_slice(&buf);
    unsafe { list_get_real_into(h_mat, "IDS Thickness", &idx, true, &mut found, &mut buf) };
    st.ids_thick[..nn].copy_from_slice(&buf);
    for j in 0..nn {
        st.epl_storing[j] = st.epl_thick[j] * st.gravity[j] * st.epl_porosity[j] * st.density[j]
            * (watcomp + st.epl_comp[j] / st.epl_porosity[j]);
        st.ids_storing[j] = st.ids_thick[j] * st.gravity[j] * st.ids_porosity[j] * st.density[j]
            * (watcomp + st.ids_comp[j] / st.ids_porosity[j]);
    }
    unsafe { list_get_real_into(h_mat, "IDS Transmitivity", &idx, true, &mut found, &mut buf) };
    st.ids_trans[..nn].copy_from_slice(&buf);
    unsafe { list_get_real_into(h_mat, "IDSHead Upper Limit", &idx, false, &mut found, &mut buf) };
    st.upper[..nn].copy_from_slice(&buf);
    if !found {
        unsafe { info("WaterTransfer", "No upper limit of solution for element no.", 10) };
    }
    unsafe { list_get_real_into(h_mat, "Leakage Factor", &idx, true, &mut found, &mut buf) };
    st.leak[..nn].copy_from_slice(&buf);
    // N.B. source spells the bound `n` in the last two calls; identical value.
    for i in 0..nn {
        let k = idx[i];
        if node == k {
            if ids_at >= st.upper[i] && epl_at > ids_at {
                st.transfer = 0.0;
            } else {
                st.transfer =
                    st.ids_trans[i] * (epl_at - ids_at) / (st.ids_thick[i] * st.leak[i]);
                if st.transfer > 0.0 {
                    st.transfer *= st.epl_storing[i];
                } else if st.transfer < 0.0 {
                    st.transfer *= st.ids_storing[i];
                }
            }
            break;
        }
    }
    st.transfer
}
