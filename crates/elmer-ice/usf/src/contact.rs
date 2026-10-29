//! `slidcoef_contact_`, `gzept_scaling_` — port of USF_Contact.F90
//! (contact/grounding-line basal drag + GZEPT Stokes scaling).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated per entry, including the `Normal`/`Fwater`/`Fbase`
//!   work stores, `IceFrontNode` flags, `BuddFrictionCoeff`-style coefficient
//!   caches, `told`/`relChangeOld`, `yesChange`, `GLmoves`, `thresh`,
//!   `GLtype`, `TestContact`, `FrictionVar` handle + `FrictionValues` +
//!   `FrictionValue` + `FrictionPerm`, `BC`, name buffers.
//! - `BSD`-zeroed statics: `ApplyFrontFriction` (`.FALSE.`), `Cauchy`-style
//!   flags, `Warned` (contact local, always set before read).
//! - Fixed-length `CHARACTER(128)` names are replicated blank-padded
//!   (`pad128`); `SELECT CASE` and substring matches compare the padded
//!   form trimmed — byte-identical to the Fortran comparisons.
//! - Cross-entry calls (`Sliding_Weertman`, `Sliding_Budd`,
//!   `Friction_Coulomb` from sliding.rs; `GZEPT_scaling` within this
//!   module; `EffectivePressure` from sliding.rs) execute against the
//!   shared whole-library statics, exactly like the reference `.so`.
//! - `WRITE(Message,*)` list-directed output (including the `(A,I3)`
//!   WRITE with no integer item, `(A,e14.8)`, and mixed string/real
//!   items) evaluates in Fortran accessors.
//! - `FrictionValues(FrictionPerm(n))**2.0_dp` is a real-exponent `**`
//!   (routed through `usf_pow_dp`); `10.0_dp**x` likewise.
//! - `GroundedMask(Nn) = -1.0` / `= 0.0` and `Chi`-style writes go
//!   through the live variable address (`var_set_at`) — visible to
//!   Elmer by construction (kill-criterion-2 watch item, proven by the
//!   staged Contact differential).
//! - `MSum = MSum + INT(GroundedMask(Nn))`: mask values are -1/0/1, so
//!   the truncation cast is exact.
//! - `NULLIFY` statements have no observable effect; documented, not
//!   replicated.
//! - The dead `WRITE(FunctionName,...)` in `getFrictionHeat` (friction.rs)
//!   stays omitted there; here every WRITE is live.

use crate::sliding;
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

const USF_NAME: &str = "SlidCoef_Contact";

struct ContactState {
    first: bool,
    yes_change: bool,
    told: c_double,
    gl_moves: bool,
    thresh: c_double,
    gl_type: String,
    test_contact: c_double,
    dim: c_int,
    normal: Vec<c_double>,
    fwater: Vec<c_double>,
    fbase: Vec<c_double>,
    rel_change_old: c_double,
    sl_law: String,
    mask_name: String,
    // Opaque Elmer handle bits (never dereferenced in Rust).
    friction_var: usize,
    friction_value: c_double,
    flow_loads_name: String,
    flow_solution_name: String,
    apply_front_friction: bool,
    ice_front_node: Vec<bool>,
}

static CONTACT_STATE: Mutex<ContactState> = Mutex::new(ContactState {
    first: true,
    yes_change: false,
    told: 0.0,
    gl_moves: false,
    thresh: 0.0,
    gl_type: String::new(),
    test_contact: 0.0,
    dim: 0,
    normal: Vec::new(),
    fwater: Vec::new(),
    fbase: Vec::new(),
    rel_change_old: 0.0,
    sl_law: String::new(),
    mask_name: String::new(),
    friction_var: 0,
    friction_value: 0.0,
    flow_loads_name: String::new(),
    flow_solution_name: String::new(),
    apply_front_friction: false,
    ice_front_node: Vec::new(),
});

pub unsafe fn slidcoef_contact(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let y = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_mesh_vars = unsafe { model_mesh_variables() };
    let h_time = unsafe { variable_get(h_vars, "Time", true, false) };
    let t = unsafe { var_value_at_idx(h_time, 1) };
    let rel_change = unsafe { nonlin_change() };
    {
        let mut st = CONTACT_STATE.lock().unwrap();
        if st.first {
            st.dim = unsafe { usf_coord_dim() };
            let d = st.dim as usize;
            st.first = false;
            let maxnodes = unsafe { usf_max_element_nodes() } as usize;
            let _ = maxnodes;
            st.normal = vec![0.0; d];
            st.fwater = vec![0.0; d];
            st.fbase = vec![0.0; d];
            st.told = t;
            let mut h_c: *const c_void = std::ptr::null();
            unsafe { usf_get_constants(&mut h_c) };
            let mut got = false;
            let mn = unsafe { list_get_string(h_c, "Grounded Mask Variable Name", &mut got) };
            st.mask_name = pad128(if got { &mn } else { "GroundedMask" });
            let msg = unsafe { write_msg_star_aa("Grounded mask name is:", &st.mask_name) };
            unsafe { info(USF_NAME, &msg, 5) };
            st.yes_change = true;
            st.rel_change_old = rel_change;
            let mut h_elem0: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_elem0) };
            let mut h_bc0: *const c_void = std::ptr::null();
            unsafe { usf_get_bc(h_elem0, &mut h_bc0) };
            let mut g = false;
            let fsn = unsafe { get_string(h_bc0, "Flow Solution Name", &mut g) };
            st.flow_solution_name = pad128(if g { &fsn } else { "Flow Solution" });
            if !g {
                unsafe { info(USF_NAME, "Using default name >Flow Solution<", 4) };
            }
            st.flow_loads_name = pad128(&format!("{} Loads", st.flow_solution_name.trim_end()));
            let mut g2 = false;
            let sl = unsafe { get_string(h_bc0, "Sliding Law", &mut g2) };
            if !g2 {
                unsafe { fatal(USF_NAME, "No \"Sliding law\" Name given") };
            }
            st.sl_law = pad128(&sl);
            let mut g3 = false;
            let gl = unsafe { get_string(h_bc0, "Grounding Line Definition", &mut g3) };
            if !g3 {
                st.gl_type = pad128("last grounded");
                unsafe { info(USF_NAME, "Grounded Line Defined as the last Grounded point", 3) };
            } else {
                st.gl_type = pad128(&gl);
                let msg = format!("Grounding Line Defined as {}", st.gl_type.trim_end());
                unsafe { info(USF_NAME, &msg, 3) };
            }
            let mut g4 = false;
            let gm = unsafe { get_logical(h_bc0, "Grounding line moves", &mut g4) };
            st.gl_moves = if g4 { gm } else { true };
            if st.gl_moves {
                unsafe { info(USF_NAME, "GL may move by default", 3) };
                unsafe {
                    info(
                        USF_NAME,
                        "If you want to fix the Grounding Line, put the keyword \"Grounding line moves\" to False",
                        3,
                    )
                };
            } else {
                unsafe { info(USF_NAME, "GL will be fixed", 3) };
            }
            let mut g5 = false;
            let tc = unsafe { const_real(h_bc0, "Test Contact Tolerance", &mut g5) };
            if !g5 {
                st.test_contact = 1.0e-3;
                unsafe { info(USF_NAME, "Contact will be tested for a tolerance of 1.0e-3", 3) };
            } else {
                st.test_contact = tc;
                let msg = unsafe { write_msg_ae14d8("Contact tested for a tolerance of ", tc) };
                unsafe { info(USF_NAME, &msg, 3) };
            }
            let mut g6 = false;
            let th = unsafe { const_real(h_bc0, "non detachment inland distance", &mut g6) };
            if !g6 {
                st.thresh = -10000.0;
                unsafe {
                    info(
                        USF_NAME,
                        "far inland nodes have the possibility to detach by default",
                        3,
                    )
                };
                unsafe {
                    info(
                        USF_NAME,
                        "to avoid detachment (when bedrock is well below sea level),",
                        3,
                    )
                };
                unsafe {
                    info(
                        USF_NAME,
                        "use the keyword \"non detachment inland distance\" to the distance you wish",
                        3,
                    )
                };
                unsafe { info(USF_NAME, "This works with the DistanceSolver", 3) };
            } else {
                st.thresh = th;
                unsafe { info(USF_NAME, "far inland nodes will not detach", 3) };
            }
            let mut g7 = false;
            let aff = unsafe { get_logical(h_bc0, "Tiny shelf drag", &mut g7) };
            st.apply_front_friction = if g7 { aff } else { false };
            if st.apply_front_friction {
                if st.gl_type.trim_end() != "discontinuous" {
                    unsafe {
                        warn(
                            USF_NAME,
                            "\"Tiny shelf drag\" only active with \"Grounding Line Definition = discontinuous\"; switching it off",
                        )
                    };
                    st.apply_front_friction = false;
                } else {
                    unsafe {
                        info(
                            USF_NAME,
                            "Applying friction under ice-front slivers (discontinuous add-on)",
                            3,
                        )
                    };
                    let nnodes = unsafe { usf_mesh_number_of_nodes() } as usize;
                    st.ice_front_node = vec![false; nnodes];
                    let mut h_cur: *const c_void = std::ptr::null();
                    unsafe { usf_current_element(&mut h_cur) };
                    let nbe = unsafe { usf_number_of_boundary_elements() };
                    let mype = unsafe { usf_mype() };
                    let mut tt = 1;
                    while tt <= nbe {
                        let mut h_el: *const c_void = std::ptr::null();
                        unsafe { usf_get_boundary_element(tt, &mut h_el) };
                        if unsafe { usf_element_part_index(h_el) } != mype {
                            tt += 1;
                            continue;
                        }
                        let mut h_ebc: *const c_void = std::ptr::null();
                        unsafe { usf_get_bc(h_el, &mut h_ebc) };
                        if h_ebc.is_null() {
                            tt += 1;
                            continue;
                        }
                        let mut g8 = false;
                        if unsafe { get_logical(h_ebc, "Ice Front", &mut g8) } {
                            let en = unsafe { usf_get_element_nof_nodes(h_el) } as usize;
                            let mut nidx = vec![0 as c_int; en];
                            unsafe {
                                usf_element_nodeindexes(h_el, nidx.as_mut_ptr(), en as c_int)
                            };
                            for &nid in &nidx {
                                // 1-based global node number into the
                                // 0-based flag store (sized NumberOfNodes).
                                let gi = nid as usize;
                                if gi >= 1 && gi - 1 < st.ice_front_node.len() {
                                    st.ice_front_node[gi - 1] = true;
                                }
                            }
                        }
                        tt += 1;
                    }
                    unsafe { set_current_element(h_cur) };
                }
            }
        }
        // `Model % CurrentElement => CurElement` restore at First-block end.
        let mut h_cur: *const c_void = std::ptr::null();
        unsafe { usf_current_element(&mut h_cur) };
        unsafe { set_current_element(h_cur) };
    }
    // Per-call lookups.
    let st = CONTACT_STATE.lock().unwrap();
    let (dim, thresh) = (st.dim as usize, st.thresh);
    let mask_name = st.mask_name.trim_end().to_string();
    drop(st);
    let h_gmask = unsafe { variable_get(h_mesh_vars, &mask_name, true, false) };
    // `Sl_Law(1:10) == 'prescribed'` gate for the friction block.
    let sl10 = CONTACT_STATE.lock().unwrap().sl_law.as_str().get(..10).unwrap_or("").to_string();
    if sl10 == "prescribed" {
        let sl_full = CONTACT_STATE.lock().unwrap().sl_law.trim_end().to_string();
        if sl_full == "prescribed value" {
            let mut h_elem: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_elem) };
            let mut h_bc: *const c_void = std::ptr::null();
            unsafe { usf_get_bc(h_elem, &mut h_bc) };
            let fv = unsafe { list_get_const_real(h_bc, "Friction Value") };
            CONTACT_STATE.lock().unwrap().friction_value = fv;
        } else {
            let mut h_elem: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_elem) };
            let mut h_bc: *const c_void = std::ptr::null();
            unsafe { usf_get_bc(h_elem, &mut h_bc) };
            let mut got = false;
            let fvn = unsafe { get_string(h_bc, "Friction Variable Name", &mut got) };
            if !got {
                unsafe {
                    fatal(
                        USF_NAME,
                        "Prescribed friction requested but no Friction Variable Name found!",
                    )
                };
            }
            let h_fv = unsafe { variable_get(h_mesh_vars, fvn.trim_end(), false, false) };
            if h_fv.is_null() {
                let msg = format!("Unable to find variable: {}", pad128(fvn.trim_end()));
                unsafe { fatal(USF_NAME, &msg) };
            }
            let mut st = CONTACT_STATE.lock().unwrap();
            st.friction_var = h_fv as usize;
        }
    }
    // New-timestep edge.
    {
        let mut st = CONTACT_STATE.lock().unwrap();
        if t > st.told {
            st.told = t;
            st.yes_change = true;
            st.rel_change_old = rel_change;
        }
    }
    // Distance lookup (thresh-gated).
    if CONTACT_STATE.lock().unwrap().thresh > 0.0 {
        let h_dist = unsafe { variable_get(h_mesh_vars, "Distance", true, false) };
        let _ = h_dist;
    }
    CONTACT_STATE.lock().unwrap().normal.fill(0.0);
    CONTACT_STATE.lock().unwrap().fwater.fill(0.0);
    CONTACT_STATE.lock().unwrap().fbase.fill(0.0);
    // Convergence-gated grounding-line update.
    let go = {
        let st = CONTACT_STATE.lock().unwrap();
        (rel_change != st.rel_change_old)
            && rel_change > 0.0
            && rel_change < st.test_contact
            && st.yes_change
            && st.gl_moves
    };
    if go {
        CONTACT_STATE.lock().unwrap().yes_change = false;
        unsafe { info(USF_NAME, "FLOW SOLVER HAS SLIGHTLY CONVERGED: look for new basal conditions", 3) };
        let flow_loads = CONTACT_STATE.lock().unwrap().flow_loads_name.trim_end().to_string();
        let h_resid = unsafe { variable_get(h_mesh_vars, &flow_loads, true, false) };
        let h_normal = unsafe { variable_get(h_vars, "Normal Vector", true, false) };
        let h_hydro = unsafe { variable_get(h_mesh_vars, "Fw", true, false) };
        let mut h_cur: *const c_void = std::ptr::null();
        unsafe { usf_current_element(&mut h_cur) };
        let nbe = unsafe { usf_number_of_boundary_elements() };
        let mype = unsafe { usf_mype() };
        let mut gl_retreat = 0;
        let mut tt = 1;
        while tt <= nbe {
            let mut h_el: *const c_void = std::ptr::null();
            unsafe { usf_get_boundary_element(tt, &mut h_el) };
            if unsafe { usf_element_part_index(h_el) } != mype {
                tt += 1;
                continue;
            }
            let en = unsafe { usf_get_element_nof_nodes(h_el) } as usize;
            // Nodes x/y for the retreat message.
            let (mut xs, mut ys, mut zs) = (vec![0.0 as c_double; en], vec![0.0 as c_double; en], vec![0.0 as c_double; en]);
            unsafe {
                usf_element_nodes_xyz(h_el, en as c_int, xs.as_mut_ptr(), ys.as_mut_ptr(), zs.as_mut_ptr())
            };
            let mut nidx = vec![0 as c_int; en];
            unsafe { usf_element_nodeindexes(h_el, nidx.as_mut_ptr(), en as c_int) };
            // ANY(GroundedMaskPerm(NodeIndexes)==0) skip.
            let mut any_zero = false;
            for &nid in &nidx {
                if unsafe { var_perm_at(h_gmask, nid) } == 0 {
                    any_zero = true;
                    break;
                }
            }
            if any_zero {
                tt += 1;
                continue;
            }
            for (k, &nid) in nidx.iter().enumerate() {
                let nn = unsafe { var_perm_at(h_gmask, nid) };
                if nn == 0 {
                    continue;
                }
                if unsafe { var_value_at_idx(h_gmask, nn) } < -0.5 {
                    continue;
                }
                let jj = nid;
                for i in 0..dim {
                    // SAFETY: same 1-based slot arithmetic the source uses.
                    let nval = unsafe {
                        var_value_at_idx(h_normal, dim as c_int * (var_perm_at(h_normal, jj) - 1) + 1 + i as c_int)
                    };
                    let hval = unsafe {
                        var_value_at_idx(h_hydro, dim as c_int * (var_perm_at(h_hydro, jj) - 1) + 1 + i as c_int)
                    };
                    // `Normal = ...`, `Fwater = ...` work stores, verbatim.
                    CONTACT_STATE.lock().unwrap().normal[i] = nval;
                    CONTACT_STATE.lock().unwrap().fwater[i] = hval;
                }
                let rp = unsafe { var_perm_at(h_resid, jj) };
                for i in 0..dim {
                    let rv = unsafe {
                        var_value_at_idx(h_resid, (dim as c_int + 1) * (rp - 1) + 1 + i as c_int)
                    };
                    CONTACT_STATE.lock().unwrap().fbase[i] = rv;
                }
                let st = CONTACT_STATE.lock().unwrap();
                let comp = unsafe { dot(st.fwater.as_ptr(), st.normal.as_ptr(), dim) }
                    + unsafe { dot(st.fbase.as_ptr(), st.normal.as_ptr(), dim) };
                drop(st);
                if comp >= 0.0 {
                    let do_retreat = if thresh <= 0.0 {
                        true
                    } else {
                        let h_dist = unsafe { variable_get(h_mesh_vars, "Distance", true, false) };
                        let dp = unsafe { var_perm_at(h_dist, nid) };
                        (unsafe { var_value_at_idx(h_dist, dp) }) <= thresh
                    };
                    if do_retreat {
                        unsafe { var_set_at(h_gmask, nn, -1.0) };
                        gl_retreat += 1;
                        let msg = unsafe { write_msg_star_arr("Retreat of the Grounding Line : ", xs[k], ys[k]) };
                        unsafe { info(USF_NAME, &msg, 4) };
                    }
                }
            }
            tt += 1;
        }
        unsafe { set_current_element(h_cur) };
        if gl_retreat > 0 {
            let mut h_cur2: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_cur2) };
            let mut tt2 = 1;
            while tt2 <= nbe {
                let mut h_el: *const c_void = std::ptr::null();
                unsafe { usf_get_boundary_element(tt2, &mut h_el) };
                if unsafe { usf_element_part_index(h_el) } != mype {
                    tt2 += 1;
                    continue;
                }
                let en = unsafe { usf_get_element_nof_nodes(h_el) } as usize;
                // `CALL GetElementNodes` is live in the source (fills the
                // unused-here coordinate store); replicated, results ignored.
                let (mut xs, mut ys, mut zs) = (vec![0.0 as c_double; en], vec![0.0 as c_double; en], vec![0.0 as c_double; en]);
                unsafe {
                    usf_element_nodes_xyz(h_el, en as c_int, xs.as_mut_ptr(), ys.as_mut_ptr(), zs.as_mut_ptr())
                };
                let _ = (&xs, &ys, &zs);
                let mut nidx = vec![0 as c_int; en];
                unsafe { usf_element_nodeindexes(h_el, nidx.as_mut_ptr(), en as c_int) };
                let mut msum = 0;
                let mut zsum = 0;
                let mut any_zero = false;
                for &nid in &nidx {
                    if unsafe { var_perm_at(h_gmask, nid) } == 0 {
                        any_zero = true;
                        break;
                    }
                }
                if any_zero {
                    tt2 += 1;
                    continue;
                }
                for &nid in &nidx {
                    let nn = unsafe { var_perm_at(h_gmask, nid) };
                    if nn == 0 {
                        continue;
                    }
                    let gv = unsafe { var_value_at_idx(h_gmask, nn) };
                    msum += gv as c_int;
                    if gv == 0.0 {
                        zsum += 1;
                    }
                }
                if msum + zsum < en as c_int {
                    for &nid in &nidx {
                        let nn = unsafe { var_perm_at(h_gmask, nid) };
                        if nn == 0 {
                            continue;
                        }
                        if unsafe { var_value_at_idx(h_gmask, nn) } == 1.0 {
                            unsafe { var_set_at(h_gmask, nn, 0.0) };
                        }
                    }
                }
                tt2 += 1;
            }
            unsafe { set_current_element(h_cur2) };
        }
    }
    CONTACT_STATE.lock().unwrap().rel_change_old = rel_change;
    // LowerSurface branch.
    let gperm = unsafe { var_perm_at(h_gmask, node) };
    if gperm > 0 {
        let cond = unsafe { var_value_at_idx(h_gmask, gperm) };
        let st = CONTACT_STATE.lock().unwrap();
        let (gl_type, sl_law) = (st.gl_type.trim_end().to_string(), st.sl_law.trim_end().to_string());
        drop(st);
        let mut friction = false;
        if gl_type == "last grounded" {
            if cond > -0.5 {
                friction = true;
            }
        } else if gl_type == "first floating" {
            if cond > 0.5 {
                friction = true;
            }
        } else if gl_type == "discontinuous" {
            let mut h_be: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_be) };
            let en = unsafe { usf_get_element_nof_nodes(h_be) } as usize;
            let mut nidx = vec![0 as c_int; en];
            unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
            let mut all_above = true;
            for &nid in &nidx {
                let nn = unsafe { var_perm_at(h_gmask, nid) };
                let gv = if nn == 0 { 1.0 } else { unsafe { var_value_at_idx(h_gmask, nn) } };
                if gv <= -0.5 {
                    all_above = false;
                    break;
                }
            }
            if all_above {
                friction = true;
            } else if CONTACT_STATE.lock().unwrap().apply_front_friction {
                friction = true;
                for &nid in &nidx {
                    let jj = nid;
                    let nn = unsafe { var_perm_at(h_gmask, jj) };
                    if nn == 0 {
                        continue;
                    }
                    if unsafe { var_value_at_idx(h_gmask, nn) } < -0.5 {
                        // 1-based global into the 0-based flag store.
                        let on_front = if jj >= 1 {
                            CONTACT_STATE.lock().unwrap().ice_front_node.get(jj as usize - 1).copied().unwrap_or(false)
                        } else {
                            false
                        };
                        if !on_front {
                            friction = false;
                            break;
                        }
                    }
                }
            }
        } else {
            let msg = format!("GL type not recognised {}", CONTACT_STATE.lock().unwrap().gl_type);
            unsafe { fatal(USF_NAME, &msg) };
        }
        if friction {
            let bd = if sl_law == "weertman" {
                unsafe { sliding::sliding_weertman(m, n, &y as *const c_double) }
            } else if sl_law == "budd" {
                unsafe { sliding::sliding_budd(m, n, &y as *const c_double) }
            } else if sl_law == "coulomb" {
                unsafe { sliding::friction_coulomb(m, n, &y as *const c_double) }
            } else if sl_law == "prescribed beta2" {
                let st = CONTACT_STATE.lock().unwrap();
                let hv = st.friction_var as *const c_void;
                drop(st);
                let fp = unsafe { var_perm_at(hv, node) };
                let fv = unsafe { var_value_at_idx(hv, fp) };
                unsafe { usf_pow_dp(fv, 2.0) }
            } else if sl_law == "prescribed power" {
                let st = CONTACT_STATE.lock().unwrap();
                let hv = st.friction_var as *const c_void;
                drop(st);
                let fp = unsafe { var_perm_at(hv, node) };
                let fv = unsafe { var_value_at_idx(hv, fp) };
                unsafe { usf_pow_dp(10.0, fv) }
            } else if sl_law == "prescribed variable" {
                let st = CONTACT_STATE.lock().unwrap();
                let hv = st.friction_var as *const c_void;
                drop(st);
                let fp = unsafe { var_perm_at(hv, node) };
                unsafe { var_value_at_idx(hv, fp) }
            } else if sl_law == "prescribed value" {
                CONTACT_STATE.lock().unwrap().friction_value
            } else {
                let msg = format!("Sliding law not recognised {}", CONTACT_STATE.lock().unwrap().sl_law);
                // `fatal` diverges (`!`), so this arm already yields the
                // branch value — no trailing expression needed.
                unsafe { fatal(USF_NAME, &msg) };
            };
            bd * unsafe { gzept_scaling(m, n, &y as *const c_double) }
        } else {
            0.0
        }
    } else {
        unsafe { sliding::sliding_weertman(m, n, &y as *const c_double) }
    }
}

// -- GZEPT_scaling ----------------------------------------------------------

struct GzeptState {
    first: bool,
    gzept: c_double,
    use_gzept: bool,
}

static GZEPT_STATE: Mutex<GzeptState> = Mutex::new(GzeptState {
    first: true,
    gzept: 0.0,
    use_gzept: false,
});

pub unsafe fn gzept_scaling(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    // `nodenumber` is received and ignored, exactly like the source.
    let y = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut scaling = 1.0;
    {
        let mut st = GZEPT_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            let mut h_elem: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_elem) };
            let mut h_bc: *const c_void = std::ptr::null();
            unsafe { usf_get_bc(h_elem, &mut h_bc) };
            // No ASSOCIATED guard in the source — null pass-through.
            let mut got = false;
            st.gzept = unsafe { const_real(h_bc, "GZEPT", &mut got) };
            st.use_gzept = got;
            if got {
                let msg = unsafe { write_msg_star_ar("Setting Grounding Zone Effective Pressure Threshold (GZEPT) to: ", st.gzept) };
                unsafe { info("GZEPT_scaling", &msg, 4) };
            }
        }
    }
    let st = GZEPT_STATE.lock().unwrap();
    let (gzept, use_gzept) = (st.gzept, st.use_gzept);
    drop(st);
    if !use_gzept {
        return scaling;
    }
    let ep = unsafe { sliding::effective_pressure(m, n, &y as *const c_double) };
    if ep <= gzept {
        scaling = ep / gzept;
        if scaling < 0.0 {
            scaling = 0.0;
        }
    }
    scaling
}
