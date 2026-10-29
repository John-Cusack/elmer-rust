//! `sliding_weertman_`, `friction_coulomb_`, `sliding_budd_`,
//! `effectivepressure_`, `freeslipshelves_` — port of USF_Sliding.F90
//! (Weertman/Budd/Coulomb basal drag, effective pressure, shelf mask).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated per entry (`normal`/`velo`/`Sig`/`Sn`/`DIM`/`Ind`
//!   stores, `FirstTime`, `FlowSolverName`, `t0`, `Cauchy`,
//!   `BuddFrictionCoeff`, `UseFloatation`, `WeertCoefName/Exp/Form`,
//!   `ConvertWeertman`, `GZEPT`-style `UseGZEPT` lives in contact.rs).
//!   Gfortran keeps non-SAVE locals static too; every local the source
//!   reads before writing on some path is replicated in the state
//!   struct: `elementNodeNumber` (Budd, stale on no-match), `Snn`
//!   (Coulomb, stale when `NVariable` is associated), `UseMask`-gated
//!   `frictionLoad` (contact.rs).
//! - `UnFoundFatal` in `Friction_Coulomb` is never assigned and not
//!   SAVE'd: BSS zero (`.FALSE.`) in the reference, literal `false` here.
//! - Cross-entry calls (`SlidCoef_Contact` and `GZEPT_scaling` call into
//!   this module) execute against the same statics — whole-library
//!   state sharing, exactly like the reference `.so`.
//! - Integer-exponent `**` (`Cw**2`) is repeated multiplication, while
//!   real-exponent `**` (`ut**(m-1.0)`, `Zab**q`, `10**Cw`) goes through
//!   `usf_pow_dp` — gfortran lowers those to `pow`, not multiplies.
//! - `ut**(m-1.0)` evaluates `m-1.0` in dp first, verbatim.
//! - Full-`SELECT CASE` string matches use the blank-padded 128-char
//!   form (`pad128`), then compare trimmed — byte-identical to the
//!   Fortran padded comparison.
//! - `10**Cw`: default-integer base with a real exponent promotes to dp
//!   `pow(10, Cw)` — routed through `usf_pow_dp` like every real `**`.
//! - `WRITE(Message,*)` list-directed output (including the `(A,I3)`
//!   WRITE with no integer item) evaluates in Fortran accessors.
//! - `Sig`/`Sn` row copies into contiguous temps precede `usf_dot`
//!   (copies are exact; strided slices never cross).
//! - The second `DGEEV` call in the reference runs on the matrix the
//!   first call destroyed (LAPACK overwrites `A` in all cases) —
//!   replicated call-for-call, bug-compatible by construction.
//! - `NULLIFY` statements have no observable effect; documented, not
//!   replicated.

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

// -- Sliding_Weertman -----------------------------------------------------

struct WeertState {
    first: bool,
    dim: c_int,
    ssa: bool,
    flow_name: String,
    normal: Vec<c_double>,
    velo: Vec<c_double>,
}

static WEERT_STATE: Mutex<WeertState> = Mutex::new(WeertState {
    first: true,
    dim: 0,
    ssa: false,
    flow_name: String::new(),
    normal: Vec::new(),
    velo: Vec::new(),
});

pub unsafe fn sliding_weertman(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _x = unsafe { *v };
    unsafe { usf_begin(m) };
    {
        let mut st = WEERT_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            st.dim = unsafe { usf_coord_dim() };
            let d = st.dim;
            if d == 2 || d == 3 {
                st.normal = vec![0.0; d as usize];
                st.velo = vec![0.0; d as usize];
            } else {
                unsafe { fatal("USF_sliding", "Bad dimension of the problem") };
            }
            let h_sv = unsafe { current_solver_values() };
            let mut got = false;
            let mut name = unsafe { get_string(h_sv, "Flow Solver Name", &mut got) };
            if !got {
                name = "Flow Solution".to_string();
            }
            st.flow_name = pad128(&name);
            // SELECT CASE (FlowSolverName): full padded compare.
            if st.flow_name.trim_end() == "ssabasalflow" {
                st.ssa = true;
            }
            let msg = unsafe { write_msg_star_asal("Flow Solver Name:", name.trim_end(), " SSA:", st.ssa) };
            unsafe { info("Sliding_Weertman", &msg, 3) };
        }
    }
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_bc: *const c_void = std::ptr::null();
    unsafe { usf_get_bc(h_elem, &mut h_bc) };
    if h_bc.is_null() {
        unsafe { fatal("Sliding_Weertman", "No BC Found") };
    }
    let nn = unsafe { usf_get_element_nof_nodes(h_elem) } as usize;
    let mut idx = vec![0 as c_int; nn];
    unsafe { usf_element_nodeindexes(h_elem, idx.as_mut_ptr(), nn as c_int) };
    let mut aux = vec![0.0 as c_double; nn];
    let mut found = false;
    unsafe { get_real_into(h_bc, "Weertman Friction Coefficient", &mut found, &mut aux) };
    if !found {
        unsafe { fatal("USF_sliding", "Need a Friction Coefficient for the Weertman sliding law") };
    }
    let mut ci = nn; // stale-safe: reference reads aux(n+1) on no-match (UB);
    for (k, &id) in idx.iter().enumerate() {
        if node == id {
            ci = k;
            break;
        }
    }
    // N.B. `ci == nn` (no match) would read out of bounds in the
    // reference too; the invariant (node in element) always holds.
    let c = aux[ci];
    let m_exp = unsafe {
        let mm = const_real(h_bc, "Weertman Exponent", &mut found);
        if !found {
            fatal("USF_sliding", "Need an Exponent for the Weertman sliding law");
        }
        mm
    };
    let ut0 = unsafe {
        let u = const_real(h_bc, "Weertman Linear Velocity", &mut found);
        if !found {
            fatal("USF_sliding", "Need a Linear Velocity for the Weertman sliding law");
        }
        u
    };
    let st = WEERT_STATE.lock().unwrap();
    let (dim, ssa) = (st.dim as usize, st.ssa);
    let flow_name = st.flow_name.trim_end().to_string();
    drop(st);
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_flow = unsafe { variable_get(h_vars, &flow_name, true, false) };
    let fp = unsafe { var_perm_at(h_flow, node) };
    let (un, mut ut);
    if !ssa {
        let h_normal = unsafe { variable_get(h_vars, "Normal Vector", true, false) };
        let np = unsafe { var_perm_at(h_normal, node) };
        let mut st = WEERT_STATE.lock().unwrap();
        for i in 0..dim {
            // SAFETY: same 1-based slot arithmetic the source uses.
            st.normal[i] = -unsafe { var_value_at_idx(h_normal, dim as c_int * (np - 1) + 1 + i as c_int) };
            st.velo[i] = unsafe {
                var_value_at_idx(h_flow, (dim as c_int + 1) * (fp - 1) + 1 + i as c_int)
            };
        }
        un = unsafe { dot(st.normal.as_ptr(), st.velo.as_ptr(), dim) };
        // (velo-un*normal)**2.0: diff in Rust, pow-sum in Fortran.
        let mut diff = vec![0.0 as c_double; dim];
        for i in 0..dim {
            diff[i] = st.velo[i] - un * st.normal[i];
        }
        drop(st);
        ut = unsafe { sum_pow2(diff.as_ptr(), dim) }.sqrt();
    } else {
        let mut st = WEERT_STATE.lock().unwrap();
        for (i, vv) in st.velo.iter_mut().enumerate().take(dim - 1) {
            *vv = unsafe {
                var_value_at_idx(h_flow, (dim as c_int - 1) * (fp - 1) + 1 + i as c_int)
            };
        }
        let vv = st.velo.clone();
        drop(st);
        ut = unsafe { sum_pow2(vv.as_ptr(), dim - 1) }.sqrt();
        un = 0.0;
    }
    let _ = un;
    ut = ut.max(ut0);
    (c * unsafe { usf_pow_dp(ut, m_exp - 1.0) }).min(1.0e20)
}

// -- Friction_Coulomb -------------------------------------------------------

struct CoulombState {
    first: bool,
    dim: c_int,
    ind: [[c_int; 3]; 3],
    t0: c_double,
    sig: Vec<c_double>,
    normal: Vec<c_double>,
    velo: Vec<c_double>,
    sn: Vec<c_double>,
    // `Snn` is a non-SAVE local: gfortran static storage, stale when
    // `NVariable` is associated (the final `t==t0` guard reads it).
    snn: c_double,
}

static COULOMB_STATE: Mutex<CoulombState> = Mutex::new(CoulombState {
    first: true,
    dim: 0,
    ind: [[0; 3]; 3],
    t0: 0.0,
    sig: Vec::new(),
    normal: Vec::new(),
    velo: Vec::new(),
    sn: Vec::new(),
    snn: 0.0,
});

pub unsafe fn friction_coulomb(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _y = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_time = unsafe { variable_get(h_vars, "Time", true, false) };
    let t = unsafe { var_value_at_idx(h_time, 1) };
    {
        let mut st = COULOMB_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            st.t0 = t;
            st.dim = unsafe { usf_coord_dim() };
            let d = st.dim;
            if d == 2 || d == 3 {
                st.sig = vec![0.0; (d * d) as usize];
                st.normal = vec![0.0; d as usize];
                st.velo = vec![0.0; d as usize];
                st.sn = vec![0.0; d as usize];
            } else {
                unsafe { fatal("Friction_Coulomb", "Bad dimension of the problem") };
            }
            st.ind[0][0] = 1;
            st.ind[1][1] = 2;
            st.ind[2][2] = 3;
            st.ind[0][1] = 4;
            st.ind[1][0] = 4;
            st.ind[1][2] = 5;
            st.ind[2][1] = 5;
            st.ind[0][2] = 6;
            st.ind[2][0] = 6;
        }
    }
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let nn0 = unsafe { usf_get_element_nof_nodes(h_elem) } as usize;
    let mut idx = vec![0 as c_int; nn0];
    unsafe { usf_element_nodeindexes(h_elem, idx.as_mut_ptr(), nn0 as c_int) };
    let mut h_bc: *const c_void = std::ptr::null();
    unsafe { usf_get_bc(h_elem, &mut h_bc) };
    if h_bc.is_null() {
        unsafe { fatal("Friction_Coulomb", "No BC Found") };
    }
    let n = nn0;
    let mut aux = vec![0.0 as c_double; n];
    let mut found = false;
    // EXIT-scan: node position in the element (`n` on no-match; the
    // reference reads one past the end there — invariant never hit).
    let scan = |idx: &[c_int], node: c_int| -> usize {
        for (k, &id) in idx.iter().enumerate() {
            if node == id {
                return k;
            }
        }
        idx.len()
    };
    unsafe { get_real_into(h_bc, "Friction Law Sliding Coefficient", &mut found, &mut aux) };
    if !found {
        unsafe {
            fatal(
                "Friction_Coulomb",
                "Need a Friction Law Sliding Coefficient for the Coulomb Friction Law",
            )
        };
    }
    let as_ = aux[scan(&idx, node)];
    unsafe { get_real_into(h_bc, "Friction Law Post-Peak Exponent", &mut found, &mut aux) };
    if !found {
        unsafe {
            fatal(
                "Friction_Coulomb",
                "Need a Friction Law Post-Peak Exponent    (>= 1) for the Coulomb Friction Law",
            )
        };
    }
    let q = aux[scan(&idx, node)];
    let a = unsafe { usf_pow_dp(q - 1.0, q - 1.0) / usf_pow_dp(q, q) };
    unsafe { get_real_into(h_bc, "Friction Law Maximum Value", &mut found, &mut aux) };
    if !found {
        unsafe {
            fatal(
                "Friction_Coulomb",
                "Need a Friction Law Maximum Value     (~ Max Bed Slope) for the Coulomb Friction Law",
            )
        };
    }
    let c = aux[scan(&idx, node)];
    let ut0 = unsafe {
        let u = const_real(h_bc, "Friction Law Linear Velocity", &mut found);
        if !found {
            fatal(
                "Friction_Coulomb",
                "Need a Friction Law Linear Velocity for the Coulomb Friction Law ",
            );
        }
        u
    };
    let m_exp = unsafe {
        let mm = const_real(h_bc, "Friction Law PowerLaw Exponent", &mut found);
        if !found {
            fatal(
                "Friction_Coulomb",
                "Need a Friction Law PowerLaw Exponent       (= n Glen law) for the Coulomb Friction Law",
            );
        }
        mm
    };
    let h_nvar = unsafe { variable_get(h_vars, "Effective Pressure", false, false) };
    let n_associated = !h_nvar.is_null();
    // Initial values required for definite assignment (every path
    // overwrites before any read); the inits themselves never survive.
    #[allow(unused_assignments)]
    let (mut ne, mut pw, mut cauchy) = (0.0, 0.0, false);
    if n_associated {
        let np = unsafe { var_perm_at(h_nvar, node) };
        ne = unsafe { var_value_at_idx(h_nvar, np) };
        pw = 0.0;
        cauchy = false;
    } else {
        unsafe { get_real_into(h_bc, "External Pressure", &mut found, &mut aux) };
        if !found {
            unsafe {
                fatal(
                    "Friction_Coulomb",
                    "Need External Pressure       Or Variable Effective Pressure",
                )
            };
        }
        pw = -aux[scan(&idx, node)];
        let h_parent = unsafe { boundary_parent(h_elem) };
        let mut h_mat: *const c_void = std::ptr::null();
        unsafe { usf_get_material(h_parent, &mut h_mat) };
        let mut got = false;
        cauchy = unsafe { list_get_logical(h_mat, "Cauchy", &mut got) };
        if !got {
            cauchy = false;
        }
    }
    drop(aux);
    let h_flow = unsafe { variable_get(h_vars, "Flow Solution", true, false) };
    let h_normal = unsafe { variable_get(h_vars, "Normal Vector", true, false) };
    let fp = unsafe { var_perm_at(h_flow, node) };
    let np = unsafe { var_perm_at(h_normal, node) };
    let mut st = COULOMB_STATE.lock().unwrap();
    let dim = st.dim as usize;
    let ind = st.ind;
    for i in 0..dim {
        // SAFETY: same 1-based slot arithmetic the source uses.
        st.normal[i] =
            -unsafe { var_value_at_idx(h_normal, dim as c_int * (np - 1) + 1 + i as c_int) };
        st.velo[i] = unsafe {
            var_value_at_idx(h_flow, (dim as c_int + 1) * (fp - 1) + 1 + i as c_int)
        };
    }
    let un = unsafe { dot(st.velo.as_ptr(), st.normal.as_ptr(), dim) };
    let mut diff = vec![0.0 as c_double; dim];
    for i in 0..dim {
        diff[i] = st.velo[i] - un * st.normal[i];
    }
    let mut ut = unsafe { sum_pow2(diff.as_ptr(), dim) }.sqrt();
    if !n_associated {
        let h_stress = unsafe { variable_get(h_vars, "Stress", true, false) };
        let spp = unsafe { var_perm_at(h_stress, node) };
        for i in 0..dim {
            for j in 0..dim {
                // SAFETY: column-major Sig(i,j) like the source: flat i + j*dim.
                let slot = 2 * dim as c_int * (spp - 1) + ind[i][j];
                st.sig[i + j * dim] = unsafe { var_value_at_idx(h_stress, slot) };
            }
        }
        let mut row = [0.0 as c_double; 3];
        for i in 0..dim {
            for j in 0..dim {
                row[j] = st.sig[i + j * dim];
            }
            st.sn[i] = unsafe { dot(row.as_ptr(), st.normal.as_ptr(), dim) };
        }
        st.snn = unsafe { dot(st.sn.as_ptr(), st.normal.as_ptr(), dim) };
        let pice = unsafe { var_value_at_idx(h_flow, (dim as c_int + 1) * fp) };
        if cauchy {
            if st.snn.abs() < 1.0e-10 * pice.abs() {
                st.snn = -pice;
            }
        } else {
            st.snn -= pice;
        }
        ne = -st.snn - pw;
    }
    let snn_now = st.snn;
    let t0 = st.t0;
    drop(st);
    let mut xi;
    let mut bdrag;
    if ne > 0.0 {
        ut = ut.max(ut0);
        xi = ut / (as_ * unsafe { usf_pow_dp(c * ne, m_exp) });
        xi = xi.min(1.0e20);
    } else {
        xi = 1.0e20;
        let msg = unsafe { write_msg_star_air("!!! Ne <=0, nodenumber", node, ne) };
        unsafe { info("Friction_Coulomb", &msg, 9) };
        ne = 0.0;
    }
    // Bdrag = C*Ne * ((Xi * ut**(-m)) / (1 + a*Xi**q))**(1/m).
    let t1 = unsafe { usf_pow_dp(ut, -m_exp) };
    let t2 = xi * t1;
    let t3 = unsafe { usf_pow_dp(xi, q) };
    let t5 = t2 / (1.0 + a * t3);
    bdrag = c * ne * unsafe { usf_pow_dp(t5, 1.0 / m_exp) };
    bdrag = bdrag.min(1.0e20);
    if (t - t0) == 0.0 && !n_associated && snn_now >= 0.0 {
        bdrag = 1.0e20;
    }
    bdrag
}

// -- Sliding_Budd -----------------------------------------------------------

struct BuddState {
    first: bool,
    dim: c_int,
    ssa: bool,
    use_floatation: bool,
    convert_weertman: bool,
    flow_name: String,
    weert_coef_name: String,
    weert_exp: c_double,
    weert_form: String,
    budd_friction_coeff: Vec<c_double>,
    normal: Vec<c_double>,
    velo: Vec<c_double>,
    // `H` (ice depth) is a non-SAVE local set only in the floatation
    // branch: gfortran static storage, stale (BSS 0.0 first) when the
    // `H_scaling` path runs without floatation.
    h: c_double,
}

static BUDD_STATE: Mutex<BuddState> = Mutex::new(BuddState {
    first: true,
    dim: 0,
    ssa: false,
    use_floatation: false,
    convert_weertman: false,
    flow_name: String::new(),
    weert_coef_name: String::new(),
    weert_exp: 0.0,
    weert_form: String::new(),
    budd_friction_coeff: Vec::new(),
    normal: Vec::new(),
    velo: Vec::new(),
    h: 0.0,
});

pub unsafe fn sliding_budd(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let z = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_bc: *const c_void = std::ptr::null();
    unsafe { usf_get_bc(h_elem, &mut h_bc) };
    if h_bc.is_null() {
        unsafe { fatal("Sliding_Budd", "No BC Found") };
    }
    {
        let mut st = BUDD_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            st.dim = unsafe { usf_coord_dim() };
            let d = st.dim;
            if d == 2 || d == 3 {
                st.normal = vec![0.0; d as usize];
                st.velo = vec![0.0; d as usize];
            } else {
                unsafe { fatal("Sliding_Budd", "Bad dimension of the problem") };
            }
            let maxnodes = unsafe { usf_max_element_nodes() } as usize;
            st.budd_friction_coeff = vec![0.0; maxnodes];
            let h_sv = unsafe { current_solver_values() };
            let mut got = false;
            let mut name = unsafe { get_string(h_sv, "Flow Solver Name", &mut got) };
            if !got {
                name = "Flow Solution".to_string();
            }
            st.flow_name = pad128(&name);
            if st.flow_name.trim_end() == "ssabasalflow" {
                st.ssa = true;
            }
            let mut g = false;
            let wcn = unsafe { get_string(h_bc, "Budd Conv Weertman coef name", &mut g) };
            if g {
                st.convert_weertman = true;
                st.weert_coef_name = pad128(&wcn);
                st.weert_exp = unsafe { const_real(h_bc, "Budd Conv Weertman exponent", &mut g) };
                if !g {
                    unsafe {
                        fatal(
                            "Sliding_Budd",
                            "Converting coef, need >Budd Conv Weertman exponent<",
                        )
                    };
                }
                let mut g2 = false;
                let wf = unsafe { get_string(h_bc, "Budd Conv Weertman formulation", &mut g2) };
                if !g2 {
                    unsafe {
                        fatal(
                            "Sliding_Budd",
                            "Converting coef, need >Budd Conv Weertman formulation<",
                        )
                    };
                }
                st.weert_form = pad128(&wf);
            } else {
                st.convert_weertman = false;
            }
        }
    }
    let element_nb_nodes = unsafe { usf_get_element_nof_nodes(h_elem) } as usize;
    let st = BUDD_STATE.lock().unwrap();
    let dim = st.dim as usize;
    // EXIT-scan over the boundary element's NodeIndexes.
    let mut idx = vec![0 as c_int; element_nb_nodes];
    unsafe { usf_element_nodeindexes(h_elem, idx.as_mut_ptr(), element_nb_nodes as c_int) };
    // EXIT-scan over the boundary element's NodeIndexes (1-based position;
    // normal DO completion leaves N+1, matching the source).
    let mut enn = element_nb_nodes + 1;
    for (k, &id) in idx.iter().enumerate() {
        if id == node {
            enn = k + 1;
            break;
        }
    }
    if h_elem.is_null() {
        unsafe { fatal("Sliding_Budd", "No boundary element found") };
    }
    let h_parent = unsafe { boundary_parent(h_elem) };
    let body_id = unsafe { usf_element_body_id(h_parent) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    drop(st);
    // Bodies/materials navigated through accessors (bodies are 1-based).
    let h_bodies = unsafe {
        let mut h_bv: *const c_void = std::ptr::null();
        usf_model_body_values(body_id, &mut h_bv);
        h_bv
    };
    let mut found = false;
    let material_id = unsafe { list_get_integer(h_bodies, "Material", true, &mut found) };
    let h_parent_mat = unsafe {
        let mut h_pm: *const c_void = std::ptr::null();
        usf_model_material_values(material_id, &mut h_pm);
        h_pm
    };
    if h_parent_mat.is_null() {
        // `WRITE(Message,'(A,I10,A,I10)')` then Fatal.
        let msg = format!(
            "No material values found for body no {:>10} under material id {:>10}",
            body_id, material_id
        );
        unsafe { fatal("Sliding_Budd", &msg) };
    }
    let rhoi = unsafe {
        let r = const_real(h_parent_mat, "Density", &mut found);
        if !found {
            fatal("Sliding_Budd", "Material property Density not found.");
        }
        r
    };
    // BuddFrictionCoeff(1:elementNbNodes) = ListGetReal(BC, ...).
    let mut buf = vec![0.0 as c_double; element_nb_nodes];
    unsafe { list_get_real_into(h_bc, "Budd Friction Coefficient", &idx, false, &mut found, &mut buf) };
    if !found {
        unsafe { fatal("Sliding_Budd", "Need a Friction Coefficient for the Budd sliding law") };
    }
    let mut st = BUDD_STATE.lock().unwrap();
    st.budd_friction_coeff[..element_nb_nodes].copy_from_slice(&buf);
    // `C = BuddFrictionCoeff(elementNodeNumber)` (1-based; `enn == N+1`
    // on no-match panics here where the reference reads out of bounds —
    // the invariant (node in element) always holds.
    let c = st.budd_friction_coeff[enn - 1];
    drop(st);
    let m_exp = unsafe {
        let mm = const_real(h_bc, "Budd Velocity Exponent", &mut found);
        if !found {
            fatal("Sliding_Budd", "Need a velocity Exponent for the Budd sliding law");
        }
        mm
    };
    let q = unsafe {
        let qq = const_real(h_bc, "Budd Zab Exponent", &mut found);
        if !found {
            fatal("Sliding_Budd", "Need a Zab Exponent for the Budd sliding law");
        }
        qq
    };
    let mut zab_offset = unsafe { const_real(h_bc, "Budd Zab Offset", &mut found) };
    if !found {
        zab_offset = 0.0;
    }
    let ut0 = unsafe {
        let u = const_real(h_bc, "Budd Linear Velocity", &mut found);
        if !found {
            fatal("Sliding_Budd", "Need a Linear Velocity for the Budd sliding law");
        }
        u
    };
    let g = unsafe {
        let gg = const_real(h_bc, "Budd Gravity", &mut found);
        if !found {
            fatal("Sliding_Budd", "Need Gravity for the Budd sliding law");
        }
        gg
    };
    let use_floatation = unsafe {
        let u = get_logical(h_bc, "Budd Floatation", &mut found);
        if !found {
            fatal("Sliding_Budd", "Need Floatation for the Budd sliding law");
        }
        u
    };
    BUDD_STATE.lock().unwrap().use_floatation = use_floatation;
    let h_scaling = unsafe {
        let h = get_logical(h_bc, "Budd Thickness Scaling", &mut found);
        if !found { false } else { h }
    };
    let st = BUDD_STATE.lock().unwrap();
    let (ssa2, flow_name) = (st.ssa, st.flow_name.trim_end().to_string());
    let (convert, weert_name, weert_exp, weert_form) = (
        st.convert_weertman,
        st.weert_coef_name.trim_end().to_string(),
        st.weert_exp,
        st.weert_form.trim_end().to_string(),
    );
    drop(st);
    let h_flow = unsafe { variable_get(h_vars, &flow_name, true, false) };
    let fp = unsafe { var_perm_at(h_flow, node) };
    let mut st = BUDD_STATE.lock().unwrap();
    let (un, mut ut);
    if !ssa2 {
        let h_normal = unsafe { variable_get(h_vars, "Normal Vector", true, false) };
        let np = unsafe { var_perm_at(h_normal, node) };
        for i in 0..dim {
            // SAFETY: same 1-based slot arithmetic the source uses.
            st.normal[i] =
                -unsafe { var_value_at_idx(h_normal, dim as c_int * (np - 1) + 1 + i as c_int) };
            st.velo[i] = unsafe {
                var_value_at_idx(h_flow, (dim as c_int + 1) * (fp - 1) + 1 + i as c_int)
            };
        }
        un = unsafe { dot(st.velo.as_ptr(), st.normal.as_ptr(), dim) };
        let mut diff = vec![0.0 as c_double; dim];
        for i in 0..dim {
            diff[i] = st.velo[i] - un * st.normal[i];
        }
        drop(st);
        ut = unsafe { sum_pow2(diff.as_ptr(), dim) }.sqrt();
    } else {
        let mut st = BUDD_STATE.lock().unwrap();
        for (i, vv) in st.velo.iter_mut().enumerate().take(dim - 1) {
            *vv = unsafe {
                var_value_at_idx(h_flow, (dim as c_int - 1) * (fp - 1) + 1 + i as c_int)
            };
        }
        let vv = st.velo.clone();
        drop(st);
        ut = unsafe { sum_pow2(vv.as_ptr(), dim - 1) }.sqrt();
        un = 0.0;
    }
    let _ = un;
    let zab = if use_floatation {
        let h_h = unsafe { variable_get(h_vars, "Depth", true, false) };
        let hp = unsafe { var_perm_at(h_h, node) };
        let h = unsafe { var_value_at_idx(h_h, hp) };
        BUDD_STATE.lock().unwrap().h = h;
        let rhow = unsafe {
            let r = const_real(h_bc, "Budd Ocean Density", &mut found);
            if !found {
                fatal("Sliding_Budd", "Need Ocean Density for the Budd sliding law");
            }
            r
        };
        let sl = unsafe {
            let s = get_creal(h_parent_mat, "Sea level", &mut found);
            if !found {
                fatal("Sliding_Budd", "Material property Sea level not found.");
            }
            s
        };
        let mut zab = if z < sl { h - (sl - z) * rhow / rhoi } else { h };
        if zab < 0.0 {
            zab = 0.0;
        }
        zab + zab_offset
    } else {
        let ep = unsafe { effective_pressure(m, n, &z as *const c_double) };
        -ep / (g * rhoi)
    };
    ut = ut.max(ut0);
    // `H` below is the stale static: set only in the floatation branch.
    let h_stale = BUDD_STATE.lock().unwrap().h;
    let mut zab = zab;
    if h_scaling {
        zab /= h_stale;
    }
    let mut c = c;
    if convert {
        let h_wc = unsafe { variable_get(h_vars, &weert_name, true, false) };
        let wp = unsafe { var_perm_at(h_wc, node) };
        let mut cw = unsafe { var_value_at_idx(h_wc, wp) };
        // Full-padded SELECT CASE on WeertForm.
        if weert_form == "power" {
            cw = unsafe { usf_pow_dp(10.0, cw) };
        } else if weert_form == "beta2" {
            cw *= cw;
        } else if weert_form == "none" {
        } else {
            unsafe { fatal("Sliding_Budd", "Unrecognised >Budd Conv Weertman formulation<") };
        }
        if weert_exp != 1.0 || m_exp != 1.0 {
            unsafe { fatal("Sliding_Budd", "Currently only works for velocity exponents = 1.0") };
        }
        c = cw / unsafe { usf_pow_dp(zab, q) };
    }
    c * unsafe { usf_pow_dp(ut, m_exp - 1.0) } * unsafe { usf_pow_dp(zab, q) }
}

// -- EffectivePressure ------------------------------------------------------

struct EPressState {
    first: bool,
    dim: c_int,
    ind: [[c_int; 3]; 3],
    sig: Vec<c_double>,
    normal: Vec<c_double>,
    sn: Vec<c_double>,
}

static EPRESS_STATE: Mutex<EPressState> = Mutex::new(EPressState {
    first: true,
    dim: 0,
    ind: [[0; 3]; 3],
    sig: Vec::new(),
    normal: Vec::new(),
    sn: Vec::new(),
});

pub unsafe fn effective_pressure(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _y = unsafe { *v };
    unsafe { usf_begin(m) };
    {
        let mut st = EPRESS_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            let msg = unsafe { write_msg_star_s("Calculating Effective Pressure") };
            unsafe { info("EffectivePressure", &msg, 5) };
            st.dim = unsafe { usf_coord_dim() };
            let d = st.dim;
            if d == 2 || d == 3 {
                st.sig = vec![0.0; (d * d) as usize];
                st.normal = vec![0.0; d as usize];
                st.sn = vec![0.0; d as usize];
            } else {
                unsafe { fatal("EffectivePressure", "Bad dimension of the problem") };
            }
            st.ind[0][0] = 1;
            st.ind[1][1] = 2;
            st.ind[2][2] = 3;
            st.ind[0][1] = 4;
            st.ind[1][0] = 4;
            st.ind[1][2] = 5;
            st.ind[2][1] = 5;
            st.ind[0][2] = 6;
            st.ind[2][0] = 6;
        }
    }
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_bc: *const c_void = std::ptr::null();
    unsafe { usf_get_bc(h_elem, &mut h_bc) };
    if h_bc.is_null() {
        unsafe { fatal("EffectivePressure", "No BC Found") };
    }
    let nn = unsafe { usf_get_element_nof_nodes(h_elem) } as usize;
    let mut idx = vec![0 as c_int; nn];
    unsafe { usf_element_nodeindexes(h_elem, idx.as_mut_ptr(), nn as c_int) };
    let mut aux = vec![0.0 as c_double; nn];
    let mut found = false;
    unsafe { get_real_into(h_bc, "External Pressure", &mut found, &mut aux) };
    let mut pos = nn; // EXIT-scan; `n+1` on no-match like the source.
    for (k, &id) in idx.iter().enumerate() {
        if id == node {
            pos = k;
            break;
        }
    }
    let pext = -aux[pos];
    drop(aux);
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_normal = unsafe { variable_get(h_vars, "Normal Vector", true, false) };
    let h_stress = unsafe { variable_get(h_vars, "Stress", true, false) };
    let h_flow = unsafe { variable_get(h_vars, "Flow Solution", true, false) };
    let np = unsafe { var_perm_at(h_normal, node) };
    let spp = unsafe { var_perm_at(h_stress, node) };
    let fp = unsafe { var_perm_at(h_flow, node) };
    let h_parent = unsafe { boundary_parent(h_elem) };
    let mut h_mat: *const c_void = std::ptr::null();
    unsafe { usf_get_material(h_parent, &mut h_mat) };
    let mut got = false;
    let mut cauchy = unsafe { list_get_logical(h_mat, "Cauchy", &mut got) };
    if !got {
        cauchy = false;
    }
    let mut st = EPRESS_STATE.lock().unwrap();
    let (dim, ind) = (st.dim as usize, st.ind);
    for i in 0..dim {
        for j in 0..dim {
            // SAFETY: column-major Sig(i,j) like the source: flat i + j*dim.
            st.sig[i + j * dim] = unsafe { var_value_at_idx(h_stress, 2 * dim as c_int * (spp - 1) + ind[i][j]) };
        }
        if !cauchy {
            let f = unsafe { var_value_at_idx(h_flow, (dim as c_int + 1) * fp) };
            st.sig[i + i * dim] -= f;
        }
    }
    for (i, nrm) in st.normal.iter_mut().enumerate().take(dim) {
        // SAFETY: normal carries the source's leading minus.
        *nrm =
            -unsafe { var_value_at_idx(h_normal, dim as c_int * (np - 1) + 1 + i as c_int) };
    }
    for i in 0..dim {
        let mut row = [0.0 as c_double; 3];
        for j in 0..dim {
            row[j] = st.sig[i + j * dim];
        }
        st.sn[i] = unsafe { dot(row.as_ptr(), st.normal.as_ptr(), dim) };
    }
    let snn = unsafe { dot(st.sn.as_ptr(), st.normal.as_ptr(), dim) };
    drop(st);
    -snn - pext
}

// -- FreeSlipShelves --------------------------------------------------------

struct ShelfState {
    first: bool,
    mask_name: String,
    beta_form: String,
}

static SHELF_STATE: Mutex<ShelfState> = Mutex::new(ShelfState {
    first: true,
    mask_name: String::new(),
    beta_form: String::new(),
});

pub unsafe fn free_slip_shelves(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let beta_in = unsafe { *v };
    unsafe { usf_begin(m) };
    {
        let mut st = SHELF_STATE.lock().unwrap();
        if st.first {
            let mut h_elem: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_elem) };
            let mut h_bc: *const c_void = std::ptr::null();
            unsafe { usf_get_bc(h_elem, &mut h_bc) };
            if h_bc.is_null() {
                unsafe { fatal("FreeSlipShelves", "No BC Found") };
            }
            let mut got = false;
            let mut mn = unsafe { get_string(h_bc, "FreeSlipShelves mask name", &mut got) };
            if !got {
                unsafe { warn("FreeSlipShelves", "Keyword >FreeSlipShelves mask name< not found in boundary condition") };
                unsafe { warn("FreeSlipShelves", "Taking default value >GroundedMask<") };
                mn = "GroundedMask".to_string();
            }
            st.mask_name = pad128(&mn);
            let mut got2 = false;
            let bf = unsafe { get_string(h_bc, "FreeSlipShelves beta formulation", &mut got2) };
            if !got2 {
                unsafe {
                    fatal(
                        "FreeSlipShelves",
                        "Need >FreeSlipShelves beta formulation< (e.g. Power or Beta2)",
                    )
                };
            }
            st.beta_form = pad128(&bf);
            st.first = false;
        }
    }
    let st = SHELF_STATE.lock().unwrap();
    let (mask_name, beta_form) = (st.mask_name.trim_end().to_string(), st.beta_form.trim_end().to_string());
    drop(st);
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    // `UnFoundFatal=.TRUE.` — missing mask variable Fatals inside.
    let h_mask = unsafe { variable_get(h_vars, &mask_name, true, false) };
    let mask = unsafe { var_value_at(h_mask, node) };
    if mask < 0.0 {
        0.0
    } else if beta_form == "Power" || beta_form == "power" {
        unsafe { usf_pow_dp(10.0, beta_in) }
    } else if beta_form == "Beta2" || beta_form == "beta2" {
        beta_in * beta_in
    } else {
        // `fatal` diverges (`!`), so this arm already yields the
        // function result type — no trailing value needed.
        unsafe { fatal("FreeSlipShelves", "beta formulation not recognised") };
    }
}
