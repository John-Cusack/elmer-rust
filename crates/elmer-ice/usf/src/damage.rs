//! `enhancementfactor_`, `sourcedamage_` — port of USF_Damage.F90
//! (continuum damage: enhancement factor + damage source).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated: `EnhancementFactor` keeps `FirstTime` + `n`
//!   (REAL, like the source); `SourceDamage` keeps `Ind(3,3)` (Voigt map),
//!   `DIM`, `sigmath`, `B`, `FirstTime`, `Cauchy` (BSS-false init).
//! - `n` is a REAL local in the source (assigned from `GetConstReal`);
//!   dead declarations (`Solver`, `tmp`) are documented, not replicated.
//! - `WRITE(Message, ...)` F10.4/L1 formatting evaluates in the Fortran
//!   accessor (`usf_write_msg_af104`/`usf_write_msg_al1`); INFO levels and
//!   call sites match exactly.
//! - `EvenRandom()` calls the in-process libgfortran stream through the
//!   accessor — Elmer seeds it deterministically at startup
//!   (`InitializeRandomSeed`), so both arms draw the identical sequence
//!   in the identical call order (`u` then `v`).
//! - `Sig` is column-major `[f64; 9]`; `DGEEV` (fixed 3x3, `Work(24)`)
//!   and `SortD` evaluate in Fortran wrappers; `Work`/`Dumy` scratch are
//!   zeroed statics (LAPACK write-before-read scratch — proven by the
//!   port-vs-port rerun in the proof script).
//! - The `Chi` write-through (`ChiValues(ChiPerm(n)) = Chi`) uses the
//!   live variable address — visible to Elmer by construction. This is
//!   the entry the SaveData-style visibility watch (kill criterion 2)
//!   is staked on; the staged Damage differential proves it.
//! - `SigDev` is replicated explicitly even though `DGEEV` consumes only
//!   `Sig` (mirrors the source, zero cleverness).
//! - The `ELSE IF (FirstTime)` PSeaD WARN arm is dead in the source
//!   (`FirstTime` is already `.FALSE.` there); kept verbatim.

use crate::support::*;

struct EnhState {
    first: bool,
    n: c_double,
}

static ENH_STATE: Mutex<EnhState> = Mutex::new(EnhState { first: true, n: 0.0 });

pub unsafe fn enhancement_factor(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let d = unsafe { *v };
    let _node = unsafe { *n };
    unsafe { usf_begin(m) };
    {
        let mut st = ENH_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            let mut h_mat: *const c_void = std::ptr::null();
            unsafe { usf_get_material_current(&mut h_mat) };
            let mut got = false;
            let g = unsafe { const_real(h_mat, "Glen Exponent", &mut got) };
            if !got {
                unsafe {
                    info(
                        "Damage EnhancementFactor",
                        "Variable Glen Exponent not found. Setting to 3.0",
                        2,
                    )
                };
                st.n = 3.0;
            } else {
                st.n = g;
                unsafe {
                    info(
                        "Damage EnhancementFactor",
                        &write_msg_af104("n = ", g),
                        2,
                    )
                };
            }
        }
    }
    let n = ENH_STATE.lock().unwrap().n;
    unsafe { usf_pow_dp(1.0 - d, -n) }
}

struct DmgState {
    first: bool,
    ind: [[c_int; 3]; 3],
    dim: c_int,
    sigmath: c_double,
    b: c_double,
    cauchy: bool,
}

static DMG_STATE: Mutex<DmgState> = Mutex::new(DmgState {
    first: true,
    ind: [[0; 3]; 3],
    dim: 0,
    sigmath: 0.0,
    b: 0.0,
    cauchy: false,
});

pub unsafe fn source_damage(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let d = unsafe { *v };
    let node = unsafe { *n };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    {
        let mut st = DMG_STATE.lock().unwrap();
        if st.first {
            st.dim = unsafe { usf_coord_dim() };
            st.ind[0][0] = 1;
            st.ind[1][1] = 2;
            st.ind[2][2] = 3;
            st.ind[0][1] = 4;
            st.ind[1][0] = 4;
            st.ind[1][2] = 5;
            st.ind[2][1] = 5;
            st.ind[2][0] = 6;
            st.ind[0][2] = 6;
            let mut h_mat: *const c_void = std::ptr::null();
            unsafe { usf_get_material_current(&mut h_mat) };
            let mut got = false;
            st.b = unsafe { const_real(h_mat, "Damage Enhancement Factor", &mut got) };
            if !got {
                unsafe { fatal("Damage Source", "Damage Enhancement Factor B not Found") };
            } else {
                unsafe {
                    info(
                        "Damage Source",
                        &write_msg_af104("Damage Enhancement Factor = ", st.b),
                        2,
                    )
                };
            }
            st.sigmath = unsafe { const_real(h_mat, "Damage Parameter sigmath", &mut got) };
            if !got {
                unsafe { fatal("Damage Source", "Damage Parameter Sigmath not Found") };
            } else {
                unsafe {
                    info(
                        "Damage Source",
                        &write_msg_af104("Damage Parameter sigmath = ", st.sigmath),
                        2,
                    )
                };
            }
            st.cauchy = unsafe { list_get_logical(h_mat, "Cauchy", &mut got) };
            unsafe {
                info(
                    "Damage Source",
                    &write_msg_al1("Cauchy stress tensor computed ? ", st.cauchy),
                    2,
                )
            };
        }
    }
    let h_stress = unsafe { variable_get(h_vars, "Stress", true, false) };
    let h_chi = unsafe { variable_get(h_vars, "Chi", true, false) };
    let h_psead = unsafe { variable_get(h_vars, "PSeaD", false, false) };
    if h_psead.is_null() {
        // Dead in the source (`FirstTime` already `.FALSE.`); kept verbatim.
        let first = DMG_STATE.lock().unwrap().first;
        if first {
            unsafe { warn("Damage Source", "PSeaD not associated, basal pressure not taken into account in damage formation") };
            unsafe { warn("Damage Source", "Taking default value PSeaD=0.0") };
        }
    }
    let h_flow = unsafe { variable_get(h_vars, "Flow Solution", true, false) };
    {
        DMG_STATE.lock().unwrap().first = false;
    }
    let st = DMG_STATE.lock().unwrap();
    let (dim, sigmath, b, cauchy) = (st.dim as usize, st.sigmath, st.b, st.cauchy);
    let ind = st.ind;
    drop(st);
    // Column-major Sig/SigDev/EigVect, like the source's (3,3) arrays.
    let mut sig = [0.0 as c_double; 9];
    let mut sigdev = [0.0 as c_double; 9];
    let sp = unsafe { var_perm_at(h_stress, node) };
    for i in 0..dim {
        for j in 0..dim {
            // SAFETY: same 1-based slot arithmetic the source uses.
            let slot = 2 * dim as c_int * (sp - 1) + ind[i][j];
            sig[j * 3 + i] = unsafe { var_value_at_idx(h_stress, slot) };
        }
    }
    if dim == 2 {
        let slot = 2 * dim as c_int * (sp - 1) + ind[2][2];
        sig[2 * 3 + 2] = unsafe { var_value_at_idx(h_stress, slot) };
    }
    let fp = unsafe { var_perm_at(h_flow, node) };
    // SAFETY: `(DIM+1)*FlowPerm` 1-based slot, verbatim.
    let fval = unsafe { var_value_at_idx(h_flow, (dim as c_int + 1) * fp) };
    if !cauchy {
        sigdev.copy_from_slice(&sig);
        for i in 0..3 {
            sig[i * 3 + i] -= fval;
        }
    } else {
        for i in 0..3 {
            sigdev[i * 3 + i] = sig[i * 3 + i] + fval;
        }
    }
    let mut eig = [0.0 as c_double; 3];
    let mut ei = [0.0 as c_double; 3];
    let mut work = [0.0 as c_double; 24];
    let mut dumy = [0.0 as c_double; 1];
    let mut infor: c_int = 0;
    unsafe {
        usf_dgeev(
            0,
            sig.as_mut_ptr(),
            eig.as_mut_ptr(),
            ei.as_mut_ptr(),
            dumy.as_mut_ptr(),
            1,
            work.as_mut_ptr(),
            &mut infor,
        )
    };
    if infor != 0 {
        unsafe { fatal("Compute EigenValues", "Failed to compute EigenValues") };
    }
    let mut eigvect = [0.0 as c_double; 9];
    unsafe {
        usf_dgeev(
            1,
            sig.as_mut_ptr(),
            eig.as_mut_ptr(),
            ei.as_mut_ptr(),
            eigvect.as_mut_ptr(),
            3,
            work.as_mut_ptr(),
            &mut infor,
        )
    };
    if infor != 0 {
        unsafe { fatal("Compute EigenVectors", "Failed to compute EigenVectors") };
    }
    let _ = eigvect;
    let mut indice = [1 as c_int, 2, 3];
    unsafe { usf_sortd(3, eig.as_mut_ptr(), indice.as_mut_ptr()) };
    let sigma_i = eig[2];
    let sigma_ii = eig[1];
    let _ = sigma_ii;
    // N.B. NOT `PI`: the source literal is `3.141592_dp` (7 digits), and
    // the value feeds `COS(2*nbrPi*v)` — full-precision PI would diverge.
    #[allow(clippy::approx_constant)]
    let nbr_pi: c_double = 3.141592;
    let u = unsafe { even_random() };
    let vv = unsafe { even_random() };
    let mut h_c: *const c_void = std::ptr::null();
    unsafe { usf_get_constants(&mut h_c) };
    let mut got = false;
    let s = unsafe { const_real(h_c, "Dev Tensile Strength Modifier", &mut got) };
    if !got {
        unsafe {
            fatal(
                "USF_Damage",
                "No \"Dev tensile strength modifier\" given, set to 0.05",
            )
        };
    }
    // N.B. the Fatal above diverts, so `s` is always defined here —
    // exactly like the source (whose `s` would otherwise be stale).
    let sigmath_var = (0.0 + s * (-2.0 * (1.0 - u).ln()).sqrt() * (2.0 * nbr_pi * vv).cos()).abs();
    let stress_threshold = sigmath * (1.0 + sigmath_var);
    let pwater = if !h_psead.is_null() {
        unsafe { var_value_at_idx(h_psead, var_perm_at(h_psead, node)) }
    } else {
        0.0
    };
    let chi = 1.0 / (1.0 - d) * (sigma_i + pwater) - stress_threshold;
    unsafe { var_set_at(h_chi, node, chi) };
    b * chi.max(0.0)
}
