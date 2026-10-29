//! `caffegetviscosity_` — port of CaffeFlow.F90 (anisotropic CAFFE
//! ice-flow viscosity factor).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated: `ElementNodes` scratch, `DIM`, `FirstTime`,
//!   `gasconst`, the `(2,nMax)` coefficient stores
//!   (`activationEnergy`, `arrheniusFactor`, `aniEnhancementFact`,
//!   `viscosityExponent`, `PressureMeltingPoint`, `LimitTemp`) sized by
//!   `MaxElementNodes`, `Hwrk` dims. `stat`/`istat` are plain locals.
//! - `Hwrk` (a Fortran 3-D POINTER re-pointed per call) is copied into a
//!   Rust `[3 x n]` store by the accessor; `SIZE(Hwrk,2) == 1` and the
//!   `MIN(3,SIZE(Hwrk,1))` row loop are replicated verbatim.
//! - `stat = ElementInfo(...)` is assigned and never checked in the
//!   source — the call is kept (side effects on `Basis`/`dBasisdx`),
//!   the value dropped, documented.
//! - `LGrad = MATMUL(...)`, `Dij = 0.5*(LGrad+TRANSPOSE)` evaluate in a
//!   Fortran wrapper (intrinsic evaluation order); `**` with a real
//!   exponent goes through `usf_pow_dp`, integer `**2` is multiplication.
//! - `Eps > 0 .OR. Eps < 0` is written out (not `!=`): NaN skips the
//!   adjustment in both arms.
//! - `1.0D00/3.0D00`, `2.7315D02` (= 273.15) and all decimal literals
//!   use the identical values the source text denotes.
//! - `WRITE(Message,...)` with explicit `(a,I2,...)`, `(a,e10.4,a)`,
//!   `(A,A,A,i3,A)` formats evaluates in Fortran accessors; plain
//!   `(A)` literals pass through (`Fatal`/`Info` TRIM).
//! - The EXIT-scan over the element's nodes leaves `n+1` on no-match
//!   (like the source's completed DO loop); the panicking index there
//!   marks the invariant (node in element) the reference would violate
//!   with an out-of-bounds read.
//! - `IBOF` (internal closure) is a pure transcription with identical
//!   literal text (rewritten `e`-notation only) and identical
//!   accumulation order.

use crate::support::*;

struct CaffeState {
    first: bool,
    dim: c_int,
    gasconst: c_double,
    activation_energy: Vec<c_double>,
    arrhenius_factor: Vec<c_double>,
    ani_enhance: Vec<c_double>,
    viscosity_exponent: Vec<c_double>,
    pressure_melting_point: Vec<c_double>,
    limit_temp: Vec<c_double>,
    maxnodes: usize,
    // Last computed value (RESULT local: gfortran static storage).
    last: c_double,
}

static CAFFE_STATE: Mutex<CaffeState> = Mutex::new(CaffeState {
    first: true,
    dim: 0,
    gasconst: 0.0,
    activation_energy: Vec::new(),
    arrhenius_factor: Vec::new(),
    ani_enhance: Vec::new(),
    viscosity_exponent: Vec::new(),
    pressure_melting_point: Vec::new(),
    limit_temp: Vec::new(),
    maxnodes: 0,
    last: 0.0,
});

/// `IBOF(a2)` closure (Chung 2002): 6-vector to 9-vector. Pure
/// transcription; accumulation order and literals mirror the source.
fn ibof(a2: &[c_double; 6]) -> [c_double; 9] {
    let (a11, a22, a33, a12, a23, a13) = (a2[0], a2[1], a2[2], a2[3], a2[4], a2[5]);
    #[rustfmt::skip]
    let mat: [[c_double; 21]; 3] = [
        [0.217774509809788e+02, -0.297570854171128e+03, 0.188686077307885e+04, -0.272941724578513e+03, 0.417148493642195e+03, 0.152038182241196e+04, -0.137643852992708e+04, -0.628895857556395e+03, -0.526081007711996e+04, -0.266096234984017e+03, -0.196278098216953e+04, -0.505266963449819e+03, -0.110483041928547e+03, 0.430488193758786e+04, -0.139197970442470e+02, -0.144351781922013e+04, -0.265701301773249e+03, -0.428821699139210e+02, -0.443236656693991e+01, 0.309742340203200e+04, 0.386473912295113e+00],
        [-0.514850598717222e+00, 0.213316362570669e+02, -0.302865564916568e+03, -0.198569416607029e+02, -0.460306750911640e+02, 0.270825710321281e+01, 0.184510695601404e+03, 0.156537424620061e+03, 0.190613131168980e+04, 0.277006550460850e+03, -0.568117055198608e+02, 0.428921546783467e+03, 0.142494945404341e+03, -0.541945228489881e+04, 0.233351898912768e+02, 0.104183218654671e+04, 0.331489412844667e+03, 0.660002154209991e+02, 0.997500770521877e+01, 0.560508628472486e+04, 0.209909225990756e+01],
        [0.203814051719994e+02, -0.283958093739548e+03, 0.173908241235198e+04, -0.195566197110461e+03, -0.138012943339611e+03, 0.523629892715050e+03, 0.859266451736379e+03, -0.805606471979730e+02, -0.468711180560599e+04, 0.889580760829066e+01, -0.782994158054881e+02, -0.437214580089117e+02, 0.112996386047623e+01, 0.401746416262936e+04, 0.104927789918320e+01, -0.139340154288711e+03, -0.170995948015951e+02, 0.545784716783902e+00, 0.971126767581517e+00, 0.141909512967882e+04, 0.994142892628410e+00],
    ];
    let inv2 = 0.5 * (1.0 - (a11 * a11 + a22 * a22 + a33 * a33 + 2.0 * (a12 * a12 + a13 * a13 + a23 * a23)));
    let inv3 = a11 * (a22 * a33 - a23 * a23) + a12 * (a23 * a13 - a12 * a33) + a13 * (a12 * a23 - a22 * a13);
    let mut vec = [0.0 as c_double; 21];
    vec[0] = 1.0;
    vec[1] = inv2;
    vec[2] = vec[1] * vec[1];
    vec[3] = inv3;
    vec[4] = vec[3] * vec[3];
    vec[5] = vec[1] * vec[3];
    vec[6] = vec[2] * vec[3];
    vec[7] = vec[1] * vec[4];
    vec[8] = vec[1] * vec[2];
    vec[9] = vec[4] * vec[3];
    vec[10] = vec[8] * vec[3];
    vec[11] = vec[2] * vec[4];
    vec[12] = vec[1] * vec[9];
    vec[13] = vec[2] * vec[2];
    vec[14] = vec[4] * vec[4];
    vec[15] = vec[13] * vec[3];
    vec[16] = vec[11] * vec[1];
    vec[17] = vec[11] * vec[3];
    vec[18] = vec[1] * vec[14];
    vec[19] = vec[13] * vec[1];
    vec[20] = vec[14] * vec[3];
    let mut beta = [0.0 as c_double; 6];
    for i in 0..3 {
        for j in 0..21 {
            beta[i] += mat[i][j] * vec[j];
        }
    }
    beta[3] = 3.0 * (-1.0 / 7.0 + beta[0] * (1.0 / 7.0 + 4.0 * inv2 / 7.0 + 8.0 * inv3 / 3.0) / 5.0
        - beta[1] * (0.2 - 8.0 * inv2 / 15.0 - 14.0 * inv3 / 15.0)
        - beta[2] * (1.0 / 35.0 - 24.0 * inv3 / 105.0 - 4.0 * inv2 / 35.0 + 16.0 * inv2 * inv3 / 15.0 + 8.0 * inv2 * inv2 / 35.0))
        / 5.0;
    beta[4] = 6.0
        * (1.0 - 0.2 * beta[0] * (1.0 + 4.0 * inv2)
            + 7.0 * beta[1] * (1.0 / 6.0 - inv2) / 5.0
            - beta[2] * (-0.2 + 2.0 * inv3 / 3.0 + 4.0 * inv2 / 5.0 - 8.0 * inv2 * inv2 / 5.0))
        / 7.0;
    beta[5] =
        -4.0 * beta[0] / 5.0 - 7.0 * beta[1] / 5.0 - 6.0 * beta[2] * (1.0 - 4.0 * inv2 / 3.0) / 5.0;
    for b in beta.iter_mut() {
        *b /= 3.0;
    }
    beta[1] /= 2.0;
    beta[4] /= 2.0;
    beta[5] /= 2.0;
    let b11 = a11 * a11 + a12 * a12 + a13 * a13;
    let b22 = a22 * a22 + a12 * a12 + a23 * a23;
    let b12 = a11 * a12 + a12 * a22 + a13 * a23;
    let b13 = a11 * a13 + a12 * a23 + a13 * a33;
    let b23 = a12 * a13 + a22 * a23 + a23 * a33;
    let mut a4 = [0.0 as c_double; 9];
    a4[0] = 3.0 * beta[3] + 6.0 * beta[4] * a11 + 3.0 * beta[0] * a11 * a11 + 6.0 * beta[1] * b11
        + 6.0 * beta[5] * a11 * b11
        + 3.0 * beta[2] * b11 * b11;
    a4[1] = 3.0 * beta[3] + 6.0 * beta[4] * a22 + 3.0 * beta[0] * a22 * a22 + 6.0 * beta[1] * b22
        + 6.0 * beta[5] * a22 * b22
        + 3.0 * beta[2] * b22 * b22;
    a4[2] = beta[3] + beta[4] * (a22 + a11) + beta[0] * (a11 * a22 + 2.0 * a12 * a12)
        + beta[1] * (b22 + b11)
        + beta[5] * (a11 * b22 + a22 * b11 + 4.0 * a12 * b12)
        + beta[2] * (b11 * b22 + 2.0 * b12 * b12);
    a4[3] = beta[4] * a23 + beta[0] * (a11 * a23 + 2.0 * a12 * a13) + beta[1] * b23
        + beta[5] * (a11 * b23 + a23 * b11 + 2.0 * (a12 * b13 + a13 * b12))
        + beta[2] * (b11 * b23 + 2.0 * b12 * b13);
    a4[4] = beta[4] * a13 + beta[0] * (a22 * a13 + 2.0 * a12 * a23) + beta[1] * b13
        + beta[5] * (a22 * b13 + a13 * b22 + 2.0 * (a12 * b23 + a23 * b12))
        + beta[2] * (b22 * b13 + 2.0 * b12 * b23);
    a4[5] = 3.0 * beta[4] * a13 + 3.0 * beta[0] * a11 * a13 + 3.0 * beta[1] * b13
        + 3.0 * beta[5] * (a11 * b13 + a13 * b11)
        + 3.0 * beta[2] * b11 * b13;
    a4[6] = 3.0 * beta[4] * a12 + 3.0 * beta[0] * a11 * a12 + 3.0 * beta[1] * b12
        + 3.0 * beta[5] * (a11 * b12 + a12 * b11)
        + 3.0 * beta[2] * b11 * b12;
    a4[7] = 3.0 * beta[4] * a23 + 3.0 * beta[0] * a22 * a23 + 3.0 * beta[1] * b23
        + 3.0 * beta[5] * (a22 * b23 + a23 * b22)
        + 3.0 * beta[2] * b22 * b23;
    a4[8] = 3.0 * beta[4] * a12 + 3.0 * beta[0] * a22 * a12 + 3.0 * beta[1] * b12
        + 3.0 * beta[5] * (a22 * b12 + a12 * b22)
        + 3.0 * beta[2] * b22 * b12;
    a4
}

pub unsafe fn caffe_get_viscosity(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let temperature = unsafe { *v };
    unsafe { usf_begin(m) };
    {
        let mut st = CAFFE_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            st.dim = unsafe { usf_coord_dim() };
            let maxnodes = unsafe { usf_max_element_nodes() } as usize;
            st.maxnodes = maxnodes;
            let mut h_c: *const c_void = std::ptr::null();
            unsafe { usf_get_constants(&mut h_c) };
            let mut got = false;
            st.gasconst = unsafe { list_get_const_real_found(h_c, "Gas Constant", &mut got) };
            if !got {
                st.gasconst = 8.314;
                let msg = unsafe {
                    write_msg_aea104a(
                        "No entry for Gas Constant (Constants) in input file found. Setting to ",
                        st.gasconst,
                        " (J/mol)",
                    )
                };
                unsafe { info("CAFFE (CAFFEViscosity)", &msg, 4) };
            }
            let nmax = maxnodes;
            st.activation_energy = vec![0.0; 2 * nmax];
            st.arrhenius_factor = vec![0.0; 2 * nmax];
            st.limit_temp = vec![0.0; nmax];
            st.ani_enhance = vec![0.0; nmax];
            st.pressure_melting_point = vec![0.0; nmax];
            st.viscosity_exponent = vec![0.0; nmax];
            unsafe { info("CAFFE (CAFFEViscosity)", "Memory allocations done", 3) };
        }
    }
    let h_svars = unsafe { current_solver_mesh_variables() };
    let h_flow = unsafe { variable_get(h_svars, "flow solution", false, false) };
    if h_flow.is_null() {
        unsafe { info("CAFFE", "No variable for velocity associated.", 4) };
    }
    // Null handles fault identically below (pass-through, never guarded).
    let stdofs = unsafe { variable_dofs(h_flow) } as usize;
    let h_fabric = unsafe { variable_get(h_svars, "fabric", false, false) };
    if h_fabric.is_null() {
        unsafe { info("CAFFE", "No variable for fabric associated.", 4) };
    }
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let body_id = unsafe { usf_element_body_id(h_elem) };
    let mut h_b1: *const c_void = std::ptr::null();
    unsafe { usf_model_body_values(body_id, &mut h_b1) };
    let mut got = false;
    let material_id = unsafe { list_get_integer(h_b1, "Material", true, &mut got) };
    if !got {
        unsafe { fatal("CAFFE (CAFFEViscosity)", "No Material ID found") };
    }
    let element_nb_nodes = unsafe { usf_get_element_nof_nodes(h_elem) } as usize;
    // `IF (.NOT. GotIt)` reuses the material-lookup `Found` — dead when
    // the lookup failed (Fatal above), always-false here; kept verbatim.
    if !got {
        let msg = unsafe { write_msg_ai2ai2a("No material id for current element of node ", node, ", body ", body_id, " found") };
        unsafe { fatal("CAFFE (CAFFEViscosity)", &msg) };
    }
    let mut idx = vec![0 as c_int; element_nb_nodes];
    unsafe { usf_element_nodeindexes(h_elem, idx.as_mut_ptr(), element_nb_nodes as c_int) };
    let mut h_mat: *const c_void = std::ptr::null();
    unsafe { usf_model_material_values(material_id, &mut h_mat) };
    if h_mat.is_null() {
        let msg = unsafe { write_msg_ai2ai2a("No Material for current element of node ", node, ", body ", body_id, " found") };
        unsafe { fatal("CAFFE (CAFFEViscosity)", &msg) };
    }
    let maxnodes = CAFFE_STATE.lock().unwrap().maxnodes;
    // `ElementNodes` scratch + `LocalVelo`/`Basis`/`dBasisdx` automatics.
    let mut enodes_x = vec![0.0 as c_double; maxnodes];
    let mut enodes_y = vec![0.0 as c_double; maxnodes];
    let mut enodes_z = vec![0.0 as c_double; maxnodes];
    for k in 0..element_nb_nodes {
        let j = idx[k] as usize;
        enodes_x[k] = unsafe { model_node_x(j as c_int) };
        enodes_y[k] = unsafe { model_node_y(j as c_int) };
        enodes_z[k] = unsafe { model_node_z(j as c_int) };
    }
    // LocalVelo(3, MaxElementNodes), zeroed; rows 1..STDOFs-1 filled.
    let mut local_velo = vec![0.0 as c_double; 3 * maxnodes];
    for r in 0..stdofs.saturating_sub(1) {
        for k in 0..element_nb_nodes {
            let j = idx[k];
            let p = unsafe { var_perm_at(h_flow, j) };
            // SAFETY: same 1-based slot arithmetic the source uses.
            local_velo[r + 3 * k] = unsafe { var_value_at_idx(h_flow, stdofs as c_int * (p - 1) + 1 + r as c_int) };
        }
    }
    // EXIT-scan over the element's nodes for `nodenumber`.
    let mut pos = element_nb_nodes + 1;
    for (k, &id) in idx.iter().enumerate() {
        if node == id {
            pos = k;
            break;
        }
    }
    if pos >= element_nb_nodes {
        // No match: the RESULT local keeps its static storage.
        return CAFFE_STATE.lock().unwrap().last;
    }
    let i = pos; // 0-based element-node position (source `i`, 1-based).
    let mut st = CAFFE_STATE.lock().unwrap();
    // Activation energies.
    let mut buf = vec![0.0 as c_double; 3 * element_nb_nodes];
    let mut g2 = false;
    let (d1, d2) = unsafe { list_get_real_array_into(h_mat, "Activation Energies", &idx, &mut g2, &mut buf) };
    if !g2 {
        let msg = unsafe { write_msg_ai2ai2("No Value for Activation Energy  found in Material ", material_id, " for node ", node) };
        unsafe { fatal(" CAFFE (CAFFEViscosity)", &msg) };
    }
    if d2 == 1 {
        for j in 0..(3.min(d1 as usize)) {
            for k in 0..element_nb_nodes {
                st.activation_energy[j + 2 * k] = buf[j + 3 * k];
            }
        }
    } else {
        let msg = unsafe { write_msg_ai2ai2("Incorrect array size for Activation Energy in Material ", material_id, " for node ", node) };
        unsafe { fatal("CAFFE (CAFFEViscosity)", &msg) };
    }
    // Arrhenius factors.
    let mut g3 = false;
    let (e1, e2) = unsafe { list_get_real_array_into(h_mat, "Arrhenius Factors", &idx, &mut g3, &mut buf) };
    if !g3 {
        let msg = unsafe { write_msg_ai2ai2("No Value for Arrhenius Factors  found in Material ", material_id, " for node ", node) };
        unsafe { fatal("CAFFE (CAFFEViscosity)", &msg) };
    }
    if e2 == 1 {
        for j in 0..(3.min(e1 as usize)) {
            for k in 0..element_nb_nodes {
                st.arrhenius_factor[j + 2 * k] = buf[j + 3 * k];
            }
        }
    } else {
        let msg = unsafe { write_msg_ai2ai2("Incorrect array size for Arrhenius Factors in Material ", material_id, " for node ", node) };
        unsafe { fatal("CAFFE (CAFFEViscosity)", &msg) };
    }
    // Threshold temperature.
    let mut g4 = false;
    let mut tmp = vec![0.0 as c_double; element_nb_nodes];
    unsafe { list_get_real_into(h_mat, "Limit Temperature", &idx, false, &mut g4, &mut tmp) };
    st.limit_temp[..element_nb_nodes].copy_from_slice(&tmp);
    if !g4 {
        st.limit_temp[..element_nb_nodes].fill(-10.0);
        let msg = unsafe { write_msg_ai2ai2a("No keyword >Limit Temperature< found in Material ", material_id, " for node ", node, ".setting to -10") };
        unsafe { info("CAFFE (CAFFEViscosity)", &msg, 4) };
    }
    // Viscosity exponent.
    let mut g5 = false;
    unsafe { list_get_real_into(h_mat, "Viscosity Exponent", &idx, false, &mut g5, &mut tmp) };
    st.viscosity_exponent[..element_nb_nodes].copy_from_slice(&tmp);
    if !g5 {
        st.viscosity_exponent[..element_nb_nodes].fill(1.0 / 3.0);
        let msg = unsafe { write_msg_ai2ai2a("No Viscosity Exponent found in Material ", material_id, " for node ", node, ".setting k=1/3") };
        unsafe { info("CAFFE (CAFFEViscosity)", &msg, 4) };
    }
    // Anisotropic enhancement factor.
    let mut g6 = false;
    unsafe { list_get_real_into(h_mat, "Anisotropic Enhancement Factor", &idx, false, &mut g6, &mut tmp) };
    st.ani_enhance[..element_nb_nodes].copy_from_slice(&tmp);
    if !g6 {
        let msg = unsafe { write_msg_ai2ai2a("Value for Maximum Enhancement factor not found in Matrial", material_id, " for node ", node, "") };
        unsafe { fatal("CAFFE (CAFFEViscosity)", &msg) };
    }
    // Critical enhancement factor (singularities guard).
    let mut g7 = false;
    let mut emin = unsafe { const_real(h_mat, "Critical Enhancement factor", &mut g7) };
    if !g7 {
        emin = 0.0001;
        let msg = unsafe { write_msg_ai2ai2a("Value for Critical Enhancement factor not found in Matrial", material_id, " for node ", node, ".setting Emin=1.Oe-4") };
        unsafe { info("CAFFE (CAFFEViscosity)", &msg, 3) };
    }
    // Critical strain rate.
    let mut g8 = false;
    let mut min_gamma = unsafe { const_real(h_mat, "Critical Strain Rate", &mut g8) };
    if !g8 {
        min_gamma = 1.0e-15;
        let msg = unsafe { write_msg_ai2ai2a("Critical Strain rate not found in Matrial", material_id, " for node ", node, ".setting MinGamma=1.0e-15") };
        unsafe { info("CAFFE (CAFFEViscosity)", &msg, 3) };
    }
    // Temperature name + pressure melting point.
    let mut g9 = false;
    let temp_name = unsafe { get_string(h_mat, "Temperature Name", &mut g9) };
    if !g9 {
        unsafe { fatal("CAFFE (CAFFEViscosity)", "No Temperature Name found") };
    }
    let pkey = format!("{} Upper Limit", temp_name.trim_end());
    let mut g10 = false;
    unsafe { list_get_real_into(h_mat, &pkey, &idx, false, &mut g10, &mut tmp) };
    st.pressure_melting_point[..element_nb_nodes].copy_from_slice(&tmp);
    let temphom = if !g10 {
        let msg = unsafe { write_msg_aaai3a("No entry for ", &pkey, " found in material no. ", material_id, ". Using 273.16 K.") };
        unsafe { warn("iceproperties (getViscosityFactor)", &msg) };
        0.0
    } else {
        (temperature - st.pressure_melting_point[i]).min(0.0)
    };
    let k = if temphom < st.limit_temp[i] { 1 } else { 2 };
    // `EXP` lowers to libm `exp` — Rust `.exp()` is the same call.
    // `-(a/b)` is bit-identical to the source's `(-1.0*a)/b` (single
    // rounding of the same real value; negation is exact).
    let rate_factor = st.arrhenius_factor[(k - 1) + 2 * i]
        * (-(st.activation_energy[(k - 1) + 2 * i] / (st.gasconst * (2.7315e2 + temphom)))).exp();
    drop(st);
    // `u, v, w` local coordinates of node `i`.
    let u = unsafe { element_type_nodeu(h_elem, (i + 1) as c_int) };
    let v = unsafe { element_type_nodev(h_elem, (i + 1) as c_int) };
    let w = unsafe { element_type_nodew(h_elem, (i + 1) as c_int) };
    let mut basis = vec![0.0 as c_double; maxnodes];
    let mut dbasis = vec![0.0 as c_double; 3 * maxnodes];
    // `stat` is assigned and never checked in the source — kept verbatim.
    let (detj, _stat) = unsafe {
        element_info(
            h_elem,
            u,
            v,
            w,
            element_nb_nodes,
            &mut basis[..element_nb_nodes],
            &mut dbasis[..3 * element_nb_nodes],
        )
    };
    let _ = detj;
    let mut dij = [0.0 as c_double; 9];
    unsafe { lgrad_dij(&local_velo[..3 * element_nb_nodes], &dbasis[..3 * element_nb_nodes], element_nb_nodes, &mut dij) };
    // Fabric second-order tensor (null handle faults identically).
    let mut a2 = [0.0 as c_double; 6];
    if !h_fabric.is_null() {
        let fp = unsafe { var_perm_at(h_fabric, node) };
        a2[0] = unsafe { var_value_at_idx(h_fabric, 5 * (fp - 1) + 1) };
        a2[1] = unsafe { var_value_at_idx(h_fabric, 5 * (fp - 1) + 2) };
        a2[2] = 1.0 - a2[0] - a2[1];
        a2[3] = unsafe { var_value_at_idx(h_fabric, 5 * (fp - 1) + 3) };
        a2[4] = unsafe { var_value_at_idx(h_fabric, 5 * (fp - 1) + 4) };
        a2[5] = unsafe { var_value_at_idx(h_fabric, 5 * (fp - 1) + 5) };
    } else {
        a2[3] = 0.0;
        a2[4] = 0.0;
        a2[5] = 0.0;
        a2[0] = 1.0 / 3.0;
        a2[1] = 1.0 / 3.0;
        a2[2] = 1.0 / 3.0;
    }
    let a4 = ibof(&a2);
    let a1133 = -a4[0] + a2[0] - a4[2];
    let a2233 = -a4[1] + a2[1] - a4[2];
    let a3333 = a2[2] - a1133 - a2233;
    let a3331 = a2[5] - a4[5] - a4[4];
    let a3332 = a2[4] - a4[3] - a4[7];
    let a3312 = a2[3] - a4[6] - a4[8];
    let mut eps = dij[0] + dij[4] + dij[8];
    // N.B. NOT `eps != 0.0`: a NaN trace must skip the adjustment (both
    // comparisons false), while `!=` would apply it. Verbatim source form.
    #[allow(clippy::double_comparisons)]
    if eps > 0.0 || eps < 0.0 {
        eps /= 3.0;
        dij[0] -= eps;
        dij[4] -= eps;
        dij[8] -= eps;
    }
    // AfirstArg / AsecondArg — `**2.0` default-real exponents are `pow`.
    let d11 = dij[0];
    let d22 = dij[4];
    let d33 = dij[8];
    let d12 = dij[3];
    let d13 = dij[6];
    let d23 = dij[7];
    let p2 = unsafe { usf_pow_dp(d11, 2.0) };
    let _ = p2;
    let afirst = a2[0]
        * (unsafe { usf_pow_dp(d11, 2.0) } + unsafe { usf_pow_dp(d12, 2.0) } + unsafe { usf_pow_dp(d13, 2.0) })
        + a2[1]
            * (unsafe { usf_pow_dp(d12, 2.0) } + unsafe { usf_pow_dp(d22, 2.0) } + unsafe { usf_pow_dp(d23, 2.0) })
        + a2[2]
            * (unsafe { usf_pow_dp(d13, 2.0) } + unsafe { usf_pow_dp(d23, 2.0) } + unsafe { usf_pow_dp(d33, 2.0) })
        + a2[3] * (2.0 * d11 * d12 + 2.0 * d22 * d12 + 2.0 * d23 * d13)
        + a2[4] * (2.0 * d13 * d12 + 2.0 * d22 * d23 + 2.0 * d33 * d23)
        + a2[5] * (2.0 * d11 * d13 + 2.0 * d23 * d12 + 2.0 * d33 * d13);
    let asecond = a4[0] * unsafe { usf_pow_dp(d11, 2.0) } + a4[1] * unsafe { usf_pow_dp(d22, 2.0) } + a3333 * unsafe { usf_pow_dp(d33, 2.0) }
        + a4[2] * (2.0 * d11 * d22 + 4.0 * unsafe { usf_pow_dp(d12, 2.0) })
        + a1133 * (2.0 * d11 * d33 + 4.0 * unsafe { usf_pow_dp(d13, 2.0) })
        + a2233 * (2.0 * d22 * d33 + 4.0 * unsafe { usf_pow_dp(d23, 2.0) })
        + 4.0 * a4[3] * d11 * d23
        + 8.0 * a4[3] * d12 * d13
        + 4.0 * a4[4] * d22 * d13
        + 8.0 * a4[4] * d12 * d23
        + 4.0 * a3312 * d33 * d12
        + 8.0 * a3312 * d13 * d23
        + 4.0 * a4[6] * d11 * d12
        + 4.0 * a4[5] * d11 * d13
        + 4.0 * a4[8] * d22 * dij[3]
        + 4.0 * a4[7] * d22 * d23
        + 4.0 * a3331 * d33 * d13
        + 4.0 * a3332 * d33 * d23;
    let st = CAFFE_STATE.lock().unwrap();
    let mut aconst = unsafe { usf_pow_dp(d11, 2.0) }
        + unsafe { usf_pow_dp(d22, 2.0) }
        + unsafe { usf_pow_dp(d33, 2.0) }
        + 2.0 * unsafe { usf_pow_dp(d12, 2.0) }
        + 2.0 * unsafe { usf_pow_dp(d13, 2.0) }
        + 2.0 * unsafe { usf_pow_dp(d23, 2.0) };
    if aconst < min_gamma {
        aconst = min_gamma;
    }
    let a = 5.0 * (afirst - asecond) / aconst;
    let enh = st.ani_enhance[i];
    let expo = st.viscosity_exponent[i];
    drop(st);
    let ee = if a >= 1.0 {
        (4.0 * unsafe { usf_pow_dp(a, 2.0) } * (enh - 1.0) + 25.0 - 4.0 * enh) / 21.0
    } else {
        (1.0 - emin) * unsafe { usf_pow_dp(a, (8.0 / 21.0) * ((enh - 1.0) / (1.0 - emin))) } + emin
    };
    let visc = 1.0 / unsafe { usf_pow_dp((2.0 * ee) * rate_factor, expo) };
    CAFFE_STATE.lock().unwrap().last = visc;
    visc
}
