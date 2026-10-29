//! `seapressure_`, `seaspring_` — port of Buoyancy.F90 (sea pressure
//! and sea-spring viscous term for floating ice).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated per entry. Fortran implicit-SAVE rules apply:
//!   `NormalFlux` (init `.TRUE.`), `UnFoundFatal` (init `.TRUE.`),
//!   `MeshChanged` (init `.FALSE.`) keep their initializers as static
//!   inits; plain non-SAVE locals (`NewTime`, `S`, `Ns`, …) start BSS
//!   (false/zero) — `NewTime` is set `.TRUE.` in every First block
//!   before any read, so the BSS start is unobservable.
//! - Dead declarations (`Nn`, `Nmax`, `NBoundary`, `NParent`,
//!   `BoundaryElementNode`, `ParentElementNode`, unread `material`,
//!   SeaPressure's `BodyForce`) are documented, not replicated.
//! - `GotIt`-gated gravity reads use only the SECOND call's `Found`
//!   (source quirk), verbatim.
//! - `WRITE(Message,*)` list-directed output (including the `(A,I3)`
//!   WRITE with no integer item, `(A,e14.8)`, `(A,F10.4)`, and mixed
//!   items) evaluates in Fortran accessors.
//! - `WRITE(Message,'(A)')` + `INFO`/`FATAL` paths pass the literal:
//!   `Fatal`/`Info` TRIM trailing blanks, so the write-into-`Message`
//!   step is unobservable — except list-directed writes (leading
//!   blanks significant), which go through the accessor.
//! - `10.0_dp**x`, `(a/b)**2.0_dp` use real-exponent `pow` through
//!   `usf_pow_dp`; integer-exponent `**2` in SeaSpring's `Nsi` is
//!   explicit multiplication (`NodeNormal(1:dim-1)**2` sums squares).
//! - The `NodeOnBoundary` map is 1-based (0 = unmapped) over global
//!   node numbers; the Rust store is 0-based with the mapping shifted
//!   by one. An unmapped node (mapping 0) panics loudly where the
//!   reference reads out of bounds — the invariant (called for boundary
//!   nodes) always holds.
//! - `S`/`Ns`/`a_perp`/`SourceFunc`/`Normal` stores are column-major
//!   flat `Vec`s with the source's `(row, col)` indexing.
//! - `TYPE(Nodes_t), SAVE` scratch lives inside the `usf_normal_vector`
//!   accessor (re-filled per call — no staleness); `GetElementNodes`
//!   calls with unused-here results are kept verbatim.
//! - `Normal(DIM,i)` zero-check uses `ABS(...) > 1.0e-20` with the
//!   `-999.0` fallback and its per-node `WARN`, verbatim.

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

// -- SeaPressure --------------------------------------------------------

struct PressureState {
    first: bool,
    new_time: bool,
    told: c_double,
    dt: c_double,
    ns: Vec<c_double>,
    dim: c_int,
    s: Vec<c_double>,
    rhow: c_double,
    gravity: c_double,
    zsl: c_double,
    normal_flux: bool,
    a_perp: Vec<c_double>,
    source_func: Vec<c_double>,
    number_on_boundary: c_int,
    node_on_boundary: Vec<c_int>,
    normal: Vec<c_double>,
    bottom_surface_name: String,
    bf_id_fs: c_int,
    old_mesh_tag: c_int,
    mesh_changed: bool,
}

static PRESSURE_STATE: Mutex<PressureState> = Mutex::new(PressureState {
    first: true,
    new_time: false,
    told: 0.0,
    dt: 0.0,
    ns: Vec::new(),
    dim: 0,
    s: Vec::new(),
    rhow: 0.0,
    gravity: 0.0,
    zsl: 0.0,
    normal_flux: true,
    a_perp: Vec::new(),
    source_func: Vec::new(),
    number_on_boundary: 0,
    node_on_boundary: Vec::new(),
    normal: Vec::new(),
    bottom_surface_name: String::new(),
    bf_id_fs: 0,
    old_mesh_tag: 0,
    mesh_changed: false,
});

pub unsafe fn sea_pressure(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _y = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_time = unsafe { variable_get(h_vars, "Time", true, false) };
    let t = unsafe { var_value_at_idx(h_time, 1) };
    let dt_now = unsafe { solver_dt() };
    let mut h_elem0: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem0) };
    if unsafe { get_element_family(h_elem0) } < 2 {
        return 0.0;
    }
    {
        let mut st = PRESSURE_STATE.lock().unwrap();
        st.dt = dt_now;
        if st.first {
            st.new_time = true;
            st.told = t;
            st.old_mesh_tag = unsafe { mesh_tag() };
            st.dim = unsafe { usf_coord_dim() };
            let mut h_c: *const c_void = std::ptr::null();
            unsafe { usf_get_constants(&mut h_c) };
            let mut got = false;
            st.rhow = unsafe { const_real(h_c, "Water Density", &mut got) };
            if !got {
                unsafe { info("SeaPressure", "Variable Water Density not found. Setting to 1.03225e-18", 20) };
                st.rhow = 1.03225e-18;
            } else {
                let msg = unsafe { write_msg_af104("Water Density = ", st.rhow) };
                unsafe { info("SeaPressure", &msg, 20) };
            }
            let mut g2 = false;
            st.normal_flux = unsafe { get_logical(h_c, "Buoyancy Use Basal Melt", &mut g2) };
            if !g2 {
                unsafe { info("SeaPressure", "Variable Buoyancy Use Basal Melt not found. Setting to FALSE", 20) };
                st.normal_flux = false;
            } else {
                let msg = unsafe { write_msg_star_al("Buoyancy Use Basal Melt = ", st.normal_flux) };
                unsafe { info("SeaPressure", &msg, 20) };
            }
            if st.normal_flux {
                let mut g3 = false;
                let bsn = unsafe { list_get_string(h_c, "Bottom Surface Name", &mut g3) };
                if !g3 {
                    unsafe { info("SeaPressure", "Variable Bottom SUrface Name not found. Setting to DyBottom", 20) };
                    st.bottom_surface_name = pad128("DyBottom");
                } else {
                    st.bottom_surface_name = pad128(&bsn);
                    let msg = unsafe { write_msg_star_aa("Bottom Surface Name = ", st.bottom_surface_name.trim_end()) };
                    unsafe { info("SeaPressure", &msg, 20) };
                }
                let mut found_id = -1;
                let nbf = unsafe { number_of_bodyforces() };
                let mut bf = 1;
                while bf <= nbf {
                    let mut h_bfv: *const c_void = std::ptr::null();
                    unsafe { usf_bodyforce_values(bf, &mut h_bfv) };
                    let key = format!("{} Accumulation", st.bottom_surface_name.trim_end());
                    if unsafe { list_check_present(h_bfv, &key) } {
                        found_id = bf;
                    }
                    bf += 1;
                }
                st.bf_id_fs = found_id;
                if found_id < 0 {
                    unsafe { fatal("Sea Pressure", "No Basal Melt found in any body Force") };
                }
            }
            let mut h_b1: *const c_void = std::ptr::null();
            unsafe { usf_model_body_values(1, &mut h_b1) };
            let nmat = unsafe { number_of_materials() };
            let mut gf = false;
            let bf_id = unsafe { list_get_integer_mm(h_b1, "Body Force", false, &mut gf, Some((1, nmat))) };
            let mut h_bfv: *const c_void = std::ptr::null();
            unsafe { usf_bodyforce_values(bf_id, &mut h_bfv) };
            let dim = st.dim;
            if dim == 2 {
                let mut g = false;
                st.gravity = -unsafe { const_real(h_bfv, "Lateral Friction Gravity 2", &mut g) };
                if !g {
                    st.gravity = -unsafe { const_real(h_bfv, "Flow BodyForce 2", &mut g) };
                }
            } else {
                let mut g = false;
                st.gravity = -unsafe { const_real(h_bfv, "Flow BodyForce 3", &mut g) };
            }
        } else if t > st.told {
            st.new_time = true;
            st.told = t;
        }
    }
    {
        let mut st = PRESSURE_STATE.lock().unwrap();
        if unsafe { mesh_tag() } != st.old_mesh_tag {
            st.old_mesh_tag = unsafe { mesh_tag() };
            st.mesh_changed = true;
        }
        if unsafe { mesh_changed() } {
            st.mesh_changed = true;
        }
        if st.first || (st.new_time && st.mesh_changed) {
            st.mesh_changed = false;
            let nnodes = unsafe { usf_mesh_number_of_nodes() } as usize;
            st.node_on_boundary = vec![0; nnodes];
            st.number_on_boundary = 0;
            let maxnodes = unsafe { usf_max_element_nodes() } as usize;
            st.source_func = vec![0.0; maxnodes];
            st.first = false;
            let mut h_cur: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_cur) };
            let nbe = unsafe { usf_number_of_boundary_elements() };
            let mut p = 1;
            while p <= nbe {
                let mut h_be: *const c_void = std::ptr::null();
                unsafe { usf_get_boundary_element(p, &mut h_be) };
                let mut h_bc: *const c_void = std::ptr::null();
                unsafe { usf_get_bc(h_be, &mut h_bc) };
                if unsafe { get_element_family(h_be) } == 1 {
                    p += 1;
                    continue;
                }
                let mut got = false;
                let cs = unsafe { get_logical(h_bc, "Compute Sea Pressure", &mut got) };
                let compute = if got { cs } else { false };
                if compute {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    let mut nidx = vec![0 as c_int; en];
                    unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                    for &nid in &nidx {
                        let j = nid as usize;
                        if st.node_on_boundary[j - 1] == 0 {
                            st.number_on_boundary += 1;
                            st.node_on_boundary[j - 1] = st.number_on_boundary;
                        }
                    }
                }
                p += 1;
            }
            unsafe { set_current_element(h_cur) };
            let count = st.number_on_boundary as usize;
            st.s = vec![0.0; count];
            st.ns = vec![0.0; count];
            if st.normal_flux {
                st.a_perp = vec![0.0; count];
            }
        }
    }
    {
        let mut st = PRESSURE_STATE.lock().unwrap();
        if st.new_time {
            st.new_time = false;
            let mut h_be0: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_be0) };
            if h_be0.is_null() {
                unsafe { fatal("Sea Pressure", "No boundary element found") };
            }
            let h_parent = unsafe { boundary_parent(h_be0) };
            let body_id = unsafe { usf_element_body_id(h_parent) };
            let mut h_bodies: *const c_void = std::ptr::null();
            unsafe { usf_model_body_values(body_id, &mut h_bodies) };
            let mut got = false;
            let material_id = unsafe { list_get_integer(h_bodies, "Material", true, &mut got) };
            let mut h_pm: *const c_void = std::ptr::null();
            unsafe { usf_model_material_values(material_id, &mut h_pm) };
            if h_pm.is_null() {
                let msg = format!(
                    "No material values found for body no {:>10} under material id {:>10}",
                    body_id, material_id
                );
                unsafe { fatal("Sea Pressure", &msg) };
            }
            let mut g2 = false;
            st.zsl = unsafe { get_creal(h_pm, "Sea level", &mut g2) };
            if !g2 {
                unsafe { info("SeaPressure", "Variable Sea level not found. Setting to 0.0", 2) };
                st.zsl = 0.0;
            } else {
                let msg = unsafe { write_msg_af104("Sea level = ", st.zsl) };
                unsafe { info("SeaPressure", &msg, 2) };
            }
            let mut h_cur: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_cur) };
            st.s.fill(0.0);
            let dim = st.dim as usize;
            let nbe = unsafe { usf_number_of_boundary_elements() };
            let mut p = 1;
            while p <= nbe {
                let mut h_be: *const c_void = std::ptr::null();
                unsafe { usf_get_boundary_element(p, &mut h_be) };
                let mut h_bc: *const c_void = std::ptr::null();
                unsafe { usf_get_bc(h_be, &mut h_bc) };
                if unsafe { get_element_family(h_be) } == 1 {
                    p += 1;
                    continue;
                }
                let mut got2 = false;
                let cs = unsafe { get_logical(h_bc, "Compute Sea Pressure", &mut got2) };
                let compute = if got2 { cs } else { false };
                if compute {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    let mut nidx = vec![0 as c_int; en];
                    unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                    for &nid in &nidx {
                        let j = nid as usize;
                        let map = st.node_on_boundary[j - 1];
                        if dim == 2 {
                            st.s[map as usize - 1] = unsafe { model_node_y(j as c_int) };
                        } else {
                            st.s[map as usize - 1] = unsafe { model_node_z(j as c_int) };
                        }
                    }
                }
                p += 1;
            }
            unsafe { set_current_element(h_cur) };
            st.ns.fill(0.0);
            if st.normal_flux {
                st.a_perp.fill(0.0);
                let count = st.number_on_boundary as usize;
                st.normal = vec![0.0; 3 * count];
                let mut h_cur2: *const c_void = std::ptr::null();
                unsafe { usf_current_element(&mut h_cur2) };
                let mut p2 = 1;
                while p2 <= nbe {
                    let mut h_be: *const c_void = std::ptr::null();
                    unsafe { usf_get_boundary_element(p2, &mut h_be) };
                    let mut h_bc: *const c_void = std::ptr::null();
                    unsafe { usf_get_bc(h_be, &mut h_bc) };
                    if unsafe { get_element_family(h_be) } == 1 {
                        p2 += 1;
                        continue;
                    }
                    let mut got3 = false;
                    let cs = unsafe { get_logical(h_bc, "Compute Sea Spring", &mut got3) };
                    let compute = if got3 { cs } else { false };
                    // `SourceFunc = 0.0` runs whether or not ComputeS.
                    st.source_func.fill(0.0);
                    if compute {
                        let en = unsafe { element_type_nodes(h_be) } as usize;
                        let mut h_bfv: *const c_void = std::ptr::null();
                        unsafe { usf_bodyforce_values(st.bf_id_fs, &mut h_bfv) };
                        let key = format!("{} Accumulation", st.bottom_surface_name.trim_end());
                        let mut g4 = false;
                        let mut buf = vec![0.0 as c_double; en];
                        unsafe {
                            get_real_into(h_bfv, &key, &mut g4, &mut buf);
                        }
                        st.source_func[..en].copy_from_slice(&buf);
                        let mut nidx = vec![0 as c_int; en];
                        unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                        for k in 0..en {
                            let j = nidx[k] as usize;
                            let bu = unsafe { element_type_nodeu(h_be, (k + 1) as c_int) };
                            let bv = if unsafe { element_type_dim(h_be) } > 1 {
                                unsafe { element_type_nodev(h_be, (k + 1) as c_int) }
                            } else {
                                0.0
                            };
                            let mut nv = [0.0 as c_double; 3];
                            unsafe { normal_vector(h_be, bu, bv, &mut nv) };
                            let map = st.node_on_boundary[j - 1];
                            for r in 0..3 {
                                st.normal[(map as usize - 1) * 3 + r] += nv[r];
                            }
                            st.a_perp[map as usize - 1] = st.source_func[k];
                        }
                    }
                    p2 += 1;
                }
                for i in 0..count {
                    let nd = st.normal[(i * 3) + (dim - 1)];
                    if nd.abs() > 1.0e-20 {
                        let mut nsi = 1.0 + unsafe { usf_pow_dp(st.normal[i * 3] / nd, 2.0) };
                        if dim > 2 {
                            nsi += unsafe { usf_pow_dp(st.normal[i * 3 + 1] / st.normal[i * 3 + 2], 2.0) };
                        }
                        st.ns[i] = nsi.sqrt();
                    } else {
                        st.ns[i] = -999.0;
                    }
                }
                unsafe { set_current_element(h_cur2) };
            }
        }
    }
    let st = PRESSURE_STATE.lock().unwrap();
    // `j = NodeOnBoundary(nodenumber)` (1-based map; 0 = unmapped).
    let map = st.node_on_boundary[node as usize - 1];
    let (s_j, ns_j, ap_j) = (
        st.s[map as usize - 1],
        st.ns[map as usize - 1],
        if st.normal_flux { st.a_perp[map as usize - 1] } else { 0.0 },
    );
    let (gravity, rhow, zsl, dt, normal_flux) = (st.gravity, st.rhow, st.zsl, st.dt, st.normal_flux);
    drop(st);
    let mut pw = if normal_flux && ns_j > 0.0 {
        -gravity * rhow * (zsl - s_j - ap_j * dt * ns_j)
    } else {
        -gravity * rhow * (zsl - s_j)
    };
    if pw > 0.0 {
        pw = 0.0;
    }
    pw
}

// -- SeaSpring ------------------------------------------------------------

struct SpringState {
    first: bool,
    new_time: bool,
    told: c_double,
    dt: c_double,
    ns: Vec<c_double>,
    dim: c_int,
    rhow: c_double,
    gravity: c_double,
    number_on_boundary: c_int,
    node_on_boundary: Vec<c_int>,
    normal: Vec<c_double>,
    // `NormalSolution` handle (opaque bits; fetched in the mesh block,
    // consumed by the per-call `ConsistentNormalVector`).
    normal_solution: usize,
    old_mesh_tag: c_int,
    mesh_changed: bool,
}

static SPRING_STATE: Mutex<SpringState> = Mutex::new(SpringState {
    first: true,
    new_time: false,
    told: 0.0,
    dt: 0.0,
    ns: Vec::new(),
    dim: 0,
    rhow: 0.0,
    gravity: 0.0,
    number_on_boundary: 0,
    node_on_boundary: Vec::new(),
    normal: Vec::new(),
    normal_solution: 0,
    old_mesh_tag: 0,
    mesh_changed: false,
});

pub unsafe fn sea_spring(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let _y = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_time = unsafe { variable_get(h_vars, "Time", true, false) };
    let t = unsafe { var_value_at_idx(h_time, 1) };
    let mut h_c0: *const c_void = std::ptr::null();
    unsafe { usf_get_constants(&mut h_c0) };
    {
        let mut st = SPRING_STATE.lock().unwrap();
        let mut got = false;
        let aux = unsafe { const_real(h_c0, "Sea Spring Timestep Size", &mut got) };
        if got {
            st.dt = aux;
        } else {
            st.dt = unsafe { solver_dt() };
        }
    }
    let mut h_elem0: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem0) };
    if unsafe { get_element_family(h_elem0) } < 2 {
        return 1.0e20;
    }
    {
        let mut st = SPRING_STATE.lock().unwrap();
        if st.first {
            st.new_time = true;
            st.old_mesh_tag = unsafe { mesh_tag() };
        } else if t > st.told {
            st.new_time = true;
            st.told = t;
        }
        if unsafe { mesh_tag() } != st.old_mesh_tag {
            st.old_mesh_tag = unsafe { mesh_tag() };
            st.mesh_changed = true;
        }
        if unsafe { mesh_changed() } {
            st.mesh_changed = true;
        }
        if st.first || (st.new_time && st.mesh_changed) {
            st.normal_solution = unsafe { variable_get(h_vars, "Normal Vector", false, false) } as usize;
            if !st.first {
                st.node_on_boundary.clear();
                st.ns.clear();
            }
            st.first = false;
            st.new_time = false;
            let nnodes = unsafe { usf_mesh_number_of_nodes() } as usize;
            st.node_on_boundary = vec![0; nnodes];
            st.number_on_boundary = 0;
            let _n = unsafe { usf_max_element_nodes() };
            st.dim = unsafe { usf_coord_dim() };
            let dim = st.dim as usize;
            let mut got = false;
            st.rhow = unsafe { const_real(h_c0, "Water Density", &mut got) };
            if !got {
                unsafe { fatal("SeaSpring", "Variable \"Water Density\" not found in Constants.") };
            } else {
                let msg = unsafe { write_msg_af104("Water Density = ", st.rhow) };
                unsafe { info("SeaSpring", &msg, 20) };
            }
            let mut h_be0: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_be0) };
            if h_be0.is_null() {
                unsafe { fatal("SeaSpring", "No boundary element found") };
            }
            let h_parent = unsafe { boundary_parent(h_be0) };
            let mut h_bf: *const c_void = std::ptr::null();
            unsafe { usf_get_bodyforce_of(h_parent, &mut h_bf) };
            let mut h_b1: *const c_void = std::ptr::null();
            unsafe { usf_model_body_values(1, &mut h_b1) };
            let nmat = unsafe { number_of_materials() };
            let mut gf = false;
            let bf_id = unsafe { list_get_integer_mm(h_b1, "Body Force", false, &mut gf, Some((1, nmat))) };
            let mut h_bfv: *const c_void = std::ptr::null();
            unsafe { usf_bodyforce_values(bf_id, &mut h_bfv) };
            if dim == 2 {
                let mut g = false;
                st.gravity = -unsafe { const_real(h_bfv, "Lateral Friction Gravity 2", &mut g) };
                if !g {
                    st.gravity = -unsafe { const_real(h_bfv, "Flow BodyForce 2", &mut g) };
                }
            } else {
                let mut g = false;
                st.gravity = -unsafe { const_real(h_bfv, "Flow BodyForce 3", &mut g) };
            }
            let mut h_cur: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_cur) };
            let nbe = unsafe { usf_number_of_boundary_elements() };
            let mut p = 1;
            while p <= nbe {
                let mut h_be: *const c_void = std::ptr::null();
                unsafe { usf_get_boundary_element(p, &mut h_be) };
                let mut h_bc: *const c_void = std::ptr::null();
                unsafe { usf_get_bc(h_be, &mut h_bc) };
                if unsafe { get_element_family(h_be) } == 1 {
                    p += 1;
                    continue;
                }
                let mut got2 = false;
                let cs = unsafe { get_logical(h_bc, "Compute Sea Spring", &mut got2) };
                let compute = if got2 { cs } else { false };
                if compute {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    let mut nidx = vec![0 as c_int; en];
                    unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                    for &nid in &nidx {
                        let j = nid as usize;
                        if st.node_on_boundary[j - 1] == 0 {
                            st.number_on_boundary += 1;
                            st.node_on_boundary[j - 1] = st.number_on_boundary;
                        }
                    }
                }
                p += 1;
            }
            unsafe { set_current_element(h_cur) };
            let count = st.number_on_boundary as usize;
            st.ns = vec![0.0; count];
            st.normal = vec![0.0; 3 * count];
            let mut h_cur2: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_cur2) };
            let mut p2 = 1;
            while p2 <= nbe {
                let mut h_be: *const c_void = std::ptr::null();
                unsafe { usf_get_boundary_element(p2, &mut h_be) };
                let mut h_bc: *const c_void = std::ptr::null();
                unsafe { usf_get_bc(h_be, &mut h_bc) };
                if unsafe { get_element_family(h_be) } == 1 {
                    p2 += 1;
                    continue;
                }
                let mut got3 = false;
                let cs = unsafe { get_logical(h_bc, "Compute Sea Spring", &mut got3) };
                let compute = if got3 { cs } else { false };
                if compute {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    let mut nidx = vec![0 as c_int; en];
                    unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                    for (k, &nid) in nidx.iter().enumerate() {
                        let j = nid as usize;
                        let bu = unsafe { element_type_nodeu(h_be, (k + 1) as c_int) };
                        let bv = if unsafe { element_type_dim(h_be) } > 1 {
                            unsafe { element_type_nodev(h_be, (k + 1) as c_int) }
                        } else {
                            0.0
                        };
                        let mut nv = [0.0 as c_double; 3];
                        unsafe { normal_vector(h_be, bu, bv, &mut nv) };
                        let map = st.node_on_boundary[j - 1];
                        for r in 0..3 {
                            st.normal[(map as usize - 1) * 3 + r] += nv[r];
                        }
                    }
                }
                p2 += 1;
            }
            for i in 0..count {
                let nd = st.normal[i * 3 + (dim - 1)];
                if nd.abs() > 1.0e-20 {
                    let mut nsi = 1.0 + unsafe { usf_pow_dp(st.normal[i * 3] / nd, 2.0) };
                    if dim > 2 {
                        nsi += unsafe { usf_pow_dp(st.normal[i * 3 + 1] / st.normal[i * 3 + 2], 2.0) };
                    }
                    st.ns[i] = nsi.sqrt();
                } else {
                    st.ns[i] = -999.0;
                    unsafe { warn("SeaSpring", "Lower surface almost is vertically aligned") };
                }
            }
            unsafe { set_current_element(h_cur2) };
        }
    }
    let mut got = false;
    let mut nsi = -999.0;
    let h_solver = unsafe { current_solver() };
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    // `NormalSolution` fetched in the mesh block, like the source.
    let h_normal = SPRING_STATE.lock().unwrap().normal_solution as *const c_void;
    let mut node_normal = [0.0 as c_double; 3];
    let dim = SPRING_STATE.lock().unwrap().dim as usize;
    unsafe {
        consistent_normal_vector(h_solver, h_normal, h_elem, node, &mut node_normal, &mut got)
    };
    if got && node_normal[dim - 1].abs() > 1.0e-20 {
        // Integer-exponent squares: explicit multiplication.
        let mut s = 0.0;
        for &nn in node_normal.iter().take(dim - 1) {
            s += nn * nn;
        }
        // `NodeNormal(dim)**2.0_dp` is real-exponent: pow.
        nsi = (1.0 + s / unsafe { usf_pow_dp(node_normal[dim - 1], 2.0) }).sqrt();
        got = true;
    }
    if !got {
        let st = SPRING_STATE.lock().unwrap();
        let map = st.node_on_boundary[node as usize - 1];
        drop(st);
        if map > 0 {
            nsi = SPRING_STATE.lock().unwrap().ns[map as usize - 1];
        }
    }
    if nsi > 0.0 {
        let st = SPRING_STATE.lock().unwrap();
        st.gravity * st.rhow * st.dt * nsi
    } else {
        1.0e20
    }
}
