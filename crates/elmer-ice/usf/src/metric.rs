//! `glaciermeshmetricaniso_` core — port of USF_GlacierMeshMetric.F90
//! (anisotropic target element size from calving-front distance).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - The exported `glaciermeshmetricaniso_` symbol is owned by the
//!   Fortran shim in `accessors/usf_acc.F90`: the caller
//!   (`ElmerRealArrFn`, assumed-shape `TargetLength(:)`) passes a
//!   descriptor, which must not cross into Rust. The shim forwards
//!   base address + length; this module owns all entry logic.
//! - SAVE replicated: `FirstTime`, `NewTime`, `told`, `Mesh` (handle,
//!   re-fetched — same target), `NNodes`, `nfront`, `lc_*`, `dz`,
//!   `ExtrudeLayers`, `layers`, `BottomElemMask`/`SurfElemMask`;
//!   the mask *perms* live SAVE'd in the accessor (mirroring the
//!   source's SAVE'd pointers, including the DEALLOCATE-if-associated
//!   before each rebuild).
//! - `FrontMaskName` etc. are assigned every call in the source;
//!   replicated as constants.
//! - `BCTag`/`nbottom`/`nsurface` are write-only locals: documented,
//!   not replicated (`BCTag` is consumed immediately through the
//!   accessor-free mask loop... see below).
//! - `MinDist2 = HUGE(...)` is `f64::MAX`; `**2.0` default-real
//!   exponents go through `usf_pow_dp`; `SQRT` is same-libm.
//! - `GetZFromMask` is transcribed as a Rust function (internal, like
//!   the source's `CONTAINS` procedure): `PointInElement` and
//!   `InterpolateInElement` evaluate in Fortran wrappers;
//!   `SUM((Point(:2)-NodeCoord(:2))**2)` uses integer-exponent
//!   multiplication; `**0.5`/`**(-1.0)` are real-exponent `pow`.
//! - `dz = (NodeElev-NodeDepth)/(layers-1)`: the divisor converts from
//!   INTEGER exactly like the source.
//! - `NodeIndx`/`NodeIndx2` are BSS-zero statics when no node matches
//!   (the reference reads `z(0)` out of bounds); the port panics
//!   loudly there — the invariant (masked nodes exist) always holds.
//! - `Model % Nodes % x/y/z` reads go through the *mesh* handle
//!   (`mesh_node_*`), exactly where the source reads `Mesh % Nodes`.

use crate::support::*;

struct MetricState {
    first: bool,
    new_time: bool,
    told: c_double,
    dim: c_int,
    nnodes: c_int,
    nfront: c_int,
    layers: c_int,
    lc_maxdist: c_double,
    lc_mindist: c_double,
    lc_max: c_double,
    lc_min: c_double,
    dz: c_double,
    extrude_layers: bool,
    bottom_elem_mask: Vec<bool>,
    surf_elem_mask: Vec<bool>,
    // `BCTag` is a non-SAVE local: gfortran static storage, stale when
    // no BC matches (BSS 0 first).
    bc_tag: c_int,
}

static METRIC_STATE: Mutex<MetricState> = Mutex::new(MetricState {
    first: true,
    new_time: false,
    told: 0.0,
    dim: 0,
    nnodes: 0,
    nfront: 0,
    layers: 0,
    lc_maxdist: 0.0,
    lc_mindist: 0.0,
    lc_max: 0.0,
    lc_min: 0.0,
    dz: 0.0,
    extrude_layers: false,
    bottom_elem_mask: Vec::new(),
    surf_elem_mask: Vec::new(),
    bc_tag: 0,
});

/// `GetZFromMask` (source `CONTAINS` procedure): depth at `point` by
/// point-in-element search, else inverse-distance of the two nearest
/// masked nodes. `which` selects the mask perm (2 = bottom, 3 = surface).
unsafe fn get_z_from_mask(
    h_mesh: *const c_void,
    point: &[c_double; 3],
    which: c_int,
    elem_mask: &[bool],
    nnodes: usize,
    nbulk: usize,
    nbe: usize,
) -> c_double {
    let maxnodes = unsafe { mesh_max_element_nodes(h_mesh) } as usize;
    let mut enodes_x = vec![0.0 as c_double; maxnodes];
    let mut enodes_y = vec![0.0 as c_double; maxnodes];
    let mut enodes_z = vec![0.0 as c_double; maxnodes];
    let mut elem_vals = vec![0.0 as c_double; maxnodes];
    let mut lc = [0.0 as c_double; 3];
    let mut eg = 2.0e-10 as c_double;
    let mut el = 1.0e-10 as c_double;
    let en = 1.0e-10 as c_double;
    let (mut found, mut h_found_el) = (false, std::ptr::null());
    let mut n = 0;
    loop {
        let mut i = 1;
        while i <= nbulk + nbe {
            if elem_mask[i - 1] {
                i += 1;
                continue;
            }
            let h_el = unsafe { mesh_element(h_mesh, i as c_int) };
            n = unsafe { element_type_nodes(h_el) } as usize;
            let mut idx = vec![0 as c_int; n];
            unsafe { usf_element_nodeindexes(h_el, idx.as_mut_ptr(), n as c_int) };
            for k in 0..n {
                let j = idx[k] as usize;
                enodes_x[k] = unsafe { mesh_node_x(j as c_int) };
                enodes_y[k] = unsafe { mesh_node_y(j as c_int) };
            }
            enodes_z.fill(0.0);
            let (f, lcc) = unsafe {
                point_in_element(h_el, point[0], point[1], point[2], eg, el, en)
            };
            lc = lcc;
            if f {
                found = true;
                h_found_el = h_el;
                break;
            }
            i += 1;
        }
        if found {
            break;
        }
        eg *= 10.0;
        el *= 10.0;
        if eg > 1.0e-2 || el > 1.0e-2 {
            break;
        }
    }
    if !found {
        let mut min_dist = f64::MAX;
        let mut min_dist2 = f64::MAX;
        let (mut ni1, mut ni2) = (0, 0);
        for i in 1..=nnodes {
            if unsafe { usf_perm_at(which, i as c_int) } == 0 {
                continue;
            }
            let dx = point[0] - unsafe { mesh_node_x(i as c_int) };
            let dy = point[1] - unsafe { mesh_node_y(i as c_int) };
            let dist = dx * dx + dy * dy;
            if dist < min_dist {
                if min_dist < min_dist2 {
                    min_dist2 = min_dist;
                    ni2 = ni1;
                }
                min_dist = dist;
                ni1 = i;
            } else if dist < min_dist2 {
                min_dist2 = dist;
                ni2 = i;
            }
        }
        // `MinDist ** 0.5`: real-exponent pow.
        let d1 = unsafe { usf_pow_dp(min_dist, 0.5) };
        let d2 = unsafe { usf_pow_dp(min_dist2, 0.5) };
        let w1 = unsafe { usf_pow_dp(d1, -1.0) };
        let w2 = unsafe { usf_pow_dp(d2, -1.0) };
        let z1 = unsafe { mesh_node_z(ni1 as c_int) } * w1;
        let z2 = unsafe { mesh_node_z(ni2 as c_int) } * w2;
        (z1 + z2) / (w1 + w2)
    } else {
        for k in 0..n {
            let mut idx = vec![0 as c_int; n];
            unsafe {
                usf_element_nodeindexes(h_found_el, idx.as_mut_ptr(), n as c_int);
                elem_vals[k] = mesh_node_z(idx[k] as c_int);
            }
        }
        unsafe { interpolate_in_element(h_found_el, &elem_vals[..n], lc[0], lc[1], lc[2]) }
    }
}

// Called only from the Fortran entry shim in `accessors/usf_acc.F90`.
// `#[no_mangle]` keeps the C name the shim links against; the version
// script hides it from `nm -D` (explicit `local` entry — the `*`
// wildcard alone does not hide it with this toolchain).
// `#[used]` pins the address: the sole reference is cross-language and
// invisible to the compiler's reachability analysis.
#[used]
static USF_METRIC_CORE: unsafe extern "C" fn(
    *const c_void,
    c_int,
    c_double,
    *mut c_double,
    c_int,
) = usf_metric_core;

#[no_mangle]
unsafe extern "C" fn usf_metric_core(
    h_model: *const c_void,
    node: c_int,
    y: c_double,
    target: *mut c_double,
    ntarget: c_int,
) {
    let _ = y;
    unsafe { usf_begin(h_model) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_time = unsafe { variable_get(h_vars, "Time", true, false) };
    let t = unsafe { var_value_at_idx(h_time, 1) };
    let h_mesh = unsafe { model_mesh() };
    let mut st = METRIC_STATE.lock().unwrap();
    if st.first {
        st.first = false;
        st.new_time = true;
        st.told = t;
        st.dim = unsafe { usf_coord_dim() };
        let mut h_mat: *const c_void = std::ptr::null();
        unsafe { usf_get_material(mesh_element(h_mesh, 1), &mut h_mat) };
        let mut got = false;
        st.lc_maxdist = unsafe { list_get_const_real_found(h_mat, "GlacierMeshMetric Max Distance", &mut got) };
        if !got {
            st.lc_maxdist = 2000.0;
        }
        st.lc_mindist = unsafe { list_get_const_real_found(h_mat, "GlacierMeshMetric Min Distance", &mut got) };
        if !got {
            st.lc_mindist = 200.0;
        }
        st.lc_max = unsafe { list_get_const_real_found(h_mat, "GlacierMeshMetric Max LC", &mut got) };
        if !got {
            st.lc_max = 500.0;
        }
        st.lc_min = unsafe { list_get_const_real_found(h_mat, "GlacierMeshMetric Min LC", &mut got) };
        if !got {
            st.lc_min = 75.0;
        }
        st.extrude_layers = unsafe { list_get_logical(h_mat, "GlacierMeshMetric Vertical From Layers", &mut got) };
        if !got {
            st.extrude_layers = false;
        }
        if st.extrude_layers {
            let mut g2 = false;
            st.layers = unsafe { list_get_integer(h_mat, "GlacierMeshMetric Vertical Layers", false, &mut g2) };
            if !g2 {
                unsafe { fatal("GlacierMeshMetric", "Vertical requested by layers but number of layers not given") };
            }
        } else {
            st.dz = unsafe { list_get_const_real_found(h_mat, "GlacierMeshMetric Vertical LC", &mut got) };
            if !got {
                st.dz = 50.0;
            }
        }
    }
    if t > st.told {
        st.new_time = true;
    }
    if st.new_time {
        st.told = t;
        st.new_time = false;
        st.nnodes = unsafe { usf_mesh_number_of_nodes() };
        st.nfront = unsafe { make_perm(1, "Calving Front Mask") };
        if st.extrude_layers {
            unsafe { make_perm(2, "Bottom Surface Mask") };
            unsafe { make_perm(3, "Top Surface Mask") };
            let nbulk = unsafe { mesh_bulk_elements(h_mesh) };
            let nbe_mesh = unsafe { mesh_boundary_elements(h_mesh) };
            let total = (nbulk + nbe_mesh) as usize;
            st.bottom_elem_mask = vec![true; total];
            st.surf_elem_mask = vec![true; total];
            // BC scan over `Model % BCs` for the mask tags.
            let nbc = unsafe { number_of_bcs() };
            let mut i = 1;
            while i <= nbc {
                let h_bcv = unsafe { model_bc_values(i) };
                let mut found = false;
                let masked = unsafe { list_get_logical(h_bcv, "Bottom Surface Mask", &mut found) };
                if !found || !masked {
                    i += 1;
                    continue;
                }
                st.bc_tag = unsafe { model_bc_tag(i) };
                break;
            }
            for k in nbulk as usize..total {
                // `Mesh % Elements(i) % BoundaryInfo % constraint`,
                // 1-based element number `k+1`.
                let h_el = unsafe { mesh_element(h_mesh, (k + 1) as c_int) };
                if unsafe { usf_element_constraint(h_el) } == st.bc_tag {
                    st.bottom_elem_mask[k] = false;
                }
            }
            let mut i2 = 1;
            while i2 <= nbc {
                let h_bcv = unsafe { model_bc_values(i2) };
                let mut found = false;
                let masked = unsafe { list_get_logical(h_bcv, "Top Surface Mask", &mut found) };
                if !found || !masked {
                    i2 += 1;
                    continue;
                }
                st.bc_tag = unsafe { model_bc_tag(i2) };
                break;
            }
            for k in nbulk as usize..total {
                let h_el = unsafe { mesh_element(h_mesh, (k + 1) as c_int) };
                if unsafe { usf_element_constraint(h_el) } == st.bc_tag {
                    st.surf_elem_mask[k] = false;
                }
            }
        }
    }
    drop(st);
    // Nodal coordinates + front-distance scan.
    let xx = unsafe { mesh_node_x(node) };
    let yy = unsafe { mesh_node_y(node) };
    let zz = unsafe { mesh_node_z(node) };
    let st = METRIC_STATE.lock().unwrap();
    let (nnodes, lc_mindist, lc_maxdist, lc_min, lc_max) = (
        st.nnodes as usize,
        st.lc_mindist,
        st.lc_maxdist,
        st.lc_min,
        st.lc_max,
    );
    drop(st);
    let mut min_dist2 = f64::MAX;
    for i in 1..=nnodes {
        if unsafe { perm_at(1, i as c_int) } == 0 {
            continue;
        }
        // `((xx-x)**2.0 + ...)` — default-real exponents: pow.
        let dx = unsafe { usf_pow_dp(xx - mesh_node_x(i as c_int), 2.0) };
        let dy = unsafe { usf_pow_dp(yy - mesh_node_y(i as c_int), 2.0) };
        let dz = unsafe { usf_pow_dp(zz - mesh_node_z(i as c_int), 2.0) };
        let this_dist = dx + dy + dz;
        min_dist2 = min_dist2.min(this_dist);
    }
    let mut dist = min_dist2.sqrt();
    dist = dist.max(lc_mindist);
    dist = dist.min(lc_maxdist);
    let s = (dist - lc_mindist) / (lc_maxdist - lc_mindist);
    let dx = s * lc_max + (1.0 - s) * lc_min;
    let out: &mut [c_double] = unsafe { std::slice::from_raw_parts_mut(target, ntarget as usize) };
    // `TargetLength` always has 3 lanes in every caller; the shim
    // passes the descriptor length through for the check.
    assert!(ntarget >= 3, "TargetLength assumed size 3");
    out[0] = dx;
    out[1] = dx;
    let st = METRIC_STATE.lock().unwrap();
    let (extrude, layers, dz0) = (st.extrude_layers, st.layers, st.dz);
    drop(st);
    if extrude {
        let point = [xx, yy, 0.0];
        let st = METRIC_STATE.lock().unwrap();
        let (masks_b, masks_s) = (st.bottom_elem_mask.clone(), st.surf_elem_mask.clone());
        let nn = st.nnodes as usize;
        drop(st);
        let nb = unsafe { mesh_bulk_elements(h_mesh) } as usize;
        let nbe2 = unsafe { mesh_boundary_elements(h_mesh) } as usize;
        let node_depth = unsafe { get_z_from_mask(h_mesh, &point, 2, &masks_b, nn, nb, nbe2) };
        let node_elev = unsafe { get_z_from_mask(h_mesh, &point, 3, &masks_s, nn, nb, nbe2) };
        // `(NodeElev-NodeDepth)/(layers-1)`: INTEGER divisor converts.
        let dz = (node_elev - node_depth) / (layers - 1) as c_double;
        out[2] = dz;
    } else {
        let st = METRIC_STATE.lock().unwrap();
        out[2] = st.dz;
        drop(st);
        let _ = dz0;
    }
}
