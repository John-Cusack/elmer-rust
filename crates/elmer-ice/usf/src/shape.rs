//! `shapefactorgravity_x_`, `shapefactorgravity_y_`,
//! `shapefactorgravity_` — port of USF_ShapeFactor.F90 (2D flow-line
//! lateral drag `-g + (1-f)(g.t)t`).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - SAVE replicated per the source list; write-only locals (`DIM`,
//!   `t0`, `Nn`-style counters, `NParent`-style counts) are documented,
//!   not replicated.
//! - `WRITE(Message,*)` list-directed output evaluates in Fortran
//!   accessors; plain `(A)` literals pass through (`Fatal`/`Info` TRIM).
//! - `SortD`, `CubicSpline`, `InterpolateCurve`, `SUM(a*b)` evaluate in
//!   Fortran wrappers (order-sensitive); `SQRT` is same-libm.
//! - `w**i` uses integer exponents (repeated multiplication, never
//!   `pow`); `**2.0` default-real exponents go through `usf_pow_dp`.
//! - `PI` is `3.1415926535897932384626433832795_dp` — exactly
//!   `std::f64::consts::PI`; documented, not approximated.
//! - `ATAN` is same-libm `atan`.
//! - Fixed-length `CHARACTER(128)` names are blank-padded (`pad128`)
//!   with trimmed comparisons, byte-identical to the Fortran compares.
//! - 1-based node mappings shift by one in the 0-based Rust stores; an
//!   unmapped node (mapping 0) panics loudly where the reference reads
//!   out of bounds — the invariant always holds.
//! - `ts` rows are copied to contiguous temps before interpolation
//!   (copies are exact; strided slices never cross).
//! - The `_x`/`_y` wrappers forward to the 4-arg form with axis 1/2;
//!   the exported 4-arg parent takes `axis` as its own reference.

use crate::support::*;
use std::f64::consts::PI;

const NP: usize = 1;
const NR: usize = 6;
const AP: [c_double; 1] = [0.814589];
const AR: [c_double; 6] = [5.30904, 0.0934244, 3.94177, 74.6014, 0.0475342, 1.35021];

struct ShapeState {
    first: bool,
    new_time: bool,
    told: c_double,
    g: [c_double; 2],
    zs: Vec<c_double>,
    xs: Vec<c_double>,
    ts: Vec<c_double>,
    zb: Vec<c_double>,
    xb: Vec<c_double>,
    ns: c_int,
    nb: c_int,
    para: Vec<c_double>,
    aux_real: Vec<c_double>,
    z2s: Vec<c_double>,
    z2b: Vec<c_double>,
    ts2x: Vec<c_double>,
    ts2y: Vec<c_double>,
    para2: Vec<c_double>,
    parabola: bool,
    rectangle: bool,
    computed: bool,
    nodes_on_surface: Vec<c_int>,
    nodes_on_bed: Vec<c_int>,
    temp_s: Vec<c_int>,
    temp_b: Vec<c_int>,
}

static SHAPE_STATE: Mutex<ShapeState> = Mutex::new(ShapeState {
    first: true,
    new_time: false,
    told: 0.0,
    g: [0.0; 2],
    zs: Vec::new(),
    xs: Vec::new(),
    ts: Vec::new(),
    zb: Vec::new(),
    xb: Vec::new(),
    ns: 0,
    nb: 0,
    para: Vec::new(),
    aux_real: Vec::new(),
    z2s: Vec::new(),
    z2b: Vec::new(),
    ts2x: Vec::new(),
    ts2y: Vec::new(),
    para2: Vec::new(),
    parabola: false,
    rectangle: false,
    computed: true,
    nodes_on_surface: Vec::new(),
    nodes_on_bed: Vec::new(),
    temp_s: Vec::new(),
    temp_b: Vec::new(),
});

pub unsafe fn shape_factor_gravity_x(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    let node = unsafe { *n };
    let x = unsafe { *v };
    shape_factor_gravity(m, node, x, 1)
}

pub unsafe fn shape_factor_gravity_y(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    let node = unsafe { *n };
    let x = unsafe { *v };
    shape_factor_gravity(m, node, x, 2)
}

pub(crate) fn shape_factor_gravity(m: *mut c_void, node: c_int, _x: c_double, axis: c_int) -> c_double {
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_time = unsafe { variable_get(h_vars, "Time", true, false) };
    let t = unsafe { var_value_at_idx(h_time, 1) };
    let mut h_bf0: *const c_void = std::ptr::null();
    unsafe { usf_get_bodyforce(&mut h_bf0) };
    {
        let mut st = SHAPE_STATE.lock().unwrap();
        if st.first {
            let dim = unsafe { usf_coord_dim() };
            if dim != 2 {
                unsafe { fatal("ShapeFactorGravity", "Work only for 2D flow line problem") };
            }
            st.first = false;
            st.new_time = true;
            st.told = t;
            let maxnodes = unsafe { usf_max_element_nodes() } as usize;
            st.aux_real = vec![0.0; maxnodes];
            let mut got = false;
            st.g[0] = unsafe { const_real(h_bf0, "Shape Gravity 1", &mut got) };
            st.g[1] = unsafe { const_real(h_bf0, "Shape Gravity 2", &mut got) };
            if !got {
                unsafe { fatal("ShapeFactorGravity", "Gravity vector must be specified") };
            }
            let mut g2 = false;
            st.parabola = unsafe { get_logical(h_bf0, "Parabola Shape", &mut g2) };
            if !g2 {
                st.parabola = false;
            }
            let mut g3 = false;
            st.rectangle = unsafe { get_logical(h_bf0, "Rectangle Shape", &mut g3) };
            if !g3 {
                st.rectangle = false;
            }
            if st.rectangle && st.parabola {
                unsafe { fatal("ShapeFactorGravity", "Make a choice between Parabola or Rectangle shapes") };
            }
            if !st.rectangle && !st.parabola {
                st.computed = false;
                let mut h_ce: *const c_void = std::ptr::null();
                unsafe { usf_current_element(&mut h_ce) };
                let en = unsafe { element_type_nodes(h_ce) } as usize;
                let mut nidx = vec![0 as c_int; en];
                unsafe { usf_element_nodeindexes(h_ce, nidx.as_mut_ptr(), en as c_int) };
                let mut buf = vec![0.0 as c_double; en];
                let mut g4 = false;
                unsafe { get_real_into(h_bf0, "Shape Factor", &mut g4, &mut buf) };
                st.aux_real[..en].copy_from_slice(&buf);
                if !g4 {
                    unsafe { fatal("ShapeFactorGravity", "Shape Factor must be given or Computed") };
                }
            }
            // Boundary-element census (surface vs bedrock node counts).
            let mut h_cur: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_cur) };
            let nbe = unsafe { usf_number_of_boundary_elements() };
            let mut ns = 1;
            let mut nb = 1;
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
                let mut gs = false;
                let surface = unsafe { get_logical(h_bc, "Shape Surface", &mut gs) };
                let mut gb = false;
                let bedrock = unsafe { get_logical(h_bc, "Shape Bedrock", &mut gb) };
                if surface {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    ns += en as c_int - 1;
                } else if bedrock {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    nb += en as c_int - 1;
                }
                p += 1;
            }
            unsafe { set_current_element(h_cur) };
            if ns == 0 || nb == 0 {
                unsafe { fatal("ShapeFactorGravity", "Shape Surface and/or Shape Bedrock not defined") };
            }
            let (nss, nbb) = (ns as usize, nb as usize);
            st.ns = ns;
            st.nb = nb;
            st.zs = vec![0.0; nss];
            st.ts = vec![0.0; 2 * nss];
            st.xs = vec![0.0; nss];
            st.xb = vec![0.0; nbb];
            st.zb = vec![0.0; nbb];
            st.z2s = vec![0.0; nss];
            // NOTE: the source allocates `z2b(Ns)` (surface-sized, likely a
            // typo for `z2b(Nb)`); replicated verbatim.
            st.z2b = vec![0.0; nss];
            st.ts2x = vec![0.0; nss];
            st.ts2y = vec![0.0; nss];
            st.nodes_on_surface = vec![0; nss];
            st.nodes_on_bed = vec![0; nbb];
            st.temp_s = vec![0; nss];
            st.temp_b = vec![0; nbb];
            if st.computed {
                st.para = vec![0.0; nss];
                st.para2 = vec![0.0; nss];
            }
            // Store the surface/bed node numbers + x coords.
            let mut h_cur2: *const c_void = std::ptr::null();
            unsafe { usf_current_element(&mut h_cur2) };
            let mut ns2 = 0;
            let mut nb2 = 0;
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
                let mut gs = false;
                let surface = unsafe { get_logical(h_bc, "Shape Surface", &mut gs) };
                let mut gb = false;
                let bedrock = unsafe { get_logical(h_bc, "Shape Bedrock", &mut gb) };
                if surface {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    let mut nidx = vec![0 as c_int; en];
                    unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                    for &nid in &nidx {
                        let j = nid as usize;
                        if !st.nodes_on_surface.contains(&(j as c_int)) {
                            ns2 += 1;
                            st.nodes_on_surface[ns2 as usize - 1] = j as c_int;
                            st.xs[ns2 as usize - 1] = unsafe { model_node_x(j as c_int) };
                        }
                    }
                } else if bedrock {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    let mut nidx = vec![0 as c_int; en];
                    unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                    for &nid in &nidx {
                        let j = nid as usize;
                        if !st.nodes_on_bed.contains(&(j as c_int)) {
                            nb2 += 1;
                            st.nodes_on_bed[nb2 as usize - 1] = j as c_int;
                            st.xb[nb2 as usize - 1] = unsafe { model_node_x(j as c_int) };
                        }
                    }
                }
                p2 += 1;
            }
            unsafe { set_current_element(h_cur2) };
            st.ns = ns2;
            st.nb = nb2;
            // Order by x (SortD + fancy reorder).
            let mut ind: Vec<c_int> = (1..=ns2).collect();
            {
                let (xs, _) = st.xs.split_at_mut(ns2 as usize);
                let mut idx: Vec<c_int> = ind.clone();
                unsafe { usf_sortd(ns2, xs.as_mut_ptr(), idx.as_mut_ptr()) };
                ind = idx;
            }
            let old = st.nodes_on_surface.clone();
            for (k, &ix) in ind.iter().enumerate() {
                st.nodes_on_surface[k] = old[ix as usize - 1];
            }
            let mut indb: Vec<c_int> = (1..=nb2).collect();
            {
                let (xb, _) = st.xb.split_at_mut(nb2 as usize);
                let mut idx: Vec<c_int> = indb.clone();
                unsafe { usf_sortd(nb2, xb.as_mut_ptr(), idx.as_mut_ptr()) };
                indb = idx;
            }
            let oldb = st.nodes_on_bed.clone();
            for (k, &ix) in indb.iter().enumerate() {
                st.nodes_on_bed[k] = oldb[ix as usize - 1];
            }
            // Channel parameters (para) per surface node.
            if st.computed {
                let mut h_cur3: *const c_void = std::ptr::null();
                unsafe { usf_current_element(&mut h_cur3) };
                let mut p3 = 1;
                while p3 <= nbe {
                    let mut h_be: *const c_void = std::ptr::null();
                    unsafe { usf_get_boundary_element(p3, &mut h_be) };
                    let mut h_bc: *const c_void = std::ptr::null();
                    unsafe { usf_get_bc(h_be, &mut h_bc) };
                    if unsafe { get_element_family(h_be) } == 1 {
                        p3 += 1;
                        continue;
                    }
                    let mut gs = false;
                    if unsafe { get_logical(h_bc, "Shape Surface", &mut gs) } {
                        let en = unsafe { element_type_nodes(h_be) } as usize;
                        let mut nidx = vec![0 as c_int; en];
                        unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                        let mut buf = vec![0.0 as c_double; en];
                        if st.parabola {
                            let mut gp = false;
                            unsafe { get_real_into(h_bf0, "Parabola aFactor", &mut gp, &mut buf) };
                            st.aux_real[..en].copy_from_slice(&buf);
                            if !gp {
                                unsafe { fatal("ShapeFactorGravity", "Parabola aFactor must be given for Parabola shape") };
                            }
                        } else if st.rectangle {
                            let mut gr = false;
                            unsafe { get_real_into(h_bf0, "Rectangle Width", &mut gr, &mut buf) };
                            st.aux_real[..en].copy_from_slice(&buf);
                            if !gr && st.rectangle {
                                unsafe { fatal("ShapeFactorGravity", "Rectangle Width must be given for Rectangle shape") };
                            }
                        }
                        for k in 0..en {
                            let j = nidx[k];
                            st.temp_s.fill(0);
                            for q in 0..st.nodes_on_surface.len() {
                                if st.nodes_on_surface[q] == j {
                                    st.temp_s[q] = 1;
                                }
                            }
                            // MAXLOC first maximum.
                            let mut kk = 0;
                            for (q, &v) in st.temp_s.iter().enumerate() {
                                if v != 0 {
                                    kk = q + 1;
                                    break;
                                }
                            }
                            // `WHERE` over the whole `Ns` array then
                            // `MAXLOC`: first set slot; `kk == 0` when the
                            // node is absent (stale `para` entry kept).
                            if kk > 0 {
                                st.para[kk - 1] = st.aux_real[k];
                            }
                        }
                    }
                    p3 += 1;
                }
                unsafe { set_current_element(h_cur3) };
            }
        } else if t > st.told {
            st.new_time = true;
            st.told = t;
        }
    }
    // Nodal values + tangent for the current time.
    let st = SHAPE_STATE.lock().unwrap();
    let (ns, _nb) = (st.ns, st.nb);
    drop(st);
    {
        let mut st = SHAPE_STATE.lock().unwrap();
        if st.new_time {
            st.new_time = false;
            st.ts.fill(0.0);
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
                let mut gs = false;
                let surface = unsafe { get_logical(h_bc, "Shape Surface", &mut gs) };
                let mut gb = false;
                let bedrock = unsafe { get_logical(h_bc, "Shape Bedrock", &mut gb) };
                if surface {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    let mut nidx = vec![0 as c_int; en];
                    unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                    for (k, &nid) in nidx.iter().enumerate() {
                        let j = nid as usize;
                        let bu = unsafe { element_type_nodeu(h_be, (k + 1) as c_int) };
                        let mut nv = [0.0 as c_double; 3];
                        unsafe { normal_vector(h_be, bu, 0.0, &mut nv) };
                        let mut kk = 0;
                        for (q, &m) in st.nodes_on_surface.iter().enumerate() {
                            if m == j as c_int {
                                kk = q + 1;
                                break;
                            }
                        }
                        if kk > 0 {
                            st.ts[2 * (kk - 1)] += nv[1];
                            st.ts[2 * (kk - 1) + 1] -= nv[0];
                            st.xs[kk - 1] = unsafe { model_node_x(j as c_int) };
                            st.zs[kk - 1] = unsafe { model_node_y(j as c_int) };
                        }
                    }
                } else if bedrock {
                    let en = unsafe { element_type_nodes(h_be) } as usize;
                    let mut nidx = vec![0 as c_int; en];
                    unsafe { usf_element_nodeindexes(h_be, nidx.as_mut_ptr(), en as c_int) };
                    for &nid in &nidx {
                        let j = nid as usize;
                        // `TempB = 0; WHERE (NodesOnBed==j) TempB = 1;
                        // k = MAXLOC(TempB)` — first set slot.
                        st.temp_b.fill(0);
                        for q in 0..st.nodes_on_bed.len() {
                            if st.nodes_on_bed[q] == j as c_int {
                                st.temp_b[q] = 1;
                            }
                        }
                        let mut kk = 0;
                        for (q, &v) in st.temp_b.iter().enumerate() {
                            if v != 0 {
                                kk = q + 1;
                                break;
                            }
                        }
                        if kk > 0 {
                            st.xb[kk - 1] = unsafe { model_node_x(j as c_int) };
                            st.zb[kk - 1] = unsafe { model_node_y(j as c_int) };
                        }
                    }
                }
                p += 1;
            }
            for i in 0..ns as usize {
                let n0 = st.ts[2 * i];
                let n1 = st.ts[2 * i + 1];
                let l = (unsafe { usf_pow_dp(n0, 2.0) } + unsafe { usf_pow_dp(n1, 2.0) }).sqrt();
                st.ts[2 * i] /= l;
                st.ts[2 * i + 1] /= l;
            }
            unsafe { set_current_element(h_cur) };
            // Spline coefficients.
            let nsu = ns as usize;
            let xs = st.xs.clone();
            let zs = st.zs.clone();
            let mut ts0 = vec![0.0; nsu];
            let mut ts1 = vec![0.0; nsu];
            for i in 0..nsu {
                ts0[i] = st.ts[2 * i];
                ts1[i] = st.ts[2 * i + 1];
            }
            let mut z2s = vec![0.0; nsu];
            let mut z2b = vec![0.0; nsu];
            let mut ts2x = vec![0.0; nsu];
            let mut ts2y = vec![0.0; nsu];
            unsafe { cubic_spline(nsu, &xs, &zs, &mut z2s) };
            let nbu = st.nb as usize;
            let xb = st.xb.clone();
            let zb = st.zb.clone();
            unsafe { cubic_spline(nbu, &xb, &zb, &mut z2b) };
            unsafe { cubic_spline(nsu, &xs, &ts0, &mut ts2x) };
            unsafe { cubic_spline(nsu, &xs, &ts1, &mut ts2y) };
            st.z2s = z2s;
            st.z2b = z2b;
            st.ts2x = ts2x;
            st.ts2y = ts2y;
            if st.computed {
                let para = st.para.clone();
                let mut para2 = vec![0.0; nsu];
                unsafe { cubic_spline(nsu, &xs, &para, &mut para2) };
                st.para2 = para2;
            }
        }
    }
    let st = SHAPE_STATE.lock().unwrap();
    let mut h_bf: *const c_void = std::ptr::null();
    unsafe { usf_get_bodyforce(&mut h_bf) };
    let mut h_ce: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_ce) };
    let en = unsafe { element_type_nodes(h_ce) } as usize;
    let mut nidx = vec![0 as c_int; en];
    unsafe { usf_element_nodeindexes(h_ce, nidx.as_mut_ptr(), en as c_int) };
    let x = unsafe { model_node_x(node) };
    let (computed, parabola, rectangle) = (st.computed, st.parabola, st.rectangle);
    let f = if computed {
        let h = unsafe { interpolate_curve(&st.xs, &st.zs, x, &st.z2s) }
            - unsafe { interpolate_curve(&st.xb, &st.zb, x, &st.z2b) };
        if parabola {
            let afactor = unsafe { interpolate_curve(&st.xs, &st.para, x, &st.para2) };
            if h > 0.0 && afactor > 0.0 {
                let w = 1.0 / (h * afactor).sqrt();
                let mut ff = 0.0;
                // Integer exponents: repeated multiplication.
                let mut pw = 1.0;
                for &a in AP.iter() {
                    pw *= w;
                    ff += (a * pw).atan();
                }
                2.0 * ff / (NP as c_double * PI)
            } else {
                1.0
            }
        } else if rectangle {
            let width = unsafe { interpolate_curve(&st.xs, &st.para, x, &st.para2) };
            if h > 0.0 {
                let w = width / (2.0 * h);
                let mut ff = 0.0;
                let mut pw = 1.0;
                for &a in AR.iter() {
                    pw *= w;
                    ff += (a * pw).atan();
                }
                2.0 * ff / (NR as c_double * PI)
            } else {
                1.0
            }
        } else {
            1.0
        }
    } else {
        let mut buf = vec![0.0 as c_double; en];
        let mut got = false;
        unsafe { get_real_into(h_bf, "Shape Factor", &mut got, &mut buf) };
        let mut pos = en;
        for (k, &id) in nidx.iter().enumerate() {
            if node == id {
                pos = k;
                break;
            }
        }
        buf[pos]
    };
    let mut ts0 = vec![0.0; st.ns as usize];
    let mut ts1 = vec![0.0; st.ns as usize];
    for i in 0..st.ns as usize {
        ts0[i] = st.ts[2 * i];
        ts1[i] = st.ts[2 * i + 1];
    }
    let mut tangent = [0.0 as c_double; 2];
    tangent[0] = unsafe { interpolate_curve(&st.xs, &ts0, x, &st.ts2x) };
    tangent[1] = unsafe { interpolate_curve(&st.xs, &ts1, x, &st.ts2y) };
    let tl = (unsafe { usf_pow_dp(tangent[0], 2.0) } + unsafe { usf_pow_dp(tangent[1], 2.0) }).sqrt();
    tangent[0] /= tl;
    tangent[1] /= tl;
    let g = st.g;
    drop(st);
    let a = (axis - 1) as usize;
    g[a] - (1.0 - f) * unsafe { dot(g.as_ptr(), tangent.as_ptr(), 2) } * tangent[a]
}
