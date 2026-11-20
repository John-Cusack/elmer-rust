//! `AddProjectorWeakGeneric` (MortarUtils.F90:4346-5567) for nodal dofs: the
//! Galerkin level projector of a 3D mesh, whose interface meshes have been
//! mapped to a plane `(x, y)`.
//!
//! Each slave element is intersected with every master element it may
//! overlap: the polygon of edge-edge crossings and of the corners of either
//! element inside the other, sorted by angle about its centre (`SortR`),
//! with nearby points merged, is cut into a fan of triangles. Each triangle
//! is integrated with the `GaussPoints2` rule of the triangle (or `Mortar BC
//! Gauss Points`): row `nrow` of slave node `j` gets `N_j N_i w` on the
//! slave nodes and `-sgn NodeScale N_j M_k w` on the master nodes.
//!
//! With a repeating (rotational) interface the master element is shifted by
//! whole periods `ArcRange` until it meets the slave element, and an
//! antisymmetric one changes sign on odd periods. Ported: nodal dofs,
//! `Projector Node Multiplier`, `Projector Max Distance`, full circles.
//! Refused (by the caller): edge dofs, p-elements, biorthogonal bases,
//! `Axial Projector Periods`, `Projector Symmetry X/Y`.

use std::collections::BTreeMap;

use elmer_core::{BasisEval, ElemId, ElementType, Logger, MeshView};
use elmer_fe::{element_info, element_type, gauss_points, MetricStatus};
use elmer_interp::geom::{global_to_local, point_in_element, PointTol};
use elmer_interp::shape::Shape;

use crate::interface::InterfaceMesh;
use crate::MortarError;

fn fe(e: elmer_fe::FeError) -> MortarError {
    MortarError::Fatal(format!("AddProjectorWeakGeneric: {e}"))
}

/// `SortR(n, a, b)` (GeneralUtils.F90): `b` descending (the heap is built
/// for the minimum), `a` permuted with it.
/// Elmer's heapsort, including where it stores the sifted entry, so that
/// ties end up where the Fortran puts them.
pub(crate) fn sort_r(a: &mut [usize], b: &mut [f64]) {
    let n = b.len();
    if n <= 1 {
        return;
    }
    // 1-based indices as in the Fortran.
    let (mut l, mut ir) = (n / 2 + 1, n);
    loop {
        let (ra, rb);
        if l > 1 {
            l -= 1;
            ra = a[l - 1];
            rb = b[l - 1];
        } else {
            ra = a[ir - 1];
            rb = b[ir - 1];
            a[ir - 1] = a[0];
            b[ir - 1] = b[0];
            ir -= 1;
            if ir == 1 {
                a[0] = ra;
                b[0] = rb;
                return;
            }
        }
        let mut i = l;
        let mut j = l + l;
        while j <= ir {
            if j < ir && b[j - 1] > b[j] {
                j += 1;
            }
            if rb > b[j - 1] {
                a[i - 1] = a[j - 1];
                b[i - 1] = b[j - 1];
                i = j;
                j += i;
            } else {
                j = ir + 1;
            }
            a[i - 1] = ra;
            b[i - 1] = rb;
        }
    }
}

/// The geometry `LevelProjector` sets up for the integrator.
pub(crate) struct GenericSetup {
    pub node_coeff: f64,
    pub node_scale: f64,
    pub max_distance: Option<f64>,
    /// `ArcCoeff`: `2 pi Radius / 360` for a rotational interface, else 1.
    pub arc_coeff: f64,
    pub xtol: f64,
    pub ytol: f64,
    pub xrange: f64,
    pub rel_tol_y: f64,
    pub repeating: bool,
    pub anti_repeating: bool,
    pub full_circle: bool,
    pub gauss_points: Option<usize>,
    pub coord_dim: usize,
}

/// The coordinates of an element, `x` scaled by `arc`, `z` zero.
fn plane_nodes(m: &InterfaceMesh, e: usize, arc: f64) -> Vec<[f64; 3]> {
    m.element(ElemId(e as u32)).nodes.iter().map(|&n| {
        let c = m.coords[n as usize];
        [arc * c[0], c[1], 0.0]
    })
    .collect()
}

fn min_of(v: impl Iterator<Item = f64>) -> f64 {
    v.fold(f64::INFINITY, f64::min)
}
fn max_of(v: impl Iterator<Item = f64>) -> f64 {
    v.fold(f64::NEG_INFINITY, f64::max)
}

/// `Element % TYPE % NumberOfEdges` of a face element.
fn n_edges(code: u16) -> Result<usize, MortarError> {
    match code / 100 {
        3 => Ok(3),
        4 => Ok(4),
        _ => Err(MortarError::Unsupported(format!("generic level projector on boundary elements of type {code}"))),
    }
}

/// The shape `GlobalToLocal` uses: the linear element of the corners
/// (`LinCode = 101 * ne`) for a higher-order one.
fn linear_shape(code: u16, ne: usize) -> Result<&'static Shape, MortarError> {
    let lin = (101 * ne) as u16;
    let c = if code != lin { lin } else { code };
    Shape::get(c).ok_or_else(|| MortarError::Unsupported(format!("level projector on element type {code}")))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn weak_generic(
    b1: &InterfaceMesh,
    b2: &InterfaceMesh,
    s: &GenericSetup,
    node_perm: &[usize],
    inv_perm: &mut [u32],
    list: &mut [BTreeMap<u32, f64>],
    log: &dyn Logger,
) -> Result<(), MortarError> {
    const CALLER: &str = "AddProjectorWeakGeneric";
    log.info(CALLER, "Creating weak constraints using a generic integrator", 8);
    let arc_tol = s.arc_coeff * s.xtol;
    let arc_range = s.arc_coeff * s.xrange;
    let dist_tol = arc_tol * arc_tol + s.ytol * s.ytol;
    let tri = element_type(303).map_err(fe)?;
    let ngp = s.gauss_points.unwrap_or(tri.gauss_points2);
    let ip = gauss_points(ElementType(303), Some(ngp)).map_err(fe)?;
    let tol = PointTol { local: Some(1.0e-8), ..Default::default() };
    let cd = s.coord_dim;
    let mut st = MetricStatus::default();
    let (mut bt, mut bs, mut bm) = (BasisEval::default(), BasisEval::default(), BasisEval::default());
    let (mut tot_cands, mut tot_hits, mut edge_hits, mut corner_hits) = (0usize, 0usize, 0usize, 0usize);
    let (mut initial_hits, mut active_hits, mut max_sub, mut max_sub_elem) = (0usize, 0usize, 0usize, 0usize);
    let (mut tot_ref, mut tot_sum) = (0.0f64, 0.0f64);
    let (mut max_err, mut min_err) = (0.0f64, f64::MAX);
    let (mut max_ind, mut min_ind) = (0usize, 0usize);
    let (mut nslave, mut nmaster) = (0usize, 0usize);
    let z_range = |m: &InterfaceMesh, e: usize, ne: usize| -> (f64, f64) {
        let z = m.element(ElemId(e as u32)).nodes[..ne].iter().map(|&n| m.coords[n as usize][2]);
        (min_of(z.clone()), max_of(z))
    };
    for ind in 0..b1.n_bulk() {
        let el = b1.element(ElemId(ind as u32));
        let code = el.etype.0;
        let def = element_type(code).map_err(fe)?;
        let n = el.nodes.len();
        let ne = n_edges(code)?;
        let shape = linear_shape(code, ne)?;
        let mut nodes = plane_nodes(b1, ind, s.arc_coeff);
        let mut left_circle = false;
        if s.full_circle {
            left_circle = nodes[..ne].iter().all(|c| c[0].abs() > s.arc_coeff * 90.0);
            if left_circle {
                for c in nodes.iter_mut().take(n) {
                    if c[0] < 0.0 {
                        c[0] += s.arc_coeff * 360.0;
                    }
                }
            }
        }
        let xmin = min_of(nodes[..ne].iter().map(|c| c[0]));
        let xmax = max_of(nodes[..ne].iter().map(|c| c[0]));
        let ymin = min_of(nodes[..ne].iter().map(|c| c[1]));
        let ymax = max_of(nodes[..ne].iter().map(|c| c[1]));
        let (zmin, zmax) = z_range(b1, ind, ne);
        element_info(def, &nodes, [0.0; 3], cd, &mut st, &mut bs).map_err(fe)?;
        let ref_area = bs.det_j * gauss_points(el.etype, None).map_err(fe)?.s.iter().sum::<f64>();
        let mut sum_area = 0.0f64;
        // List_AddMatrixIndex: the rows of the slave nodes exist.
        for &k in el.nodes {
            let j = b1.inv_perm[k as usize];
            let nrow = node_perm[j as usize];
            if nrow != 0 {
                list[nrow - 1].entry(j).or_insert(0.0);
            }
        }
        let (mut elem_cands, mut elem_hits) = (0usize, 0usize);
        for indm in 0..b2.n_bulk() {
            let elm = b2.element(ElemId(indm as u32));
            let code_m = elm.etype.0;
            let def_m = element_type(code_m).map_err(fe)?;
            let n_m = elm.nodes.len();
            let ne_corner = (code_m / 100) as usize;
            if let Some(md) = s.max_distance {
                let (zminm, zmaxm) = z_range(b2, indm, ne_corner);
                if zmaxm < zmin - md || zminm > zmax + md {
                    continue;
                }
            }
            let mut nodes_m = plane_nodes(b2, indm, s.arc_coeff);
            let yminm = min_of(nodes_m[..ne_corner].iter().map(|c| c[1]));
            if yminm > ymax {
                continue;
            }
            let ymaxm = max_of(nodes_m[..ne_corner].iter().map(|c| c[1]));
            if ymaxm < ymin {
                continue;
            }
            if left_circle {
                if nodes_m[..ne_corner].iter().all(|c| c[0].abs() - s.arc_coeff * 90.0 < arc_tol) {
                    continue;
                }
                for c in nodes_m[..ne_corner].iter_mut() {
                    if c[0] < 0.0 {
                        c[0] += s.arc_coeff * 360.0;
                    }
                }
            }
            let (mut nrange, mut nrange1, mut nrange2) = (0i64, 0i64, 0i64);
            if s.repeating {
                let xminm = min_of(nodes_m.iter().map(|c| c[0]));
                let xmaxm = max_of(nodes_m.iter().map(|c| c[0]));
                nrange1 = ((xmaxm - xmin + arc_tol) / arc_range).floor() as i64;
                nrange2 = ((xmax - xminm + arc_tol) / arc_range).floor() as i64;
                if nrange1 != 0 {
                    for c in nodes_m.iter_mut() {
                        c[0] -= nrange1 as f64 * arc_range;
                    }
                }
                nrange = nrange1;
            }
            let mut xminm = min_of(nodes_m[..ne_corner].iter().map(|c| c[0]));
            let mut xmaxm = max_of(nodes_m[..ne_corner].iter().map(|c| c[0]));
            if s.full_circle && !left_circle && xmaxm - xminm > s.arc_coeff * 180.0 {
                continue;
            }
            let ne_m = n_edges(code_m)?;
            let shape_m = linear_shape(code_m, ne_m)?;
            // Label 200; label 100 shifts a repeating master once more.
            loop {
                if !(xminm > xmax || xmaxm < xmin) {
                    elem_cands += 1;
                    let mut x: Vec<f64> = Vec::new();
                    let mut y: Vec<f64> = Vec::new();
                    let mut corner = vec![false; ne];
                    let mut corner_m = vec![false; ne_m];
                    for i in 0..ne {
                        let i2 = if i + 1 >= ne { 0 } else { i + 1 };
                        let (x1, y1, x2, y2) = (nodes[i][0], nodes[i][1], nodes[i2][0], nodes[i2][1]);
                        for im in 0..ne_m {
                            let i2m = if im + 1 >= ne_m { 0 } else { im + 1 };
                            let (x1m, y1m, x2m, y2m) = (nodes_m[im][0], nodes_m[im][1], nodes_m[i2m][0], nodes_m[i2m][1]);
                            // A(2,2) column-major.
                            let mut a = [x2 - x1, y2 - y1, x1m - x2m, y1m - y2m];
                            let det_a = a[0] * a[3] - a[2] * a[1];
                            let abs_a = (a[0].abs() + a[2].abs()) * (a[1].abs() + a[3].abs());
                            if det_a.abs() < 1.0e-8 * abs_a + 1.0e-20 {
                                continue;
                            }
                            let b = [x1m - x1, y1m - y1];
                            elmer_fe::invert_matrix(&mut a, 2)
                                .ok_or_else(|| MortarError::Fatal(format!("{CALLER}: singular edge system")))?;
                            let c = [a[0] * b[0] + a[2] * b[1], a[1] * b[0] + a[3] * b[1]];
                            if c.iter().any(|&v| v < 0.0) || c.iter().any(|&v| v > 1.0) {
                                continue;
                            }
                            x.push(x1 + c[0] * (x2 - x1));
                            y.push(y1 + c[0] * (y2 - y1));
                            if c[0].abs() < 1.0e-6 {
                                corner[i] = true;
                            } else if (c[0] - 1.0).abs() < 1.0e-6 {
                                corner[i2] = true;
                            }
                            if c[1].abs() < 1.0e-6 {
                                corner_m[im] = true;
                            } else if (c[1] - 1.0).abs() < 1.0e-6 {
                                corner_m[i2m] = true;
                            }
                            edge_hits += 1;
                        }
                    }
                    for i in 0..ne {
                        if corner[i] {
                            continue;
                        }
                        let p = [nodes[i][0], nodes[i][1], 0.0];
                        if p[0] < xminm - arc_tol || p[0] > xmaxm + arc_tol || p[1] < yminm - s.ytol || p[1] > ymaxm + s.ytol {
                            continue;
                        }
                        if point_in_element(Shape::get(code_m).unwrap_or(shape_m), &nodes_m, p, tol, false, cd).inside {
                            x.push(p[0]);
                            y.push(p[1]);
                            corner_hits += 1;
                        }
                    }
                    for i in 0..ne_m {
                        if corner_m[i] {
                            continue;
                        }
                        let p = [nodes_m[i][0], nodes_m[i][1], 0.0];
                        if p[0] < xmin - arc_tol || p[0] > xmax + arc_tol || p[1] < ymin - s.ytol || p[1] > ymax + s.ytol {
                            continue;
                        }
                        if point_in_element(Shape::get(code).unwrap_or(shape), &nodes, p, tol, false, cd).inside {
                            x.push(p[0]);
                            y.push(p[1]);
                            corner_hits += 1;
                        }
                    }
                    let mut kmax = x.len();
                    if kmax >= 3 {
                        let sgn0: f64 = if s.anti_repeating && nrange.rem_euclid(2) != 0 { -1.0 } else { 1.0 };
                        initial_hits += kmax;
                        let xt = x.iter().sum::<f64>() / kmax as f64;
                        let yt = y.iter().sum::<f64>() / kmax as f64;
                        let mut phi: Vec<f64> = (0..kmax).map(|k| (y[k] - yt).atan2(x[k] - xt)).collect();
                        let mut inds: Vec<usize> = (0..kmax).collect();
                        sort_r(&mut inds, &mut phi);
                        let xs: Vec<f64> = inds.iter().map(|&k| x[k]).collect();
                        let ys: Vec<f64> = inds.iter().map(|&k| y[k]).collect();
                        let (mut x, mut y) = (xs, ys);
                        let mut j = 0;
                        for k in 1..kmax {
                            let dist = (x[j] - x[k]).powi(2) + (y[j] - y[k]).powi(2);
                            if dist > dist_tol {
                                j += 1;
                                if j != k {
                                    x[j] = x[k];
                                    y[j] = y[k];
                                }
                            }
                        }
                        kmax = j + 1;
                        if kmax >= 3 {
                            elem_hits += 1;
                            active_hits += kmax;
                            if kmax > max_sub {
                                max_sub = kmax;
                                max_sub_elem = ind + 1;
                            }
                            for k in 0..kmax - 2 {
                                let d_area = 0.5
                                    * ((x[k + 1] - x[0]) * (y[k + 2] - y[0]) - (x[k + 2] - x[0]) * (y[k + 1] - y[0])).abs();
                                if d_area < s.rel_tol_y * s.rel_tol_y * ref_area {
                                    continue;
                                }
                                let tnodes = [[x[0], y[0], 0.0], [x[k + 1], y[k + 1], 0.0], [x[k + 2], y[k + 2], 0.0]];
                                for g in 0..ip.s.len() {
                                    if element_info(tri, &tnodes, [ip.u[g], ip.v[g], ip.w[g]], cd, &mut st, &mut bt).is_err() {
                                        break;
                                    }
                                    let xq = bt.basis[0] * tnodes[0][0] + bt.basis[1] * tnodes[1][0] + bt.basis[2] * tnodes[2][0];
                                    let yq = bt.basis[0] * tnodes[0][1] + bt.basis[1] * tnodes[1][1] + bt.basis[2] * tnodes[2][1];
                                    let wtemp = bt.det_j * ip.s[g];
                                    sum_area += wtemp;
                                    let uvw = global_to_local(shape, &nodes[..ne], [xq, yq, 0.0], cd);
                                    element_info(def, &nodes, uvw, cd, &mut st, &mut bs).map_err(fe)?;
                                    let uvw_m = global_to_local(shape_m, &nodes_m[..ne_m], [xq, yq, 0.0], cd);
                                    if element_info(def_m, &nodes_m, uvw_m, cd, &mut st, &mut bm).is_err() {
                                        continue;
                                    }
                                    for jn in 0..n {
                                        let jj = b1.inv_perm[el.nodes[jn] as usize];
                                        let nrow = node_perm[jj as usize];
                                        if nrow == 0 {
                                            continue;
                                        }
                                        inv_perm[nrow - 1] = jj;
                                        let val = bs.basis[jn] * wtemp;
                                        let row = &mut list[nrow - 1];
                                        for i in 0..n {
                                            nslave += 1;
                                            *row.entry(b1.inv_perm[el.nodes[i] as usize]).or_insert(0.0) +=
                                                s.node_coeff * bs.basis[i] * val;
                                        }
                                        for i in 0..n_m {
                                            nmaster += 1;
                                            *row.entry(b2.inv_perm[elm.nodes[i] as usize]).or_insert(0.0) +=
                                                -sgn0 * s.node_scale * s.node_coeff * bm.basis[i] * val;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                // Label 100.
                if s.repeating && nrange != nrange2 {
                    nrange = nrange2;
                    for c in nodes_m.iter_mut() {
                        c[0] += arc_range * (nrange2 - nrange1) as f64;
                    }
                    xminm = min_of(nodes_m[..ne_m].iter().map(|c| c[0]));
                    xmaxm = max_of(nodes_m[..ne_m].iter().map(|c| c[0]));
                    continue;
                }
                break;
            }
        }
        tot_cands += elem_cands;
        tot_hits += elem_hits;
        tot_sum += sum_area;
        tot_ref += ref_area;
        let err = sum_area / ref_area;
        if err > max_err {
            max_err = err;
            max_ind = ind + 1;
        }
        if err < min_err {
            min_err = err;
            min_ind = ind + 1;
        }
    }
    let es = elmer_solver::fortran_fmt::es;
    log.info(CALLER, &format!("Number of integration pair candidates: {tot_cands}"), 10);
    log.info(CALLER, &format!("Number of integration pairs: {tot_hits}"), 10);
    log.info(CALLER, &format!("Number of edge intersections: {edge_hits}"), 10);
    log.info(CALLER, &format!("Number of corners inside element: {corner_hits}"), 10);
    log.info(CALLER, &format!("Number of initial corners: {initial_hits}"), 10);
    log.info(CALLER, &format!("Number of active corners: {active_hits}"), 10);
    log.info(CALLER, &format!("Number of most subelement corners: {max_sub}"), 10);
    log.info(CALLER, &format!("Element of most subelement corners: {max_sub_elem}"), 10);
    log.info(CALLER, &format!("Total reference area:{}", es(tot_ref, 12, 5)), 8);
    log.info(CALLER, &format!("Total integrated area:{}", es(tot_sum, 12, 5)), 8);
    log.info(CALLER, &format!("Average ratio in area integration:{}", es(tot_sum / tot_ref, 15, 6)), 8);
    log.info(
        CALLER,
        &format!("Maximum relative discrepancy in areas (element: {max_ind}):{}", es(max_err - 1.0, 12, 4)),
        8,
    );
    log.info(
        CALLER,
        &format!("Minimum relative discrepancy in areas (element: {min_ind}):{}", es(min_err - 1.0, 12, 4)),
        8,
    );
    log.info(CALLER, &format!("Number of slave entries: {nslave}"), 10);
    log.info(CALLER, &format!("Number of master entries: {nmaster}"), 10);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sort_r_sorts_descending_and_permutes() {
        let mut b = vec![0.3, -1.0, 2.0, 0.1, -0.5];
        let mut a: Vec<usize> = (0..b.len()).collect();
        sort_r(&mut a, &mut b);
        // As the Fortran sorts: descending.
        assert_eq!(b, vec![2.0, 0.3, 0.1, -0.5, -1.0]);
        assert_eq!(a, vec![2, 0, 3, 4, 1]);
    }
}
