//! The 2-D face (edge) table of `FindMeshEdges2D` (MeshBasics.F90) as the
//! particle march uses it, and the robust segment tests of ParticleUtils.F90
//! (`SegmentElementIntersection2`, `SegmentElementInside`) on top of
//! `LineFaceIntersection2` (ElementGeometry.F90).

use std::collections::HashMap;

use elmer_core::{ElemId, MeshView};

pub(crate) type V3 = [f64; 3];

/// One edge: its two nodes in the order of the first element that met it,
/// and its bulk parents (1-based element indexes, 0 for none).
#[derive(Clone, Debug)]
pub(crate) struct Face {
    pub nodes: [u32; 2],
    pub left: usize,
    pub right: usize,
}

pub(crate) struct Faces {
    pub faces: Vec<Face>,
    /// `Element % EdgeIndexes` of each bulk element (0-based face ids).
    pub of_element: Vec<Vec<usize>>,
}

impl Faces {
    /// `FindMeshEdges2D`: element `k`'s edge `i` joins its nodes `i` and
    /// `i + 1` (cyclically); an edge keeps the node order of the first
    /// element with it, its `Left` is the first bulk element and its `Right`
    /// the second. Bulk elements must be linear triangles or quadrilaterals.
    pub fn build(mesh: &dyn MeshView) -> Result<Faces, String> {
        let nb = mesh.n_bulk();
        let mut faces: Vec<Face> = Vec::new();
        let mut hash: HashMap<(u32, u32), usize> = HashMap::new();
        let mut of_element = Vec::with_capacity(nb);
        for t in 0..nb {
            let el = mesh.element(ElemId(t as u32));
            let n = match el.etype.0 {
                303 => 3,
                404 => 4,
                c => return Err(format!("bulk element type {c} (303 and 404 only)")),
            };
            let mut idx = Vec::with_capacity(n);
            for k in 0..n {
                let a = el.nodes[k];
                let b = el.nodes[(k + 1) % n];
                let key = (a.min(b), a.max(b));
                let f = *hash.entry(key).or_insert_with(|| {
                    faces.push(Face {
                        nodes: [a, b],
                        left: 0,
                        right: 0,
                    });
                    faces.len() - 1
                });
                if faces[f].left != 0 {
                    faces[f].right = t + 1;
                } else {
                    faces[f].left = t + 1;
                }
                idx.push(f);
            }
            of_element.push(idx);
        }
        Ok(Faces { faces, of_element })
    }
}

/// The node coordinates of bulk element `e1` (1-based).
pub(crate) fn element_nodes(mesh: &dyn MeshView, e1: usize) -> Vec<V3> {
    mesh.element(ElemId((e1 - 1) as u32)).nodes.iter().map(|&k| mesh.coord(k)).collect()
}

/// `SUM( x(1:n) ) / n` per coordinate.
pub(crate) fn center(x: &[V3]) -> V3 {
    let n = x.len() as f64;
    std::array::from_fn(|d| {
        let mut s = 0.0;
        for p in x {
            s += p[d];
        }
        s / n
    })
}

/// `LineFaceIntersection2` for a two-node face: the segment parameter of
/// the crossing of `rinit -> rfin` with the face's line, `None` when they
/// are parallel or the crossing is off the face.
pub(crate) fn line_face_intersection2(face: [V3; 2], rinit: V3, rfin: V3) -> Option<f64> {
    let eps = f64::EPSILON;
    let eps2 = f64::MIN_POSITIVE.sqrt();
    let eps3 = 1.0e-12;
    // A(row, col), column-major for InvertMatrix.
    let a11 = rfin[0] - rinit[0];
    let a21 = rfin[1] - rinit[1];
    let a12 = face[0][0] - face[1][0];
    let a22 = face[0][1] - face[1][1];
    let det = a11 * a22 - a12 * a21;
    let mut r1 = 0.0;
    for v in [a11, a12] {
        r1 += v.abs();
    }
    let mut r2 = 0.0;
    for v in [a21, a22] {
        r2 += v.abs();
    }
    let abs_a = r1 * r2;
    if det.abs() <= eps * abs_a + eps2 {
        return None;
    }
    let b = [face[0][0] - rinit[0], face[0][1] - rinit[1]];
    let mut m = [a11, a21, a12, a22];
    elmer_fe::invert_matrix(&mut m, 2)?;
    let c1 = 0.0 + m[0] * b[0] + m[2] * b[1];
    let c2 = 0.0 + m[1] * b[0] + m[3] * b[1];
    if c2 < -eps3 || c2 > 1.0 + eps3 {
        return None;
    }
    Some(c1)
}

fn face_coords(mesh: &dyn MeshView, f: &Face) -> [V3; 2] {
    [mesh.coord(f.nodes[0]), mesh.coord(f.nodes[1])]
}

/// `SegmentElementIntersection2`: the first face of element `e1` crossed
/// at a parameter above `1e-10` (else the largest), and the parameter;
/// `None` when the largest is below `-1e-7`.
pub(crate) fn segment_element_intersection2(
    mesh: &dyn MeshView,
    faces: &Faces,
    e1: usize,
    rinit: V3,
    rfin: V3,
) -> (Option<usize>, f64) {
    let pos_eps = 1.0e-10;
    let neg_eps = -1.0e-7;
    let mut min_lambda = -f64::MAX;
    let mut face = None;
    for &f in &faces.of_element[e1 - 1] {
        let Some(lambda) = line_face_intersection2(face_coords(mesh, &faces.faces[f]), rinit, rfin) else {
            continue;
        };
        if lambda > min_lambda {
            min_lambda = lambda;
            face = Some(f);
            if min_lambda > pos_eps {
                break;
            }
        }
    }
    if min_lambda < neg_eps {
        (None, min_lambda)
    } else {
        (face, min_lambda.max(0.0))
    }
}

/// `SegmentElementInside`: `rfin` is inside element `e1` when it passes the
/// bounding box and the segment from the centre reaches no face before it.
pub(crate) fn segment_element_inside(mesh: &dyn MeshView, faces: &Faces, e1: usize, rfin: V3) -> bool {
    let eps = f64::EPSILON;
    let x = element_nodes(mesh, e1);
    let edim = 2;
    for (i, &r) in rfin.iter().enumerate().take(edim) {
        let mut maxx = x[0][i];
        let mut minx = x[0][i];
        for p in &x[1..] {
            maxx = maxx.max(p[i]);
            minx = minx.min(p[i]);
        }
        let dx = eps * (maxx - minx);
        if r < minx - dx || r > maxx + dx {
            return false;
        }
    }
    let rinit = center(&x);
    let mut hits = 0;
    let mut min_lambda = f64::MAX;
    for &f in &faces.of_element[e1 - 1] {
        let Some(lambda) = line_face_intersection2(face_coords(mesh, &faces.faces[f]), rinit, rfin) else {
            continue;
        };
        hits += 1;
        if lambda > 0.0 && min_lambda > lambda {
            min_lambda = lambda;
        }
    }
    hits > 0 && min_lambda > 1.0 - eps && min_lambda > -eps
}
