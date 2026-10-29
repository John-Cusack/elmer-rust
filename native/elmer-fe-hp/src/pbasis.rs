//! Assembly of the complete p-hierarchic H1 basis of one element (nodal, then
//! edge, face and bubble functions), after
//! `ElementInfoVec_ComputePElementBasis` in `ElemInfo.F90`, and its
//! `ElementBasis` front end.

use crate::brick::*;
use crate::line::*;
use crate::maps::{bubble_dofs, edge_map, face_map, tetra_edge_direction, tetra_face_direction};
use crate::pyramid::*;
use crate::quad::*;
use crate::tetra::*;
use crate::triangle::*;
use crate::types::{EdgeDir, FaceDir, MAX_EDGES, MAX_FACES};
use crate::wedge::*;
use elmer_core::fe::{BasisEval, ElementBasis};
use elmer_core::mesh::ElementType;

/// Everything `ElementInfoVec_ComputePElementBasis` takes from the element,
/// the solver's `Def_Dofs` and the mesh edge/face definitions.
#[derive(Clone, Debug, PartialEq)]
pub struct PElementSpec {
    /// Element code: 202, 303, 404, 504, 605, 706 or 808.
    pub ecode: u16,
    /// `EdgeDegree(i)`: polynomial degree of each edge (`BDOFs+1`); `0` or
    /// `1` means no edge functions on that edge.
    pub edge_degree: [i32; MAX_EDGES],
    /// `FaceDegree(i)`: polynomial degree of each face.
    pub face_degree: [i32; MAX_FACES],
    /// Polynomial degree `p` handed to the bubble routine (after
    /// `getEffectiveBubbleP`); `0` for no bubbles.
    pub bubble_p: i32,
    /// `EdgeDirection` (1-based local nodes, from [`crate::maps::edge_direction`]
    /// or the tetra tables).
    pub edge_dir: EdgeDir,
    /// `FaceDirection`, likewise.
    pub face_dir: FaceDir,
    /// `PDefs % Serendipity`.
    pub serendipity: bool,
    /// Line element: `invertBubble` (`PDefs % isEdge` with decreasing global
    /// node numbers).
    pub invert_bubble: bool,
    /// 2D element living on a 3D face (`PDefs % isEdge`): its
    /// `FaceDirection(1:n,1)`, which directs the interior bubble.
    pub bubble_face_dir: Option<[usize; 4]>,
}

fn n_nodes(ecode: u16) -> Result<usize, String> {
    Ok(match ecode {
        202 => 2,
        303 => 3,
        404 => 4,
        504 => 4,
        605 => 5,
        706 => 6,
        808 => 8,
        _ => return Err(format!("p basis not available for element {ecode}")),
    })
}

fn n_edges(ecode: u16) -> usize {
    match ecode {
        202 => 1,
        303 => 3,
        404 => 4,
        504 => 6,
        605 => 8,
        706 => 9,
        808 => 12,
        _ => 0,
    }
}

fn n_faces(ecode: u16) -> usize {
    match ecode {
        504 => 4,
        605 | 706 => 5,
        808 => 6,
        _ => 0,
    }
}

impl PElementSpec {
    /// Uniform degree `p` on every edge, face and the bubble, with the
    /// reference edge/face maps (tetrahedra use type 1) unless explicit
    /// directions are given. Not serendipity.
    pub fn uniform(
        ecode: u16,
        p: u8,
        edge_dir: Option<&[[usize; 2]]>,
        face_dir: Option<&[[usize; 4]]>,
    ) -> Result<Self, String> {
        n_nodes(ecode)?;
        let p = i32::from(p);
        let mut edge_degree = [0; MAX_EDGES];
        edge_degree[..n_edges(ecode)].fill(p);
        let mut face_degree = [0; MAX_FACES];
        face_degree[..n_faces(ecode)].fill(p);
        let ed = match edge_dir {
            Some(d) => {
                let mut a = [[0usize; 2]; MAX_EDGES];
                a[..d.len().min(MAX_EDGES)].copy_from_slice(&d[..d.len().min(MAX_EDGES)]);
                a
            }
            None if ecode == 504 => tetra_edge_direction(1),
            None => edge_map(ecode).unwrap_or([[0; 2]; MAX_EDGES]),
        };
        let fd = match face_dir {
            Some(d) => {
                let mut a = [[0usize; 4]; MAX_FACES];
                a[..d.len().min(MAX_FACES)].copy_from_slice(&d[..d.len().min(MAX_FACES)]);
                a
            }
            None if ecode == 504 => tetra_face_direction(1),
            None => face_map(ecode).unwrap_or([[0; 4]; MAX_FACES]),
        };
        Ok(PElementSpec {
            ecode,
            edge_degree,
            face_degree,
            bubble_p: p,
            edge_dir: ed,
            face_dir: fd,
            serendipity: false,
            invert_bubble: false,
            bubble_face_dir: None,
        })
    }

    fn edge_max(&self) -> i32 {
        self.edge_degree[..n_edges(self.ecode)]
            .iter()
            .copied()
            .max()
            .unwrap_or(0)
    }

    fn face_max(&self) -> i32 {
        self.face_degree[..n_faces(self.ecode)]
            .iter()
            .copied()
            .max()
            .unwrap_or(0)
    }

    /// Upper bound of the number of functions the element can produce.
    fn capacity(&self) -> usize {
        let q = self
            .edge_max()
            .max(self.face_max())
            .max(self.bubble_p)
            .max(2) as usize;
        4 + 12 * q + 6 * q * q + q * q * q + 16
    }

    /// Number of basis functions (nodal + edge + face + bubble).
    pub fn n_basis(&self) -> Result<usize, String> {
        let mut val = vec![0.0; self.capacity()];
        let mut grad = vec![[0.0; 3]; self.capacity()];
        h1_p_basis_spec(self, [0.1, 0.1, 0.1], &mut val, &mut grad)
    }
}

fn need(len: usize, n: usize) -> Result<(), String> {
    if len < n {
        Err(format!("output holds {len} functions, {n} needed"))
    } else {
        Ok(())
    }
}

/// Values and local gradients of the complete p basis at `uvw`, in the order
/// nodal, edges, faces, bubble. Returns the number of functions written.
pub fn h1_p_basis_spec(
    s: &PElementSpec,
    uvw: [f64; 3],
    val: &mut [f64],
    grad: &mut [[f64; 3]],
) -> Result<usize, String> {
    let [u, v, w] = uvw;
    need(val.len().min(grad.len()), n_nodes(s.ecode)?)?;
    for g in grad.iter_mut() {
        *g = [0.0; 3];
    }
    let mut nbp = 0usize;
    let mut nbd = 0usize;
    let ne = n_edges(s.ecode);
    let nf = n_faces(s.ecode);
    let sd = s.serendipity;
    let ed: &[[usize; 2]] = &s.edge_dir[..ne];
    let fd: &[[usize; 4]] = &s.face_dir[..nf];
    let ep = &s.edge_degree[..ne];
    let fp = &s.face_degree[..nf];
    let edge_extra: usize = ep.iter().map(|&d| (d - 1).max(0) as usize).sum();
    let bdofs = bubble_dofs(s.ecode, s.bubble_p, sd);

    match s.ecode {
        202 => {
            line_nodal(u, val, &mut nbp);
            d_line_nodal(u, grad, &mut nbd);
            if bdofs > 0 {
                need(val.len().min(grad.len()), nbp + bdofs)?;
                line_bubble_p(u, s.bubble_p, val, &mut nbp, Some(s.invert_bubble));
                d_line_bubble_p(u, s.bubble_p, grad, &mut nbd, Some(s.invert_bubble));
            }
        }
        303 => {
            triangle_nodal_p(u, v, val, &mut nbp);
            d_triangle_nodal_p(u, v, grad, &mut nbd);
            if s.edge_max() > 1 {
                need(val.len().min(grad.len()), nbp + edge_extra)?;
                triangle_edge_p(u, v, ep, val, &mut nbp, ed);
                d_triangle_edge_p(u, v, ep, grad, &mut nbd, ed);
            }
            if bdofs > 0 {
                need(val.len().min(grad.len()), nbp + bdofs)?;
                let ln = s.bubble_face_dir.map(|d| [d[0], d[1], d[2]]);
                triangle_bubble_p(u, v, s.bubble_p, val, &mut nbp, ln);
                d_triangle_bubble_p(u, v, s.bubble_p, grad, &mut nbd, ln);
            }
        }
        404 => {
            quad_nodal(u, v, val, &mut nbp);
            d_quad_nodal(u, v, grad, &mut nbd);
            if s.edge_max() > 1 {
                need(val.len().min(grad.len()), nbp + edge_extra)?;
                if sd {
                    sd_quad_edge_p(u, v, ep, val, &mut nbp, ed);
                    sd_d_quad_edge_p(u, v, ep, grad, &mut nbd, ed);
                } else {
                    quad_edge_p(u, v, ep, val, &mut nbp, ed);
                    d_quad_edge_p(u, v, ep, grad, &mut nbd, ed);
                }
            }
            if bdofs > 0 {
                need(val.len().min(grad.len()), nbp + bdofs)?;
                let ln = s.bubble_face_dir;
                if sd {
                    sd_quad_bubble_p(u, v, s.bubble_p, val, &mut nbp, ln);
                    sd_d_quad_bubble_p(u, v, s.bubble_p, grad, &mut nbd, ln);
                } else {
                    quad_bubble_p(u, v, s.bubble_p, val, &mut nbp, ln);
                    d_quad_bubble_p(u, v, s.bubble_p, grad, &mut nbd, ln);
                }
            }
        }
        504 => {
            tetra_nodal_p(u, v, w, val, &mut nbp);
            d_tetra_nodal_p(grad, &mut nbd);
            if s.edge_max() > 1 {
                need(val.len().min(grad.len()), nbp + edge_extra)?;
                tetra_edge_p(u, v, w, ep, val, &mut nbp, ed);
                d_tetra_edge_p(u, v, w, ep, grad, &mut nbd, ed);
            }
            if s.face_max() > 1 {
                let extra: usize = fp
                    .iter()
                    .map(|&d| (0..=d).map(|j| (d - j - 2).max(0) as usize).sum::<usize>())
                    .sum();
                need(val.len().min(grad.len()), nbp + extra)?;
                tetra_face_p(u, v, w, fp, val, &mut nbp, fd);
                d_tetra_face_p(u, v, w, fp, grad, &mut nbd, fd);
            }
            if bdofs > 0 {
                need(val.len().min(grad.len()), nbp + bdofs)?;
                tetra_bubble_p(u, v, w, s.bubble_p, val, &mut nbp);
                d_tetra_bubble_p(u, v, w, s.bubble_p, grad, &mut nbd);
            }
        }
        605 => {
            if sd {
                return Err("p-Pyramid not available for serendipity scheme".into());
            }
            pyramid_nodal_p(u, v, w, val, &mut nbp);
            d_pyramid_nodal_p(u, v, w, grad, &mut nbd);
            if s.edge_max() > 1 {
                need(val.len().min(grad.len()), nbp + edge_extra)?;
                pyramid_edge_p(u, v, w, ep, val, &mut nbp, ed);
                d_pyramid_edge_p(u, v, w, ep, grad, &mut nbd, ed);
            }
            if s.face_max() > 1 {
                let mut extra = 0usize;
                extra += (0..=fp[0] - 2).map(|_| (fp[0] - 1).max(0) as usize).sum::<usize>();
                for &d in &fp[1..5] {
                    extra += (0..=d - 3).map(|j| (d - j - 2).max(0) as usize).sum::<usize>();
                }
                need(val.len().min(grad.len()), nbp + extra)?;
                pyramid_face_p(u, v, w, fp, val, &mut nbp, fd);
                d_pyramid_face_p(u, v, w, fp, grad, &mut nbd, fd);
            }
            if bdofs > 0 {
                need(val.len().min(grad.len()), nbp + bdofs)?;
                pyramid_bubble_p(u, v, w, s.bubble_p, val, &mut nbp);
                d_pyramid_bubble_p(u, v, w, s.bubble_p, grad, &mut nbd);
            }
        }
        706 => {
            wedge_nodal_p(u, v, w, val, &mut nbp);
            d_wedge_nodal_p(u, v, w, grad, &mut nbd);
            if s.edge_max() > 1 {
                need(val.len().min(grad.len()), nbp + edge_extra)?;
                if sd {
                    sd_wedge_edge_p(u, v, w, ep, val, &mut nbp, ed);
                    sd_d_wedge_edge_p(u, v, w, ep, grad, &mut nbd, ed);
                } else {
                    wedge_edge_p(u, v, w, ep, val, &mut nbp, ed);
                    d_wedge_edge_p(u, v, w, ep, grad, &mut nbd, ed);
                }
            }
            if s.face_max() > 1 {
                let mut extra = 0usize;
                for &d in &fp[..2] {
                    extra += (0..=d - 3).map(|j| (d - j - 2).max(0) as usize).sum::<usize>();
                }
                for &d in &fp[2..5] {
                    if sd {
                        extra += (2..=d - 2).map(|j| (d - j - 1).max(0) as usize).sum::<usize>();
                    } else {
                        extra += (0..=d - 2).map(|_| (d - 1).max(0) as usize).sum::<usize>();
                    }
                }
                need(val.len().min(grad.len()), nbp + extra)?;
                if sd {
                    sd_wedge_face_p(u, v, w, fp, val, &mut nbp, fd);
                    sd_d_wedge_face_p(u, v, w, fp, grad, &mut nbd, fd);
                } else {
                    wedge_face_p(u, v, w, fp, val, &mut nbp, fd);
                    d_wedge_face_p(u, v, w, fp, grad, &mut nbd, fd);
                }
            }
            if bdofs > 0 {
                need(val.len().min(grad.len()), nbp + bdofs)?;
                if sd {
                    sd_wedge_bubble_p(u, v, w, s.bubble_p, val, &mut nbp);
                    sd_d_wedge_bubble_p(u, v, w, s.bubble_p, grad, &mut nbd);
                } else {
                    wedge_bubble_p(u, v, w, s.bubble_p, val, &mut nbp);
                    d_wedge_bubble_p(u, v, w, s.bubble_p, grad, &mut nbd);
                }
            }
        }
        808 => {
            brick_nodal(u, v, w, val, &mut nbp);
            d_brick_nodal(u, v, w, grad, &mut nbd);
            if s.edge_max() > 1 {
                need(val.len().min(grad.len()), nbp + edge_extra)?;
                if sd {
                    sd_brick_edge_p(u, v, w, ep, val, &mut nbp, ed);
                    d_sd_brick_edge_p(u, v, w, ep, grad, &mut nbd, ed);
                } else {
                    brick_edge_p(u, v, w, ep, val, &mut nbp, ed);
                    d_brick_edge_p(u, v, w, ep, grad, &mut nbd, ed);
                }
            }
            if s.face_max() > 1 {
                let extra: usize = fp
                    .iter()
                    .map(|&d| (2..=d).map(|j| (d - j - 1).max(0) as usize).sum::<usize>())
                    .sum();
                need(val.len().min(grad.len()), nbp + extra)?;
                if sd {
                    sd_brick_face_p(u, v, w, fp, val, &mut nbp, fd);
                    d_sd_brick_face_p(u, v, w, fp, grad, &mut nbd, fd);
                } else {
                    brick_face_p(u, v, w, fp, val, &mut nbp, fd);
                    d_brick_face_p(u, v, w, fp, grad, &mut nbd, fd);
                }
            }
            if bdofs > 0 {
                need(val.len().min(grad.len()), nbp + bdofs)?;
                if sd {
                    sd_brick_bubble_p(u, v, w, s.bubble_p, val, &mut nbp);
                    d_sd_brick_bubble_p(u, v, w, s.bubble_p, grad, &mut nbd);
                } else {
                    brick_bubble_p(u, v, w, s.bubble_p, val, &mut nbp);
                    d_brick_bubble_p(u, v, w, s.bubble_p, grad, &mut nbd);
                }
            }
        }
        other => return Err(format!("p basis not available for element {other}")),
    }
    debug_assert_eq!(nbp, nbd);
    Ok(nbp)
}

/// Spec'd entry point: uniform degree `p` on every edge, face and bubble,
/// with the given edge/face directions (`EdgeDir`/`FaceDir` rows,
/// 1-based local nodes). Returns the number of functions written.
pub fn h1_p_basis(
    etype: ElementType,
    p: u8,
    uvw: [f64; 3],
    edge_dirs: &[[usize; 2]],
    face_dirs: &[[usize; 4]],
    out_val: &mut [f64],
    out_grad: &mut [[f64; 3]],
) -> Result<usize, String> {
    let spec = PElementSpec::uniform(etype.0, p, Some(edge_dirs), Some(face_dirs))?;
    h1_p_basis_spec(&spec, uvw, out_val, out_grad)
}

/// `(nodal, edge, face, bubble)` counts for degree `p`; see
/// [`crate::maps::p_dof_counts`].
pub fn p_dof_counts(etype: ElementType, p: u8) -> Option<(usize, usize, usize, usize)> {
    crate::maps::p_dof_counts(etype.0, p)
}

/// The p-hierarchic basis of one element as an `ElementBasis`: values and
/// global gradients with the Jacobian determinant from the element's nodal
/// geometry (`ElementMetricVec`, using only the first `n_nodes` local
/// derivatives).
#[derive(Clone, Debug)]
pub struct PElementBasis {
    pub spec: PElementSpec,
}

impl PElementBasis {
    pub fn new(spec: PElementSpec) -> Self {
        PElementBasis { spec }
    }
}

impl ElementBasis for PElementBasis {
    fn eval(
        &self,
        etype: ElementType,
        nodal_xyz: &[[f64; 3]],
        uvw: [f64; 3],
        out: &mut BasisEval,
    ) -> Result<(), String> {
        if etype.0 != self.spec.ecode {
            return Err(format!(
                "PElementBasis built for element {} evaluated as {}",
                self.spec.ecode, etype.0
            ));
        }
        let nn = n_nodes(etype.0)?;
        if nodal_xyz.len() < nn {
            return Err(format!("{} nodal coordinates given, {nn} needed", nodal_xyz.len()));
        }
        let cap = self.spec.capacity();
        let mut val = vec![0.0; cap];
        let mut grad = vec![[0.0; 3]; cap];
        let nb = h1_p_basis_spec(&self.spec, uvw, &mut val, &mut grad)?;
        let dim = match etype.family() {
            2 => 1,
            3 | 4 => 2,
            _ => 3,
        };
        // dx(j,l) = sum_k dL(k,l) * xyz(k,j), nodal functions only.
        let mut dx = [[0.0f64; 3]; 3];
        for l in 0..dim {
            for j in 0..3 {
                let mut acc = 0.0;
                for k in 0..nn {
                    acc += grad[k][l] * nodal_xyz[k][j];
                }
                dx[j][l] = acc;
            }
        }
        // Covariant metric G = dx^T dx (upper triangle, dim x dim).
        let mut g = [[0.0f64; 3]; 3];
        for j in 0..dim {
            for i in 0..=j {
                g[i][j] = dx[0][i] * dx[0][j] + dx[1][i] * dx[1][j] + dx[2][i] * dx[2][j];
                g[j][i] = g[i][j];
            }
        }
        let tiny = f64::MIN_POSITIVE;
        let (det, metric): (f64, [[f64; 3]; 3]) = match dim {
            1 => {
                let d = g[0][0];
                if d <= tiny {
                    return Err("degenerate 1D element".into());
                }
                (d, [[1.0 / d, 0.0, 0.0], [0.0; 3], [0.0; 3]])
            }
            2 => {
                let d = g[0][0] * g[1][1] - g[0][1] * g[0][1];
                if d <= tiny {
                    return Err("degenerate 2D element".into());
                }
                let s = 1.0 / d;
                (
                    d,
                    [[s * g[1][1], -s * g[0][1], 0.0], [-s * g[0][1], s * g[0][0], 0.0], [0.0; 3]],
                )
            }
            _ => {
                let d = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[1][2])
                    + g[0][1] * (g[1][2] * g[0][2] - g[0][1] * g[2][2])
                    + g[0][2] * (g[0][1] * g[1][2] - g[1][1] * g[0][2]);
                if d <= tiny {
                    return Err("degenerate 3D element".into());
                }
                let s = 1.0 / d;
                let m11 = s * (g[1][1] * g[2][2] - g[1][2] * g[1][2]);
                let m12 = -s * (g[0][1] * g[2][2] - g[0][2] * g[1][2]);
                let m22 = s * (g[0][0] * g[2][2] - g[0][2] * g[0][2]);
                let m13 = s * (g[0][1] * g[1][2] - g[1][1] * g[0][2]);
                let m23 = -s * (g[0][0] * g[1][2] - g[0][1] * g[0][2]);
                let m33 = s * (g[0][0] * g[1][1] - g[0][1] * g[0][1]);
                (d, [[m11, m12, m13], [m12, m22, m23], [m13, m23, m33]])
            }
        };
        // LtoGMap(i,l) = sum_m dx(i,m) * Metric(m,l)
        let mut l2g = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for l in 0..dim {
                let mut acc = 0.0;
                for m in 0..dim {
                    acc += dx[i][m] * metric[m][l];
                }
                l2g[i][l] = acc;
            }
        }
        out.det_j = det.sqrt();
        out.basis.clear();
        out.basis.extend_from_slice(&val[..nb]);
        out.dbasis_dx.clear();
        for gl in &grad[..nb] {
            let mut r = [0.0; 3];
            for (j, rj) in r.iter_mut().enumerate() {
                let mut acc = 0.0;
                for l in 0..dim {
                    acc += gl[l] * l2g[j][l];
                }
                *rj = acc;
            }
            out.dbasis_dx.push(r);
        }
        Ok(())
    }
}
