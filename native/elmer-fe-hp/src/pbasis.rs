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
    /// Order the functions of permuted wedge square faces as the scalar
    /// `ElementInfo` does (see [`h1_p_basis_spec`]); `false` gives the order
    /// of the vectorised `ElementInfoVec` path that `H1Basis_WedgeFaceP`
    /// implements. Only matters for full (non-serendipity) wedges.
    pub scalar_face_order: bool,
    /// Evaluate as the scalar `ElementInfo` does (PElementBase.F90 and
    /// `ElementMetric`) instead of `ElementInfoVec` (H1Basis.F90 and
    /// `ElementMetricVec`); see `scalar.rs`. Triangles and quads
    /// (both schemes) have the scalar path so far; other elements ignore it.
    pub scalar_element_info: bool,
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
            scalar_face_order: true,
            scalar_element_info: false,
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
                extra += (0..=fp[0] - 2)
                    .map(|_| (fp[0] - 1).max(0) as usize)
                    .sum::<usize>();
                for &d in &fp[1..5] {
                    extra += (0..=d - 3)
                        .map(|j| (d - j - 2).max(0) as usize)
                        .sum::<usize>();
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
                    extra += (0..=d - 3)
                        .map(|j| (d - j - 2).max(0) as usize)
                        .sum::<usize>();
                }
                for &d in &fp[2..5] {
                    if sd {
                        extra += (2..=d - 2)
                            .map(|j| (d - j - 1).max(0) as usize)
                            .sum::<usize>();
                    } else {
                        extra += (0..=d - 2).map(|_| (d - 1).max(0) as usize).sum::<usize>();
                    }
                }
                need(val.len().min(grad.len()), nbp + extra)?;
                if sd {
                    sd_wedge_face_p(u, v, w, fp, val, &mut nbp, fd);
                    sd_d_wedge_face_p(u, v, w, fp, grad, &mut nbd, fd);
                } else {
                    let start = nbp;
                    wedge_face_p(u, v, w, fp, val, &mut nbp, fd);
                    d_wedge_face_p(u, v, w, fp, grad, &mut nbd, fd);
                    if s.scalar_face_order {
                        transpose_permuted_wedge_faces(start, fp, fd, val, grad);
                    }
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

/// The scalar `ElementInfo` lists the functions of a wedge square face whose
/// first two local nodes do not lie in one triangle (`wedgeOrdering` false) with
/// the two face indices exchanged relative to `H1Basis_WedgeFaceP`; undo that
/// here by transposing the `(p-1) x (p-1)` block of such faces.
fn transpose_permuted_wedge_faces(
    start: usize,
    fp: &[i32],
    fd: &[[usize; 4]],
    val: &mut [f64],
    grad: &mut [[f64; 3]],
) {
    let mut off = start;
    for (i, &d) in fp.iter().enumerate().take(2) {
        for j in 0..=d - 3 {
            off += (d - j - 2).max(0) as usize;
        }
        let _ = i;
    }
    for i in 2..5 {
        let d = fp[i];
        let m = (d - 1).max(0) as usize;
        let nodes = fd[i];
        let nonpermuted = (nodes[0] <= 3 && nodes[1] <= 3) || (nodes[0] >= 4 && nodes[1] >= 4);
        if !nonpermuted {
            for a in 0..m {
                for b in (a + 1)..m {
                    val.swap(off + a * m + b, off + b * m + a);
                    grad.swap(off + a * m + b, off + b * m + a);
                }
            }
        }
        off += m * m;
    }
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

impl PElementBasis {
    /// The scalar `ElementInfo` of a triangle or a quad: the
    /// functions `scalar.rs` ports from PElementBase.F90, H1Basis values for
    /// the rest, geometry from `ElementMetric`.
    #[allow(clippy::too_many_arguments)]
    fn eval_scalar_2d(
        &self,
        ecode: u16,
        nodal_xyz: &[[f64; 3]],
        uvw: [f64; 3],
        mut val: Vec<f64>,
        mut grad: Vec<[f64; 3]>,
        nb: usize,
        out: &mut BasisEval,
    ) -> Result<(), String> {
        let (u, v) = (uvw[0], uvw[1]);
        let tri = ecode == 303;
        // The nodal functions: the triangle's differ from H1Basis's; the
        // quad's QuadNodalPBasisAll are H1Basis's QuadNodal already.
        let nn = if tri { 3 } else { 4 };
        if tri {
            crate::scalar::triangle_nodal(u, v, &mut val, &mut grad);
        }
        // Edge functions in H1Basis order: edge by edge, degree 2 up.
        // ElemInfo.F90:854-860 / the 404 case: invert when the edge's first
        // map node has the larger global index; edge_dir holds the map pair
        // reordered by global index.
        let maps: &[[usize; 2]] = if tri {
            &[[1, 2], [2, 3], [3, 1]]
        } else {
            &[[1, 2], [2, 3], [4, 3], [1, 4]]
        };
        let mut q = nn;
        for (i, map) in maps.iter().enumerate() {
            let invert = self.spec.edge_dir[i][0] != map[0];
            for k in 2..=self.spec.edge_degree[i].max(1) as usize {
                let f = if tri {
                    crate::scalar::triangle_edge(i + 1, k, u, v, invert)
                } else if self.spec.serendipity {
                    crate::scalar::sd_quad_edge(i + 1, k, u, v, invert)
                } else {
                    crate::scalar::quad_edge(i + 1, k, u, v, invert)
                };
                if let Some((value, g)) = f {
                    val[q] = value;
                    grad[q] = [g[0], g[1], 0.0];
                }
                q += 1;
            }
        }
        // Serendipity quad bubbles, i = 2..p-2 outer, j = 2..p-i inner, as
        // `SD_bubbles_quad`.
        if !tri
            && self.spec.serendipity
            && self.spec.bubble_p > 0
            && self.spec.bubble_face_dir.is_none()
        {
            let bp = self.spec.bubble_p as usize;
            'sd: for i in 2..=bp.saturating_sub(2) {
                for j in 2..=bp.saturating_sub(i) {
                    if q >= nb {
                        break 'sd;
                    }
                    if let Some((value, g)) = crate::scalar::sd_quad_bubble(i, j, u, v) {
                        val[q] = value;
                        grad[q] = [g[0], g[1], 0.0];
                    }
                    q += 1;
                }
            }
        }
        // Triangle bubbles, i = 0..p-3 outer, j = 0..p-i-3 inner, as
        // `bubbles_triangle` (ElemInfo.F90:901).
        if tri && self.spec.bubble_p > 0 && self.spec.bubble_face_dir.is_none() {
            let bp = self.spec.bubble_p as usize;
            'tb: for i in 0..=bp.saturating_sub(3) {
                for j in 0..=(bp.saturating_sub(3)).saturating_sub(i) {
                    if bp < 3 || q >= nb {
                        break 'tb;
                    }
                    let (value, g) = crate::scalar::triangle_bubble(i, j, u, v);
                    val[q] = value;
                    grad[q] = [g[0], g[1], 0.0];
                    q += 1;
                }
            }
        }
        // Quad interior bubbles, i (u) outer, j (v) inner, as `bubbles_quad`.
        if !tri
            && !self.spec.serendipity
            && self.spec.bubble_p > 0
            && self.spec.bubble_face_dir.is_none()
        {
            let bp = self.spec.bubble_p as usize;
            'bubbles: for i in 0..=bp.saturating_sub(2) {
                for j in 0..=bp.saturating_sub(2) {
                    if q >= nb {
                        break 'bubbles;
                    }
                    if let Some((value, g)) = crate::scalar::quad_bubble(i, j, u, v) {
                        val[q] = value;
                        grad[q] = [g[0], g[1], 0.0];
                    }
                    q += 1;
                }
            }
        }
        let cdim = if nodal_xyz[..nn].iter().all(|x| x[2] == 0.0) {
            2
        } else {
            3
        };
        // ElementMetric(q, ...) sums over all q functions with `Nodes % x(1:q)`:
        // the coordinates beyond the nodes are zero unless the caller gives
        // them (a p mesh whose `Mesh % Nodes` hold every dof, as one whose
        // geometry a solver has curved).
        let ng = nodal_xyz.len().clamp(nn, nb);
        let (det_j, l2g) = crate::scalar::element_metric(2, cdim, &nodal_xyz[..ng], &grad[..nb])?;
        out.det_j = det_j;
        out.basis.clear();
        out.basis.extend_from_slice(&val[..nb]);
        out.dbasis_dx.clear();
        for gl in &grad[..nb] {
            let mut r = [0.0; 3];
            for (j, rj) in r.iter_mut().enumerate() {
                let mut acc = 0.0;
                for l in 0..2 {
                    acc += gl[l] * l2g[j][l];
                }
                *rj = acc;
            }
            out.dbasis_dx.push(r);
        }
        Ok(())
    }
}

impl PElementBasis {
    /// Second derivatives in the reference coordinates (`ddLBasisddx`):
    /// along each reference axis a basis function's gradient is a
    /// polynomial of degree at most its p, so the derivative of its
    /// interpolant through `p + 3` Chebyshev-Lobatto points around `uvw`
    /// (half width 1/8) is the exact one up to rounding; the result is
    /// symmetrized.
    fn reference_second_derivatives(&self, uvw: [f64; 3], dim: usize, nb: usize) -> Result<Vec<[[f64; 3]; 3]>, String> {
        let s = &self.spec;
        let pmax = s
            .edge_degree
            .iter()
            .chain(&s.face_degree)
            .copied()
            .chain([s.bubble_p, 1])
            .max()
            .unwrap_or(1)
            .max(1) as usize;
        let m = pmax + 3;
        let h = 0.125;
        let t: Vec<f64> = (0..m).map(|j| h * (std::f64::consts::PI * j as f64 / (m - 1) as f64).cos()).collect();
        // L_j'(0) of the Lagrange basis on the points t.
        let w: Vec<f64> = (0..m)
            .map(|j| {
                let mut sum = 0.0;
                for i in (0..m).filter(|&i| i != j) {
                    let mut prod = 1.0 / (t[j] - t[i]);
                    for l in (0..m).filter(|&l| l != i && l != j) {
                        prod *= (0.0 - t[l]) / (t[j] - t[l]);
                    }
                    sum += prod;
                }
                sum
            })
            .collect();
        let cap = s.capacity();
        let mut val = vec![0.0; cap];
        let mut grad = vec![[0.0; 3]; cap];
        let mut dd = vec![[[0.0f64; 3]; 3]; nb];
        for k in 0..dim {
            for (j, &tj) in t.iter().enumerate() {
                let mut p = uvw;
                p[k] += tj;
                let n2 = h1_p_basis_spec(s, p, &mut val, &mut grad)?;
                if n2 != nb {
                    return Err("p basis size changed along the axis".into());
                }
                for q in 0..nb {
                    for l in 0..dim {
                        dd[q][l][k] += w[j] * grad[q][l];
                    }
                }
            }
        }
        for d in dd.iter_mut() {
            for l in 0..dim {
                for k in l + 1..dim {
                    let a = 0.5 * (d[l][k] + d[k][l]);
                    d[l][k] = a;
                    d[k][l] = a;
                }
            }
        }
        Ok(dd)
    }
}

impl ElementBasis for PElementBasis {
    /// `ddBasisddx`: the reference second derivatives
    /// ([`Self::reference_second_derivatives`]) taken to the global
    /// coordinates by `GlobalSecondDerivatives` (ElemInfo.F90:10663), with
    /// the element's geometry from all the coordinates given (as `eval`).
    fn second_derivatives(&self, etype: ElementType, nodal_xyz: &[[f64; 3]], uvw: [f64; 3]) -> Result<Vec<[[f64; 3]; 3]>, String> {
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
        let ddl = self.reference_second_derivatives(uvw, dim, nb)?;
        let ng = nodal_xyz.len().clamp(nn, nb);
        let x = &nodal_xyz[..ng];
        let mut dx = [[0.0f64; 3]; 3];
        for r in 0..3 {
            for i in 0..dim {
                dx[r][i] = (0..ng).fold(0.0, |a, q| a + x[q][r] * grad[q][i]);
            }
        }
        // Metric = the inverse of G = dx^T dx.
        let mut g = [[0.0f64; 3]; 3];
        for i in 0..dim {
            for j in 0..dim {
                g[i][j] = (0..3).fold(0.0, |a, r| a + dx[r][i] * dx[r][j]);
            }
        }
        let metric: [[f64; 3]; 3] = match dim {
            1 => [[1.0 / g[0][0], 0.0, 0.0], [0.0; 3], [0.0; 3]],
            2 => {
                let s = 1.0 / (g[0][0] * g[1][1] - g[0][1] * g[1][0]);
                [[s * g[1][1], -s * g[0][1], 0.0], [-s * g[1][0], s * g[0][0], 0.0], [0.0; 3]]
            }
            _ => {
                let d = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
                    + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
                let s = 1.0 / d;
                let c = |a: usize, b: usize, c2: usize, d2: usize| g[a][b] * g[c2][d2] - g[a][d2] * g[c2][b];
                [
                    [s * c(1, 1, 2, 2), -s * c(0, 1, 2, 2), s * c(0, 1, 1, 2)],
                    [-s * c(1, 0, 2, 2), s * c(0, 0, 2, 2), -s * c(0, 0, 1, 2)],
                    [s * c(1, 0, 2, 1), -s * c(0, 0, 2, 1), s * c(0, 0, 1, 1)],
                ]
            }
        };
        // ddx(r,i,j), C2(i,j,k) = ddx . dx, C1(i,j,k) = Metric(k,:) . C2(i,j,:).
        let mut ddx = [[[0.0f64; 3]; 3]; 3];
        for r in 0..3 {
            for i in 0..dim {
                for j in 0..dim {
                    ddx[r][i][j] = (0..ng).fold(0.0, |a, q| a + ddl[q][i][j] * x[q][r]);
                }
            }
        }
        let mut c1 = [[[0.0f64; 3]; 3]; 3];
        for i in 0..dim {
            for j in 0..dim {
                let mut c2 = [0.0f64; 3];
                for (k, c) in c2.iter_mut().enumerate().take(dim) {
                    *c = (0..3).fold(0.0, |a, r| a + ddx[r][i][j] * dx[r][k]);
                }
                for k in 0..dim {
                    c1[i][j][k] = (0..dim).fold(0.0, |a, l| a + metric[k][l] * c2[l]);
                }
            }
        }
        let mut out = Vec::with_capacity(nb);
        for q in 0..nb {
            let df = grad[q];
            let mut ddf = ddl[q];
            for i in 0..dim {
                for j in 0..dim {
                    ddf[i][j] -= (0..dim).fold(0.0, |a, k| a + c1[i][j][k] * df[k]);
                }
            }
            let mut cddf = [[0.0f64; 3]; 3];
            for i in 0..dim {
                for j in 0..dim {
                    let mut s = 0.0;
                    for k in 0..dim {
                        for l in 0..dim {
                            s += metric[i][k] * metric[j][l] * ddf[k][l];
                        }
                    }
                    cddf[i][j] = s;
                }
            }
            let mut v = [[0.0f64; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    let mut s = 0.0;
                    for k in 0..dim {
                        for l in 0..dim {
                            s += dx[i][k] * dx[j][l] * cddf[k][l];
                        }
                    }
                    v[i][j] = s;
                }
            }
            out.push(v);
        }
        Ok(out)
    }


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
            return Err(format!(
                "{} nodal coordinates given, {nn} needed",
                nodal_xyz.len()
            ));
        }
        let cap = self.spec.capacity();
        let mut val = vec![0.0; cap];
        let mut grad = vec![[0.0; 3]; cap];
        let nb = h1_p_basis_spec(&self.spec, uvw, &mut val, &mut grad)?;
        if self.spec.scalar_element_info && (etype.0 == 303 || etype.0 == 404) {
            return self.eval_scalar_2d(etype.0, nodal_xyz, uvw, val, grad, nb, out);
        }
        let dim = match etype.family() {
            2 => 1,
            3 | 4 => 2,
            _ => 3,
        };
        // dx(j,l) = sum_k dL(k,l) * xyz(k,j) over `MIN(SIZE(Nodes % x),
        // ndof)` functions (ElementMetricVec): the nodal ones, and the p
        // dofs too when their coordinates are given (a p mesh's `Mesh %
        // Nodes` past the nodes, e.g. from `Follow P Curvature`).
        let ng = nodal_xyz.len().clamp(nn, nb);
        let mut dx = [[0.0f64; 3]; 3];
        for l in 0..dim {
            for j in 0..3 {
                let mut acc = 0.0;
                for k in 0..ng {
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
                    [
                        [s * g[1][1], -s * g[0][1], 0.0],
                        [-s * g[0][1], s * g[0][0], 0.0],
                        [0.0; 3],
                    ],
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
