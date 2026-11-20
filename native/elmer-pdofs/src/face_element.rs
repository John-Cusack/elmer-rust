//! The face (div-conforming) element pieces that `ModelMixedPoisson` and
//! the `FaceFEInterpolationTest` fixture share: the bulk element with its
//! face maps and Gauss rule ([`FaceCell`]), the slots of a face's dofs, and
//! the `<Variable> {f}` Dirichlet values of `DefaultDirichletBCs`
//! (DefUtils.F90, the `{f}` branch): `VectorElementEdgeDOFs(...,
//! FaceElement = .TRUE.)` on the edges of triangles (RT(0), BDM(1)) and
//! `FaceElementDOFs` on the triangular faces of tetrahedra (the Nedelec
//! tetrahedra of the first and second kind). Of the 48-dof brick's
//! (`LUSolve` of the face mass matrix) only zero values are ported.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{ElemId, ElementBasis, BasisEval, ElementType, ParamMap, NO_DOF};
use elmer_fe::quadrature::GaussRule;
use elmer_fe_hp::face::{
    brick_face_orientation, face_element_info_brick, face_element_info_tetra, face_element_info_triangle,
    tetra_face_basis_ordering, tetra_face_orientation, triangle_face_orientation, FaceInfo, BRICK_FACE_MAP,
    TETRA_FACE_MAP, TRIANGLE_EDGE_MAP,
};
use elmer_solver::{SolverCtx, SolverError};

use crate::quadrature::{gauss_points_tetra_p, gauss_points_triangle_p};

const CALLER: &str = "SetDefaultDirichlet";

/// A bulk element of a face-element (div-conforming) solver: a triangle
/// (303), a tetrahedron (504) or a brick (808).
pub struct FaceCell {
    pub code: u16,
    pub x: [f64; 8],
    pub y: [f64; 8],
    pub z: [f64; 8],
    pub g: [u32; 8],
}

impl FaceCell {
    pub fn of(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Self, SolverError> {
        let el = ctx.mesh.element(e);
        let code = el.etype.0;
        if !matches!(code, 303 | 504 | 808) {
            return Err(SolverError::Unsupported(format!("face elements on element type {code}")));
        }
        let mut c = FaceCell { code, x: [0.0; 8], y: [0.0; 8], z: [0.0; 8], g: [0; 8] };
        for i in 0..c.n() {
            let p = ctx.mesh.coord(el.nodes[i]);
            (c.x[i], c.y[i], c.z[i]) = (p[0], p[1], p[2]);
            c.g[i] = el.nodes[i];
        }
        Ok(c)
    }

    pub fn brick(&self) -> bool {
        self.code == 808
    }

    pub fn n(&self) -> usize {
        match self.code {
            303 => 3,
            504 => 4,
            _ => 8,
        }
    }

    pub fn faces(&self) -> usize {
        match self.code {
            303 => 3,
            504 => 4,
            _ => 6,
        }
    }

    pub fn dim(&self) -> usize {
        if self.code == 303 { 2 } else { 3 }
    }

    /// Dofs per face: RT(0) one; BDM(1) two on an edge, three on a face;
    /// four on a brick's face.
    pub fn per_face(&self, second: bool) -> usize {
        match (self.code, second) {
            (808, _) => 4,
            (_, false) => 1,
            (303, true) => 2,
            _ => 3,
        }
    }

    /// The flux basis functions: the face ones, and a brick's 24 interior
    /// ones (elementwise bubbles before the potential's).
    pub fn flux_dofs(&self, second: bool) -> usize {
        self.faces() * self.per_face(second) + if self.brick() { 24 } else { 0 }
    }

    /// The Gauss rule `LocalMatrix` and `GetFlux` pick for the element.
    pub fn rule(&self) -> Result<GaussRule, SolverError> {
        let r = match self.code {
            303 => gauss_points_triangle_p(3),
            504 => gauss_points_tetra_p(4),
            _ => elmer_fe::quadrature::gauss_points_brick(64).map_err(|e| crate::Error::Unsupported(e.to_string())),
        };
        r.map_err(|m| SolverError::fatal(CALLER, m.to_string()))
    }

    pub fn info(&self, rule: &GaussRule, t: usize, second: bool) -> FaceInfo {
        match self.code {
            303 => {
                let (x, y, z, g) = ([self.x[0], self.x[1], self.x[2]], [self.y[0], self.y[1], self.y[2]], [self.z[0], self.z[1], self.z[2]], [self.g[0], self.g[1], self.g[2]]);
                face_element_info_triangle(&x, &y, &z, &g, rule.u[t], rule.v[t], second, true)
            }
            504 => {
                let (x, y, z, g) = (first4(&self.x), first4(&self.y), first4(&self.z), [self.g[0], self.g[1], self.g[2], self.g[3]]);
                face_element_info_tetra(&x, &y, &z, &g, rule.u[t], rule.v[t], rule.w[t], second, true)
            }
            _ => face_element_info_brick(&self.x, &self.y, &self.z, &self.g, rule.u[t], rule.v[t], rule.w[t], true),
        }
    }

    /// The parent face's nodes in the order of its local face map.
    pub fn face_nodes(&self, face: usize) -> Vec<u32> {
        match self.code {
            303 => TRIANGLE_EDGE_MAP[face - 1].iter().map(|&k| self.g[k - 1]).collect(),
            504 => TETRA_FACE_MAP[face - 1].iter().map(|&k| self.g[k - 1]).collect(),
            _ => BRICK_FACE_MAP[face - 1].iter().map(|&k| self.g[k - 1]).collect(),
        }
    }

    /// `PickActiveFace`: the 1-based local face (an edge of a triangle) the
    /// boundary element with nodes `b` lies on.
    pub fn face_of(&self, b: &[u32]) -> Option<usize> {
        (1..=self.faces()).find(|&q| {
            let f = self.face_nodes(q);
            f.len() == b.len() && f.iter().all(|i| b.contains(i))
        })
    }

    /// `FaceElementOrientation(Parent, ReverseSign, ActiveFaceId)`.
    pub fn reverse(&self, face: usize) -> bool {
        match self.code {
            303 => triangle_face_orientation(&[self.g[0], self.g[1], self.g[2]])[face - 1],
            504 => tetra_face_orientation(&[self.g[0], self.g[1], self.g[2], self.g[3]])[face - 1],
            _ => brick_face_orientation(&self.g)[face - 1],
        }
    }
}

pub fn first4(a: &[f64; 8]) -> [f64; 4] {
    [a[0], a[1], a[2], a[3]]
}


/// The parent of a boundary element and the face it lies on.
pub fn active_face(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Option<(ElemId, FaceCell, usize)>, SolverError> {
    let el = ctx.mesh.element(e);
    let Some(parent) = el.parents[0].or(el.parents[1]) else { return Ok(None) };
    let cell = FaceCell::of(ctx, parent)?;
    Ok(cell.face_of(el.nodes).map(|q| (parent, cell, q)))
}

/// `np`: the nodal dofs of an element with `nd` dofs (`n` with an `n:1`
/// element definition, else none), before its face dofs and its
/// `bubbles` elementwise dofs.
pub fn nodal_dofs(cell: &FaceCell, nd: usize, second: bool, bubbles: usize) -> Result<usize, SolverError> {
    let flux = cell.flux_dofs(second);
    match nd.checked_sub(flux + bubbles) {
        Some(np) if np == 0 || np == cell.n() => Ok(np),
        _ => Err(SolverError::fatal(CALLER, format!("an element with {nd} dofs ({flux} face dofs expected)"))),
    }
}

/// The slots of the face dofs of local face `face` (1-based) of bulk
/// element `parent`, `per_face` of them.
pub fn face_slots(ctx: &SolverCtx<'_>, parent: ElemId, cell: &FaceCell, face: usize, second: bool, bubbles: usize) -> Result<Vec<u32>, SolverError> {
    let mut s = Vec::new();
    ctx.element_slots(parent, &mut s);
    let np = nodal_dofs(cell, s.len(), second, bubbles)?;
    let per_face = cell.per_face(second);
    let first = np + (face - 1) * per_face;
    Ok(s[first..first + per_face].to_vec())
}

/// Elmer's `CrossProduct`.
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], -a[0] * b[2] + a[2] * b[0], a[0] * b[1] - a[1] * b[0]]
}

/// `v / SQRT(SUM(v**2))`.
fn unit(v: [f64; 3]) -> [f64; 3] {
    let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / r, v[1] / r, v[2] / r]
}

/// The values of `key` (and of `key 1..3`) at `nodes`, which are boundary
/// element `e`'s in another order; zero where the condition lacks them.
fn bc_values(ctx: &SolverCtx<'_>, bc: &ParamMap, key: &str, e: ElemId, nodes: &[u32]) -> Result<(Vec<f64>, Vec<[f64; 3]>), SolverError> {
    let orig = ctx.mesh.element(e).nodes;
    let at = |k: &str| -> Result<Vec<f64>, SolverError> {
        let v = ctx.real_at_nodes(bc, k, e)?.unwrap_or_else(|| vec![0.0; orig.len()]);
        Ok(nodes.iter().map(|n| orig.iter().position(|m| m == n).map_or(0.0, |i| v[i])).collect())
    };
    let load = at(key)?;
    let mut vload = vec![[0.0; 3]; nodes.len()];
    for c in 0..3 {
        let v = at(&format!("{key} {}", c + 1))?;
        for i in 0..nodes.len() {
            vload[i][c] = v[i];
        }
    }
    Ok((load, vload))
}

/// `ElementInfo` of the linear boundary element `etype` on `xyz`.
fn element_info(etype: u16, xyz: &[[f64; 3]], u: f64, v: f64, ev: &mut BasisEval) -> Result<(), SolverError> {
    elmer_fe::LagrangeBasis::new(3).eval(ElementType(etype), xyz, [u, v, 0.0], ev).map_err(|m| SolverError::fatal("ElementInfo", m))
}

/// `VectorElementEdgeDOFs(BC, Edge, n, Parent, np, Name, Integral, EDOFs,
/// SecondFamily, FaceElement = .TRUE.)`: the dofs of edge `face` of the
/// triangle `cell` from the condition's normal component. The mesh edge
/// of a boundary edge has its triangle's node order (`FindMeshEdges2D`),
/// so the edge map matches it unreversed.
fn edge_face_dofs(cell: &FaceCell, face: usize, second: bool, load: &[f64], vload: &[[f64; 3]]) -> Result<Vec<f64>, SolverError> {
    let [ia, ib] = TRIANGLE_EDGE_MAP[face - 1];
    let (j, k) = (ia - 1, ib - 1);
    let p = |i: usize| [cell.x[i], cell.y[i], cell.z[i]];
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let mut e = unit(sub(p(k), p(j)));
    // DivConforming: the in-plane normal of the edge.
    let (i1, i2, i3) = (TRIANGLE_EDGE_MAP[0][0] - 1, TRIANGLE_EDGE_MAP[0][1] - 1, TRIANGLE_EDGE_MAP[1][1] - 1);
    let d = unit(cross(sub(p(i2), p(i1)), sub(p(i3), p(i2))));
    e = cross(e, d);
    let rule = elmer_fe::quadrature::gauss_points(ElementType(202), None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let xyz = [p(j), p(k)];
    let s3 = 3f64.sqrt();
    let mut ev = BasisEval::default();
    let mut integral = vec![0.0; if second { 2 } else { 1 }];
    for t in 0..rule.len() {
        element_info(202, &xyz, rule.u[t], 0.0, &mut ev)?;
        let s = rule.s[t] * ev.det_j;
        let l = load[0] * ev.basis[0] + load[1] * ev.basis[1];
        let mut vl = [0.0; 3];
        for c in 0..3 {
            vl[c] = vload[0][c] * ev.basis[0] + vload[1][c] * ev.basis[1];
        }
        let ve = vl[0] * e[0] + vl[1] * e[1] + vl[2] * e[2];
        let u = rule.u[t];
        if second {
            let v = 0.5 * (1.0 - s3 * u);
            integral[0] += s * (l + ve) * v;
            let v = 0.5 * (1.0 + s3 * u);
            integral[1] += s * (l + ve) * v;
        } else {
            integral[0] += s * (l + ve);
        }
    }
    if cell.g[k] < cell.g[j] {
        for x in integral.iter_mut() {
            *x = -*x;
        }
    }
    Ok(integral)
}

/// `NormalVector(Face, Nodes, u, v, .TRUE.)` of a linear triangle: the
/// unit normal of its parametrization, turned away from the tetrahedron
/// `cell` (`CheckNormalDirection` for a face with one parent).
fn outward_triangle_normal(xyz: &[[f64; 3]], cell: &FaceCell) -> [f64; 3] {
    let d = |c: usize, w: [f64; 3]| (w[0] * xyz[0][c] + w[1] * xyz[1][c]) + w[2] * xyz[2][c];
    let (du, dv) = ([-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]);
    let (dxdu, dydu, dzdu) = (d(0, du), d(1, du), d(2, du));
    let (dxdv, dydv, dzdv) = (d(0, dv), d(1, dv), d(2, dv));
    let auu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    let auv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
    let avv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    let det = 1.0 / (auu * avv - auv * auv).sqrt();
    let mut nrm = [(dydu * dzdv - dydv * dzdu) * det, (dxdv * dzdu - dxdu * dzdv) * det, (dxdu * dydv - dxdv * dydu) * det];
    let third = 1.0 / 3.0;
    let centre = |c: usize| (third * xyz[0][c] + third * xyz[1][c]) + third * xyz[2][c];
    let pc = |v: &[f64; 8]| (0..4).fold(0.0, |s, i| s + 0.25 * v[i]);
    let dc = [pc(&cell.x) - centre(0), pc(&cell.y) - centre(1), pc(&cell.z) - centre(2)];
    if dc[0] * nrm[0] + dc[1] * nrm[1] + dc[2] * nrm[2] > 0.0 {
        nrm = nrm.map(|x| -x);
    }
    nrm
}

/// `FaceElementDOFs(BC, Face, n, Parent, FaceId, Name, Integral, FDOFs,
/// SecondFamily)` for a triangular face of a tetrahedron, then the
/// orientation (and, of the second kind, the ordering) the assembly uses.
/// The face's nodes are taken in the parent's face-map order (the second
/// kind's `ElementCopy`; of the first kind only the summation order of a
/// nodal value could tell them from the mesh face's).
fn tetra_face_dofs(cell: &FaceCell, face: usize, second: bool, load: &[f64], vload: &[[f64; 3]], xyz: &[[f64; 3]]) -> Result<Vec<f64>, SolverError> {
    let rule = elmer_fe::quadrature::gauss_points(ElementType(303), Some(3)).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let normal = outward_triangle_normal(xyz, cell);
    let s3 = 3f64.sqrt();
    let mut ev = BasisEval::default();
    let mut integral = vec![0.0; if second { 3 } else { 1 }];
    for t in 0..rule.len() {
        element_info(303, xyz, rule.u[t], rule.v[t], &mut ev)?;
        let mut vl = [0.0; 3];
        for c in 0..3 {
            vl[c] = (vload[0][c] * ev.basis[0] + vload[1][c] * ev.basis[1]) + vload[2][c] * ev.basis[2];
        }
        let l = ((load[0] * ev.basis[0] + load[1] * ev.basis[1]) + load[2] * ev.basis[2])
            + ((vl[0] * normal[0] + vl[1] * normal[1]) + vl[2] * normal[2]);
        let s = rule.s[t] * ev.det_j;
        if second {
            let u = -1.0 + 2.0 * rule.u[t] + rule.v[t];
            let v = s3 * rule.v[t];
            let f = [
                s3 * 0.5 * (1.0 - 2.0 * u + 1.0 / 3.0 - 2.0 / s3 * v),
                s3 * 0.5 * (1.0 + 2.0 * u + 1.0 / 3.0 - 2.0 / s3 * v),
                s3 * (-1.0 / 3.0 + 2.0 / s3 * v),
            ];
            for i in 0..3 {
                integral[i] += l * f[i] * s;
            }
        } else {
            integral[0] += l * s;
        }
    }
    let g4 = [cell.g[0], cell.g[1], cell.g[2], cell.g[3]];
    let sign = if tetra_face_orientation(&g4)[face - 1] { -1.0 } else { 1.0 };
    if second {
        let map = tetra_face_basis_ordering(&g4)[face - 1];
        Ok(map.iter().map(|&k| sign * integral[k - 1]).collect())
    } else {
        Ok(integral.iter().map(|&x| if sign < 0.0 { -x } else { x }).collect())
    }
}

/// `DefaultDirichletBCs` for `<name> {f}` (DefUtils.F90): the face dofs of
/// each boundary element's parent face, set from the condition's `<name>
/// {f}` (along the normal) and `<name> {f} 1..3` (a vector, its normal
/// component) by [`edge_face_dofs`] or [`tetra_face_dofs`]. A boundary with
/// a positive `<name> Condition`, or whose parent the solver does not
/// cover, is skipped. `bubbles` is the elementwise dofs after the face
/// ones in each element.
pub fn face_dirichlet(ctx: &SolverCtx<'_>, name: &str, second: bool, bubbles: usize) -> Result<Vec<(u32, f64)>, SolverError> {
    let prefix = format!("{name} {{f}}");
    let condition = format!("{name} condition");
    let mut fixed = Vec::new();
    let mut slots = Vec::new();
    let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
    for &e in &boundary {
        let Some(bc) = ctx.bc_params(e) else { continue };
        if !bc.0.iter().any(|(k, _)| k.starts_with(&prefix)) {
            continue;
        }
        if let Some(c) = ctx.real_at_nodes(bc, &condition, e)? {
            if c.iter().sum::<f64>() > 0.0 {
                continue;
            }
        }
        let Some((parent, cell, face)) = active_face(ctx, e)? else { continue };
        ctx.element_slots(parent, &mut slots);
        if slots.is_empty() || slots.contains(&NO_DOF) {
            continue;
        }
        let nodes = cell.face_nodes(face);
        let (load, vload) = bc_values(ctx, bc, &prefix, e, &nodes)?;
        let work = match (cell.code, ctx.mesh.element(e).etype.0) {
            (303, 202) => edge_face_dofs(&cell, face, second, &load, &vload)?,
            (504, 303) => {
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                tetra_face_dofs(&cell, face, second, &load, &vload, &xyz)?
            }
            // FaceElementDOFs of a brick's face solves the face mass matrix
            // (`LUSolve`); zero values give zero dofs whatever the solve.
            (808, 404) if load.iter().chain(vload.iter().flatten()).all(|&v| v == 0.0) => vec![0.0; cell.per_face(second)],
            (808, 404) => return Err(SolverError::Unsupported(format!("a nonzero \"{prefix}\" on the faces of bricks"))),
            (p, b) => return Err(SolverError::Unsupported(format!("\"{prefix}\" on a {b} boundary of a {p} element"))),
        };
        for (s, w) in face_slots(ctx, parent, &cell, face, second, bubbles)?.into_iter().zip(work) {
            fixed.push((s, w));
        }
    }
    Ok(fixed)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(code: u16, xyz: &[[f64; 3]]) -> FaceCell {
        let mut c = FaceCell { code, x: [0.0; 8], y: [0.0; 8], z: [0.0; 8], g: [0; 8] };
        for (i, p) in xyz.iter().enumerate() {
            (c.x[i], c.y[i], c.z[i]) = (p[0], p[1], p[2]);
            c.g[i] = i as u32;
        }
        c
    }

    /// RT(0) on the unit triangle: the dof of edge 2-3 is the flux of
    /// `(1, 1)` through it, `|e| (1, 1).n = 2`; BDM(1) splits a constant
    /// flux evenly between the edge's two dofs.
    #[test]
    fn edge_dofs_are_the_normal_flux() {
        let c = cell(303, &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let v = [[1.0, 1.0, 0.0]; 2];
        let rt = edge_face_dofs(&c, 2, false, &[0.0; 2], &v).unwrap();
        assert!((rt[0] - 2.0).abs() < 1e-14, "{rt:?}");
        let bdm = edge_face_dofs(&c, 2, true, &[0.0; 2], &v).unwrap();
        assert!((bdm[0] - 1.0).abs() < 1e-14 && (bdm[1] - 1.0).abs() < 1e-14, "{bdm:?}");
        // The scalar `<name> {f}` is the normal component itself.
        let s = edge_face_dofs(&c, 1, false, &[3.0; 2], &[[0.0; 3]; 2]).unwrap();
        assert!((s[0] - 3.0).abs() < 1e-14, "{s:?}");
    }

    /// The first kind on the unit tetrahedron: face 2-3-4 (area sqrt(3)/2,
    /// outward normal (1,1,1)/sqrt(3)) carries the flux 3/2 of (1, 1, 1),
    /// signed by the face's orientation.
    #[test]
    fn tetra_face_dof_is_the_outward_flux() {
        let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let c = cell(504, &xyz);
        let face = 3;
        let nodes: Vec<[f64; 3]> = TETRA_FACE_MAP[face - 1].iter().map(|&k| xyz[k - 1]).collect();
        let w = tetra_face_dofs(&c, face, false, &[0.0; 3], &[[1.0; 3]; 3], &nodes).unwrap();
        let sign = if tetra_face_orientation(&[0, 1, 2, 3])[face - 1] { -1.0 } else { 1.0 };
        assert!((w[0] - sign * 1.5).abs() < 1e-14, "{w:?}");
        // The second kind's weights live on the p reference triangle (area
        // sqrt(3)): a constant flux gives each of the three dofs flux/sqrt(3).
        let w2 = tetra_face_dofs(&c, face, true, &[0.0; 3], &[[1.0; 3]; 3], &nodes).unwrap();
        assert!(w2.iter().all(|&x| (x - sign * 1.5 / 3f64.sqrt()).abs() < 1e-14), "{w2:?}");
    }
}
