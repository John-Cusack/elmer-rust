//! `FollowCurvedBoundary` (MeshBasics.F90:2893) with `Follow Polynom
//! Boundary = Integer ndeg`: the nodes of the boundary are moved onto a
//! polynomial patch, once when the mesh is loaded (LoadMesh2, MeshLoad.F90:
//! 692) and again after each `Mesh Levels` split (SplitMeshEqual,
//! MeshSplit.F90:1690), so the new boundary nodes lie on the surface.
//!
//! The patch (`PolynomBoundaryFit`, MeshBasics.F90:2580) is fitted on the
//! first call and kept in the boundary condition (`Patch Corners x/y/z`,
//! `Patch Height Basis`), as Elmer stores it in the BC's list; a given
//! `Patch Height Basis` is used as is. Its four corners are the boundary
//! nodes with the smallest sums of element angles, the farthest pair made
//! corners 1 and 3; the height over the corners' plane is fitted in the
//! serendipity p basis of degree `ndeg` on the corners' quadrilateral,
//! weighted by the nodes' boundary measure, with the corner heights imposed.
//! A node is then projected onto the plane and lifted by the patch height
//! there (`PatchElementApply`).
//!
//! `Follow Circle Boundary` (2D boundaries) fits a circle (`CylinderFit`,
//! GeometryFitting.F90:68, kept in the BC as `Cylinder Center` and
//! `Cylinder Radius`, or taken from them) through three boundary nodes: the
//! extremes along the tangent `Tangent1` of the z axis and the node farthest
//! from both; each node is then scaled radially onto it.
//!
//! Refused: the cylinder, sphere, torus and function boundaries, and circles
//! on 3D boundaries.

use elmer_core::{ElementBasis, MeshView, ParamMap, Params, Value};
use elmer_mesh_levels::Mesh as LevelsMesh;
use elmer_solver::{Model, SolverError};

const CALLER: &str = "FollowCurvedBoundary";

/// Whether any boundary condition has a `Follow ...` keyword
/// (`ListCheckPrefixAnyBC(Model, 'Follow')`).
pub(crate) fn any_follow(bcs: &[ParamMap]) -> bool {
    bcs.iter().any(|bc| bc.0.keys().any(|k| k.starts_with("follow")))
}

fn fatal(m: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, m)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `NormalDirection(t1, t2)`: the normalized cross product.
fn normal_direction(t1: [f64; 3], t2: [f64; 3]) -> [f64; 3] {
    let n = cross(t1, t2);
    let s = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    [n[0] / s, n[1] / s, n[2] / s]
}

/// The patch: its corners and the heights of its basis.
struct Patch {
    corners: [[f64; 3]; 4],
    height: Vec<f64>,
    ndeg: usize,
}

impl Patch {
    fn normal(&self) -> [f64; 3] {
        let [c1, c2, _, c4] = self.corners;
        normal_direction(sub(c2, c1), sub(c4, c1))
    }

    /// The serendipity p basis of the corners' quadrilateral at the point
    /// `v2` of its plane: `QuadNodalPBasisAll`, then per edge `j` the
    /// `SD_QuadEdgePBasis(j, k+1, u, v, j == 4)`, `k = 1..ndeg-1`.
    fn basis(&self, v2: [f64; 3]) -> Result<Vec<f64>, SolverError> {
        let shape = elmer_interp::shape::Shape::get(404).ok_or_else(|| fatal("no 404 shape"))?;
        let uvw = elmer_interp::geom::global_to_local(shape, &self.corners, v2, 3);
        let (u, v) = (uvw[0], uvw[1]);
        let edofs = self.ndeg.saturating_sub(1);
        let mut b = vec![0.0; 4 * (1 + edofs)];
        let mut nb = 0;
        elmer_pdofs::elmer_fe_hp::quad::quad_nodal(u, v, &mut b, &mut nb);
        let mut q = 4;
        for j in 1..=4 {
            for k in 1..=edofs {
                let (f, _) = elmer_pdofs::elmer_fe_hp::sd_quad_edge(j, k + 1, u, v, j == 4)
                    .ok_or_else(|| fatal(format!("no SD_QuadEdgePBasis of degree {}", k + 1)))?;
                b[q] = f;
                q += 1;
            }
        }
        Ok(b)
    }

    /// `PatchElementApply`.
    fn apply(&self, v1: [f64; 3]) -> Result<[f64; 3], SolverError> {
        let c1 = self.corners[0];
        let normal = self.normal();
        let norm_proj = dot(sub(v1, c1), normal);
        let mut v2 = [v1[0] - norm_proj * normal[0], v1[1] - norm_proj * normal[1], v1[2] - norm_proj * normal[2]];
        let b = self.basis(v2)?;
        let h: f64 = b.iter().zip(&self.height).map(|(b, h)| b * h).sum();
        for (x, n) in v2.iter_mut().zip(normal) {
            *x += h * n;
        }
        Ok(v2)
    }
}

/// The boundary elements of BC `bc` (1-based).
fn bc_elements<'a>(mesh: &'a LevelsMesh, model: &'a Model, bc: usize) -> impl Iterator<Item = &'a elmer_mesh_levels::Element> {
    mesh.boundary().iter().filter(move |el| model.bc_of_tag(el.tag) == Some(bc))
}

/// `PolynomBoundaryFit` without a stored patch: FindBoundaryCorners,
/// SetBoundaryWeights, FitBoundaryPatch.
fn fit(mesh: &LevelsMesh, model: &Model, bc: usize, ndeg: usize) -> Result<Patch, SolverError> {
    let x = &mesh.coords;
    let n_nodes = x.len();
    let mut active = vec![false; n_nodes];
    let mut angle_sum = vec![0.0f64; n_nodes];
    let mut weights = vec![0.0f64; n_nodes];
    // FindBoundaryCorners.
    for el in bc_elements(mesh, model, bc) {
        let n = (el.etype.0 % 100) as usize;
        if !(3..=4).contains(&n) {
            return Err(SolverError::fatal("PolynomBoundaryFit", "2D polynom can only bet fitted on 2D elements!"));
        }
        for i in 0..n {
            let (i1, i2) = ((i + n - 1) % n, (i + 1) % n);
            let j = el.nodes[i] as usize;
            let mut v1 = sub(x[el.nodes[i1] as usize], x[j]);
            let mut v2 = sub(x[el.nodes[i2] as usize], x[j]);
            let n1 = dot(v1, v1).sqrt();
            let n2 = dot(v2, v2).sqrt();
            v1 = v1.map(|c| c / n1);
            v2 = v2.map(|c| c / n2);
            angle_sum[j] += dot(v1, v2).acos();
            active[j] = true;
        }
    }
    let mut corners = [0usize; 4];
    for c in corners.iter_mut() {
        // MINLOC(AngleSum, Mask = ActiveNode): the first minimum.
        let k = (0..n_nodes)
            .filter(|&i| active[i])
            .fold(None, |m: Option<usize>, i| match m {
                Some(j) if angle_sum[j] <= angle_sum[i] => Some(j),
                _ => Some(i),
            })
            .ok_or_else(|| fatal(format!("Boundary Condition {bc}: no boundary nodes to fit")))?;
        *c = k;
        angle_sum[k] = 3.0 * std::f64::consts::PI;
    }
    let (mut i1, mut i2, mut maxdist) = (0usize, 0usize, 0.0f64);
    for i in 0..4 {
        for j in i + 1..4 {
            let v = sub(x[corners[j]], x[corners[i]]);
            let dist = dot(v, v).sqrt();
            if dist > maxdist {
                (i1, i2, maxdist) = (i, j, dist);
            }
        }
    }
    if i1 != 0 {
        corners.swap(0, i1);
    }
    if i2 != 2 {
        corners.swap(2, i2);
    }
    // SetBoundaryWeights.
    let lb = elmer_fe::LagrangeBasis::new(3);
    let mut ev = elmer_core::BasisEval::default();
    for el in bc_elements(mesh, model, bc) {
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| x[k as usize]).collect();
        let def = elmer_fe::element_type(el.etype.0).map_err(|e| fatal(e.to_string()))?;
        let rule = elmer_fe::quadrature::gauss_points(el.etype, Some(def.gauss_points)).map_err(|e| fatal(e.to_string()))?;
        for i in 0..rule.len() {
            lb.eval(el.etype, &xyz, [rule.u[i], rule.v[i], rule.w[i]], &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
            for (k, &nd) in el.nodes.iter().enumerate() {
                weights[nd as usize] += rule.s[i] * ev.det_j * ev.basis[k];
            }
        }
    }
    // FitBoundaryPatch.
    let patch0 = Patch { corners: corners.map(|c| x[c]), height: Vec::new(), ndeg };
    let edofs = ndeg.saturating_sub(1);
    let nd = 4 * (1 + edofs);
    let c1 = patch0.corners[0];
    let normal = patch0.normal();
    let mut mass = vec![0.0f64; nd * nd]; // column-major MASS(i, j)
    let mut force = vec![0.0f64; nd];
    let mut dir = [0.0f64; 4];
    for i in 0..n_nodes {
        if !active[i] {
            continue;
        }
        let v1 = x[i];
        let norm_proj = dot(sub(v1, c1), normal);
        let v2 = [v1[0] - norm_proj * normal[0], v1[1] - norm_proj * normal[1], v1[2] - norm_proj * normal[2]];
        let basis = patch0.basis(v2)?;
        let weight = weights[i];
        for q in 0..nd {
            for r in 0..nd {
                mass[r + q * nd] += weight * basis[r] * basis[q];
            }
        }
        for r in 0..nd {
            force[r] += weight * basis[r] * norm_proj;
        }
        for (j, &c) in corners.iter().enumerate() {
            if c == i {
                dir[j] = norm_proj;
            }
        }
    }
    for j in 0..4 {
        for c in 0..nd {
            mass[j + c * nd] = 0.0;
        }
        mass[j + j * nd] = 1.0;
        force[j] = dir[j];
    }
    let mut pivot = vec![0usize; nd];
    elmer_fe::lu_decomp(&mut mass, nd, &mut pivot).ok_or_else(|| SolverError::fatal("FitBoundaryPatch", "LU-decomposition fails"))?;
    elmer_fe::lu_solve_pivot(&mass, nd, &mut force, &pivot);
    Ok(Patch { height: force, ..patch0 })
}

fn reals(p: &ParamMap, key: &str) -> Option<Vec<f64>> {
    match p.get(key)? {
        Value::Real { data, .. } => Some(data.clone()),
        _ => None,
    }
}

fn column(v: &[f64]) -> Value {
    Value::Real { rows: v.len(), cols: 1, data: v.to_vec() }
}

/// `FollowCurvedBoundary(Model, Mesh, SetP = .FALSE.)` on `mesh`: each BC
/// with `Follow Polynom Boundary` moves its nodes onto its patch, fitted on
/// the first call and kept in the BC.
pub(crate) fn follow_curved_boundary(mesh: &mut LevelsMesh, model: &mut Model) -> Result<(), SolverError> {
    if !any_follow(&model.boundary_conditions) {
        return Ok(());
    }
    for b in 0..model.boundary_conditions.len() {
        let bc = &model.boundary_conditions[b];
        if bc.logical("follow circle boundary") == Some(true) {
            follow_circle(mesh, model, b + 1)?;
            continue;
        }
        for mode in ["follow cylinder boundary", "follow sphere boundary", "follow function boundary", "follow toroid boundary"] {
            if bc.logical(mode) == Some(true) {
                return Err(SolverError::Unsupported(format!("Boundary Condition {}: {mode}", b + 1)));
            }
        }
        if bc.get("follow polynom boundary").is_none() {
            continue;
        }
        let ndeg = bc.integer("follow polynom boundary").unwrap_or(0).max(0) as usize;
        if reals(bc, "patch height basis").is_none() {
            let patch = fit(mesh, model, b + 1, ndeg)?;
            let bc = &mut model.boundary_conditions[b];
            for (c, key) in ["Patch Corners x", "Patch Corners y", "Patch Corners z"].iter().enumerate() {
                bc.insert(key, column(&patch.corners.map(|p| p[c])));
            }
            bc.insert("Patch Height Basis", column(&patch.height));
        }
        let bc = &model.boundary_conditions[b];
        let height = reals(bc, "patch height basis").expect("set above");
        let get = |k: &str| reals(bc, k).filter(|v| v.len() >= 4).ok_or_else(|| fatal(format!("{k} not found")));
        let (cx, cy, cz) = (get("patch corners x")?, get("patch corners y")?, get("patch corners z")?);
        let patch = Patch { corners: std::array::from_fn(|i| [cx[i], cy[i], cz[i]]), height, ndeg };
        // SetCurvedBoundary, `.NOT. SetP`: every node of the BC's elements.
        let mut done = vec![false; mesh.coords.len()];
        for el in bc_elements(mesh, model, b + 1) {
            for &k in &el.nodes {
                done[k as usize] = true;
            }
        }
        for (j, d) in done.iter().enumerate() {
            if *d {
                mesh.coords[j] = patch.apply(mesh.coords[j])?;
            }
        }
    }
    Ok(())
}

/// `Det3x3` (LinearAlgebra.F90:883), `a[row][col]`.
fn det3x3(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// `TangentDirections(Normal, Tangent1, Tangent2)` (ElementUtils.F90:2634).
fn tangent_directions(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut t1 = if n[0].abs() <= n[2].abs() && n[1].abs() <= n[2].abs() { [0.0, -n[2], n[1]] } else { [-n[1], n[0], 0.0] };
    let s = dot(t1, t1).sqrt();
    t1 = t1.map(|c| c / s);
    let mut t2 = [n[1] * t1[2] - n[2] * t1[1], n[2] * t1[0] - n[0] * t1[2], n[0] * t1[1] - n[1] * t1[0]];
    let s = dot(t2, t2).sqrt();
    t2 = t2.map(|c| c / s);
    (t1, t2)
}

/// `CylinderFit(Mesh, BC, bc, 2, FitParams)` then `SetCurvedBoundary` in
/// mode 1 (circle).
fn follow_circle(mesh: &mut LevelsMesh, model: &mut Model, bc: usize) -> Result<(), SolverError> {
    let x = &mesh.coords;
    // cdim from the first element of the BC.
    let Some(first) = bc_elements(mesh, model, bc).next() else { return Ok(()) };
    if first.etype.0 >= 300 {
        return Err(SolverError::Unsupported(format!("Boundary Condition {bc}: Follow Circle Boundary on a 3D boundary")));
    }
    let mut active = vec![false; x.len()];
    for el in bc_elements(mesh, model, bc) {
        for &k in &el.nodes {
            active[k as usize] = true;
        }
    }
    let params = &model.boundary_conditions[bc - 1];
    let center = reals(params, "cylinder center");
    let radius = params.real("cylinder radius");
    let (orig, rad) = match (center, radius) {
        (Some(c), Some(r)) if c.len() >= 2 => ([c[0], c[1]], r),
        _ => {
            let axis = [0.0, 0.0, 1.0];
            let (t1, t2) = tangent_directions(axis);
            let (mut min_d, mut max_d) = (f64::MAX, -f64::MAX);
            let mut ind = [0usize; 3];
            for (i, p) in x.iter().enumerate() {
                if !active[i] {
                    continue;
                }
                let d1 = dot(t1, *p);
                if d1 < min_d {
                    min_d = d1;
                    ind[0] = i;
                }
                if d1 > max_d {
                    max_d = d1;
                    ind[1] = i;
                }
            }
            let mut cc = [[0.0f64; 3]; 3];
            for j in 0..2 {
                let p = x[ind[j]];
                cc[j] = [dot(t1, p), dot(t2, p), dot(axis, p)];
            }
            let mut maxdist = 0.0;
            for (i, p) in x.iter().enumerate() {
                if !active[i] {
                    continue;
                }
                let (d1, d2) = (dot(t1, *p), dot(t2, *p));
                let mut mind = f64::MAX;
                for c in cc.iter().take(2) {
                    mind = mind.min((d1 - c[0]).powi(2) + (d2 - c[1]).powi(2));
                }
                if maxdist < mind {
                    maxdist = mind;
                    ind[2] = i;
                }
            }
            let p = x[ind[2]];
            cc[2] = [dot(t1, p), dot(t2, p), dot(axis, p)];
            let sq = |k: usize| cc[k][0] * cc[k][0] + cc[k][1] * cc[k][1];
            let m11 = det3x3(std::array::from_fn(|k| [cc[k][0], cc[k][1], 1.0]));
            let m12 = det3x3(std::array::from_fn(|k| [sq(k), cc[k][1], 1.0]));
            let m13 = det3x3(std::array::from_fn(|k| [sq(k), cc[k][0], 1.0]));
            let m14 = det3x3(std::array::from_fn(|k| [sq(k), cc[k][0], cc[k][1]]));
            if m11.abs() < f64::EPSILON {
                return Err(SolverError::fatal("CylinderFit", "Points cannot be an a circle"));
            }
            let x0 = 0.5 * m12 / m11;
            let y0 = -0.5 * m13 / m11;
            let rad = (x0 * x0 + y0 * y0 + m14 / m11).sqrt();
            let c: [f64; 3] = std::array::from_fn(|a| x0 * t1[a] + y0 * t2[a]);
            let bcp = &mut model.boundary_conditions[bc - 1];
            bcp.insert("Cylinder Center", column(&c));
            bcp.insert("Cylinder Radius", Value::Real { rows: 1, cols: 1, data: vec![rad] });
            ([c[0], c[1]], rad)
        }
    };
    for (j, p) in mesh.coords.iter_mut().enumerate() {
        if !active[j] {
            continue;
        }
        let c = [p[0] - orig[0], p[1] - orig[1]];
        let rat = rad / (c[0] * c[0] + c[1] * c[1]).sqrt();
        *p = [rat * c[0] + orig[0], rat * c[1] + orig[1], p[2]];
    }
    Ok(())
}

/// The loaded mesh as a [`LevelsMesh`] followed onto its curved boundaries,
/// then `Mesh Levels - 1` equal splits, each followed again.
pub(crate) fn follow_levels(mesh: &dyn MeshView, model: &mut Model) -> Result<LevelsMesh, SolverError> {
    let err = |e: elmer_mesh_levels::MeshLevelsError| match e {
        elmer_mesh_levels::MeshLevelsError::Unsupported(m) => SolverError::Unsupported(m),
        elmer_mesh_levels::MeshLevelsError::Invalid(m) => SolverError::fatal("LoadModel", m),
    };
    let mut cur = LevelsMesh::from_view(mesh);
    follow_curved_boundary(&mut cur, model)?;
    let levels = model.simulation.integer("mesh levels").unwrap_or(1);
    if levels > 1 {
        // `mesh_levels`' own checks, one split at a time, on BCs without the
        // `Follow` keywords it refuses.
        let mut sim = model.simulation.clone();
        sim.insert("Mesh Levels", Value::Integer(vec![2]));
        let bcs: Vec<ParamMap> = model
            .boundary_conditions
            .iter()
            .map(|bc| ParamMap(bc.0.iter().filter(|(k, _)| !k.starts_with("follow")).map(|(k, v)| (k.clone(), v.clone())).collect()))
            .collect();
        for _ in 1..levels {
            cur = elmer_mesh_levels::mesh_levels(&cur, &sim, &bcs).map_err(err)?.ok_or_else(|| fatal("no split"))?;
            follow_curved_boundary(&mut cur, model)?;
        }
    }
    Ok(cur)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_flat_patch_projects_onto_the_corners_plane_and_a_lift_follows_the_basis() {
        let corners = [[0.0, 0.0, 1.0], [2.0, 0.0, 1.0], [2.0, 2.0, 1.0], [0.0, 2.0, 1.0]];
        let flat = Patch { corners, height: vec![0.0; 8], ndeg: 2 };
        assert_eq!(flat.normal(), [0.0, 0.0, 1.0]);
        let p = flat.apply([0.5, 1.5, 7.0]).unwrap();
        assert!((p[0] - 0.5).abs() < 1e-12 && (p[1] - 1.5).abs() < 1e-12 && (p[2] - 1.0).abs() < 1e-12);
        // Corner heights lift by the bilinear nodal basis: the centre by the mean.
        let lifted = Patch { corners, height: vec![1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], ndeg: 2 };
        let c = lifted.apply([1.0, 1.0, 0.0]).unwrap();
        assert!((c[2] - 2.0).abs() < 1e-12, "{c:?}");
    }
}
