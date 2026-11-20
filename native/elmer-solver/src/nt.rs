//! Normal-tangential systems (`Normal-Tangential <var>`, `Normal-Tangential
//! Velocity` for a flow solution): the vector variable is solved in a
//! coordinate system rotated, at the nodes of such boundaries, to the
//! averaged boundary normal and its tangents.
//!
//! - [`NtSystem::build`]: `CheckNormalTangentialBoundary` and
//!   `AverageBoundaryNormals` (BoundaryConditionUtils.F90:3420, :3699), run
//!   by `DefaultInitialize` (DefUtils.F90:3910-3943).
//! - [`rotate_local`]: `RotateMatrix` (ElementUtils.F90:2539) as
//!   `UpdateGlobalEquations` applies it (SolverBasics.F90:700-727).
//! - [`NtSystem::rotate_vector`]: `RotateNTSystem`, for the Dirichlet rows
//!   of Cartesian components on rotated nodes.
//! - [`NtSystem::back_rotate`]: `BackRotateNTSystem` after the solve.
//!
//! Refused, not approximated: conditional n-t boundaries (`<name>
//! Condition`), `Rotational Normals`, `Mass Consistent Normals`, a given
//! `Normals` variable, `Normal-Tangential Priority`, `Use p Normals`,
//! normal-tangential systems on mortar or contact boundaries in 3D (the
//! left-hand tangent system), and periodic projectors that are not nodal.

use elmer_core::{ElemId, LocalSystem, MeshView, ParamMap, Params, Value};

use crate::backend::NodalProjectorView;
use crate::error::SolverError;
use crate::model::Model;

/// `Solver % NormalTangential`.
#[derive(Clone, Debug)]
pub(crate) struct NtSystem {
    /// `NormalTangentialName`.
    pub name: String,
    /// `BoundaryReorder`: mesh node -> 1-based n-t index, 0 for none.
    pub reorder: Vec<u32>,
    pub normals: Vec<[f64; 3]>,
    pub tangent1: Vec<[f64; 3]>,
    pub tangent2: Vec<[f64; 3]>,
}

fn unsupported(what: impl Into<String>) -> SolverError {
    SolverError::Unsupported(what.into())
}

/// The n-t keyword of a variable (DefUtils.F90:3911-3916).
pub(crate) fn nt_name(var: &str) -> String {
    if var.eq_ignore_ascii_case("flow solution") {
        "normal-tangential velocity".into()
    } else {
        format!("normal-tangential {}", elmer_core::normalize_key(var))
    }
}

/// `TangentDirections(Normal, Tangent1, Tangent2)`.
pub fn tangent_directions(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let (n1, n2, n3) = (n[0].abs(), n[1].abs(), n[2].abs());
    let mut t1 = if n1 <= n3 && n2 <= n3 { [0.0, -n[2], n[1]] } else { [-n[1], n[0], 0.0] };
    let s = (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]).sqrt();
    t1 = t1.map(|x| x / s);
    let mut t2 = [
        n[1] * t1[2] - n[2] * t1[1],
        n[2] * t1[0] - n[0] * t1[2],
        n[0] * t1[1] - n[1] * t1[0],
    ];
    let s = (t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2]).sqrt();
    t2 = t2.map(|x| x / s);
    (t1, t2)
}

/// `NormalVector(Element, Nodes, u, v, .TRUE.)` of a linear boundary line
/// (2D), triangle or quadrilateral (3D) at reference point `(u, v)`, turned
/// away from its parent (`CheckNormalDirection`, ElementGeometry.F90): the
/// parent with the larger body id of two, the only one otherwise; turned
/// when the vector from the point to the parent's centre has a positive
/// component along it. `None`: no unit normal.
pub fn boundary_normal(
    mesh: &dyn MeshView,
    model: &Model,
    e: ElemId,
    uv: [f64; 2],
) -> Result<[f64; 3], SolverError> {
    let el = mesh.element(e);
    let x: Vec<[f64; 3]> = el.nodes.iter().map(|&n| mesh.coord(n)).collect();
    let (u, v) = (uv[0], uv[1]);
    // NodalFirstDerivatives of the linear element maps, summed in node order
    // (DOT_PRODUCT).
    let dot = |d: &[f64], a: usize| d.iter().zip(&x).fold(0.0, |s, (di, xi)| s + di * xi[a]);
    let mut n = match el.etype.0 {
        202 | 203 => {
            // NodalFirstDerivatives of the line: (u^2 -+ u)/2 and 1 - u^2.
            let d: Vec<f64> = if el.etype.0 == 202 { vec![-0.5, 0.5] } else { vec![u - 0.5, u + 0.5, -2.0 * u] };
            let (dxdu, dydu) = (dot(&d, 0), dot(&d, 1));
            let det = dxdu * dxdu + dydu * dydu;
            if det <= 0.0 {
                return Ok([0.0; 3]);
            }
            let det = 1.0 / det.sqrt();
            [-dydu * det, dxdu * det, 0.0]
        }
        303 | 404 => {
            let (du, dv): (Vec<f64>, Vec<f64>) = if el.etype.0 == 303 {
                (vec![-1.0, 1.0, 0.0], vec![-1.0, 0.0, 1.0])
            } else {
                (
                    vec![-(1.0 - v) / 4.0, (1.0 - v) / 4.0, (1.0 + v) / 4.0, -(1.0 + v) / 4.0],
                    vec![-(1.0 - u) / 4.0, -(1.0 + u) / 4.0, (1.0 + u) / 4.0, (1.0 - u) / 4.0],
                )
            };
            let (xu, yu, zu) = (dot(&du, 0), dot(&du, 1), dot(&du, 2));
            let (xv, yv, zv) = (dot(&dv, 0), dot(&dv, 1), dot(&dv, 2));
            let auu = xu * xu + yu * yu + zu * zu;
            let auv = xu * xv + yu * yv + zu * zv;
            let avv = xv * xv + yv * yv + zv * zv;
            let det = 1.0 / (auu * avv - auv * auv).sqrt();
            [(yu * zv - yv * zu) * det, (xv * zu - xu * zv) * det, (xu * yv - xv * yu) * det]
        }
        t => return Err(unsupported(format!("normal-tangential boundaries on boundary elements of type {t}"))),
    };
    // CheckNormalDirection from the element's centre.
    let nn = x.len() as f64;
    let point: [f64; 3] = std::array::from_fn(|a| x.iter().fold(0.0, |s, xi| s + xi[a] / nn));
    let [left, right] = el.parents;
    let parent = match (left, right) {
        (Some(l), Some(r)) => {
            let body = |p: ElemId| model.body_of_tag(mesh.element(p).tag);
            match body(l).cmp(&body(r)) {
                std::cmp::Ordering::Greater => Some(l),
                std::cmp::Ordering::Less => Some(r),
                std::cmp::Ordering::Equal => {
                    return Err(unsupported("normal direction on a boundary between two elements of one body"))
                }
            }
        }
        (Some(l), None) => Some(l),
        (None, r) => r,
    };
    if let Some(p) = parent {
        let pe = mesh.element(p);
        let corners = match pe.etype.family() {
            3 => 3,
            4 => 4,
            5 => 4,
            // The wedge at (1/3, 1/3, 0): the mean of its six nodes.
            7 => 6,
            8 => 8,
            f => return Err(unsupported(format!("normal direction from a parent element of family {f}"))),
        };
        let c: [f64; 3] = std::array::from_fn(|a| {
            pe.nodes[..corners].iter().map(|&nn| mesh.coord(nn)[a]).sum::<f64>() / corners as f64
        });
        let d: f64 = (0..3).map(|a| (c[a] - point[a]) * n[a]).sum();
        if d > 0.0 {
            n = n.map(|c| -c);
        }
    }
    Ok(n)
}

/// Reference coordinates of the nodes of a linear boundary element.
fn node_uv(etype: u16, j: usize) -> Result<[f64; 2], SolverError> {
    Ok(match etype {
        202 => [[-1.0, 0.0], [1.0, 0.0]][j],
        203 => [[-1.0, 0.0], [1.0, 0.0], [0.0, 0.0]][j],
        303 => [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]][j],
        404 => [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]][j],
        t => return Err(unsupported(format!("normal-tangential boundaries on boundary elements of type {t}"))),
    })
}

impl NtSystem {
    /// `CheckNormalTangentialBoundary` + `AverageBoundaryNormals` for the
    /// n-t keyword `name`; `None` when no boundary rotates.
    pub(crate) fn build(
        mesh: &dyn MeshView,
        model: &Model,
        bc_of: &dyn Fn(ElemId) -> Option<usize>,
        name: &str,
        dim: usize,
        projectors: &[NodalProjectorView],
        fe: Option<&crate::backend::FeBackend>,
    ) -> Result<Option<Self>, SolverError> {
        let bcs = &model.boundary_conditions;
        let rotate = |bc: &ParamMap| bc.logical(&format!("{name} rotate")).unwrap_or(true);
        let active = |bc: &ParamMap| bc.logical(name) == Some(true) && rotate(bc);
        if !bcs.iter().any(active) {
            return Ok(None);
        }
        for bc in bcs.iter().filter(|b| b.logical(name) == Some(true)) {
            for (key, what) in [
                (format!("{name} condition"), "conditional normal-tangential boundaries"),
                ("rotational normals".to_string(), "Rotational Normals"),
                ("normal-tangential priority".to_string(), "Normal-Tangential Priority"),
            ] {
                if bc.get(&key).is_some() {
                    return Err(unsupported(what));
                }
            }
            if dim == 3 && (bc.get("mortar bc").is_some() || bc.get("contact bc").is_some()) {
                return Err(unsupported("normal-tangential mortar or contact boundaries in 3D"));
            }
        }
        if model.solvers.iter().any(|s| s.logical("use p normals") == Some(true)) {
            return Err(unsupported("Use p Normals"));
        }
        if model.simulation.logical("use cylinder system").is_some() || model.simulation.logical("use lhs system") == Some(true) {
            return Err(unsupported("Use Cylinder System / Use Lhs System"));
        }
        let nb = mesh.n_bulk();
        let boundary: Vec<ElemId> = (nb..nb + mesh.n_boundary()).map(|t| ElemId(t as u32)).collect();
        let bc_params = |e: ElemId| bc_of(e).and_then(|i| bcs.get(i - 1));

        // CheckNormalTangentialBoundary: the n-t index of each node, in
        // boundary-element order.
        let mut reorder = vec![0u32; mesh.n_nodes()];
        let mut count = 0u32;
        for &e in &boundary {
            let el = mesh.element(e);
            if el.etype.0 == 101 {
                continue;
            }
            let Some(bc) = bc_params(e) else { continue };
            if !active(bc) {
                continue;
            }
            for &n in el.nodes {
                if reorder[n as usize] == 0 {
                    count += 1;
                    reorder[n as usize] = count;
                }
            }
        }
        if count == 0 {
            return Ok(None);
        }
        let m = count as usize;
        let mut normals = vec![[0.0f64; 3]; m];
        // AverageBoundaryNormals: the sum of the element normals at the nodes.
        for &e in &boundary {
            let el = mesh.element(e);
            if el.etype.0 < 200 {
                continue;
            }
            let Some(bc) = bc_params(e) else { continue };
            if !active(bc) {
                continue;
            }
            // `Mass Consistent Normals`: IntegMassConsistent
            // (BoundaryConditionUtils.F90:4374), the element normal weighted
            // by the node's basis function over the element's Gauss points.
            let mass_consistent = bc.logical("mass consistent normals") == Some(true);
            let ips = if mass_consistent {
                let fe = fe.ok_or_else(|| unsupported("Mass Consistent Normals without an element basis"))?;
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| mesh.coord(n)).collect();
                let mut out = Vec::new();
                let mut be = elmer_core::BasisEval::default();
                for qp in fe.quadrature.gauss_points(el.etype, None) {
                    fe.basis
                        .eval(el.etype, &xyz, qp.uvw, &mut be)
                        .map_err(|m| SolverError::fatal("IntegMassConsistent", m.to_string()))?;
                    let lnrm = boundary_normal(mesh, model, e, [qp.uvw[0], qp.uvw[1]])?;
                    out.push((qp.weight, lnrm, be.det_j, be.basis.clone()));
                }
                out
            } else {
                Vec::new()
            };
            for (j, &n) in el.nodes.iter().enumerate() {
                let k = reorder[n as usize];
                if k == 0 {
                    continue;
                }
                let nrm = if mass_consistent {
                    let mut nrm = [0.0f64; 3];
                    for (s, lnrm, det_j, basis) in &ips {
                        for a in 0..3 {
                            nrm[a] += s * lnrm[a] * det_j * basis[j];
                        }
                    }
                    nrm
                } else {
                    boundary_normal(mesh, model, e, node_uv(el.etype.0, j)?)?
                };
                let s = &mut normals[k as usize - 1];
                for a in 0..3 {
                    s[a] += nrm[a];
                }
            }
        }
        // The nodal periodic projectors: rotated ones flip the normals of
        // their target nodes; the others average the normals over the pair.
        let ri = |node: u32| -> usize { reorder.get(node as usize).copied().unwrap_or(0) as usize };
        for p in projectors {
            let bc = &bcs[p.bc - 1];
            if bc.logical("rotational normals") == Some(true) {
                continue;
            }
            if periodic_rotate(bc)? {
                let mut done = vec![false; m];
                for i in 0..p.rows.len() - 1 {
                    if ri(p.inv_perm[i]) == 0 {
                        continue;
                    }
                    for &col in &p.cols[p.rows[i]..p.rows[i + 1]] {
                        if col == elmer_core::NO_DOF {
                            continue;
                        }
                        let mm = ri(col);
                        if mm > 0 && !done[mm - 1] {
                            done[mm - 1] = true;
                            normals[mm - 1] = normals[mm - 1].map(|c| -c);
                        }
                    }
                }
                continue;
            }
            if bc.logical("one sided normals") != Some(true) {
                for i in 0..p.rows.len() - 1 {
                    let k = ri(p.inv_perm[i]);
                    if k == 0 {
                        continue;
                    }
                    for l in p.rows[i]..p.rows[i + 1] {
                        if p.cols[l] == elmer_core::NO_DOF {
                            continue;
                        }
                        let mm = ri(p.cols[l]);
                        if mm > 0 {
                            let add = normals[mm - 1];
                            for a in 0..3 {
                                normals[k - 1][a] += p.values[l] * add[a];
                            }
                        }
                    }
                }
            }
            for i in 0..p.rows.len() - 1 {
                if ri(p.inv_perm[i]) == 0 {
                    continue;
                }
                for &col in &p.cols[p.rows[i]..p.rows[i + 1]] {
                    if col == elmer_core::NO_DOF {
                        continue;
                    }
                    let mm = ri(col);
                    if mm > 0 {
                        normals[mm - 1] = [0.0; 3];
                    }
                }
            }
        }
        // The master side's values on the slave side (hierarchical).
        for p in projectors {
            let bc = &bcs[p.bc - 1];
            if bc.logical("rotational normals") == Some(true) || periodic_rotate(bc)? {
                continue;
            }
            for i in 0..p.rows.len() - 1 {
                let k = ri(p.inv_perm[i]);
                if k == 0 {
                    continue;
                }
                for l in p.rows[i]..p.rows[i + 1] {
                    if p.cols[l] == elmer_core::NO_DOF {
                        continue;
                    }
                    let mm = ri(p.cols[l]);
                    if mm > 0 {
                        let add = normals[k - 1];
                        for a in 0..3 {
                            normals[mm - 1][a] += p.values[l] * add[a];
                        }
                    }
                }
            }
        }
        // Normalize, and the tangents in 3D.
        let mut tangent1 = vec![[0.0f64; 3]; m];
        let mut tangent2 = vec![[0.0f64; 3]; m];
        for &k in &reorder {
            let k = k as usize;
            if k == 0 {
                continue;
            }
            let nk = normals[k - 1];
            let s = (nk[0] * nk[0] + nk[1] * nk[1] + nk[2] * nk[2]).sqrt();
            if s > f64::MIN_POSITIVE {
                normals[k - 1] = nk.map(|c| c / s);
                if dim > 2 {
                    let (t1, t2) = tangent_directions(normals[k - 1]);
                    tangent1[k - 1] = t1;
                    tangent2[k - 1] = t2;
                }
            }
        }
        if model.simulation.logical("save averaged normals") == Some(true) {
            return Err(unsupported("Save Averaged Normals"));
        }
        Ok(Some(NtSystem { name: name.to_string(), reorder, normals, tangent1, tangent2 }))
    }

    /// `RotateNTSystem(Vec, node)`: a Cartesian vector in the n-t system of
    /// mesh node `node` (unchanged off the n-t boundary).
    pub(crate) fn rotate_vector(&self, v: [f64; 3], node: u32, dim: usize) -> [f64; 3] {
        let k = self.reorder.get(node as usize).copied().unwrap_or(0) as usize;
        if k == 0 {
            return v;
        }
        let n = self.normals[k - 1];
        if dim < 3 {
            let (bu, bv) = (v[0], v[1]);
            [n[0] * bu + n[1] * bv, -n[1] * bu + n[0] * bv, v[2]]
        } else {
            let (t1, t2) = (self.tangent1[k - 1], self.tangent2[k - 1]);
            let (bu, bv, bw) = (v[0], v[1], v[2]);
            [
                n[0] * bu + n[1] * bv + n[2] * bw,
                t1[0] * bu + t1[1] * bv + t1[2] * bw,
                t2[0] * bu + t2[1] * bv + t2[2] * bw,
            ]
        }
    }

    /// `BackRotateNTSystem(Solution, Perm, NDOFs)`.
    pub(crate) fn back_rotate(&self, x: &mut [f64], perm: &[u32], ndofs: usize, dim: usize) {
        if ndofs < dim {
            return;
        }
        for (i, &k) in self.reorder.iter().enumerate() {
            if k == 0 {
                continue;
            }
            let Some(&j) = perm.get(i) else { continue };
            if j == elmer_core::NO_DOF {
                continue;
            }
            let (k, j) = (k as usize - 1, j as usize);
            let n = self.normals[k];
            if dim < 3 {
                let (bu, bv) = (x[ndofs * j], x[ndofs * j + 1]);
                x[ndofs * j] = n[0] * bu - n[1] * bv;
                x[ndofs * j + 1] = n[1] * bu + n[0] * bv;
            } else {
                let (t1, t2) = (self.tangent1[k], self.tangent2[k]);
                let (bu, bv, bw) = (x[ndofs * j], x[ndofs * j + 1], x[ndofs * j + 2]);
                // RM(1,:) = N, RM(2,:) = T1, RM(3,:) = T2; x = RM^T b.
                x[ndofs * j] = n[0] * bu + t1[0] * bv + t2[0] * bw;
                x[ndofs * j + 1] = n[1] * bu + t1[1] * bv + t2[1] * bw;
                x[ndofs * j + 2] = n[2] * bu + t1[2] * bv + t2[2] * bw;
            }
        }
    }
}

/// `Periodic BC Rotate` given with a nonzero angle.
fn periodic_rotate(bc: &ParamMap) -> Result<bool, SolverError> {
    match bc.get("periodic bc rotate") {
        None => Ok(false),
        Some(Value::Real { data, .. }) => Ok(data.iter().any(|&x| x != 0.0)),
        Some(_) => Err(unsupported("Periodic BC Rotate other than constant reals with normal-tangential boundaries")),
    }
}

/// `RotateMatrix(Matrix, Vector, n, dim, DOFs, NormalIndexes, ...)`: for
/// each local node `i` with n-t index `normal[i] > 0`, `M = R M R^T`, `F =
/// R F`, R the identity but for that node's block, the rows N (and T1, T2 in
/// 3D) or (n1, n2; -n2, n1) in 2D. Dense products in the Fortran's order.
pub(crate) fn rotate_local(local: &mut LocalSystem, normal: &[u32], nt: &NtSystem, dim: usize, dofs: usize) {
    let size = local.n;
    let nodes = (size / dofs).min(normal.len());
    let mut r = vec![0.0f64; size * size];
    let mut q = vec![0.0f64; size * size];
    for (i, &ii) in normal.iter().enumerate().take(nodes) {
        let ii = ii as usize;
        if ii == 0 || ii > nt.normals.len() {
            continue;
        }
        r.iter_mut().for_each(|x| *x = 0.0);
        for j in 0..size {
            r[j * size + j] = 1.0;
        }
        let n1 = nt.normals[ii - 1];
        let b = dofs * i;
        match dim {
            2 => {
                r[b * size + b] = n1[0];
                r[b * size + b + 1] = n1[1];
                r[(b + 1) * size + b] = -n1[1];
                r[(b + 1) * size + b + 1] = n1[0];
            }
            _ => {
                let (t1, t2) = (nt.tangent1[ii - 1], nt.tangent2[ii - 1]);
                for (row, vecr) in [n1, t1, t2].iter().enumerate() {
                    for c in 0..3 {
                        r[(b + row) * size + b + c] = vecr[c];
                    }
                }
            }
        }
        let m = &mut local.stiff;
        for j in 0..size {
            for k in 0..size {
                let mut s = 0.0;
                for l in 0..size {
                    s += r[j * size + l] * m[l * size + k];
                }
                q[j * size + k] = s;
            }
        }
        for j in 0..size {
            for k in 0..size {
                let mut s = 0.0;
                for l in 0..size {
                    s += q[j * size + l] * r[k * size + l];
                }
                m[j * size + k] = s;
            }
        }
        let f = &mut local.force;
        let mut qv = vec![0.0; size];
        for (j, qj) in qv.iter_mut().enumerate() {
            let mut s = 0.0;
            for k in 0..size {
                s += r[j * size + k] * f[k];
            }
            *qj = s;
        }
        f.copy_from_slice(&qv);
    }
}
