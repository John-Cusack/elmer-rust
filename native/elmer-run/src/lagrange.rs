//! `Lagrange Element Degree` of the VTU output (VtuOutputSolver.F90 `LagN`):
//! p fields written on order-`LagN` Lagrange cells.
//!
//! * `GetLagrangeIndexes` (MeshBasics.F90): each bulk element's Lagrange
//!   points are its corners, then `LagN - 1` points per edge numbered from
//!   the global edge (reversed when the orientations do not match, with
//!   the Fortran's own test), then the interior points. The points used by
//!   bulk elements are numbered in that global order (`GenerateSavePermutation`).
//! * `HierarchicPToLagrange` (PElementVisual.F90:51): the VTK point tables
//!   in the p reference element, where the p solution is evaluated with the
//!   element's p basis.
//! * `p2LagrangeSwapper` (SolverBasics.F90): element by element, each point
//!   takes the value of the last element that visits it. Coordinates go
//!   through the same swapper with the p definitions of the first p solver;
//!   only their nodal values are known, so they are bilinear (linear) maps.
//!
//! Only 2D meshes (triangles and quadrilaterals, boundary lines) and nodal
//! p fields of a p solver are supported.

use std::sync::Arc;

use elmer_core::{BasisEval, ElemId, ElementRef, ElementType, FieldKind, FieldRef, FieldSource, MeshView, NO_DOF};
use elmer_solver::{SolverDofs, SolverError, VarHandle, VariableStore};

const CALLER: &str = "VtuOutputSolver";
const MAX_DEGREE: usize = 8;

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("VTU output: Lagrange Element Degree with {what}"))
}

/// `NodesOnEdge`: `k - 1` points from corner `a` to corner `b` of `corners`,
/// by repeated steps of `2/k` along the unit direction.
fn nodes_on_edge(corners: &[(f64, f64)], a: usize, b: usize, k: usize, out: &mut Vec<(f64, f64)>) {
    let mut r = [corners[b].0 - corners[a].0, corners[b].1 - corners[a].1];
    let len = (r[0] * r[0] + r[1] * r[1]).sqrt();
    r = [r[0] / len, r[1] / len];
    let delta = 2.0 / k as f64;
    let (mut ut, mut vt) = corners[a];
    for _ in 1..k {
        ut += delta * r[0];
        vt += delta * r[1];
        out.push((ut, vt));
    }
}

/// `GetTriangleEdgeMap` / `GetQuadEdgeMap` (PElementMaps.F90), 0-based:
/// the p-element edge directions (quad edges 3 and 4 run 4 -> 3 and 1 -> 4).
const TRI_EDGES: [(usize, usize); 3] = [(0, 1), (1, 2), (2, 0)];
const QUAD_EDGES: [(usize, usize); 4] = [(0, 1), (1, 2), (3, 2), (0, 3)];

/// `VTKTriangleU/V(k, :)`: corners, edge points, interior points.
fn triangle_points(k: usize) -> Vec<(f64, f64)> {
    let s3 = 3.0f64.sqrt();
    let corners = [(-1.0, 0.0), (1.0, 0.0), (0.0, s3)];
    let mut p = corners.to_vec();
    for (a, b) in TRI_EDGES {
        nodes_on_edge(&corners, a, b, k, &mut p);
    }
    if k == 3 {
        p.push((0.0, s3 / 3.0));
    } else if k >= 4 {
        let inner = triangle_points(k - 3);
        let f = (k as f64 - 3.0) / k as f64;
        let vt = s3 / k as f64;
        for &(u, v) in &inner {
            p.push((f * u, vt + f * v));
        }
    }
    p
}

/// `VTKQuadU/V(k, :)`: corners, edge points, the interior grid row by row.
fn quad_points(k: usize) -> Vec<(f64, f64)> {
    let corners = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let mut p = corners.to_vec();
    for (a, b) in QUAD_EDGES {
        nodes_on_edge(&corners, a, b, k, &mut p);
    }
    let delta = 2.0 / k as f64;
    let mut vt = -1.0;
    for _ in 1..k {
        vt += delta;
        let mut ut = -1.0;
        for _ in 1..k {
            ut += delta;
            p.push((ut, vt));
        }
    }
    p
}

/// `LagrangeDOFCount`: (points per edge, interior points) of a 2D family.
fn dof_count(family: u16, lag_n: usize) -> (usize, usize) {
    let m = lag_n - 1;
    match family {
        3 => (m, m * m.saturating_sub(1) / 2),
        _ => (m, m * m),
    }
}

/// The Lagrange points of the bulk mesh.
pub(crate) struct LagrangeMesh<'a> {
    base: &'a dyn MeshView,
    /// Per element (bulk, then boundary), its points (compact, 0-based).
    elems: Vec<Vec<u32>>,
    etypes: Vec<ElementType>,
    coords: Vec<[f64; 3]>,
}

impl MeshView for LagrangeMesh<'_> {
    fn dim(&self) -> usize {
        self.base.dim()
    }
    fn n_nodes(&self) -> usize {
        self.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.base.n_bulk()
    }
    fn n_boundary(&self) -> usize {
        self.base.n_boundary()
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let mut el = self.base.element(e);
        if let Some(nodes) = self.elems.get(e.0 as usize) {
            el.etype = self.etypes[e.0 as usize];
            el.nodes = nodes;
        }
        el
    }
}

/// The saved fields on the Lagrange points (no permutation).
pub(crate) struct LagrangeFields<'a> {
    inner: &'a dyn FieldSource,
    fields: Vec<(String, usize, Vec<f64>)>,
}

impl FieldSource for LagrangeFields<'_> {
    fn fields(&self) -> Vec<FieldRef<'_>> {
        self.fields
            .iter()
            .map(|(name, dofs, values)| FieldRef { name, kind: FieldKind::Nodal, dofs: *dofs, perm: None, values })
            .collect()
    }
    fn time(&self) -> f64 {
        self.inner.time()
    }
    fn step(&self) -> usize {
        self.inner.step()
    }
}

/// Indexes and point tables shared by the swaps of one write.
struct Layout {
    lag_n: usize,
    /// `GetLagrangeIndexes` of each bulk element, compact 0-based.
    elems: Vec<Vec<u32>>,
    n_points: usize,
    tri: Vec<(f64, f64)>,
    quad: Vec<(f64, f64)>,
}

fn layout(mesh: &dyn MeshView, pdofs: &dyn SolverDofs, lag_n: usize) -> Result<Layout, SolverError> {
    if lag_n > MAX_DEGREE {
        return Err(SolverError::fatal("HierarchicPToLagrange", "Lagrange Element Degree exceeds the supported maximum"));
    }
    let nb = mesh.n_bulk();
    let (mut nedge_max, mut nelem_max) = (0, 0);
    for e in 0..nb {
        let fam = mesh.element(ElemId(e as u32)).etype.0 / 100;
        if !(fam == 3 || fam == 4) || mesh.dim() != 2 {
            return Err(unsupported("elements other than 2D triangles and quadrilaterals"));
        }
        let (ne, ni) = dof_count(fam, lag_n);
        nedge_max = nedge_max.max(ne);
        nelem_max = nelem_max.max(ni);
    }
    let n_nodes = mesh.n_nodes();
    let edges_active = (0..nb).any(|e| pdofs.element_edges(ElemId(e as u32)).is_some());
    let n_edges = pdofs.n_edges();
    let mut global = Vec::with_capacity(nb);
    for e in 0..nb {
        let el = mesh.element(ElemId(e as u32));
        let fam = el.etype.0 / 100;
        let nc = fam as usize;
        // 1-based indexes, as the Fortran.
        let mut ix: Vec<usize> = el.nodes[..nc].iter().map(|&n| n as usize + 1).collect();
        let mut n0 = n_nodes;
        let (nedge, nelem) = dof_count(fam, lag_n);
        if edges_active {
            let edges = pdofs.element_edges(ElemId(e as u32)).ok_or_else(|| unsupported("an element without edges"))?;
            let map: &[(usize, usize)] = if fam == 3 { &TRI_EDGES } else { &QUAD_EDGES };
            for (i, &(a, b)) in map.iter().enumerate() {
                let (g, [e1, e2]) = edges[i];
                let matches = if e2 < e1 { e1 == el.nodes[b] } else { e1 == el.nodes[a] };
                let at = |j: usize| n0 + nedge_max * g as usize + j;
                if matches {
                    ix.extend((1..=nedge).map(at));
                } else {
                    ix.extend((1..=nedge).rev().map(at));
                }
            }
            n0 += n_edges * nedge_max;
        }
        ix.extend((1..=nelem).map(|j| n0 + nelem_max * e + j));
        global.push(ix);
    }
    // `DgPerm`: the indexes of the bulk elements, numbered in index order.
    let size = global.iter().flatten().copied().max().unwrap_or(0) + 1;
    let mut dg = vec![0u32; size];
    for &i in global.iter().flatten() {
        dg[i] = 1;
    }
    let mut m = 0u32;
    for d in dg.iter_mut() {
        if *d > 0 {
            m += 1;
            *d = m;
        }
    }
    // Boundary lines (`ElemDim == 1` with `BoundaryInfo`): the corners, then
    // the parent's edge points (`Find_Edge`) in the line's direction. The
    // nodes are in the p order (`AssignLocalNumber`).
    for e in nb..nb + mesh.n_boundary() {
        let el = mesh.element(ElemId(e as u32));
        if el.etype.0 / 100 != 2 {
            return Err(unsupported("boundary elements other than lines"));
        }
        let nodes = pdofs.element_nodes(ElemId(e as u32));
        let mut ix: Vec<usize> = nodes[..2].iter().map(|&n| n as usize + 1).collect();
        if edges_active {
            if let Some(parent) = el.parents[0].or(el.parents[1]) {
                let edges = pdofs.element_edges(parent).ok_or_else(|| unsupported("a parent without edges"))?;
                let (g, [e1, e2]) = *edges
                    .iter()
                    .find(|(_, en)| (en[0] == nodes[0] && en[1] == nodes[1]) || (en[0] == nodes[1] && en[1] == nodes[0]))
                    .ok_or_else(|| SolverError::fatal(CALLER, "boundary line is not an edge of its parent"))?;
                let matches = if e2 < e1 { e1 == nodes[1] } else { e1 == nodes[0] };
                let nedge = lag_n - 1;
                let at = |j: usize| n_nodes + nedge_max * g as usize + j;
                if matches {
                    ix.extend((1..=nedge).map(at));
                } else {
                    ix.extend((1..=nedge).rev().map(at));
                }
            }
        }
        global.push(ix);
    }
    let elems = global
        .into_iter()
        .map(|ix| {
            ix.into_iter()
                .map(|i| match dg.get(i) {
                    Some(&d) if d > 0 => Ok(d - 1),
                    _ => Err(SolverError::fatal(CALLER, "a boundary Lagrange point outside the bulk points")),
                })
                .collect()
        })
        .collect::<Result<Vec<Vec<u32>>, _>>()?;
    Ok(Layout { lag_n, elems, n_points: m as usize, tri: triangle_points(lag_n), quad: quad_points(lag_n) })
}

/// Where the p solution of one field lives.
enum Source<'s> {
    /// A field of a p solver: element slots and `dofs` components per slot.
    P { dofs_of: &'s dyn SolverDofs, values: &'s [f64], perm: bool },
    /// A field of a solver without p elements: `HierarchicPToLagrange`
    /// returns zeros.
    NotP { perm: Option<&'s [u32]> },
    /// The coordinates: node values, zero for the other p dofs.
    Coord { dofs_of: &'s dyn SolverDofs, axis: usize },
}

/// `p2LagrangeSwapper` of one `dofs`-component field.
fn swap(mesh: &dyn MeshView, lay: &Layout, src: &Source<'_>, dofs: usize) -> Result<Vec<f64>, SolverError> {
    let mut out = vec![0.0; dofs * lay.n_points];
    // `pSol` keeps the entries a smaller element does not overwrite.
    let mut psol: Vec<f64> = Vec::new();
    let mut slots = Vec::new();
    let mut eval = BasisEval::default();
    for e in 0..mesh.n_bulk() {
        let id = ElemId(e as u32);
        let el = mesh.element(id);
        let nc = (el.etype.0 / 100) as usize;
        let pts = if nc == 3 { &lay.tri } else { &lay.quad };
        let mut lsol = vec![0.0; dofs * pts.len()];
        let dofs_of = match src {
            Source::NotP { perm } => {
                if let Some(p) = perm {
                    if el.nodes[..nc].iter().any(|&n| p.get(n as usize).is_none_or(|&s| s == NO_DOF)) {
                        continue;
                    }
                }
                None
            }
            Source::P { dofs_of, values, perm } => {
                dofs_of.element_slots(id, &mut slots);
                if *perm && slots[..nc].contains(&NO_DOF) {
                    continue;
                }
                if psol.len() < dofs * slots.len() {
                    psol.resize(dofs * slots.len(), 0.0);
                }
                for (i, &s) in slots.iter().enumerate() {
                    for k in 0..dofs {
                        psol[dofs * i + k] = if s == NO_DOF { 0.0 } else { values[dofs * s as usize + k] };
                    }
                }
                Some(*dofs_of)
            }
            Source::Coord { dofs_of, axis } => {
                dofs_of.element_slots(id, &mut slots);
                if psol.len() < slots.len() {
                    psol.resize(slots.len(), 0.0);
                }
                let nodes = dofs_of.element_nodes(id);
                for i in 0..slots.len() {
                    psol[i] = if i < nodes.len() { mesh.coord(nodes[i])[*axis] } else { 0.0 };
                }
                Some(*dofs_of)
            }
        };
        if let Some(d) = dofs_of.filter(|d| d.is_p_element(id)) {
            if lay.lag_n == 1 {
                lsol[..dofs * nc].copy_from_slice(&psol[..dofs * nc]);
            } else {
                let basis = d.element_basis(id).map_err(|err| SolverError::fatal(CALLER, err))?;
                let nodes = d.element_nodes(id);
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&n| mesh.coord(n)).collect();
                for (t, &(u, v)) in pts.iter().enumerate() {
                    basis.eval(el.etype, &xyz, [u, v, 0.0], &mut eval).map_err(|err| SolverError::fatal(CALLER, err))?;
                    let nd = eval.basis.len();
                    if psol.len() < dofs * nd {
                        psol.resize(dofs * nd, 0.0);
                    }
                    for k in 0..dofs {
                        let mut s = 0.0;
                        for (i, b) in eval.basis.iter().enumerate() {
                            s += b * psol[dofs * i + k];
                        }
                        lsol[dofs * t + k] = s;
                    }
                }
            }
        }
        for (i, &j) in lay.elems[e].iter().enumerate().take(pts.len()) {
            out[dofs * j as usize..dofs * (j as usize + 1)].copy_from_slice(&lsol[dofs * i..dofs * (i + 1)]);
        }
    }
    Ok(out)
}

/// The p solver of each variable, and the first one (`PDefsSolver`).
fn owners(sim: &elmer_solver::Simulation) -> Vec<(elmer_solver::VarId, Option<Arc<dyn SolverDofs>>)> {
    sim.slots
        .iter()
        .enumerate()
        .filter_map(|(ix, s)| Some((s.variable?, sim.slot_dofs(ix))))
        .collect()
}

/// The Lagrange mesh and fields of one write of `fields`' requested names.
pub(crate) fn build<'a>(
    sim: &'a elmer_solver::Simulation,
    mesh: &'a dyn MeshView,
    store: &'a VariableStore,
    fields: &'a dyn FieldSource,
    opts: &elmer_vtu::VtuOptions,
) -> Result<(LagrangeMesh<'a>, LagrangeFields<'a>), SolverError> {
    let owners = owners(sim);
    let is_p = |d: &Arc<dyn SolverDofs>| (0..mesh.n_bulk()).any(|e| d.is_p_element(ElemId(e as u32)));
    let pdefs = owners
        .iter()
        .find_map(|(_, d)| d.as_ref().filter(|d| is_p(d)))
        .ok_or_else(|| unsupported("no p solver"))?;
    let lay = layout(mesh, pdefs.as_ref(), opts.lagrange_degree)?;

    let mut coords = vec![[0.0; 3]; lay.n_points];
    for axis in 0..3 {
        let c = swap(mesh, &lay, &Source::Coord { dofs_of: pdefs.as_ref(), axis }, 1)?;
        for (p, x) in coords.iter_mut().zip(c) {
            p[axis] = x;
        }
    }

    // The requested names, or their `name k` components.
    let wanted = |name: &str| {
        let name = name.trim();
        opts.scalar_fields.iter().chain(&opts.vector_fields).any(|w| {
            let w = w.trim();
            w.eq_ignore_ascii_case(name)
                || name
                    .rsplit_once(' ')
                    .is_some_and(|(b, k)| b.trim().eq_ignore_ascii_case(w) && matches!(k, "1" | "2" | "3"))
        })
    };
    let mut out = Vec::new();
    for f in fields.fields() {
        if !wanted(f.name) {
            continue;
        }
        if f.kind != FieldKind::Nodal {
            return Err(unsupported(&format!("the non-nodal field {}", f.name)));
        }
        let parent = match store.resolve(f.name) {
            Some(VarHandle::Primary(id)) | Some(VarHandle::Component { parent: id, .. }) => id,
            None => return Err(unsupported(&format!("the field {} without a variable", f.name))),
        };
        let owner = owners.iter().find(|(v, _)| *v == parent).map(|(_, d)| d.as_ref());
        let src = match owner {
            Some(Some(d)) => Source::P { dofs_of: d.as_ref(), values: f.values, perm: f.perm.is_some() },
            Some(None) => Source::NotP { perm: f.perm },
            None => return Err(unsupported(&format!("the field {} of no solver", f.name))),
        };
        let values = swap(mesh, &lay, &src, f.dofs)?;
        out.push((f.name.to_string(), f.dofs, values));
    }
    let etypes = lay
        .elems
        .iter()
        .enumerate()
        .map(|(e, ix)| {
            let fam = mesh.element(ElemId(e as u32)).etype.0 / 100;
            let n = ix.len() as u16;
            // `j = 100 * (code / 100) + n` (1000 * .. above 100 nodes).
            ElementType(if n > 100 { 1000 * fam + n } else { 100 * fam + n })
        })
        .collect();
    let lmesh = LagrangeMesh { base: mesh, elems: lay.elems, etypes, coords };
    Ok((lmesh, LagrangeFields { inner: fields, fields: out }))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Point counts of `VTKTriangleU` / `VTKQuadU` (`TriangleNodeCounts`,
    /// `QuadNodeCounts`): the full Lagrange triangle and quadrilateral.
    #[test]
    fn point_tables_have_the_lagrange_counts() {
        for k in 1..=MAX_DEGREE {
            assert_eq!(triangle_points(k).len(), (k + 1) * (k + 2) / 2, "triangle {k}");
            assert_eq!(quad_points(k).len(), (k + 1) * (k + 1), "quad {k}");
            let (ne, ni) = dof_count(3, k);
            assert_eq!(3 + 3 * ne + ni, (k + 1) * (k + 2) / 2);
            let (ne, ni) = dof_count(4, k);
            assert_eq!(4 + 4 * ne + ni, (k + 1) * (k + 1));
        }
        // Degree 3's interior point is the centroid; degree 4's interior is
        // the degree-1 corners scaled by 1/4 and lifted by sqrt(3)/4.
        assert_eq!(triangle_points(3)[9], (0.0, 3.0f64.sqrt() / 3.0));
        let t4 = triangle_points(4);
        assert_eq!(t4[12], (-0.25, 3.0f64.sqrt() / 4.0));
    }
}
