//! `ModelMixedPoisson` (fem/src/modules/ModelMixedPoisson.F90): the mixed
//! formulation of the Poisson equation, `q = -grad u` with the flux `q` in
//! a div-conforming face element space and the potential `u` piecewise
//! constant (one bubble per element), on triangles, tetrahedra and
//! bricks (the 48-function space of `-brick b:25 -quad_face b:4`): the
//! Raviart-Thomas space RT(0) (`Element = "n:0 e:1 -tri b:1"`, `"n:0
//! -tetra b:1 -tri_face b:1"`) and, with `Second Kind Basis`, BDM(1)
//! (`"n:0 e:2 b:1"`, `-tri_face b:3`). `MixedPoisson_Init0` sets the element
//! and the exported variables; `MixedPoisson_post` (its PostSolver) writes
//! the elemental potential `mixedpot` and the L2-projected DG flux
//! `mixedflux` after the solve.
//!
//! The dof layout (edge dofs, then the element's bubble, element after
//! element) is `elmer-pdofs`' `InitialPermutation`; the face basis is
//! `elmer_fe_hp::face`. Refused: quadrilaterals, wedges and pyramids,
//! transient runs, a vector `Material Tensor`, `Convection Velocity`, a `<Variable>
//! Condition`, and `{f}` Dirichlet values other than zero.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_core::{BasisEval, ElementBasis};
use elmer_fe_hp::face::{
    brick_face_basis_ordering, brick_face_orientation, face_element_info_brick, face_element_info_tetra,
    face_element_info_triangle, tetra_face_basis_ordering, tetra_face_orientation, triangle_face_orientation,
    triangle_nodal_p_basis, FaceInfo, BRICK_FACE_MAP, TETRA_FACE_MAP, TRIANGLE_EDGE_MAP,
};
use elmer_solver::{DgIndexes, Solver, SolverCtx, SolverError, SolverRegistry};

use elmer_fe::quadrature::GaussRule;

use crate::quadrature::{gauss_points_tetra_p, gauss_points_triangle_p};

/// `Procedure = "ModelMixedPoisson" "MixedPoisson"`.
pub const PROCEDURE: (&str, &str) = ("ModelMixedPoisson", "MixedPoisson");

const CALLER: &str = "MixedPoisson";

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| Ok(Box::new(MixedPoisson)));
}

fn add_new(p: &mut ParamMap, key: &str, v: Value) {
    if p.get(key).is_none() {
        p.insert(key, v);
    }
}

/// `NextFreeKeyword('Exported Variable', Params)`.
fn next_exported(p: &ParamMap) -> String {
    let mut k = 1;
    while p.get(&format!("exported variable {k}")).is_some() {
        k += 1;
    }
    format!("Exported Variable {k}")
}

/// `MixedPoisson_Init0`: the keywords it adds to the solver section
/// (`ListAddNew*`: a keyword already there wins). `csys` is the
/// Simulation's `Coordinate System`. Of the two `Exported Variable`s it adds
/// per field, the plain-named second one names a variable that already
/// exists and adds nothing, so it is left out.
pub fn solver_defaults(p: &mut ParamMap, csys: Option<&str>) {
    let second = p.logical("second kind basis").unwrap_or(false);
    match csys.map(|s| s.trim().to_ascii_lowercase()).as_deref() {
        Some("cartesian 2d") => {
            let e = if second { "n:0 e:2 b:1" } else { "n:0 e:1 -tri b:1 -quad b:3" };
            add_new(p, "Element", Value::Str(e.into()));
        }
        Some("cartesian 3d") => {
            let e = if second {
                "n:0 -tetra b:1 -brick b:25 -quad_face b:4 -tri_face b:3"
            } else {
                "n:0 -tetra b:1 -brick b:25 -quad_face b:4 -tri_face b:1"
            };
            add_new(p, "Element", Value::Str(e.into()));
        }
        _ => {}
    }
    add_new(p, "Bubbles in Global System", Value::Logical(true));
    add_new(p, "Potential Variable", Value::Str("mixedpot".into()));
    let pot = p.string("potential variable").unwrap_or("mixedpot").to_string();
    let k = next_exported(p);
    p.insert(&k, Value::Str(format!("-elem {pot}")));
    add_new(p, "Flux Variable", Value::Str("mixedflux".into()));
    let flux = p.string("flux variable").unwrap_or("mixedflux").to_string();
    let k = next_exported(p);
    p.insert(&k, Value::Str(format!("-dofs 3 -dg {flux}")));
}

struct MixedPoisson;

/// `ListGetElementReal(h, Basis, Element)`: a constant as it is, else the
/// nodal values interpolated with `Basis`. `None` when the section lacks
/// the keyword.
fn element_real(ctx: &SolverCtx<'_>, section: Option<&ParamMap>, key: &str, e: ElemId) -> Result<Option<Elemental>, SolverError> {
    let Some(sec) = section else { return Ok(None) };
    match sec.get(key) {
        None => Ok(None),
        Some(Value::Real { data, .. }) if data.len() == 1 => Ok(Some(Elemental::Constant(data[0]))),
        Some(_) => Ok(ctx.real_at_nodes(sec, key, e)?.map(Elemental::Nodal)),
    }
}

enum Elemental {
    Constant(f64),
    Nodal(Vec<f64>),
}

impl Elemental {
    fn at(&self, basis: &[f64]) -> f64 {
        match self {
            Elemental::Constant(c) => *c,
            Elemental::Nodal(v) => {
                let mut s = 0.0;
                for i in 0..v.len() {
                    s += basis[i] * v[i];
                }
                s
            }
        }
    }
}

/// A bulk element of the solver: a triangle (303), a tetrahedron (504)
/// or a brick (808).
struct Cell {
    code: u16,
    x: [f64; 8],
    y: [f64; 8],
    z: [f64; 8],
    g: [u32; 8],
}

impl Cell {
    fn of(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Self, SolverError> {
        let el = ctx.mesh.element(e);
        let code = el.etype.0;
        if !matches!(code, 303 | 504 | 808) {
            return Err(SolverError::Unsupported(format!("ModelMixedPoisson on element type {code}")));
        }
        let mut c = Cell { code, x: [0.0; 8], y: [0.0; 8], z: [0.0; 8], g: [0; 8] };
        for i in 0..c.n() {
            let p = ctx.mesh.coord(el.nodes[i]);
            (c.x[i], c.y[i], c.z[i]) = (p[0], p[1], p[2]);
            c.g[i] = el.nodes[i];
        }
        Ok(c)
    }

    fn brick(&self) -> bool {
        self.code == 808
    }

    fn n(&self) -> usize {
        match self.code {
            303 => 3,
            504 => 4,
            _ => 8,
        }
    }

    fn faces(&self) -> usize {
        match self.code {
            303 => 3,
            504 => 4,
            _ => 6,
        }
    }

    fn dim(&self) -> usize {
        if self.code == 303 { 2 } else { 3 }
    }

    /// Dofs per face: RT(0) one; BDM(1) two on an edge, three on a face;
    /// four on a brick's face.
    fn per_face(&self, second: bool) -> usize {
        match (self.code, second) {
            (808, _) => 4,
            (_, false) => 1,
            (303, true) => 2,
            _ => 3,
        }
    }

    /// The flux basis functions: the face ones, and a brick's 24 interior
    /// ones (elementwise bubbles before the potential's).
    fn flux_dofs(&self, second: bool) -> usize {
        self.faces() * self.per_face(second) + if self.brick() { 24 } else { 0 }
    }

    /// The Gauss rule `LocalMatrix` and `GetFlux` pick for the element.
    fn rule(&self) -> Result<GaussRule, SolverError> {
        let r = match self.code {
            303 => gauss_points_triangle_p(3),
            504 => gauss_points_tetra_p(4),
            _ => elmer_fe::quadrature::gauss_points_brick(64).map_err(|e| crate::Error::Unsupported(e.to_string())),
        };
        r.map_err(|m| SolverError::fatal(CALLER, m.to_string()))
    }

    fn info(&self, rule: &GaussRule, t: usize, second: bool) -> FaceInfo {
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
    fn face_nodes(&self, face: usize) -> Vec<u32> {
        match self.code {
            303 => TRIANGLE_EDGE_MAP[face - 1].iter().map(|&k| self.g[k - 1]).collect(),
            504 => TETRA_FACE_MAP[face - 1].iter().map(|&k| self.g[k - 1]).collect(),
            _ => BRICK_FACE_MAP[face - 1].iter().map(|&k| self.g[k - 1]).collect(),
        }
    }

    /// `PickActiveFace`: the 1-based local face (an edge of a triangle) the
    /// boundary element with nodes `b` lies on.
    fn face_of(&self, b: &[u32]) -> Option<usize> {
        (1..=self.faces()).find(|&q| {
            let f = self.face_nodes(q);
            f.len() == b.len() && f.iter().all(|i| b.contains(i))
        })
    }

    /// `FaceElementOrientation(Parent, ReverseSign, ActiveFaceId)`.
    fn reverse(&self, face: usize) -> bool {
        match self.code {
            303 => triangle_face_orientation(&[self.g[0], self.g[1], self.g[2]])[face - 1],
            504 => tetra_face_orientation(&[self.g[0], self.g[1], self.g[2], self.g[3]])[face - 1],
            _ => brick_face_orientation(&self.g)[face - 1],
        }
    }
}

fn first4(a: &[f64; 8]) -> [f64; 4] {
    [a[0], a[1], a[2], a[3]]
}

/// A constant rank-2 `Material Tensor` (`mu_tensor(i,j)`), `None` when the
/// material has none; other shapes are refused.
fn material_tensor(material: Option<&ParamMap>) -> Result<Option<[[f64; 3]; 3]>, SolverError> {
    let Some(v) = material.and_then(|m| m.get("material tensor")) else { return Ok(None) };
    match v {
        Value::Real { rows: 3, cols: 3, data } => {
            let mut t = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    t[i][j] = data[i * 3 + j];
                }
            }
            Ok(Some(t))
        }
        _ => Err(SolverError::Unsupported("ModelMixedPoisson with a Material Tensor that is not a constant 3x3 matrix".into())),
    }
}

/// The parent of a boundary element and the face it lies on.
fn active_face(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Option<(ElemId, Cell, usize)>, SolverError> {
    let el = ctx.mesh.element(e);
    let Some(parent) = el.parents[0].or(el.parents[1]) else { return Ok(None) };
    let cell = Cell::of(ctx, parent)?;
    Ok(cell.face_of(el.nodes).map(|q| (parent, cell, q)))
}

/// `np`: the nodal dofs of an element with `nd` dofs (`n` with an `n:1`
/// element definition, else none), before its face dofs and its bubble.
fn nodal_dofs(cell: &Cell, nd: usize, second: bool) -> Result<usize, SolverError> {
    let flux = cell.flux_dofs(second);
    match nd.checked_sub(flux + 1) {
        Some(np) if np == 0 || np == cell.n() => Ok(np),
        _ => Err(SolverError::fatal(CALLER, format!("an element with {nd} dofs ({flux} face dofs expected)"))),
    }
}

/// The slots of the face dofs of local face `face` (1-based) of bulk
/// element `parent`, `per_face` of them.
fn face_slots(ctx: &SolverCtx<'_>, parent: ElemId, cell: &Cell, face: usize, second: bool) -> Result<Vec<u32>, SolverError> {
    let mut s = Vec::new();
    ctx.element_slots(parent, &mut s);
    let np = nodal_dofs(cell, s.len(), second)?;
    let per_face = cell.per_face(second);
    let first = np + (face - 1) * per_face;
    Ok(s[first..first + per_face].to_vec())
}

impl Solver for MixedPoisson {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(SolverError::Unsupported("transient ModelMixedPoisson".into()));
        }
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::fatal(CALLER, "the face element dofs (no DofsBackend)"));
        }
        let second = ctx.params().logical("second kind basis").unwrap_or(false);
        let name = ctx.variable_name().unwrap_or("").to_ascii_lowercase();
        ctx.default_initialize()?;

        // Bulk assembly: LocalMatrix.
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut slots = Vec::new();
        for &e in &active {
            let cell = Cell::of(ctx, e)?;
            let dim = cell.dim();
            let rule = cell.rule()?;
            ctx.element_slots(e, &mut slots);
            let nd = slots.len();
            let np = nodal_dofs(&cell, nd, second)?;
            if slots.contains(&NO_DOF) {
                return Err(SolverError::fatal(CALLER, format!("element {} lacks dofs", e.0)));
            }
            let bp = ctx.body_params(e);
            let (material, body_force) = (bp.as_ref().and_then(|b| b.material), bp.as_ref().and_then(|b| b.body_force));
            if material.is_some_and(|m| m.get("convection velocity").is_some()) {
                return Err(SolverError::Unsupported("ModelMixedPoisson with Convection Velocity".into()));
            }
            let mu = material_tensor(material)?;
            let matpar = element_real(ctx, material, "material parameter", e)?;
            let source = element_real(ctx, body_force, "source field", e)?;
            let mut local = LocalSystem::zeros(nd);
            let m = nd - 1; // the bubble (`nd` in the Fortran)
            let nf = nd - np - 1; // the flux dofs, after the nodal ones
            let n = cell.n();
            for t in 0..rule.len() {
                let fi = cell.info(&rule, t, second);
                let s = fi.det_f * rule.s[t];
                // The Galerkin projection of the potential onto the nodal dofs.
                if np > 0 {
                    for p in 0..n {
                        for q in 0..n {
                            local.stiff[p * nd + q] += fi.basis[p] * fi.basis[q] * s;
                        }
                        local.stiff[p * nd + m] -= fi.basis[p] * 1.0 * s;
                    }
                }
                if let Some(mu) = &mu {
                    // tmp = MATMUL(mu_tensor(1:dim,1:dim), FaceBasis(q,1:dim)).
                    for p in 0..nf {
                        for q in 0..nf {
                            let mut tmp = [0.0f64; 3];
                            for i in 0..dim {
                                let mut t = 0.0;
                                for j in 0..dim {
                                    t += mu[i][j] * fi.fbasis[q][j];
                                }
                                tmp[i] = t;
                            }
                            let mut d = 0.0;
                            for k in 0..dim {
                                d += tmp[k] * fi.fbasis[p][k];
                            }
                            local.stiff[(np + p) * nd + np + q] += d * s;
                        }
                    }
                } else {
                    let a = matpar.as_ref().map_or(1.0, |h| h.at(&fi.basis));
                    for p in 0..nf {
                        for q in 0..nf {
                            let mut d = 0.0;
                            for k in 0..dim {
                                d += fi.fbasis[q][k] * fi.fbasis[p][k];
                            }
                            local.stiff[(np + p) * nd + np + q] += 1.0 / a * d * s;
                        }
                    }
                }
                for p in 0..nf {
                    local.stiff[(np + p) * nd + m] += 1.0 * fi.div[p] * s;
                }
                for q in 0..nf {
                    local.stiff[m * nd + np + q] += 1.0 * fi.div[q] * s;
                }
                if let Some(h) = &source {
                    let f = h.at(&fi.basis);
                    local.force[m] -= f * 1.0 * s;
                }
            }
            ctx.default_update(&local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;

        // Boundary assembly: LocalMatrixBC, the `Scalar Field` of the
        // potential through the trace of the flux basis.
        let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
        let condition = format!("{name} condition");
        let s3 = 3f64.sqrt();
        for &e in &boundary {
            let Some(bc) = ctx.bc_params(e).cloned() else { continue };
            if bc.get(&condition).is_some() {
                return Err(SolverError::Unsupported(format!("\"{name} Condition\" in ModelMixedPoisson")));
            }
            let Some(g) = element_real(ctx, Some(&bc), "scalar field", e)? else { continue };
            let el = ctx.mesh.element(e);
            let mut bnodes = el.nodes.to_vec();
            let Some((parent, cell, face)) = active_face(ctx, e)? else { continue };
            let per_face = cell.per_face(second);
            let rev = cell.reverse(face);
            let sign = if rev { -1.0 } else { 1.0 };
            let mut force = vec![0.0; per_face];
            match (cell.code, el.etype.0) {
                (303, 202) => {
                    let [ia, ib] = TRIANGLE_EDGE_MAP[face - 1];
                    let matches = if rev { bnodes[0] == cell.g[ib - 1] } else { bnodes[0] == cell.g[ia - 1] };
                    let rule = elmer_fe::quadrature::gauss_points(el.etype, None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
                    for t in 0..rule.len() {
                        let u = rule.u[t];
                        let basis = [(1.0 - u) / 2.0, (1.0 + u) / 2.0];
                        let trace: Vec<f64> = match (second, matches) {
                            (true, true) => vec![sign * 0.5 * (1.0 - s3 * u), sign * 0.5 * (1.0 + s3 * u)],
                            (true, false) => vec![sign * 0.5 * (1.0 + s3 * u), sign * 0.5 * (1.0 - s3 * u)],
                            (false, _) => vec![sign * 0.5],
                        };
                        let w = rule.s[t];
                        let gv = g.at(&basis);
                        for p in 0..per_face {
                            force[p] += gv * trace[p] * w;
                        }
                    }
                }
                (504, 303) => {
                    // BDM: the parent's ordering of the face dofs, with the
                    // element's nodes taken in the parent face's order.
                    let g4 = [cell.g[0], cell.g[1], cell.g[2], cell.g[3]];
                    let fmap = if second { Some(tetra_face_basis_ordering(&g4)[face - 1]) } else { None };
                    if second {
                        bnodes = cell.face_nodes(face);
                    }
                    let g = match g {
                        Elemental::Nodal(_) if bnodes != el.nodes => element_real_at(ctx, &bc, "scalar field", e, &bnodes)?,
                        other => other,
                    };
                    let rule = gauss_points_triangle_p(3).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
                    for t in 0..rule.len() {
                        let (u, v) = (rule.u[t], rule.v[t]);
                        let basis: Vec<f64> = (1..=3).map(|k| triangle_nodal_p_basis(k, u, v)).collect();
                        let trace: Vec<f64> = match fmap {
                            Some(map) => {
                                let wt = [
                                    -(-8.0 - 12.0 * u + 4.0 * s3 * v) / 12.0,
                                    -(u + (-8.0 + 4.0 * s3 * v) / 12.0),
                                    -(4.0 - 8.0 * s3 * v) / 12.0,
                                ];
                                map.iter().map(|&k| sign * wt[k - 1]).collect()
                            }
                            None => vec![sign / s3],
                        };
                        let w = rule.s[t];
                        let gv = g.at(&basis);
                        for p in 0..per_face {
                            force[p] += gv * trace[p] * w;
                        }
                    }
                }
                (808, 404) => {
                    // The parent's ordering of the face dofs; the element's
                    // nodes in the parent face's order; the trace is the
                    // element's own nodal basis.
                    let fmap = brick_face_basis_ordering(&cell.g)[face - 1];
                    bnodes = cell.face_nodes(face);
                    let g = match g {
                        Elemental::Nodal(_) if bnodes != el.nodes => element_real_at(ctx, &bc, "scalar field", e, &bnodes)?,
                        other => other,
                    };
                    let xyz: Vec<[f64; 3]> = bnodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                    let lagrange = elmer_fe::LagrangeBasis::new(3);
                    let rule = elmer_fe::quadrature::gauss_points(el.etype, None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
                    let mut ev = BasisEval::default();
                    for t in 0..rule.len() {
                        lagrange
                            .eval(el.etype, &xyz, [rule.u[t], rule.v[t], rule.w[t]], &mut ev)
                            .map_err(|m| SolverError::fatal("ElementInfo", m))?;
                        let w = rule.s[t];
                        let gv = g.at(&ev.basis);
                        for p in 0..per_face {
                            force[p] += gv * (sign * ev.basis[fmap[p] - 1]) * w;
                        }
                    }
                }
                _ => continue,
            }
            // DefaultUpdateForce on the boundary element: its dofs are the
            // parent face's.
            let fs = face_slots(ctx, parent, &cell, face, second)?;
            ctx.element_slots(e, &mut slots);
            let mut local = LocalSystem::zeros(slots.len());
            for (p, s) in fs.iter().enumerate() {
                let Some(i) = slots.iter().position(|x| x == s) else {
                    return Err(SolverError::fatal(CALLER, format!("boundary element {} lacks its face dofs", e.0)));
                };
                local.force[i] = force[p];
            }
            ctx.default_update(&local, e)?;
        }

        // DefaultDirichletBCs for `<name> {f} j`: the face dofs of the
        // boundary element's parent face (only zero values are ported).
        let prefix = format!("{name} {{f}}");
        let mut fixed: Vec<(u32, f64)> = Vec::new();
        for &e in &boundary {
            let Some(bc) = ctx.bc_params(e) else { continue };
            let mut given = false;
            for (k, v) in bc.0.iter() {
                if k.starts_with(&prefix) {
                    given = true;
                    match v {
                        Value::Real { data, .. } if data.iter().all(|&x| x == 0.0) => {}
                        _ => return Err(SolverError::Unsupported(format!("a nonzero or non-constant \"{k}\" in ModelMixedPoisson"))),
                    }
                }
            }
            if !given {
                continue;
            }
            let Some((parent, cell, face)) = active_face(ctx, e)? else { continue };
            for s in face_slots(ctx, parent, &cell, face, second)? {
                fixed.push((s, 0.0));
            }
        }
        ctx.default_dirichlet_after(&fixed)?;
        ctx.default_solve()?;
        post(ctx, &active, second)
    }
}

/// [`element_real`] at nodes `nodes` (the element's, reordered).
fn element_real_at(ctx: &SolverCtx<'_>, sec: &ParamMap, key: &str, e: ElemId, nodes: &[u32]) -> Result<Elemental, SolverError> {
    let vals = ctx.real_at_nodes(sec, key, e)?.unwrap_or_default();
    let orig = ctx.mesh.element(e).nodes;
    let re: Vec<f64> = nodes.iter().map(|n| orig.iter().position(|m| m == n).map_or(0.0, |i| vals[i])).collect();
    Ok(Elemental::Nodal(re))
}

/// `MixedPoisson_post`: the bubble value of each element into the
/// elemental `Potential Variable`, and the L2 projection of the flux onto
/// the element's nodal basis into the DG `Flux Variable`.
fn post(ctx: &mut SolverCtx<'_>, active: &[ElemId], second: bool) -> Result<(), SolverError> {
    let pot = ctx.params().string("potential variable").unwrap_or("mixedpot").to_ascii_lowercase();
    let flux = ctx.params().string("flux variable").unwrap_or("mixedflux").to_ascii_lowercase();
    let dg = DgIndexes::create(ctx.mesh);
    let mut slots = Vec::new();
    let mut out: Vec<(ElemId, f64, Vec<[f64; 3]>)> = Vec::with_capacity(active.len());
    {
        let sol = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no primary variable"))?.values.to_vec();
        for &e in active {
            let cell = Cell::of(ctx, e)?;
            let rule = cell.rule()?;
            ctx.element_slots(e, &mut slots);
            let nd = slots.len();
            let np = nodal_dofs(&cell, nd, second)?;
            let val = sol[slots[nd - 1] as usize];
            let vals: Vec<f64> = slots[np..nd - 1].iter().map(|&s| sol[s as usize]).collect();
            let n = cell.n();
            let mut mass = vec![0.0f64; n * n]; // column-major, as InvertMatrix takes it
            let mut force = [[0.0f64; 8]; 3]; // [component][node]
            for t in 0..rule.len() {
                let fi = cell.info(&rule, t, second);
                let s = rule.s[t] * fi.det_f;
                for i in 0..n {
                    for j in 0..n {
                        mass[i + j * n] += s * fi.basis[i] * fi.basis[j];
                    }
                    for c in 0..3 {
                        let mut d = 0.0;
                        for k in 0..vals.len() {
                            d += fi.fbasis[k][c] * vals[k];
                        }
                        force[c][i] += s * d * fi.basis[i];
                    }
                }
            }
            // InvertMatrix(MASS, n); Flux = MATMUL(MASS, FORCE), which the
            // optimised build that made the reference results inlines as
            // `Flux(i) = Flux(i) + MASS(i,j) * FORCE(j)` in `j` order.
            elmer_fe::invert_matrix(&mut mass, n)
                .ok_or_else(|| SolverError::fatal("InvertMatrix", "inversion needs successful LU-decomposition"))?;
            let mut q = vec![[0.0f64; 3]; n]; // [node][component]
            for c in 0..3 {
                for i in 0..n {
                    let mut s = 0.0f64;
                    for j in 0..n {
                        s += mass[i + j * n] * force[c][j];
                    }
                    q[i][c] = s;
                }
            }
            out.push((e, val, q));
        }
    }
    {
        let pv = ctx
            .variable_mut(&pot)
            .ok_or_else(|| SolverError::fatal(CALLER, format!("no variable {pot}")))?;
        let perm = pv.perm.map(<[u32]>::to_vec);
        for (e, val, _) in &out {
            let k = perm.as_ref().map_or(e.0, |p| p[e.0 as usize]);
            if k != NO_DOF {
                pv.values[k as usize] = *val;
            }
        }
    }
    let fv = ctx
        .variable_mut(&flux)
        .ok_or_else(|| SolverError::fatal(CALLER, format!("no variable {flux}")))?;
    let perm = fv.perm.map(<[u32]>::to_vec);
    for (e, _, q) in &out {
        for (i, &d) in dg.index[e.0 as usize].iter().enumerate() {
            let k = perm.as_ref().map_or(d, |p| p[d as usize]);
            if k == NO_DOF {
                continue;
            }
            for c in 0..3 {
                fv.values[3 * k as usize + c] = q[i][c];
            }
        }
    }
    Ok(())
}
