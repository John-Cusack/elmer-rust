//! `ModelMixedPoisson` (fem/src/modules/ModelMixedPoisson.F90): the mixed
//! formulation of the Poisson equation, `q = -grad u` with the flux `q` in
//! a div-conforming face element space and the potential `u` piecewise
//! constant (one bubble per element), on triangles: the Raviart-Thomas
//! space RT(0) (`Element = "n:0 e:1 -tri b:1"`) and, with `Second Kind
//! Basis`, BDM(1) (`"n:0 e:2 b:1"`). `MixedPoisson_Init0` sets the element
//! and the exported variables; `MixedPoisson_post` (its PostSolver) writes
//! the elemental potential `mixedpot` and the L2-projected DG flux
//! `mixedflux` after the solve.
//!
//! The dof layout (edge dofs, then the element's bubble, element after
//! element) is `elmer-pdofs`' `InitialPermutation`; the face basis is
//! `elmer_fe_hp::face`. Refused: 3D (face dofs), quadrilaterals, transient
//! runs, `Material Tensor`, `Convection Velocity`, a `<Variable>
//! Condition`, and `{f}` Dirichlet values other than zero.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_fe_hp::face::{face_element_info_triangle, triangle_face_orientation, FaceInfo, TRIANGLE_EDGE_MAP};
use elmer_solver::{DgIndexes, Solver, SolverCtx, SolverError, SolverRegistry};

use crate::quadrature::gauss_points_triangle_p;

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

/// A triangle's `x`, `y`, `z` and global node indices.
type Triangle = ([f64; 3], [f64; 3], [f64; 3], [u32; 3]);

/// The triangle's coordinates and global node indices.
fn triangle(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Triangle, SolverError> {
    let el = ctx.mesh.element(e);
    if el.etype.0 != 303 {
        return Err(SolverError::Unsupported(format!("ModelMixedPoisson on element type {}", el.etype.0)));
    }
    let (mut x, mut y, mut z, mut g) = ([0.0; 3], [0.0; 3], [0.0; 3], [0u32; 3]);
    for i in 0..3 {
        let c = ctx.mesh.coord(el.nodes[i]);
        (x[i], y[i], z[i]) = (c[0], c[1], c[2]);
        g[i] = el.nodes[i];
    }
    Ok((x, y, z, g))
}

/// The parent edge (1-based local index) a boundary line element lies on
/// (`PickActiveFace`), its parent, and the parent's global nodes.
fn active_face(ctx: &SolverCtx<'_>, e: ElemId) -> Option<(ElemId, usize)> {
    let el = ctx.mesh.element(e);
    let parent = el.parents[0].or(el.parents[1])?;
    let pn = ctx.mesh.element(parent).nodes;
    for (q, [a, b]) in TRIANGLE_EDGE_MAP.iter().enumerate() {
        let (na, nb) = (pn[a - 1], pn[b - 1]);
        if (el.nodes[0] == na && el.nodes[1] == nb) || (el.nodes[0] == nb && el.nodes[1] == na) {
            return Some((parent, q + 1));
        }
    }
    None
}

/// The slots of the face (edge) dofs of local edge `face` (1-based) of
/// bulk element `parent`, `per_face` of them.
fn face_slots(ctx: &SolverCtx<'_>, parent: ElemId, face: usize, per_face: usize) -> Vec<u32> {
    let mut s = Vec::new();
    ctx.element_slots(parent, &mut s);
    let first = (face - 1) * per_face;
    s[first..first + per_face].to_vec()
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
        let per_face = if second { 2 } else { 1 };
        let name = ctx.variable_name().unwrap_or("").to_ascii_lowercase();
        ctx.default_initialize()?;
        let rule = gauss_points_triangle_p(3).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;

        // Bulk assembly: LocalMatrix.
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut slots = Vec::new();
        for &e in &active {
            let (x, y, z, g) = triangle(ctx, e)?;
            ctx.element_slots(e, &mut slots);
            let nd = slots.len();
            if nd != 3 * per_face + 1 || slots.contains(&NO_DOF) {
                return Err(SolverError::fatal(CALLER, format!("element {} has {nd} dofs", e.0)));
            }
            let bp = ctx.body_params(e);
            let (material, body_force) = (bp.as_ref().and_then(|b| b.material), bp.as_ref().and_then(|b| b.body_force));
            if material.is_some_and(|m| m.get("material tensor").is_some() || m.get("convection velocity").is_some()) {
                return Err(SolverError::Unsupported("ModelMixedPoisson with Material Tensor or Convection Velocity".into()));
            }
            let matpar = element_real(ctx, material, "material parameter", e)?;
            let source = element_real(ctx, body_force, "source field", e)?;
            let mut local = LocalSystem::zeros(nd);
            let m = nd - 1; // the flux dofs; the bubble is `nd`
            for t in 0..rule.len() {
                let fi = face_element_info_triangle(&x, &y, &z, &g, rule.u[t], rule.v[t], second, true);
                let s = fi.det_f * rule.s[t];
                let a = matpar.as_ref().map_or(1.0, |h| h.at(&fi.basis));
                for p in 0..m {
                    for q in 0..m {
                        let mut d = 0.0;
                        for k in 0..2 {
                            d += fi.fbasis[q][k] * fi.fbasis[p][k];
                        }
                        local.stiff[p * nd + q] += 1.0 / a * d * s;
                    }
                }
                for p in 0..m {
                    local.stiff[p * nd + m] += 1.0 * fi.div[p] * s;
                }
                for q in 0..m {
                    local.stiff[m * nd + q] += 1.0 * fi.div[q] * s;
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
        for &e in &boundary {
            let Some(bc) = ctx.bc_params(e).cloned() else { continue };
            if bc.get(&condition).is_some() {
                return Err(SolverError::Unsupported(format!("\"{name} Condition\" in ModelMixedPoisson")));
            }
            let Some(g) = element_real(ctx, Some(&bc), "scalar field", e)? else { continue };
            let el = ctx.mesh.element(e);
            if el.etype.0 != 202 {
                continue;
            }
            let bnodes = el.nodes.to_vec();
            let Some((parent, face)) = active_face(ctx, e) else { continue };
            let (_, _, _, pg) = triangle(ctx, parent)?;
            let rev = triangle_face_orientation(&pg)[face - 1];
            let sign = if rev { -1.0 } else { 1.0 };
            let [ia, ib] = TRIANGLE_EDGE_MAP[face - 1];
            let matches = if rev { bnodes[0] == pg[ib - 1] } else { bnodes[0] == pg[ia - 1] };
            let rule = elmer_fe::quadrature::gauss_points(el.etype, None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
            let mut force = vec![0.0; per_face];
            let s3 = 3f64.sqrt();
            for t in 0..rule.len() {
                let u = rule.u[t];
                let basis = [(1.0 - u) / 2.0, (1.0 + u) / 2.0];
                let trace: Vec<f64> = if second {
                    if matches {
                        vec![sign * 0.5 * (1.0 - s3 * u), sign * 0.5 * (1.0 + s3 * u)]
                    } else {
                        vec![sign * 0.5 * (1.0 + s3 * u), sign * 0.5 * (1.0 - s3 * u)]
                    }
                } else {
                    vec![sign * 0.5]
                };
                let w = rule.s[t];
                let gv = g.at(&basis);
                for p in 0..per_face {
                    force[p] += gv * trace[p] * w;
                }
            }
            // DefaultUpdateForce on the boundary element: its dofs are the
            // parent edge's.
            let fs = face_slots(ctx, parent, face, per_face);
            ctx.element_slots(e, &mut slots);
            let mut local = LocalSystem::zeros(slots.len());
            for (p, s) in fs.iter().enumerate() {
                let Some(i) = slots.iter().position(|x| x == s) else {
                    return Err(SolverError::fatal(CALLER, format!("boundary element {} lacks its edge dofs", e.0)));
                };
                local.force[i] = force[p];
            }
            ctx.default_update(&local, e)?;
        }

        // DefaultDirichletBCs for `<name> {f} j`: the face dofs of the
        // boundary element's parent edge (only zero values are ported).
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
            let Some((parent, face)) = active_face(ctx, e) else { continue };
            for s in face_slots(ctx, parent, face, per_face) {
                fixed.push((s, 0.0));
            }
        }
        ctx.default_dirichlet_after(&fixed)?;
        ctx.default_solve()?;
        post(ctx, &active, second)
    }
}

/// `MixedPoisson_post`: the bubble value of each element into the
/// elemental `Potential Variable`, and the L2 projection of the flux onto
/// the element's nodal basis into the DG `Flux Variable`.
fn post(ctx: &mut SolverCtx<'_>, active: &[ElemId], second: bool) -> Result<(), SolverError> {
    let pot = ctx.params().string("potential variable").unwrap_or("mixedpot").to_ascii_lowercase();
    let flux = ctx.params().string("flux variable").unwrap_or("mixedflux").to_ascii_lowercase();
    let rule = gauss_points_triangle_p(3).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let dg = DgIndexes::create(ctx.mesh);
    let mut slots = Vec::new();
    let mut out: Vec<(ElemId, f64, [[f64; 3]; 3])> = Vec::with_capacity(active.len());
    {
        let sol = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no primary variable"))?.values.to_vec();
        for &e in active {
            let (x, y, z, g) = triangle(ctx, e)?;
            ctx.element_slots(e, &mut slots);
            let nd = slots.len();
            let val = sol[slots[nd - 1] as usize];
            let vals: Vec<f64> = slots[..nd - 1].iter().map(|&s| sol[s as usize]).collect();
            let n = 3;
            let mut mass = [[0.0f64; 3]; 3];
            let mut force = [[0.0f64; 3]; 3]; // [component][node]
            for t in 0..rule.len() {
                let fi: FaceInfo = face_element_info_triangle(&x, &y, &z, &g, rule.u[t], rule.v[t], second, true);
                let s = rule.s[t] * fi.det_f;
                for i in 0..n {
                    for j in 0..n {
                        mass[i][j] += s * fi.basis[i] * fi.basis[j];
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
            let mut a = vec![0.0f64; n * n];
            for i in 0..n {
                for j in 0..n {
                    a[i + j * n] = mass[i][j];
                }
            }
            elmer_fe::invert_matrix(&mut a, n)
                .ok_or_else(|| SolverError::fatal("InvertMatrix", "inversion needs successful LU-decomposition"))?;
            let mut q = [[0.0f64; 3]; 3]; // [node][component]
            for c in 0..3 {
                for i in 0..n {
                    let mut s = 0.0f64;
                    for j in 0..n {
                        s += a[i + j * n] * force[c][j];
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
