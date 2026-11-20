//! Steady-state `WhitneyAVSolver` with the Piola-transformed edge basis
//! of the first kind (`Use Piola Transform = True`, and of degree 2 with
//! `Quadratic Approximation`): `DoSolve` through the
//! solver framework, whose `elmer-pdofs` layout of the `Element`
//! definition (`n:1 e:1 -brick b:3 -quad_face b:2`, or the `n:0` one
//! without the scalar potential) numbers the nodal, edge, face and bubble
//! DOFs as Elmer's `GetElementDOFs` does.
//!
//! `LocalMatrix` is the lowest-order kernel ([`crate::local::add_bulk_point`])
//! with `ElementInfo(..., EdgeBasis, RotBasis, USolver)` taking the Piola
//! branch: `EdgeElementInfo(ApplyPiolaTransform=.TRUE.)`, whose nodal
//! `Basis` / `dBasisdx` are the p-reference element's. The rule is
//! `GaussPointsAdapt(EdgeBasis=.TRUE.)`: `ElementalNp = [1,2,3,9,11,27,18,27]`
//! on the p-reference element.
//!
//! `DefaultDirichletBCs` sets the nodal values (`P`) and, for `P {e}`, the
//! DOFs of the boundary face's edges and of the face itself
//! (`VectorElementEdgeDOFs` / `SolveLocalFaceDOFs`); only zero `{e}` values
//! are ported, for which both give zero. `ConstrainUnused` follows.
//!
//! Bricks only: the 8-node one, and the 27-node one whose geometry and
//! nodal (scalar potential) basis are its own 27-node functions.
//!
//! Refused: the second-kind basis, other element types, the boundary terms of
//! `LocalMatrixBC` / air gaps / thin lines, coils, velocities, nonlinear
//! materials and the options [`crate::steady`] refuses.

use elmer_core::{ElemId, ElementType, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_fe_hp::edge::{edge_element_info_geometry, edge_element_info_nodal};
use elmer_fe_hp::face::{BRICK_FACE_MAP, TETRA_FACE_MAP};
use elmer_solver::{SolverCtx, SolverError};

use crate::case::has_prefix;
use crate::local::{add_bulk_point, BulkCoeffs, Local};
use crate::material::{conductivity, const_real, const_vector, reluctivity};
use crate::topology::edge_map;

const CALLER: &str = "WhitneyAVSolver";

/// Whether the section asks for a Piola-transformed basis
/// (`EdgeElementStyle`).
pub fn is_piola(sp: &ParamMap) -> bool {
    ["use piola transform", "quadratic approximation", "second kind basis"]
        .iter()
        .any(|k| sp.logical(k) == Some(true))
}

/// `WhitneyAVSolver_Init0`'s `Element` for the Piola basis of the first
/// kind (`second_order`: `Quadratic Approximation`) in steady state.
pub fn init_element(static_cond: bool, second_order: bool) -> &'static str {
    match (static_cond, second_order) {
        (true, false) => "n:1 e:1 -brick b:3 -quad_face b:2",
        (false, false) => "n:0 e:1 -brick b:3 -quad_face b:2",
        (true, true) => "n:1 e:2 -brick b:6 -prism b:2 -pyramid b:3 -quad_face b:4 -tri_face b:2",
        (false, true) => "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2",
    }
}

/// `EdgeBasisDegree`: 2 with `Quadratic Approximation`.
fn degree(sp: &ParamMap) -> u8 {
    if sp.logical("quadratic approximation") == Some(true) {
        2
    } else {
        1
    }
}

fn unsupported(m: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER} (Piola): {}", m.into()))
}

fn magneto(e: crate::error::MagnetoError) -> SolverError {
    e.into()
}

const BC_TERMS: &[&str] = &[
    "magnetic field strength",
    "magnetic transfer coefficient",
    "electric current density",
    "electric transfer coefficient",
    "perfect electric conductor",
    "outflow bc",
    "thin line crossection area",
    "air gap length",
    "thin sheet thickness",
    "layer thickness",
];

fn check(ctx: &SolverCtx<'_>, sp: &ParamMap) -> Result<(), SolverError> {
    for k in crate::steady::REFUSED_SOLVER_KEYS {
        if matches!(*k, "use piola transform" | "quadratic approximation") {
            continue;
        }
        match sp.get(k) {
            None | Some(Value::Logical(false)) => {}
            Some(_) => return Err(unsupported(format!("option '{k}'"))),
        }
    }
    if sp.logical("bubbles in global system") == Some(false) {
        return Err(unsupported("Bubbles in Global System = False"));
    }
    let model = &ctx.sim.model;
    let sim_type = model.simulation.string("simulation type").unwrap_or("steady state");
    if !sim_type.trim().to_ascii_lowercase().starts_with("steady") {
        return Err(unsupported(format!("Simulation Type = {sim_type}")));
    }
    match sp.string("linear system solver").map(elmer_core::normalize_key).as_deref() {
        Some("direct") | Some("block") => return Err(unsupported("a direct or block linear solver")),
        _ => {}
    }
    let fix = sp
        .logical("fix input current density")
        .unwrap_or_else(|| model.body_forces.iter().any(|bf| has_prefix(bf, "current density")));
    if fix {
        return Err(unsupported("Fix Input Current Density (JFixPotentialSolver)"));
    }
    if !model.components.is_empty() {
        return Err(unsupported("Component sections (coils)"));
    }
    for bf in &model.body_forces {
        for k in ["angular velocity", "lorentz velocity", "local lagrange gauge", "current density multiplier"] {
            if has_prefix(bf, k) {
                return Err(unsupported(format!("body force '{k}'")));
            }
        }
    }
    for m in &model.materials {
        for k in ["laminate stack model", "magnetization", "electric permittivity"] {
            if has_prefix(m, k) {
                return Err(unsupported(format!("material '{k}'")));
            }
        }
    }
    for bc in &model.boundary_conditions {
        for k in crate::steady::REFUSED_BC_PREFIXES.iter().chain(BC_TERMS) {
            if has_prefix(bc, k) {
                return Err(unsupported(format!("boundary condition '{k}'")));
            }
        }
    }
    Ok(())
}

/// `GaussPointsAdapt(EdgeBasis=.TRUE.)` of the Piola basis: 27 points on
/// the brick for either degree (the p-reference brick is the standard
/// one). Other element types are refused.
fn rule(code: u16) -> Result<elmer_fe::quadrature::GaussRule, SolverError> {
    match code {
        808 | 827 => elmer_fe::quadrature::gauss_points(ElementType(808), Some(27)).map_err(|e| SolverError::fatal(CALLER, e.to_string())),
        c => Err(unsupported(format!("element type {c}"))),
    }
}

/// The bulk element's local system.
fn local_matrix(ctx: &SolverCtx<'_>, e: ElemId, nd: usize, degree: u8) -> Result<(Local, usize), SolverError> {
    let el = ctx.mesh.element(e);
    let code = el.etype.0;
    let n = el.nodes.len();
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let gindex: Vec<i64> = el.nodes.iter().map(|&k| k as i64 + 1).collect();
    let bp = ctx.body_params(e).ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
    let material = bp.material.ok_or_else(|| SolverError::fatal(CALLER, format!("body {} has no material", bp.body_id)))?;
    let nu = reluctivity(material, &ctx.sim.model.constants).map_err(magneto)?;
    let sigma = conductivity(material).map_err(magneto)?;
    let mut load = [0.0; 7];
    if let Some(bf) = bp.body_force {
        if let Some(j) = const_vector(bf, "current density").map_err(magneto)? {
            load[..3].copy_from_slice(&j);
        }
        if let Some(m) = const_vector(bf, "magnetization").map_err(magneto)? {
            load[3..6].copy_from_slice(&m);
        }
        if let Some(v) = const_real(bf, "electric potential").map_err(magneto)? {
            load[6] = v;
        }
    }
    let co = BulkCoeffs { nu: vec![nu; n], sigma: vec![sigma; n], load: vec![load; n] };
    let rule = rule(code)?;
    let mut loc = Local::new(nd);
    let mut np = usize::MAX;
    for t in 0..rule.len() {
        let uvw = [rule.u[t], rule.v[t], rule.w[t]];
        let (ev, basis, dbasis) = edge_info(el.etype, &xyz, &gindex, uvw, degree)?;
        // np = n * Def_Dofs(family, body, 1).
        np = match nd.checked_sub(ev.dofs) {
            Some(p) if p == 0 || p == n => p,
            _ => return Err(SolverError::fatal(CALLER, format!("{nd} element DOFs for {} edge basis functions", ev.dofs))),
        };
        add_bulk_point(&mut loc, &basis, &dbasis, &ev.edge, &ev.curl, ev.det, rule.s[t], np, &co);
    }
    Ok((loc, np))
}

/// The edge basis, nodal `Basis` and `dBasisdx` at a point.
pub(crate) type EdgeInfo = (elmer_fe_hp::edge::EdgeEval, Vec<f64>, Vec<[f64; 3]>);

/// `ElementInfo(..., EdgeBasis, RotBasis, dBasisdx, USolver)` with the Piola
/// basis: `EdgeElementInfo`, whose nodal basis is the p-reference brick's
/// for 8 nodes and `NodalBasisFunctions3D` / `NodalFirstDerivatives` for
/// the 27-node brick (`n > 8`).
pub(crate) fn edge_info(
    etype: ElementType,
    xyz: &[[f64; 3]],
    gindex: &[i64],
    uvw: [f64; 3],
    degree: u8,
) -> Result<EdgeInfo, SolverError> {
    let fatal = |m: String| SolverError::fatal(CALLER, m);
    match etype.0 {
        827 => {
            let def = elmer_fe::element_type(827).map_err(|e| fatal(e.to_string()))?;
            let mut b = vec![0.0; 27];
            let mut d = vec![[0.0; 3]; 27];
            elmer_fe::basis::nodal_basis(def, uvw[0], uvw[1], uvw[2], &mut b);
            elmer_fe::basis::nodal_first_derivatives(def, uvw[0], uvw[1], uvw[2], &mut d);
            edge_element_info_geometry(ElementType(808), xyz, &gindex[..8], uvw, degree, 3, b, d).map_err(fatal)
        }
        _ => edge_element_info_nodal(etype, xyz, gindex, uvw, degree, 3).map_err(fatal),
    }
}

/// `DefaultDirichletBCs` for `<var> {e}`: the slots of the boundary face's
/// edges and face DOFs in the parent, all zero.
fn edge_dirichlet(
    ctx: &SolverCtx<'_>,
    var: &str,
    degree: u8,
    n_parent_nodal: &dyn Fn(ElemId) -> usize,
) -> Result<Vec<(u32, f64)>, SolverError> {
    // DOFs per edge and per square face of the first-kind Piola basis.
    let (edofs, fdofs) = if degree == 2 { (2, 4) } else { (1, 2) };
    let prefix = format!("{var} {{e}}");
    let nb = ctx.mesh.n_bulk();
    let active: std::collections::HashSet<u32> = ctx.active_elements().map(|e| e.0).collect();
    let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
    let mut fixed = Vec::new();
    let mut slots = Vec::new();
    for &b in &boundary {
        let Some(bc) = ctx.bc_params(b) else { continue };
        if !has_prefix(bc, &prefix) {
            continue;
        }
        for (k, _) in bc.0.iter().filter(|(k, _)| k.starts_with(&prefix)) {
            if const_real(bc, k).map_err(magneto)? != Some(0.0) {
                return Err(unsupported(format!("a nonzero edge Dirichlet value '{k}'")));
            }
        }
        let el = ctx.mesh.element(b);
        let fam = el.etype.0 / 100;
        if fam != 3 && fam != 4 {
            return Err(unsupported(format!("'{prefix}' on element type {}", el.etype.0)));
        }
        let Some(parent) = el.parents.iter().flatten().copied().find(|p| (p.0 as usize) < nb && active.contains(&p.0)) else {
            continue;
        };
        let pel = ctx.mesh.element(parent);
        let pnodes = pel.nodes;
        let bnodes = el.nodes;
        let faces: Vec<Vec<usize>> = match pel.etype.0 {
            504 => TETRA_FACE_MAP.iter().map(|f| f.iter().map(|&k| k - 1).collect()).collect(),
            808 | 827 => BRICK_FACE_MAP.iter().map(|f| f.iter().map(|&k| k - 1).collect()).collect(),
            c => return Err(unsupported(format!("'{prefix}' on a parent of type {c}"))),
        };
        let Some(q) = faces.iter().position(|f| f.len() == bnodes.len() && f.iter().all(|&i| bnodes.contains(&pnodes[i]))) else {
            continue;
        };
        let emap = edge_map(pel.etype.0 / 100);
        let np = n_parent_nodal(parent);
        ctx.element_slots(parent, &mut slots);
        for (l, &[i, j]) in emap.iter().enumerate() {
            if bnodes.contains(&pnodes[i]) && bnodes.contains(&pnodes[j]) {
                fixed.extend_from_slice(&slots[np + edofs * l..np + edofs * (l + 1)]);
            }
        }
        // The square face's own DOFs (a tetrahedron's faces have none in
        // the first-order basis).
        if pel.etype.0 / 100 == 8 {
            let first = np + edofs * emap.len() + fdofs * q;
            fixed.extend_from_slice(&slots[first..first + fdofs]);
        } else if degree == 2 {
            return Err(unsupported("second-order '{e}' conditions on tetrahedra"));
        }
    }
    Ok(fixed.into_iter().filter(|&s| s != NO_DOF).map(|s| (s, 0.0)).collect())
}

/// `ConstrainUnused`: rows with a zero diagonal become identity rows.
fn constrain_unused(ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
    let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no linear system"))?;
    let a = &mut sys.matrix;
    let diag = a.diag.clone().ok_or_else(|| SolverError::fatal(CALLER, "matrix without a diagonal"))?;
    for r in 0..a.nrows {
        if a.values[diag[r]] == 0.0 {
            a.values[a.row_ptr[r]..a.row_ptr[r + 1]].fill(0.0);
            a.values[diag[r]] = 1.0;
            sys.rhs[r] = 0.0;
        }
    }
    Ok(())
}

/// The steady Piola `WhitneyAVSolver`.
pub fn solve(ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
    let sp = ctx.params().clone();
    if sp.logical("second kind basis") == Some(true) {
        return Err(unsupported("the second-kind basis"));
    }
    let degree = degree(&sp);
    check(ctx, &sp)?;
    if ctx.solver_dofs().is_none() {
        return Err(SolverError::fatal(CALLER, "the edge element DOFs (no DofsBackend)"));
    }
    if ctx.mesh.dim() != 3 {
        return Err(SolverError::fatal(CALLER, "Solver requires 3D mesh!"));
    }
    let var = ctx.variable_name().unwrap_or("av").to_ascii_lowercase();
    let n_nodes = ctx.mesh.n_nodes();
    let id = ctx.solver_id();
    let ix = ctx
        .sim
        .slots
        .iter()
        .position(|s| s.id == id)
        .ok_or_else(|| SolverError::fatal(CALLER, "solver slot not found"))?;
    ctx.sim.slots[ix].params.insert("Norm Permutation", Value::Integer(vec![n_nodes as i32 + 1]));

    ctx.default_initialize()?;
    let nbulk = ctx.mesh.n_bulk();
    let active: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nbulk).collect();
    let mut slots = Vec::new();
    let mut nodal = std::collections::HashMap::new();
    for &e in &active {
        ctx.element_slots(e, &mut slots);
        if slots.contains(&NO_DOF) {
            return Err(SolverError::fatal(CALLER, format!("element {} lacks DOFs", e.0)));
        }
        let nd = slots.len();
        let (loc, np) = local_matrix(ctx, e, nd, degree)?;
        nodal.insert(e.0, np);
        let ls = LocalSystem { n: nd, stiff: loc.stiff, mass: None, damp: None, force: loc.force };
        ctx.default_update(&ls, e)?;
    }
    ctx.default_finish_bulk_assembly()?;
    ctx.default_finish_boundary_assembly()?;
    let fixed = edge_dirichlet(ctx, &var, degree, &|p| nodal.get(&p.0).copied().unwrap_or(0))?;
    ctx.default_dirichlet_after(&fixed)?;
    constrain_unused(ctx)?;
    ctx.default_solve()?;
    if nodal.values().all(|&np| np == 0) {
        // WhitneyAVSolver_Init0: `Variable Output = False` without the
        // scalar potential in steady state.
        if let Some(var_id) = ctx.sim.slots[ix].variable {
            ctx.sim.vars.by_id_mut(var_id).output = false;
        }
    }
    Ok(())
}
