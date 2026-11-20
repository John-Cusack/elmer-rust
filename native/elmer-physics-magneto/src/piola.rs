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
    let st = sim_type.trim().to_ascii_lowercase();
    if !st.starts_with("steady") && st != "transient" {
        return Err(unsupported(format!("Simulation Type = {sim_type}")));
    }
    if st == "transient" {
        if sp.logical("use global mass matrix") != Some(true) {
            return Err(SolverError::fatal(CALLER, "a transient run without the global mass matrix of WhitneyAVSolver_Init0"));
        }
        if ctx.params().integer("time derivative order").is_some_and(|o| o != 1) {
            return Err(unsupported("a second-order time derivative"));
        }
        if degree(sp) == 2 {
            return Err(unsupported("a transient second-order basis"));
        }
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
    // `Apply Mortar BCs`: the framework's mortar constraints (elmer-run's
    // MortarSolver) take the `Mortar BC` boundaries.
    let mortar = sp.logical("apply mortar bcs") == Some(true);
    for bc in &model.boundary_conditions {
        for k in crate::steady::REFUSED_BC_PREFIXES.iter().chain(BC_TERMS) {
            if mortar && *k == "mortar bc" {
                continue;
            }
            if has_prefix(bc, k) {
                return Err(unsupported(format!("boundary condition '{k}'")));
            }
        }
    }
    Ok(())
}

/// `GaussPointsAdapt(EdgeBasis=.TRUE.)` of the Piola basis: 27 points on
/// the brick for either degree (the p-reference brick is the standard
/// one), 11 on the p-reference tetrahedron (first order). Other element
/// types are refused.
fn rule(code: u16, degree: u8) -> Result<elmer_fe::quadrature::GaussRule, SolverError> {
    let fatal = |e: elmer_fe::FeError| SolverError::fatal(CALLER, e.to_string());
    match code {
        808 | 827 => elmer_fe::quadrature::gauss_points(ElementType(808), Some(27)).map_err(fatal),
        504 if degree == 1 => {
            // GaussPointsTetra(11, PReferenceElement=.TRUE.).
            let mut p = elmer_fe::quadrature::gauss_points_tetra(11).map_err(fatal)?;
            for i in 0..p.len() {
                let (uq, vq, wq, sq) = (p.u[i], p.v[i], p.w[i], p.s[i]);
                p.u[i] = -1.0 + 2.0 * uq + vq + wq;
                p.v[i] = 3.0f64.sqrt() * vq + 1.0 / 3.0f64.sqrt() * wq;
                p.w[i] = 8.0f64.sqrt() / 3.0f64.sqrt() * wq;
                p.s[i] = 8.0f64.sqrt() * 2.0 * sq;
            }
            Ok(p)
        }
        c => Err(unsupported(format!("element type {c}"))),
    }
}

/// `VectorElementEdgeDOFs(List, Edge, 2, Parent, np, Name, Integral)` of the
/// lowest-order basis (DefUtils.F90): the integral over the edge of `L +
/// V . e`, `L` the scalar `Name` and `V` the components `Name 1..3` at the
/// edge's nodes, `e` the unit tangent from node `a` to node `b` (the
/// parent's edge map), negated when `b`'s global id is the smaller.
fn edge_dof_value(ctx: &SolverCtx<'_>, list: &ParamMap, name: &str, a: u32, b: u32) -> Result<f64, SolverError> {
    let nodes = [a, b];
    let zero = || vec![0.0; 2];
    let load = ctx.real_at_node_list(list, name, &nodes)?.unwrap_or_else(zero);
    let mut vload = [zero(), zero(), zero()];
    for (c, v) in vload.iter_mut().enumerate() {
        if let Some(x) = ctx.real_at_node_list(list, &format!("{name} {}", c + 1), &nodes)? {
            *v = x;
        }
    }
    let xyz = [ctx.mesh.coord(a), ctx.mesh.coord(b)];
    let mut e = [xyz[1][0] - xyz[0][0], xyz[1][1] - xyz[0][1], xyz[1][2] - xyz[0][2]];
    let len = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
    for x in e.iter_mut() {
        *x /= len;
    }
    let fatal = |m: String| SolverError::fatal("VectorElementEdgeDOFs", m);
    let def = elmer_fe::element_type(202).map_err(|x| fatal(x.to_string()))?;
    let ip = elmer_fe::quadrature::gauss_points(ElementType(202), None).map_err(|x| fatal(x.to_string()))?;
    let mut st = elmer_fe::MetricStatus::default();
    let mut be = elmer_core::BasisEval::default();
    let mut integral = 0.0;
    for p in 0..ip.len() {
        elmer_fe::element_info(def, &xyz, [ip.u[p], ip.v[p], ip.w[p]], 3, &mut st, &mut be).map_err(|x| fatal(x.to_string()))?;
        let s = ip.s[p] * be.det_j;
        let l = load[0] * be.basis[0] + load[1] * be.basis[1];
        let vl: [f64; 3] = std::array::from_fn(|c| vload[c][0] * be.basis[0] + vload[c][1] * be.basis[1]);
        integral += s * (l + (vl[0] * e[0] + vl[1] * e[1] + vl[2] * e[2]));
    }
    if b < a {
        integral = -integral;
    }
    Ok(integral)
}

/// `InitCond` for `<var> {e}` (ElmerSolver.F90:2332-2347): each active bulk
/// element with a body `Initial Condition` giving it sets its edges' DOFs
/// with [`edge_dof_value`], element by element. First-order basis only.
pub(crate) fn edge_initial_conditions(ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
    let var = ctx.variable_name().unwrap_or("av").to_ascii_lowercase();
    let key = format!("{var} {{e}}");
    if !ctx.sim.model.initial_conditions.iter().any(|ic| has_prefix(ic, &key)) {
        return Ok(());
    }
    if degree(ctx.params()) != 1 {
        return Err(unsupported("'{e}' initial conditions of the second-order basis"));
    }
    if ctx.solver_dofs().is_none() {
        return Err(SolverError::fatal(CALLER, "the edge element DOFs (no DofsBackend)"));
    }
    let nb = ctx.mesh.n_bulk();
    let active: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nb).collect();
    let mut slots = Vec::new();
    let mut set: Vec<(u32, f64)> = Vec::new();
    for &e in &active {
        let Some(ic) = ctx.body_params(e).and_then(|b| b.initial_condition).cloned() else { continue };
        if !has_prefix(&ic, &key) {
            continue;
        }
        let el = ctx.mesh.element(e);
        let (code, nodes) = (el.etype.0, el.nodes.to_vec());
        let emap = edge_map(code / 100);
        ctx.element_slots(e, &mut slots);
        let np = slots.len() - emap.len();
        for (l, &[i, j]) in emap.iter().enumerate() {
            let slot = slots[np + l];
            if slot == NO_DOF {
                continue;
            }
            set.push((slot, edge_dof_value(ctx, &ic, &key, nodes[i], nodes[j])?));
        }
    }
    let id = ctx.solver_id();
    let var_id = ctx
        .sim
        .slots
        .iter()
        .find(|s| s.id == id)
        .and_then(|s| s.variable)
        .ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    let v = ctx.sim.vars.by_id_mut(var_id);
    for (s, x) in set {
        v.values[s as usize] = x;
    }
    v.values_changed = true;
    Ok(())
}

/// The bulk element's local system.
fn local_matrix(
    ctx: &SolverCtx<'_>,
    e: ElemId,
    nd: usize,
    degree: u8,
    transient: bool,
) -> Result<(Local, usize, Option<Vec<f64>>), SolverError> {
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
    let rule = rule(code, degree)?;
    let mut loc = Local::new(nd);
    let mut mass = transient.then(|| vec![0.0; nd * nd]);
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
        if let Some(m) = mass.as_mut() {
            crate::local::add_mass_point(m, nd, &basis, &dbasis, &ev.edge, ev.det, rule.s[t], np, &co);
        }
    }
    Ok((loc, np, mass))
}

/// `ElementInfo( Face, Nodes, u, v, w, detJ, Basis, EdgeBasis = WBasis,
/// USolver )` of the lowest-order Piola basis on a boundary face, at a
/// p-reference point (the mortar projector's `AddProjectorWeakGeneric`):
/// the edge basis and the nodal basis. `cdim` is 3.
pub fn face_edge_basis(
    code: u16,
    nodes: &[[f64; 3]],
    gindex: &[i64],
    uvw: [f64; 3],
) -> Result<(Vec<[f64; 3]>, Vec<f64>), String> {
    let (ev, basis, _) = edge_element_info_nodal(ElementType(code), nodes, gindex, uvw, 1, 3)?;
    Ok((ev.edge[..ev.dofs].to_vec(), basis))
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
    edge_info_cdim(etype, xyz, gindex, uvw, degree, 3)
}

/// [`edge_info`] in a `cdim`-dimensional coordinate system
/// (`CoordinateSystemDimension()`): a 2D mesh's elements take
/// `EdgeElementInfo`'s `cdim == dim` branch.
pub(crate) fn edge_info_cdim(
    etype: ElementType,
    xyz: &[[f64; 3]],
    gindex: &[i64],
    uvw: [f64; 3],
    degree: u8,
    cdim: usize,
) -> Result<EdgeInfo, SolverError> {
    let fatal = |m: String| SolverError::fatal(CALLER, m);
    match etype.0 {
        827 => {
            let def = elmer_fe::element_type(827).map_err(|e| fatal(e.to_string()))?;
            let mut b = vec![0.0; 27];
            let mut d = vec![[0.0; 3]; 27];
            elmer_fe::basis::nodal_basis(def, uvw[0], uvw[1], uvw[2], &mut b);
            elmer_fe::basis::nodal_first_derivatives(def, uvw[0], uvw[1], uvw[2], &mut d);
            edge_element_info_geometry(ElementType(808), xyz, &gindex[..8], uvw, degree, cdim, b, d).map_err(fatal)
        }
        _ => edge_element_info_nodal(etype, xyz, gindex, uvw, degree, cdim).map_err(fatal),
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
    let mut fixed: Vec<(u32, f64)> = Vec::new();
    let mut slots = Vec::new();
    for &b in &boundary {
        let Some(bc) = ctx.bc_params(b) else { continue };
        if !has_prefix(bc, &prefix) {
            continue;
        }
        // Nonzero values: VectorElementEdgeDOFs, first-order basis only.
        let mut nonzero = false;
        for (k, _) in bc.0.iter().filter(|(k, _)| k.starts_with(&prefix)) {
            if const_real(bc, k).ok().flatten() != Some(0.0) {
                nonzero = true;
            }
        }
        if nonzero && degree != 1 {
            return Err(unsupported("a nonzero edge Dirichlet value of the second-order basis"));
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
        // The face map holds corners; a quadratic boundary element lists its
        // corners first.
        let Some(q) = faces.iter().position(|f| bnodes.len() >= f.len() && f.iter().all(|&i| bnodes[..f.len()].contains(&pnodes[i]))) else {
            continue;
        };
        let emap = edge_map(pel.etype.0 / 100);
        let np = n_parent_nodal(parent);
        ctx.element_slots(parent, &mut slots);
        for (l, &[i, j]) in emap.iter().enumerate() {
            if bnodes.contains(&pnodes[i]) && bnodes.contains(&pnodes[j]) {
                if nonzero {
                    let v = edge_dof_value(ctx, bc, &prefix, pnodes[i], pnodes[j])?;
                    fixed.push((slots[np + l], v));
                } else {
                    fixed.extend(slots[np + edofs * l..np + edofs * (l + 1)].iter().map(|&s| (s, 0.0)));
                }
            }
        }
        // The square face's own DOFs (a tetrahedron's faces have none in
        // the first-order basis).
        if pel.etype.0 / 100 == 8 {
            if nonzero {
                return Err(unsupported("a nonzero edge Dirichlet value on bricks (SolveLocalFaceDOFs)"));
            }
            let first = np + edofs * emap.len() + fdofs * q;
            fixed.extend(slots[first..first + fdofs].iter().map(|&s| (s, 0.0)));
        } else if degree == 2 {
            return Err(unsupported("second-order '{e}' conditions on tetrahedra"));
        }
    }
    Ok(fixed.into_iter().filter(|&(s, _)| s != NO_DOF).collect())
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

/// The Piola `WhitneyAVSolver`, steady or (first order) transient.
pub fn solve(ctx: &mut SolverCtx<'_>, transient: bool) -> Result<(), SolverError> {
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
        let (loc, np, mass) = local_matrix(ctx, e, nd, degree, transient)?;
        nodal.insert(e.0, np);
        // DefaultUpdateMass: the global mass matrix (`Use Global Mass
        // Matrix`), whose time terms DefaultFinishAssembly adds.
        if let Some(m) = mass {
            ctx.default_update_mass(&m, e)?;
        }
        let ls = LocalSystem { n: nd, stiff: loc.stiff, mass: None, damp: None, force: loc.force };
        ctx.default_update(&ls, e)?;
    }
    ctx.default_finish_bulk_assembly()?;
    ctx.default_finish_boundary_assembly()?;
    let fixed = edge_dirichlet(ctx, &var, degree, &|p| nodal.get(&p.0).copied().unwrap_or(0))?;
    ctx.default_dirichlet_after(&fixed)?;
    constrain_unused(ctx)?;
    ctx.default_solve()?;
    if !transient && sp.get("variable output").is_none() && nodal.values().all(|&np| np == 0) {
        // WhitneyAVSolver_Init0: `Variable Output = False` without the
        // scalar potential in steady state.
        if let Some(var_id) = ctx.sim.slots[ix].variable {
            ctx.sim.vars.by_id_mut(var_id).output = false;
        }
    }
    Ok(())
}
