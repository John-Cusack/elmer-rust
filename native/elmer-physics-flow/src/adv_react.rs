//! `AdvectionReactionSolver` (fem/src/modules/AdvectionReaction.F90): the
//! steady advection-reaction equation `u . grad c + gamma c = f` with
//! discontinuous Galerkin, upwind fluxes on the interior edges (faces in
//! 3D) and a weak inflow condition on the boundary.
//!
//! - Bulk (`LocalMatrix`): `(gamma + div v_mesh) c w - c (u . grad w)`, the
//!   source `f w`.
//! - Every interior edge (`LocalJumps`): `u.n {c} [w] + |u.n|/2 [c] [w]`,
//!   the normal out of the left parent, the parents' bases at the edge's
//!   integration points (`FindParentUVW`), glued at both parents' DG
//!   indexes.
//! - Every boundary element (`LocalMatrixBoundary`): on inflow (the
//!   variable given on the boundary and the mean `u.n < 0`) the load
//!   `-u.n c_in w`, else the outflow term `u.n c w`, on the parent's basis.
//!
//! The velocity is the material's `Convection Velocity i` (`Convection =
//! constant` in the equation) or `none`. Ported for 2D meshes on the DG
//! layouts of `elmer-pdofs`. `Calculate Nodal Average` adds the nodal
//! `<var>_nodal` and, after each solve, sets it to the mean of the DG
//! values at each node (`AdvectionReactionSolver_post`). Refused:
//! `computed` convection, a mesh velocity, `Limit Solution`, transient
//! runs, 3D.

use std::sync::Arc;

use elmer_core::{BasisEval, ElemId, ElementType, FieldKind, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::expr::eval_at_nodes_with;
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "AdvectionReaction";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(AdvectionReaction))
}

pub struct AdvectionReaction;

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("AdvectionReactionSolver: {what}"))
}

/// `GetReal(section, key, Found, Element)` at `nodes`; zeros without it.
fn at_nodes(ctx: &SolverCtx<'_>, section: Option<&ParamMap>, key: &str, nodes: &[u32]) -> Result<Vec<f64>, SolverError> {
    let Some(v) = section.and_then(|s| s.get(key)) else {
        return Ok(vec![0.0; nodes.len()]);
    };
    let b = &ctx.sim.backends;
    eval_at_nodes_with(v, nodes, &ctx.fields(), b.matc.as_deref(), &b.udfs)
}

/// `GetLocalALEVelocity` for `Convection = constant | none`: the velocity
/// at `nodes`, from the material of `parent`.
fn velocity(ctx: &SolverCtx<'_>, parent: ElemId, nodes: &[u32]) -> Result<Vec<[f64; 3]>, SolverError> {
    let bp = ctx.body_params(parent).ok_or_else(|| SolverError::fatal(CALLER, "an element without a body"))?;
    let flag = bp
        .equation
        .and_then(|q| q.string("convection"))
        .map(elmer_core::normalize_key)
        .ok_or_else(|| SolverError::fatal(CALLER, "No string for keyword >Convection< found in Equation"))?;
    let mut v = vec![[0.0; 3]; nodes.len()];
    match flag.as_str() {
        "constant" => {
            for i in 0..3 {
                for (o, x) in v.iter_mut().zip(at_nodes(ctx, bp.material, &format!("convection velocity {}", i + 1), nodes)?) {
                    o[i] = x;
                }
            }
        }
        "none" => {}
        "computed" => return Err(unsupported("Convection = computed")),
        other => return Err(SolverError::fatal(CALLER, format!("Convection flag >{other}< not recognised"))),
    }
    Ok(v)
}

/// `AdvectionReactionSolver_post`: `<var>_nodal` (or `Nodal Result`) as
/// the mean of the DG values at each node, over the active elements.
fn nodal_average(ctx: &mut SolverCtx<'_>, var: &str) -> Result<(), SolverError> {
    let target = [format!("{var}_nodal"), "nodal result".to_string()]
        .into_iter()
        .find(|n| ctx.variable(n).is_some());
    let Some(target) = target else { return Ok(()) };
    let (sperm, svals) = {
        let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
        (v.perm.map(|p| p.to_vec()), v.values.to_vec())
    };
    let ix = ctx.solver_dofs().and_then(|d| d.dg_layout()).cloned().ok_or_else(|| SolverError::fatal(CALLER, "no DG layout"))?;
    let tperm = ctx.variable(&target).and_then(|v| v.perm.map(|p| p.to_vec())).ok_or_else(|| SolverError::fatal(CALLER, "nodal average without a Perm"))?;
    let n = ctx.variable(&target).map_or(0, |v| v.values.len());
    let mut sum = vec![0.0; n];
    let mut cnt = vec![0u32; n];
    let elems: Vec<ElemId> = ctx.active_elements().collect();
    for e in elems {
        for (i, &node) in ctx.mesh.element(e).nodes.iter().enumerate() {
            let d = ix.index[e.0 as usize][i];
            let j = match &sperm {
                Some(p) if d != NO_DOF => p[d as usize],
                _ => NO_DOF,
            };
            if j == NO_DOF {
                continue;
            }
            let k = tperm[node as usize];
            if k == NO_DOF {
                continue;
            }
            sum[k as usize] += svals[j as usize];
            cnt[k as usize] += 1;
        }
    }
    if let Some(v) = ctx.variable_mut(&target) {
        for (i, x) in v.values.iter_mut().enumerate() {
            *x = if cnt[i] > 0 { sum[i] / f64::from(cnt[i]) } else { 0.0 };
        }
    }
    Ok(())
}

/// `FindParentUVW`: the parent's reference coordinates at a point of
/// `face` (whose basis values there are `fb`).
fn parent_uvw(face: &[u32], parent: &[u32], ptype: ElementType, fb: &[f64]) -> Result<[f64; 3], SolverError> {
    let def = elmer_fe::element_type(ptype.0).map_err(|e| SolverError::fatal("FindParentUVW", e.to_string()))?;
    let mut nu = [vec![0.0; face.len()], vec![0.0; face.len()], vec![0.0; face.len()]];
    for (i, f) in face.iter().enumerate() {
        let j = parent
            .iter()
            .position(|p| p == f)
            .ok_or_else(|| SolverError::fatal("FindParentUVW", "Could not find all nodes in parent!"))?;
        nu[0][i] = def.node_u[j];
        nu[1][i] = def.node_v[j];
        nu[2][i] = def.node_w[j];
    }
    Ok(std::array::from_fn(|a| (0..face.len()).fold(0.0, |acc, i| acc + fb[i] * nu[a][i])))
}

/// `NormalVector(Face, Nodes, u, v, .FALSE.)` of a line in 2D: the
/// unoriented unit normal `(-y', x')`.
fn line_normal(xyz: &[[f64; 3]], u: f64) -> [f64; 3] {
    let d: Vec<f64> = match xyz.len() {
        2 => vec![-0.5, 0.5],
        _ => vec![u - 0.5, u + 0.5, -2.0 * u],
    };
    let dot = |a: usize| d.iter().zip(xyz).fold(0.0, |s, (di, x)| s + di * x[a]);
    let (dxdu, dydu) = (dot(0), dot(1));
    let det = 1.0 / (dxdu * dxdu + dydu * dydu).sqrt();
    [-dydu * det, dxdu * det, 0.0]
}

impl Solver for AdvectionReaction {
    /// `AdvectionReactionSolver_init`: `Discontinuous Galerkin = True`.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        if p.logical("discontinuous galerkin") != Some(true) {
            return Err(unsupported("a run without Discontinuous Galerkin = True in the sif (the init adds it)"));
        }
        if p.logical("calculate nodal average") == Some(true) {
            // `Exported Variable = -nodal <var>_nodal` (CreateNodalPerm: the
            // nodes of the active elements in node order).
            let name = format!("{}_nodal", ctx.variable_name().unwrap_or("tracer"));
            if ctx.variable(&name).is_none() {
                let mut mark = vec![false; ctx.mesh.n_nodes()];
                for e in ctx.active_elements() {
                    for &k in ctx.mesh.element(e).nodes {
                        mark[k as usize] = true;
                    }
                }
                let mut j = 0u32;
                let perm: Vec<u32> = mark
                    .iter()
                    .map(|&m| {
                        if m {
                            j += 1;
                            j - 1
                        } else {
                            NO_DOF
                        }
                    })
                    .collect();
                let v = elmer_solver::Variable::new(&name, FieldKind::Nodal, 1, Some(Arc::from(perm)), vec![0.0; j as usize]);
                ctx.sim.vars.add(v)?;
            }
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("transient runs"));
        }
        if ctx.mesh.dim() != 2 {
            return Err(unsupported("3D meshes (faces)"));
        }
        let p = ctx.params().clone();
        if p.logical("limit solution") == Some(true) {
            return Err(unsupported("Limit Solution"));
        }
        if ctx.variable("mesh velocity").is_some() {
            return Err(unsupported("a Mesh Velocity"));
        }
        if ctx.solver_dofs().and_then(|d| d.dg_layout()).is_none() {
            return Err(unsupported("a solver without a DG layout"));
        }
        let name = ctx.variable_name().unwrap_or_default().to_string();
        let iter_max = p.integer("nonlinear system max iterations").unwrap_or(1).max(1);
        let fatal = |m: &str| SolverError::fatal(CALLER, m);
        let mut be = BasisEval::default();
        let mut pb = BasisEval::default();
        let mut slots = Vec::new();
        for _iter in 1..=iter_max {
            ctx.default_initialize()?;
            // Bulk.
            let elems: Vec<ElemId> = ctx.active_elements().collect();
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let nodes = el.nodes.to_vec();
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let bp = ctx.body_params(e).ok_or_else(|| fatal("an element without a body"))?;
                let load = at_nodes(ctx, bp.body_force, &format!("{name} source"), &nodes)?;
                let gamma = at_nodes(ctx, bp.material, &format!("{name} gamma"), &nodes)?;
                let velo = velocity(ctx, e, &nodes)?;
                let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| fatal("no FE backend"))?;
                let mut local = LocalSystem::zeros(n);
                for qp in fe.quadrature.gauss_points(etype, None) {
                    fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let (b, g) = (&be.basis[..n], &be.dbasis_dx[..n]);
                    let s = qp.weight * be.det_j;
                    let l = (0..n).fold(0.0, |acc, k| acc + load[k] * b[k]);
                    let gg = (0..n).fold(0.0, |acc, k| acc + b[k] * gamma[k]);
                    let div_mvelo = 0.0;
                    let cu: [f64; 2] = std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + b[k] * velo[k][i]));
                    for pp in 0..n {
                        for q in 0..n {
                            let a = &mut local.stiff[pp * n + q];
                            *a += s * (gg + div_mvelo) * b[q] * b[pp];
                            for (i, c) in cu.iter().enumerate() {
                                *a -= s * c * b[q] * g[pp][i];
                            }
                        }
                    }
                    for (pp, f) in local.force.iter_mut().enumerate() {
                        *f += s * l * b[pp];
                    }
                }
                ctx.default_update(&local, e)?;
            }

            // Interior edges, in Mesh % Edges order, both parents active.
            let active: std::collections::HashSet<ElemId> = elems.iter().copied().collect();
            for (lp, rp, face) in elmer_solver::dg_interfaces(ctx.mesh)? {
                if !active.contains(&lp) || !active.contains(&rp) {
                    continue;
                }
                let n = face.len();
                let ftype = ElementType(match n {
                    2 => 202,
                    3 => 203,
                    _ => return Err(unsupported("edges of more than 3 nodes")),
                });
                let (ln, rn) = (ctx.mesh.element(lp).nodes.to_vec(), ctx.mesh.element(rp).nodes.to_vec());
                let (lt, rt) = (ctx.mesh.element(lp).etype, ctx.mesh.element(rp).etype);
                let (n1, n2) = (ln.len(), rn.len());
                let fxyz: Vec<[f64; 3]> = face.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let lxyz: Vec<[f64; 3]> = ln.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let rxyz: Vec<[f64; 3]> = rn.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let velo = velocity(ctx, lp, &face)?;
                // LeftOut: face centre minus left parent centre.
                let left_out: [f64; 3] = std::array::from_fn(|a| {
                    fxyz.iter().map(|x| x[a]).sum::<f64>() / n as f64 - lxyz.iter().map(|x| x[a]).sum::<f64>() / n1 as f64
                });
                let m = n1 + n2;
                let mut local = LocalSystem::zeros(m);
                let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| fatal("no FE backend"))?;
                for qp in fe.quadrature.gauss_points(ftype, None) {
                    fe.basis.eval(ftype, &fxyz, qp.uvw, &mut be).map_err(|e| SolverError::fatal("ElementInfo", e))?;
                    let fb = be.basis[..n].to_vec();
                    let s = qp.weight * be.det_j;
                    let mut normal = line_normal(&fxyz, qp.uvw[0]);
                    if left_out.iter().zip(&normal).map(|(a, b)| a * b).sum::<f64>() < 0.0 {
                        normal = normal.map(|c| -c);
                    }
                    let mut jump = vec![0.0; m];
                    let mut avg = vec![0.0; m];
                    let uvw = parent_uvw(&face, &ln, lt, &fb)?;
                    fe.basis.eval(lt, &lxyz, uvw, &mut pb).map_err(|e| SolverError::fatal("ElementInfo", e))?;
                    for k in 0..n1 {
                        jump[k] = pb.basis[k];
                        avg[k] = pb.basis[k] / 2.0;
                    }
                    let uvw = parent_uvw(&face, &rn, rt, &fb)?;
                    fe.basis.eval(rt, &rxyz, uvw, &mut pb).map_err(|e| SolverError::fatal("ElementInfo", e))?;
                    for k in 0..n2 {
                        jump[n1 + k] = -pb.basis[k];
                        avg[n1 + k] = pb.basis[k] / 2.0;
                    }
                    let cu: [f64; 3] = std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + velo[k][i] * fb[k]));
                    let udotn = normal.iter().zip(&cu).map(|(a, b)| a * b).sum::<f64>();
                    for pp in 0..m {
                        for q in 0..m {
                            let a = &mut local.stiff[pp * m + q];
                            *a += s * udotn * avg[q] * jump[pp];
                            *a += s * udotn.abs() / 2.0 * jump[q] * jump[pp];
                        }
                    }
                }
                // GetElementDOFs of the face: the left, then the right parent's DG indexes.
                let d = ctx.solver_dofs().expect("checked above");
                let mut s2 = Vec::new();
                d.element_slots(lp, &mut slots);
                d.element_slots(rp, &mut s2);
                slots.extend_from_slice(&s2);
                if slots.contains(&NO_DOF) {
                    continue;
                }
                ctx.update_global_equations(&local, &slots)?;
            }
            ctx.default_finish_bulk_assembly()?;

            // Every boundary element whose dofs are all in the system.
            let nb = ctx.mesh.n_bulk();
            for t in nb..nb + ctx.mesh.n_boundary() {
                let e = ElemId(t as u32);
                let d = ctx.solver_dofs().expect("checked above");
                d.element_slots(e, &mut slots);
                if slots.is_empty() || slots.contains(&NO_DOF) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let parent = el.parents[0].or(el.parents[1]).ok_or_else(|| fatal("a boundary element without a parent"))?;
                let (enodes, etype) = (el.nodes.to_vec(), el.etype);
                let n = enodes.len();
                let pel = ctx.mesh.element(parent);
                let (pnodes, ptype) = (pel.nodes.to_vec(), pel.etype);
                let np = pnodes.len();
                let xyz: Vec<[f64; 3]> = enodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let pxyz: Vec<[f64; 3]> = pnodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let velo = velocity(ctx, parent, &enodes)?;
                let bc = ctx.bc_params(e).cloned();
                let (load, inflow_bc) = match bc.as_ref().and_then(|b| b.get(&name)) {
                    Some(_) => (at_nodes(ctx, bc.as_ref(), &name, &enodes)?, true),
                    None => (vec![0.0; n], false),
                };
                if matches!(bc.as_ref().and_then(|b| b.get(&name)), Some(Value::Logical(_))) {
                    return Err(unsupported("a logical boundary value"));
                }
                let normal0 = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [0.0, 0.0])?;
                let cu0: [f64; 3] = std::array::from_fn(|i| velo.iter().map(|v| v[i]).sum::<f64>() / n as f64);
                let inflow = inflow_bc && normal0.iter().zip(&cu0).map(|(a, b)| a * b).sum::<f64>() < 0.0;
                // The local system on the boundary element's dofs: its
                // (left) parent's, then a right parent's (left untouched).
                let m = slots.len();
                let mut local = LocalSystem::zeros(m);
                let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| fatal("no FE backend"))?;
                for qp in fe.quadrature.gauss_points(etype, None) {
                    let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [qp.uvw[0], qp.uvw[1]])?;
                    fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let b = be.basis[..n].to_vec();
                    let s = qp.weight * be.det_j;
                    let uvw = parent_uvw(&enodes, &pnodes, ptype, &b)?;
                    fe.basis.eval(ptype, &pxyz, uvw, &mut pb).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let pbas = &pb.basis[..np];
                    let l = (0..n).fold(0.0, |acc, k| acc + load[k] * b[k]);
                    let cu: [f64; 2] = std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + velo[k][i] * b[k]));
                    let udotn = normal[0] * cu[0] + normal[1] * cu[1] + normal[2] * 0.0;
                    for pp in 0..np {
                        if inflow {
                            local.force[pp] -= s * udotn * l * pbas[pp];
                        } else {
                            for q in 0..np {
                                local.stiff[pp * m + q] += s * udotn * pbas[q] * pbas[pp];
                            }
                        }
                    }
                }
                ctx.update_global_equations(&local, &slots)?;
            }
            ctx.default_finish_boundary_assembly()?;
            ctx.default_solve()?;
            let converged = ctx.convergence().is_some_and(|c| c.nonlin_converged == 1);
            if converged {
                break;
            }
        }
        if p.logical("calculate nodal average") == Some(true) {
            nodal_average(ctx, &name)?;
        }
        Ok(())
    }
}
