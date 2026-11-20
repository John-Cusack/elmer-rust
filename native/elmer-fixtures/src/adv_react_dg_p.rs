//! `fem/tests/AdvReactDG_P/AdvReactDG_P.F90`, the case's own copy of an
//! older `AdvectionReactionSolver` (strategy D12: a test-local plugin, not
//! product code): the steady advection-reaction equation with
//! discontinuous Galerkin on p elements (`Element = "p:N d:M"`), where
//! every loop runs over `GetElementNOFDOFs` basis functions of the p basis
//! while the coefficients are interpolated with the first (nodal) ones.
//!
//! - Bulk (`LocalMatrix`): `gamma c w - c (u . grad w)`, the source `f w`.
//! - Interior edges (`LocalJumps`): `u.n {c} [w] + |u.n|/2 [c] [w]`, on
//!   the edge's p integration rule, the parents' p bases at the edge's
//!   points (`FindParentUVW`).
//! - Boundary elements (`LocalMatrixBoundary`): inflow `-u.n c_in w` when
//!   the field is given and the mean `u.n < 0`, else `u.n c w`.
//! - Afterwards `Exported Variable 1` becomes a nodal field over every
//!   mesh node: the mean of the DG values at each node of the active
//!   elements.
//!
//! Refused: `Convection = computed`, a mesh velocity, transient runs, 3D,
//! elements other than quadrilaterals (`FindParentUVW` on the p reference
//! triangle is not ported) and a layout without p elements.

use std::sync::Arc;

use elmer_core::{BasisEval, ElemId, FieldKind, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::expr::eval_at_nodes_with;
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "AdvectionReaction";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(AdvReactDgP))
}

pub struct AdvReactDgP;

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("AdvReactDG_P AdvectionReactionSolver: {what}"))
}

fn fatal(m: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, m)
}

/// `GetReal(section, key, Found, Element)` at `nodes`; zeros without it.
fn at_nodes(ctx: &SolverCtx<'_>, section: Option<&ParamMap>, key: &str, nodes: &[u32]) -> Result<Vec<f64>, SolverError> {
    let Some(v) = section.and_then(|s| s.get(key)) else {
        return Ok(vec![0.0; nodes.len()]);
    };
    let b = &ctx.sim.backends;
    eval_at_nodes_with(v, nodes, &ctx.fields(), b.matc.as_deref(), &b.udfs)
}

/// `GetLocalALEVelocity` for `Convection = constant | none`.
fn velocity(ctx: &SolverCtx<'_>, parent: ElemId, nodes: &[u32]) -> Result<Vec<[f64; 3]>, SolverError> {
    let bp = ctx.body_params(parent).ok_or_else(|| fatal("an element without a body"))?;
    let flag = bp
        .equation
        .and_then(|q| q.string("convection"))
        .map(elmer_core::normalize_key)
        .ok_or_else(|| fatal("No string for keyword >Convection< found in Equation"))?;
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
        other => return Err(fatal(format!("Convection flag >{other}< not recognised"))),
    }
    Ok(v)
}

/// `FindParentUVW` on a quadrilateral parent: the face's nodal basis
/// times the parent's local node coordinates.
fn parent_uvw(face: &[u32], parent: &[u32], fb: &[f64]) -> Result<[f64; 3], SolverError> {
    let def = elmer_fe::element_type(404).map_err(|e| SolverError::fatal("FindParentUVW", e.to_string()))?;
    let mut uvw = [0.0; 3];
    for (i, f) in face.iter().enumerate() {
        let j = parent
            .iter()
            .position(|p| p == f)
            .ok_or_else(|| SolverError::fatal("FindParentUVW", "Could not find all nodes in parent!"))?;
        uvw[0] += fb[i] * def.node_u[j];
        uvw[1] += fb[i] * def.node_v[j];
    }
    Ok(uvw)
}

/// `NormalVector(Face, Nodes, u, v, .FALSE.)` of a linear line in 2D.
fn line_normal(xyz: &[[f64; 3]]) -> [f64; 3] {
    let dxdu = 0.5 * (xyz[1][0] - xyz[0][0]);
    let dydu = 0.5 * (xyz[1][1] - xyz[0][1]);
    let det = 1.0 / (dxdu * dxdu + dydu * dydu).sqrt();
    [-dydu * det, dxdu * det, 0.0]
}

/// The linear line's nodal basis and `detJ` at `u` (the first two p basis
/// functions of an edge are these).
fn line_basis(xyz: &[[f64; 3]], u: f64) -> ([f64; 2], f64) {
    let dx = 0.5 * (xyz[1][0] - xyz[0][0]);
    let dy = 0.5 * (xyz[1][1] - xyz[0][1]);
    ([(1.0 - u) / 2.0, (1.0 + u) / 2.0], (dx * dx + dy * dy).sqrt())
}

impl Solver for AdvReactDgP {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("transient runs"));
        }
        if ctx.mesh.dim() != 2 {
            return Err(unsupported("3D meshes"));
        }
        if ctx.variable("mesh velocity").is_some() {
            return Err(unsupported("a Mesh Velocity"));
        }
        let p = ctx.params().clone();
        let exp_name = p
            .string("exported variable 1")
            .map(str::to_string)
            .ok_or_else(|| fatal("No value > Exported Variable 1 < found in Solver"))?;
        let d = ctx.solver_dofs().ok_or_else(|| unsupported("a solver without p elements"))?;
        if d.dg_layout().is_none() {
            return Err(unsupported("a solver without a DG layout"));
        }
        let name = ctx.variable_name().unwrap_or_default().to_string();
        let iter_max = p.integer("nonlinear system max iterations").unwrap_or(1).max(1);
        let mut be = BasisEval::default();
        let mut pb = BasisEval::default();
        let mut slots = Vec::new();
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        for &e in &elems {
            if ctx.mesh.element(e).etype.0 != 404 || !ctx.solver_dofs().is_some_and(|d| d.is_p_element(e)) {
                return Err(unsupported("elements other than p quadrilaterals"));
            }
        }
        for _iter in 1..=iter_max {
            ctx.default_initialize()?;
            for &e in &elems {
                let d = ctx.solver_dofs().expect("checked");
                d.element_slots(e, &mut slots);
                let nd = slots.len();
                let nodes = d.element_nodes(e);
                let n = nodes.len();
                let etype = ctx.mesh.element(e).etype;
                let basis = d.element_basis(e).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let bp = ctx.body_params(e).ok_or_else(|| fatal("an element without a body"))?;
                let load = at_nodes(ctx, bp.body_force, &format!("{name} source"), &nodes)?;
                let gamma = at_nodes(ctx, bp.material, &format!("{name} gamma"), &nodes)?;
                let velo = velocity(ctx, e, &nodes)?;
                let mut local = LocalSystem::zeros(nd);
                for qp in &rule {
                    basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let (b, g) = (&be.basis[..nd], &be.dbasis_dx[..nd]);
                    let s = qp.weight * be.det_j;
                    let l = (0..n).fold(0.0, |acc, k| acc + load[k] * b[k]);
                    let gg = (0..n).fold(0.0, |acc, k| acc + b[k] * gamma[k]);
                    let cu: [f64; 2] = std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + b[k] * velo[k][i]));
                    for pp in 0..nd {
                        for q in 0..nd {
                            let a = &mut local.stiff[pp * nd + q];
                            *a += s * gg * b[q] * b[pp];
                            for (i, c) in cu.iter().enumerate() {
                                *a -= s * c * b[q] * g[pp][i];
                            }
                        }
                    }
                    for (pp, f) in local.force.iter_mut().enumerate() {
                        *f += s * l * b[pp];
                    }
                }
                ctx.update_global_equations(&local, &slots)?;
            }

            for (lp, rp, face) in elmer_solver::dg_interfaces(ctx.mesh)? {
                if face.len() != 2 {
                    return Err(unsupported("edges of more than 2 nodes"));
                }
                let d = ctx.solver_dofs().expect("checked");
                let mut s2 = Vec::new();
                d.element_slots(lp, &mut slots);
                d.element_slots(rp, &mut s2);
                let (nd1, nd2) = (slots.len(), s2.len());
                slots.extend_from_slice(&s2);
                if slots.contains(&NO_DOF) {
                    continue;
                }
                let rule = match d.edge_rule([face[0], face[1]]) {
                    Some(r) => r.map_err(|m| SolverError::fatal("GaussPoints", m))?,
                    None => return Err(unsupported("an edge without p definitions")),
                };
                let (ln, rn) = (d.element_nodes(lp), d.element_nodes(rp));
                let (lbas, rbas) = (
                    d.element_basis(lp).map_err(|m| SolverError::fatal("ElementInfo", m))?,
                    d.element_basis(rp).map_err(|m| SolverError::fatal("ElementInfo", m))?,
                );
                let (lt, rt) = (ctx.mesh.element(lp).etype, ctx.mesh.element(rp).etype);
                let fxyz: Vec<[f64; 3]> = face.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let lxyz: Vec<[f64; 3]> = ln.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let rxyz: Vec<[f64; 3]> = rn.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let velo = velocity(ctx, lp, &face)?;
                let left_out: [f64; 3] = std::array::from_fn(|a| {
                    fxyz.iter().map(|x| x[a]).sum::<f64>() / 2.0 - lxyz.iter().map(|x| x[a]).sum::<f64>() / ln.len() as f64
                });
                let m = nd1 + nd2;
                let mut local = LocalSystem::zeros(m);
                for qp in &rule {
                    let (fb, det_j) = line_basis(&fxyz, qp.uvw[0]);
                    let s = qp.weight * det_j;
                    let mut normal = line_normal(&fxyz);
                    if left_out.iter().zip(&normal).map(|(a, b)| a * b).sum::<f64>() < 0.0 {
                        normal = normal.map(|c| -c);
                    }
                    let mut jump = vec![0.0; m];
                    let mut avg = vec![0.0; m];
                    lbas.eval(lt, &lxyz, parent_uvw(&face, &ln, &fb)?, &mut pb).map_err(|e| SolverError::fatal("ElementInfo", e))?;
                    for k in 0..nd1 {
                        jump[k] = pb.basis[k];
                        avg[k] = pb.basis[k] / 2.0;
                    }
                    rbas.eval(rt, &rxyz, parent_uvw(&face, &rn, &fb)?, &mut pb).map_err(|e| SolverError::fatal("ElementInfo", e))?;
                    for k in 0..nd2 {
                        jump[nd1 + k] = -pb.basis[k];
                        avg[nd1 + k] = pb.basis[k] / 2.0;
                    }
                    let cu: [f64; 3] = std::array::from_fn(|i| (0..2).fold(0.0, |acc, k| acc + velo[k][i] * fb[k]));
                    let udotn = normal.iter().zip(&cu).map(|(a, b)| a * b).sum::<f64>();
                    for pp in 0..m {
                        for q in 0..m {
                            let a = &mut local.stiff[pp * m + q];
                            *a += s * udotn * avg[q] * jump[pp];
                            *a += s * udotn.abs() / 2.0 * jump[q] * jump[pp];
                        }
                    }
                }
                ctx.update_global_equations(&local, &slots)?;
            }
            ctx.default_finish_bulk_assembly()?;

            let nb = ctx.mesh.n_bulk();
            for t in nb..nb + ctx.mesh.n_boundary() {
                let e = ElemId(t as u32);
                let d = ctx.solver_dofs().expect("checked");
                d.element_slots(e, &mut slots);
                if slots.is_empty() || slots.contains(&NO_DOF) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let parent = el.parents[0].or(el.parents[1]).ok_or_else(|| fatal("a boundary element without a parent"))?;
                let enodes = d.element_nodes(e);
                let n = enodes.len();
                if n != 2 {
                    return Err(unsupported("boundary elements of more than 2 nodes"));
                }
                let pnodes = d.element_nodes(parent);
                let ptype = ctx.mesh.element(parent).etype;
                let pbas = d.element_basis(parent).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
                let xyz: Vec<[f64; 3]> = enodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let pxyz: Vec<[f64; 3]> = pnodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let velo = velocity(ctx, parent, &enodes)?;
                let bc = ctx.bc_params(e).cloned();
                if matches!(bc.as_ref().and_then(|b| b.get(&name)), Some(Value::Logical(_))) {
                    return Err(unsupported("a logical boundary value"));
                }
                let (load, inflow_bc) = match bc.as_ref().and_then(|b| b.get(&name)) {
                    Some(_) => (at_nodes(ctx, bc.as_ref(), &name, &enodes)?, true),
                    None => (vec![0.0; n], false),
                };
                let normal0 = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [0.0, 0.0])?;
                let cu0: [f64; 3] = std::array::from_fn(|i| velo.iter().map(|v| v[i]).sum::<f64>() / n as f64);
                let inflow = inflow_bc && normal0.iter().zip(&cu0).map(|(a, b)| a * b).sum::<f64>() < 0.0;
                let npd = slots.len();
                let mut local = LocalSystem::zeros(npd);
                for qp in &rule {
                    let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [qp.uvw[0], qp.uvw[1]])?;
                    let (b, det_j) = line_basis(&xyz, qp.uvw[0]);
                    let s = qp.weight * det_j;
                    pbas.eval(ptype, &pxyz, parent_uvw(&enodes, &pnodes, &b)?, &mut pb)
                        .map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let l = (0..n).fold(0.0, |acc, k| acc + load[k] * b[k]);
                    let cu: [f64; 2] = std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + velo[k][i] * b[k]));
                    let udotn = normal[0] * cu[0] + normal[1] * cu[1];
                    for pp in 0..npd {
                        if inflow {
                            local.force[pp] -= s * udotn * l * pb.basis[pp];
                        } else {
                            for q in 0..npd {
                                local.stiff[pp * npd + q] += s * udotn * pb.basis[q] * pb.basis[pp];
                            }
                        }
                    }
                }
                ctx.update_global_equations(&local, &slots)?;
            }
            ctx.default_finish_boundary_assembly()?;
            ctx.default_solve()?;
            if ctx.convergence().is_some_and(|c| c.nonlin_converged == 1) {
                break;
            }
        }
        nodal_result(ctx, &exp_name, &elems)
    }
}

/// The plugin's tail: `Exported Variable 1` turned into a nodal field over
/// every mesh node (Perm `i -> i`), the mean of the DG values at each node
/// of the active elements.
fn nodal_result(ctx: &mut SolverCtx<'_>, var: &str, elems: &[ElemId]) -> Result<(), SolverError> {
    if ctx.variable(var).is_none() {
        return Err(fatal(format!("Exported Variable >{var}< not found")));
    }
    let (sperm, svals) = {
        let v = ctx.primary().ok_or_else(|| fatal("no variable"))?;
        (v.perm.map(|p| p.to_vec()).ok_or_else(|| fatal("no Perm"))?, v.values.to_vec())
    };
    let ix = ctx.solver_dofs().and_then(|d| d.dg_layout()).cloned().ok_or_else(|| fatal("no DG layout"))?;
    let n1 = ctx.mesh.n_nodes();
    let mut sum = vec![0.0; n1];
    let mut cnt = vec![0u32; n1];
    for &e in elems {
        for (i, &node) in ctx.mesh.element(e).nodes.iter().enumerate() {
            let dg = ix.index[e.0 as usize][i];
            if dg == NO_DOF || sperm[dg as usize] == NO_DOF {
                continue;
            }
            sum[node as usize] += svals[sperm[dg as usize] as usize];
            cnt[node as usize] += 1;
        }
    }
    let id = ctx.sim.vars.id(var).expect("checked");
    let v = ctx.sim.vars.by_id_mut(id);
    let len = v.values.len().max(n1);
    let mut perm = vec![NO_DOF; len];
    for (i, p) in perm.iter_mut().enumerate().take(n1) {
        *p = i as u32;
    }
    let mut values = vec![0.0; len];
    for i in 0..n1 {
        values[i] = if cnt[i] > 0 { sum[i] / f64::from(cnt[i]) } else { 0.0 };
    }
    v.kind = FieldKind::Nodal;
    v.perm = Some(Arc::from(perm));
    v.values = values;
    v.values_changed = true;
    Ok(())
}
