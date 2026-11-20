//! `MeshSolver` (fem/src/modules/MeshSolve.F90): moves the mesh nodes by
//! solving linear elasticity for the `Mesh Update` displacement (`Mesh
//! Elastic Modulus` / `Youngs Modulus`, default 1, and `Mesh Poisson Ratio`
//! / `Poisson Ratio`, default 0.25), with `Mesh Force i`, `Mesh Normal
//! Force` and `Mesh Penalty Factor i` boundary terms. The previous update is
//! taken back from the mesh before the assembly and the new one added after
//! the solve (`DisplaceMesh`), unless `Skip Displace Mesh`.
//!
//! Transient runs compute the first-order `Mesh Velocity` of
//! `MeshSolver_Init`. Refused: `Compute Mesh Velocity Order`, p
//! elements, and a `Displacement` field to combine with (the `TPerm` path).

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::util::{apply_coords, coords_of, has_prefix};

const CALLER: &str = "MeshSolve";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(MeshSolver::default()))
}

#[derive(Default)]
pub struct MeshSolver {
    /// `AllocationsDone`: a solve has displaced the mesh before.
    visited: bool,
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

/// `DisplaceMesh(Mesh, Update, sgn, Perm, DOFs)`: `x += sgn * u` at every
/// node with a slot, the first `dofs` coordinates.
fn displace(ctx: &SolverCtx<'_>, sgn: f64) -> Result<(), SolverError> {
    let f = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
    let dofs = f.dofs;
    let old = coords_of(ctx.mesh);
    let mut new = old.clone();
    for (i, c) in new.iter_mut().enumerate() {
        let k = match f.perm {
            Some(p) => {
                if p[i] == NO_DOF {
                    continue;
                }
                p[i] as usize
            }
            None => i,
        };
        for (d, x) in c.iter_mut().enumerate().take(dofs.min(3)) {
            *x += sgn * f.values[dofs * k + d];
        }
    }
    apply_coords(CALLER, &old, &new)
}

/// `Mesh Velocity` of a transient run (MeshSolve.F90:370-430), first
/// order: `(MeshUpdate - PrevValues(:,1)) / dt` on the velocity's nodes,
/// unless `Compute Mesh Velocity = False` (zero) or `Skip First Mesh
/// Velocity` on the first visit (left as it is).
fn mesh_velocity(ctx: &mut SolverCtx<'_>, dt: f64, first: bool) -> Result<(), SolverError> {
    let p = ctx.params();
    let compute = p.logical("compute mesh velocity").unwrap_or(true);
    let skip_first = first && compute && p.logical("skip first mesh velocity") == Some(true);
    if skip_first {
        ctx.info("MeshSolve", "Skipping computation of initial Mesh Velocity", 3);
        return Ok(());
    }
    let Some(vperm) = ctx.variable("mesh velocity").map(|v| v.perm.map(|q| q.to_vec())) else { return Ok(()) };
    if !compute {
        ctx.variable_mut("mesh velocity").expect("checked").values.iter_mut().for_each(|x| *x = 0.0);
        return Ok(());
    }
    let (upd, perm, dofs, prev) = {
        let name = ctx.variable_name().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?.to_string();
        let v = ctx.sim.vars.get(&name).ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
        let prev = v.prev_values.first().cloned().ok_or_else(|| unsupported("a mesh velocity without the previous mesh update"))?;
        (v.values.clone(), v.perm.as_deref().map(|q| q.to_vec()), v.dofs, prev)
    };
    let n = ctx.mesh.n_nodes();
    let vel = ctx.variable_mut("mesh velocity").expect("checked");
    if vel.dofs != dofs {
        return Err(unsupported("a Mesh Velocity with other dofs than the mesh update"));
    }
    for i in 0..n {
        let l = vperm.as_ref().map_or(i as u32, |q| q[i]);
        if l == NO_DOF {
            continue;
        }
        let j = perm.as_ref().map_or(i as u32, |q| q[i]);
        if j == NO_DOF {
            return Err(SolverError::fatal(CALLER, "a mesh velocity node without a mesh update"));
        }
        for m in 0..dofs {
            vel.values[l as usize * dofs + m] = (upd[j as usize * dofs + m] - prev[j as usize * dofs + m]) / dt;
        }
    }
    Ok(())
}

impl Solver for MeshSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient && ctx.params().get("compute mesh velocity order").is_some() {
            return Err(unsupported("Compute Mesh Velocity Order (BDF mesh velocities)"));
        }
        if ctx.solver_dofs().is_some() {
            return Err(unsupported("p elements"));
        }
        let p = ctx.params().clone();
        let stdofs = ctx.primary().map(|f| f.dofs).ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
        let mut dim = ctx.mesh.dim();
        if stdofs < dim {
            ctx.warn("MeshSolve", "Displacement dofs smaller than dim, assuming reduced dimension!");
            dim = stdofs;
        }
        let skip = p.logical("skip displace mesh") == Some(true);
        let displace_first = p.logical("first time non-zero") == Some(true);
        let ignore = p.logical("ignore displacement") == Some(true);
        if !skip && !ignore && ctx.sim.vars.get("displacement").is_some() {
            return Err(unsupported("a Displacement field to combine with (TPerm)"));
        }
        if skip {
            ctx.info("MeshSolver", "Skipping the displacement of mesh!", 5);
        } else if self.visited || displace_first {
            displace(ctx, -1.0)?;
        }
        let fe = ctx
            .sim
            .backends
            .fe
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
        let _ = fe;
        ctx.info("MeshSolve", "MESH UPDATE SOLVER:", 4);
        ctx.info("MeshSolve", "Starting assembly...", 4);
        let maxiter = p.integer("nonlinear system max iterations").unwrap_or(1).max(1);
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let do_bc = ctx.sim.model.boundary_conditions.iter().any(|b| b.get("mesh normal force").is_some() || has_prefix(b, "mesh force"));
        let bnd: Vec<ElemId> = if do_bc { ctx.boundary_elements().collect() } else { Vec::new() };
        let mut dofs = Vec::new();
        for _ in 0..maxiter {
            ctx.default_initialize()?;
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let mat = ctx.body_params(e).and_then(|b| b.material).cloned().unwrap_or_default();
                let nodal = |ctx: &SolverCtx<'_>, keys: [&str; 2], default: f64| -> Result<Vec<f64>, SolverError> {
                    for k in keys {
                        if let Some(v) = ctx.real_at_nodes(&mat, k, e)? {
                            return Ok(v);
                        }
                    }
                    Ok(vec![default; n])
                };
                let young = nodal(ctx, ["Mesh Elastic Modulus", "Youngs Modulus"], 1.0)?;
                let poisson = nodal(ctx, ["Mesh Poisson Ratio", "Poisson Ratio"], 0.25)?;
                let lame1: Vec<f64> =
                    young.iter().zip(&poisson).map(|(&y, &nu)| y * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))).collect();
                let lame2: Vec<f64> = young.iter().zip(&poisson).map(|(&y, &nu)| y / (2.0 * (1.0 + nu))).collect();
                let m = dim * n;
                let mut local = LocalSystem::zeros(m);
                let fe = ctx.sim.backends.fe.as_ref().expect("checked above");
                let mut be = BasisEval::default();
                for qp in fe.quadrature.gauss_points(el.etype, None) {
                    fe.basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(|msg| SolverError::fatal("ElementInfo", msg))?;
                    let s = be.det_j * qp.weight;
                    let b = &be.basis[..n];
                    let db = &be.dbasis_dx[..n];
                    let interp = |v: &[f64]| v.iter().zip(b).fold(0.0, |a, (&x, &bi)| a + x * bi);
                    let (l1, l2) = (interp(&lame1), interp(&lame2));
                    for pp in 0..n {
                        for qq in pp..n {
                            for i in 0..dim {
                                for j in 0..dim {
                                    let at = |r: usize, c: usize| (dim * pp + r) * m + dim * qq + c;
                                    local.stiff[at(i, j)] += s * l1 * db[qq][j] * db[pp][i];
                                    local.stiff[at(i, i)] += s * l2 * db[qq][j] * db[pp][j];
                                    local.stiff[at(i, j)] += s * l2 * db[qq][i] * db[pp][j];
                                }
                            }
                        }
                    }
                }
                // The lower triangle from the upper one.
                for r in 0..m {
                    for c in 0..r {
                        local.stiff[r * m + c] = local.stiff[c * m + r];
                    }
                }
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let Some(bc) = ctx.bc_params(e).cloned() else { continue };
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let mut alpha = vec![vec![0.0; n]; dim];
                for (i, a) in alpha.iter_mut().enumerate() {
                    if let Some(v) = ctx.real_at_nodes(&bc, &format!("Mesh Penalty Factor {}", i + 1), e)? {
                        *a = v;
                    }
                }
                let mut got = false;
                let mut load = vec![vec![0.0; n]; 3];
                for (i, l) in load.iter_mut().enumerate() {
                    if let Some(v) = ctx.real_at_nodes(&bc, &format!("Mesh Force {}", i + 1), e)? {
                        *l = v;
                        got = true;
                    }
                }
                let beta = match ctx.real_at_nodes(&bc, "Mesh Normal Force", e)? {
                    Some(v) => {
                        got = true;
                        v
                    }
                    None => vec![0.0; n],
                };
                if !got {
                    continue;
                }
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let m = dim * n;
                let mut local = LocalSystem::zeros(m);
                let fe = ctx.sim.backends.fe.as_ref().expect("checked above");
                let mut be = BasisEval::default();
                for qp in fe.quadrature.gauss_points(el.etype, None) {
                    fe.basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(|msg| SolverError::fatal("ElementInfo", msg))?;
                    let s = be.det_j * qp.weight;
                    let b = &be.basis[..n];
                    let interp = |v: &[f64]| v.iter().zip(b).fold(0.0, |a, (&x, &bi)| a + x * bi);
                    let mut load_ip = [0.0; 3];
                    let mut alpha_ip = [0.0; 3];
                    for i in 0..dim {
                        load_ip[i] = interp(&load[i]);
                        alpha_ip[i] = interp(&alpha[i]);
                    }
                    let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [qp.uvw[0], qp.uvw[1]])?;
                    let bn = interp(&beta);
                    for (l, nv) in load_ip.iter_mut().zip(normal) {
                        *l += bn * nv;
                    }
                    for pp in 0..n {
                        for qq in 0..n {
                            for (i, &a) in alpha_ip.iter().enumerate().take(dim) {
                                local.stiff[(dim * pp + i) * m + dim * qq + i] += s * a * b[qq] * b[pp];
                            }
                        }
                    }
                    for qq in 0..n {
                        for (i, &l) in load_ip.iter().enumerate().take(dim) {
                            local.force[dim * qq + i] += s * b[qq] * l;
                        }
                    }
                }
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_boundary_assembly()?;
            ctx.info("MeshSolve", "Assembly done", 4);
            ctx.default_dirichlet()?;
            ctx.info("MeshSolve", "Set boundaries done", 4);
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        if transient {
            mesh_velocity(ctx, dt, !self.visited)?;
        }
        if !skip {
            displace(ctx, 1.0)?;
        }
        self.visited = true;
        Ok(())
    }
}
