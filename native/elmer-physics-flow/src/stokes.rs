//! `StokesSolver` (fem/src/modules/Stokes.F90): the Stokes equations on p
//! elements (`Element = "p:k"`), the velocity components and the pressure
//! in one `dim + 1` dof variable. The velocity uses all the element's dofs,
//! the pressure only the vertex functions: the pressure rows of the other
//! dofs are the identity (P(k)/P1 pairs). Dirichlet conditions are the
//! solver's own `SetBoundaryConditions`: the row of each boundary node a
//! condition gives becomes the identity with the value on the right-hand
//! side, the edge dofs of the element get zero.
//!
//! Ported: `Convective = False` (the plain Stokes solve) without `Block
//! Preconditioning`. Refused: convection (Picard, Newton, `hybrid`), the
//! block-preconditioned GCR outer iteration, grad-div stabilization and
//! transient runs.

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, NO_DOF};
use elmer_solver::backend::DirectRow;
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "NavierStokesSolver";

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("Stokes StokesSolver: {what}"))
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(StokesSolver))
}

pub struct StokesSolver;

impl Solver for StokesSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        if transient {
            return Err(unsupported("transient runs"));
        }
        if p.logical("convective").unwrap_or(true) {
            return Err(unsupported("Convective (the Picard/Newton iteration)"));
        }
        if p.logical("block preconditioning") == Some(true) {
            return Err(unsupported("Block Preconditioning (GCR outer iteration)"));
        }
        if p.logical("grad-div stabilization") == Some(true) {
            return Err(unsupported("Grad-Div Stabilization"));
        }
        if p.get("nonlinear system relaxation factor").is_some() {
            return Err(unsupported("Nonlinear System Relaxation Factor"));
        }
        let dim = ctx.mesh.dim();
        let c = dim + 1;
        let (dofs, name) = {
            let f = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
            (f.dofs, ctx.variable_name().unwrap_or_default().to_string())
        };
        if dofs != c {
            return Err(unsupported("a variable without dim + 1 dofs"));
        }
        ctx.default_initialize()?;
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut slots = Vec::new();
        for &e in &elems {
            let d = ctx.solver_dofs().ok_or_else(|| unsupported("nodal elements (it needs Element = p:k)"))?;
            d.element_slots(e, &mut slots);
            let nd = slots.len();
            if d.local_bubbles(e) > 0 {
                return Err(unsupported("local bubbles"));
            }
            let order = d.element_nodes(e);
            let basis = d.element_basis_scalar(e).map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
            let el = ctx.mesh.element(e);
            if order != el.nodes {
                return Err(unsupported("p elements whose dofs reorder the nodes"));
            }
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let bp = ctx.body_params(e);
            let (mat, bf) = (bp.and_then(|b| b.material).cloned(), bp.and_then(|b| b.body_force).cloned());
            let nodal = |ctx: &SolverCtx<'_>, s: &Option<ParamMap>, key: &str| -> Result<Vec<f64>, SolverError> {
                match s {
                    Some(s) => Ok(ctx.real_at_nodes(s, key, e)?.unwrap_or_else(|| vec![0.0; n])),
                    None => Ok(vec![0.0; n]),
                }
            };
            let mut load = vec![vec![0.0; n]; 3];
            if bf.is_some() {
                for (k, l) in load.iter_mut().enumerate() {
                    *l = nodal(ctx, &bf, &format!("Force {}", k + 1))?;
                }
            }
            // `Density` only enters the mass matrix of a transient run.
            let mu = nodal(ctx, &mat, "Viscosity")?;
            let m = c * nd;
            let mut local = LocalSystem::zeros(m);
            let mut be = BasisEval::default();
            let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).fold(0.0, |s, (x, y)| s + x * y);
            for qp in &rule {
                basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(|msg| SolverError::fatal("ElementInfo", msg))?;
                if be.basis.len() < nd {
                    return Err(SolverError::fatal("ElementInfo", "fewer basis functions than element dofs"));
                }
                let b = &be.basis[..nd];
                let g = &be.dbasis_dx[..nd];
                let s = qp.weight * be.det_j;
                let mu_ip = dot(&b[..n], &mu);
                // MATMUL(Basis(1:n), LOAD(1:n,1:dim+1)).
                let load_ip: Vec<f64> = (0..c).map(|j| dot(&b[..n], &load[j])).collect();
                for pp in 0..nd {
                    for qq in 0..nd {
                        let at = |i: usize, j: usize| (c * pp + i) * m + c * qq + j;
                        for i in 0..dim {
                            for j in 0..dim {
                                local.stiff[at(i, i)] += s * mu_ip * g[qq][j] * g[pp][j];
                                local.stiff[at(i, j)] += s * mu_ip * g[qq][i] * g[pp][j];
                            }
                            if qq < n {
                                local.stiff[at(i, dim)] -= s * b[qq] * g[pp][i];
                            }
                            if pp < n {
                                local.stiff[at(dim, i)] -= s * g[qq][i] * b[pp];
                            }
                        }
                    }
                    for (k, &l) in load_ip.iter().enumerate() {
                        local.force[c * pp + k] += s * l * b[pp];
                    }
                }
            }
            // The pressure rows of the edge (and face) dofs.
            for pp in n..nd {
                let i = c * (pp + 1) - 1;
                local.force[i] = 0.0;
                for j in 0..m {
                    local.stiff[i * m + j] = 0.0;
                }
                local.stiff[i * m + i] = 1.0;
            }
            ctx.default_update(&local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;

        // SetBoundaryConditions for components 1 .. dim+1: the boundary
        // elements in order, each condition with the element's tag.
        let mut rows: Vec<DirectRow> = Vec::new();
        let nb = ctx.mesh.n_bulk();
        let (perm, comps) = {
            let f = ctx.primary().expect("checked above");
            (f.perm.map(|p| p.to_vec()), ctx.sim.vars.component_names(ctx.sim.slots.iter().find(|s| s.id == ctx.solver_id()).and_then(|s| s.variable).expect("variable")))
        };
        let slot_of = |node: u32| -> Option<usize> {
            match &perm {
                Some(p) => (p[node as usize] != NO_DOF).then(|| p[node as usize] as usize),
                None => Some(node as usize),
            }
        };
        for dof in 0..c {
            let key = comps.get(dof).cloned().unwrap_or_else(|| format!("{name} {}", dof + 1));
            for t in nb..nb + ctx.mesh.n_boundary() {
                let e = ElemId(t as u32);
                let Some(bc) = ctx.bc_params(e).cloned() else { continue };
                let Some(vals) = ctx.real_at_nodes(&bc, &key, e)? else { continue };
                let nodes = ctx.mesh.element(e).nodes.to_vec();
                for (&node, &x) in nodes.iter().zip(&vals) {
                    if let Some(sl) = slot_of(node) {
                        let k = (c * sl + dof) as u32;
                        rows.push(DirectRow { row: k, entries: vec![(k, 1.0)], rhs: x });
                    }
                }
                if let Some(d) = ctx.solver_dofs() {
                    d.element_slots(e, &mut slots);
                    for &sl in slots.iter().skip(nodes.len()) {
                        if sl == NO_DOF {
                            return Err(SolverError::fatal(CALLER, "a boundary edge dof outside the system"));
                        }
                        let k = (c * sl as usize + dof) as u32;
                        rows.push(DirectRow { row: k, entries: vec![(k, 1.0)], rhs: 0.0 });
                    }
                }
            }
        }
        ctx.apply_dirichlet_values_with_rows(&[], &rows)?;
        ctx.default_solve()?;
        Ok(())
    }
}
