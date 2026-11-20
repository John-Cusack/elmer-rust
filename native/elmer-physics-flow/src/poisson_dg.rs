//! `PoissonSolver` of the library module `PoissonDG` (fem/src/modules/
//! PoissonDG.F90): the Poisson equation `-lap u = f` with symmetric
//! interior-penalty discontinuous Galerkin.
//!
//! - Bulk (`LocalMatrix`): `grad u . grad w` and the `Source`.
//! - Every interior edge (`LocalJumps`, or `LocalJumpsDiscontBC` on a
//!   boundary with `Discont BC`): `gamma/hE [u][w] - {dw/dn}[u] + [w]{du/dn}`,
//!   with `[w] = w_right - w_left`, the normal out of the left parent and
//!   `gamma = gamma1` (default 1000) or `beta` (default `gamma`), glued at
//!   the left, then the right parent's DG indexes.
//! - Dirichlet conditions through `DefaultDirichletBCs`.
//!
//! The module's `LocalMatrixBoundary` is never called. Ported for 2D
//! meshes. Refused: transient runs, 3D.

use elmer_core::{BasisEval, ElemId, ElementType, LocalSystem, ParamMap, Params, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::adv_react::{line_normal, parent_uvw};

const CALLER: &str = "PoissonDG";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(PoissonDg))
}

pub struct PoissonDg;

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("PoissonDG PoissonSolver: {what}"))
}

impl Solver for PoissonDg {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("transient runs"));
        }
        if ctx.mesh.dim() != 2 {
            return Err(unsupported("3D meshes (faces)"));
        }
        if ctx.solver_dofs().and_then(|d| d.dg_layout()).is_none() {
            return Err(unsupported("a solver without Discontinuous Galerkin"));
        }
        let p = ctx.params().clone();
        let fatal = |m: &str| SolverError::fatal(CALLER, m);
        let mut be = BasisEval::default();
        let mut pb = BasisEval::default();

        ctx.default_initialize()?;
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        for &e in &elems {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let etype = el.etype;
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let bf = ctx.body_params(e).and_then(|b| b.body_force).cloned();
            let load = match &bf {
                Some(b) => ctx.real_at_nodes(b, "source", e)?.unwrap_or_else(|| vec![0.0; n]),
                None => vec![0.0; n],
            };
            let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| fatal("no FE backend"))?;
            let mut local = LocalSystem::zeros(n);
            for qp in fe.quadrature.gauss_points(etype, None) {
                fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let (b, g) = (&be.basis[..n], &be.dbasis_dx[..n]);
                let s = qp.weight * be.det_j;
                let l = (0..n).fold(0.0, |acc, k| acc + load[k] * b[k]);
                for i in 0..n {
                    for j in 0..n {
                        // MATMUL(dBasisdx, TRANSPOSE(dBasisdx)) over the 3 columns.
                        let m = (0..3).fold(0.0, |acc, a| acc + g[i][a] * g[j][a]);
                        local.stiff[i * n + j] += s * m;
                    }
                }
                for (i, f) in local.force.iter_mut().enumerate() {
                    *f += s * l * b[i];
                }
            }
            ctx.default_update(&local, e)?;
        }

        let gamma = p.real("gamma1").unwrap_or(1.0e3);
        let beta = p.real("beta").unwrap_or(gamma);
        // The boundary condition of an edge: that of the boundary element on
        // the same nodes (AssignConstraints).
        let nb = ctx.mesh.n_bulk();
        let mut bc_of: std::collections::HashMap<Vec<u32>, ElemId> = std::collections::HashMap::new();
        for t in nb..nb + ctx.mesh.n_boundary() {
            let mut key = ctx.mesh.element(ElemId(t as u32)).nodes.to_vec();
            key.sort_unstable();
            bc_of.insert(key, ElemId(t as u32));
        }
        let mut lslots = Vec::new();
        let mut rslots = Vec::new();
        for (lp, rp, face) in elmer_solver::dg_edges(ctx.mesh)? {
            let Some(rp) = rp else { continue };
            // ActiveBoundaryElement(elm, DGBoundary = .TRUE.): all of the left
            // or all of the right parent's DG indexes in the system.
            {
                let d = ctx.solver_dofs().expect("checked above");
                d.element_slots(lp, &mut lslots);
                d.element_slots(rp, &mut rslots);
            }
            if lslots.contains(&NO_DOF) && rslots.contains(&NO_DOF) {
                continue;
            }
            let mut key = face.clone();
            key.sort_unstable();
            let discont = bc_of
                .get(&key)
                .and_then(|&b| ctx.bc_params(b))
                .and_then(|bc| bc.logical("discont bc"))
                == Some(true);
            let g = if discont { beta } else { gamma };
            let n = face.len();
            let ftype = ElementType(match n {
                2 => 202,
                3 => 203,
                _ => return Err(unsupported("edges of more than 3 nodes")),
            });
            let (ln, rn) = (ctx.mesh.element(lp).nodes.to_vec(), ctx.mesh.element(rp).nodes.to_vec());
            let (lt, rt) = (ctx.mesh.element(lp).etype, ctx.mesh.element(rp).etype);
            let (nl, nr) = (ln.len(), rn.len());
            let fxyz: Vec<[f64; 3]> = face.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let lxyz: Vec<[f64; 3]> = ln.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let rxyz: Vec<[f64; 3]> = rn.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let he = elmer_fe::element_diameter(ftype, &fxyz, false);
            let left_out: [f64; 3] = std::array::from_fn(|a| {
                fxyz.iter().map(|x| x[a]).sum::<f64>() / n as f64 - lxyz.iter().map(|x| x[a]).sum::<f64>() / nl as f64
            });
            let m = nl + nr;
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
                let mut flux = vec![0.0; m];
                let dn = |g: &[f64; 3]| g.iter().zip(&normal).map(|(a, b)| a * b).sum::<f64>();
                let uvw = parent_uvw(&face, &ln, lt, &fb)?;
                fe.basis.eval(lt, &lxyz, uvw, &mut pb).map_err(|e| SolverError::fatal("ElementInfo", e))?;
                for k in 0..nl {
                    jump[k] = -pb.basis[k];
                    flux[k] = dn(&pb.dbasis_dx[k]) / 2.0;
                }
                let uvw = parent_uvw(&face, &rn, rt, &fb)?;
                fe.basis.eval(rt, &rxyz, uvw, &mut pb).map_err(|e| SolverError::fatal("ElementInfo", e))?;
                for k in 0..nr {
                    jump[nl + k] = pb.basis[k];
                    flux[nl + k] = dn(&pb.dbasis_dx[k]) / 2.0;
                }
                for pp in 0..m {
                    for q in 0..m {
                        let a = &mut local.stiff[pp * m + q];
                        *a += (g / he) * jump[pp] * jump[q] * s;
                        *a -= flux[pp] * jump[q] * s;
                        *a += jump[pp] * flux[q] * s;
                    }
                }
            }
            // In parallel a face between partitions is half on each.
            let Some(share) = ctx.partition_share([Some(lp), Some(rp)]) else { continue };
            if share != 1.0 {
                local.stiff.iter_mut().for_each(|v| *v /= 2.0);
                local.force.iter_mut().for_each(|v| *v /= 2.0);
            }
            // GetElementDOFs of the edge: the left, then the right parent's.
            let mut slots = lslots.clone();
            slots.extend_from_slice(&rslots);
            ctx.update_global_equations(&local, &slots)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}
