//! `Q1Q0/Q1Q0.f90`: `SUBROUTINE StokesSolver`, steady Navier-Stokes on
//! Q1 velocity with an elementwise constant pressure carried by the
//! element's one bubble dof (`Element = "b:1"`, `Bubbles In Global System =
//! False`). The variable has only the `dim` velocity components; the
//! pressure row (`l = dim*n + 1`, with a `1e-8` diagonal) is condensed out
//! locally (`LCondensateP`). Before that, the elementwise pressure of the
//! previous velocity, `-SUM(u * K(l,:)) / K(l,l)`, is averaged into the
//! nodes of the exported `aPressure`. Convection (`Convect`, default on) is
//! Picard, Newton for good once successive norms are within 1 % (`SAVE
//! Newton, Norm`).
//!
//! Ported: the path the case takes, the bubble outside the global system
//! (`n == nd`). Refused: transient runs, global bubbles.

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::nodal;

const CALLER: &str = "StokesSolve";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Q1Q0 { newton: false, norm: 0.0 }))
}

pub struct Q1Q0 {
    newton: bool,
    norm: f64,
}

impl Solver for Q1Q0 {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        let fatal = |m: String| SolverError::fatal(CALLER, m);
        if transient {
            return Err(SolverError::Unsupported("Q1Q0 StokesSolver: transient runs".into()));
        }
        let dim = crate::hooks::coordinate_system(ctx)?.dim;
        let dofs = ctx.primary().ok_or_else(|| fatal("no solver variable".into()))?.dofs;
        if dofs != dim {
            return Err(SolverError::Unsupported(format!("Q1Q0 StokesSolver: a variable of {dofs} dofs in {dim}D")));
        }
        let convect = ctx.params().logical("convect").unwrap_or(true);
        // P % Values = 0 and the counts, on aPressure's own permutation.
        let p_perm: Option<Vec<u32>> = ctx.variable("apressure").and_then(|v| (v.dofs == 1).then(|| v.perm.map(|p| p.to_vec()))).flatten();
        let p_size = ctx.variable("apressure").map_or(0, |v| v.values.len());
        let mut p_sum = vec![0.0; p_size];
        let mut p_count = vec![0u32; p_size];

        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        if ctx.sim.backends.fe.is_none() {
            return Err(fatal("no FE backend".into()));
        }
        let mut slots = Vec::new();
        for e in active {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let etype = el.etype;
            let nodes = el.nodes.to_vec();
            let (nd, nb) = match ctx.solver_dofs() {
                Some(d) => {
                    d.element_slots(e, &mut slots);
                    (slots.len(), d.local_bubbles(e))
                }
                None => (n, 1),
            };
            if nd != n || nb != 1 {
                return Err(SolverError::Unsupported(format!(
                    "Q1Q0 StokesSolver: {nd} element dofs and {nb} bubbles (only the local b:1 pressure)"
                )));
            }
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let (m, bf) = ctx.body_params(e).map_or((None, None), |b| (b.material, b.body_force));
            let mut load = vec![[0.0; 3]; n];
            for (k, key) in ["source x", "source y", "source z"].into_iter().enumerate() {
                for (l, x) in load.iter_mut().zip(nodal(ctx, bf, key, e, &nodes)?) {
                    l[k] = x;
                }
            }
            let rho = nodal(ctx, m, "density", e, &nodes)?;
            let mu = nodal(ctx, m, "viscosity", e, &nodes)?;
            // GetVectorLocalSolution.
            ctx.element_dofs(e, &mut slots);
            let velo: Vec<[f64; 3]> = {
                let v = ctx.primary().ok_or_else(|| fatal("no solver variable".into()))?;
                (0..n)
                    .map(|i| std::array::from_fn(|k| if k < dim && slots[i * dofs + k] != NO_DOF { v.values[slots[i * dofs + k] as usize] } else { 0.0 }))
                    .collect()
            };

            // LocalMatrix: velocity rows dim*(p-1)+i, the pressure row l.
            let l = dim * n;
            let m_sz = l + 1;
            let mut stiff = vec![0.0; m_sz * m_sz];
            let mut force = vec![0.0; m_sz];
            let fe = ctx.sim.backends.fe.as_ref().expect("checked above");
            let mut be = BasisEval::default();
            for qp in fe.quadrature.gauss_points(etype, None) {
                fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let (b, g) = (&be.basis[..n], &be.dbasis_dx[..n]);
                let s = qp.weight * be.det_j;
                let mu_ip = (0..n).fold(0.0, |acc, k| acc + b[k] * mu[k]);
                let rho_ip = (0..n).fold(0.0, |acc, k| acc + b[k] * rho[k]);
                let mut vel = [0.0; 3];
                let mut grad = [[0.0; 3]; 3];
                for i in 0..dim {
                    vel[i] = (0..n).fold(0.0, |acc, k| acc + velo[k][i] * b[k]);
                    for j in 0..dim {
                        grad[i][j] = (0..n).fold(0.0, |acc, k| acc + velo[k][i] * g[k][j]);
                    }
                }
                let mut load_ip = [0.0; 3];
                for (j, li) in load_ip.iter_mut().enumerate().take(dim) {
                    *li = (0..n).fold(0.0, |acc, k| acc + b[k] * load[k][j]);
                }
                if convect && self.newton {
                    for i in 0..dim {
                        load_ip[i] += rho_ip * (0..dim).fold(0.0, |acc, j| acc + grad[i][j] * vel[j]);
                    }
                }
                for p in 0..n {
                    for q in 0..n {
                        let at = |i: usize, j: usize| (dim * p + i) * m_sz + dim * q + j;
                        for i in 0..dim {
                            for j in 0..dim {
                                stiff[at(i, i)] += s * mu_ip * g[q][j] * g[p][j];
                                stiff[at(i, j)] += s * mu_ip * g[q][i] * g[p][j];
                                if convect {
                                    stiff[at(i, i)] += s * rho_ip * vel[j] * g[q][j] * b[p];
                                    if self.newton {
                                        stiff[at(i, j)] += s * rho_ip * grad[i][j] * b[q] * b[p];
                                    }
                                }
                            }
                        }
                    }
                    for i in 0..dim {
                        stiff[(dim * p + i) * m_sz + l] -= s * g[p][i];
                        stiff[l * m_sz + dim * p + i] += s * rho_ip * g[p][i];
                    }
                    for i in 0..dim {
                        force[dim * p + i] += s * load_ip[i] * b[p];
                    }
                }
                stiff[l * m_sz + l] += s * 1.0e-8;
            }

            // The elementwise pressure of the current velocity.
            if p_size > 0 {
                let mut pressure = 0.0;
                for j in 0..dofs {
                    for (i, v) in velo.iter().enumerate() {
                        pressure -= v[j] * stiff[l * m_sz + dofs * i + j];
                    }
                }
                pressure /= stiff[l * m_sz + l];
                for &k in &nodes {
                    let s = p_perm.as_ref().map_or(k, |p| p[k as usize]);
                    if s != NO_DOF {
                        p_sum[s as usize] += pressure;
                        p_count[s as usize] += 1;
                    }
                }
            }

            // LCondensateP(n*DOFs, STIFF).
            let kbb = stiff[l * m_sz + l];
            let mut local = LocalSystem::zeros(l);
            for i in 0..l {
                for j in 0..l {
                    local.stiff[i * l + j] = stiff[i * m_sz + j] - stiff[i * m_sz + l] * stiff[l * m_sz + j] / kbb;
                }
                local.force[i] = force[i];
            }
            ctx.default_update(&local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        let prev = self.norm;
        self.norm = ctx.default_solve()?;
        if p_size > 0 {
            if let Some(v) = ctx.variable_mut("apressure") {
                for (i, x) in v.values.iter_mut().enumerate() {
                    *x = if p_count[i] > 0 { p_sum[i] / f64::from(p_count[i]) } else { 0.0 };
                }
            }
        }
        if 2.0 * (prev - self.norm).abs() / (prev + self.norm) < 1.0e-2 {
            self.newton = true;
        }
        Ok(())
    }
}
