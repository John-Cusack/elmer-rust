//! `1dtests/Poisson.f90` as a test-local solver.

use elmer_core::{BasisEval, ElemId, ElementBasis, LocalSystem, ParamMap, Quadrature};
use elmer_solver::{Solver, SolverCtx, SolverError};

/// `1dtests/Poisson.f90`: `-u'' = Source` with `Source` from the body force.
pub struct Poisson;

pub fn poisson(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Poisson))
}

impl Solver for Poisson {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut ev = BasisEval::default();
        for e in active {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let load = match ctx.body_params(e).and_then(|b| b.body_force) {
                Some(bf) => ctx.real_at_nodes(bf, "source", e)?.unwrap_or_else(|| vec![0.0; n]),
                None => vec![0.0; n],
            };
            let mut local = LocalSystem::zeros(n);
            for qp in super::quad::TestQuadrature.gauss_points(el.etype, None) {
                super::basis::TestBasis
                    .eval(el.etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|m| SolverError::fatal("PoissonSolve", m))?;
                let w = qp.weight * ev.det_j;
                let f: f64 = ev.basis.iter().zip(&load).map(|(b, l)| b * l).sum();
                for p in 0..n {
                    for q in 0..n {
                        let g = (ev.dbasis_dx[p], ev.dbasis_dx[q]);
                        local.stiff[p * n + q] += w * (g.0[0] * g.1[0] + g.0[1] * g.1[1] + g.0[2] * g.1[2]);
                    }
                    local.force[p] += w * f * ev.basis[p];
                }
            }
            ctx.default_update(&local, e)?;
        }
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}
