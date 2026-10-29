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

/// Transient `T' - T'' = 0` with unit conductivity, density and heat
/// capacity: HeatSolver's diffusion terms with consistent mass, through
/// `Default1stOrderTime`.
pub struct Heat;

pub fn heat(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Heat))
}

impl Solver for Heat {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut ev = BasisEval::default();
        for e in active {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let mut local = LocalSystem::zeros(n);
            let mut mass = vec![0.0; n * n];
            for qp in super::quad::TestQuadrature.gauss_points(el.etype, None) {
                super::basis::TestBasis
                    .eval(el.etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|m| SolverError::fatal("HeatSolve", m))?;
                let w = qp.weight * ev.det_j;
                for p in 0..n {
                    for q in 0..n {
                        let g = (ev.dbasis_dx[p], ev.dbasis_dx[q]);
                        local.stiff[p * n + q] += w * (g.0[0] * g.1[0] + g.0[1] * g.1[1] + g.0[2] * g.1[2]);
                        mass[p * n + q] += w * ev.basis[p] * ev.basis[q];
                    }
                }
            }
            local.mass = Some(mass);
            if transient {
                ctx.default_first_order_time(&mut local, e, dt)?;
            }
            ctx.default_update(&local, e)?;
        }
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}

/// HelmholtzSolve's bulk matrix (HelmholtzSolve.F90 `LocalMatrix`) for
/// constant material values: `((M + iD) N_q N_p + grad N_q . grad N_p) /
/// rho` with `k = omega / c`, `M = -k^2`, `D = k * damping`, assembled in
/// Elmer's interleaved real form (`[[Re, -Im], [Im, Re]]` per complex entry).
pub struct Helmholtz;

pub fn helmholtz(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Helmholtz))
}

impl Solver for Helmholtz {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        use elmer_core::Params;
        ctx.default_initialize()?;
        let omega = ctx.simulation_params().real("angular frequency").expect("Angular Frequency");
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut ev = BasisEval::default();
        for e in active {
            let mat = ctx.body_params(e).and_then(|b| b.material).expect("material");
            let (c, damp, rho) = (
                mat.real("sound speed").unwrap(),
                mat.real("sound damping").unwrap_or(0.0),
                mat.real("density").unwrap_or(1.0),
            );
            let k = omega / c;
            let (m, d) = (-k * k, k * damp);
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&i| ctx.mesh.coord(i)).collect();
            let mut local = LocalSystem::zeros(2 * n);
            for qp in super::quad::TestQuadrature.gauss_points(el.etype, None) {
                super::basis::TestBasis
                    .eval(el.etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|msg| SolverError::fatal("HelmholtzSolve", msg))?;
                let s = qp.weight * ev.det_j;
                for p in 0..n {
                    for q in 0..n {
                        let g = (ev.dbasis_dx[q], ev.dbasis_dx[p]);
                        let grad = g.0[0] * g.1[0] + g.0[1] * g.1[1] + g.0[2] * g.1[2];
                        let re = s * ((m * ev.basis[q] * ev.basis[p] + grad) / rho);
                        let im = s * ((d * ev.basis[q] * ev.basis[p]) / rho);
                        let w = 2 * n;
                        local.stiff[2 * p * w + 2 * q] += re;
                        local.stiff[2 * p * w + 2 * q + 1] -= im;
                        local.stiff[(2 * p + 1) * w + 2 * q] += im;
                        local.stiff[(2 * p + 1) * w + 2 * q + 1] += re;
                    }
                }
            }
            ctx.default_update(&local, e)?;
        }
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}
