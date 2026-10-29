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

/// HeatSolve's diffusion with constant material values and a body-force
/// `Heat Source` (any keyword value, so MATC and field dependence work),
/// iterated `Nonlinear System Max Iterations` times until
/// `default_converged`: stiffness `k grad N . grad N`, mass `rho c_p N N`
/// through `Default1stOrderTime`, load `rho * source` interpolated at the
/// Gauss points, and the Robin boundary term.
pub struct HeatNl;

pub fn heat_nl(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(HeatNl))
}

impl Solver for HeatNl {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        use elmer_core::Params;
        let iters = ctx.params().integer("nonlinear system max iterations").unwrap_or(1);
        for _ in 0..iters {
            ctx.default_initialize()?;
            let active: Vec<ElemId> = ctx.active_elements().collect();
            let mut ev = BasisEval::default();
            for e in active {
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let b = ctx.body_params(e).expect("body");
                let mat = b.material.expect("material");
                let (k, rho, cp) = (
                    mat.real("heat conductivity").unwrap_or(1.0),
                    mat.real("density").unwrap_or(1.0),
                    mat.real("heat capacity").unwrap_or(1.0),
                );
                let source = match b.body_force {
                    Some(bf) => ctx.real_at_nodes(bf, "heat source", e)?.unwrap_or_else(|| vec![0.0; n]),
                    None => vec![0.0; n],
                };
                let load: Vec<f64> = source.iter().map(|s| s * rho).collect();
                let mut local = LocalSystem::zeros(n);
                let mut mass = vec![0.0; n * n];
                for qp in super::quad::TestQuadrature.gauss_points(el.etype, None) {
                    super::basis::TestBasis
                        .eval(el.etype, &xyz, qp.uvw, &mut ev)
                        .map_err(|m| SolverError::fatal("HeatSolve", m))?;
                    let w = qp.weight * ev.det_j;
                    let f: f64 = ev.basis.iter().zip(&load).map(|(a, b)| a * b).sum();
                    for p in 0..n {
                        for q in 0..n {
                            let g = (ev.dbasis_dx[p], ev.dbasis_dx[q]);
                            local.stiff[p * n + q] += w * k * (g.0[0] * g.1[0] + g.0[1] * g.1[1] + g.0[2] * g.1[2]);
                            mass[p * n + q] += w * rho * cp * ev.basis[p] * ev.basis[q];
                        }
                        local.force[p] += w * f * ev.basis[p];
                    }
                }
                local.mass = Some(mass);
                if transient {
                    ctx.default_first_order_time(&mut local, e, dt)?;
                }
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;
            // Robin: `Heat Transfer Coefficient` h and `External Temperature`
            // T_ext add h N N and h T_ext N.
            let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
            for e in boundary {
                let n = ctx.mesh.element(e).nodes.len();
                let Some(bc) = ctx.bc_params(e) else { continue };
                let Some(h) = ctx.real_at_nodes(bc, "heat transfer coefficient", e)? else { continue };
                let text = ctx.real_at_nodes(bc, "external temperature", e)?.unwrap_or_else(|| vec![0.0; n]);
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let mut local = LocalSystem::zeros(n);
                let mut ev = BasisEval::default();
                for qp in super::quad::TestQuadrature.gauss_points(el.etype, None) {
                    super::basis::TestBasis
                        .eval(el.etype, &xyz, qp.uvw, &mut ev)
                        .map_err(|m| SolverError::fatal("HeatSolve", m))?;
                    let w = qp.weight * ev.det_j;
                    let hh: f64 = ev.basis.iter().zip(&h).map(|(a, b)| a * b).sum();
                    let tt: f64 = ev.basis.iter().zip(&text).map(|(a, b)| a * b).sum();
                    for p in 0..n {
                        for q in 0..n {
                            local.stiff[p * n + q] += w * hh * ev.basis[p] * ev.basis[q];
                        }
                        local.force[p] += w * hh * tt * ev.basis[p];
                    }
                }
                ctx.default_update(&local, e)?;
            }
            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        ctx.default_finish()
    }
}
