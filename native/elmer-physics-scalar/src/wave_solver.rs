//! `WaveSolver` (fem/src/modules/WaveSolver.F90): the compressional
//! pressure wave
//!
//! ```text
//!   1/(rho c^2) d2p/dt2 + (att/(rho c^2)) grad.grad dp/dt
//!     + react/(rho c^2) dp/dt - 1/rho div grad p = -s/rho
//! ```
//!
//! with `Sound Speed` `c`, `Sound Damping` `att`, `Sound Reaction Damping`
//! `react` and the optional `Density` `rho` (Material), `Sound Source` `s`
//! (Body Force), and on boundaries `Outflow Boundary` / `Plane Wave BC`
//! (the damping `1/(rho c)`) and `Source Acceleration`.
//!
//! `WaveSolver_Init` makes the time derivative order two. Refused: the
//! harmonic and eigen modes (`DefaultUpdateMass`/`Damp`) and bubbles.

use elmer_core::{ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::common::{fe_provider, init_dirichlet, real_at_nodes, unsupported};

/// `Procedure = "WaveSolver" "WaveSolver"`.
pub const PROCEDURE: (&str, &str) = ("WaveSolver", "WaveSolver");

const CALLER: &str = "WaveSolver";

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| {
        let (basis, quad) = fe_provider()?();
        Ok(Box::new(WaveSolver { basis, quad }))
    });
}

struct WaveSolver {
    basis: Box<dyn ElementBasis + Send>,
    quad: Box<dyn Quadrature + Send>,
}

fn flag(p: &ParamMap, key: &str) -> bool {
    p.logical(key) == Some(true)
}

/// A keyword of an element at its nodes, for `ListGetElementReal`: a
/// constant as it is, anything else at the nodes (then interpolated).
enum AtIp {
    Absent,
    Constant(f64),
    Nodal(Vec<f64>),
}

impl AtIp {
    fn get(ctx: &SolverCtx<'_>, section: Option<&ParamMap>, key: &str, e: ElemId) -> Result<AtIp, SolverError> {
        let Some(p) = section else { return Ok(AtIp::Absent) };
        Ok(match p.get(key) {
            None => AtIp::Absent,
            Some(Value::Real { data, .. }) if data.len() == 1 => AtIp::Constant(data[0]),
            Some(Value::Integer(v)) if v.len() == 1 => AtIp::Constant(f64::from(v[0])),
            Some(_) => match real_at_nodes(ctx, Some(p), key, e)? {
                Some(v) => AtIp::Nodal(v),
                None => AtIp::Absent,
            },
        })
    }

    /// `ListGetElementReal(h, Basis, Element, Found)`.
    fn at(&self, basis: &[f64]) -> Option<f64> {
        match self {
            AtIp::Absent => None,
            AtIp::Constant(c) => Some(*c),
            AtIp::Nodal(v) => {
                let mut s = 0.0;
                for (b, x) in basis.iter().zip(v) {
                    s += b * x;
                }
                Some(s)
            }
        }
    }
}

/// `MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )(p,q)`.
fn grad_dot(d: &[[f64; 3]], p: usize, q: usize) -> f64 {
    let mut s = 0.0;
    for k in 0..3 {
        s += d[p][k] * d[q][k];
    }
    s
}

impl Solver for WaveSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        // WaveSolver_Init: `ListAddInteger( SolverParams, 'Time derivative order', 2 )`.
        ctx.set_time_derivative_order(2)?;
        init_dirichlet(ctx)?;
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "Solving the compressional pressure wave", 4);
        let p = ctx.params();
        let need_mass = ctx.nof_eigen() > 0 && flag(p, "eigen analysis")
            || flag(p, "harmonic mode")
            || flag(p, "harmonic analysis");
        if need_mass {
            return Err(unsupported("WaveSolver harmonic or eigen mode", "DefaultUpdateMass/Damp are not ported here"));
        }
        if flag(p, "bubbles") {
            return Err(unsupported("Solver keyword \"Bubbles\"", "bubble condensation is not ported"));
        }
        let max_iter = p.integer("nonlinear system max iterations").map_or(1, |n| n.max(1));
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut dofs = Vec::new();
        let mut ev = elmer_core::BasisEval::default();
        for _ in 1..=max_iter {
            ctx.default_initialize()?;
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let nd = xyz.len();
                let bp = ctx
                    .body_params(e)
                    .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
                let (mat, bf) = (bp.material, bp.body_force);
                let speed = AtIp::get(ctx, mat, "sound speed", e)?;
                if matches!(speed, AtIp::Absent) {
                    return Err(SolverError::fatal("ListGetElementReal", "Failed to find \"Sound speed\""));
                }
                let damping = AtIp::get(ctx, mat, "sound damping", e)?;
                let density = AtIp::get(ctx, mat, "density", e)?;
                let react = AtIp::get(ctx, mat, "sound reaction damping", e)?;
                let load = AtIp::get(ctx, bf, "sound source", e)?;
                let mut local = LocalSystem::zeros(nd);
                local.mass = Some(vec![0.0; nd * nd]);
                local.damp = Some(vec![0.0; nd * nd]);
                for qp in self.quad.gauss_points(etype, None) {
                    self.basis
                        .eval(etype, &xyz, qp.uvw, &mut ev)
                        .map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let b = &ev.basis[..nd];
                    let d = &ev.dbasis_dx[..nd];
                    let c = speed.at(b).expect("checked");
                    let att = damping.at(b);
                    let reaction = react.at(b);
                    let load_ip = load.at(b);
                    let mut weight = qp.weight * ev.det_j;
                    if let Some(rho) = density.at(b) {
                        weight /= rho;
                    }
                    for p in 0..nd {
                        for q in 0..nd {
                            local.stiff[p * nd + q] += weight * grad_dot(d, p, q);
                        }
                    }
                    if let Some(att) = att {
                        let dm = local.damp.as_mut().expect("set");
                        for p in 0..nd {
                            for q in 0..nd {
                                dm[p * nd + q] += weight * (att / c.powi(2)) * grad_dot(d, p, q);
                            }
                        }
                    }
                    for p in 0..nd {
                        for q in 0..nd {
                            if let Some(r) = reaction {
                                let dm = local.damp.as_mut().expect("set");
                                dm[p * nd + q] += weight * (r / c.powi(2)) * b[q] * b[p];
                            }
                            let m = local.mass.as_mut().expect("set");
                            m[p * nd + q] += weight * (1.0 / c.powi(2)) * b[q] * b[p];
                        }
                    }
                    if let Some(s) = load_ip {
                        for p in 0..nd {
                            local.force[p] -= weight * s * b[p];
                        }
                    }
                }
                if transient {
                    ctx.default_second_order_time(&mut local, e, dt)?;
                }
                local.mass = None;
                local.damp = None;
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let Some(bc) = ctx.bc_params(e).cloned() else { continue };
                let outflow = bc.logical("plane wave bc").or_else(|| bc.logical("outflow boundary")) == Some(true);
                let flux = AtIp::get(ctx, Some(&bc), "source acceleration", e)?;
                if !outflow && matches!(flux, AtIp::Absent) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let etype = el.etype;
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let nd = xyz.len();
                // `ListGetElementRealParent`: the parent's material at the
                // boundary nodes.
                let parent = el
                    .parents
                    .iter()
                    .flatten()
                    .next()
                    .copied()
                    .ok_or_else(|| SolverError::fatal(CALLER, format!("boundary element {} has no parent element", e.0)))?;
                let pmat = ctx.body_params(parent).and_then(|b| b.material);
                let speed = AtIp::get(ctx, pmat, "sound speed", e)?;
                let density = AtIp::get(ctx, pmat, "density", e)?;
                if matches!(density, AtIp::Absent) {
                    return Err(SolverError::fatal("ListGetElementRealParent", "Failed to find \"Density\""));
                }
                let mut local = LocalSystem::zeros(nd);
                local.mass = Some(vec![0.0; nd * nd]);
                local.damp = Some(vec![0.0; nd * nd]);
                for qp in self.quad.gauss_points(etype, None) {
                    self.basis
                        .eval(etype, &xyz, qp.uvw, &mut ev)
                        .map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let b = &ev.basis[..nd];
                    let mut weight = qp.weight * ev.det_j;
                    if let Some(rho) = density.at(b) {
                        weight /= rho;
                    }
                    if outflow {
                        let c = speed
                            .at(b)
                            .ok_or_else(|| SolverError::fatal("ListGetElementRealParent", "Failed to find \"Sound speed\""))?;
                        let dm = local.damp.as_mut().expect("set");
                        for p in 0..nd {
                            for q in 0..nd {
                                dm[p * nd + q] += weight * b[q] * b[p] / c;
                            }
                        }
                    }
                    if let Some(g) = flux.at(b) {
                        for p in 0..nd {
                            local.force[p] -= weight * g * b[p];
                        }
                    }
                }
                if transient {
                    ctx.default_second_order_time(&mut local, e, dt)?;
                }
                local.mass = None;
                local.damp = None;
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_boundary_assembly()?;
            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
            if flag(ctx.params(), "set average to zero") {
                let v = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
                let mut s = 0.0;
                for x in v.values.iter() {
                    s += x;
                }
                let ave = s / v.values.len() as f64;
                for x in v.values.iter_mut() {
                    *x -= ave;
                }
            }
        }
        Ok(())
    }
}
