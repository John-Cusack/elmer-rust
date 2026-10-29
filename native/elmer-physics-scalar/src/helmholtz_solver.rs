//! `HelmholtzSolver`: the `Procedure = "HelmholtzSolve" "HelmholtzSolver"`
//! solver, driving the kernels of [`crate::helmholtz`] through [`SolverCtx`]
//! the way `HelmholtzSolver` (HelmholtzSolve.F90:99) does.
//!
//! The pressure is a two-dof variable (`Pressure 1` real, `Pressure 2`
//! imaginary). The Fortran marks the system `Linear System Complex` and lets
//! the complex Krylov methods and the complex UMFPACK work on the interleaved
//! storage; `elmer-solver`'s `default_solve` does not take that flag yet, so
//! the interleaved real system is solved instead. It is the same linear
//! system, so the solution is the same; only the iterative methods'
//! convergence path and the matrix scaling differ.

use elmer_core::{ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};
use num_complex::Complex64;

use crate::common::{axisymmetric, fe_provider, init_dirichlet, real_at_nodes, unsupported};
use crate::helmholtz::{
    helmholtz_boundary_local, helmholtz_bulk_local, HelmholtzBcAtNodes, HelmholtzBulkAtNodes,
    HelmholtzConfig, HelmholtzLocal,
};
use crate::HELMHOLTZ_PROCEDURE;

const CALLER: &str = "HelmholtzSolve";

/// Registers `("HelmholtzSolve", "HelmholtzSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(
        HELMHOLTZ_PROCEDURE.0,
        HELMHOLTZ_PROCEDURE.1,
        helmholtz_factory,
    );
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`HelmholtzSolver`].
/// Fails with `NotWired` until [`install_fe`](crate::install_fe) has been called.
pub fn helmholtz_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    let (basis, quad) = fe_provider()?();
    Ok(Box::new(HelmholtzSolver::new(basis, quad)))
}

pub struct HelmholtzSolver {
    basis: Box<dyn ElementBasis + Send>,
    quad: Box<dyn Quadrature + Send>,
    complex: HelmholtzLocal,
    local: LocalSystem,
}

impl HelmholtzSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis,
            quad,
            complex: HelmholtzLocal::default(),
            local: LocalSystem::zeros(0),
        }
    }
}

fn flag(p: &ParamMap, key: &str) -> bool {
    p.logical(key) == Some(true)
}

/// Refuses what this port does not implement, instead of dropping it.
fn check_supported(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    if flag(ctx.params(), "bubbles") {
        return Err(unsupported(
            "Solver keyword \"Bubbles\"",
            "bubble condensation is not ported",
        ));
    }
    let model = &ctx.sim.model;
    for m in &model.materials {
        if m.0.keys().any(|k| k.starts_with("convection velocity")) {
            return Err(unsupported(
                "Material keyword \"Convection Velocity\"",
                "the boundary convection terms need NormalVector(Element, Nodes, u, v, check orientation)",
            ));
        }
    }
    for bc in &model.boundary_conditions {
        for key in ["flow interface", "structure interface"] {
            if flag(bc, key) {
                return Err(unsupported(
                    format!("Boundary Condition \"{key}\""),
                    "harmonic flow and structure coupling is not ported",
                ));
            }
        }
    }
    Ok(())
}

/// `GetAngularFrequency()` without a list or element (Lists.F90:10543):
/// `Frequency` (times 2 pi) then `Angular Frequency` in Simulation, then in
/// the solver section, then in the equation and the material of the first
/// active element.
fn angular_frequency(ctx: &SolverCtx<'_>) -> Option<f64> {
    fn from(p: &ParamMap) -> Option<f64> {
        p.real("frequency")
            .map(|f| 2.0 * std::f64::consts::PI * f)
            .or_else(|| p.real("angular frequency"))
    }
    from(ctx.simulation_params())
        .or_else(|| from(ctx.params()))
        .or_else(|| {
            let bp = ctx.body_params(ctx.active_elements().next()?)?;
            bp.equation
                .and_then(from)
                .or_else(|| bp.material.and_then(from))
        })
}

impl Solver for HelmholtzSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        check_supported(ctx)?;
        axisymmetric(ctx)?;
        let dofs = ctx.primary().map(|v| v.dofs);
        if dofs != Some(2) {
            return Err(SolverError::fatal(
                CALLER,
                format!("the pressure is complex: \"Variable DOFs\" must be 2, got {dofs:?}"),
            ));
        }
        init_dirichlet(ctx)?;
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        check_supported(ctx)?;
        let axi = axisymmetric(ctx)?;
        let dim = ctx.mesh.dim();
        let omega = angular_frequency(ctx)
            .ok_or_else(|| SolverError::fatal(CALLER, "Could not figure out Frequency!"))?;
        ctx.sim.model.simulation.insert(
            "res: Frequency",
            Value::real(omega / (2.0 * std::f64::consts::PI)),
        );
        if (flag(
            ctx.params(),
            "linear system preconditioning damp coefficient",
        ) || ctx
            .params()
            .get("linear system preconditioning damp coefficient")
            .is_some()
            || ctx
                .params()
                .get("linear system preconditioning damp coefficient im")
                .is_some())
            && ctx
                .params()
                .string("linear system solver")
                .map(str::to_ascii_lowercase)
                .as_deref()
                == Some("iterative")
        {
            ctx.warn(
                CALLER,
                "Damped preconditioning matrix is not built; the solution is unaffected, the iteration may differ",
            );
        }

        let mut cfg = HelmholtzConfig::new(dim, omega);
        cfg.axisymmetric = axi;
        let maxiter = ctx
            .params()
            .integer("nonlinear system max iterations")
            .map_or(1, |n| n.max(1));
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut dofs = Vec::new();

        for iter in 1..=maxiter {
            ctx.info(CALLER, " ", 4);
            ctx.info(CALLER, "-------------------------------------", 4);
            ctx.info(CALLER, &format!(" Helmholtz iteration {iter}"), 4);
            ctx.info(
                CALLER,
                &format!(" Frequency (Hz): {}", omega / (2.0 * std::f64::consts::PI)),
                4,
            );
            ctx.info(CALLER, "-------------------------------------", 4);
            ctx.default_initialize()?;

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let mat = {
                    let bp = ctx.body_params(e).ok_or_else(|| {
                        SolverError::fatal(CALLER, format!("element {} has no body", e.0))
                    })?;
                    let required = |key: &str| -> Result<Vec<f64>, SolverError> {
                        real_at_nodes(ctx, bp.material, key, e)?.ok_or_else(|| {
                            SolverError::fatal("ListGetReal", format!("Failed to find \"{key}\""))
                        })
                    };
                    HelmholtzBulkAtNodes {
                        sound_speed: required("Sound Speed")?,
                        density: required("Density")?,
                        damping: real_at_nodes(ctx, bp.material, "Sound Damping", e)?,
                        load: complex_at_nodes(ctx, bp.body_force, "Pressure Source", e, n)?,
                    }
                };
                helmholtz_bulk_local(
                    &*self.basis,
                    &*self.quad,
                    el.etype,
                    &xyz,
                    &mat,
                    &cfg,
                    &mut self.complex,
                )?;
                self.complex.to_real(&mut self.local);
                ctx.default_update(&self.local, e)?;
            }

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let bc = {
                    let p = ctx.bc_params(e);
                    let parent = el.parents.iter().flatten().next().copied().ok_or_else(|| {
                        SolverError::fatal(
                            CALLER,
                            format!("boundary element {} has no parent element", e.0),
                        )
                    })?;
                    let pmat = ctx.body_params(parent).and_then(|b| b.material);
                    let parent_prop = |key: &str| -> Result<Vec<f64>, SolverError> {
                        real_at_nodes(ctx, pmat, key, e)?.ok_or_else(|| {
                            SolverError::fatal(
                                "GetParentMatProp",
                                format!("Failed to find \"{key}\""),
                            )
                        })
                    };
                    let impedance = if p.is_some_and(|p| flag(p, "plane wave bc")) {
                        Some(
                            parent_prop("Sound Speed")?
                                .into_iter()
                                .map(|c| Complex64::new(c, 0.0))
                                .collect(),
                        )
                    } else {
                        complex_at_nodes(ctx, p, "Wave Impedance", e, n)?
                    };
                    HelmholtzBcAtNodes {
                        density: parent_prop("Density")?,
                        impedance,
                        wave_flux: complex_at_nodes(ctx, p, "Wave Flux", e, n)?,
                    }
                };
                helmholtz_boundary_local(
                    &*self.basis,
                    &*self.quad,
                    el.etype,
                    &xyz,
                    &bc,
                    &cfg,
                    &mut self.complex,
                )?;
                self.complex.to_real(&mut self.local);
                ctx.default_update(&self.local, e)?;
            }

            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        Ok(())
    }
}

/// `GetReal(section, "<base> 1")` + i `GetReal(section, "<base> 2")`: `None`
/// when neither part is given, a missing part is zero.
fn complex_at_nodes(
    ctx: &SolverCtx<'_>,
    p: Option<&ParamMap>,
    base: &str,
    e: ElemId,
    n: usize,
) -> Result<Option<Vec<Complex64>>, SolverError> {
    let re = real_at_nodes(ctx, p, &format!("{base} 1"), e)?;
    let im = real_at_nodes(ctx, p, &format!("{base} 2"), e)?;
    if re.is_none() && im.is_none() {
        return Ok(None);
    }
    let get = |v: &Option<Vec<f64>>, i: usize| v.as_ref().map_or(0.0, |v| v[i]);
    Ok(Some(
        (0..n)
            .map(|i| Complex64::new(get(&re, i), get(&im, i)))
            .collect(),
    ))
}
