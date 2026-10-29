//! `HeatSolver`: the `Procedure = "HeatSolve" "HeatSolver"` solver, driving the
//! element kernels through [`SolverCtx`] the way `HeatSolver` in
//! HeatSolve.F90:254 drives `LocalMatrixVec` and `LocalMatrixBC`.

use std::sync::OnceLock;

use elmer_core::{ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::material::{
    Conductivity, HeatBcAtNodes, HeatBoundaryConfig, HeatBulkConfig, Stabilization,
};
use crate::preflight::stabilized;
use crate::{
    check_supported, heat_boundary_local, heat_bulk_local, HeatError, HeatMaterialAtNodes,
    HeatSections, PROCEDURE,
};

/// Builds the element basis and quadrature the solver integrates with. The
/// registry's factories are plain `fn` pointers, so the front end (`elmer-fe`)
/// hands its implementations over through [`install_fe`] instead.
pub type FeProvider = fn() -> (Box<dyn ElementBasis + Send>, Box<dyn Quadrature + Send>);

static FE: OnceLock<FeProvider> = OnceLock::new();

/// Installs the basis and quadrature provider used by [`factory`]. Returns
/// `false` if one is already installed.
pub fn install_fe(provider: FeProvider) -> bool {
    FE.set(provider).is_ok()
}

/// Registers `("HeatSolve", "HeatSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`HeatSolver`]. Fails
/// with `NotWired` until [`install_fe`] has been called.
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    let provider = FE
        .get()
        .ok_or_else(|| SolverError::NotWired("element basis and quadrature (elmer-fe)".into()))?;
    let (basis, quad) = provider();
    Ok(Box::new(HeatSolver::new(basis, quad)))
}

impl From<HeatError> for SolverError {
    fn from(e: HeatError) -> Self {
        match e {
            HeatError::Unsupported { .. } => SolverError::Unsupported(e.to_string()),
            e => SolverError::fatal("HeatSolver", e.to_string()),
        }
    }
}

pub struct HeatSolver {
    basis: Box<dyn ElementBasis + Send>,
    quad: Box<dyn Quadrature + Send>,
    local: LocalSystem,
}

impl HeatSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis,
            quad,
            local: LocalSystem::zeros(0),
        }
    }
}

/// `ListGetElementReal` of a scalar keyword at the nodes of element `e`:
/// constants, tables, MATC and procedures through `elmer-expr`, evaluated at
/// the current iterate of the fields they depend on.
fn real_at_nodes(
    ctx: &SolverCtx<'_>,
    p: Option<&ParamMap>,
    key: &str,
    e: ElemId,
) -> Result<Option<Vec<f64>>, SolverError> {
    match p {
        Some(p) => ctx.real_at_nodes(p, key, e),
        None => Ok(None),
    }
}

/// `Heat Conductivity` by shape: 1x1 scalar, n x 1 diagonal, n x n full.
/// A scalar may depend on fields; a tensor must be constant.
fn conductivity(
    ctx: &SolverCtx<'_>,
    p: &ParamMap,
    e: ElemId,
    n: usize,
) -> Result<Option<Conductivity>, SolverError> {
    let key = "Heat Conductivity";
    match p.get(key) {
        None => Ok(None),
        Some(Value::Real { rows, cols, data }) if rows * cols > 1 => {
            let (r, c) = (*rows, *cols);
            if (2..=3).contains(&r) && c == 1 {
                let mut d = [0.0; 3];
                d[..r].copy_from_slice(&data[..r]);
                Ok(Some(Conductivity::Diagonal(vec![d; n])))
            } else if (2..=3).contains(&r) && r == c {
                let mut k = [[0.0; 3]; 3];
                for i in 0..r {
                    for j in 0..r {
                        k[i][j] = data[i * c + j];
                    }
                }
                Ok(Some(Conductivity::Full(vec![k; n])))
            } else {
                Err(SolverError::Unsupported(format!(
                    "\"{key}\" with shape {r}x{c}"
                )))
            }
        }
        Some(Value::Dependent(d)) if d.rows * d.cols > 1 => Err(SolverError::Unsupported(format!(
            "\"{key}\": a non-constant tensor"
        ))),
        Some(_) => Ok(ctx.real_at_nodes(p, key, e)?.map(Conductivity::Scalar)),
    }
}

/// `Convection = Constant`: `Convection Velocity i` of the material at the
/// element nodes; a component without its own keyword is zero.
fn convection_velocity(
    ctx: &SolverCtx<'_>,
    m: Option<&ParamMap>,
    e: ElemId,
    n: usize,
) -> Result<Vec<[f64; 3]>, SolverError> {
    let mut v = vec![[0.0; 3]; n];
    for i in 0..3 {
        if let Some(c) = real_at_nodes(ctx, m, &format!("Convection Velocity {}", i + 1), e)? {
            for (vn, ci) in v.iter_mut().zip(c) {
                vn[i] = ci;
            }
        }
    }
    Ok(v)
}

fn axisymmetric(ctx: &SolverCtx<'_>) -> Result<bool, SolverError> {
    match ctx.simulation_params().string("coordinate system") {
        None => Ok(false),
        Some(s) => {
            let s = elmer_core::normalize_key(s);
            if s.starts_with("cartesian") {
                Ok(false)
            } else if s == "axi symmetric" || s == "cylindric symmetric" {
                Ok(true)
            } else {
                Err(SolverError::Unsupported(format!("Coordinate System = {s}")))
            }
        }
    }
}

impl Solver for HeatSolver {
    /// Refuses every keyword in the sections this solver reads that the
    /// kernels would otherwise ignore (see [`check_supported`]), and
    /// `Convection = Computed`, which needs the flow solution field.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        axisymmetric(ctx)?;
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut bodies: Vec<usize> = Vec::new();
        for e in elems {
            if let Some(b) = ctx.body_params(e) {
                if !bodies.contains(&b.body_id) {
                    bodies.push(b.body_id);
                }
            }
        }
        let model = &ctx.sim.model;
        let bcs: Vec<&dyn Params> = model
            .boundary_conditions
            .iter()
            .map(|p| p as &dyn Params)
            .collect();
        for b in bodies {
            let Some(bp) = model.body_params(b) else {
                continue;
            };
            let sections = HeatSections {
                solver: ctx.params(),
                equation: bp.equation.map(|p| p as &dyn Params),
                material: bp.material.map(|p| p as &dyn Params),
                body_force: bp.body_force.map(|p| p as &dyn Params),
                boundaries: bcs.clone(),
            };
            check_supported(&sections)?;
            if bp
                .equation
                .and_then(|e| e.string("Convection"))
                .is_some_and(|c| c.eq_ignore_ascii_case("computed"))
            {
                return Err(SolverError::Unsupported(
                    "Convection = Computed: reading the flow solution as the convection \
                     field is not implemented"
                        .into(),
                ));
            }
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        dt: f64,
        transient: bool,
    ) -> Result<(), SolverError> {
        const CALLER: &str = "HeatSolver";
        let axisym = axisymmetric(ctx)?;
        let dim = ctx.mesh.dim();
        let maxiter = ctx
            .params()
            .integer("nonlinear system max iterations")
            .map_or(1, |n| n.max(1));
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();

        ctx.info(
            CALLER,
            "------------------------------------------------",
            4,
        );
        ctx.info(CALLER, "Solving energy equation for temperature", 4);

        let mut bulk_cfg = HeatBulkConfig::new(dim);
        bulk_cfg.axisymmetric = axisym;
        bulk_cfg.transient = transient;
        let bnd_cfg = HeatBoundaryConfig {
            axisymmetric: axisym,
            quad_points: None,
        };
        let mut dofs = Vec::new();
        // SUPG (`Stabilize`) takes `Element % hK` and `Element %
        // StabilizationMK`, which `MeshStabParams` fills at mesh load with
        // `StabParam` in the model dimension.
        let supg = stabilized(ctx.params());
        let long_edge = ctx
            .simulation_params()
            .logical("Stabilization Use Longest Element Edge")
            .unwrap_or(false);

        for iter in 1..=maxiter {
            ctx.info(CALLER, &format!("Heat solver iteration: {iter}"), 4);
            ctx.default_initialize()?;

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let mat = {
                    let bp = ctx.body_params(e).ok_or_else(|| {
                        SolverError::fatal(CALLER, format!("element {} has no body", e.0))
                    })?;
                    let (m, bf) = (bp.material, bp.body_force);
                    let convected = bp
                        .equation
                        .and_then(|q| q.string("Convection"))
                        .is_some_and(|c| c.eq_ignore_ascii_case("constant"));
                    HeatMaterialAtNodes {
                        conductivity: match m {
                            Some(m) => conductivity(ctx, m, e, n)?,
                            None => None,
                        },
                        density: real_at_nodes(ctx, m, "Density", e)?,
                        heat_capacity: real_at_nodes(ctx, m, "Heat Capacity", e)?,
                        heat_transfer_multiplier: real_at_nodes(
                            ctx,
                            m,
                            "Heat Transfer Multiplier",
                            e,
                        )?,
                        heat_source: real_at_nodes(ctx, bf, "Heat Source", e)?,
                        volumetric_heat_source: real_at_nodes(
                            ctx,
                            bf,
                            "Volumetric Heat Source",
                            e,
                        )?,
                        convection_velocity: if convected {
                            Some(convection_velocity(ctx, m, e, n)?)
                        } else {
                            None
                        },
                    }
                };
                bulk_cfg.stabilization = if supg && mat.convection_velocity.is_some() {
                    let (m_k, h_k) = elmer_fe::stab_param(el.etype, &xyz, dim, long_edge)
                        .map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
                    Stabilization::Supg { h_k, m_k }
                } else {
                    Stabilization::None
                };
                heat_bulk_local(
                    &*self.basis,
                    &*self.quad,
                    el.etype,
                    &xyz,
                    &mat,
                    &bulk_cfg,
                    &mut self.local,
                )?;
                if transient {
                    ctx.default_first_order_time(&mut self.local, e, dt)?;
                }
                ctx.default_update(&self.local, e)?;
            }

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let bc = {
                    let p = ctx.bc_params(e);
                    HeatBcAtNodes {
                        heat_flux: real_at_nodes(ctx, p, "Heat Flux", e)?,
                        heat_transfer_coefficient: real_at_nodes(
                            ctx,
                            p,
                            "Heat Transfer Coefficient",
                            e,
                        )?,
                        external_temperature: real_at_nodes(ctx, p, "External Temperature", e)?,
                    }
                };
                if bc.heat_flux.is_none() && bc.heat_transfer_coefficient.is_none() {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                heat_boundary_local(
                    &*self.basis,
                    &*self.quad,
                    el.etype,
                    &xyz,
                    &bc,
                    &bnd_cfg,
                    &mut self.local,
                )?;
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
