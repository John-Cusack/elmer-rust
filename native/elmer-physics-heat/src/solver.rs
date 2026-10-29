//! `HeatSolver`: the `Procedure = "HeatSolve" "HeatSolver"` solver, driving the
//! element kernels through [`SolverCtx`] the way `HeatSolver` in
//! HeatSolve.F90:254 drives `LocalMatrixVec` and `LocalMatrixBC`.

use std::sync::OnceLock;

use elmer_core::{ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::material::{Conductivity, HeatBcAtNodes, HeatBoundaryConfig, HeatBulkConfig};
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

/// A keyword that must be a constant scalar.
fn const_real(p: &dyn Params, key: &str) -> Result<Option<f64>, SolverError> {
    match p.get(key) {
        None => Ok(None),
        Some(Value::Real {
            rows: 1,
            cols: 1,
            data,
        }) => Ok(Some(data[0])),
        Some(_) => Err(SolverError::Unsupported(format!(
            "\"{key}\" must be a constant scalar (property evaluation is not wired)"
        ))),
    }
}

/// `Heat Conductivity` by shape: 1x1 scalar, n x 1 diagonal, n x n full.
fn conductivity(p: &dyn Params, n: usize) -> Result<Option<Conductivity>, SolverError> {
    let key = "Heat Conductivity";
    match p.get(key) {
        None => Ok(None),
        Some(Value::Real { rows, cols, data }) => {
            let (r, c) = (*rows, *cols);
            if r == 1 && c == 1 {
                Ok(Some(Conductivity::Scalar(vec![data[0]; n])))
            } else if (2..=3).contains(&r) && c == 1 {
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
        Some(_) => Err(SolverError::Unsupported(format!(
            "\"{key}\" must be constant (property evaluation is not wired)"
        ))),
    }
}

fn nodal(v: Option<f64>, n: usize) -> Option<Vec<f64>> {
    v.map(|x| vec![x; n])
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
    /// convection, whose SUPG parameters need element sizes that
    /// `MeshView` does not provide.
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
                .is_some_and(|c| !c.eq_ignore_ascii_case("none"))
            {
                return Err(SolverError::Unsupported(
                    "Convection: element size hK for SUPG is not available from MeshView".into(),
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
                    let m = bp.material;
                    let bf = bp.body_force;
                    let get = |p: Option<&ParamMap>, key: &str| match p {
                        Some(p) => const_real(p, key),
                        None => Ok(None),
                    };
                    HeatMaterialAtNodes {
                        conductivity: match m {
                            Some(m) => conductivity(m, n)?,
                            None => None,
                        },
                        density: nodal(get(m, "Density")?, n),
                        heat_capacity: nodal(get(m, "Heat Capacity")?, n),
                        heat_transfer_multiplier: nodal(get(m, "Heat Transfer Multiplier")?, n),
                        heat_source: nodal(get(bf, "Heat Source")?, n),
                        volumetric_heat_source: nodal(get(bf, "Volumetric Heat Source")?, n),
                        convection_velocity: None,
                    }
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
                    let p = ctx.bc_params(e).expect("boundary_elements filters on BC");
                    let n = dofs.len();
                    HeatBcAtNodes {
                        heat_flux: nodal(const_real(p, "Heat Flux")?, n),
                        heat_transfer_coefficient: nodal(
                            const_real(p, "Heat Transfer Coefficient")?,
                            n,
                        ),
                        external_temperature: nodal(const_real(p, "External Temperature")?, n),
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
