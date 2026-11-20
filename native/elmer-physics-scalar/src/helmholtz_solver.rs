//! `HelmholtzSolver`: the `Procedure = "HelmholtzSolve" "HelmholtzSolver"`
//! solver, driving the kernels of [`crate::helmholtz`] through [`SolverCtx`]
//! the way `HelmholtzSolver` (HelmholtzSolve.F90:99) does.
//!
//! The pressure is a two-dof variable (`Pressure 1` real, `Pressure 2`
//! imaginary). As the Fortran's `_Init`, the system is marked `Linear
//! System Complex` (unless the case says otherwise), so `default_solve`
//! solves the interleaved storage as the complex system: the Krylov methods
//! in complex arithmetic (CG with the unconjugated product, which a complex
//! symmetric Helmholtz matrix needs) and UMFPACK on the interleaved form.

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
        if flag(bc, "flow interface") {
            return Err(unsupported(
                "Boundary Condition \"flow interface\"",
                "harmonic flow coupling is not ported",
            ));
        }
    }
    if model.boundary_conditions.iter().any(|bc| flag(bc, "structure interface")) {
        let p = ctx.params();
        if p.integer("displacement variable eigenmode").is_none_or(|m| m <= 0) {
            return Err(unsupported(
                "Boundary Condition \"structure interface\" without \"Displacement Variable EigenMode\"",
                "only the eigenmode (DispMode 1) displacement is ported",
            ));
        }
        if flag(p, "displacement variable frequency") {
            return Err(unsupported("Solver keyword \"Displacement Variable Frequency\"", "not ported"));
        }
        if axisymmetric(ctx)? {
            return Err(unsupported("a structure interface in axisymmetric coordinates", "the metric factor is not ported"));
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
        // HelmholtzSolver_Init (HelmholtzSolve.F90:51): the interleaved
        // system is solved as the complex one, so the Krylov methods run in
        // complex arithmetic (complex CG with the unconjugated product).
        if ctx.params().get("linear system complex").is_none() {
            ctx.params_mut().insert("Linear System Complex", Value::Logical(true));
        }
        if ctx.params().0.keys().any(|k| k.starts_with("linear system preconditioning damp coefficient"))
            && ctx.params().get("allocate preconditioning matrix").is_none()
        {
            ctx.params_mut().insert("Allocate Preconditioning Matrix", Value::Logical(true));
        }
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
        // HelmholtzSolve.F90:268-280: the damped preconditioning matrix of an
        // iterative solve, `ShiftCoeff = Damp Coefficient + i Damp
        // Coefficient Im`.
        let (re, im) = (
            ctx.params().real("linear system preconditioning damp coefficient"),
            ctx.params().real("linear system preconditioning damp coefficient im"),
        );
        let mut prec_shift = (re.is_some() || im.is_some()).then(|| Complex64::new(re.unwrap_or(0.0), im.unwrap_or(0.0)));
        if prec_shift.is_some()
            && ctx.params().string("linear system solver").map(str::to_ascii_lowercase).as_deref() != Some("iterative")
        {
            ctx.warn(CALLER, "Damped preconditioning makes sence only for iterative methods, canceling!");
            prec_shift = None;
        }

        let mut cfg = HelmholtzConfig::new(dim, omega);
        cfg.axisymmetric = axi;
        cfg.prec_shift = prec_shift;
        let mut prec_local = LocalSystem::zeros(0);
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
                if cfg.prec_shift.is_some() {
                    self.complex.prec_to_real(&mut prec_local);
                    ctx.default_update_prec(&prec_local, e)?;
                }
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
                if cfg.prec_shift.is_some() {
                    self.complex.prec_to_real(&mut prec_local);
                    ctx.default_update_prec(&prec_local, e)?;
                }
            }

            self.structure_interface(ctx, &bnd, omega)?;

            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        Ok(())
    }
}

impl HelmholtzSolver {
    /// HelmholtzSolve.F90:410-460 for `Structure Interface` boundaries: the
    /// wall velocity `i omega u` of the displacement eigenmode
    /// (`Displacement Variable EigenMode`, `DispMode = 1`) enters the
    /// right-hand side through `LocalInterfaceMatrix`: `FORCE(p) -= s
    /// Basis(p) i omega (n . v)`, with the checked `NormalVector`.
    fn structure_interface(&mut self, ctx: &mut SolverCtx<'_>, bnd: &[ElemId], omega: f64) -> Result<(), SolverError> {
        if !ctx.sim.model.boundary_conditions.iter().any(|bc| flag(bc, "structure interface")) {
            return Ok(());
        }
        let p = ctx.params().clone();
        let vname = p.string("displacement variable name").unwrap_or("displacement").to_string();
        let mode = p.integer("displacement variable eigenmode").unwrap_or(0) as usize;
        let (eig, dperm, ddofs) = {
            let dv = ctx
                .sim
                .vars
                .get(&vname)
                .ok_or_else(|| SolverError::fatal(CALLER, format!("No displacement variable associated:{vname}")))?;
            let eig = dv
                .eigen_vectors
                .get(mode - 1)
                .cloned()
                .ok_or_else(|| SolverError::fatal(CALLER, "Requested eigenmode does not exist"))?;
            (eig, dv.perm.as_deref().map(<[u32]>::to_vec), dv.dofs)
        };
        let dim = ctx.mesh.dim();
        if ddofs < dim {
            return Err(SolverError::fatal(CALLER, "Eigenmode displacement field should have at least 1*dim components"));
        }
        let iu = Complex64::new(0.0, 1.0);
        let mut dofs = Vec::new();
        let mut be = elmer_core::BasisEval::default();
        for &e in bnd {
            ctx.element_dofs(e, &mut dofs);
            if dofs.contains(&NO_DOF) {
                continue;
            }
            let el = ctx.mesh.element(e);
            if el.etype.family() == 1 {
                continue;
            }
            if !ctx.bc_params(e).is_some_and(|bc| flag(bc, "structure interface")) {
                continue;
            }
            let nodes = el.nodes.to_vec();
            let etype = el.etype;
            let n = nodes.len();
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            // WallVelocity = ImUnit * AngularFrequency * u.
            let mut wall = vec![[Complex64::new(0.0, 0.0); 3]; n];
            for (j, &node) in nodes.iter().enumerate() {
                let k = match &dperm {
                    Some(pm) if pm[node as usize] == NO_DOF => {
                        return Err(SolverError::fatal("HelmholtzSolve", "Displacement solution is not available on boundary"))
                    }
                    Some(pm) => pm[node as usize] as usize,
                    None => node as usize,
                };
                for l in 0..dim {
                    wall[j][l] = eig[k * ddofs + l];
                }
            }
            let iw = iu * omega;
            for w in wall.iter_mut() {
                for c in w.iter_mut() {
                    *c = iw * *c;
                }
            }
            self.complex.reset(n);
            let ips = self.quad.gauss_points(etype, None);
            for qp in &ips {
                self.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m.to_string()))?;
                let s = qp.weight * be.det_j;
                let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [qp.uvw[0], qp.uvw[1]])?;
                let mut normvelo = Complex64::new(0.0, 0.0);
                for i in 0..dim {
                    let mut sum = Complex64::new(0.0, 0.0);
                    for j in 0..n {
                        sum += wall[j][i] * be.basis[j];
                    }
                    normvelo += normal[i] * sum;
                }
                for pp in 0..n {
                    let t = s * be.basis[pp] * iu * omega * normvelo;
                    self.complex.force[pp] -= t;
                }
            }
            self.complex.to_real(&mut self.local);
            ctx.default_update(&self.local, e)?;
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
