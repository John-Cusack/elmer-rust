//! `FlowSolver`: the `Procedure = "FlowSolve" "FlowSolver"` solver, driving
//! the element kernels through [`SolverCtx`] the way `FlowSolver`
//! (FlowSolve.F90:117) drives `NavierStokesCompose` and
//! `NavierStokesBoundary`, with its Picard/Newton switching.

use std::sync::OnceLock;

use elmer_core::{ElemId, ElementBasis, LocalSystem, MeshView, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::sim::parse_variable_spec;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::boundary::{element_center, flow_boundary_local, FlowBoundaryInput};
use crate::bulk::{flow_bulk_local, Discretization, FlowBulkInput, PowerLaw};
use crate::error::FlowError;
use crate::bulk::StabParams;
use crate::cylindric::{cylindric_bulk_local, ns_condensate, CylindricBulkInput};

const CALLER: &str = "FlowSolver";

/// Builds the element basis and quadrature the solver integrates with, in
/// place of the `elmer-fe` defaults ([`elmer_fe::LagrangeBasis`] for the
/// mesh dimension and [`elmer_fe::ElmerQuadrature`]).
pub type FeProvider = fn() -> (Box<dyn ElementBasis + Send>, Box<dyn Quadrature + Send>);

static FE: OnceLock<FeProvider> = OnceLock::new();

/// Installs a basis and quadrature provider for [`factory`]-built solvers.
/// Returns `false` if one is already installed.
pub fn install_fe(provider: FeProvider) -> bool {
    FE.set(provider).is_ok()
}

/// Registers `("FlowSolve", "FlowSolver")` and `("FixTangentVelo", "FixTangentVelo")`.
pub fn register(registry: &mut SolverRegistry) {
    let (lib, entry) = crate::REGISTRY_ENTRY;
    registry.register(lib, entry, factory);
    registry.register("Stokes", "StokesSolver", crate::stokes::factory);
    registry.register("StreamSolver", "StreamSolver", crate::stream::factory);
    registry.register("FixTangentVelo", "FixTangentVelo", crate::fix_tangent::factory);
    registry.register("AdvectionReaction", "AdvectionReactionSolver", crate::adv_react::factory);
    registry.register("PoissonDG", "PoissonSolver", crate::poisson_dg::factory);
    registry.register("ParStokes", "StokesSolver", crate::parstokes::factory);
    registry.register("FreeSurfaceSolver", "FreeSurfaceSolver", crate::fs_solver::factory);
    registry.register("ReynoldsSolver", "ReynoldsSolver", crate::reynolds::solver_factory);
    registry.register("ReynoldsSolver", "ReynoldsPostprocess", crate::reynolds::post_factory);
    registry.register("VelocityPrecond", "VelocityPrecond", crate::parstokes::velocity_precond_factory);
    registry.register("PressurePrecond", "PressurePrecond", crate::parstokes::pressure_precond_factory);
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`FlowSolver`]: the
/// installed provider's basis and quadrature, else the `elmer-fe` ones.
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(match FE.get() {
        Some(provider) => {
            let (basis, quad) = provider();
            FlowSolver::new(basis, quad)
        }
        None => FlowSolver::default(),
    }))
}

impl From<FlowError> for SolverError {
    fn from(e: FlowError) -> Self {
        match e {
            FlowError::Unsupported(_) => SolverError::Unsupported(format!("{CALLER}: {e}")),
            e => SolverError::fatal(CALLER, e.to_string()),
        }
    }
}

fn unsupported(what: impl std::fmt::Display) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

fn flag(p: Option<&ParamMap>, key: &str) -> bool {
    p.and_then(|p| p.logical(key)).unwrap_or(false)
}

/// The incompressible stabilized P1/P1 solver.
pub struct FlowSolver {
    /// `None` until `init` picks the `elmer-fe` defaults.
    basis: Option<Box<dyn ElementBasis + Send>>,
    quad: Option<Box<dyn Quadrature + Send>>,
    local: LocalSystem,
    /// FlowSolve's `NewtonLinearization` and `RecheckNewton`: initialized
    /// locals, so Fortran keeps them across calls (and, unlike here, across
    /// every FlowSolve instance).
    newton: bool,
    recheck_newton: bool,
    /// FlowSolve's SAVEd `ComputeFree`: free surfaces move from the first
    /// iteration past `Free Surface After Iterations` / `Tolerance` on.
    compute_free: bool,
}

impl Default for FlowSolver {
    /// A solver on the `elmer-fe` basis and Gauss rules.
    fn default() -> Self {
        Self {
            basis: None,
            quad: None,
            local: LocalSystem::zeros(0),
            newton: false,
            recheck_newton: false,
            compute_free: false,
        }
    }
}

impl FlowSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis: Some(basis),
            quad: Some(quad),
            ..Self::default()
        }
    }
}

/// The `Stabilization Method` FlowSolver resolves (FlowSolve.F90:505-547).
/// Only `stabilized` is ported.
/// `StabilizeFlag` of FlowSolve.F90:509-524.
fn stabilization_method(p: &ParamMap) -> String {
    match p.string("stabilization method") {
        Some(m) => elmer_core::normalize_key(m),
        None if flag(Some(p), "stabilize") => "stabilized".into(),
        None if flag(Some(p), "bubbles") || p.get("element").is_some() => "bubbles".into(),
        None => "stabilized".into(),
    }
}

fn check_stabilization(p: &ParamMap) -> Result<(), SolverError> {
    let method = stabilization_method(p);
    match method.as_str() {
        "stabilized" => {
            if flag(Some(p), "bubbles") {
                return Err(SolverError::fatal(CALLER, "You cant have stabilization and bubbles both!"));
            }
            Ok(())
        }
        // NavierStokesCylindricalCompose's p bubbles; the solver checks the
        // coordinates and the elements.
        "pbubbles" => Ok(()),
        // Legacy bubbles; the solver checks the coordinates and elements.
        "bubbles" => Ok(()),
        "p2/p1" | "p2p1" | "vms" => {
            Err(unsupported(format!("Stabilization Method = {method}")))
        }
        other => Err(SolverError::fatal(CALLER, format!("Unknown \"stabilization method\": {other}"))),
    }
}

/// Solver-section options this crate does not port.
fn check_solver_params(p: &ParamMap) -> Result<(), SolverError> {
    check_stabilization(p)?;
    if p.get("solver coordinate system").is_some() {
        return Err(unsupported("Solver Coordinate System"));
    }
    for (key, what) in [
        ("out of plane flow", "out of plane flow"),
        ("adaptive mesh refinement", "adaptive mesh refinement"),
        ("nonlinear system linesearch", "nonlinear line search"),
    ] {
        if flag(Some(p), key) {
            return Err(unsupported(what));
        }
    }
    if p.get("implicit friction direction vector").is_some() {
        return Err(unsupported("implicit friction"));
    }
    Ok(())
}

/// Material, body force and equation options this crate does not port.
fn check_body(eq: Option<&ParamMap>, mat: Option<&ParamMap>, bf: Option<&ParamMap>) -> Result<(), SolverError> {
    if let Some(m) = mat {
        if let Some(c) = m.string("compressibility model") {
            let c = elmer_core::normalize_key(c);
            if matches!(
                c.as_str(),
                "perfect gas" | "perfect gas equation 1" | "user defined" | "pressure dependent"
                    | "artificial compressible"
            ) {
                return Err(unsupported(format!("Compressibility Model = {c}")));
            }
        }
        if m.get("viscosity model").is_some()
            && m.string("viscosity model").map(elmer_core::normalize_key).as_deref() != Some("power law")
        {
            return Err(unsupported("Viscosity Model other than power law (non-newtonian or anisotropic viscosity)"));
        }
    }
    if let Some(b) = bf {
        for (key, what) in [
            ("lorentz force", "Lorentz force"),
            ("boussinesq", "Boussinesq body force"),
            ("potential force", "potential force"),
        ] {
            if flag(Some(b), key) {
                return Err(unsupported(what));
            }
        }
    }
    let _ = eq;
    Ok(())
}

/// FlowSolve.F90:1279-1302: `Alpha`, the `Surface Tension Coefficient`, and
/// `Beta`, zero without a `Temperature` variable, else the coefficient times
/// `1 - Surface Tension Expansion Coefficient * Temperature` (the expansion
/// coefficient alone at a node outside the temperature's perm), or the
/// coefficient when no expansion coefficient is given.
fn surface_tension(ctx: &SolverCtx<'_>, bc: &ParamMap, e: ElemId, n: usize) -> Result<(Vec<f64>, Vec<f64>), SolverError> {
    let coeff = ctx.real_at_nodes(bc, "surface tension coefficient", e)?;
    let alpha = coeff.clone().unwrap_or_else(|| vec![0.0; n]);
    let Some(temp) = ctx.variable("temperature") else {
        return Ok((alpha, vec![0.0; n]));
    };
    let beta = match ctx.real_at_nodes(bc, "surface tension expansion coefficient", e)? {
        Some(mut b) => {
            for (j, &node) in ctx.mesh.element(e).nodes.iter().enumerate() {
                let k = temp.perm.map_or(node, |p| p[node as usize]);
                if k != NO_DOF {
                    b[j] = 1.0 - b[j] * temp.values[k as usize * temp.dofs];
                }
            }
            let c = coeff.unwrap_or_else(|| vec![0.0; n]);
            b.iter().zip(&c).map(|(b, c)| b * c).collect()
        }
        None => alpha.clone(),
    };
    Ok((alpha, beta))
}

/// FlowSolve.F90:916-932: `Density * (1 - Heat Expansion Coefficient *
/// (T - Reference Temperature))`, `T` the `Temperature` at the nodes when
/// all of them have one, else zero.
fn thermal_density(ctx: &SolverCtx<'_>, mat: Option<&ParamMap>, e: ElemId, n: usize, density: &mut [f64]) -> Result<(), SolverError> {
    let beta = nodal(ctx, mat, "heat expansion coefficient", e, n)?;
    let t0 = nodal(ctx, mat, "reference temperature", e, n)?;
    let mut temp = vec![0.0; n];
    if let Some(v) = ctx.variable("temperature") {
        let nodes = ctx.mesh.element(e).nodes;
        let slot = |k: u32| v.perm.map_or(Some(k), |p| (p[k as usize] != NO_DOF).then(|| p[k as usize]));
        if nodes.iter().all(|&k| slot(k).is_some()) {
            for (t, &k) in temp.iter_mut().zip(nodes) {
                *t = v.values[slot(k).expect("checked") as usize * v.dofs];
            }
        }
    }
    for k in 0..n {
        density[k] *= 1.0 - beta[k] * (temp[k] - t0[k]);
    }
    Ok(())
}

/// FlowSolve.F90:940-958: with `Porous Media`, `Porous Resistivity` at the
/// nodes: a constant array (one value for every component, else one per
/// component), or `Porous Resistivity 1..dim`.
fn porous_resistivity(
    ctx: &SolverCtx<'_>,
    mat: Option<&ParamMap>,
    e: ElemId,
    n: usize,
    vdim: usize,
) -> Result<Option<Vec<[f64; 3]>>, SolverError> {
    let Some(m) = mat.filter(|m| flag(Some(m), "porous media")) else {
        return Ok(None);
    };
    let mut drag = vec![[0.0; 3]; n];
    match m.get("porous resistivity") {
        Some(Value::Real { rows, cols: 1, data }) => {
            for d in drag.iter_mut() {
                for (i, x) in d.iter_mut().enumerate().take(vdim) {
                    *x = if *rows == 1 { data[0] } else if i < *rows { data[i] } else { 0.0 };
                }
            }
        }
        Some(_) => return Err(unsupported("a non-constant Porous Resistivity array")),
        None => {
            for i in 0..vdim {
                let v = nodal(ctx, Some(m), &format!("porous resistivity {}", i + 1), e, n)?;
                for (d, x) in drag.iter_mut().zip(v) {
                    d[i] = x;
                }
            }
        }
    }
    Ok(Some(drag))
}

/// Boundary-condition options this crate does not port.
fn check_bc(bc: &ParamMap, _var_name: &str) -> Result<(), SolverError> {
    for (key, what) in [
        ("wall law", "wall law"),
        ("vms wall", "VMS wall"),
    ] {
        if flag(Some(bc), key) {
            return Err(unsupported(what));
        }
    }
    Ok(())
}

/// `Mesh % Stabilize` (MeshStabParams, MeshBasics.F90:2338): only when some
/// solver asks for stabilization by keyword does Elmer fill
/// `Element % StabilizationMK`; otherwise mK stays zero and FlowSolve's
/// default "stabilized" method assembles with Tau = Delta = 0.
fn mesh_stabilize(ctx: &SolverCtx<'_>) -> bool {
    ctx.sim.model.solvers.iter().any(|p| {
        let method = p.string("stabilization method").map(elmer_core::normalize_key);
        flag(Some(p), "stabilize")
            || matches!(method.as_deref(), Some("vms") | Some("stabilized"))
            || flag(Some(p), "pressure stabilization")
    })
}

/// The coordinate system FlowSolve assembles in: `false` for cartesian,
/// `true` for `Cylindric Symmetric` (NavierStokesCylindricalCompose, a
/// swirl velocity as the third component).
fn cylindric(ctx: &SolverCtx<'_>) -> Result<bool, SolverError> {
    match ctx.simulation_params().string("coordinate system") {
        None => Ok(false),
        Some(s) if elmer_core::normalize_key(s).starts_with("cartesian") => Ok(false),
        Some(s) if elmer_core::normalize_key(s) == "cylindric symmetric" => Ok(true),
        Some(s) => Err(unsupported(format!("Coordinate System = {s}"))),
    }
}

/// The component names of the solver variable (`velocity 1`, ...,
/// `pressure`), as DefaultDirichletBCs looks them up with `ComponentName`.
fn component_names(ctx: &SolverCtx<'_>) -> Result<Vec<String>, SolverError> {
    let spec = ctx
        .params()
        .string("variable")
        .ok_or_else(|| SolverError::fatal(CALLER, "solver has no 'Variable'"))?;
    Ok(parse_variable_spec(spec, ctx.params().integer("variable dofs"))?.components)
}

fn slot_params_mut<'b>(ctx: &'b mut SolverCtx<'_>) -> &'b mut ParamMap {
    let id = ctx.solver_id();
    &mut ctx
        .sim
        .slots
        .iter_mut()
        .find(|s| s.id == id)
        .expect("the running solver has a slot")
        .params
}

/// Parent of a boundary element whose centre orients the outward normal
/// (`CheckNormalDirection`, ElementGeometry.F90): the only parent, or, at a
/// body-body boundary, the one with the higher body id.
fn normal_parent(ctx: &SolverCtx<'_>, e: ElemId) -> Option<ElemId> {
    let [left, right] = ctx.mesh.element(e).parents;
    match (left, right) {
        (Some(l), Some(r)) => {
            let body = |x: ElemId| ctx.sim.model.body_of_tag(ctx.mesh.element(x).tag);
            match body(l).cmp(&body(r)) {
                std::cmp::Ordering::Greater => Some(l),
                std::cmp::Ordering::Less => Some(r),
                std::cmp::Ordering::Equal => None,
            }
        }
        (Some(l), None) => Some(l),
        (None, r) => r,
    }
}

fn node_xyz(mesh: &dyn MeshView, nodes: &[u32]) -> Vec<[f64; 3]> {
    nodes.iter().map(|&n| mesh.coord(n)).collect()
}

/// Topological dimension of an element (`GetElementDim`).
fn element_dim(family: u16) -> usize {
    match family {
        1 => 0,
        2 => 1,
        3 | 4 => 2,
        _ => 3,
    }
}

/// FlowSolve's noslip wall code (FlowSolve.F90:1416-1433): after the
/// boundary assembly a `Noslip Wall BC` boundary gets `Velocity i = 0`
/// (the velocity components of the flow variable) in its section, which
/// `DefaultDirichletBCs` then applies like any other Dirichlet value.
fn noslip_walls(ctx: &mut SolverCtx<'_>, comps: &[String], dim: usize) {
    for bc in ctx.sim.model.boundary_conditions.iter_mut() {
        if flag(Some(bc), "noslip wall bc") {
            for key in comps.iter().take(dim) {
                bc.insert(key, Value::real(0.0));
            }
        }
    }
}

/// `AddEquationBasics` (MainUtils.F90:1241, 1827) creates the variable of
/// an `Equation = Navier-Stokes` solver filled with 1e-6, which initial
/// conditions then overwrite. elmer-solver creates it zeroed, so the values
/// no `Initial Condition` sets (the same element walk as its
/// `apply_initial_conditions`) get Elmer's start value here.
fn navier_stokes_init_value(ctx: &mut SolverCtx<'_>, comps: &[String]) -> Result<(), SolverError> {
    const INIT_VALUE: f64 = 1.0e-6;
    let is_ns = ctx
        .params()
        .string("equation")
        .is_some_and(|e| elmer_core::normalize_key(e) == "navier-stokes");
    if !is_ns {
        return Ok(());
    }
    let dofs = comps.len();
    let n = ctx.primary().map_or(0, |f| f.values.len());
    let perm = ctx.primary().and_then(|f| f.perm.map(|p| p.to_vec()));
    let mut set = vec![false; n];
    let mut slots = Vec::new();
    let elems: Vec<ElemId> = ctx.active_elements().collect();
    for e in elems {
        let Some(ic) = ctx.body_params(e).and_then(|b| b.initial_condition) else { continue };
        // Initial conditions set the nodes only; the edge, face and bubble
        // dofs of a p element keep the start value.
        slots.clear();
        let nodes = ctx.mesh.element(e).nodes;
        match &perm {
            Some(p) => slots.extend(nodes.iter().map(|&nn| p.get(nn as usize).copied().unwrap_or(NO_DOF))),
            None => slots.extend_from_slice(nodes),
        }
        for (k, key) in comps.iter().enumerate() {
            if ic.get(key).is_some() {
                for &s in &slots {
                    if s != NO_DOF {
                        set[s as usize * dofs + k] = true;
                    }
                }
            }
        }
    }
    let var = ctx
        .primary_mut()
        .ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    for (x, done) in var.values.iter_mut().zip(set) {
        if !done {
            *x = INIT_VALUE;
        }
    }
    Ok(())
}


fn nodal(ctx: &SolverCtx<'_>, p: Option<&ParamMap>, key: &str, e: ElemId, n: usize) -> Result<Vec<f64>, SolverError> {
    Ok(match p {
        Some(p) => ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0; n]),
        None => vec![0.0; n],
    })
}

/// `Gravity` from the Constants section as direction times magnitude, or
/// FlowSolve's default (0, -9.81, 0).
fn gravity(ctx: &SolverCtx<'_>) -> [f64; 3] {
    match ctx.constants().get("gravity") {
        Some(Value::Real { data, .. }) if data.len() >= 4 => [data[0] * data[3], data[1] * data[3], data[2] * data[3]],
        _ => [0.0, -9.81, 0.0],
    }
}

impl Solver for FlowSolver {
    /// Refuses up front every option of the sections this solver reads that
    /// the kernels would otherwise ignore.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let cyl = cylindric(ctx)?;
        check_solver_params(ctx.params())?;
        let dim = ctx.mesh.dim();
        let dofs = ctx.primary().map(|f| f.dofs);
        if dofs != Some(dim + 1 + usize::from(cyl)) {
            return Err(unsupported(format!(
                "a flow variable with {dofs:?} dofs in {dim}D (expected velocity components and pressure)"
            )));
        }
        if ctx.variable("mesh velocity").is_some() {
            return Err(unsupported("Mesh Velocity"));
        }
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut bodies = Vec::new();
        for e in elems {
            let b = ctx
                .body_params(e)
                .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
            if !bodies.contains(&b.body_id) {
                bodies.push(b.body_id);
                if b.material.is_none() {
                    return Err(SolverError::fatal(CALLER, format!("body {} has no material", b.body_id)));
                }
                check_body(b.equation, b.material, b.body_force)?;
            }
        }
        let var_name = ctx.variable_name().unwrap_or_default().to_string();
        for bc in &ctx.sim.model.boundary_conditions {
            check_bc(bc, &var_name)?;
        }
        let comps = component_names(ctx)?;
        navier_stokes_init_value(ctx, &comps)?;
        if self.basis.is_none() {
            self.basis = Some(Box::new(elmer_fe::LagrangeBasis::new(dim)));
        }
        if self.quad.is_none() {
            self.quad = Some(Box::new(elmer_fe::ElmerQuadrature));
        }
        // FlowSolver_init.
        slot_params_mut(ctx).insert("time derivative order", Value::Integer(vec![1]));
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        let cyl = cylindric(ctx)?;
        check_solver_params(ctx.params())?;
        let p = ctx.params().clone();
        let (convect_model, transient) = match p.string("flow model").map(elmer_core::normalize_key).as_deref() {
            Some("no convection") => (false, transient),
            Some("stokes") => (false, false),
            _ => (true, transient),
        };
        let dim = ctx.mesh.dim();
        let comps = component_names(ctx)?;
        let nsdofs = comps.len();
        if nsdofs != dim + 1 + usize::from(cyl) {
            return Err(unsupported(format!("{nsdofs} flow dofs in {dim}D")));
        }
        // The velocity components: `dim`, and the swirl in cylindric symmetry.
        let vdim = nsdofs - 1;
        if !cyl && stabilization_method(&p) == "pbubbles" {
            return Err(unsupported("Stabilization Method = pbubbles in cartesian coordinates"));
        }
        if cyl {
            if stabilization_method(&p) != "pbubbles" {
                return Err(unsupported(format!(
                    "Cylindric Symmetric FlowSolve with Stabilization Method = {} (only pbubbles is ported)",
                    stabilization_method(&p)
                )));
            }
            if transient {
                return Err(unsupported("transient Cylindric Symmetric FlowSolve (NSCondensateTransientH)"));
            }
            if ctx.sim.model.boundary_conditions.iter().any(|bc| bc.get("implicit friction coefficient").is_some()) {
                return Err(unsupported("Implicit Friction Coefficient in Cylindric Symmetric coordinates"));
            }
        }

        let disc = Discretization {
            div: flag(Some(&p), "div discretization"),
            gradp: flag(Some(&p), "gradp discretization"),
            laplace: flag(Some(&p), "laplace discretization"),
        };
        let h_scale = p.real("h scale").unwrap_or(1.0);
        let stabilize_flag = stabilization_method(&p);
        let long_edge = flag(Some(ctx.simulation_params()), "stabilization use longest element edge");
        let stabilize_mesh = mesh_stabilize(ctx);
        let nonlinear_tol = p.real("nonlinear system convergence tolerance").unwrap_or(0.0);
        let newton_tol = p.real("nonlinear system newton after tolerance").unwrap_or(0.0);
        let newton_ubound = p.real("nonlinear system newton max tolerance");
        if newton_ubound.is_some() {
            self.recheck_newton = true;
        }
        let newton_ubound = newton_ubound.unwrap_or(0.0);
        let newton_iter = p.integer("nonlinear system newton after iterations").unwrap_or(0);
        if newton_iter == 0 {
            self.newton = true;
        }
        let newton_max_iter = p.integer("nonlinear system newton max iterations");
        self.recheck_newton |= newton_max_iter.is_some();
        let newton_max_iter = newton_max_iter.unwrap_or(0);
        if flag(Some(&p), "nonlinear system reset newton") {
            self.newton = false;
        }
        // ListGetInteger without a default: an absent keyword is 0 iterations.
        let nonlinear_iter = p.integer("nonlinear system max iterations").unwrap_or(0);
        let max_norm = p.real("nonlinear system max norm return").unwrap_or(f64::MAX);
        let free_tol = p.real("free surface after tolerance").unwrap_or(f64::MAX);
        let free_iter = p.integer("free surface after iterations").unwrap_or(0);
        let free_surface = crate::free_surface::any_free_surface(ctx);

        // "Enforcing relative pressure relaxation" when no material of the
        // model names a compressibility model other than incompressible.
        let all_incompressible = ctx.sim.model.materials.iter().all(|m| {
            m.string("compressibility model")
                .is_none_or(|c| elmer_core::normalize_key(c) == "incompressible")
        });
        {
            let sp = slot_params_mut(ctx);
            if sp.get("nonlinear system norm dofs").is_none() {
                sp.insert("nonlinear system norm dofs", Value::Integer(vec![nsdofs as i32 - 1]));
            }
            if all_incompressible && sp.get("relative pressure relaxation").is_none() {
                sp.insert("relative pressure relaxation", Value::Logical(true));
            }
        }

        let grav = gravity(ctx);
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let var_name = ctx.variable_name().unwrap_or_default().to_string();
        let mut dofs = Vec::new();
        let (Some(basis), Some(quad)) = (self.basis.as_deref(), self.quad.as_deref()) else {
            return Err(SolverError::fatal(CALLER, "solve before init"));
        };

        for iter in 1..=nonlinear_iter {
            // NavierStokesCompose: a positive `Newtonian Viscosity Condition`
            // turns a non-newtonian viscosity model off for this assembly.
            let newtonian = ctx.creal(&p, "newtonian viscosity condition")?.is_some_and(|c| c > 0.0);
            ctx.info("FlowSolve", "-------------------------------------", 4);
            ctx.info("FlowSolve", &format!(" NAVIER-STOKES ITERATION {iter:>11}"), 4);
            ctx.info("FlowSolve", "-------------------------------------", 4);
            ctx.default_initialize()?;

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let xyz = node_xyz(ctx.mesh, el.nodes);
                let bp = ctx
                    .body_params(e)
                    .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
                let (eq, mat, bf) = (bp.equation, bp.material, bp.body_force);
                check_body(eq, mat, bf)?;
                let convect = convect_model && eq.and_then(|q| q.logical("ns convect")).unwrap_or(true);

                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    return Err(SolverError::fatal(CALLER, format!("active element {} has nodes outside the system", e.0)));
                }
                // A p element: NavierStokesCompose(..., Element, n, ...) takes
                // its n vertex functions, on its own rule; the other dofs get
                // no bulk contribution.
                let pe = match ctx.solver_dofs() {
                    Some(d) if !cyl => Some(p_parts(d, e, el.etype, el.nodes, &stabilize_flag)?),
                    _ => None,
                };
                let values = ctx.primary().expect("checked in init").values;
                let velocity: Vec<[f64; 3]> = (0..n)
                    .map(|k| std::array::from_fn(|i| if i < vdim { values[dofs[k * nsdofs + i] as usize] } else { 0.0 }))
                    .collect();

                let thermal = mat
                    .and_then(|m| m.string("compressibility model"))
                    .is_some_and(|c| elmer_core::normalize_key(c) == "thermal");
                let mut density = nodal(ctx, mat, "density", e, n)?;
                if thermal {
                    if cyl {
                        return Err(unsupported("Compressibility Model = thermal in Cylindric Symmetric coordinates"));
                    }
                    if transient {
                        return Err(unsupported("transient Compressibility Model = thermal (PrevDensity)"));
                    }
                    thermal_density(ctx, mat, e, n, &mut density)?;
                }
                let viscosity = nodal(ctx, mat, "viscosity", e, n)?;
                let porous = porous_resistivity(ctx, mat, e, n, vdim)?;
                if porous.is_some() && cyl {
                    return Err(unsupported("porous media in Cylindric Symmetric coordinates"));
                }
                // NavierStokesCompose takes bubbles on a linear element of a
                // compressible model whatever the method (and FlowSolve.F90:1128
                // condenses them).
                let bubbles = !cyl && (stabilize_flag == "bubbles" || thermal);
                if bubbles {
                    if pe.is_some() {
                        return Err(unsupported("Bubbles on p elements (element bubbles)"));
                    }
                    if transient {
                        return Err(unsupported("transient legacy bubbles (NSCondensateTransientH)"));
                    }
                }
                let power_law = mat.is_some_and(|m| m.get("viscosity model").is_some()) && !newtonian;
                let (exponent, critical, nominal) = if power_law {
                    if cyl {
                        return Err(unsupported("Viscosity Model = power law in Cylindric Symmetric coordinates"));
                    }
                    if pe.is_some() {
                        return Err(unsupported("Viscosity Model = power law on p elements"));
                    }
                    let m = mat.expect("checked above");
                    let opt = |key: &str| ctx.real_at_nodes(m, key, e);
                    (opt("viscosity exponent")?.unwrap_or_else(|| vec![0.0; n]), opt("critical shear rate")?, opt("nominal shear rate")?)
                } else {
                    (Vec::new(), None, None)
                };
                let mut load = vec![[0.0; 4]; n];
                if let Some(b) = bf {
                    let hydrostatic = b
                        .logical("hydrostatic pressure")
                        .or_else(|| eq.and_then(|q| q.logical("hydrostatic pressure")))
                        .unwrap_or(false);
                    if hydrostatic {
                        for l in load.iter_mut() {
                            l[..vdim].copy_from_slice(&grav[..vdim]);
                        }
                    }
                    for i in 0..vdim {
                        if let Some(f) = ctx.real_at_nodes(b, &format!("flow bodyforce {}", i + 1), e)? {
                            for (l, v) in load.iter_mut().zip(f) {
                                l[i] += v;
                            }
                        }
                    }
                }

                if cyl {
                    if angular_velocity(bf)?.is_some() {
                        return Err(SolverError::fatal(
                            "FlowSolve",
                            "Rotating coordinate implemented only for cartesian coordinates",
                        ));
                    }
                    let d = ctx
                        .solver_dofs()
                        .ok_or_else(|| unsupported("Cylindric Symmetric pbubbles without p elements (legacy bubbles)"))?;
                    let nb = d.local_bubbles(e);
                    if nb == 0 {
                        return Err(unsupported("Cylindric Symmetric pbubbles on an element without local bubbles"));
                    }
                    if dofs.len() != n * nsdofs || d.element_nodes(e) != el.nodes {
                        return Err(unsupported("Cylindric Symmetric pbubbles on p elements with edge or face dofs"));
                    }
                    let pbasis = d.element_basis(e).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
                    let inp = CylindricBulkInput {
                        etype: el.etype,
                        xyz: &xyz,
                        density: &density,
                        viscosity: &viscosity,
                        velocity: &velocity,
                        load: &load,
                        convect,
                        newton: self.newton,
                        gradp: disc.gradp,
                        bubbles: nb,
                    };
                    cylindric_bulk_local(&*pbasis, &rule, &inp, &mut self.local)?;
                    ns_condensate(n, nb, vdim, &mut self.local)?;
                    ctx.default_update(&self.local, e)?;
                    continue;
                }

                // MeshStabParams: StabParam with Mesh % Stabilize, else only
                // ElementDiameter (mK stays zero).
                let (mk, hk) = elmer_fe::stab_param(el.etype, &xyz, dim, long_edge)
                    .map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
                let stab = StabParams {
                    hk: hk * h_scale,
                    mk: if stabilize_mesh { mk } else { 0.0 },
                };
                let inp = FlowBulkInput {
                    etype: el.etype,
                    xyz: &xyz,
                    dim,
                    density: &density,
                    viscosity: &viscosity,
                    velocity: &velocity,
                    mesh_velocity: None,
                    load: &load,
                    convect,
                    newton: self.newton,
                    stab,
                    discretization: disc,
                    want_mass: transient,
                    np: None,
                    omega: angular_velocity(bf)?,
                    bubbles,
                    porous: porous.as_deref(),
                    thermal,
                    power_law: power_law.then_some(PowerLaw {
                        exponent: &exponent,
                        critical: critical.as_deref(),
                        nominal: nominal.as_deref(),
                    }),
                };
                match &pe {
                    None if bubbles => {
                        flow_bulk_local(basis, quad, &inp, &mut self.local)?;
                        // FlowSolve.F90:1159, 1222: nb = n, NSCondensate.
                        ns_condensate(n, n, vdim, &mut self.local)?;
                    }
                    None => flow_bulk_local(basis, quad, &inp, &mut self.local)?,
                    Some(p) => {
                        flow_bulk_local(&*p.basis, &FixedRule(p.rule.clone()), &inp, &mut self.local)?;
                        widen(&mut self.local, p.nd * nsdofs);
                    }
                }
                if transient {
                    ctx.default_first_order_time(&mut self.local, e, dt)?;
                }
                ctx.default_update(&self.local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                if dim - element_dim(el.etype.family()).min(dim) > 1 {
                    continue;
                }
                let bc = ctx.bc_params(e).expect("boundary_elements filters on BC");
                check_bc(bc, &var_name)?;
                if bc.logical("flow force bc") == Some(false) {
                    continue;
                }
                let n = el.nodes.len();
                let bc_nt = bc;
                let bc = Some(bc);
                let mut ext = nodal(ctx, bc, "external pressure", e, n)?;
                if bc.is_some_and(|b| b.get("external pressure").is_none()) {
                    ext = nodal(ctx, bc, "normal pressure", e, n)?;
                }
                let mut traction = vec![[0.0; 3]; n];
                let mut slip = vec![[0.0; 3]; n];
                for i in 0..3 {
                    let t = nodal(ctx, bc, &format!("pressure {}", i + 1), e, n)?;
                    let s = nodal(ctx, bc, &format!("slip coefficient {}", i + 1), e, n)?;
                    for k in 0..n {
                        traction[k][i] = t[k];
                        slip[k][i] = s[k];
                    }
                }
                let mass_flux = nodal(ctx, bc, "mass flux", e, n)?;
                let (alpha, beta) = surface_tension(ctx, bc_nt, e, n)?;
                let any = |v: &[f64]| v.iter().any(|&x| x != 0.0);
                if !any(&ext)
                    && !any(&mass_flux)
                    && !any(&alpha)
                    && !any(&beta)
                    && !traction.iter().any(|t| any(t))
                    && !slip.iter().any(|s| any(s))
                {
                    continue;
                }
                if cyl {
                    return Err(unsupported(
                        "boundary forces or slip in Cylindric Symmetric coordinates (NavierStokesCylindricalBoundary)",
                    ));
                }
                let parent_center = normal_parent(ctx, e).map(|pe| element_center(&node_xyz(ctx.mesh, ctx.mesh.element(pe).nodes)));
                let xyz = node_xyz(ctx.mesh, el.nodes);
                let inp = FlowBoundaryInput {
                    etype: el.etype,
                    xyz: &xyz,
                    dim,
                    traction: &traction,
                    ext_pressure: &ext,
                    mass_flux: &mass_flux,
                    slip: &slip,
                    parent_center,
                    np: None,
                    normal_tangential: normal_tangential(bc_nt, &var_name),
                    surface_tension: (any(&alpha) || any(&beta)).then_some((&alpha[..], &beta[..])),
                };
                match ctx.solver_dofs() {
                    None => flow_boundary_local(basis, quad, &inp, &mut self.local)?,
                    Some(d) => {
                        // NavierStokesBoundary: GaussPoints(Element) and the n
                        // vertex functions of the p boundary element.
                        let p = p_parts(d, e, el.etype, el.nodes, &stabilize_flag)?;
                        flow_boundary_local(&*p.basis, &FixedRule(p.rule.clone()), &inp, &mut self.local)?;
                        widen(&mut self.local, p.nd * nsdofs);
                    }
                }
                ctx.default_update(&self.local, e)?;
            }

            ctx.default_finish_boundary_assembly()?;
            // Implicit friction boundaries, on the matrix level.
            set_implicit_friction(ctx, dim)?;

            noslip_walls(ctx, &comps, vdim);
            // DefaultDirichletBCs: element values (in BC order with `Set
            // Dirichlet BCs by BC Numbering`), zero on p edge dofs, then the
            // point values of `Target Nodes` / `Target Coordinates`.
            let (fixed, rows) = ctx.dirichlet_values_with_rows()?;
            ctx.apply_dirichlet_values_with_rows(&fixed, &rows)?;
            let norm = ctx.default_solve()?;

            let change = nonlin_change(ctx);
            ctx.info("FlowSolve", &format!(" Result Norm     : {norm:>25.17E}"), 4);
            ctx.info("FlowSolve", &format!(" Relative Change : {change:>25.17E}"), 4);

            if change < newton_tol || iter > newton_iter {
                self.newton = true;
            }
            if self.recheck_newton && self.newton {
                if change > newton_ubound {
                    self.newton = false;
                    ctx.info("FlowSolve", "Newton tolerance exceeded, switching back to picard", 6);
                }
                if iter >= newton_max_iter {
                    self.newton = false;
                    ctx.info("FlowSolve", "Newton iteration limit exceeded, switching back to picard", 6);
                }
            }
            if change < nonlinear_tol && iter < nonlinear_iter {
                break;
            }
            if norm > max_norm {
                ctx.warn("FlowSolve", "Exiting as nonlinear norm is above allowed maximum!");
                break;
            }
            // FlowSolve.F90:1517-1529: move the free surfaces.
            if free_surface {
                if change < free_tol || iter > free_iter {
                    self.compute_free = true;
                }
                if self.compute_free && p.logical("internal move boundary").unwrap_or(true) {
                    if cyl {
                        return Err(unsupported("free surface in Cylindric Symmetric coordinates"));
                    }
                    let relax = ctx.creal(&p, "free surface relaxation factor")?.unwrap_or(1.0);
                    crate::free_surface::move_boundary(ctx, relax, basis, quad)?;
                }
            }
        }
        Ok(())
    }
}

/// `SetImplicitFriction(Model, Solver, 'Implicit Friction Coefficient')`
/// (SolverBasics.F90:10074): at each node of a boundary with an `Implicit
/// Friction Coefficient` c (once per node, in boundary-element order), the
/// row of the velocity component normal to the boundary (the component the
/// element normal points along most) is added, times c times the unit
/// direction of the current tangential velocity, to the tangential rows:
/// the tangential traction becomes c times the normal one. `Implicit
/// Friction Exclude Pressure` leaves the pressure columns out. Not ported:
/// `Implicit Friction Direction Vector` (refused) and normal-tangential
/// systems (FlowSolve refuses them).
fn set_implicit_friction(ctx: &mut SolverCtx<'_>, dim: usize) -> Result<(), SolverError> {
    const KEY: &str = "implicit friction coefficient";
    if !ctx.sim.model.boundary_conditions.iter().any(|bc| bc.get(KEY).is_some()) {
        return Ok(());
    }
    ctx.info("SetImplicitFriction", "Setting fluid friction for boundaries on matrix level!", 7);
    let exclude_pressure = flag(Some(ctx.params()), "implicit friction exclude pressure");
    let (dofs, perm, values) = {
        let f = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
        (f.dofs, f.perm.map(|p| p.to_vec()), f.values.to_vec())
    };
    let perm = perm.ok_or_else(|| SolverError::fatal("SetImplicitFriction", "flow variable has no permutation"))?;
    let mut done = vec![false; perm.len()];
    let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
    for e in bnd {
        let Some(bc) = ctx.bc_params(e) else { continue };
        if bc.get(KEY).is_none() {
            continue;
        }
        let bc = bc.clone();
        let el = ctx.mesh.element(e);
        let nodes = el.nodes.to_vec();
        let normal = element_normal(ctx.mesh, el.etype, &nodes)?;
        let coeff = ctx
            .real_at_nodes(&bc, KEY, e)?
            .ok_or_else(|| SolverError::fatal("SetImplicitFriction", "Implicit Friction Coefficient not found"))?;
        for (i, &j) in nodes.iter().enumerate() {
            let j = j as usize;
            if done[j] || perm[j] == NO_DOF {
                continue;
            }
            done[j] = true;
            let c = coeff[i];
            if c.abs() < 1.0e-10 {
                continue;
            }
            // The Cartesian component the normal is associated to.
            let mut dof_n = 0;
            for k in 1..dim {
                if normal[k].abs() > normal[dof_n].abs() {
                    dof_n = k;
                }
            }
            if normal[dof_n].abs() < 0.99 {
                ctx.warn("SetImplicitFriction", "No normal-tangential system for implicit friction!");
            }
            let dof_t1 = if dof_n == 0 { 1 } else { 0 };
            let dof_t2 = 3 - dof_t1 - dof_n;
            let slot = perm[j] as usize;
            let mut velo = [0.0f64; 3];
            velo[dof_t1] = values[dofs * slot + dof_t1];
            if dim == 3 {
                velo[dof_t2] = values[dofs * slot + dof_t2];
            }
            let abs = velo.iter().map(|v| v * v).sum::<f64>().sqrt();
            if abs <= f64::MIN_POSITIVE {
                continue;
            }
            for v in velo.iter_mut() {
                *v /= abs;
            }
            for v in velo.iter_mut() {
                *v *= c;
            }
            let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal("SetImplicitFriction", "no linear system"))?;
            let a = &mut sys.matrix;
            let k = dofs * slot + dof_n;
            let k2 = dofs * slot + dof_t1;
            let k3 = dofs * slot + dof_t2;
            sys.rhs[k2] += velo[dof_t1] * sys.rhs[k];
            if dim == 3 {
                sys.rhs[k3] += velo[dof_t2] * sys.rhs[k];
            }
            for l in a.row_ptr[k]..a.row_ptr[k + 1] {
                let col = a.col_idx[l] as usize;
                if dofs > dim && exclude_pressure && col % dofs == dofs - 1 {
                    continue;
                }
                let find = |row: usize| (a.row_ptr[row]..a.row_ptr[row + 1]).find(|&m| a.col_idx[m] as usize == col);
                let v = a.values[l];
                let l2 = find(k2).ok_or_else(|| SolverError::fatal("SetImplicitFriction", "matrix pattern"))?;
                a.values[l2] += velo[dof_t1] * v;
                if dim == 3 {
                    let l3 = find(k3).ok_or_else(|| SolverError::fatal("SetImplicitFriction", "matrix pattern"))?;
                    a.values[l3] += velo[dof_t2] * v;
                }
            }
        }
    }
    let n = done.iter().filter(|&&d| d).count();
    ctx.info("SetImplicitFriction", &format!("Number of friction nodes set: {n}"), 10);
    Ok(())
}

/// `NormalVector(Element, Nodes)` at the element centre without the
/// orientation check, for the direction it is associated to only: a line's
/// `(dy, -dx)`, a face's `dx/du x dx/dv` of its first three nodes.
fn element_normal(mesh: &dyn MeshView, etype: elmer_core::ElementType, nodes: &[u32]) -> Result<[f64; 3], SolverError> {
    let x: Vec<[f64; 3]> = nodes.iter().map(|&n| mesh.coord(n)).collect();
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let n = match etype.family() {
        2 => {
            let t = sub(x[1], x[0]);
            [t[1], -t[0], 0.0]
        }
        3 | 4 => {
            let (u, v) = (sub(x[1], x[0]), sub(x[2], x[0]));
            [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        }
        f => return Err(unsupported(format!("Implicit Friction on boundary elements of family {f}"))),
    };
    let s = n.iter().map(|v| v * v).sum::<f64>().sqrt();
    Ok(n.map(|v| v / s))
}

/// `Angular Velocity(3)` of the element's body force (FlowSolve.F90:745):
/// a rotating frame, Cartesian coordinates only (FlowSolve's own Fatal
/// otherwise; FlowSolve here refuses other coordinate systems anyway).
fn angular_velocity(bf: Option<&ParamMap>) -> Result<Option<[f64; 3]>, SolverError> {
    let Some(v) = bf.and_then(|b| b.get("angular velocity")) else { return Ok(None) };
    match v {
        Value::Real { data, .. } if data.len() >= 3 => Ok(Some([data[0], data[1], data[2]])),
        _ => Err(unsupported("Angular Velocity other than three constant reals")),
    }
}

/// FlowSolve.F90:1322-1328: `Normal-Tangential Velocity`, else
/// `Normal-Tangential <variable>`.
fn normal_tangential(bc: &ParamMap, var_name: &str) -> bool {
    match bc.logical("normal-tangential velocity") {
        Some(v) => v,
        None => flag(Some(bc), &format!("normal-tangential {var_name}")),
    }
}

/// A p element's basis, rule and dof count.
struct PParts {
    basis: Box<dyn ElementBasis + Send + Sync>,
    rule: Vec<elmer_core::QuadPoint>,
    nd: usize,
}

/// The p dofs of element `e` as FlowSolve's `Stabilization Method =
/// Stabilized` uses them: only the quadrilateral and brick families (and
/// their boundary lines and quadrilaterals), whose p vertex functions are
/// the nodal ones on the same reference element, with the nodes in mesh
/// order.
fn p_parts(
    d: &dyn elmer_solver::SolverDofs,
    e: ElemId,
    etype: elmer_core::ElementType,
    nodes: &[u32],
    method: &str,
) -> Result<PParts, SolverError> {
    if method != "stabilized" {
        return Err(unsupported(format!("p elements with Stabilization Method = {method}")));
    }
    if ![2, 4, 8].contains(&etype.family()) {
        return Err(unsupported(format!("p elements of type {} (only line, quadrilateral, brick)", etype.0)));
    }
    if d.element_nodes(e) != nodes {
        return Err(unsupported("p elements whose dofs reorder the nodes"));
    }
    let mut slots = Vec::new();
    d.element_slots(e, &mut slots);
    Ok(PParts {
        basis: d.element_basis(e).map_err(|m| SolverError::fatal("ElementInfo", m))?,
        rule: d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?,
        nd: slots.len(),
    })
}

/// A fixed integration rule (`GaussPoints(Element)` of a p element).
struct FixedRule(Vec<elmer_core::QuadPoint>);

impl Quadrature for FixedRule {
    fn gauss_points(&self, _etype: elmer_core::ElementType, _np: Option<usize>) -> Vec<elmer_core::QuadPoint> {
        self.0.clone()
    }
}

/// The local system of the first `local.n` rows widened to `m` with zero
/// rows and columns (the element's other dofs).
fn widen(local: &mut LocalSystem, m: usize) {
    let n = local.n;
    if m == n {
        return;
    }
    fn grow(a: &[f64], n: usize, m: usize) -> Vec<f64> {
        let mut out = vec![0.0; m * m];
        for r in 0..n {
            out[r * m..r * m + n].copy_from_slice(&a[r * n..r * n + n]);
        }
        out
    }
    local.stiff = grow(&local.stiff, n, m);
    local.mass = local.mass.as_deref().map(|a| grow(a, n, m));
    local.damp = local.damp.as_deref().map(|a| grow(a, n, m));
    local.force.resize(m, 0.0);
    local.n = m;
}

/// `Solver % Variable % NonlinChange` after the last solve.
fn nonlin_change(ctx: &SolverCtx<'_>) -> f64 {
    let id = ctx.solver_id();
    ctx.sim
        .slots
        .iter()
        .find(|s| s.id == id)
        .and_then(|s| s.variable)
        .map_or(0.0, |v| ctx.sim.vars.by_id(v).conv.nonlin_change)
}
