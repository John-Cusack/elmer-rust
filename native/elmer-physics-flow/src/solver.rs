//! `FlowSolver`: the `Procedure = "FlowSolve" "FlowSolver"` solver, driving
//! the element kernels through [`SolverCtx`] the way `FlowSolver`
//! (FlowSolve.F90:117) drives `NavierStokesCompose` and
//! `NavierStokesBoundary`, with its Picard/Newton switching.

use std::sync::OnceLock;

use elmer_core::{ElemId, ElementBasis, LocalSystem, MeshView, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::sim::parse_variable_spec;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::boundary::{element_center, flow_boundary_local, FlowBoundaryInput};
use crate::bulk::{flow_bulk_local, Discretization, FlowBulkInput};
use crate::error::FlowError;
use crate::bulk::StabParams;

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

/// Registers `("FlowSolve", "FlowSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    let (lib, entry) = crate::REGISTRY_ENTRY;
    registry.register(lib, entry, factory);
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

fn has_prefix(p: &ParamMap, prefix: &str) -> bool {
    p.0.keys().any(|k| k.starts_with(prefix))
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
fn check_stabilization(p: &ParamMap) -> Result<(), SolverError> {
    let method = match p.string("stabilization method") {
        Some(m) => elmer_core::normalize_key(m),
        None if flag(Some(p), "stabilize") => "stabilized".into(),
        None if flag(Some(p), "bubbles") || p.get("element").is_some() => "bubbles".into(),
        None => "stabilized".into(),
    };
    match method.as_str() {
        "stabilized" => {
            if flag(Some(p), "bubbles") {
                return Err(SolverError::fatal(CALLER, "You cant have stabilization and bubbles both!"));
            }
            Ok(())
        }
        "bubbles" | "pbubbles" | "p2/p1" | "p2p1" | "vms" => {
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
        ("calculate loads", "load computation"),
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
                "perfect gas" | "perfect gas equation 1" | "thermal" | "user defined" | "pressure dependent"
                    | "artificial compressible"
            ) {
                return Err(unsupported(format!("Compressibility Model = {c}")));
            }
        }
        if m.get("viscosity model").is_some() {
            return Err(unsupported("Viscosity Model (non-newtonian or anisotropic viscosity)"));
        }
        if flag(Some(m), "porous media") {
            return Err(unsupported("porous media"));
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
        if has_prefix(b, "angular velocity") {
            return Err(unsupported("rotating frame (Angular Velocity)"));
        }
    }
    let _ = eq;
    Ok(())
}

/// Boundary-condition options this crate does not port.
fn check_bc(bc: &ParamMap, var_name: &str) -> Result<(), SolverError> {
    for (key, what) in [
        ("free surface", "free surface"),
        ("normal-tangential velocity", "normal-tangential velocity"),
        ("wall law", "wall law"),
        ("vms wall", "VMS wall"),
    ] {
        if flag(Some(bc), key) {
            return Err(unsupported(what));
        }
    }
    if flag(Some(bc), &format!("normal-tangential {var_name}")) {
        return Err(unsupported("normal-tangential velocity"));
    }
    for (prefix, what) in [
        ("surface tension", "surface tension"),
        ("implicit friction", "implicit friction"),
        ("periodic bc", "periodic boundary conditions"),
        ("target coordinates", "point boundary conditions (Target Coordinates)"),
        ("target nodes", "point boundary conditions (Target Nodes)"),
    ] {
        if has_prefix(bc, prefix) {
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

fn cartesian(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    match ctx.simulation_params().string("coordinate system") {
        None => Ok(()),
        Some(s) if elmer_core::normalize_key(s).starts_with("cartesian") => Ok(()),
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

/// The Dirichlet rows of the solver variable from its component names: for
/// each boundary element in mesh order, `(row, value)` for every component
/// its BC sets, a later element overwriting an earlier one
/// (DefaultDirichletBCs). `noslip` also fixes the velocity components to
/// zero on `Noslip Wall BC` boundaries, as FlowSolve adds `Velocity i = 0`
/// to those BCs after its first boundary assembly.
fn dirichlet_rows(ctx: &SolverCtx<'_>, comps: &[String], dim: usize, noslip: bool) -> Result<Vec<(u32, f64)>, SolverError> {
    let dofs = comps.len();
    let mut fixed = std::collections::BTreeMap::new();
    let mut idx = Vec::new();
    let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
    for e in bnd {
        let bc = ctx.bc_params(e).expect("boundary_elements filters on BC");
        ctx.element_slots(e, &mut idx);
        let wall = noslip && flag(Some(bc), "noslip wall bc");
        for (k, key) in comps.iter().enumerate() {
            if bc.get(&format!("{key} condition")).is_some() {
                return Err(unsupported(format!("Dirichlet mask > {key} Condition <")));
            }
            let vals = if wall && k < dim {
                Some(vec![0.0; idx.len()])
            } else {
                ctx.real_at_nodes(bc, key, e)?
            };
            let Some(vals) = vals else { continue };
            for (&s, &x) in idx.iter().zip(&vals) {
                if s != NO_DOF {
                    fixed.insert(s * dofs as u32 + k as u32, x);
                }
            }
        }
    }
    Ok(fixed.into_iter().collect())
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
    let mut set = vec![false; n];
    let mut slots = Vec::new();
    let elems: Vec<ElemId> = ctx.active_elements().collect();
    for e in elems {
        let Some(ic) = ctx.body_params(e).and_then(|b| b.initial_condition) else { continue };
        ctx.element_slots(e, &mut slots);
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

/// `DefaultDirichletBCs` with the component names of the flow variable.
fn apply_dirichlet(ctx: &mut SolverCtx<'_>, fixed: &[(u32, f64)]) -> Result<(), SolverError> {
    let id = ctx.solver_id();
    let sim = &mut *ctx.sim;
    let backend = sim
        .backends
        .assembly
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
    let slot = sim
        .slots
        .iter_mut()
        .find(|s| s.id == id)
        .expect("the running solver has a slot");
    let sys = slot
        .system
        .as_mut()
        .ok_or_else(|| SolverError::fatal("DefaultDirichletBCs", "call default_initialize first"))?;
    backend.try_apply_dirichlet(&mut sys.matrix, &mut sys.rhs, fixed, &slot.params)
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
        cartesian(ctx)?;
        check_solver_params(ctx.params())?;
        let dim = ctx.mesh.dim();
        let dofs = ctx.primary().map(|f| f.dofs);
        if dofs != Some(dim + 1) {
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
        cartesian(ctx)?;
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
        if nsdofs != dim + 1 {
            return Err(unsupported(format!("{nsdofs} flow dofs in {dim}D")));
        }

        let disc = Discretization {
            div: flag(Some(&p), "div discretization"),
            gradp: flag(Some(&p), "gradp discretization"),
            laplace: flag(Some(&p), "laplace discretization"),
        };
        let h_scale = p.real("h scale").unwrap_or(1.0);
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
                let values = ctx.primary().expect("checked in init").values;
                let velocity: Vec<[f64; 3]> = (0..n)
                    .map(|k| std::array::from_fn(|i| if i < dim { values[dofs[k * nsdofs + i] as usize] } else { 0.0 }))
                    .collect();

                let density = nodal(ctx, mat, "density", e, n)?;
                let viscosity = nodal(ctx, mat, "viscosity", e, n)?;
                let mut load = vec![[0.0; 4]; n];
                if let Some(b) = bf {
                    let hydrostatic = b
                        .logical("hydrostatic pressure")
                        .or_else(|| eq.and_then(|q| q.logical("hydrostatic pressure")))
                        .unwrap_or(false);
                    if hydrostatic {
                        for l in load.iter_mut() {
                            l[..dim].copy_from_slice(&grav[..dim]);
                        }
                    }
                    for i in 0..dim {
                        if let Some(f) = ctx.real_at_nodes(b, &format!("flow bodyforce {}", i + 1), e)? {
                            for (l, v) in load.iter_mut().zip(f) {
                                l[i] += v;
                            }
                        }
                    }
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
                };
                flow_bulk_local(basis, quad, &inp, &mut self.local)?;
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
                let any = |v: &[f64]| v.iter().any(|&x| x != 0.0);
                if !any(&ext)
                    && !any(&mass_flux)
                    && !traction.iter().any(|t| any(t))
                    && !slip.iter().any(|s| any(s))
                {
                    continue;
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
                };
                flow_boundary_local(basis, quad, &inp, &mut self.local)?;
                ctx.default_update(&self.local, e)?;
            }

            let fixed = dirichlet_rows(ctx, &comps, dim, true)?;
            apply_dirichlet(ctx, &fixed)?;
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
        }
        Ok(())
    }
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
