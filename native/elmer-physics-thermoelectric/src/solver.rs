//! `ThermoElectricSolver` (ThermoElectricSolver.F90): the
//! `Procedure = "ThermoElectricSolver" "ThermoElectricSolver"` solver, driving
//! the element kernels through [`SolverCtx`] with the Fortran's nonlinear
//! loop (Picard, with Newton switched on for good once the change or the
//! iteration count says so).

use std::collections::BTreeMap;
use std::sync::OnceLock;

use elmer_core::{ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, NO_DOF};
use elmer_solver::sim::parse_variable_spec;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::boundary::{thermoelectric_boundary_local, ThermoElectricBoundaryInput};
use crate::bubble::{parse_element, BubbleConfig};
use crate::condense::{add_backward_euler, condense_bubbles, BubbleRecovery};
use crate::bulk::{thermoelectric_bulk_local, Stabilization, Tensor, ThermoElectricBulkInput};
use crate::error::ThermoElectricError;

const CALLER: &str = "ThermoElectricSolver";

/// Builds the element basis and quadrature the solver integrates with, in
/// place of the `elmer-fe` defaults.
pub type FeProvider = fn() -> (Box<dyn ElementBasis + Send>, Box<dyn Quadrature + Send>);

static FE: OnceLock<FeProvider> = OnceLock::new();

/// Installs a basis and quadrature provider for [`factory`]-built solvers.
/// Returns `false` if one is already installed.
pub fn install_fe(provider: FeProvider) -> bool {
    FE.set(provider).is_ok()
}

/// Registers `("ThermoElectricSolver", "ThermoElectricSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    for (lib, entry) in crate::REGISTRY_ENTRIES {
        registry.register(lib, entry, factory);
    }
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`ThermoElectricSolver`].
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(match FE.get() {
        Some(provider) => {
            let (basis, quad) = provider();
            ThermoElectricSolver::new(basis, quad)
        }
        None => ThermoElectricSolver::default(),
    }))
}

impl From<ThermoElectricError> for SolverError {
    fn from(e: ThermoElectricError) -> Self {
        match e {
            ThermoElectricError::Unsupported(_) => SolverError::Unsupported(format!("{CALLER}: {e}")),
            e => SolverError::fatal(CALLER, e.to_string()),
        }
    }
}

fn unsupported(what: impl std::fmt::Display) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

fn flag(p: &ParamMap, key: &str) -> bool {
    p.logical(key).unwrap_or(false)
}

/// `GetStabilizeFlag`: `Stabilize`, overridden by `Stabilization Method`.
fn stabilize_flag(p: &ParamMap) -> bool {
    let mut s = flag(p, "stabilize");
    if let Some(m) = p.string("stabilization method") {
        match elmer_core::normalize_key(m).as_str() {
            "stabilized" => s = true,
            "bubbles" | "none" => s = false,
            _ => {}
        }
    }
    s
}

/// The thermoelectric solver.
pub struct ThermoElectricSolver {
    basis: Option<Box<dyn ElementBasis + Send>>,
    quad: Option<Box<dyn Quadrature + Send>>,
    local: LocalSystem,
    /// `NewtonLinearization`, a saved local of the Fortran: once on, on for
    /// every later call.
    newton: bool,
    /// Interior bubble values per element, `[T, P]` per bubble function: the
    /// latest recovery (`bx`) and the last timestep's (`bxprev`), the
    /// per-element history `CondensatePTransient` keeps.
    bubbles: BTreeMap<u32, Vec<f64>>,
    bubbles_prev: BTreeMap<u32, Vec<f64>>,
    /// `stimestep`: the timestep `bubbles_prev` belongs to.
    bubble_step: Option<i64>,
}

impl Default for ThermoElectricSolver {
    /// A solver on the `elmer-fe` basis and Gauss rules.
    fn default() -> Self {
        Self {
            basis: None,
            quad: None,
            local: LocalSystem::zeros(0),
            newton: false,
            bubbles: BTreeMap::new(),
            bubbles_prev: BTreeMap::new(),
            bubble_step: None,
        }
    }
}

impl ThermoElectricSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis: Some(basis),
            quad: Some(quad),
            ..Self::default()
        }
    }
}

/// Solver-section options this crate does not port. Returns the bubbles per
/// element of `Element = "p:1 b:N"` (zero without `b:`).
fn check_solver_params(p: &ParamMap) -> Result<usize, SolverError> {
    let mut nb = 0;
    if let Some(e) = p.string("element") {
        nb = parse_element(e)?;
    }
    if nb > 0 && p.logical("bubbles in global system") == Some(false) {
        return Err(unsupported("Bubbles in Global System = False (the bubbles are always condensed here)"));
    }
    Ok(nb)
}

/// Bubble-element setup of one solve: the `Element` definition, plus the
/// time scheme check (the condensation here is backward Euler).
fn bubble_config(ctx: &SolverCtx<'_>, nb: usize, transient: bool) -> Result<Option<BubbleConfig>, SolverError> {
    if nb == 0 {
        return Ok(None);
    }
    if transient {
        let sim = ctx.simulation_params();
        let method = ctx
            .params()
            .string("timestepping method")
            .or_else(|| sim.string("timestepping method"))
            .map(elmer_core::normalize_key)
            .unwrap_or_default();
        let order = ctx
            .params()
            .integer("bdf order")
            .or_else(|| sim.integer("bdf order"))
            .unwrap_or(2);
        let euler = method == "implicit euler" || (method == "bdf" && order == 1);
        if !euler {
            return Err(unsupported(format!(
                "transient bubble condensation with Timestepping Method = {method} (BDF order {order}); only backward Euler"
            )));
        }
    }
    let serendipity = ctx.simulation_params().logical("serendipity p elements").unwrap_or(true);
    Ok(Some(BubbleConfig { nb, serendipity }))
}

fn component_names(ctx: &SolverCtx<'_>) -> Result<Vec<String>, SolverError> {
    let spec = ctx
        .params()
        .string("variable")
        .ok_or_else(|| SolverError::fatal(CALLER, "solver has no 'Variable'"))?;
    Ok(parse_variable_spec(spec, ctx.params().integer("variable dofs"))?.components)
}

/// The Dirichlet rows of the two components (`DefaultDirichletBCs`).
fn dirichlet_rows(ctx: &SolverCtx<'_>, comps: &[String]) -> Result<Vec<(u32, f64)>, SolverError> {
    let dofs = comps.len();
    let mut fixed = std::collections::BTreeMap::new();
    let mut extra = Vec::new();
    let perm = ctx.primary().and_then(|f| f.perm);
    let slot = |n: u32| perm.map_or(n, |p| p[n as usize]);
    let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
    for e in bnd {
        let bc = ctx.bc_params(e).expect("boundary_elements filters on BC");
        let nodes = ctx.mesh.element(e).nodes;
        if let Some(d) = ctx.solver_dofs() {
            d.boundary_extra_slots(e, &mut extra);
        }
        for (k, key) in comps.iter().enumerate() {
            if bc.get(&format!("{key} condition")).is_some() {
                return Err(unsupported(format!("Dirichlet mask > {key} Condition <")));
            }
            let Some(vals) = ctx.real_at_nodes(bc, key, e)? else { continue };
            if !extra.is_empty() {
                return Err(unsupported(format!("Dirichlet > {key} < on the edge or face dofs of p elements")));
            }
            for (&n, &x) in nodes.iter().zip(&vals) {
                let s = slot(n);
                if s != NO_DOF {
                    fixed.insert(s * dofs as u32 + k as u32, x);
                }
            }
        }
    }
    Ok(fixed.into_iter().collect())
}

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

/// The basis, rule and node order of element `e` from the solver's p dofs
/// (`ElementInfo`, `GaussPoints`), with the bubbles in the global system.
struct PElement {
    basis: Box<dyn ElementBasis + Send + Sync>,
    rule: FixedRule,
    nodes: Vec<u32>,
}

fn p_element(d: &dyn elmer_solver::SolverDofs, e: ElemId) -> Result<PElement, SolverError> {
    if d.local_bubbles(e) > 0 {
        return Err(unsupported("Bubbles in Global System = False (the bubbles are always condensed here)"));
    }
    Ok(PElement {
        basis: d.element_basis(e).map_err(|m| SolverError::fatal("ElementInfo", m))?,
        rule: FixedRule(d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?),
        nodes: d.element_nodes(e),
    })
}

/// One element's rule as a [`Quadrature`].
struct FixedRule(Vec<elmer_core::QuadPoint>);

impl Quadrature for FixedRule {
    fn gauss_points(&self, _etype: elmer_core::ElementType, _np: Option<usize>) -> Vec<elmer_core::QuadPoint> {
        self.0.clone()
    }
}

/// `values` (at `mesh_nodes`) at `nodes` instead: a p element lists its
/// nodes in the order of its dofs.
fn reorder<T: Copy>(values: Vec<T>, mesh_nodes: &[u32], nodes: &[u32]) -> Vec<T> {
    if mesh_nodes == nodes {
        return values;
    }
    nodes.iter().map(|n| values[mesh_nodes.iter().position(|m| m == n).expect("same node set")]).collect()
}

fn node_xyz(mesh: &dyn elmer_core::MeshView, nodes: &[u32]) -> Vec<[f64; 3]> {
    nodes.iter().map(|&n| mesh.coord(n)).collect()
}

/// `GetReal(section, key, ..)`: the keyword at the nodes, zero where the
/// section lacks it.
fn nodal(ctx: &SolverCtx<'_>, p: Option<&ParamMap>, key: &str, e: ElemId, n: usize) -> Result<Vec<f64>, SolverError> {
    Ok(match p {
        Some(p) => ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0; n]),
        None => vec![0.0; n],
    })
}

/// `InputTensor`: `Name` as a tensor at the nodes, all zero when the section
/// lacks it. A scalar becomes the diagonal, a row or column vector the
/// diagonal of its first three entries, a matrix its leading 3x3 block.
fn input_tensor(p: Option<&ParamMap>, key: &str, ctx: &SolverCtx<'_>, e: ElemId, n: usize) -> Result<Vec<Tensor>, SolverError> {
    use elmer_core::Value;
    let zero = vec![[[0.0; 3]; 3]; n];
    let Some(v) = p.and_then(|p| p.get(key)) else { return Ok(zero) };
    match v {
        Value::Real { rows, cols, data } => {
            let mut t = [[0.0; 3]; 3];
            if *rows == 1 && *cols == 1 {
                for (i, row) in t.iter_mut().enumerate() {
                    row[i] = data[0];
                }
            } else if *rows == 1 || *cols == 1 {
                for (i, row) in t.iter_mut().enumerate().take(data.len().min(3)) {
                    row[i] = data[i];
                }
            } else {
                for i in 0..(*rows).min(3) {
                    for j in 0..(*cols).min(3) {
                        t[i][j] = data[i * cols + j];
                    }
                }
            }
            Ok(vec![t; n])
        }
        Value::Dependent(d) if d.rows * d.cols <= 1 => {
            let vals = ctx.real_at_nodes(p.expect("checked above"), key, e)?.unwrap_or_else(|| vec![0.0; n]);
            Ok(vals
                .iter()
                .map(|&x| std::array::from_fn(|i| std::array::from_fn(|j| if i == j { x } else { 0.0 })))
                .collect())
        }
        _ => Err(unsupported(format!("{key} as a tensor that varies with a field"))),
    }
}

/// `GetRealVector(BC, x, 'Electric Flux Vector', Found)`: the components
/// `<name> 1..3` at the nodes, `None` when none is given.
fn real_vector(ctx: &SolverCtx<'_>, p: &ParamMap, key: &str, e: ElemId, n: usize) -> Result<Option<Vec<[f64; 3]>>, SolverError> {
    let mut out = vec![[0.0; 3]; n];
    let mut found = false;
    for c in 0..3 {
        if let Some(v) = ctx.real_at_nodes(p, &format!("{key} {}", c + 1), e)? {
            found = true;
            for (o, x) in out.iter_mut().zip(v) {
                o[c] = x;
            }
        }
    }
    Ok(found.then_some(out))
}

impl Solver for ThermoElectricSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if let Some(s) = ctx.simulation_params().string("coordinate system") {
            let s = elmer_core::normalize_key(s);
            if !s.starts_with("cartesian") {
                return Err(unsupported(format!("Coordinate System = {s} (only Cartesian)")));
            }
        }
        check_solver_params(ctx.params())?;
        let dofs = ctx.primary().map(|f| f.dofs);
        if dofs != Some(2) {
            return Err(unsupported(format!(
                "a variable with {dofs:?} dofs (expected temperature and potential)"
            )));
        }
        let dim = ctx.mesh.dim();
        if self.basis.is_none() {
            self.basis = Some(Box::new(elmer_fe::LagrangeBasis::new(dim)));
        }
        if self.quad.is_none() {
            self.quad = Some(Box::new(elmer_fe::ElmerQuadrature));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        let nb = check_solver_params(ctx.params())?;
        // With p dofs (`SolverDofs`) the bubbles are in the global system, as
        // in Elmer; without, this crate condenses them.
        let pdofs = ctx.solver_dofs().is_some();
        let bubble = if pdofs { None } else { bubble_config(ctx, nb, transient)? };
        let p = ctx.params().clone();
        let dim = ctx.mesh.dim();
        let comps = component_names(ctx)?;
        let stabilize = stabilize_flag(&p);
        if pdofs && stabilize {
            return Err(unsupported("Stabilize on p elements"));
        }
        let long_edge = flag(ctx.simulation_params(), "stabilization use longest element edge");

        let newton_iter = p.integer("nonlinear system newton after iterations").unwrap_or(0);
        let newton_tol = p.real("nonlinear system newton after tolerance").unwrap_or(0.0);
        let max_iter = p.integer("nonlinear system max iterations").unwrap_or(1);
        let nonlin_tol = p.real("nonlinear system convergence tolerance").unwrap_or(0.0);

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut dofs = Vec::new();
        if bubble.is_some() && transient {
            let step = ctx.scalar("timestep").map_or(-1, |t| t.round() as i64);
            if self.bubble_step != Some(step) {
                self.bubbles_prev = self.bubbles.clone();
                self.bubble_step = Some(step);
            }
        }
        let var_name = ctx.variable_name().map(str::to_owned);
        let (Some(basis), Some(quad)) = (self.basis.as_deref(), self.quad.as_deref()) else {
            return Err(SolverError::fatal(CALLER, "solve before init"));
        };

        for iter in 1..=max_iter {
            ctx.default_initialize()?;
            let mut recoveries: Vec<(ElemId, Vec<usize>, BubbleRecovery)> = Vec::new();

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let pe = match ctx.solver_dofs() {
                    Some(d) => Some(p_element(d, e)?),
                    None => None,
                };
                let nodes = pe.as_ref().map_or(el.nodes, |p| &p.nodes[..]);
                let ord = |v| reorder(v, el.nodes, nodes);
                let ordt = |v| reorder(v, el.nodes, nodes);
                let xyz = node_xyz(ctx.mesh, nodes);
                let bp = ctx
                    .body_params(e)
                    .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
                let (mat, bf) = (bp.material, bp.body_force);

                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    return Err(SolverError::fatal(CALLER, format!("active element {} has nodes outside the system", e.0)));
                }
                let values = ctx.primary().expect("checked in init").values;
                let nf = dofs.len() / 2;
                let mut solution: Vec<[f64; 2]> =
                    (0..nf).map(|k| [values[dofs[2 * k] as usize], values[dofs[2 * k + 1] as usize]]).collect();
                let bubble_basis = match &bubble {
                    Some(cfg) => {
                        let b = self.bubbles.get(&e.0).cloned().unwrap_or_else(|| vec![0.0; 2 * nb]);
                        solution.extend(b.chunks(2).map(|c| [c[0], c[1]]));
                        Some(cfg.basis(el.etype)?)
                    }
                    None => None,
                };

                let hcond = ordt(input_tensor(mat, "heat conductivity", ctx, e, n)?);
                let alpha = ordt(input_tensor(mat, "seebeck coefficient", ctx, e, n)?);
                let sigma = ordt(input_tensor(mat, "electric conductivity", ctx, e, n)?);
                let epsil0 = if mat.is_some_and(|m| m.get("relative permittivity").is_some()) {
                    ctx.simulation_params().real("permittivity of vacuum").ok_or_else(|| {
                        SolverError::fatal(CALLER, "> Relative Permittivity < requires > Permittivity of Vacuum < !")
                    })?
                } else {
                    1.0
                };
                let epsil = ordt(input_tensor(mat, "electric permittivity", ctx, e, n)?);
                let rho = nodal(ctx, mat, "density", e, n)?;
                let cp = ord(nodal(ctx, mat, "heat capacity", e, n)?);

                let mut heat_source = vec![0.0; n];
                let mut current_source = vec![0.0; n];
                if let Some(b) = bf {
                    match ctx.real_at_nodes(b, "heat source", e)? {
                        Some(h) => {
                            for ((l, h), r) in heat_source.iter_mut().zip(h).zip(&rho) {
                                *l = r * h;
                            }
                        }
                        None => heat_source = nodal(ctx, bf, "volumetric heat source", e, n)?,
                    }
                    current_source = nodal(ctx, bf, "current source", e, n)?;
                }
                let (rho, heat_source, current_source) = (ord(rho), ord(heat_source), ord(current_source));

                let stab = if stabilize {
                    let (mk, hk) = elmer_fe::stab_param(el.etype, &xyz, dim, long_edge)
                        .map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
                    Some(Stabilization { mk, hk })
                } else {
                    None
                };
                let inp = ThermoElectricBulkInput {
                    etype: el.etype,
                    xyz: &xyz,
                    hcond: &hcond,
                    alpha: &alpha,
                    sigma: &sigma,
                    epsil: &epsil,
                    epsil0,
                    rho: &rho,
                    cp: &cp,
                    heat_source: &heat_source,
                    current_source: &current_source,
                    solution: &solution,
                    newton: self.newton,
                    stabilize: stab,
                    want_mass: transient,
                    np: bubble.as_ref().map(|c| c.gauss_points(el.etype.0)),
                };
                if let Some(p) = &pe {
                    thermoelectric_bulk_local(&*p.basis, &p.rule, &inp, &mut self.local)?;
                    if transient {
                        ctx.default_first_order_time(&mut self.local, e, dt)?;
                    }
                } else if let Some(bb) = &bubble_basis {
                    thermoelectric_bulk_local(bb, quad, &inp, &mut self.local)?;
                    if transient {
                        let name = var_name.as_deref().expect("a solver with a variable");
                        let prev = &ctx.sim.vars.get(name).expect("primary variable").prev_values;
                        let Some(prev) = prev.first() else {
                            return Err(SolverError::fatal(CALLER, "transient run without a previous timestep"));
                        };
                        let mut xprev: Vec<f64> = dofs.iter().map(|&d| prev[d as usize]).collect();
                        match self.bubbles_prev.get(&e.0) {
                            Some(b) => xprev.extend_from_slice(b),
                            None => xprev.resize(2 * (n + nb), 0.0),
                        }
                        add_backward_euler(&mut self.local, dt, &xprev);
                    }
                    let rec = condense_bubbles(&mut self.local, 2 * n)?;
                    recoveries.push((e, dofs.iter().map(|&d| d as usize).collect(), rec));
                } else {
                    thermoelectric_bulk_local(basis, quad, &inp, &mut self.local)?;
                    if transient {
                        ctx.default_first_order_time(&mut self.local, e, dt)?;
                    }
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
                let n = el.nodes.len();
                if dofs.len() != 2 * n {
                    return Err(unsupported("boundary elements with edge or face dofs (p elements)"));
                }
                let bc = ctx.bc_params(e).expect("boundary_elements filters on BC");
                // A p boundary element lists its nodes in its parent's edge or
                // face order; its values follow that order.
                let nodes = ctx.solver_dofs().map_or_else(|| el.nodes.to_vec(), |d| d.element_nodes(e));
                let ord = |v| reorder(v, el.nodes, &nodes);
                let xyz = node_xyz(ctx.mesh, &nodes);
                let heat_flux = ord(nodal(ctx, Some(bc), "heat flux", e, n)?);
                let ext_temperature = ord(nodal(ctx, Some(bc), "external temperature", e, n)?);
                let heat_transfer = ord(nodal(ctx, Some(bc), "heat transfer coefficient", e, n)?);
                let electric_flux = ord(nodal(ctx, Some(bc), "electric flux", e, n)?);
                let efv = real_vector(ctx, bc, "electric flux vector", e, n)?.map(|v| reorder(v, el.nodes, &nodes));
                let ext_potential = ord(nodal(ctx, Some(bc), "external potential", e, n)?);
                let etc = if bc.get("electric transfer coefficient").is_some() {
                    Some(reorder(input_tensor(Some(bc), "electric transfer coefficient", ctx, e, n)?, el.nodes, &nodes))
                } else {
                    None
                };
                let parent_center = [el.parents[0], el.parents[1]]
                    .into_iter()
                    .flatten()
                    .next()
                    .map(|pe| elmer_physics_flow::boundary::element_center(&node_xyz(ctx.mesh, ctx.mesh.element(pe).nodes)));
                let inp = ThermoElectricBoundaryInput {
                    etype: el.etype,
                    xyz: &xyz,
                    heat_flux: &heat_flux,
                    ext_temperature: &ext_temperature,
                    heat_transfer: &heat_transfer,
                    electric_flux: &electric_flux,
                    electric_flux_vector: efv.as_deref(),
                    ext_potential: &ext_potential,
                    transfer_tensor: etc.as_deref(),
                    parent_center,
                    np: None,
                };
                thermoelectric_boundary_local(basis, quad, &inp, &mut self.local)?;
                if transient {
                    // The boundary system has no mass matrix; Default1stOrderTime
                    // with none leaves it unchanged.
                }
                ctx.default_update(&self.local, e)?;
            }
            ctx.default_finish_boundary_assembly()?;

            let fixed = dirichlet_rows(ctx, &comps)?;
            apply_dirichlet(ctx, &fixed)?;
            let norm = ctx.default_solve()?;
            if !recoveries.is_empty() {
                let values = ctx.primary().expect("checked in init").values;
                for (e, idx, rec) in &recoveries {
                    let x: Vec<f64> = idx.iter().map(|&d| values[d]).collect();
                    self.bubbles.insert(e.0, rec.bubbles(&x));
                }
            }
            let change = ctx.convergence().map_or(0.0, |c| c.nonlin_change);
            ctx.info(CALLER, &format!(" Result Norm   : {norm}"), 4);
            ctx.info(CALLER, &format!(" Relative Change : {change}"), 4);

            if change < newton_tol || iter >= newton_iter {
                self.newton = true;
            }
            if change < nonlin_tol {
                break;
            }
        }
        Ok(())
    }
}
