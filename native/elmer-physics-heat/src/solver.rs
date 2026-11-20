//! `HeatSolver`: the `Procedure = "HeatSolve" "HeatSolver"` solver, driving the
//! element kernels through [`SolverCtx`] the way `HeatSolver` in
//! HeatSolve.F90:254 drives `LocalMatrixVec` and `LocalMatrixBC`.

use std::sync::OnceLock;

use std::collections::BTreeMap;

use elmer_core::{
    ElemId, ElementBasis, ElementType, LocalSystem, ParamMap, Params, QuadPoint, Quadrature, Value,
    NO_DOF,
};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::condense::{condensate_p, condensate_p_transient};
use crate::info_vec::ElementInfoVec;
use crate::material::{
    Conductivity, HeatBcAtNodes, HeatBoundaryConfig, HeatBulkConfig, Stabilization,
};
use crate::preflight::stabilized;
use crate::radiation::{is_diffuse_gray, temps4, DiffuseGray};
use crate::radiators::{is_radiator_bc, Radiators};
use crate::{
    check_supported, heat_boundary_local, heat_bulk_local_vec, CompressModel, Compressibility,
    HeatError, HeatMaterialAtNodes, HeatSections, JouleHeat, PROCEDURE,
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
    pub(crate) basis: Box<dyn ElementBasis + Send>,
    pub(crate) quad: Box<dyn Quadrature + Send>,
    pub(crate) local: LocalSystem,
    /// The diffuse gray radiation, found and factored on the first solve.
    diffuse: Option<Option<DiffuseGray>>,
    /// The point radiators, likewise.
    radiators: Option<Option<Radiators>>,
    /// The condensed bubble values per bulk element (`BubbleValues`), and
    /// those of the previous timestep (`BubblePrevValues`).
    bubbles: BTreeMap<u32, Vec<f64>>,
    bubbles_prev: BTreeMap<u32, Vec<f64>>,
    /// The timestep `bubbles_prev` was shifted at (`BubbleTimestep`).
    bubble_step: Option<i64>,
    /// RefineMesh's SAVEd state (`Adaptive Mesh Refinement`).
    refine: elmer_adaptive::RefineState,
}

impl HeatSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis,
            quad,
            local: LocalSystem::zeros(0),
            diffuse: None,
            radiators: None,
            bubbles: BTreeMap::new(),
            bubbles_prev: BTreeMap::new(),
            bubble_step: None,
            refine: Default::default(),
        }
    }
}

/// `ListGetElementReal` of a scalar keyword at the nodes of element `e`:
/// constants, tables, MATC and procedures through `elmer-expr`, evaluated at
/// the current iterate of the fields they depend on.
pub(crate) fn real_at_nodes(
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
pub(crate) fn conductivity(
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
pub(crate) fn convection_velocity(
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

/// The kernel inputs of bulk element `e` with `n` nodes: material and body
/// force evaluated at the element nodes. `legacy` forms the heat source as
/// HeatSolveLegacy does, the nodal product `Density * Heat Source`, instead
/// of `HeatSolveVec`'s product of interpolants. `dt` is the step of a
/// transient run (`None` when steady), for the `dp/dt` pressure work.
pub(crate) fn bulk_inputs(
    ctx: &SolverCtx<'_>,
    e: ElemId,
    n: usize,
    legacy: bool,
    dt: Option<f64>,
) -> Result<HeatMaterialAtNodes, SolverError> {
    const CALLER: &str = "HeatSolver";
    let bp = ctx
        .body_params(e)
        .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
    let (m, bf) = (bp.material, bp.body_force);
    let convection = bp
        .equation
        .and_then(|q| q.string("Convection"))
        .map(elmer_core::normalize_key);
    let computed = convection.as_deref() == Some("computed");
    let convected = computed || convection.as_deref() == Some("constant");
    // HeatSolveLegacy reads the flow field and the compressibility through
    // its own paths, which are not ported. Elmer's `init` accepts these, so
    // the refusal waits for the solve (a solver with `Exec Solver = Never`
    // never gets here).
    if legacy && computed {
        return Err(SolverError::Unsupported(
            "Convection = Computed with Legacy Assembly: reading the flow solution as the \
             convection field is not implemented"
                .into(),
        ));
    }
    if legacy && m.is_some_and(|m| m.get("Compressibility Model").is_some()) {
        return Err(SolverError::Unsupported(
            "material keyword \"Compressibility Model\" with Legacy Assembly: compressible \
             heat is not implemented"
                .into(),
        ));
    }
    let flow = if computed {
        Some(flow_at_nodes(ctx, e, dt.is_some())?)
    } else {
        None
    };
    let density = real_at_nodes(ctx, m, "Density", e)?;
    let mut heat_source = real_at_nodes(ctx, bf, "Heat Source", e)?;
    let mut volumetric_heat_source = real_at_nodes(ctx, bf, "Volumetric Heat Source", e)?;
    if legacy && volumetric_heat_source.is_none() {
        if let Some(src) = heat_source.take() {
            let rho = density.as_deref();
            volumetric_heat_source = Some(
                src.iter()
                    .enumerate()
                    .map(|(i, s)| s * rho.map_or(0.0, |r| r[i]))
                    .collect(),
            );
        }
    }
    Ok(HeatMaterialAtNodes {
        conductivity: match m {
            Some(m) => conductivity(ctx, m, e, n)?,
            None => None,
        },
        density,
        heat_capacity: real_at_nodes(ctx, m, "Heat Capacity", e)?,
        heat_transfer_multiplier: real_at_nodes(ctx, m, "Heat Transfer Multiplier", e)?,
        heat_source,
        volumetric_heat_source,
        convection_velocity: match (&flow, convected) {
            (Some(f), _) => Some(f.velocity.clone()),
            (None, true) => Some(convection_velocity(ctx, m, e, n)?),
            (None, false) => None,
        },
        convection_computed: flow.is_some(),
        mesh_velocity: match (&flow, convected) {
            (None, false) => None,
            _ => mesh_velocity(ctx, e),
        },
        compressibility: match m {
            Some(m) => compressibility(ctx, m, e, flow, dt)?,
            None => None,
        },
        joule_heat: match bf.is_some_and(|b| b.logical("Joule Heat") == Some(true)) {
            true => joule_heat(ctx, m, e, n)?,
            false => None,
        },
    })
}

/// `JouleHeat` (Differentials.F90:284): the field Joule heating takes, the
/// first of `Joule Heating e`, `Joule Field`, `Potential` and `Magnetic
/// Field 1` that exists. `None` (no heat) when none does, or when a node of
/// `e` is outside the field.
fn joule_heat(
    ctx: &SolverCtx<'_>,
    m: Option<&ParamMap>,
    e: ElemId,
    n: usize,
) -> Result<Option<JouleHeat>, SolverError> {
    let unsupported = |what: &str| {
        Err(SolverError::Unsupported(format!(
            "HeatSolver: Joule Heat from the {what} variable"
        )))
    };
    let vars = &ctx.sim.vars;
    let name = if vars.get("Joule Heating e").is_some() {
        return unsupported("Joule Heating e");
    } else if vars.get("Joule Field").is_some() {
        "Joule Field"
    } else if vars.get("Potential").is_some() {
        "Potential"
    } else if vars.get("Magnetic Field 1").is_some() {
        return unsupported("Magnetic Field 1");
    } else {
        ctx.warn(
            "JouleHeat",
            "Joule heating requested but no field to compute it!",
        );
        return Ok(None);
    };
    let var = vars.get(name).expect("checked above");
    let mut values = Vec::with_capacity(n);
    for &nn in ctx.mesh.element(e).nodes {
        match var.perm.as_deref().map_or(nn, |p| p[nn as usize]) {
            NO_DOF => return Ok(None),
            k => values.push(var.values[var.dofs * k as usize]),
        }
    }
    let electric_conductivity =
        real_at_nodes(ctx, m, "Electric Conductivity", e)?.unwrap_or_else(|| vec![0.0; n]);
    Ok(Some(if name == "Potential" {
        JouleHeat::Potential {
            electric_conductivity,
            potential: values,
        }
    } else {
        JouleHeat::Field {
            electric_conductivity,
            field: values,
        }
    }))
}

/// The convection field (`Temperature Convection Field`, else `Flow
/// Solution`) at the nodes of element `e`: its first `dim` components as the
/// velocity, its last as the pressure (`ConvField_h`, `PressureField_h` and,
/// transient, `PrevFlowField_h` of `LocalMatrixVec`).
struct FlowAtNodes {
    velocity: Vec<[f64; 3]>,
    pressure: Vec<f64>,
    prev_pressure: Option<Vec<f64>>,
}

fn flow_at_nodes(
    ctx: &SolverCtx<'_>,
    e: ElemId,
    transient: bool,
) -> Result<FlowAtNodes, SolverError> {
    const CALLER: &str = "HeatSolver";
    let name = ctx
        .params()
        .string("Temperature Convection Field")
        .unwrap_or("Flow Solution");
    let Some(var) = ctx.sim.vars.get(name) else {
        return Err(SolverError::fatal(
            CALLER,
            format!("Convection = Computed: no convection field \"{name}\""),
        ));
    };
    let dim = ctx.mesh.dim();
    let dofs = var.dofs;
    let nodes = ctx.mesh.element(e).nodes;
    let mut slots = Vec::with_capacity(nodes.len());
    for &nn in nodes {
        // `ListGetElementVectorSolutionVec` leaves an element with a node
        // outside the field's Perm without values (Fortran then reads a
        // null pointer).
        match var.perm.as_deref().map_or(nn, |p| p[nn as usize]) {
            NO_DOF => {
                return Err(SolverError::fatal(
                    CALLER,
                    format!(
                        "Convection = Computed: element {} has a node outside the \"{name}\" field",
                        e.0
                    ),
                ))
            }
            k => slots.push(k as usize),
        }
    }
    let velocity = slots
        .iter()
        .map(|&k| {
            let mut v = [0.0; 3];
            for (i, vi) in v.iter_mut().enumerate().take(dofs.min(dim)) {
                *vi = var.values[dofs * k + i];
            }
            v
        })
        .collect();
    // The pressure is component `dim + 1` (1-based) of the field.
    let at = |vals: &[f64]| -> Vec<f64> {
        slots
            .iter()
            .map(|&k| {
                if dim < dofs {
                    vals[dofs * k + dim]
                } else {
                    0.0
                }
            })
            .collect()
    };
    let pressure = at(&var.values);
    let prev_pressure = match transient {
        false => None,
        true => match var.prev_values.first() {
            Some(pv) => Some(at(pv)),
            None => {
                return Err(SolverError::fatal(
                    CALLER,
                    format!("Convection = Computed: \"{name}\" has no previous time step"),
                ))
            }
        },
    };
    Ok(FlowAtNodes {
        velocity,
        pressure,
        prev_pressure,
    })
}

/// The `Mesh Velocity` field at the nodes of `e`, which convection is
/// relative to (`MeshVeloField_h`); `None` unless it covers every node.
fn mesh_velocity(ctx: &SolverCtx<'_>, e: ElemId) -> Option<Vec<[f64; 3]>> {
    let mv = ctx.sim.vars.get("Mesh Velocity")?;
    let dim = ctx.mesh.dim();
    let nodes = ctx.mesh.element(e).nodes;
    nodes
        .iter()
        .map(
            |&nn| match mv.perm.as_deref().map_or(nn, |p| p[nn as usize]) {
                NO_DOF => None,
                k => {
                    let mut v = [0.0; 3];
                    for (i, vi) in v.iter_mut().enumerate().take(mv.dofs.min(dim)) {
                        *vi = mv.values[mv.dofs * k as usize + i];
                    }
                    Some(v)
                }
            },
        )
        .collect()
}

/// The material's `Compressibility Model` (`Thermal`, `Perfect Gas` or
/// `Perfect Gas Equation 1`) with the solver's temperature and, with
/// computed convection, the flow pressure. Any other value is incompressible,
/// as in `LocalMatrixVec`.
fn compressibility(
    ctx: &SolverCtx<'_>,
    m: &ParamMap,
    e: ElemId,
    flow: Option<FlowAtNodes>,
    dt: Option<f64>,
) -> Result<Option<Compressibility>, SolverError> {
    const CALLER: &str = "HeatSolver";
    let Some(name) = m
        .string("Compressibility Model")
        .map(elmer_core::normalize_key)
    else {
        return Ok(None);
    };
    let model = match name.as_str() {
        "thermal" => CompressModel::Thermal {
            heat_expansion: ctx.real_at_nodes(m, "Heat Expansion Coefficient", e)?,
            reference_temperature: ctx.real_at_nodes(m, "Reference Temperature", e)?,
        },
        "perfect gas" | "perfect gas equation 1" => {
            if flow.is_none() {
                return Err(SolverError::fatal(
                    CALLER,
                    "Compressibility Model \"Perfect Gas\" requires computed convection",
                ));
            }
            CompressModel::PerfectGas {
                specific_heat_ratio: ctx.real_at_nodes(m, "Specific Heat Ratio", e)?,
                reference_pressure: ctx.real_at_nodes(m, "Reference Pressure", e)?,
            }
        }
        "user defined" => {
            return Err(SolverError::Unsupported(
                "Compressibility Model = User Defined".into(),
            ))
        }
        _ => return Ok(None),
    };
    // `TempField_h`: the solver's own variable at the nodes.
    let Some(tv) = ctx.primary() else {
        return Err(SolverError::fatal(
            CALLER,
            "the heat solver has no Temperature variable",
        ));
    };
    let mut temperature = Vec::new();
    for &nn in ctx.mesh.element(e).nodes {
        match tv.perm.map_or(nn, |p| p[nn as usize]) {
            NO_DOF => {
                return Err(SolverError::fatal(
                    CALLER,
                    format!("element {} has a node without a temperature", e.0),
                ))
            }
            k => temperature.push(tv.values[k as usize * tv.dofs]),
        }
    }
    let (pressure, prev_pressure) = match flow {
        Some(f) => (Some(f.pressure), f.prev_pressure.zip(dt)),
        None => (None, None),
    };
    Ok(Some(Compressibility {
        model,
        temperature,
        pressure,
        prev_pressure,
    }))
}

pub(crate) fn axisymmetric(ctx: &SolverCtx<'_>) -> Result<bool, SolverError> {
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

/// The keywords `HeatSolver_Init0` and `HeatSolver_Init` (HeatSolve.F90)
/// add to the solver section (`ListAddNew*`: a keyword already there wins).
/// The `Element` must be in the section before the mesh and the solver
/// variable are created, so a caller building the model calls this before
/// [`Simulation::setup`](elmer_solver::Simulation::setup).
///
/// Without `Legacy Assembly`, convection (`Convection` or `Draw Velocity` in
/// any equation) or `Bubbles = True` makes the element a p element with
/// residual-free bubbles, condensed locally, unless the solver is
/// stabilized (`Stabilize`), when Elmer gives it `Element = n:1`, the nodal
/// element this crate uses anyway. A transient run with condensed bubbles
/// takes two nonlinear iterations.
pub fn solver_defaults(
    p: &mut ParamMap,
    equations: &[ParamMap],
    serendipity: bool,
    transient: bool,
) {
    if p.logical("Legacy Assembly") == Some(true) {
        return;
    }
    // HeatSolver_init: `ListAddLogical` (it overwrites).
    if p.logical("Discontinuous Galerkin") == Some(true) || p.logical("DG Reduced Basis") == Some(true) {
        p.insert("DG Indirect Connections", Value::Logical(true));
    }
    let convection = equations
        .iter()
        .any(|e| e.get("Convection").is_some() || e.get("Draw Velocity").is_some());
    if (convection || p.logical("Bubbles") == Some(true))
        && p.get("Element").is_none()
        && !stabilized(&*p)
    {
        p.insert("Element", Value::Str(bubble_element(serendipity).into()));
        add_new(p, "Bubbles in Global System", Value::Logical(false));
    }
    let bubbles = p.string("Element").is_some_and(|e| e.contains("b:"));
    if bubbles && transient && p.logical("Bubbles in Global System") != Some(true) {
        add_new(
            p,
            "Nonlinear System Min Iterations",
            Value::Integer(vec![2]),
        );
        add_new(
            p,
            "Nonlinear System Max Iterations",
            Value::Integer(vec![2]),
        );
    }
}

/// The `Element` `HeatSolver_Init0` adds for a residual-free bubble.
fn bubble_element(serendipity: bool) -> &'static str {
    if serendipity {
        "p:1 -tri b:1 -tetra b:1 -quad b:3 -brick b:4 -prism b:4 -pyramid b:4"
    } else {
        "p:1 -tri b:1 -tetra b:1 -quad b:4 -brick b:8 -prism b:4 -pyramid b:4"
    }
}

fn add_new(p: &mut ParamMap, key: &str, v: Value) {
    if p.get(key).is_none() {
        p.insert(key, v);
    }
}

/// A p element of the solver's dofs: `ElementInfoVec`'s basis,
/// `GaussPointsAdapt`'s rule, `GetElementNOFDOFs` (retained) and
/// `GetElementNOFBDOFs` (condensed bubbles).
struct PElement {
    basis: Box<dyn ElementBasis + Send + Sync>,
    rule: FixedRule,
    nd: usize,
    nb: usize,
}

fn p_element(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Option<PElement>, SolverError> {
    let Some(d) = ctx.solver_dofs() else {
        return Ok(None);
    };
    // A DG solver's elements are the nodal ones on DG slots.
    if d.dg_layout().is_some() {
        return Ok(None);
    }
    if d.element_nodes(e) != ctx.mesh.element(e).nodes {
        return Err(SolverError::Unsupported(
            "HeatSolver: a p element whose basis reorders its nodes".into(),
        ));
    }
    let mut slots = Vec::new();
    d.element_slots(e, &mut slots);
    // The edge DOFs of a 2D p element are fixed on Dirichlet boundaries by
    // the framework; 3D face DOFs are not.
    if slots.len() != ctx.mesh.element(e).nodes.len() && ctx.mesh.dim() == 3 {
        return Err(SolverError::Unsupported(
            "HeatSolver: 3D p elements beyond p:1 (face dofs)".into(),
        ));
    }
    Ok(Some(PElement {
        basis: d
            .element_basis(e)
            .map_err(|m| SolverError::fatal("ElementInfoVec", m))?,
        rule: FixedRule(
            d.element_rule(e)
                .map_err(|m| SolverError::fatal("GaussPoints", m))?,
        ),
        nd: slots.len(),
        nb: d.local_bubbles(e),
    }))
}

/// What a DG solver (`Discontinuous Galerkin` with `DG Reduced Basis`)
/// cannot do here: keywords that depend on a DG variable would be evaluated
/// at node numbers (`ListGetReal` reads a DG field through the element's
/// DG indexes), so a model whose material, body force or boundary keywords
/// depend on one is refused; so are `Heat Gap Parent` body jumps and `Jump
/// Integration Order`.
fn check_dg(ctx: &SolverCtx<'_>, var_name: &str) -> Result<(), SolverError> {
    let m = &ctx.sim.model;
    let dg_vars: Vec<String> = ctx
        .sim
        .vars
        .iter()
        .filter(|(_, v)| v.kind == elmer_core::FieldKind::Dg)
        .map(|(_, v)| v.name.clone())
        .chain(std::iter::once(elmer_core::normalize_key(var_name)))
        .collect();
    let sections = m.materials.iter().chain(&m.body_forces).chain(&m.boundary_conditions);
    for s in sections {
        for (key, v) in &s.0 {
            if let Value::Dependent(d) = v {
                if d.variables.iter().any(|x| dg_vars.contains(&elmer_core::normalize_key(x))) {
                    return Err(SolverError::Unsupported(format!(
                        "HeatSolver: > {key} < depends on a DG variable (evaluated through the DG indexes)"
                    )));
                }
            }
        }
    }
    if m.materials.iter().any(|mat| mat.get("Heat Gap Parent").is_some()) {
        return Err(SolverError::Unsupported("HeatSolver: Heat Gap Parent (body jumps of a DG solver)".into()));
    }
    if ctx.params().get("Jump Integration Order").is_some() {
        return Err(SolverError::Unsupported("HeatSolver: Jump Integration Order".into()));
    }
    Ok(())
}

/// The body id of bulk element `e`.
fn body_id(ctx: &SolverCtx<'_>, e: ElemId) -> usize {
    ctx.sim.model.body_of_tag(ctx.mesh.element(e).tag).unwrap_or(0)
}

/// The solver slot of the DG index `parent` gives `node`.
fn dg_slot(ctx: &SolverCtx<'_>, ix: &elmer_solver::DgIndexes, parent: ElemId, node: u32) -> Result<u32, SolverError> {
    let pos = ctx
        .mesh
        .element(parent)
        .nodes
        .iter()
        .position(|&m| m == node)
        .ok_or_else(|| SolverError::fatal("HeatSolver", "a DG boundary node is not in its parent"))?;
    let d = ix.index[parent.0 as usize][pos];
    let perm = ctx.primary().and_then(|f| f.perm).ok_or_else(|| SolverError::fatal("HeatSolver", "no permutation"))?;
    Ok(if d == NO_DOF { NO_DOF } else { perm[d as usize] })
}

/// `TempPerm(DgRadiationIndexes(Element, n, Indexes, .FALSE.))`
/// (ElementUtils.F90:1120): the DG indexes of the boundary element's nodes
/// in its parent, the one with the larger body id when there are two.
/// `TempPerm(DgRadiationIndexes(Element, n, Indexes, .TRUE.))`
/// (ElementUtils.F90, `DiffuseGray`): the DG indexes of the boundary
/// element's nodes in the first parent (left, then right) whose material
/// has an `Emissivity`.
pub(crate) fn dg_diffuse_gray_slots(ctx: &SolverCtx<'_>, ix: &elmer_solver::DgIndexes, e: ElemId) -> Result<Vec<u32>, SolverError> {
    const C: &str = "DGRadiationIndexes";
    let el = ctx.mesh.element(e);
    let parent = match el.parents {
        [Some(l), Some(r)] => {
            let mut found = None;
            for p in [l, r] {
                let b = ctx
                    .body_params(p)
                    .ok_or_else(|| SolverError::fatal(C, format!("Body {} has no Material associated!", body_id(ctx, p))))?;
                let m = b
                    .material
                    .ok_or_else(|| SolverError::fatal(C, format!("Body {} has no Material associated!", body_id(ctx, p))))?;
                if m.get("emissivity").is_some() {
                    found = Some(p);
                    break;
                }
            }
            found.ok_or_else(|| SolverError::fatal(C, "DG boundary parents should have emissivity!"))?
        }
        [Some(p), None] | [None, Some(p)] => p,
        [None, None] => return Err(SolverError::fatal(C, "DG boundary should have parents!")),
    };
    el.nodes.iter().map(|&n| dg_slot(ctx, ix, parent, n)).collect()
}

pub(crate) fn dg_radiation_slots(ctx: &SolverCtx<'_>, ix: &elmer_solver::DgIndexes, e: ElemId) -> Result<Vec<u32>, SolverError> {
    let el = ctx.mesh.element(e);
    let parent = match el.parents {
        [Some(l), Some(r)] => {
            if body_id(ctx, l) > body_id(ctx, r) {
                l
            } else {
                r
            }
        }
        [Some(p), None] | [None, Some(p)] => p,
        [None, None] => return Err(SolverError::fatal("DGRadiationIndexes", "DG boundary should have parents!")),
    };
    el.nodes.iter().map(|&n| dg_slot(ctx, ix, parent, n)).collect()
}

/// The DG jumps of `HeatSolver` (HeatSolve.F90:664-695): over the edges
/// (2D) or faces (3D) between two active elements, in mesh order, those a
/// boundary condition with `Heat Gap = True` covers get
/// `LocalJumpsDiscontBC` (HeatSolve.F90:3303): `alpha (T_l - T_r)` glued
/// with the standard Galerkin (or `Diagonal Jump Glue`) projector, `Heat
/// Gap Coefficient` alpha and half the `Heat Gap Flux` on each side, at the
/// left (first) and right parent's DG indexes. With `DG Reduced Basis` no
/// other jump terms are assembled.
fn dg_jumps(
    ctx: &mut SolverCtx<'_>,
    ix: &elmer_solver::DgIndexes,
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    axisym: bool,
    local: &mut LocalSystem,
) -> Result<(), SolverError> {
    let diag = ctx.params().logical("Diagonal Jump Glue") == Some(true);
    let alpha0 = ctx.params().real("Heat Gap Coefficient").unwrap_or(0.0);
    let beta0 = ctx.params().real("Heat Gap Flux").unwrap_or(0.0);
    // AssignConstraints: an edge (face) carries the BC of the boundary
    // element on the same nodes (the last one).
    let nb = ctx.mesh.n_bulk();
    let mut bc_of: std::collections::HashMap<Vec<u32>, ElemId> = std::collections::HashMap::new();
    for t in nb..nb + ctx.mesh.n_boundary() {
        let mut key = ctx.mesh.element(ElemId(t as u32)).nodes.to_vec();
        key.sort_unstable();
        bc_of.insert(key, ElemId(t as u32));
    }
    let active: std::collections::HashSet<ElemId> = ctx.active_elements().collect();
    for (l, r, nodes) in elmer_solver::dg_interfaces(ctx.mesh)? {
        if !active.contains(&l) || !active.contains(&r) {
            continue;
        }
        let mut key = nodes.clone();
        key.sort_unstable();
        let Some(&be) = bc_of.get(&key) else { continue };
        let Some(bc) = ctx.bc_params(be).cloned() else { continue };
        if bc.logical("Heat Gap") != Some(true) {
            continue;
        }
        let n = nodes.len();
        let etype = ElementType(match n {
            2 => 202,
            3 => 303,
            4 => 404,
            _ => return Err(SolverError::Unsupported(format!("DG jump over a {n}-node interface"))),
        });
        let mut slots = Vec::with_capacity(2 * n);
        for &nd in &nodes {
            slots.push(dg_slot(ctx, ix, l, nd)?);
        }
        for &nd in &nodes {
            slots.push(dg_slot(ctx, ix, r, nd)?);
        }
        if slots.contains(&NO_DOF) {
            return Err(SolverError::fatal("HeatSolver", "There should not be zero DG indexes!"));
        }
        // ListGetElementReal of the edge: a constant as it is, nodal values
        // (on the boundary element's nodes, in the edge's order) interpolated.
        let bel = ctx.mesh.element(be).nodes.to_vec();
        let nodal = |ctx: &SolverCtx<'_>, key: &str| -> Result<Option<Vec<f64>>, SolverError> {
            let Some(vals) = ctx.real_at_nodes(&bc, key, be)? else { return Ok(None) };
            Ok(Some(nodes.iter().map(|nd| vals[bel.iter().position(|m| m == nd).unwrap_or(0)]).collect()))
        };
        let alpha_n = match bc.get("Heat Gap Coefficient") {
            Some(Value::Dependent(_)) => nodal(ctx, "Heat Gap Coefficient")?,
            _ => None,
        };
        let beta_n = match bc.get("Heat Gap Flux") {
            Some(Value::Dependent(_)) => nodal(ctx, "Heat Gap Flux")?,
            _ => None,
        };
        let alpha_c = bc.real("Heat Gap Coefficient").unwrap_or(alpha0);
        let beta_c = bc.real("Heat Gap Flux").unwrap_or(beta0);
        let xyz: Vec<[f64; 3]> = nodes.iter().map(|&nd| ctx.mesh.coord(nd)).collect();
        let m = 2 * n;
        *local = LocalSystem::zeros(m);
        let mut be_eval = elmer_core::BasisEval::default();
        for qp in quad.gauss_points(etype, None) {
            basis
                .eval(etype, &xyz, qp.uvw, &mut be_eval)
                .map_err(|msg| SolverError::fatal("ElementInfo", msg))?;
            let b = &be_eval.basis[..n];
            let mut s = qp.weight * be_eval.det_j;
            if axisym {
                s *= xyz.iter().zip(b).fold(0.0, |acc, (x, &bi)| acc + x[0] * bi);
            }
            let interp = |v: &[f64]| v.iter().zip(b).fold(0.0, |acc, (&x, &bi)| acc + x * bi);
            let alpha = alpha_n.as_deref().map_or(alpha_c, interp);
            let beta = beta_n.as_deref().map_or(beta_c, interp);
            for p in 0..n {
                if diag {
                    local.stiff[p * m + p] += alpha * b[p] * s;
                    local.stiff[p * m + n + p] -= alpha * b[p] * s;
                    local.stiff[(p + n) * m + p] -= alpha * b[p] * s;
                    local.stiff[(p + n) * m + n + p] += alpha * b[p] * s;
                } else {
                    for q in 0..n {
                        local.stiff[p * m + q] += alpha * b[p] * b[q] * s;
                        local.stiff[p * m + n + q] -= alpha * b[p] * b[q] * s;
                        local.stiff[(p + n) * m + q] -= alpha * b[p] * b[q] * s;
                        local.stiff[(p + n) * m + n + q] += alpha * b[p] * b[q] * s;
                    }
                }
                local.force[p] += beta / 2.0 * b[p] * s;
                local.force[p + n] += beta / 2.0 * b[p] * s;
            }
        }
        ctx.update_global_equations(local, &slots)?;
    }
    Ok(())
}

/// One element's integration rule as a [`Quadrature`].
struct FixedRule(Vec<QuadPoint>);

impl Quadrature for FixedRule {
    fn gauss_points(&self, _: ElementType, _: Option<usize>) -> Vec<QuadPoint> {
        self.0.clone()
    }
}

/// `Default1stOrderTimeBubble`: a condensed bubble's history is BDF(1).
fn check_bubble_time_integration(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    let method = ctx
        .params()
        .string("Timestepping Method")
        .or_else(|| ctx.simulation_params().string("Timestepping Method"))
        .map(elmer_core::normalize_key);
    let order = ctx
        .params()
        .integer("BDF Order")
        .or_else(|| ctx.simulation_params().integer("BDF Order"))
        .unwrap_or(1);
    let method_ok = method
        .as_deref()
        .is_none_or(|m| m == "implicit euler" || m == "bdf");
    if order != 1 || !method_ok {
        return Err(SolverError::fatal(
            "Default1stOrderTime",
            format!(
                "a locally condensed bubble only supports Implicit Euler / BDF order 1, got \
                 Timestepping Method = {:?}, order {order}",
                method.unwrap_or_default()
            ),
        ));
    }
    Ok(())
}

impl Solver for HeatSolver {
    /// Refuses every keyword in the sections this solver reads that the
    /// kernels would otherwise ignore (see [`check_supported`]).
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        // HeatSolver_init (HeatSolve.F90:201): `Calculate Conductance
        // Matrix` is the library's constraint-mode lumping, with the
        // `Conductivity Body` of the BCs as the modes.
        elmer_solver::constraint_modes_lumping(ctx, &elmer_solver::CONDUCTANCE);
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
        }
        // The bubble element of `HeatSolver_Init0` (see [`solver_defaults`])
        // must have been added before the solver's dofs were numbered.
        let legacy = ctx.params().logical("Legacy Assembly") == Some(true);
        let convection = ctx
            .sim
            .model
            .equations
            .iter()
            .any(|e| e.get("Convection").is_some() || e.get("Draw Velocity").is_some())
            || ctx.params().logical("Bubbles") == Some(true);
        if legacy && ctx.solver_dofs().is_some() {
            return Err(SolverError::Unsupported(
                "HeatSolver: Legacy Assembly with a p element (Element)".into(),
            ));
        }
        if !legacy
            && convection
            && ctx.params().get("Element").is_none()
            && !stabilized(ctx.params())
        {
            return Err(SolverError::fatal(
                "HeatSolver_Init0",
                "convection without SUPG needs the residual-free bubble element; \
                 call elmer_physics_heat::solver_defaults before the setup",
            ));
        }
        // `AddEquationBasics` computes the Gebhart factors of a `Radiation
        // Solver` before the permutation (MainUtils.F90:1712-1721), so the
        // bandwidth optimisation sees their couplings.
        if ctx.params().logical("Legacy Assembly") != Some(true) && self.diffuse.is_none() {
            let axisym = axisymmetric(ctx)?;
            let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
            let dg = DiffuseGray::setup(ctx, &bnd, axisym)?;
            if let Some(dg) = &dg {
                let extra = dg.couplings(ctx.mesh);
                ctx.reorder_perm(&extra)?;
            }
            self.diffuse = Some(dg);
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
        if ctx.params().logical("Legacy Assembly") == Some(true) {
            self.solve_legacy(ctx, transient)?;
            // `DefaultFinish` (HeatSolve.F90): `Apply Explicit Control`.
            return ctx.default_finish();
        }
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
        // DefaultBubbleHistoryUpdate: the bubbles of the last solve become
        // the previous timestep's, once per timestep.
        if transient {
            let step = ctx.scalar("timestep").map_or(-1, |t| t.round() as i64);
            if self.bubble_step != Some(step) {
                self.bubbles_prev = self.bubbles.clone();
                self.bubble_step = Some(step);
            }
        }
        let var_name = ctx.variable_name().unwrap_or_default().to_string();
        let info_vec = ElementInfoVec { cdim: dim };
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

        if self.diffuse.is_none() {
            self.diffuse = Some(DiffuseGray::setup(ctx, &bnd, axisym)?);
        }
        if self.radiators.is_none() {
            self.radiators = Some(Radiators::setup(ctx, &bnd)?);
        }
        let stefan = ctx.constants().real("Stefan Boltzmann");
        // A DG solver (`Discontinuous Galerkin` with `DG Reduced Basis`).
        let dg_ix = ctx.solver_dofs().and_then(|d| d.dg_layout()).cloned();
        if dg_ix.is_some() {
            check_dg(ctx, &var_name)?;
            if self.radiators.as_ref().is_some_and(Option::is_some) {
                return Err(SolverError::Unsupported("HeatSolver: radiators of a DG solver".into()));
            }
        }
        if self.diffuse.as_ref().is_some_and(Option::is_some) && stefan.is_none() {
            return Err(SolverError::fatal(
                CALLER,
                "Stefan Boltzmann is not defined".to_string(),
            ));
        }

        for iter in 1..=maxiter {
            ctx.info(CALLER, &format!("Heat solver iteration: {iter}"), 4);
            ctx.default_initialize()?;
            if let Some(Some(dg)) = self.diffuse.as_mut() {
                dg.extend_topology(ctx)?;
            }
            let newton = {
                let id = ctx.solver_id();
                ctx.sim
                    .slots
                    .iter()
                    .find(|s| s.id == id)
                    .is_some_and(|s| s.state.newton_active)
            };
            let t4 = match self.diffuse {
                // The imported surfaces of a partitioned run too.
                Some(Some(_)) => {
                    let all: Vec<ElemId> = bnd.iter().copied().chain(ctx.foreign_boundary_elements()).collect();
                    temps4(ctx, &all)?
                }
                _ => Default::default(),
            };

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let n = xyz.len();
                let mat = bulk_inputs(ctx, e, n, false, transient.then_some(dt))?;
                bulk_cfg.stabilization = if supg && mat.convection_velocity.is_some() {
                    let (m_k, h_k) = elmer_fe::stab_param(el.etype, &xyz, dim, long_edge)
                        .map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
                    Stabilization::Supg { h_k, m_k }
                } else {
                    Stabilization::None
                };
                let pe = p_element(ctx, e)?;
                let (nd, nb) = pe.as_ref().map_or((n, 0), |p| (p.nd, p.nb));
                bulk_cfg.extra_functions = nd + nb - n;
                // LocalMatrixVec evaluates the basis with ElementInfoVec.
                let (basis, quad): (&dyn ElementBasis, &dyn Quadrature) = match &pe {
                    Some(p) => (&*p.basis, &p.rule),
                    None => (&info_vec, &*self.quad),
                };
                heat_bulk_local_vec(
                    basis,
                    quad,
                    el.etype,
                    &xyz,
                    &mat,
                    &bulk_cfg,
                    &mut self.local,
                )?;
                if transient && nb > 0 {
                    // Default1stOrderTime with Nb: implicit Euler on the
                    // whole element, the bubbles' previous values from the
                    // history, then CondensatePTransient.
                    check_bubble_time_integration(ctx)?;
                    ctx.element_dofs(e, &mut dofs);
                    let var = ctx.sim.vars.get(&var_name).ok_or_else(|| {
                        SolverError::fatal(CALLER, "the solver variable is missing")
                    })?;
                    let Some(prev) = var.prev_values.first() else {
                        return Err(SolverError::fatal(
                            CALLER,
                            "transient run without a previous timestep",
                        ));
                    };
                    let xprev: Vec<f64> = dofs.iter().map(|&d| prev[d as usize]).collect();
                    let x: Vec<f64> = dofs.iter().map(|&d| var.values[d as usize]).collect();
                    let zero = vec![0.0; nb];
                    let bprev = self.bubbles_prev.get(&e.0).unwrap_or(&zero);
                    let cur =
                        condensate_p_transient(&mut self.local, nd, nb, dt, &xprev, &x, bprev)?;
                    self.bubbles.insert(e.0, cur);
                } else {
                    if transient {
                        ctx.default_first_order_time(&mut self.local, e, dt)?;
                    }
                    condensate_p(&mut self.local, nd, nb)?;
                }
                ctx.default_update(&self.local, e)?;
            }

            ctx.default_finish_bulk_assembly()?;

            for &e in &bnd {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                let mut bc = {
                    let p = ctx.bc_params(e);
                    // A flux of an integration-point variable, at the points.
                    let flux_at_ips = match p {
                        Some(p) => {
                            let el = ctx.mesh.element(e);
                            let uvw: Vec<[f64; 3]> =
                                self.quad.gauss_points(el.etype, bnd_cfg.quad_points).iter().map(|q| q.uvw).collect();
                            ctx.boundary_real_at_ips(p, "heat flux", e, &uvw)?
                        }
                        None => None,
                    };
                    HeatBcAtNodes {
                        heat_flux: match flux_at_ips {
                            Some(_) => None,
                            None => real_at_nodes(ctx, p, "Heat Flux", e)?,
                        },
                        heat_flux_at_ips: flux_at_ips,
                        heat_transfer_coefficient: real_at_nodes(
                            ctx,
                            p,
                            "Heat Transfer Coefficient",
                            e,
                        )?,
                        external_temperature: real_at_nodes(ctx, p, "External Temperature", e)?,
                    }
                };
                // Radiator BC: the radiators' power adds to the flux, unless
                // the element is a diffuse gray one.
                if let Some(Some(r)) = &self.radiators {
                    if is_radiator_bc(ctx, e) && !is_diffuse_gray(ctx, e) {
                        if let Some(q) = r.flux(e) {
                            if bc.heat_flux_at_ips.is_some() {
                                return Err(SolverError::Unsupported(
                                    "HeatSolver: a radiator BC with a heat flux of integration-point variables".into(),
                                ));
                            }
                            let n = ctx.mesh.element(e).nodes.len();
                            let flux = bc.heat_flux.get_or_insert_with(|| vec![0.0; n]);
                            flux.iter_mut().for_each(|x| *x += q);
                        }
                    }
                }
                if bc.heat_flux.is_some() || bc.heat_flux_at_ips.is_some() || bc.heat_transfer_coefficient.is_some() {
                    let el = ctx.mesh.element(e);
                    // LocalMatrixBC spans the boundary's p DOFs, which this
                    // nodal boundary assembly leaves out.
                    if let (Some(d), Some(parent)) = (ctx.solver_dofs(), el.parents[0].or(el.parents[1])) {
                        let mut ps = Vec::new();
                        d.element_slots(parent, &mut ps);
                        if d.dg_layout().is_none() && ps.len() > ctx.mesh.element(parent).nodes.len() {
                            return Err(SolverError::Unsupported(
                                "HeatSolver: a heat flux or transfer BC on a p element beyond p:1".into(),
                            ));
                        }
                    }
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
                    match &dg_ix {
                        // LocalMatrixBC of a DG solver: UpdateGlobalEquations at
                        // TempPerm(DgRadiationIndexes).
                        Some(ix) => {
                            let slots = dg_radiation_slots(ctx, ix, e)?;
                            ctx.update_global_equations(&self.local, &slots)?;
                        }
                        None => ctx.default_update(&self.local, e)?,
                    }
                }
                if let Some(Some(dg)) = &self.diffuse {
                    if is_diffuse_gray(ctx, e) {
                        let sigma = stefan.expect("checked above");
                        dg.assemble(ctx, &*self.basis, &*self.quad, e, &t4, newton, sigma)?;
                    }
                }
            }

            if let Some(ix) = &dg_ix {
                dg_jumps(ctx, ix, &*self.basis, &*self.quad, axisym, &mut self.local)?;
            }

            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        // `DefaultFinish` (HeatSolve.F90): `Apply Explicit Control`.
        ctx.default_finish()?;
        // HeatSolve.F90:789: RefineMesh with HeatSolve's residuals.
        let p = ctx.params();
        if p.logical("adaptive mesh refinement") == Some(true) && p.logical("library adaptivity") != Some(true) {
            crate::adaptive::HeatEstimator::check(ctx)?;
            elmer_adaptive::refine_mesh(ctx, &mut self.refine, &mut crate::adaptive::HeatEstimator)?;
        }
        Ok(())
    }
}
