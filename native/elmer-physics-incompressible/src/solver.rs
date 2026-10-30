//! `IncompressibleNSSolver` (IncompressibleNS.F90:2509): the
//! `Procedure = "IncompressibleNS" "IncompressibleNSSolver"` solver, driving
//! the element kernels through [`SolverCtx`] with the nonlinear loop of the
//! Fortran (Picard, switching to Newton as `ComputeChange` decides).

use std::sync::OnceLock;

use elmer_core::{ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::sim::parse_variable_spec;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::boundary::{incompressible_boundary_local, incompressible_boundary_local_with, IncompressibleBoundaryInput};
use crate::bulk::{
    incompressible_bulk_local_sums, lcondensate, lcondensate_bdf1, pin_pressure, HMode, ViscositySums,
    Compressibility, CompressModel, GlenFields, GlenParams, IncompressibleBulkInput, PressureStab, ViscosityModel,
};
use crate::error::IncompressibleError;

const CALLER: &str = "IncompressibleNSSolver";

/// Builds the element basis and quadrature the solver integrates with, in
/// place of the `elmer-fe` defaults.
pub type FeProvider = fn() -> (Box<dyn ElementBasis + Send>, Box<dyn Quadrature + Send>);

static FE: OnceLock<FeProvider> = OnceLock::new();

/// Installs a basis and quadrature provider for [`factory`]-built solvers.
/// Returns `false` if one is already installed.
pub fn install_fe(provider: FeProvider) -> bool {
    FE.set(provider).is_ok()
}

/// Registers `("IncompressibleNS", "IncompressibleNSSolver")` and its
/// `IncompressibleNSVec` front.
pub fn register(registry: &mut SolverRegistry) {
    for (lib, entry) in crate::REGISTRY_ENTRIES {
        registry.register(lib, entry, factory);
    }
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`IncompressibleNsSolver`].
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(match FE.get() {
        Some(provider) => {
            let (basis, quad) = provider();
            IncompressibleNsSolver::new(basis, quad)
        }
        None => IncompressibleNsSolver::default(),
    }))
}

impl From<IncompressibleError> for SolverError {
    fn from(e: IncompressibleError) -> Self {
        match e {
            IncompressibleError::Unsupported(_) => SolverError::Unsupported(format!("{CALLER}: {e}")),
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

fn has_prefix(p: &ParamMap, prefix: &str) -> bool {
    p.0.keys().any(|k| k.starts_with(prefix))
}

/// `ListAddString(NextFreeKeyword('Exported Variable'), spec)`, unless the
/// section exports that variable already.
fn add_exported(p: &mut ParamMap, spec: &str) {
    let name = spec.split_whitespace().filter(|w| !w.starts_with('-')).collect::<Vec<_>>().join(" ").to_lowercase();
    for l in 1.. {
        let key = format!("exported variable {l}");
        match p.string(&key) {
            Some(s) => {
                let have = s.split_whitespace().filter(|w| !w.starts_with('-')).collect::<Vec<_>>().join(" ");
                if have.to_lowercase() == name {
                    return;
                }
            }
            None => {
                p.insert(&key, Value::Str(spec.into()));
                return;
            }
        }
    }
}

/// `Save Viscosity`'s nodal fields (`Viscosity`, `Shearrate`, `Viscosity
/// Weight`): zeroed before each assembly, accumulated element by element.
const SAVED_VISCOSITY: [&str; 3] = ["viscosity", "shearrate", "viscosity weight"];

fn add_new(p: &mut ParamMap, key: &str, v: Value) {
    if p.get(key).is_none() {
        p.insert(key, v);
    }
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

/// `Pressure Stabilization`, or its synonym.
fn pressure_stabilization(p: &ParamMap) -> bool {
    flag(p, "pressure stabilization") || stabilize_flag(p)
}

/// The keywords `IncompressibleNSSolver_Init0` and `_init` add to the solver
/// section (`ListAddNew*`: a keyword already there wins). The `Variable`
/// keyword and the `Element` string must be in the section before the mesh
/// and the solver variable are created, so a caller building the model calls
/// this before [`Simulation::setup`](elmer_solver::Simulation::setup);
/// `init` repeats it for the keywords that still take effect afterwards.
pub fn solver_defaults(p: &mut ParamMap, dim: usize, transient: bool) {
    solver_defaults_with(p, dim, transient, true);
}

/// [`solver_defaults`] with the simulation's `Serendipity P Elements`
/// (default true), which picks the bubble count of the MINI element
/// `_Init0` adds without stabilization.
pub fn solver_defaults_with(p: &mut ParamMap, dim: usize, transient: bool, serendipity: bool) {
    if pressure_stabilization(p) {
        add_new(p, "Element", Value::Str("n:1".into()));
    } else {
        let mini = if serendipity {
            "p:1 -tri b:1 -tetra b:1 -quad b:3 -brick b:4 -prism b:4 -pyramid b:4"
        } else {
            "p:1 -tri b:1 -tetra b:1 -quad b:4 -brick b:8 -prism b:9 -pyramid b:4"
        };
        add_new(p, "Element", Value::Str(mini.into()));
    }
    // _Init: a bubble element gets its tetra and prism rules.
    if p.string("element").is_some_and(|e| e.contains("b:")) {
        add_new(p, "Element Integration Points", Value::Str("-tetra 24 -prism 85".into()));
    }
    add_new(
        p,
        "Variable",
        Value::Str(format!("Flow Solution[Velocity:{dim} Pressure:1]")),
    );
    add_new(p, "Nonlinear System Norm DOFs", Value::Integer(vec![dim as i32]));
    add_new(p, "Relative Pressure Relaxation", Value::Logical(true));
    add_new(p, "GradP Discretization", Value::Logical(false));
    add_new(p, "Div-Curl Discretization", Value::Logical(false));
    add_new(p, "Bubbles in Global System", Value::Logical(false));
    if let Some(m) = p.string("flow model").map(elmer_core::normalize_key) {
        if m == "stokes" {
            add_new(p, "Stokes Flow", Value::Logical(true));
        }
    }
    // _Init: Save Viscosity exports the viscosity, the shear rate and the
    // weight that normalizes the former.
    if flag(p, "save viscosity") {
        for name in ["Viscosity", "Shearrate", "-nooutput Viscosity Weight"] {
            add_exported(p, name);
        }
    }
    // _Init: a transient Navier-Stokes run with condensed bubbles takes at
    // least two nonlinear iterations.
    let bubbles = p.string("element").is_some_and(|e| e.contains("b:"));
    if bubbles && transient && !flag(p, "stokes flow") && !flag(p, "bubbles in global system") {
        add_new(p, "Nonlinear System Min Iterations", Value::Integer(vec![2]));
        add_new(p, "Nonlinear System Max Iterations", Value::Integer(vec![2]));
    }
}

/// The incompressible stabilized equal-order solver.
pub struct IncompressibleNsSolver {
    basis: Option<Box<dyn ElementBasis + Send>>,
    quad: Option<Box<dyn Quadrature + Send>>,
    local: LocalSystem,
    /// The condensed bubble velocities per element (`BubbleHistorySetCur`),
    /// `dim` per bubble, and those of the previous timestep.
    bubbles: std::collections::BTreeMap<u32, Vec<f64>>,
    bubbles_prev: std::collections::BTreeMap<u32, Vec<f64>>,
    /// The timestep `bubbles_prev` belongs to.
    bubble_step: Option<i64>,
}

impl Default for IncompressibleNsSolver {
    /// A solver on the `elmer-fe` basis and Gauss rules.
    fn default() -> Self {
        Self {
            basis: None,
            quad: None,
            local: LocalSystem::zeros(0),
            bubbles: Default::default(),
            bubbles_prev: Default::default(),
            bubble_step: None,
        }
    }
}

impl IncompressibleNsSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis: Some(basis),
            quad: Some(quad),
            ..Self::default()
        }
    }
}

/// Solver-section options this crate does not port. `pdofs`: the solver
/// has p element dofs (`SolverCtx::solver_dofs`).
fn check_solver_params(p: &ParamMap, pdofs: bool) -> Result<(), SolverError> {
    let element = p.string("element").map(|e| e.trim().to_ascii_lowercase());
    if pressure_stabilization(p) {
        // `p:1` from SolverDofs is the nodal element as well.
        if element.as_deref().is_some_and(|e| e != "n:1" && !(pdofs && e == "p:1")) {
            return Err(unsupported(format!(
                "Element = \"{}\" with Pressure Stabilization (only the nodal \"n:1\" or \"p:1\" element)",
                element.unwrap_or_default()
            )));
        }
    } else {
        // p elements from SolverDofs: Taylor-Hood (the pressure pinned to the
        // nodal one) and the MINI element (velocity bubbles, condensed).
        let Some(e) = element.filter(|e| pdofs && !e.contains("n:1")) else {
            return Err(unsupported(
                "p elements without p dofs (SolverDofs): the default MINI element (condensed velocity \
                 bubbles), Taylor-Hood; only Pressure Stabilization (Element = \"n:1\") is ported then",
            ));
        };
        if !e.starts_with("p:") || (e.split_whitespace().next() == Some("p:1") && !e.contains("b:")) {
            return Err(unsupported(format!("Element = \"{e}\" without Pressure Stabilization")));
        }
    }
    for (key, what) in [
        ("newtonian viscosity condition", "Newtonian Viscosity Condition"),
        ("constant-viscosity start", "Constant-Viscosity Start"),
        ("block preconditioner", "the block preconditioner"),
        ("nonlinear system linesearch", "nonlinear line search"),
        ("save slip", "Save Slip"),
        ("bubbles in global system", "bubbles in the global system"),
        ("create schur matrix approximation", "Schur matrix approximation"),
    ] {
        if flag(p, key) {
            return Err(unsupported(what));
        }
    }
    if has_prefix(p, "schur operator") {
        return Err(unsupported("Schur operators"));
    }
    if p.get("implicit friction coefficient").is_some() {
        return Err(unsupported("implicit friction"));
    }
    if let Some(s) = p.string("pressure stabilization element size") {
        parse_hmode(s)?;
    }
    IntegrationPoints::from_params(p)?;
    Ok(())
}

/// The point counts `GaussPointsAdapt` takes from the solver section for
/// bulk elements: `Element Integration Points` (one count per element
/// family, `ParseElementalGaussRules`) over `Number of Integration Points`.
/// Zero leaves the element its own rule.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct IntegrationPoints {
    base: usize,
    family: [usize; 9],
}

impl IntegrationPoints {
    fn from_params(p: &ParamMap) -> Result<Self, SolverError> {
        for (key, what) in [
            ("relative integration order", "Relative Integration Order"),
            ("bulk relative integration order", "Bulk Relative Integration Order"),
            ("element relative integration order", "Element Relative Integration Order"),
            ("adaptive integration variable", "adaptive integration"),
        ] {
            if p.get(key).is_some() {
                return Err(unsupported(what));
            }
        }
        let mut out = Self {
            base: p.integer("number of integration points").filter(|&n| n > 0).map_or(0, |n| n as usize),
            ..Self::default()
        };
        if let Some(def) = p.string("element integration points") {
            let def = def.to_ascii_lowercase();
            let names = ["-line", "-tri", "-quad", "-tetra", "-pyramid", "-prism", "-brick"];
            for (k, name) in names.iter().enumerate() {
                let Some(j) = def.find(name) else { continue };
                let rest = def[j + name.len()..].trim_start();
                let num: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
                out.family[k + 2] = num.parse().map_err(|_| {
                    SolverError::fatal("ParseElementalGaussRules", format!("Problems reading integer from: {rest}"))
                })?;
            }
        }
        Ok(out)
    }

    /// The count for an element of type `etype`, `None` for its own rule.
    fn for_type(&self, etype: elmer_core::ElementType) -> Option<usize> {
        let n = match self.family.get(etype.0 as usize / 100) {
            Some(&n) if n > 0 => n,
            _ => self.base,
        };
        (n > 0).then_some(n)
    }
}

fn parse_hmode(s: &str) -> Result<HMode, SolverError> {
    match elmer_core::normalize_key(s).as_str() {
        "shortest edge" => Ok(HMode::Default),
        "longest edge" => Ok(HMode::LongestEdge),
        "volume" => Ok(HMode::Volume),
        other => Err(SolverError::fatal(
            CALLER,
            format!("Unknown \"Pressure Stabilization Element Size\": {other}. Use \"shortest edge\", \"longest edge\" or \"volume\""),
        )),
    }
}

/// Material and body force options this crate does not port.
fn check_body(mat: Option<&ParamMap>, bf: Option<&ParamMap>) -> Result<(), SolverError> {
    if let Some(m) = mat {
        if let Some(c) = m.string("compressibility model") {
            let c = elmer_core::normalize_key(c);
            if !["incompressible", "thermal", "perfect gas", "perfect gas equation 1"].contains(&c.as_str()) {
                return Err(unsupported(format!("Compressibility Model = {c}")));
            }
        }
        if let Some(v) = m.string("viscosity model") {
            let v = elmer_core::normalize_key(v);
            if !["power law", "carreau", "cross", "powell eyring", "glen"].contains(&v.as_str()) {
                return Err(unsupported(format!(
                    "Viscosity Model = {v} (only power law, carreau, cross, powell eyring, glen)"
                )));
            }
            // Fatal in the pristine solver ("Parameter > Basis < is required").
            for key in ["nominal shear rate", "yasuda exponent"] {
                if m.get(key).is_some() {
                    return Err(unsupported(format!("{key} of a Viscosity Model")));
                }
            }
        } else if m.get("viscosity model").is_some() {
            return Err(unsupported("Viscosity Model given as other than a string"));
        }
        if m.get("turbulence model").is_some() {
            return Err(unsupported("turbulence model viscosity"));
        }
    }
    if let Some(b) = bf {
        if flag(b, "boussinesq") {
            return Err(unsupported("Boussinesq body force"));
        }
    }
    Ok(())
}

/// The material's `Compressibility Model`, `"thermal"` or `"perfect gas"`
/// (also `Perfect Gas Equation 1`); `None` for an incompressible one.
fn compress_model(mat: Option<&ParamMap>) -> Option<&'static str> {
    match mat.and_then(|m| m.string("compressibility model")).map(elmer_core::normalize_key).as_deref() {
        Some("thermal") => Some("thermal"),
        Some("perfect gas" | "perfect gas equation 1") => Some("perfect gas"),
        _ => None,
    }
}

/// The material's non-newtonian `Viscosity Model` (`EffectiveViscosityVec`),
/// its parameters element constants.
fn viscosity_model(ctx: &SolverCtx<'_>, m: &ParamMap, e: ElemId) -> Result<Option<ViscosityModel>, SolverError> {
    let Some(v) = m.string("viscosity model").map(elmer_core::normalize_key) else { return Ok(None) };
    let c = |key: &str| element_constant(ctx, m, key, e);
    Ok(Some(match v.as_str() {
        "power law" => ViscosityModel::PowerLaw {
            exponent: c("viscosity exponent")?,
            critical: match m.get("critical shear rate") {
                Some(_) => Some(c("critical shear rate")?),
                None => None,
            },
        },
        "carreau" => ViscosityModel::Carreau {
            difference: c("viscosity difference")?,
            exponent: c("viscosity exponent")?,
            transition: c("viscosity transition")?,
        },
        "cross" => ViscosityModel::Cross {
            difference: c("viscosity difference")?,
            exponent: c("viscosity exponent")?,
            transition: c("viscosity transition")?,
        },
        "powell eyring" => ViscosityModel::PowellEyring {
            difference: c("viscosity difference")?,
            transition: c("viscosity transition")?,
        },
        "glen" => {
            if m.get("glen enhancement factor function").is_some() {
                return Err(SolverError::fatal(CALLER, "No Glen function API yet!"));
            }
            let or = |key: &str, default: f64| -> Result<f64, SolverError> {
                match m.get(key) {
                    Some(_) => c(key),
                    None => Ok(default),
                }
            };
            let set_arrhenius = m.logical("set arrhenius factor").unwrap_or(false);
            if !set_arrhenius && m.get("relative temperature").is_none() {
                return Err(SolverError::fatal(CALLER, ">Relative Temperature< not given for viscosity model \"glen\""));
            }
            ViscosityModel::Glen(GlenParams {
                exponent: or("glen exponent", 3.0)?,
                critical: match m.get("critical shear rate") {
                    Some(_) => Some(c("critical shear rate")?),
                    None => None,
                },
                arrhenius: match set_arrhenius {
                    true => Some(c("arrhenius factor")?),
                    false => None,
                },
                t_limit: or("limit temperature", -10.0)?,
                rate: [or("rate factor 1", 3.985e-13)?, or("rate factor 2", 1.916e3)?],
                energy: [or("activation energy 1", 60.0e3)?, or("activation energy 2", 139.0e3)?],
                gas_constant: ctx.sim.model.constants.real("gas constant").unwrap_or(8.314),
            })
        }
        other => return Err(unsupported(format!("Viscosity Model = {other}"))),
    }))
}

/// `ListGetElementReal` of a keyword constant over the element (zero when
/// absent).
fn element_constant(ctx: &SolverCtx<'_>, p: &ParamMap, key: &str, e: ElemId) -> Result<f64, SolverError> {
    let v = ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0]);
    if v.iter().any(|&x| x != v[0]) {
        return Err(unsupported(format!("{key} varying within an element")));
    }
    Ok(v[0])
}

/// Boundary options this crate does not port.
fn check_bc(bc: &ParamMap, var_name: &str) -> Result<(), SolverError> {
    for (key, what) in [
        ("normal-tangential velocity", "normal-tangential velocity"),
        ("wall law", "wall law"),
    ] {
        if flag(bc, key) {
            return Err(unsupported(what));
        }
    }
    if flag(bc, &format!("normal-tangential {var_name}")) {
        return Err(unsupported("normal-tangential velocity"));
    }
    for (prefix, what) in [
        ("normal slip coefficient", "Normal Slip Coefficient"),
        ("friction coefficient", "Friction Coefficient"),
        ("weertman friction coefficient", "Weertman friction"),
        ("fssa", "FSSA"),
        ("implicit friction", "implicit friction"),
    ] {
        if has_prefix(bc, prefix) {
            return Err(unsupported(what));
        }
    }
    for k in 1..=3 {
        if bc.get(&format!("pressure {k}")).is_some() {
            return Err(unsupported(format!("Pressure {k} (use Surface Traction {k})")));
        }
    }
    Ok(())
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

fn newton_active(ctx: &SolverCtx<'_>) -> bool {
    let id = ctx.solver_id();
    ctx.sim.slots.iter().find(|s| s.id == id).is_some_and(|s| s.state.newton_active)
}

/// Parent whose centre orients a boundary normal (`CheckNormalDirection`):
/// the only parent, or at a body-body boundary the one with the higher body.
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

/// The basis, rule, node order and dof count of bulk element `e` from the
/// solver's p dofs (`ElementInfo`, `GaussPointsAdapt`, `GetElementNOFDOFs`).
struct PElement {
    basis: Box<dyn ElementBasis + Send + Sync>,
    rule: Vec<elmer_core::QuadPoint>,
    nodes: Vec<u32>,
    nd: usize,
    /// Local bubbles (`GetElementNOFBDOFs`), condensed.
    nb: usize,
}

/// `values` (at `mesh_nodes`) at `nodes` instead.
fn reorder<T: Copy>(values: Vec<T>, mesh_nodes: &[u32], nodes: &[u32]) -> Vec<T> {
    if mesh_nodes == nodes {
        return values;
    }
    nodes
        .iter()
        .map(|n| values[mesh_nodes.iter().position(|m| m == n).expect("same node set")])
        .collect()
}

fn p_element(d: &dyn elmer_solver::SolverDofs, e: ElemId, np: Option<usize>) -> Result<PElement, SolverError> {
    let mut slots = Vec::new();
    d.element_slots(e, &mut slots);
    Ok(PElement {
        basis: d.element_basis(e).map_err(|m| SolverError::fatal("ElementInfo", m))?,
        rule: match np {
            Some(np) => d.element_rule_np(e, np),
            None => d.element_rule(e),
        }
        .map_err(|m| SolverError::fatal("GaussPoints", m))?,
        nodes: d.element_nodes(e),
        nd: slots.len(),
        nb: d.local_bubbles(e),
    })
}

fn node_xyz(mesh: &dyn elmer_core::MeshView, nodes: &[u32]) -> Vec<[f64; 3]> {
    nodes.iter().map(|&n| mesh.coord(n)).collect()
}

/// `GetElementFamily`.
fn family(etype: elmer_core::ElementType) -> usize {
    etype.family() as usize
}

fn component_names(ctx: &SolverCtx<'_>) -> Result<Vec<String>, SolverError> {
    let spec = ctx
        .params()
        .string("variable")
        .ok_or_else(|| SolverError::fatal(CALLER, "solver has no 'Variable'"))?;
    Ok(parse_variable_spec(spec, ctx.params().integer("variable dofs"))?.components)
}

/// The Dirichlet rows of the solver variable from its component names
/// (`DefaultDirichletBCs`); `Noslip Wall BC` fixes the velocity components
/// to zero, as the solver adds `Velocity i = 0` to those BCs. On p elements
/// a component the BC sets has its edge and face dofs fixed to zero, the
/// `LocalBcBDOFs` projection of a nodally interpolated value.
fn dirichlet_rows(ctx: &SolverCtx<'_>, comps: &[String], dim: usize) -> Result<Vec<(u32, f64)>, SolverError> {
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
        let wall = flag(bc, "noslip wall bc");
        for (k, key) in comps.iter().enumerate() {
            if bc.get(&format!("{key} condition")).is_some() {
                return Err(unsupported(format!("Dirichlet mask > {key} Condition <")));
            }
            let vals = if wall && k < dim {
                Some(vec![0.0; nodes.len()])
            } else {
                ctx.real_at_nodes(bc, key, e)?
            };
            let Some(vals) = vals else { continue };
            if !extra.is_empty() {
                if let Some(Value::Dependent(d)) = bc.get(key) {
                    if matches!(d.kind, elmer_core::DependentKind::Procedure { .. }) {
                        return Err(unsupported(format!("p-element Dirichlet value > {key} < given by a procedure")));
                    }
                }
                for &s in &extra {
                    fixed.insert(s * dofs as u32 + k as u32, 0.0);
                }
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

/// `AddEquationBasics` (MainUtils.F90:1241) fills the variable of an
/// `Equation = Navier-Stokes` solver with 1e-6, which initial conditions
/// overwrite; elmer-solver creates it zeroed.
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

fn nodal(ctx: &SolverCtx<'_>, p: Option<&ParamMap>, key: &str, e: ElemId, n: usize) -> Result<Vec<f64>, SolverError> {
    Ok(match p {
        Some(p) => ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0; n]),
        None => vec![0.0; n],
    })
}

impl Solver for IncompressibleNsSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        match ctx.simulation_params().string("coordinate system") {
            None => {}
            Some(s) => {
                let s = elmer_core::normalize_key(s);
                if !s.starts_with("cartesian") && s != "axi symmetric" {
                    return Err(unsupported(format!("Coordinate System = {s} (only Cartesian and Axi Symmetric)")));
                }
            }
        }
        let dim = ctx.mesh.dim();
        let transient = ctx.simulation_params().string("simulation type").is_some_and(|s| {
            elmer_core::normalize_key(s) == "transient"
        });
        {
            let sp = slot_params_mut(ctx);
            let mut copy = sp.clone();
            solver_defaults(&mut copy, dim, transient);
            *sp = copy;
        }
        check_solver_params(ctx.params(), ctx.solver_dofs().is_some())?;
        let dofs = ctx.primary().map(|f| f.dofs);
        if dofs != Some(dim + 1) {
            return Err(unsupported(format!(
                "a flow variable with {dofs:?} dofs in {dim}D (expected velocity components and pressure)"
            )));
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
                check_body(b.material, b.body_force)?;
            }
            if ctx.mesh.element(e).etype.n_nodes() != ctx.mesh.element(e).nodes.len() {
                return Err(unsupported("elements with more nodes than their type"));
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
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        let pdofs = ctx.solver_dofs().is_some();
        check_solver_params(ctx.params(), pdofs)?;
        let p = ctx.params().clone();
        let pstab = pressure_stabilization(&p);
        let dim = ctx.mesh.dim();
        let comps = component_names(ctx)?;
        let nsdofs = comps.len();
        if nsdofs != dim + 1 {
            return Err(unsupported(format!("{nsdofs} flow dofs in {dim}D")));
        }

        let linear_assembly = flag(&p, "linear equation");
        let stokes = flag(&p, "stokes flow");
        let gradp = flag(&p, "gradp discretization");
        let div_curl = flag(&p, "div-curl discretization");
        let specific_load = flag(&p, "specific load");
        // Without Pressure Stabilization (Taylor-Hood) LocalBulkMatrix skips the
        // stabilization terms, which a zero tau does exactly.
        let coeff = match pstab {
            true => ctx.creal(&p, "pressure stabilization coefficient")?.unwrap_or(1.0),
            false => 0.0,
        };
        let hmode = match p.string("pressure stabilization element size") {
            Some(s) => parse_hmode(s)?,
            None => HMode::Default,
        };
        let long_edge = flag(ctx.simulation_params(), "stabilization use longest element edge");
        let save_viscosity = flag(&p, "save viscosity");
        let axisymmetric = ctx
            .simulation_params()
            .string("coordinate system")
            .is_some_and(|s| elmer_core::normalize_key(s) == "axi symmetric");
        let rule_np = IntegrationPoints::from_params(&p)?;

        let max_iter = p.integer("nonlinear system max iterations").unwrap_or(1);
        let mut min_iter = p.integer("nonlinear system min iterations").unwrap_or(1);
        if min_iter > max_iter {
            ctx.warn(CALLER, &format!("\"Nonlinear System Min Iterations\" exceeds max, using: {max_iter}"));
            min_iter = max_iter;
        }

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let var_name = ctx.variable_name().unwrap_or_default().to_string();
        let mut dofs = Vec::new();
        let (Some(basis), Some(quad)) = (self.basis.as_deref(), self.quad.as_deref()) else {
            return Err(SolverError::fatal(CALLER, "solve before init"));
        };

        if transient {
            let step = ctx.scalar("timestep").map_or(-1, |t| t.round() as i64);
            if self.bubble_step != Some(step) {
                self.bubbles_prev = self.bubbles.clone();
                self.bubble_step = Some(step);
            }
        }
        for iter in 1..=max_iter {
            ctx.info(CALLER, "--------------------------------------------------------", 4);
            ctx.info(CALLER, &format!("Nonlinear iteration: {iter:>4}"), 4);
            ctx.info(CALLER, "--------------------------------------------------------", 4);
            ctx.default_initialize()?;
            if save_viscosity {
                for name in SAVED_VISCOSITY {
                    if let Some(v) = ctx.variable_mut(name) {
                        v.values.iter_mut().for_each(|x| *x = 0.0);
                    }
                }
            }
            let newton = newton_active(ctx);

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let pe = match ctx.solver_dofs() {
                    Some(d) => Some(p_element(d, e, rule_np.for_type(el.etype))?),
                    None => None,
                };
                let nf = pe.as_ref().map_or(n, |p| p.nd + p.nb);
                // A p element's basis takes its nodes in the order of its dofs
                // (p tetrahedra are reordered); the nodal data follows them.
                let nodes = pe.as_ref().map_or(el.nodes, |p| &p.nodes[..]);
                let xyz = node_xyz(ctx.mesh, nodes);
                let bp = ctx
                    .body_params(e)
                    .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
                let (mat, bf) = (bp.material, bp.body_force);
                check_body(mat, bf)?;

                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    return Err(SolverError::fatal(CALLER, format!("active element {} has nodes outside the system", e.0)));
                }
                let values = ctx.primary().expect("checked in init").values;
                // GetLocalSolution: the condensed bubbles read as zero.
                let solution: Vec<[f64; 4]> = (0..nf)
                    .map(|k| {
                        std::array::from_fn(|i| match dofs.get(k * nsdofs + i) {
                            Some(&d) if i < nsdofs => values[d as usize],
                            _ => 0.0,
                        })
                    })
                    .collect();

                let density = nodal(ctx, mat, "density", e, n)?;
                if density.iter().any(|&d| d != density[0]) {
                    return Err(unsupported("Density varying within an element"));
                }
                let viscosity = reorder(nodal(ctx, mat, "viscosity", e, n)?, el.nodes, nodes);
                let viscosity_model = match mat {
                    Some(m) => viscosity_model(ctx, m, e)?,
                    None => None,
                };
                // The glen model's rate factor from the temperature.
                let glen_fields = match (&viscosity_model, mat) {
                    (Some(ViscosityModel::Glen(g)), Some(m)) if g.arrhenius.is_none() => {
                        let t = reorder(nodal(ctx, Some(m), "relative temperature", e, n)?, el.nodes, nodes);
                        let ehf = match m.get("glen enhancement factor") {
                            Some(_) => reorder(nodal(ctx, Some(m), "glen enhancement factor", e, n)?, el.nodes, nodes),
                            None => vec![1.0; n],
                        };
                        Some((t, ehf))
                    }
                    _ => None,
                };
                if viscosity_model.is_some() && linear_assembly {
                    return Err(unsupported("Linear Equation with a non-newtonian Viscosity Model"));
                }
                let compress = match (compress_model(mat), mat) {
                    (Some(model), Some(m)) => {
                        if pstab || (transient && stokes) {
                            return Err(unsupported(format!(
                                "Compressibility Model = {model} with Pressure Stabilization or in a transient \
                                 Stokes run"
                            )));
                        }
                        let t_ref = element_constant(ctx, m, "reference temperature", e)?;
                        let model = if model == "thermal" {
                            let beta = element_constant(ctx, m, "heat expansion coefficient", e)?;
                            CompressModel::Thermal { beta, t_ref }
                        } else {
                            if linear_assembly {
                                return Err(unsupported(
                                    "Compressibility Model = Perfect Gas with Linear Equation = True",
                                ));
                            }
                            let gamma = match m.get("specific heat ratio") {
                                Some(_) => element_constant(ctx, m, "specific heat ratio", e)?,
                                None => 5.0 / 3.0,
                            };
                            let cp = element_constant(ctx, m, "heat capacity", e)?;
                            let p_ref = element_constant(ctx, m, "reference pressure", e)?;
                            CompressModel::PerfectGas { gas_constant: (gamma - 1.0) / gamma * cp, p_ref }
                        };
                        let Some(tv) = ctx.variable("temperature") else {
                            return Err(SolverError::fatal(
                                CALLER,
                                "\"Compressibility Model = Thermal\"/\"Perfect Gas\" requires a solved \
                                 \"Temperature\" field",
                            ));
                        };
                        let temp: Vec<f64> = nodes
                            .iter()
                            .map(|&nn| match tv.perm.map_or(nn, |p| p[nn as usize]) {
                                NO_DOF => t_ref,
                                k => tv.values[k as usize * tv.dofs],
                            })
                            .collect();
                        // TempVar % PrevValues(:, 1) at a transient run.
                        let prev = match transient {
                            false => None,
                            true => {
                                let tvar = ctx.sim.vars.get("temperature").expect("looked up above");
                                let Some(pv) = tvar.prev_values.first() else {
                                    return Err(SolverError::fatal(
                                        CALLER,
                                        "a transient compressible run needs the previous Temperature",
                                    ));
                                };
                                let tperm = tvar.perm.as_deref();
                                let p: Vec<f64> = nodes
                                    .iter()
                                    .map(|&nn| match tperm.map_or(nn, |p| p[nn as usize]) {
                                        NO_DOF => t_ref,
                                        k => pv[k as usize * tvar.dofs],
                                    })
                                    .collect();
                                Some(p)
                            }
                        };
                        Some((temp, model, prev))
                    }
                    _ => None,
                };
                let mut load = vec![[0.0; 3]; n];
                if let Some(b) = bf {
                    for i in 0..dim {
                        if let Some(f) = ctx.real_at_nodes(b, &format!("flow bodyforce {}", i + 1), e)? {
                            for (l, v) in load.iter_mut().zip(reorder(f, el.nodes, nodes)) {
                                l[i] = v;
                            }
                        }
                    }
                }

                let (_mk, hk) = elmer_fe::stab_param(el.etype, &xyz, dim, long_edge)
                    .map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
                let inp = IncompressibleBulkInput {
                    etype: el.etype,
                    xyz: &xyz,
                    dim,
                    density: density[0],
                    viscosity: &viscosity,
                    load: &load,
                    solution: if linear_assembly { None } else { Some(&solution) },
                    stokes,
                    newton,
                    want_mass: transient && !stokes,
                    gradp,
                    div_curl,
                    specific_load,
                    pstab: PressureStab { coeff, hmode, hk },
                    np: rule_np.for_type(el.etype),
                    axisymmetric,
                    viscosity_model,
                    glen_fields: glen_fields
                        .as_ref()
                        .map(|(t, ehf)| GlenFields { relative_temperature: t, enhancement: ehf }),
                    compress: compress.as_ref().map(|(temperature, model, prev)| Compressibility {
                        temperature,
                        model: *model,
                        transient: prev.as_deref().map(|p| (p, dt)),
                    }),
                };
                let mut sums = ViscositySums::default();
                let want_sums = save_viscosity.then_some(&mut sums);
                match &pe {
                    None => {
                        let points = quad.gauss_points(el.etype, inp.np);
                        incompressible_bulk_local_sums(basis, &points, n, &inp, &mut self.local, want_sums)?
                    }
                    Some(p) => {
                        incompressible_bulk_local_sums(&*p.basis, &p.rule, nf, &inp, &mut self.local, want_sums)?;
                        // LocalBulkMatrix: a Stokes run with bubbles condenses
                        // them and pins nothing; otherwise the pressure is
                        // pinned to the nodal one first.
                        if !(stokes && p.nb > 0) && !pstab {
                            pin_pressure(&mut self.local, n, nsdofs);
                        }
                        if p.nb > 0 && transient && !stokes {
                            // LCondensate with the bubble history (BDF1, whatever
                            // the timestepping method): in place of
                            // Default1stOrderTime.
                            let prev = &ctx.sim.vars.get(&var_name).expect("the primary variable").prev_values;
                            let Some(prev) = prev.first() else {
                                return Err(SolverError::fatal(CALLER, "transient run without a previous timestep"));
                            };
                            let retained = &dofs[..p.nd * nsdofs];
                            let xprev: Vec<f64> = retained.iter().map(|&d| prev[d as usize]).collect();
                            let xl: Vec<f64> = retained.iter().map(|&d| values[d as usize]).collect();
                            let zero = vec![0.0; dim * p.nb];
                            let bprev = self.bubbles_prev.get(&e.0).unwrap_or(&zero);
                            let cur = lcondensate_bdf1(&mut self.local, p.nd, p.nb, dim, dt, &xprev, bprev, &xl)?;
                            self.bubbles.insert(e.0, cur);
                        } else if p.nb > 0 {
                            lcondensate(&mut self.local, p.nd, p.nb, dim)?;
                        }
                    }
                }
                if save_viscosity {
                    let sums = [&sums.viscosity, &sums.shear_rate, &sums.weight];
                    for (name, vals) in SAVED_VISCOSITY.into_iter().zip(sums) {
                        let Some(v) = ctx.variable_mut(name) else { continue };
                        for (&node, x) in nodes.iter().zip(vals.iter()) {
                            let k = v.perm.map_or(node, |p| p[node as usize]);
                            if k != NO_DOF {
                                v.values[k as usize] += x;
                            }
                        }
                    }
                }
                if transient && !stokes && pe.as_ref().is_none_or(|p| p.nb == 0) {
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
                // Skip 101 elements in 2D, and additionally 202's in 3D.
                if family(el.etype) < dim {
                    continue;
                }
                let bc = ctx.bc_params(e).expect("boundary_elements filters on BC");
                check_bc(bc, &var_name)?;
                let n = el.nodes.len();
                let bc = Some(bc);
                let mut ext = nodal(ctx, bc, "normal surface traction", e, n)?;
                if bc.is_some_and(|b| b.get("normal surface traction").is_none()) {
                    ext = nodal(ctx, bc, "external pressure", e, n)?;
                }
                let mut traction = vec![[0.0; 3]; n];
                let mut slip = vec![[0.0; 3]; n];
                for i in 0..3 {
                    let t = nodal(ctx, bc, &format!("surface traction {}", i + 1), e, n)?;
                    let s = nodal(ctx, bc, &format!("slip coefficient {}", i + 1), e, n)?;
                    for k in 0..n {
                        traction[k][i] = t[k];
                        slip[k][i] = s[k];
                    }
                }
                let any = |v: &[f64]| v.iter().any(|&x| x != 0.0);
                if !any(&ext) && !traction.iter().any(|t| any(t)) && !slip.iter().any(|s| any(s)) {
                    continue;
                }
                let parent_center = normal_parent(ctx, e)
                    .map(|pe| elmer_physics_flow::boundary::element_center(&node_xyz(ctx.mesh, ctx.mesh.element(pe).nodes)));
                // A p boundary element lists its nodes in its parent's edge or
                // face order; its values follow that order.
                let p_nodes = ctx.solver_dofs().map(|d| d.element_nodes(e));
                let nodes = p_nodes.as_deref().unwrap_or(el.nodes);
                let (ext, traction, slip) = (
                    reorder(ext, el.nodes, nodes),
                    reorder(traction, el.nodes, nodes),
                    reorder(slip, el.nodes, nodes),
                );
                let xyz = node_xyz(ctx.mesh, nodes);
                let inp = IncompressibleBoundaryInput {
                    etype: el.etype,
                    xyz: &xyz,
                    dim,
                    traction: &traction,
                    ext_pressure: &ext,
                    slip: &slip,
                    parent_center,
                    np: None,
                    axisymmetric,
                };
                match ctx.solver_dofs() {
                    Some(d) => {
                        let pb = d.element_basis(e).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                        let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
                        let nd = dofs.len() / nsdofs;
                        incompressible_boundary_local_with(&*pb, &rule, nd, &inp, &mut self.local)?;
                    }
                    None => incompressible_boundary_local(basis, quad, &inp, &mut self.local)?,
                }
                ctx.default_update(&self.local, e)?;
            }
            ctx.default_finish_boundary_assembly()?;

            let fixed = dirichlet_rows(ctx, &comps, dim)?;
            apply_dirichlet(ctx, &fixed)?;
            ctx.default_solve()?;

            if ctx.convergence().is_some_and(|c| c.nonlin_converged == 1) && iter >= min_iter {
                break;
            }
        }
        // The saved viscosity as a weighted average where the weight exceeds
        // 1e-20 (the Shearrate is left as summed: the Fortran normalizes a
        // "Strainrate").
        if save_viscosity {
            let weight = ctx.variable("viscosity weight").map(|w| w.values.to_vec());
            if let (Some(w), Some(v)) = (weight, ctx.variable_mut("viscosity")) {
                for (x, &wi) in v.values.iter_mut().zip(&w) {
                    if wi > 1.0e-20 {
                        *x /= wi;
                    }
                }
            }
        }
        Ok(())
    }
}
