//! `FluxSolver`: the `Procedure = "FluxSolver" "FluxSolver"` solver, driving
//! [`flux_bulk_local`] through [`SolverCtx`] the way `FluxSolver` in
//! FluxSolver.F90:44 drives its `BulkAssembly`.

use std::sync::{Arc, OnceLock};
use std::time::Instant;

use elmer_core::{
    normalize_key, ElemId, ElementBasis, FieldKind, LocalSystem, ParamMap, Params, Quadrature,
    Value, NO_DOF,
};
use elmer_solver::sim::node_order_perm;
use elmer_solver::{Perm, Solver, SolverCtx, SolverError, SolverRegistry, SolverSlot, Variable};

use crate::{flux_bulk_local, FluxCoefficient, FluxConfig, FluxError, FluxLocal, PROCEDURE};

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

/// Registers `("FluxSolver", "FluxSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`FluxSolver`]. Fails
/// with `NotWired` until [`install_fe`] has been called.
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    let provider = FE
        .get()
        .ok_or_else(|| SolverError::NotWired("element basis and quadrature (elmer-fe)".into()))?;
    let (basis, quad) = provider();
    Ok(Box::new(FluxSolver::new(basis, quad)))
}

impl From<FluxError> for SolverError {
    fn from(e: FluxError) -> Self {
        match e {
            FluxError::Unsupported { .. } => SolverError::Unsupported(e.to_string()),
            e => SolverError::fatal("FluxSolver", e.to_string()),
        }
    }
}

const CALLER: &str = "FluxSolver";

pub struct FluxSolver {
    basis: Box<dyn ElementBasis + Send>,
    quad: Box<dyn Quadrature + Send>,
    local: FluxLocal,
}

impl FluxSolver {
    pub fn new(basis: Box<dyn ElementBasis + Send>, quad: Box<dyn Quadrature + Send>) -> Self {
        Self {
            basis,
            quad,
            local: FluxLocal::default(),
        }
    }
}

/// The four `Calculate ...` requests and whether `Calculate Flux` was given
/// (`GotIt` after the last `GetLogical`, FluxSolver.F90:90-94).
struct Requests {
    grad_mag: bool,
    grad: bool,
    flux_mag: bool,
    flux: bool,
    flux_given: bool,
}

impl Requests {
    fn read(p: &dyn Params) -> Self {
        Self {
            grad_mag: p.logical("calculate grad magnitude") == Some(true),
            grad: p.logical("calculate grad") == Some(true),
            flux_mag: p.logical("calculate flux magnitude") == Some(true),
            flux: p.logical("calculate flux") == Some(true),
            flux_given: p.logical("calculate flux").is_some(),
        }
    }

    fn any(&self) -> bool {
        self.grad_mag || self.grad || self.flux_mag || self.flux
    }

    fn config(&self, dim: usize) -> FluxConfig {
        FluxConfig {
            calculate_flux: self.flux,
            calculate_flux_magnitude: self.flux_mag,
            calculate_grad: self.grad,
            calculate_grad_magnitude: self.grad_mag,
            ..FluxConfig::new(dim)
        }
    }
}

/// What the solver differentiates: `Flux Variable`, else `Target Variable`
/// (both a field name), else `Target Expression` (a material keyword), else
/// `Temperature` with a warning (FluxSolver.F90:117-128).
struct Target {
    name: String,
    /// `Target Expression`: `name` is a material keyword, not a field.
    expression: bool,
    defaulted: bool,
}

fn target_of(p: &dyn Params) -> Target {
    let s = |k: &str| p.string(k).map(str::to_string);
    if let Some(name) = s("flux variable").or_else(|| s("target variable")) {
        return Target {
            name,
            expression: false,
            defaulted: false,
        };
    }
    if let Some(name) = s("target expression") {
        return Target {
            name,
            expression: true,
            defaulted: false,
        };
    }
    Target {
        name: "Temperature".into(),
        expression: false,
        defaulted: true,
    }
}

fn slot_mut<'a>(ctx: &'a mut SolverCtx<'_>) -> &'a mut SolverSlot {
    let id = ctx.solver_id();
    ctx.sim
        .slots
        .iter_mut()
        .find(|s| s.id == id)
        .expect("the running solver has a slot")
}

fn add_new(p: &mut ParamMap, key: &str, v: Value) {
    if p.get(key).is_none() {
        p.insert(key, v);
    }
}

fn axisymmetric(ctx: &SolverCtx<'_>) -> Result<bool, SolverError> {
    match ctx.simulation_params().string("coordinate system") {
        None => Ok(false),
        Some(s) => {
            let s = normalize_key(s);
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

/// `Flux Coefficient` of a material at the nodes of element `e`, by the shape
/// of the keyword's array: 1 row is a scalar, 1 column a diagonal, else a full
/// tensor (GetTensorRank, FluxSolver.F90:553). An absent keyword is the
/// identity.
fn coefficient(
    ctx: &SolverCtx<'_>,
    material: Option<&ParamMap>,
    key: &str,
    e: ElemId,
    n: usize,
    dim: usize,
) -> Result<FluxCoefficient, SolverError> {
    let Some(v) = material.and_then(|m| m.get(key)) else {
        return Ok(FluxCoefficient::Identity);
    };
    match v {
        Value::Real { rows, cols, data } => {
            let (r, c) = (*rows, *cols);
            if r == 1 {
                Ok(FluxCoefficient::Scalar(vec![data[0]; n]))
            } else if c == 1 {
                if r < dim {
                    return Err(SolverError::fatal(
                        CALLER,
                        format!("\"{key}\" has {r} rows, the problem has {dim} dimensions"),
                    ));
                }
                let mut d = [0.0; 3];
                d[..dim].copy_from_slice(&data[..dim]);
                Ok(FluxCoefficient::Diagonal(vec![d; n]))
            } else {
                if r < dim || c < dim {
                    return Err(SolverError::fatal(
                        CALLER,
                        format!("\"{key}\" is {r}x{c}, the problem has {dim} dimensions"),
                    ));
                }
                let mut k = [[0.0; 3]; 3];
                for (i, row) in k.iter_mut().enumerate().take(dim) {
                    for (j, kij) in row.iter_mut().enumerate().take(dim) {
                        *kij = data[i * c + j];
                    }
                }
                Ok(FluxCoefficient::Full(vec![k; n]))
            }
        }
        Value::Dependent(d) if d.rows == 1 && d.cols == 1 => {
            let m = material.expect("a value came from the material");
            let vals = ctx
                .real_at_nodes(m, key, e)?
                .expect("the keyword is present");
            Ok(FluxCoefficient::Scalar(vals))
        }
        _ => Err(SolverError::Unsupported(format!(
            "\"{key}\" as a tensor-valued dependent property"
        ))),
    }
}

impl Solver for FluxSolver {
    /// `FluxSolver_Init` (FluxSolver.F90:580): requests defaults, creates the
    /// hidden `-nooutput <Equation>_temp` solver variable when `Variable` is
    /// absent, exports the requested fields as variables, and sets the zero
    /// time derivative order and the CG/ILU0 linear-system defaults. The
    /// framework has no hook ahead of variable creation, so the variable and
    /// the time order it would have got are settled here.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let dim = ctx.mesh.dim();
        if !(2..=3).contains(&dim) {
            return Err(SolverError::fatal(
                "FluxSolver_init",
                "Flux computation makes sense only in 2D and 3D",
            ));
        }
        axisymmetric(ctx)?;
        if ctx.params().logical("discontinuous galerkin") == Some(true) {
            return Err(FluxError::Unsupported {
                what: "Discontinuous Galerkin = True".into(),
                reason: "elemental DG solver variables and face terms are not ported",
            }
            .into());
        }

        let mut req = Requests::read(ctx.params());
        if !req.any() {
            if req.flux_given {
                ctx.warn("FluxSolver_init", "No field computation requested!");
                return Ok(());
            }
            ctx.warn(
                "FluxSolver_init",
                "No field computation requested, setting > Calculate Flux = True <",
            );
            req.flux = true;
            slot_mut(ctx)
                .params
                .insert("Calculate Flux", Value::Logical(true));
        }

        let target = target_of(ctx.params());
        if target.defaulted {
            ctx.warn(CALLER, "> Target Variable < not given, using Temperature");
        }
        let var = normalize_key(&target.name);

        // The solver's own variable carries the mass-matrix system.
        if ctx.params().get("variable").is_none() {
            let eq = ctx
                .params()
                .string("equation")
                .ok_or_else(|| SolverError::fatal("ListGetString", "Failed to find \"Equation\""))?
                .to_string();
            let perm = solver_perm(ctx)?;
            let mut v = Variable::new(
                &format!("{eq}_temp"),
                FieldKind::Nodal,
                1,
                Some(Arc::from(perm.perm.clone())),
                vec![0.0; perm.n_slots],
            );
            v.output = false;
            let perm_arc = v.perm.clone();
            let id = ctx.sim.vars.add(v)?;
            let slot = slot_mut(ctx);
            slot.variable = Some(id);
            slot.perm = perm_arc;
        }
        {
            let slot = slot_mut(ctx);
            let id = slot
                .variable
                .ok_or_else(|| SolverError::fatal("FluxSolver_init", "solver has no variable"))?;
            let p = &mut slot.params;
            p.insert("Time derivative order", Value::Integer(vec![0]));
            add_new(p, "Linear System Solver", Value::Str("Iterative".into()));
            add_new(p, "Linear System Iterative Method", Value::Str("cg".into()));
            add_new(
                p,
                "Linear System Preconditioning",
                Value::Str("ILU0".into()),
            );
            add_new(p, "Linear System Max Iterations", Value::Integer(vec![500]));
            add_new(p, "Linear System Residual Output", Value::Integer(vec![10]));
            add_new(
                p,
                "Linear System Convergence Tolerance",
                Value::real(1.0e-10),
            );
            slot.time_order = 0;
            ctx.sim.vars.by_id_mut(id).prev_values.clear();
        }

        let perm = slot_mut(ctx).perm.clone();
        let n_slots = {
            let id = slot_mut(ctx).variable.expect("checked above");
            let v = ctx.sim.vars.by_id(id);
            if v.dofs != 1 {
                return Err(SolverError::Unsupported(format!(
                    "FluxSolver variable \"{}\" with {} dofs",
                    v.name, v.dofs
                )));
            }
            v.values.len()
        };

        let export = |ctx: &mut SolverCtx<'_>, name: String, dofs: usize| {
            export_variable(ctx, &name, dofs, perm.clone(), n_slots)
        };
        if req.flux {
            export(ctx, format!("{var} flux"), dim)?;
            if ctx.params().logical("calculate flux abs") == Some(true) {
                export(ctx, format!("{var} flux_abs"), 1)?;
            }
        }
        if req.flux_mag {
            export(ctx, format!("{var} flux_mag"), 1)?;
        }
        if req.grad {
            export(ctx, format!("{var} grad"), dim)?;
            if ctx.params().logical("calculate grad abs") == Some(true) {
                export(ctx, format!("{var} grad_abs"), 1)?;
            }
        }
        if req.grad_mag {
            export(ctx, format!("{var} grad_mag"), 1)?;
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "-------------------------------------", 4);
        ctx.info(CALLER, "Computing the flux and/or gradient", 4);
        ctx.info(CALLER, "-------------------------------------", 4);

        let dim = ctx.mesh.dim();
        let Some(primary) = ctx.primary() else {
            return Ok(());
        };
        // COUNT( Solver % Variable % Perm > 0 ) <= 0
        if primary.perm.is_some_and(|p| p.iter().all(|&s| s == NO_DOF)) {
            return Ok(());
        }
        let req = Requests::read(ctx.params());
        let mut cfg = req.config(dim);
        let dofs = cfg.dofs();
        if dofs == 0 {
            ctx.warn(CALLER, "No field computation requested, exiting...");
            return Ok(());
        }
        cfg.axisymmetric = axisymmetric(ctx)?;
        let first_mag = cfg.first_magnitude();

        // Only the last component is asked for a nonlinear change, so the
        // coupled system norm converges on it.
        let auto_norm = ctx.params().get("skip compute nonlinear change").is_none();
        if auto_norm {
            slot_mut(ctx)
                .params
                .insert("Skip Compute Nonlinear Change", Value::Logical(true));
        }

        let target = target_of(ctx.params());
        if target.defaulted {
            ctx.warn(CALLER, "> Target Variable < not given, using Temperature");
        }
        let var = normalize_key(&target.name);
        let cond_name = if req.flux || req.flux_mag {
            match ctx.params().string("flux coefficient") {
                Some(s) => s.to_string(),
                None => {
                    ctx.warn(
                        CALLER,
                        "> Flux Coefficient < not given, using Heat Conductivity",
                    );
                    "Heat Conductivity".into()
                }
            }
        } else {
            String::new()
        };

        // The target field, copied so the assembly loop can borrow `ctx`.
        let source = if target.expression {
            None
        } else {
            let v = ctx.sim.vars.get(&var).ok_or_else(|| {
                SolverError::fatal(CALLER, format!("Variable not found: {}", target.name))
            })?;
            if v.dofs != 1 || v.kind != FieldKind::Nodal {
                return Err(SolverError::Unsupported(format!(
                    "target variable \"{}\" must be a scalar nodal field",
                    v.name
                )));
            }
            Some((v.values.clone(), v.perm.clone()))
        };

        let t0 = Instant::now();
        ctx.default_initialize()?;
        let n_rows = ctx.system().map_or(0, |s| s.rhs.len());
        let mut forces: Vec<Vec<f64>> = vec![vec![0.0; n_rows]; dofs];

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut slots = Vec::new();
        for &e in &elems {
            let el = ctx.mesh.element(e);
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let n = xyz.len();

            let (potential, coeff) = {
                let bp = ctx.body_params(e);
                let material = bp.as_ref().and_then(|b| b.material);
                let potential = match &source {
                    None => {
                        let m = material.ok_or_else(|| {
                            SolverError::fatal(
                                CALLER,
                                format!("element {} has no material for \"{}\"", e.0, target.name),
                            )
                        })?;
                        ctx.real_at_nodes(m, &target.name, e)?.ok_or_else(|| {
                            SolverError::fatal(
                                "ListGetReal",
                                format!("Failed to find \"{}\" in the material", target.name),
                            )
                        })?
                    }
                    Some((values, perm)) => el
                        .nodes
                        .iter()
                        .map(|&nd| {
                            let s = match perm {
                                Some(p) => p.get(nd as usize).copied().filter(|&s| s != NO_DOF),
                                None => Some(nd),
                            };
                            s.map_or(0.0, |s| values[s as usize])
                        })
                        .collect(),
                };
                let coeff = if req.flux || req.flux_mag {
                    coefficient(ctx, material, &cond_name, e, n, dim)?
                } else {
                    FluxCoefficient::Identity
                };
                (potential, coeff)
            };

            flux_bulk_local(
                &*self.basis,
                &*self.quad,
                el.etype,
                &xyz,
                &potential,
                &coeff,
                &cfg,
                &mut self.local,
            )?;

            let ls = LocalSystem {
                n,
                stiff: self.local.mass.clone(),
                mass: None,
                damp: None,
                force: self.local.force[0].clone(),
            };
            ctx.default_update(&ls, e)?;
            ctx.element_slots(e, &mut slots);
            for (i, f) in forces.iter_mut().enumerate() {
                for (p, &s) in slots.iter().enumerate() {
                    if s != NO_DOF {
                        f[s as usize] += self.local.force[i][p];
                    }
                }
            }
        }
        ctx.info(
            CALLER,
            &format!("Assembly Time: {}", t0.elapsed().as_secs_f64()),
            5,
        );

        let t1 = Instant::now();
        let enforce_positive = ctx.params().logical("enforce positive magnitude") == Some(true);
        let save_rhs = ctx.system().map(|s| s.rhs.clone());
        let mut fields: Vec<Vec<f64>> = Vec::with_capacity(dofs);
        let mut norm = 0.0;
        for (i, f) in forces.into_iter().enumerate() {
            if let Some(s) = ctx.system_mut() {
                s.rhs = f;
            }
            if i + 1 == dofs && auto_norm {
                slot_mut(ctx)
                    .params
                    .0
                    .remove("skip compute nonlinear change");
            }
            norm = ctx.default_solve()?;
            let mut sol = ctx.primary().map(|p| p.values.to_vec()).unwrap_or_default();
            if enforce_positive && i >= first_mag {
                sol.iter_mut().for_each(|v| *v = v.max(0.0));
            }
            fields.push(sol);
        }
        if let (Some(s), Some(r)) = (ctx.system_mut(), save_rhs) {
            s.rhs = r;
        }

        // Only the last component is asked for a nonlinear change, so `norm`
        // is the norm of that one alone (FluxSolver.F90:241-248).
        if let Some(id) = slot_mut(ctx).variable {
            ctx.sim.vars.by_id_mut(id).conv.norm = norm;
        }

        let mut k = 0;
        if req.flux {
            store_vector(ctx, &format!("{var} flux"), dim, &fields[k..k + dim]);
            k += dim;
        }
        if req.grad {
            store_vector(ctx, &format!("{var} grad"), dim, &fields[k..k + dim]);
            k += dim;
        }
        if req.flux_mag {
            store_scalar(ctx, &format!("{var} flux_mag"), &fields[k]);
            k += 1;
        }
        if req.grad_mag {
            store_scalar(ctx, &format!("{var} grad_mag"), &fields[k]);
        }

        // Flux_abs and Grad_abs exist only when requested at init.
        if let Some(v) = ctx.sim.vars.get_mut(&format!("{var} flux_abs")) {
            absolute(&mut v.values, &fields[..dim]);
        }
        if let Some(v) = ctx.sim.vars.get_mut(&format!("{var} grad_abs")) {
            let k = if req.flux { dim } else { 0 };
            absolute(&mut v.values, &fields[k..k + dim]);
        }

        ctx.info(
            CALLER,
            &format!("Solution Time: {}", t1.elapsed().as_secs_f64()),
            5,
        );
        ctx.info(CALLER, &format!("Result Norm: {norm:?}"), 4);
        ctx.info(CALLER, "All done", 6);
        ctx.info(CALLER, "-------------------------------------", 6);
        Ok(())
    }
}

/// The permutation of the solver's own variable, as `Simulation::setup` builds
/// it for a `Variable` keyword.
fn solver_perm(ctx: &mut SolverCtx<'_>) -> Result<Perm, SolverError> {
    let mesh = ctx.mesh;
    let id = ctx.solver_id();
    let slot = ctx
        .sim
        .slots
        .iter()
        .find(|s| s.id == id)
        .expect("the running solver has a slot");
    match &ctx.sim.backends.assembly {
        Some(b) => b.build_perm_for(mesh, &slot.active, &slot.params),
        None => Ok(node_order_perm(mesh, &slot.active)),
    }
}

/// `Exported Variable = Name[Name:dofs]` on the solver's permutation. An
/// existing variable of that name is reused.
fn export_variable(
    ctx: &mut SolverCtx<'_>,
    name: &str,
    dofs: usize,
    perm: Option<Arc<[u32]>>,
    n_slots: usize,
) -> Result<(), SolverError> {
    if let Some(v) = ctx.sim.vars.get(name) {
        if v.dofs != dofs || v.values.len() != n_slots * dofs {
            return Err(SolverError::fatal(
                "FluxSolver_init",
                format!("variable \"{name}\" exists with another size"),
            ));
        }
        return Ok(());
    }
    let v = Variable::new(
        name,
        FieldKind::Nodal,
        dofs,
        perm,
        vec![0.0; n_slots * dofs],
    );
    if dofs > 1 {
        let comps: Vec<String> = (1..=dofs).map(|k| format!("{} {k}", v.name)).collect();
        ctx.sim.vars.add_with_components(v, &comps)?;
    } else {
        ctx.sim.vars.add(v)?;
    }
    Ok(())
}

fn store_vector(ctx: &mut SolverCtx<'_>, name: &str, dim: usize, comps: &[Vec<f64>]) {
    for (c, sol) in comps.iter().enumerate().take(dim) {
        if let Some(mut s) = ctx.component_mut(&format!("{name} {}", c + 1)) {
            for (j, &x) in sol.iter().enumerate().take(s.len()) {
                s.set(j, x);
            }
        }
    }
}

fn store_scalar(ctx: &mut SolverCtx<'_>, name: &str, sol: &[f64]) {
    if let Some(mut s) = ctx.component_mut(name) {
        for (j, &x) in sol.iter().enumerate().take(s.len()) {
            s.set(j, x);
        }
    }
}

/// `sqrt(sum_c field_c^2)` per slot (FluxSolver.F90:252-256).
fn absolute(out: &mut [f64], comps: &[Vec<f64>]) {
    for (j, o) in out.iter_mut().enumerate() {
        *o = comps.iter().map(|c| c[j] * c[j]).sum::<f64>().sqrt();
    }
}
