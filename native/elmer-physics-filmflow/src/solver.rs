//! The `FilmFlowSolver` procedure: `FilmFlowSolver_Init0`, `_Init` and the
//! solver loop of FilmFlowSolver.F90:154.

// Component loops follow the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{normalize_key, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::kernel::{
    local_boundary_matrix, local_bulk_matrix, ns_condensate, BulkInput, BulkOptions, Friction,
};
use crate::FilmFlowError;

/// `Procedure = "FilmFlowSolver" "FilmFlowSolver"`.
pub const PROCEDURE: (&str, &str) = ("FilmFlowSolver", "FilmFlowSolver");

const CALLER: &str = "FilmFlowSolver";

impl From<FilmFlowError> for SolverError {
    fn from(e: FilmFlowError) -> Self {
        SolverError::fatal(CALLER, e.to_string())
    }
}

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(FilmFlowSolver::default()))
}

fn add_new(p: &mut ParamMap, key: &str, v: Value) {
    if p.get(key).is_none() {
        p.insert(key, v);
    }
}

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// `NextFreeKeyword('Exported Variable ', Params)`.
fn next_exported(p: &ParamMap) -> String {
    let mut k = 1;
    while p.get(&format!("exported variable {k}")).is_some() {
        k += 1;
    }
    format!("Exported Variable {k}")
}

/// The keywords `FilmFlowSolver_Init0` and `FilmFlowSolver_Init` add to the
/// solver section (`ListAddNew*`: a keyword already there wins). The
/// `Element` and `Variable` must be in the section before the mesh and the
/// variables are created, so a caller building the model calls this before
/// [`Simulation::setup`](elmer_solver::Simulation::setup). A missing `Model
/// Dimension` is left for the solver's `init` to refuse.
pub fn solver_defaults(p: &mut ParamMap, serendipity: bool) {
    add_new(p, "Time Derivative Order", Value::Integer(vec![1]));
    let element = if serendipity {
        "p:1 -line b:1 -tri b:1 -quad b:3"
    } else {
        "p:1 -line b:1 -tri b:1 -quad b:4"
    };
    add_new(p, "Element", Value::Str(element.into()));
    let mdim = p.integer("model dimension");
    match mdim {
        Some(2) => add_new(
            p,
            "Variable",
            Value::Str("Flow[FilmVelocity:2 FilmPressure:1]".into()),
        ),
        Some(1) => add_new(
            p,
            "Variable",
            Value::Str("Flow[FilmSpeed:1 FilmPressure:1]".into()),
        ),
        _ => return,
    }
    let mdim = mdim.unwrap_or(1);
    add_new(p, "Nonlinear System Norm DOFs", Value::Integer(vec![mdim]));
    add_new(p, "Relative Pressure Relaxation", Value::Logical(true));
    add_new(p, "Bubbles in Global System", Value::Logical(false));
    add_new(p, "Global Mass Matrix", Value::Logical(true));
    let friction = flag(p, "calculate friction heating");
    let pressure = flag(p, "calculate pressure heating");
    if friction {
        let k = next_exported(p);
        p.insert(&k, Value::Str("FilmFriction HeatFlux".into()));
    }
    if pressure {
        let k = next_exported(p);
        p.insert(&k, Value::Str("FilmPressure HeatFlux".into()));
    }
    if friction || pressure {
        let k = next_exported(p);
        p.insert(&k, Value::Str("Heating Energy".into()));
    }
}

/// The state `FilmFlowSolver` keeps between calls (its `SAVE` variables).
#[derive(Default)]
pub struct FilmFlowSolver {
    /// `FirstRound`: cleared after the first nonlinear iteration.
    started: bool,
    /// `DoneWeight` and the heating normalisation weights `HeatingW`.
    heating_w: Option<Vec<f64>>,
    /// `itime0`, the timestep the heating fields were last accumulated in.
    itime0: Option<i64>,
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("FilmFlowSolver: {what}"))
}

/// `ListGetConstReal`-style scalar of a section.
fn const_real(ctx: &SolverCtx<'_>, p: &ParamMap, key: &str) -> Result<Option<f64>, SolverError> {
    ctx.creal(p, key)
}

impl Solver for FilmFlowSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        let mdim = p.integer("model dimension").ok_or_else(|| {
            SolverError::fatal(
                "FilmFlowSolver_init",
                "Give \"Model Dimension\" i.e. the dimension of N-S equation!",
            )
        })?;
        if !(1..=2).contains(&mdim) {
            return Err(SolverError::fatal(
                "FilmFlowSolver_init",
                format!("This module does not make sense in dim: {mdim}"),
            ));
        }
        for (k, what) in [
            ("calculate pressure heating", "Calculate Pressure Heating"),
            ("calculate sensitivity", "Calculate Sensitivity"),
            ("use gap average", "Use Gap Average"),
            ("skip first solution", "Skip First Solution"),
        ] {
            if flag(p, k) {
                return Err(unsupported(what));
            }
        }
        if ctx
            .sim
            .model
            .materials
            .iter()
            .any(|m| m.get("artificial compressibility").is_some())
        {
            return Err(unsupported("Artificial Compressibility"));
        }
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::fatal(
                CALLER,
                "the MINI element of FilmFlowSolver_Init0 is missing; call elmer_physics_filmflow::solver_defaults before the setup",
            ));
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        dt: f64,
        transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(
            CALLER,
            "Computing reduced dimensional Navier-Stokes equations!",
            4,
        );
        let p = ctx.params().clone();
        let mdim = p.integer("model dimension").unwrap_or(1) as usize;
        let dofs = mdim + 1;
        let gradp = flag(&p, "gradp discretization");
        let lateral_strain = flag(&p, "lateral strain");
        let mingap = ctx.creal(&p, "min gap height")?.unwrap_or(1.0e-20);
        let mut convect = flag(&p, "convect");
        let density_name = p.string("density name").unwrap_or("Density").to_string();
        let viscosity_name = p
            .string("viscosity name")
            .unwrap_or("Viscosity")
            .to_string();
        let gap_factor = ctx.creal(&p, "gap addition factor")?;
        if gap_factor.is_some_and(|g| g.abs() > 1.0) {
            ctx.warn(
                CALLER,
                "\"Gap Addition Factor\" greater to unity does not make sense!",
            );
        }
        let csymmetry = match p.logical("axi symmetric") {
            Some(b) => b,
            None => ctx
                .simulation_params()
                .string("coordinate system")
                .is_some_and(|c| {
                    let c = normalize_key(c);
                    c == "axi symmetric" || c == "cylindric symmetric"
                }),
        };
        let friction = match p.string("friction model").map(normalize_key).as_deref() {
            None | Some("laminar") => Friction::Laminar,
            Some("darcy") => Friction::Darcy,
            Some("darcy2") => Friction::Darcy2,
            Some("manning") => Friction::Manning,
            Some("manning2") => Friction::Manning2,
            Some(s) => {
                return Err(SolverError::fatal(
                    CALLER,
                    format!("Uknown friction model: {s}"),
                ))
            }
        };
        // Gravity: the absolute value of the last entry of the constant.
        let gravity = match ctx.constants().get("gravity") {
            Some(Value::Real { data, .. }) if !data.is_empty() => Some(data[data.len() - 1].abs()),
            _ => None,
        };
        let use_gravity = p.logical("use gravity");
        if gravity.is_none() && use_gravity == Some(true) {
            return Err(SolverError::fatal(
                CALLER,
                "Gravity requested but not given as constant!",
            ));
        }
        if matches!(friction, Friction::Manning | Friction::Manning2) {
            if gravity.is_none() {
                return Err(SolverError::fatal(
                    CALLER,
                    "Manning equation not possible without gravity!",
                ));
            }
            if csymmetry {
                return Err(SolverError::fatal(
                    CALLER,
                    "Manning equation not applicable to axial symmetry!",
                ));
            }
        }
        let calc_heating = flag(&p, "calculate friction heating");
        let use_heating = flag(&p, "use heating source");
        let min_speed = ctx.creal(&p, "min speed")?.unwrap_or(1.0e-6);

        // The heating fields: a repeated visit of the same timestep takes the
        // previous contribution out of the energy first.
        let itime = ctx.scalar("timestep").map_or(0, |t| t.round() as i64);
        if calc_heating {
            let same = self.itime0 == Some(itime);
            let flux: Vec<f64> = ctx
                .variable("filmfriction heatflux")
                .map(|v| v.values.to_vec())
                .ok_or_else(|| {
                    SolverError::fatal(CALLER, "the FilmFriction HeatFlux variable is missing")
                })?;
            if same {
                if let Some(e) = ctx.variable_mut("heating energy") {
                    for (x, f) in e.values.iter_mut().zip(&flux) {
                        *x -= dt * f;
                    }
                }
            }
            if let Some(f) = ctx.variable_mut("filmfriction heatflux") {
                f.values.iter_mut().for_each(|x| *x = 0.0);
            }
            self.itime0 = Some(itime);
        }

        let mut maxiter = p
            .integer("nonlinear system max iterations")
            .unwrap_or(1)
            .max(1);
        if calc_heating {
            maxiter = maxiter.max(2);
        }
        let mut fsi_mult = ctx.creal(&p, "fsi velocity multiplier")?.unwrap_or(1.0);
        if itime == 1 && flag(&p, "steady state fsi start") {
            ctx.info(CALLER, "Ignoring FSI velocity for the 1st timestep", 7);
            fsi_mult = 0.0;
        }

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut local = LocalSystem::zeros(0);
        let mut dofs_idx = Vec::new();
        // FrictionHeatFlux, zeroed once per call and summed over its
        // nonlinear iterations; HeatingW, summed until DoneWeight.
        let n_heat = ctx
            .variable("filmfriction heatflux")
            .map_or(0, |v| v.values.len());
        let mut heat_flux_acc = vec![0.0; n_heat];
        let done_weight = self.heating_w.is_some();
        let mut heat_w_acc = self.heating_w.clone().unwrap_or_else(|| vec![0.0; n_heat]);
        let lagrange = elmer_fe::LagrangeBasis::new(ctx.mesh.dim());

        for _iter in 1..=maxiter {
            ctx.default_initialize()?;
            let newton = {
                let id = ctx.solver_id();
                ctx.sim
                    .slots
                    .iter()
                    .find(|s| s.id == id)
                    .is_some_and(|s| s.state.newton_active)
            };
            if let Some(s) = ctx.creal(&p, "convect condition")? {
                convect = s > 0.0;
            }
            let heat_perm = ctx
                .variable("filmfriction heatflux")
                .and_then(|v| v.perm.map(<[u32]>::to_vec));
            let first_round = !self.started;

            for &e in &elems {
                let el = ctx.mesh.element(e);
                let nodes = el.nodes.to_vec();
                let n = nodes.len();
                let d = ctx.solver_dofs().expect("checked in init");
                if d.element_nodes(e) != nodes {
                    return Err(unsupported("a p element whose basis reorders its nodes"));
                }
                let mut slots = Vec::new();
                d.element_slots(e, &mut slots);
                if slots.len() != n {
                    return Err(unsupported(
                        "p elements beyond the MINI element (edge or face dofs)",
                    ));
                }
                let nb = d.local_bubbles(e);
                let basis = d
                    .element_basis_scalar(e)
                    .map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let rule = d
                    .element_rule(e)
                    .map_err(|m| SolverError::fatal("GaussPoints", m))?;

                let bp = ctx.body_params(e).ok_or_else(|| {
                    SolverError::fatal(CALLER, format!("element {} has no body", e.0))
                })?;
                let (material, body_force) = (bp.material, bp.body_force);
                let zero = || vec![0.0; n];
                let mut load = vec![zero(); mdim + 4];
                if let Some(bf) = body_force {
                    for i in 0..mdim {
                        if let Some(v) =
                            ctx.real_at_nodes(bf, &format!("FilmFlow Bodyforce {}", i + 1), e)?
                        {
                            load[i] = v;
                        }
                    }
                    if let Some(v) = ctx.real_at_nodes(bf, "Normal Velocity", e)? {
                        load[mdim] = v;
                    }
                    if let Some(v) = ctx.real_at_nodes(bf, "Fsi Velocity", e)? {
                        load[mdim + 1] = v.iter().map(|x| fsi_mult * x).collect();
                    }
                    if let Some(v) = ctx.real_at_nodes(bf, "Flow Admittance", e)? {
                        load[mdim + 2] = v;
                    }
                    if let Some(v) = ctx.real_at_nodes(bf, "External FilmPressure", e)? {
                        load[mdim + 3] = v;
                    }
                }
                let Some(material) = material else {
                    return Err(SolverError::fatal(
                        CALLER,
                        format!("element {} has no material", e.0),
                    ));
                };
                let rho = ctx
                    .real_at_nodes(material, &density_name, e)?
                    .unwrap_or_else(zero);
                let mu = ctx
                    .real_at_nodes(material, &viscosity_name, e)?
                    .unwrap_or_else(zero);
                let mut gap = ctx
                    .real_at_nodes(material, "Gap Height", e)?
                    .unwrap_or_else(zero);
                let height = ctx.real_at_nodes(material, "Bedrock Height", e)?;
                let roughness = match friction {
                    Friction::Darcy | Friction::Darcy2 => {
                        const_real(ctx, material, "darcy roughness")?.ok_or_else(|| {
                            SolverError::fatal(CALLER, "Darcy Roughness is required")
                        })?
                    }
                    Friction::Manning | Friction::Manning2 => {
                        const_real(ctx, material, "manning coefficient")?.ok_or_else(|| {
                            SolverError::fatal(CALLER, "Manning Coefficient is required")
                        })?
                    }
                    Friction::Laminar => 0.0,
                };
                let melt_heat = if use_heating {
                    Some(
                        const_real(ctx, material, "latent heat")?
                            .ok_or_else(|| SolverError::fatal(CALLER, "Latent Heat is required"))?,
                    )
                } else {
                    None
                };
                for g in gap.iter_mut() {
                    if *g < mingap {
                        *g = mingap;
                    }
                }
                // GetVectorLocalSolution: the previous iterate at the nodes.
                ctx.element_dofs(e, &mut dofs_idx);
                let prim = ctx
                    .primary()
                    .ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
                let comp = |c: usize| -> Vec<f64> {
                    (0..n)
                        .map(|k| match dofs_idx.get(k * dofs + c) {
                            Some(&i) if i != NO_DOF => prim.values[i as usize],
                            _ => 0.0,
                        })
                        .collect()
                };
                let velo: Vec<Vec<f64>> = (0..dofs).map(comp).collect();
                let pres = velo[mdim].clone();
                let mut xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                if let Some(gf) = gap_factor {
                    // Gap Addition Factor: the nodes move by the gap in the
                    // coordinate beyond the model dimension.
                    for (x, g) in xyz.iter_mut().zip(&gap) {
                        x[mdim] += gf * g;
                    }
                }
                let opts = BulkOptions {
                    mdim,
                    csymmetry,
                    friction,
                    roughness,
                    gravity,
                    convect,
                    newton,
                    lateral_strain,
                    gradp,
                    transient,
                    min_speed,
                    calc_friction_heating: calc_heating,
                    melt_heat,
                    first_round,
                };
                let inp = BulkInput {
                    etype: el.etype,
                    xyz: &xyz,
                    ntot: n + nb,
                    rho: &rho,
                    mu: &mu,
                    gap: &gap,
                    height: height.as_deref(),
                    pres: &pres,
                    velo: &velo,
                    load: &load,
                };
                let heat = local_bulk_matrix(&*basis, &rule, &inp, &opts, &mut local)?;
                if calc_heating {
                    if let Some(perm) = &heat_perm {
                        for (k, &node) in nodes.iter().enumerate() {
                            let s = perm[node as usize];
                            if s != NO_DOF {
                                heat_flux_acc[s as usize] += heat.flux[k];
                                if !done_weight {
                                    heat_w_acc[s as usize] += heat.weight[k];
                                }
                            }
                        }
                    }
                }
                ns_condensate(&mut local, n, nb, mdim)?;
                if transient {
                    ctx.default_first_order_time(&mut local, e, dt)?;
                }
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;

            // Boundary loads: FilmPressure i and Mass Flux.
            for &e in &bnd {
                let Some(bc) = ctx.bc_params(e) else { continue };
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let mut load = vec![vec![0.0; n]; dofs];
                let mut any = false;
                for i in 0..mdim {
                    if let Some(v) = ctx.real_at_nodes(bc, &format!("FilmPressure {}", i + 1), e)? {
                        load[i] = v;
                        any = true;
                    }
                }
                if let Some(v) = ctx.real_at_nodes(bc, "Mass Flux", e)? {
                    load[mdim] = v;
                    any = true;
                }
                if !any {
                    continue;
                }
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let rule = elmer_core::Quadrature::gauss_points(
                    &elmer_fe::ElmerQuadrature,
                    el.etype,
                    None,
                );
                local_boundary_matrix(&lagrange, &rule, el.etype, &xyz, &load, mdim, &mut local)?;
                ctx.default_update(&local, e)?;
            }

            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
            // FirstRound = .FALSE. only after an iteration that did not
            // converge.
            self.started = true;
        }
        ctx.default_finish()?;
        if calc_heating {
            if let Some(f) = ctx.variable_mut("filmfriction heatflux") {
                f.values.copy_from_slice(&heat_flux_acc);
            }
            self.heating_w = Some(heat_w_acc);
        }

        if calc_heating {
            // The nodal flux is normalised by the weight, then integrated in
            // time into the energy.
            let w = self.heating_w.clone().unwrap_or_default();
            let mut tot = 0.0;
            if let Some(f) = ctx.variable("filmfriction heatflux") {
                tot = f.values.iter().sum();
            }
            let enforce = flag(&p, "enforce zero heating");
            let mut flux = Vec::new();
            if let Some(f) = ctx.variable_mut("filmfriction heatflux") {
                for (x, &wi) in f.values.iter_mut().zip(&w) {
                    if wi > 0.0 {
                        *x /= wi;
                    }
                }
                if enforce {
                    let area: f64 = w.iter().sum();
                    let cfix = -tot / area;
                    f.values.iter_mut().for_each(|x| *x += cfix);
                }
                flux = f.values.to_vec();
            }
            if let Some(e) = ctx.variable_mut("heating energy") {
                for (x, f) in e.values.iter_mut().zip(&flux) {
                    *x += dt * f;
                }
            }
            ctx.info(CALLER, &format!("Total heating power: {tot:12.5E}"), 7);
        }
        Ok(())
    }
}
