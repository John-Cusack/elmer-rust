//! `ParticleDynamics` (ParticleDynamics.F90) and the ParticleUtils.F90
//! routines its time loop calls: `SetParticlePreliminaries`,
//! `InitializeParticles` (`box random`, `box uniform cubic`),
//! `GetParticleTimeStep`, `ParticleFieldInteraction`,
//! `ParticleAdvanceTimestep`, `ParticleBoxPeriodic`, `DeleteLostParticles`,
//! `ParticlePathIntegral` and the solver's `PseudoNorm`.

use elmer_core::{DependentKind, ElemId, FieldKind, ParamMap, Params, Value};
use elmer_interp::geom::{point_in_element, PointTol};
use elmer_interp::shape::Shape;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable};

use crate::faces::{element_nodes, Faces, V3};
use crate::particles::*;

pub const PROCEDURE: (&str, &str) = ("ParticleDynamics", "ParticleDynamics");

const CALLER: &str = "ParticleDynamics";
const MAXPARFIELDS: usize = 20;

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(ParticleDynamics::default()))
}

/// `ParticleDynamics_Init`: the solver's global pseudo variable.
pub fn solver_defaults(p: &mut ParamMap) {
    if p.get("variable").is_none() {
        p.insert(
            "Variable",
            Value::Str("-nooutput -global ParticleDynamics_var".into()),
        );
    }
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("ParticleDynamics: {what}"))
}

fn flag(p: &ParamMap, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// An integer keyword that MATC may have left as a real.
fn int_of(p: &ParamMap, k: &str) -> Option<i64> {
    match p.get(k)? {
        Value::Integer(v) => v.first().map(|&x| x as i64),
        Value::Real { data, .. } => data.first().map(|&x| x.round() as i64),
        _ => None,
    }
}

/// The real array of a keyword, row-major.
fn real_array(p: &ParamMap, k: &str) -> Option<(usize, usize, Vec<f64>)> {
    match p.get(k)? {
        Value::Real { rows, cols, data } => Some((*rows, *cols, data.clone())),
        _ => None,
    }
}

/// One particle-to-field operator (`Field i`).
#[derive(Clone)]
struct ActiveField {
    oper: u8,
    var: String,
    particle_var: Option<String>,
}

/// The solver's `SAVE`d state and `GlobalParticles`.
#[derive(Default)]
pub struct ParticleDynamics {
    visited: usize,
    particles: Particles,
    faces: Option<Faces>,
    particle_to_field: bool,
    field_reset: bool,
    nstep: usize,
    located: bool,
    tottime: f64,
    // ParticleFieldInteraction
    pfi_visited: bool,
    gravity: Option<V3>,
    pot_var: Option<String>,
    pot2_var: Option<String>,
    b_var: Option<String>,
    fields: Vec<ActiveField>,
    prev_dtime: f64,
    // GetParticleTimeStep
    tprev: f64,
}

/// `ParticleElementInfo`: the basis (and its gradient) of bulk element `e1`
/// at the global point `coord`, `None` when the point is not in it.
fn particle_element_info(
    ctx: &SolverCtx<'_>,
    e1: usize,
    coord: V3,
) -> Result<Option<elmer_core::BasisEval>, SolverError> {
    let el = ctx.mesh.element(ElemId((e1 - 1) as u32));
    let x = element_nodes(ctx.mesh, e1);
    let shape = Shape::get(el.etype.0)
        .ok_or_else(|| unsupported(&format!("element type {}", el.etype.0)))?;
    let tol = PointTol {
        global: Some(-1.0),
        local: Some(1.0e3),
        numeric: None,
    };
    let probe = point_in_element(shape, &x, coord, tol, false, ctx.mesh.dim());
    if !probe.inside {
        return Ok(None);
    }
    Ok(Some(basis_at(ctx, e1, &x, probe.uvw)?))
}

fn basis_at(
    ctx: &SolverCtx<'_>,
    e1: usize,
    x: &[V3],
    uvw: V3,
) -> Result<elmer_core::BasisEval, SolverError> {
    let el = ctx.mesh.element(ElemId((e1 - 1) as u32));
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
    let mut ev = elmer_core::BasisEval::default();
    fe.basis
        .eval(el.etype, x, uvw, &mut ev)
        .map_err(|m| SolverError::fatal("ElementInfo", m))?;
    Ok(ev)
}

/// `SUM( Basis(1:n) * f(1:n) )`.
fn sum_basis(b: &[f64], f: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in b.iter().zip(f) {
        s += x * y;
    }
    s
}

/// The nodal values of a mesh variable at `nodes` (`Var % Perm` when it has
/// one), `None` where a node has no dof.
fn nodal(var: &Variable, nodes: &[u32], dofs: usize, k: usize) -> Option<Vec<f64>> {
    nodes
        .iter()
        .map(|&n| {
            let j = match &var.perm {
                Some(p) => {
                    let j = p[n as usize];
                    if j == elmer_core::NO_DOF {
                        return None;
                    }
                    j as usize
                }
                None => n as usize,
            };
            Some(var.values[dofs * j + k])
        })
        .collect()
}

/// `ListGetElementReal` of a body-force or initial-condition keyword at a
/// point of element `e1` with `basis`, the first `dummy.len()` arguments
/// given (`DummyVals`): a constant, or a MATC function of the dummies and
/// coordinates or globals, the latter at the point.
fn element_real_at(
    ctx: &SolverCtx<'_>,
    v: &Value,
    e1: usize,
    basis: &[f64],
    dummy: &[f64],
) -> Result<f64, SolverError> {
    match v {
        Value::Real { data, .. } if !data.is_empty() => Ok(data[0]),
        Value::Dependent(d) => {
            let DependentKind::Matc(expr) = &d.kind else {
                return Err(unsupported(
                    "a particle source other than a constant or MATC",
                ));
            };
            let nodes = ctx.mesh.element(ElemId((e1 - 1) as u32)).nodes;
            let x = element_nodes(ctx.mesh, e1);
            let mut t: Vec<f64> = dummy.to_vec();
            // ListParseStrToVars: the names are separated by commas.
            let names: Vec<String> = d
                .variables
                .iter()
                .flat_map(|v| v.split(','))
                .map(elmer_core::normalize_key)
                .filter(|v| !v.is_empty())
                .collect();
            for name in names.iter().skip(dummy.len()) {
                let name = name.clone();
                if let Some(c) = name.strip_prefix("coordinate ") {
                    let k: usize = c
                        .parse()
                        .map_err(|_| unsupported(&format!("argument {name}")))?;
                    let f: Vec<f64> = x.iter().map(|p| p[k - 1]).collect();
                    t.push(sum_basis(basis, &f));
                } else if let Some(var) = ctx.sim.vars.get(&name) {
                    if var.kind == FieldKind::Global {
                        t.push(var.values[0]);
                    } else {
                        let f = nodal(var, nodes, var.dofs, 0).ok_or_else(|| {
                            unsupported(&format!("argument {name} without a dof at a node"))
                        })?;
                        t.push(sum_basis(basis, &f));
                    }
                } else {
                    return Err(unsupported(&format!("argument {name}")));
                }
            }
            let matc = ctx
                .sim
                .backends
                .matc
                .as_deref()
                .ok_or_else(|| SolverError::NotWired("MATC".into()))?;
            let r = matc
                .eval(expr, &t)
                .map_err(|m| SolverError::fatal("GetMatcReal", m))?;
            Ok(r.first().copied().unwrap_or(0.0))
        }
        _ => Err(unsupported(
            "a particle source other than a constant or MATC",
        )),
    }
}

impl Solver for ParticleDynamics {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        for k in [
            "number of particle groups",
            "particle particle collision",
            "particle particle contact",
            "box particle collision",
            "box particle contact",
            "wall particle contact",
            "particle scattering",
            "runge kutta",
            "reinitialize particles",
            "particle lift",
            "velocity gradient correction",
            "source gradient correction",
            "source time correction",
            "particle accurate at face",
            "particle stop at face",
            "interface nodes",
        ] {
            if flag(&p, k) {
                return Err(unsupported(k));
            }
        }
        if p.logical("particle dt constant") == Some(false) {
            return Err(unsupported("Particle Dt Constant = False"));
        }
        if p.logical("particle locate robust") == Some(false) {
            return Err(unsupported("Particle Locate Robust = False"));
        }
        for k in [
            "output interval",
            "velocity variable name",
            "velocity variable",
            "velocity condition variable name",
            "velocity condition variable",
            "coordinate condition variable name",
            "coordinate condition variable",
            "particle to field weight",
            "particle decay distance",
            "particle decay time",
            "field decay time",
            "particle nodal weights",
            "particle release number",
            "particle release fraction",
            "initialization condition variable",
            "initialization mask variable",
            "velocity initialization method",
            "initial velocity time",
            "particle speed constant",
            "timestep size",
            "timestep distance",
            "timestep courant number",
            "timestep unisotropic courant number",
        ] {
            if p.get(k).is_some() && !(k == "output interval" && int_of(&p, k).unwrap_or(0) <= 0) {
                return Err(unsupported(k));
            }
        }
        if flag(&p, "initial coordinate search") {
            return Err(unsupported("Initial Coordinate Search"));
        }
        let model = &ctx.sim.model;
        for bc in &model.boundary_conditions {
            for k in [
                "particle outlet",
                "particle wall",
                "particle reflect",
                "particle tangent",
                "particle interact",
            ] {
                if flag(bc, k) {
                    return Err(unsupported(&format!("the boundary condition {k}")));
                }
            }
        }
        if model
            .body_forces
            .iter()
            .any(|b| b.get("particle initialization weight").is_some())
        {
            return Err(unsupported("Particle Initialization Weight"));
        }
        if ctx.mesh.dim() != 2 {
            return Err(unsupported(
                "three-dimensional meshes (2-D face march only)",
            ));
        }
        if ctx.comm.size() > 1 {
            return Err(unsupported("more than one partition"));
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "-----------------------------------------", 4);
        ctx.info(CALLER, "Following the path of the particles", 4);
        self.visited += 1;
        let p = ctx.params().clone();
        if self.visited == 1 {
            let time_order = int_of(&p, "time order").unwrap_or(2) as i32;
            self.set_preliminaries(ctx, time_order)?;
            self.particle_to_field = flag(&p, "particle to field");
            self.nstep = int_of(&p, "max timestep intervals").unwrap_or(1).max(0) as usize;
            self.located = false;
            if self.particle_to_field {
                self.field_reset = flag(&p, "particle to field reset");
            }
            self.initialize_particles(ctx, &p)?;
            self.located = false;
            if flag(&p, "particle distance") {
                self.particles.var_create("particle distance");
            }
        }
        let mut dtime = 0.0;
        for i in 1..=self.nstep {
            dtime = self.time_step(&p, i == 1, dt)?;
            let assembly = self.particle_to_field && i > 1 && !self.field_reset;
            if dtime < f64::MIN_POSITIVE {
                break;
            }
            self.tottime += dtime;
            // ParticleInitializeTime.
            self.particles.force = vec![[0.0; 3]; self.particles.force.len()];
            if !(i == 1 && self.located) {
                self.locate(ctx, &p);
                self.located = true;
            }
            // EliminateExitingParticles returns at once.
            self.field_interaction(ctx, &p, dtime, true, assembly)?;
            self.advance(&p)?;
            self.box_periodic(&p);
            self.particles.delete_lost();
            self.path_integral(ctx, &p)?;
        }
        if self.particle_to_field {
            self.locate(ctx, &p);
            self.field_interaction(ctx, &p, dtime, false, true)?;
            self.located = true;
        }
        self.pseudo_norm(ctx)?;
        ctx.info(CALLER, "All done", 4);
        ctx.info(CALLER, "-----------------------------------------", 4);
        Ok(())
    }
}

impl ParticleDynamics {
    fn locate(&mut self, ctx: &SolverCtx<'_>, p: &ParamMap) {
        let mp = MarchParams {
            stop_at_face: false,
            max_trials: int_of(p, "max particle search trials")
                .map_or(ctx.mesh.n_bulk(), |v| v.max(0) as usize),
            eps: p.real("particle hit tolerance").unwrap_or(1.0e-8),
        };
        let faces = self.faces.as_ref().expect("set up on the first visit");
        locate_particles(
            ctx.mesh,
            faces,
            &mp,
            &mut self.particles,
            flag(p, "particle accurate always"),
        );
    }

    /// `SetParticlePreliminaries`: the mesh bounding box and its edges.
    fn set_preliminaries(
        &mut self,
        ctx: &SolverCtx<'_>,
        time_order: i32,
    ) -> Result<(), SolverError> {
        let pt = &mut self.particles;
        pt.time_order = time_order;
        pt.dim = ctx.mesh.dim();
        let n = ctx.mesh.n_nodes();
        let mut mn = [f64::MAX; 3];
        let mut mx = [-f64::MAX; 3];
        for i in 0..n as u32 {
            let c = ctx.mesh.coord(i);
            for d in 0..3 {
                mn[d] = mn[d].min(c[d]);
                mx[d] = mx[d].max(c[d]);
            }
        }
        pt.local_min = mn;
        pt.local_max = mx;
        pt.global_min = mn;
        pt.global_max = mx;
        self.faces = Some(Faces::build(ctx.mesh).map_err(|m| unsupported(&m))?);
        Ok(())
    }

    /// `InitializeParticles` without groups, masks or weights.
    fn initialize_particles(
        &mut self,
        ctx: &SolverCtx<'_>,
        p: &ParamMap,
    ) -> Result<(), SolverError> {
        const C: &str = "InitializeParticles";
        let dim = self.particles.dim;
        let method = p
            .string("coordinate initialization method")
            .map(elmer_core::normalize_key)
            .unwrap_or_default();
        if !matches!(method.as_str(), "box random" | "box uniform cubic") {
            return Err(unsupported(&format!(
                "Coordinate Initialization Method = {method}"
            )));
        }
        let pt = &self.particles;
        let eps = ctx.creal(p, "wall particle radius")?.unwrap_or(1.0e-8);
        let mut minc = [0.0f64; 3];
        let mut maxc = [0.0f64; 3];
        for k in 0..3 {
            minc[k] = match ctx.creal(p, &format!("min initial coordinate {}", k + 1))? {
                Some(v) => v,
                None if k < 2 => pt.local_min[k] + eps,
                None => pt.local_min[k],
            };
            maxc[k] = match ctx.creal(p, &format!("max initial coordinate {}", k + 1))? {
                Some(v) => v,
                None => pt.local_max[k] - eps,
            };
        }
        let mut diam = 0.0;
        let (mut nx, mut ny, mut nz) = (0i64, 0i64, 1i64);
        if method == "box uniform cubic" {
            diam = match ctx.creal(p, "particle cell radius")? {
                Some(r) => 2.0 * r,
                None => match ctx.creal(p, "particle cell size")? {
                    Some(s) => s,
                    None => match ctx.creal(p, "particle radius")? {
                        Some(r) => 2.0 * r,
                        None => return Err(SolverError::fatal(C, "Size of unit cell not given")),
                    },
                },
            };
            nx = ((maxc[0] - minc[0]) / diam).round() as i64;
            ny = ((maxc[1] - minc[1]) / diam).round() as i64;
            if dim == 3 {
                nz = ((maxc[2] - minc[2]) / diam).round() as i64;
            }
        }
        let new = if method == "box uniform cubic" {
            (nx * ny * nz) as usize
        } else {
            int_of(p, "number of particles")
                .ok_or_else(|| {
                    unsupported("the number of particles other than Number of Particles")
                })?
                .max(0) as usize
        };
        if new == 0 {
            return Err(SolverError::fatal(C, "No Particles to Initialize"));
        }
        ctx.info(C, &format!("Number of Particles: {new}"), 6);
        let pt = &mut self.particles;
        pt.allocate(new);
        pt.number = new;
        if method == "box random" {
            for k in 0..new {
                let mut c = [0.0f64; 3];
                for j in 0..dim {
                    c[j] = minc[j] + (maxc[j] - minc[j]) * elmer_itermethods::even_random();
                }
                pt.coord[k][..dim].copy_from_slice(&c[..dim]);
            }
        } else {
            ctx.info(C, "Initializing particles in a grid", 10);
            let mut k = 0;
            for ix in 1..=nx {
                for iy in 1..=ny {
                    for iz in 1..=nz {
                        pt.coord[k][0] = minc[0] + (1.0 * ix as f64 - 0.5) * diam;
                        pt.coord[k][1] = minc[1] + (1.0 * iy as f64 - 0.5) * diam;
                        if dim == 3 {
                            pt.coord[k][2] = minc[2] + (1.0 * iz as f64 - 0.5) * diam;
                        }
                        k += 1;
                    }
                }
            }
        }
        if let Some((rows, cols, data)) = real_array(p, "initial velocity") {
            if cols != dim {
                return Err(SolverError::fatal(
                    "ParticleTracker",
                    "Wrong dimension in Initial Velocity",
                ));
            }
            if rows == 1 {
                for v in pt.velo.iter_mut() {
                    v[..dim].copy_from_slice(&data[..dim]);
                }
            } else if rows != new {
                return Err(SolverError::fatal(
                    "ParticleTracker",
                    "Wrong number of particles in Initial Velocity",
                ));
            } else {
                for (i, v) in pt.velo.iter_mut().enumerate() {
                    v[..dim].copy_from_slice(&data[i * dim..(i + 1) * dim]);
                }
            }
        }
        for s in pt.status.iter_mut() {
            *s = PARTICLE_INITIATED;
        }
        pt.prev_coord = pt.coord.clone();
        let model = &ctx.sim.model;
        for name in ["particle distance integral", "particle time integral"] {
            let src = format!("{name} source");
            if model.body_forces.iter().any(|b| b.get(&src).is_some()) {
                self.particles.var_create(name);
                self.variable_initialize(ctx, p, name)?;
            }
        }
        Ok(())
    }

    /// `ParticleVariableInitialize` from the `Initial Condition` keyword of
    /// the same name (there is no mesh variable of it).
    fn variable_initialize(
        &mut self,
        ctx: &SolverCtx<'_>,
        p: &ParamMap,
        name: &str,
    ) -> Result<(), SolverError> {
        if self.particles.element[..self.particles.number].contains(&0) {
            // LocateParticles without a wall kernel.
            self.locate(ctx, p);
        }
        if ctx.sim.vars.get(name).is_some() {
            return Err(unsupported(&format!(
                "{name} initialised from a mesh variable"
            )));
        }
        if !ctx
            .sim
            .model
            .initial_conditions
            .iter()
            .any(|ic| ic.get(name).is_some())
        {
            return Ok(());
        }
        let iv = self.particles.var(name).expect("created");
        for i in 0..self.particles.number {
            let e1 = self.particles.element[i];
            if e1 == 0 {
                continue;
            }
            let Some(ev) = particle_element_info(ctx, e1, self.particles.coord[i])? else {
                continue;
            };
            let bp = ctx.body_params(ElemId((e1 - 1) as u32));
            let Some(v) = bp
                .and_then(|b| b.initial_condition)
                .and_then(|ic| ic.get(name))
            else {
                // ListGetElementReal of an absent keyword: its default 0.
                self.particles.vars[iv].values[i] = 0.0;
                continue;
            };
            let n = ctx.mesh.element(ElemId((e1 - 1) as u32)).nodes.len();
            let val = element_real_at(ctx, v, e1, &ev.basis[..n], &[])?;
            self.particles.vars[iv].values[i] = val;
        }
        Ok(())
    }

    /// `GetParticleTimeStep` with a constant step: `Simulation Timestep
    /// Sizes` (the simulation step over `Max Timestep Intervals`) or `Max
    /// Cumulative Time`.
    fn time_step(&mut self, p: &ParamMap, init: bool, sim_dt: f64) -> Result<f64, SolverError> {
        let nstep = int_of(p, "max timestep intervals").unwrap_or(1);
        let tfin = if flag(p, "simulation timestep sizes") {
            Some(sim_dt)
        } else {
            p.real("max cumulative time")
        };
        if init {
            self.tprev = 0.0;
        }
        let Some(tfin) = tfin else {
            return Err(SolverError::fatal(
                "GetParticleTimeStep",
                "Cannot determine timestep size!",
            ));
        };
        let mut dt = tfin / nstep as f64;
        if dt + self.tprev > tfin {
            dt = tfin - self.tprev;
        }
        self.tprev += dt;
        self.particles.dtime = dt;
        Ok(dt)
    }

    /// `ParticleFieldInteraction`: the forces on the particles from gravity,
    /// damping and the electric and magnetic fields, and the particle data
    /// summed to nodal fields.
    fn field_interaction(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        p: &ParamMap,
        dtime: f64,
        set_particles: bool,
        set_fields: bool,
    ) -> Result<(), SolverError> {
        const C: &str = "ParticleFieldInteraction";
        let dim = ctx.mesh.dim();
        if !self.pfi_visited {
            if flag(p, "particle gravity") {
                let g = real_array(ctx.constants(), "gravity")
                    .ok_or_else(|| SolverError::fatal(C, "Gravity and Lift requires gravity!"))?
                    .2;
                self.gravity = Some([g[3] * g[0], g[3] * g[1], g[3] * g[2]]);
            }
            let var_name = |a: &str, b: &str| {
                p.string(a)
                    .or_else(|| p.string(b))
                    .map(elmer_core::normalize_key)
            };
            for (slot, a, b, what) in [
                (
                    0,
                    "potential variable name",
                    "potential variable",
                    "Potential",
                ),
                (
                    1,
                    "secondary potential variable name",
                    "secondary potential variable",
                    "Potential",
                ),
                (
                    2,
                    "magnetic field variable name",
                    "magnetic field variable",
                    "Magnetic",
                ),
            ] {
                if let Some(n) = var_name(a, b) {
                    if ctx.sim.vars.get(&n).is_none() {
                        return Err(SolverError::fatal(
                            C,
                            format!("{what} field variable does not exist: {n}"),
                        ));
                    }
                    match slot {
                        0 => self.pot_var = Some(n),
                        1 => self.pot2_var = Some(n),
                        _ => self.b_var = Some(n),
                    }
                }
            }
            if self.particle_to_field {
                if p.get("particle to field mode").is_some() {
                    return Err(SolverError::fatal(
                        C,
                        "> Particle to Field Mode < is obsolete. Use > Field 1 < instead!",
                    ));
                }
                for i in 1..=MAXPARFIELDS {
                    let Some(mode) = p.string(&format!("field {i}")) else {
                        break;
                    };
                    let mode = elmer_core::normalize_key(mode);
                    let (oper, particle_var) = match mode.as_str() {
                        "weight" => (1, None),
                        "energy" => (2, None),
                        "kinetic energy" => (3, None),
                        "potential energy" => (4, None),
                        "electrostatic energy" => (5, None),
                        "charge" => (6, None),
                        "speed" => (7, None),
                        "force" => (8, None),
                        m if m.starts_with("particle") => match self.particles.var(m) {
                            Some(_) => (9, Some(m.to_string())),
                            None => {
                                return Err(unsupported(&format!(
                                    "field mode {m} without its particle variable"
                                )))
                            }
                        },
                        m => return Err(SolverError::fatal(C, format!("Unknown field mode: {m}"))),
                    };
                    if p.logical(&format!("field normalize {i}")) == Some(true) {
                        return Err(unsupported("Field Normalize with Particle Nodal Weights"));
                    }
                    let var = particle_var
                        .clone()
                        .unwrap_or_else(|| format!("particle {mode}"));
                    if ctx.sim.vars.get(&var).is_none() {
                        ctx.info(C, &format!("Creating variable: {var}"), 4);
                        let n = ctx.mesh.n_nodes();
                        ctx.sim.vars.add(Variable::new(
                            &var,
                            FieldKind::Nodal,
                            1,
                            None,
                            vec![0.0; n],
                        ))?;
                    }
                    self.fields.push(ActiveField {
                        oper,
                        var,
                        particle_var,
                    });
                }
                ctx.info(
                    C,
                    &format!("Number of particle fields: {}", self.fields.len()),
                    6,
                );
            }
            self.pfi_visited = true;
        }
        if self.particle_to_field && self.field_reset {
            for f in &self.fields {
                if let Some(v) = ctx.sim.vars.get_mut(&f.var) {
                    v.values.iter_mut().for_each(|x| *x = 0.0);
                }
            }
        }
        let user_coeff = ctx
            .creal(p, "particle to field coefficient")?
            .unwrap_or(1.0);
        let mass = ctx.creal(p, "particle mass")?;
        let damping = ctx.creal(p, "particle damping")?;
        let charge = ctx.creal(p, "particle charge")?;
        if self.gravity.is_some() && mass.is_none() {
            ctx.warn(C, "> Particle Mass < is needed by gravity!");
        }
        let got_pot = self.pot_var.is_some();
        if (got_pot || self.b_var.is_some()) && charge.is_none() {
            return Err(SolverError::fatal(
                C,
                "> Particle Charge < is needed by electric and magnetic fields!",
            ));
        }
        let (mass, charge) = (mass.unwrap_or(0.0), charge.unwrap_or(0.0));
        let time_order = self.particles.time_order;
        let got_field = self.particle_to_field || got_pot || self.b_var.is_some();
        let mut adds: Vec<(usize, u32, f64)> = Vec::new();
        for no in 0..self.particles.number {
            let status = self.particles.status[no];
            if status >= PARTICLE_LOST || status <= PARTICLE_INITIATED {
                continue;
            }
            let e1 = self.particles.element[no];
            if e1 < 1 {
                continue;
            }
            let coord = self.particles.coord[no];
            let velo = self.particles.velo[no];
            let mut force = [0.0f64; 3];
            if let Some(g) = self.gravity {
                for d in 0..3 {
                    force[d] += g[d] * mass;
                }
            }
            if let (Some(dc), 2) = (damping, time_order) {
                for d in 0..3 {
                    force[d] -= dc * velo[d] * mass;
                }
            }
            let mut pot_at = 0.0;
            let mut basis: Vec<f64> = Vec::new();
            if got_field {
                let Some(ev) = particle_element_info(ctx, e1, coord)? else {
                    ctx.warn(C, "Particle not in element");
                    continue;
                };
                let nodes = ctx.mesh.element(ElemId((e1 - 1) as u32)).nodes;
                let n = nodes.len();
                basis = ev.basis[..n].to_vec();
                // GetScalarFieldInMesh with the gradient: `Force - charge * grad`.
                let scalar = |name: &str, force: &mut [f64; 3]| -> f64 {
                    let var = ctx.sim.vars.get(name).expect("checked");
                    match nodal(var, nodes, var.dofs, 0) {
                        Some(f) => {
                            for d in 0..dim {
                                let mut g = 0.0;
                                for (a, b) in ev.dbasis_dx[..n].iter().zip(&f) {
                                    g += a[d] * b;
                                }
                                force[d] -= charge * g;
                            }
                            sum_basis(&basis, &f)
                        }
                        None => 0.0,
                    }
                };
                if let Some(pv) = &self.pot_var {
                    pot_at = scalar(pv, &mut force);
                }
                if let Some(bn) = &self.b_var {
                    // GetVectorFieldInMesh.
                    let var = ctx.sim.vars.get(bn).expect("checked");
                    let dofs = var.dofs.min(3);
                    let mut b = [0.0f64; 3];
                    for (k, bk) in b.iter_mut().enumerate().take(dofs) {
                        if let Some(f) = nodal(var, nodes, var.dofs, k) {
                            *bk = sum_basis(&basis, &f);
                        }
                    }
                    // CrossProduct( Velo, B ).
                    let c = [
                        velo[1] * b[2] - velo[2] * b[1],
                        -velo[0] * b[2] + velo[2] * b[0],
                        velo[0] * b[1] - velo[1] * b[0],
                    ];
                    for d in 0..3 {
                        force[d] += charge * c[d];
                    }
                    if dim < 3 {
                        force[2] = 0.0;
                    }
                }
                if let Some(pv) = &self.pot2_var {
                    pot_at = scalar(pv, &mut force);
                }
            }
            if set_fields {
                let nodes = ctx.mesh.element(ElemId((e1 - 1) as u32)).nodes;
                let val_coeff = user_coeff * 1.0;
                for f in &self.fields {
                    let mut v2 = 0.0;
                    for d in 0..dim {
                        v2 += velo[d].powi(2);
                    }
                    let mut gc = 0.0;
                    if let Some(g) = self.gravity {
                        for d in 0..dim {
                            gc += g[d] * coord[d];
                        }
                    }
                    let val = match f.oper {
                        1 => 1.0,
                        2 => 0.5 * mass * v2 + mass * gc + charge * pot_at,
                        3 => 0.5 * mass * v2,
                        4 => mass * gc,
                        5 => charge * pot_at,
                        6 => charge,
                        // SQRT( SUM( VeloAtPoint(1:dim)**2 ) ): no velocity field.
                        7 => 0.0,
                        8 => {
                            let mut s = 0.0;
                            for d in 0..dim {
                                s += force[d].powi(2);
                            }
                            s.sqrt()
                        }
                        _ => {
                            let pvn = f.particle_var.as_deref().expect("set");
                            self.particles.vars[self.particles.var(pvn).expect("set")].values[no]
                        }
                    };
                    let val = val * val_coeff;
                    let vid = ctx.sim.vars.id(&f.var).expect("created");
                    for (i, &nd) in nodes.iter().enumerate() {
                        adds.push((vid.0, nd, basis[i] * val));
                    }
                }
            }
            if set_particles {
                for d in 0..dim {
                    self.particles.force[no][d] += force[d];
                }
                self.particles.status[no] = status;
            }
        }
        for (vid, nd, x) in adds {
            let v = ctx.sim.vars.by_id_mut(elmer_solver::VarId(vid));
            let j = match &v.perm {
                Some(p) => {
                    let j = p[nd as usize];
                    if j == elmer_core::NO_DOF {
                        continue;
                    }
                    j as usize
                }
                None => nd as usize,
            };
            v.values[j] += x;
        }
        self.prev_dtime = dtime;
        Ok(())
    }

    /// `ParticleAdvanceTimestep` with a constant step, without Runge-Kutta.
    fn advance(&mut self, p: &ParamMap) -> Result<(), SolverError> {
        const C: &str = "ParticleAdvanceTimestep";
        let pt = &mut self.particles;
        let dim = pt.dim;
        let dtime = pt.dtime;
        let time_var = pt.var("particle time");
        let dist_var = pt.var("particle distance");
        let (mass, drag) = match pt.time_order {
            2 => (
                p.real("particle mass")
                    .ok_or_else(|| SolverError::fatal(C, "> Particle Mass < should be given!"))?,
                0.0,
            ),
            1 => (
                0.0,
                p.real("particle drag coefficient").ok_or_else(|| {
                    SolverError::fatal(C, "> Particle Drag Coefficient < should be given!")
                })?,
            ),
            _ => (0.0, 0.0),
        };
        for no in 0..pt.number {
            let status = pt.status[no];
            if status >= PARTICLE_LOST
                || status <= PARTICLE_INITIATED
                || status == PARTICLE_WALLBOUNDARY
            {
                continue;
            }
            if let Some(tv) = time_var {
                pt.vars[tv].values[no] += pt.dtime;
            }
            if status == PARTICLE_FIXEDCOORD {
                pt.velo[no] = [0.0; 3];
                continue;
            } else if status == PARTICLE_FIXEDVELO {
            } else if pt.time_order == 2 {
                for d in 0..dim {
                    pt.velo[no][d] += dtime * pt.force[no][d] / mass;
                }
            } else if pt.time_order == 1 {
                for d in 0..dim {
                    pt.velo[no][d] = pt.force[no][d] / drag;
                }
            } else if pt.time_order != 0 {
                return Err(SolverError::fatal(C, "Unknown time order"));
            }
            let mut dc = [0.0f64; 3];
            for d in 0..dim {
                dc[d] = dtime * pt.velo[no][d];
            }
            pt.prev_coord[no] = pt.coord[no];
            pt.force[no] = [0.0; 3];
            pt.status[no] = PARTICLE_READY;
            for d in 0..dim {
                pt.coord[no][d] += dc[d];
            }
            if let Some(dv) = dist_var {
                let mut s = 0.0;
                for x in &dc[..dim] {
                    s += x.powi(2);
                }
                pt.vars[dv].values[no] += s.sqrt();
            }
        }
        pt.time += pt.dtime;
        Ok(())
    }

    /// `ParticleBoxPeriodic` (`Box Particle Periodic` or `Box Periodic
    /// Directions`).
    fn box_periodic(&mut self, p: &ParamMap) {
        let pt = &mut self.particles;
        let dirs: Vec<usize> = match p.get("box periodic directions") {
            Some(Value::Integer(v)) => v.iter().map(|&d| d as usize).collect(),
            _ if flag(p, "box particle periodic") => (1..=pt.dim).collect(),
            _ => return,
        };
        if dirs.is_empty() {
            return;
        }
        let mut mn = pt.global_min;
        let mut mx = pt.global_max;
        let mut ext = mx[0] - mn[0];
        for d in 1..3 {
            ext = ext.max(mx[d] - mn[d]);
        }
        let eps_coord = f64::EPSILON * ext;
        for d in 0..3 {
            mn[d] += eps_coord;
            mx[d] -= eps_coord;
        }
        for no in 0..pt.number {
            let status = pt.status[no];
            if status >= PARTICLE_LOST || status < PARTICLE_INITIATED {
                continue;
            }
            for _ in 0..dirs.len() {
                let mut mapped = false;
                for &k1 in &dirs {
                    let k = k1 - 1;
                    let c = pt.coord[no][k];
                    if c < mn[k] {
                        pt.coord[no][k] = mx[k] - mn[k] + c;
                        mapped = true;
                    } else if c > mx[k] {
                        pt.coord[no][k] = mn[k] - mx[k] + c;
                        mapped = true;
                    }
                }
                if !mapped {
                    break;
                }
            }
        }
    }

    /// `ParticlePathIntegral` without the gradient and time corrections:
    /// the time and distance integrals of their body-force sources, at the
    /// particle's new position in the element it was last located in.
    fn path_integral(&mut self, ctx: &SolverCtx<'_>, p: &ParamMap) -> Result<(), SolverError> {
        const C: &str = "ParticlePathIntegral";
        let tv = self.particles.var("particle time integral");
        let dv = self.particles.var("particle distance integral");
        if tv.is_none() && dv.is_none() {
            return Ok(());
        }
        let use_dummy = flag(p, "particle integral dummy argument");
        let dim = self.particles.dim;
        let dtime = self.particles.dtime;
        for no in 0..self.particles.number {
            let status = self.particles.status[no];
            if status >= PARTICLE_LOST
                || status <= PARTICLE_INITIATED
                || status == PARTICLE_WALLBOUNDARY
                || status == PARTICLE_FIXEDCOORD
            {
                continue;
            }
            let e1 = self.particles.element[no];
            if e1 == 0 {
                continue;
            }
            let el = ctx.mesh.element(ElemId((e1 - 1) as u32));
            let Some(bf) = ctx
                .body_params(ElemId((e1 - 1) as u32))
                .and_then(|b| b.body_force)
            else {
                continue;
            };
            let x = element_nodes(ctx.mesh, e1);
            let shape = Shape::get(el.etype.0).ok_or_else(|| unsupported("element type"))?;
            let coord = self.particles.coord[no];
            let probe =
                point_in_element(shape, &x, coord, PointTol::default(), false, ctx.mesh.dim());
            if !probe.inside {
                continue;
            }
            let ev = basis_at(ctx, e1, &x, probe.uvw)?;
            let n = el.nodes.len();
            let basis = &ev.basis[..n];
            let prev = self.particles.prev_coord[no];
            if let Some(tv) = tv {
                if let Some(v) = bf.get("particle time integral source") {
                    let dummy = [self.particles.vars[tv].values[no]];
                    let source =
                        element_real_at(ctx, v, e1, basis, if use_dummy { &dummy } else { &[] })?;
                    self.particles.vars[tv].values[no] += dtime * source;
                }
            }
            if let Some(dv) = dv {
                let mut ds = 0.0;
                for d in 0..dim {
                    ds += (prev[d] - coord[d]).powi(2);
                }
                let ds = ds.sqrt();
                if let Some(v) = bf.get("particle distance integral source") {
                    let dummy = [self.particles.vars[dv].values[no]];
                    let source =
                        element_real_at(ctx, v, e1, basis, if use_dummy { &dummy } else { &[] })?;
                    self.particles.vars[dv].values[no] += ds * source;
                }
            }
        }
        let _ = C;
        Ok(())
    }

    /// `PseudoNorm`: the root mean square particle speed, as the norm and
    /// the value of the solver's global variable.
    fn pseudo_norm(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let pt = &self.particles;
        let dim = pt.dim;
        let mut v2 = 0.0;
        for i in 0..pt.number {
            let mut s = 0.0;
            for d in 0..dim {
                s += pt.velo[i][d].powi(2);
            }
            v2 += s;
        }
        let nrm = (v2 / pt.number as f64).sqrt();
        let name = ctx.variable_name().map(str::to_string);
        if let Some(var) = name.and_then(|n| ctx.sim.vars.get_mut(&n)) {
            var.conv.norm = nrm;
            var.values.iter_mut().for_each(|x| *x = nrm);
        }
        Ok(())
    }
}
