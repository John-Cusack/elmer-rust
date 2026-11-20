//! Mortar contact (`Apply Contact BCs`): `DetermineContact`
//! (ContactUtils.F90) with the projectors it works on and the constraint
//! matrix they give, as the [`elmer_solver::ContactBackend`] the library
//! calls during a solver's nonlinear iteration.
//!
//! Per iteration (DefUtils.F90):
//! 1. `DefaultInitialize`: the `Mortar BC Nonlinear` projectors are generated
//!    on the mesh displaced by the current solution (`Displace Mesh At Init`).
//! 2. `DefaultFinishBoundaryAssembly`: `DetermineContact` on the assembled
//!    system, before the Dirichlet conditions: the contact load `A x - b`, the
//!    contact pressure, the mortar distance, the normal and tangent active
//!    sets, and their transfer to the master side.
//! 3. `DefaultSolve`: `GenerateConstraintMatrix` from the projectors with
//!    `MortarBC % Active`, `Rhs` and `Perm`.
//!
//! Ported: Cartesian slip and tie contact without residual mode on a flat
//! projector of nodal dofs. Refused: normal-tangential contact, stick and
//! friction, residual mode, velocities, the plane, rotational and normal
//! projectors, p-elements, conservative and conditional active sets, linear
//! gaps and loads, Lagrange-multiplier based contact, dual projectors.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use elmer_core::{BasisEval, ElemId, ElementRef, FieldKind, LinearOperator, MeshView, NullLogger, Params, NO_DOF};
use elmer_fe::{element_info, element_type, gauss_points, MetricStatus};
use elmer_solver::{SolverCtx, SolverError, Variable};

use crate::constraint::{generate_constraint_matrix, ConstraintTarget};
use crate::interface::{create_interface_meshes, map_interface_coordinate};
use crate::projector::{generate_mortar_projectors, Projector, ProjectorEnv};

const CALLER: &str = "DetermineContact";

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

fn fatal(msg: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, msg)
}

fn mortar_err(e: crate::MortarError) -> SolverError {
    SolverError::from(e)
}

/// Replaces the mesh's node coordinates (the runner's movable mesh).
pub type NodeMover = Arc<dyn Fn(&[[f64; 3]]) + Send + Sync>;

thread_local! {
    static MOVER: std::cell::RefCell<Option<NodeMover>> = const { std::cell::RefCell::new(None) };
}

/// Installs the runner's node mover for the solvers run on this thread,
/// returning the previous one. Without one the projectors are generated on
/// a displaced view of the mesh, which leaves the coordinates untouched.
pub fn install_node_mover(m: Option<NodeMover>) -> Option<NodeMover> {
    MOVER.with(|c| std::mem::replace(&mut *c.borrow_mut(), m))
}

fn mover() -> Option<NodeMover> {
    MOVER.with(|c| c.borrow().clone())
}

/// `DisplaceMesh( Mesh, Update, sgn, Perm, DOFs )` without `UpdateDirs`
/// on a copy of the coordinates.
fn displaced(coords: &[[f64; 3]], update: &[f64], sgn: f64, perm: &[u32], dofs: usize) -> Vec<[f64; 3]> {
    let mut out = coords.to_vec();
    for (i, x) in out.iter_mut().enumerate().take(perm.len()) {
        let k = perm[i];
        if k == NO_DOF {
            continue;
        }
        let k = dofs * k as usize;
        for d in 0..dofs.min(3) {
            x[d] += sgn * update[k + d];
        }
    }
    out
}

/// The mesh with the solution added to its coordinates (`DisplaceMesh(+1,
/// Values, Perm, Dofs)` without `UpdateDirs`: the first `min(dofs, 3)`
/// components).
struct DisplacedMesh<'a> {
    mesh: &'a dyn MeshView,
    values: &'a [f64],
    perm: &'a [u32],
    dofs: usize,
}

impl MeshView for DisplacedMesh<'_> {
    fn dim(&self) -> usize {
        self.mesh.dim()
    }
    fn n_nodes(&self) -> usize {
        self.mesh.n_nodes()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        let mut x = self.mesh.coord(node);
        if let Some(&s) = self.perm.get(node as usize) {
            if s != NO_DOF {
                for (d, xd) in x.iter_mut().enumerate().take(self.dofs.min(3)) {
                    *xd += self.values[self.dofs * s as usize + d];
                }
            }
        }
        x
    }
    fn n_bulk(&self) -> usize {
        self.mesh.n_bulk()
    }
    fn n_boundary(&self) -> usize {
        self.mesh.n_boundary()
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        self.mesh.element(e)
    }
}

/// `Model % Solver % MortarBCs(bc)` of a contact boundary.
#[derive(Default)]
struct MortarBc {
    /// `Active`, `cdofs` per projector row.
    active: Vec<bool>,
    /// `Rhs`, `cdofs` per projector row.
    rhs: Vec<f64>,
    /// `Perm`: parent node -> projector row + 1 (0 when absent).
    perm: Vec<u32>,
    master_scale: f64,
}

/// One contact solver's state.
#[derive(Default)]
struct SolverContact {
    /// `MortarBCs(i) % Projector`, regenerated at each `DefaultInitialize`.
    projectors: Vec<Option<Projector>>,
    /// `Flat Projector Coordinate` as `FlatInterfaceMeshes` stores it.
    flat_dim: Vec<usize>,
    bcs: HashMap<usize, MortarBc>,
}

/// The contact backend: the state of every contact solver, by solver id.
pub struct MortarContact {
    /// Where `Save Projector` files go.
    pub case_dir: PathBuf,
    /// `CoordinateSystemDimension()`.
    pub coord_dim: usize,
    state: Mutex<HashMap<usize, SolverContact>>,
    /// `DetermineContact`'s `FirstTime` (a `SAVE`d local): true until its
    /// first call returns.
    first_time: Mutex<bool>,
}

impl MortarContact {
    pub fn new(case_dir: PathBuf, coord_dim: usize) -> Arc<Self> {
        Arc::new(Self { case_dir, coord_dim, state: Mutex::new(HashMap::new()), first_time: Mutex::new(true) })
    }
}

/// Solver-level keywords whose code is not ported.
fn check_solver(p: &dyn Params, simulation: &dyn Params) -> Result<(), SolverError> {
    for k in [
        "apply mortar bcs",
        "apply integral bcs",
        "apply contact velocity",
        "use nodal normal",
        "debug normals",
        "use lagrange multiplier for contact",
        "linear system residual mode",
        "contact linear basis",
        "eliminate linear constraints",
    ] {
        if flag(p, k) {
            return Err(unsupported(&format!("'{k}' with Apply Contact BCs")));
        }
    }
    for k in [
        "apply limiter conservative add after iterations",
        "apply limiter conservative remove after iterations",
        "nonlinear system initial friction iterations",
        "mortar bc master solver",
    ] {
        if p.get(k).is_some() {
            return Err(unsupported(&format!("'{k}' with Apply Contact BCs")));
        }
    }
    for k in ["contact bcs linear gap", "contact bcs linear loads"] {
        if flag(simulation, k) {
            return Err(unsupported(k));
        }
    }
    let st = simulation.string("simulation type").unwrap_or("steady state").to_ascii_lowercase();
    if st == "transient" {
        return Err(unsupported("contact in a transient simulation (contact velocity)"));
    }
    Ok(())
}

/// Boundary-level keywords of a contact boundary whose code is not ported.
fn check_bc(bc: &dyn Params, var: &str) -> Result<(), SolverError> {
    for k in [
        "rotational projector",
        "cylindrical projector",
        "normal projector",
        "plane projector",
        "create dual projector",
        "stick contact global",
    ] {
        if flag(bc, k) {
            return Err(unsupported(&format!("'{k}' on a contact boundary")));
        }
    }
    if flag(bc, &format!("normal-tangential {var}")) {
        return Err(unsupported("normal-tangential contact"));
    }
    for k in [
        "contact bc",
        "stick contact coefficient",
        "contact active condition",
        "contact passive condition",
        "stick active condition",
        "stick passive condition",
        "contact active set minimum",
        "dynamic friction coefficient",
        "normal sign negative",
        "normal sign positive",
        "nonlinear system initial tie contact",
        "nonlinear system initial stick contact",
        "nonlinear system initial slip contact",
        "transient system initial tie contact",
        "transient system initial stick contact",
        "transient system initial slip contact",
        "initial tie contact",
        "initial stick contact",
        "initial slip contact",
    ] {
        if bc.get(k).is_some() {
            return Err(unsupported(&format!("'{k}'")));
        }
    }
    if !flag(bc, "mortar bc nonlinear") {
        return Err(unsupported("a contact boundary without Mortar BC Nonlinear"));
    }
    Ok(())
}

/// The contact type of a boundary (CU:344-373).
#[derive(Clone, Copy, PartialEq, Eq)]
enum ContactType {
    Slip,
    Tie,
}

fn contact_type(bc: &dyn Params, ctx: &SolverCtx<'_>) -> Result<ContactType, SolverError> {
    if let Some(t) = bc.string("contact type") {
        return match t.to_ascii_lowercase().as_str() {
            "slip" | "slide" => Ok(ContactType::Slip),
            "tie" => Ok(ContactType::Tie),
            "stick" | "friction" => Err(unsupported(&format!("Contact Type = {t}"))),
            other => Err(fatal(format!("Unknown contact type: {other}"))),
        };
    }
    for (k, t) in [("stick contact", None), ("tie contact", Some(ContactType::Tie)), ("friction contact", None)] {
        if let Some(v) = bc.logical(k) {
            return match (v, t) {
                (true, Some(t)) => Ok(t),
                (true, None) => Err(unsupported(k)),
                (false, _) => Ok(ContactType::Slip),
            };
        }
    }
    for k in ["slip contact", "slide contact"] {
        if bc.logical(k).is_some() {
            return Ok(ContactType::Slip);
        }
    }
    ctx.warn(CALLER, "No contact type given, assuming > Slip Contact <");
    Ok(ContactType::Slip)
}

impl elmer_solver::ContactBackend for MortarContact {
    fn initialize(&self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        check_solver(ctx.params(), ctx.simulation_params())?;
        if ctx.solver_dofs().is_some() {
            return Err(unsupported("contact on p-elements"));
        }
        let name = ctx.variable_name().ok_or_else(|| fatal("solver has no variable"))?.to_string();
        // `<var> Contact Load` (MainUtils.F90:2437-2466): the solver
        // variable's dofs and perm.
        let load = format!("{name} contact load");
        if ctx.sim.vars.id(&load).is_none() {
            let v = ctx.sim.vars.get(&name).ok_or_else(|| fatal("no solver variable"))?;
            let var = Variable::new(&load, FieldKind::Nodal, v.dofs, v.perm.clone(), vec![0.0; v.values.len()]);
            ctx.sim.vars.add(var)?;
        }
        for (i, bc) in ctx.sim.model.boundary_conditions.iter().enumerate() {
            if bc.get("mortar bc").is_some() || bc.get("contact bc").is_some() {
                check_bc(bc, &name).map_err(|e| match e {
                    SolverError::Unsupported(m) => SolverError::Unsupported(format!("{m} (boundary condition {})", i + 1)),
                    e => e,
                })?;
            }
        }
        // DefaultInitialize (DefUtils.F90:3916-3952): `Displace Mesh At Init`
        // (ElasticSolver_Init's default) moves the mesh by the solution, the
        // nonlinear projectors are generated, and the mesh is moved back:
        // the round trip is what the rest of the iteration sees.
        let f = ctx.primary().ok_or_else(|| fatal("solver has no variable"))?;
        let perm = f.perm.ok_or_else(|| fatal("solver variable has no permutation"))?.to_vec();
        let (values, dofs) = (f.values.to_vec(), f.dofs);
        let any_proj = ctx.sim.model.boundary_conditions.iter().any(|b| flag(b, "mortar bc nonlinear"));
        let do_displace = dofs > 1 && any_proj && ctx.params().logical("displace mesh at init").unwrap_or(true);
        let mover = mover().filter(|_| do_displace);
        let base: Vec<[f64; 3]> = (0..ctx.mesh.n_nodes() as u32).map(|n| ctx.mesh.coord(n)).collect();
        let moved = do_displace.then(|| displaced(&base, &values, 1.0, &perm, dofs));
        if let (Some(m), Some(x)) = (&mover, &moved) {
            ctx.info("DefaultInitialize", "Displacing mesh for nonlinear projectors", 8);
            m(x);
        }
        let view_values = if do_displace && mover.is_none() { values.clone() } else { vec![0.0; values.len()] };
        let disp = DisplacedMesh { mesh: ctx.mesh, values: &view_values, perm: &perm, dofs };
        let model = &ctx.sim.model;
        let bcs: Vec<&dyn Params> = model.boundary_conditions.iter().map(|b| b as &dyn Params).collect();
        let bc_of = |tag: i32| model.bc_of_tag(tag);
        let env = ProjectorEnv {
            mesh: &disp,
            bcs: &bcs,
            bc_of: &bc_of,
            simulation: &model.simulation,
            solver: Some(ctx.params()),
            coord_dim: self.coord_dim,
            log: ctx.log,
            discont_perm: None,
        };
        let timestep = ctx.sim.vars.scalar("timestep").map_or(0, |t| t.round() as i32);
        let projectors = match generate_mortar_projectors(&env, &[], true, &self.case_dir, timestep).map_err(mortar_err)? {
            Some((p, false)) => p,
            Some((_, true)) => return Err(unsupported("Save Projector And Stop")),
            None => vec![None; bcs.len()],
        };
        // FlatInterfaceMeshes re-detects the coordinate at each generation
        // and stores it in the boundary condition (CU:260 reads it back).
        let mut flat_dim = vec![0usize; bcs.len()];
        for (i, p) in projectors.iter().enumerate() {
            if p.is_none() || !flag(bcs[i], "flat projector") {
                continue;
            }
            let target = bcs[i].integer("mortar bc").unwrap_or(0).max(0) as usize;
            let split = model.simulation.logical("mortar bcs split quadratic").unwrap_or(false);
            if let Some((mut b1, mut b2)) =
                create_interface_meshes(&disp, &bc_of, bcs.len(), i + 1, target, split).map_err(mortar_err)?
            {
                let map: Option<Vec<i32>> = match bcs[i].get("projector coordinate mapping") {
                    Some(elmer_core::Value::Integer(v)) => Some(v.clone()),
                    _ => None,
                };
                map_interface_coordinate(&mut b1, &mut b2, map.as_deref()).map_err(mortar_err)?;
                flat_dim[i] = crate::level::flat_interface_meshes(&mut b1, &mut b2, bcs[i], self.coord_dim, &NullLogger);
            }
        }
        if let (Some(m), Some(x)) = (&mover, &moved) {
            m(&displaced(x, &values, -1.0, &perm, dofs));
        }
        let id = ctx.solver_id();
        let mut state = self.state.lock().expect("contact state");
        let s = state.entry(id).or_default();
        s.projectors = projectors;
        s.flat_dim = flat_dim;
        Ok(())
    }

    fn finish_boundary(&self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let id = ctx.solver_id();
        let mut state = self.state.lock().expect("contact state");
        let s = state.get_mut(&id).ok_or_else(|| fatal("Cannot apply contact without projectors!"))?;
        let mut first = self.first_time.lock().expect("first time");
        determine_contact(ctx, s, *first)?;
        *first = false;
        Ok(())
    }

    fn before_solve(&self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let id = ctx.solver_id();
        let state = self.state.lock().expect("contact state");
        let Some(s) = state.get(&id) else { return Ok(()) };
        let projectors: Vec<Option<Projector>> = s
            .projectors
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let mut p = p.clone()?;
                let m = s.bcs.get(&(i + 1))?;
                p.active = Some(m.active.clone());
                p.rhs = Some(m.rhs.clone());
                p.mortar_perm = true;
                p.master_scale = m.master_scale;
                Some(p)
            })
            .collect();
        let model = &ctx.sim.model;
        let bcs: Vec<&dyn Params> = model.boundary_conditions.iter().map(|b| b as &dyn Params).collect();
        let name = ctx.variable_name().ok_or_else(|| fatal("solver has no variable"))?;
        let f = ctx.primary().ok_or_else(|| fatal("solver has no variable"))?;
        let perm = f.perm.ok_or_else(|| fatal("solver variable has no permutation"))?;
        let target = ConstraintTarget {
            solver: ctx.params(),
            name,
            dofs: f.dofs,
            dim: self.coord_dim,
            perm,
            n_rows: f.values.len(),
        };
        let c = generate_constraint_matrix(&bcs, &projectors, &target, ctx.log).map_err(mortar_err)?;
        // `Save Constraint Matrix`: the last generation's.
        if let (Some(c), true) = (&c, flag(ctx.params(), "save constraint matrix")) {
            crate::save::save_constraint_matrix(c, &self.case_dir.join("cm"))
                .map_err(|e| fatal(format!("Save Constraint Matrix: {e}")))?;
        }
        let n_cols = target.n_rows;
        drop(state);
        let c = c.map(|c| {
            let n = c.n_rows();
            elmer_solver::ConstraintMatrix {
                matrix: elmer_core::CsrMatrix {
                    nrows: n,
                    ncols: n_cols + n,
                    row_ptr: c.rows,
                    col_idx: c.cols,
                    values: c.values,
                    diag: None,
                },
                rhs: c.rhs,
                inv_perm: c.inv_perm.into_iter().map(|i| i as u32).collect(),
            }
        });
        ctx.set_generated_constraints(c);
        Ok(())
    }
}

/// The contact field names (`GetContactFields`, CU:729-748), in the order
/// they are created.
const FIELDS: [&str; 8] = [
    "contact distance",
    "contact gap",
    "contact normalload",
    "contact slipload",
    "contact weight",
    "contact active",
    "contact stick",
    "lagrange multiplier",
];

/// `DetermineContact( Solver )` for the ported subset. The loops index the
/// projector rows and several fields together, as the Fortran does.
#[allow(clippy::needless_range_loop)]
fn determine_contact(ctx: &mut SolverCtx<'_>, s: &mut SolverContact, first_time: bool) -> Result<(), SolverError> {
    ctx.info(CALLER, "Setting up contact conditions", 8);
    let name = ctx.variable_name().ok_or_else(|| fatal("solver has no variable"))?.to_string();
    let (dofs, field_perm, values) = {
        let f = ctx.primary().ok_or_else(|| fatal("solver has no variable"))?;
        (f.dofs, f.perm.ok_or_else(|| fatal("solver variable has no permutation"))?.to_vec(), f.values.to_vec())
    };
    let dim = ctx.mesh.dim();
    let cdofs = if dofs == dim {
        dim
    } else if dofs == dim + 1 {
        ctx.info(CALLER, "We seem to have mixed formulation, ignoring pressure!", 7);
        dim
    } else {
        return Err(fatal(format!("Invalid number of dofs for contact problem: {dofs}")));
    };
    let load_eps = ctx.params().real("limiter load tolerance").unwrap_or(f64::EPSILON);
    let val_eps = ctx.params().real("limiter value tolerance").unwrap_or(f64::EPSILON);
    let n_nodes = ctx.mesh.n_nodes();

    // GetContactFields: on the first call the fields over the nodes of the
    // mortar boundaries and their targets.
    if ctx.sim.vars.id(&format!("{name} {}", FIELDS[0])).is_none() {
        ctx.info(CALLER, "Creating contact fields", 8);
        let model = &ctx.sim.model;
        let mut active_bcs = vec![false; model.boundary_conditions.len() + 1];
        for (i, bc) in model.boundary_conditions.iter().enumerate() {
            let j = bc.integer("mortar bc").or_else(|| bc.integer("contact bc")).unwrap_or(0);
            if j > 0 {
                active_bcs[i + 1] = true;
                if (j as usize) < active_bcs.len() {
                    active_bcs[j as usize] = true;
                }
            }
        }
        let mut bperm = vec![0u32; field_perm.len()];
        let nb = ctx.mesh.n_bulk() as u32;
        for e in (nb..nb + ctx.mesh.n_boundary() as u32).map(ElemId) {
            if ctx.bc_id(e).is_some_and(|b| active_bcs[b]) {
                for &n in ctx.mesh.element(e).nodes {
                    bperm[n as usize] = 1;
                }
            }
        }
        let mut j = 0u32;
        for p in bperm.iter_mut() {
            if *p > 0 {
                j += 1;
                *p = j - 1;
            } else {
                *p = NO_DOF;
            }
        }
        let perm: Arc<[u32]> = bperm.into();
        for f in FIELDS {
            let var = Variable::new(&format!("{name} {f}"), FieldKind::Nodal, 1, Some(perm.clone()), vec![0.0; j as usize]);
            ctx.sim.vars.add(var)?;
        }
    }
    let bperm: Vec<u32> = ctx
        .sim
        .vars
        .get(&format!("{name} {}", FIELDS[0]))
        .and_then(|v| v.perm.as_deref().map(<[u32]>::to_vec))
        .ok_or_else(|| fatal("contact fields have no permutation"))?;
    let mut fields: Vec<Vec<f64>> = FIELDS
        .iter()
        .map(|f| ctx.sim.vars.get(&format!("{name} {f}")).map(|v| v.values.clone()).unwrap_or_default())
        .collect();
    let [dist_v, gap_v, nload_v, sload_v, weight_v, nact_v, stick_v, _lag_v] = &mut fields[..] else { unreachable!() };
    nact_v.iter_mut().for_each(|v| *v = -1.0);
    stick_v.iter_mut().for_each(|v| *v = -1.0);

    // CalculateContactLoad: A x - b of the assembled system.
    let load = {
        let sys = ctx.system().ok_or_else(|| fatal("no assembled system"))?;
        let mut r = vec![0.0; sys.matrix.nrows];
        sys.matrix.apply(&values[..sys.matrix.nrows], &mut r);
        for (ri, bi) in r.iter_mut().zip(&sys.rhs) {
            *ri -= bi;
        }
        r
    };
    if let Some(v) = ctx.sim.vars.get_mut(&format!("{name} contact load")) {
        let m = v.values.len().min(load.len());
        v.values[..m].copy_from_slice(&load[..m]);
        v.values_changed = true;
    }

    let n_bcs = ctx.sim.model.boundary_conditions.len();
    for bc_ind in 1..=n_bcs {
        let Some(projector) = s.projectors.get(bc_ind - 1).and_then(Option::as_ref) else { continue };
        ctx.info(CALLER, &format!("Set contact for boundary: {bc_ind}"), 8);
        let bc = &ctx.sim.model.boundary_conditions[bc_ind - 1];
        if !flag(bc, "flat projector") {
            return Err(fatal("Projector must be current either flat, plane, cylindrical or rotational!"));
        }
        // CU:259-261: ActiveDirection from the stored flat coordinate.
        let active_dir = match s.flat_dim.get(bc_ind - 1).copied().filter(|&d| d > 0) {
            Some(d) => d,
            None => bc.integer("flat projector coordinate").map_or(cdofs, |d| d as usize),
        };
        let ctype = contact_type(bc, ctx)?;
        let dof_n = active_dir;
        let dof_t1 = (1..=cdofs).find(|&i| i != dof_n).unwrap_or(0);
        let mut normal = [0.0; 3];
        normal[active_dir - 1] = 1.0;

        // InitializeMortarVectors.
        let rows = projector.n_rows();
        let totsize = cdofs * rows;
        let mut perm = vec![0u32; field_perm.len()];
        for (i, &j) in projector.inv_perm.iter().enumerate() {
            if j != NO_DOF {
                perm[j as usize] = i as u32 + 1;
            }
        }
        let m = s.bcs.entry(bc_ind).or_insert_with(|| {
            ctx.info(CALLER, &format!("Allocating projector mortar vectors of size: {totsize}"), 10);
            MortarBc { active: vec![false; totsize], rhs: vec![0.0; totsize], perm: perm.clone(), master_scale: 1.0 }
        });
        if m.perm != perm {
            if m.rhs.len() != totsize {
                m.rhs = vec![0.0; totsize];
            }
            let mut active = vec![false; totsize];
            for (i, &j) in perm.iter().enumerate() {
                let k = m.perm[i];
                if j == 0 || k == 0 {
                    continue;
                }
                for l in 0..cdofs {
                    active[cdofs * (j as usize - 1) + l] = m.active[cdofs * (k as usize - 1) + l];
                }
            }
            m.active = active;
            m.perm = perm;
            ctx.info(CALLER, "Copied > Active < flag to changed projector", 8);
        }
        m.master_scale = 1.0;

        // CalculateContactPressure: the slave side's normal and slip loads.
        let mut node_done = vec![false; field_perm.len()];
        let (mut sign0, mut sign, mut count) = (0i32, 0i32, 0usize);
        let mut st = MetricStatus::default();
        let mut be = BasisEval::default();
        let nb = ctx.mesh.n_bulk() as u32;
        let cartesian = !ctx
            .simulation_params()
            .string("coordinate system")
            .is_some_and(|c| c.to_ascii_lowercase().contains("symmetric") || c.to_ascii_lowercase().contains("cylindric"));
        if !cartesian {
            return Err(unsupported("contact in non-Cartesian coordinates"));
        }
        for e in (nb..nb + ctx.mesh.n_boundary() as u32).map(ElemId) {
            if ctx.bc_id(e) != Some(bc_ind) {
                continue;
            }
            let el = ctx.mesh.element(e);
            let def = element_type(el.etype.0).map_err(|x| fatal(x.to_string()))?;
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let rule = gauss_points(el.etype, None).map_err(|x| fatal(x.to_string()))?;
            for t in 0..rule.s.len() {
                let uvw = [rule.u[t], rule.v[t], rule.w[t]];
                element_info(def, &xyz, uvw, ctx.mesh.dim(), &mut st, &mut be).map_err(|x| fatal(x.to_string()))?;
                let sw = be.det_j * rule.s[t];
                let nrm = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [uvw[0], uvw[1]])?;
                let dot = nrm[0] * normal[0] + nrm[1] * normal[1] + nrm[2] * normal[2];
                sign = if dot < 0.0 { 1 } else { -1 };
                if sign0 == 0 {
                    sign0 = sign;
                } else if sign0 != sign {
                    count += 1;
                }
                for (i, &node) in el.nodes.iter().enumerate() {
                    let j = bperm[node as usize];
                    if j == NO_DOF {
                        continue;
                    }
                    let j = j as usize;
                    if !node_done[node as usize] {
                        node_done[node as usize] = true;
                        weight_v[j] = 0.0;
                        nload_v[j] = 0.0;
                        sload_v[j] = 0.0;
                    }
                    let k = field_perm[node as usize];
                    if k == NO_DOF {
                        continue;
                    }
                    let mut force = [0.0; 3];
                    for (l, fl) in force.iter_mut().enumerate().take(cdofs) {
                        *fl = load[dofs * k as usize + l];
                    }
                    let nf = force[0] * nrm[0] + force[1] * nrm[1] + force[2] * nrm[2];
                    let sf = (force[0].powi(2) + force[1].powi(2) + force[2].powi(2) - nf.powi(2)).max(0.0).sqrt();
                    nload_v[j] -= sw * be.basis[i] * nf;
                    sload_v[j] += sw * be.basis[i] * sf;
                    weight_v[j] += sw * be.basis[i];
                }
            }
        }
        for (node, done) in node_done.iter().enumerate() {
            if !done {
                continue;
            }
            let j = bperm[node];
            if j == NO_DOF {
                continue;
            }
            let w = weight_v[j as usize];
            if w.is_nan() {
                continue;
            }
            if w.abs() > f64::EPSILON {
                sload_v[j as usize] /= w.powi(2);
                nload_v[j as usize] /= w.powi(2);
            }
        }
        let dist_sign = if count == 0 {
            ctx.info(CALLER, "All normals are consistently signed", 10);
            sign
        } else {
            ctx.warn(CALLER, &format!("There are normals with conflicting signs: {count}"));
            1
        };
        ctx.info(CALLER, &format!("Normal direction for distance measure: {dist_sign}"), 8);

        // CalculateMortarDistance (one pass: no dual projector).
        let mut slave = vec![false; field_perm.len()];
        for e in (nb..nb + ctx.mesh.n_boundary() as u32).map(ElemId) {
            if ctx.bc_id(e) == Some(bc_ind) {
                for &n in ctx.mesh.element(e).nodes {
                    slave[n as usize] = true;
                }
            }
        }
        for i in 0..rows {
            let j = projector.inv_perm[i];
            if j == NO_DOF {
                continue;
            }
            let (mut wsum, mut wsum_m) = (0.0f64, 0.0f64);
            let (mut dist, mut gap) = (0.0f64, 0.0f64);
            if ctype == ContactType::Slip {
                for (k, coeff) in projector.row(i) {
                    if field_perm[k as usize] == NO_DOF {
                        continue;
                    }
                    if slave[k as usize] {
                        wsum += coeff;
                    } else {
                        wsum_m += coeff;
                    }
                }
                if wsum.abs() <= f64::MIN_POSITIVE {
                    return Err(SolverError::fatal("CalculateMortarDistance", "wsum seems to be almost zero!"));
                }
                if wsum_m.abs() <= f64::MIN_POSITIVE {
                    return Err(SolverError::fatal("CalculateMortarDistance", "wsumM seems to be almost zero!"));
                }
                let mult = (wsum / wsum_m).abs();
                for (k, coeff) in projector.row(i) {
                    let l = field_perm[k as usize];
                    if l == NO_DOF {
                        continue;
                    }
                    let coeff = if slave[k as usize] { coeff } else { mult * coeff };
                    let l = l as usize;
                    let disp = [values[dofs * l], values[dofs * l + 1], if cdofs == 2 { 0.0 } else { values[dofs * l + 2] }];
                    let mut c = ctx.mesh.coord(k);
                    gap += coeff * (normal[0] * c[0] + normal[1] * c[1] + normal[2] * c[2]);
                    for d in 0..3 {
                        c[d] += disp[d];
                    }
                    dist += coeff * (normal[0] * c[0] + normal[1] * c[1] + normal[2] * c[2]);
                }
                gap /= wsum;
                dist = dist_sign as f64 * dist / wsum;
            }
            // Label 200.
            m.rhs[cdofs * i + dof_n - 1] = -gap;
            if ctype == ContactType::Tie {
                m.rhs[cdofs * i + dof_t1 - 1] = -0.0;
            }
            let jj = bperm[j as usize] as usize;
            dist_v[jj] = dist;
            gap_v[jj] = gap;
        }

        // NormalContactSet.
        let bc = &ctx.sim.model.boundary_conditions[bc_ind - 1];
        let mut found = None;
        if first_time {
            found = bc.real("mortar bc initial contact depth").or_else(|| bc.real("contact depth offset initial"));
        }
        let offset = found.or_else(|| bc.real("contact depth offset")).unwrap_or(0.0);
        let (mut added, mut removed) = (0usize, 0usize);
        for i in 0..rows {
            let j = projector.inv_perm[i];
            if j == NO_DOF || field_perm[j as usize] == NO_DOF {
                continue;
            }
            let k = bperm[j as usize] as usize;
            let ind = cdofs * i + dof_n - 1;
            if ctype == ContactType::Tie {
                m.active[ind] = true;
                continue;
            }
            if m.active[ind] {
                if -nload_v[k] > -load_eps {
                    removed += 1;
                    m.active[ind] = false;
                }
            } else if dist_v[k] < -val_eps + offset {
                added += 1;
                m.active[ind] = true;
            }
        }
        if added > 0 {
            ctx.info(CALLER, &format!("Added {added} nodes to the set"), 6);
        }
        if removed > 0 {
            ctx.info(CALLER, &format!("Removed {removed} nodes from the set"), 6);
        }

        // TangentContactSet: slip frees the tangent, tie follows the normal.
        for i in 0..rows {
            let n_active = m.active[cdofs * i + dof_n - 1];
            for t in (1..=cdofs).filter(|&t| t != dof_n) {
                m.active[cdofs * i + t - 1] = ctype == ContactType::Tie && n_active;
            }
        }
        for i in 0..rows {
            let j = projector.inv_perm[i];
            if j == NO_DOF {
                continue;
            }
            let k = bperm[j as usize] as usize;
            nact_v[k] = if m.active[cdofs * i + dof_n - 1] { 1.0 } else { -1.0 };
            stick_v[k] = if m.active[cdofs * i + dof_t1 - 1] { 1.0 } else { -1.0 };
        }

        // ProjectFromSlaveToMaster.
        let mut coeff_table = vec![0.0f64; dist_v.len()];
        let mut done = vec![false; dist_v.len()];
        for i in 0..rows {
            let ip = projector.inv_perm[i];
            if ip == NO_DOF {
                continue;
            }
            let l = bperm[ip as usize] as usize;
            if l >= n_nodes {
                continue;
            }
            for (k, coeff) in projector.row(i) {
                if k as usize >= n_nodes || field_perm[k as usize] == NO_DOF || slave[k as usize] {
                    continue;
                }
                let l2 = bperm[k as usize] as usize;
                if !done[l2] {
                    dist_v[l2] = 0.0;
                    gap_v[l2] = 0.0;
                    nact_v[l2] = 0.0;
                    stick_v[l2] = 0.0;
                    nload_v[l2] = 0.0;
                    sload_v[l2] = 0.0;
                    done[l2] = true;
                }
                coeff_table[l2] += coeff;
                dist_v[l2] += coeff * dist_v[l];
                gap_v[l2] += coeff * gap_v[l];
                nact_v[l2] += coeff * nact_v[l];
                stick_v[l2] += coeff * stick_v[l];
                nload_v[l2] += coeff * nload_v[l];
                sload_v[l2] += coeff * sload_v[l];
            }
        }
        let eps = 1.0e-8 * coeff_table.iter().fold(0.0f64, |a, &c| a.max(c.abs()));
        for i in 0..coeff_table.len() {
            if done[i] && coeff_table[i].abs() > eps {
                let c = coeff_table[i];
                dist_v[i] /= c;
                gap_v[i] /= c;
                nact_v[i] /= c;
                stick_v[i] /= c;
                if nact_v[i] >= 0.0 {
                    nload_v[i] /= c;
                    sload_v[i] /= c;
                } else {
                    nload_v[i] = 0.0;
                    sload_v[i] = 0.0;
                }
            }
        }
    }

    for (f, vals) in FIELDS.iter().zip(fields) {
        if let Some(v) = ctx.sim.vars.get_mut(&format!("{name} {f}")) {
            v.values = vals;
            v.values_changed = true;
        }
    }
    Ok(())
}
