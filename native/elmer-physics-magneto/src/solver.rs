//! The `elmer-solver` [`Solver`] adapter: `Procedure = "MagnetoDynamics"
//! "WhitneyAVSolver"` and `"MagnetoDynamics" "WhitneyAVHarmonicSolver"`.
//!
//! The simulation creates the primary variable over the nodes; Elmer's
//! edge solver has a `Perm` over nodes and edges (`Element = "n:1 e:1"`),
//! so the adapter replaces the variable's perm and values with the
//! layout of [`AvLayout`](crate::layout::AvLayout) after each solve: perm
//! entries `0..n_nodes` are the nodal slots, entry `n_nodes + e` is edge
//! `e`'s slot. The steady solver then sets `Norm Permutation = nNodes+1`
//! as `WhitneyAVSolver` does, and the change is computed by
//! `elmer-solver` (`ComputeChange`).

use std::sync::Arc;

use elmer_core::{ElemId, MeshView, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::case::Case;
use crate::error::MagnetoError;

/// The steady solver's `(library, entry)`.
pub const STEADY: (&str, &str) = ("MagnetoDynamics", "WhitneyAVSolver");
/// The time-harmonic solver's `(library, entry)`.
pub const HARMONIC: (&str, &str) = ("MagnetoDynamics", "WhitneyAVHarmonicSolver");

/// `MagnetoDynamicsCalcFields`' `(library, entry)`.
pub const CALC_FIELDS: (&str, &str) = ("MagnetoDynamics", "MagnetoDynamicsCalcFields");

/// Registers the procedures of this crate.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(STEADY.0, STEADY.1, steady_factory);
    registry.register(HARMONIC.0, HARMONIC.1, harmonic_factory);
    registry.register(CALC_FIELDS.0, CALC_FIELDS.1, calc_fields_factory);
    crate::mgdyn2d::register(registry);
    crate::mgdyn2d_steady::register(registry);
    crate::statmag::register(registry);
}

pub fn calc_fields_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(CalcFields))
}

/// Whether `params` runs `WhitneyAVSolver` (`Some(false)`) or
/// `WhitneyAVHarmonicSolver` (`Some(true)`).
pub fn whitney_kind(params: &ParamMap) -> Option<bool> {
    let (lib, entry) = elmer_solver::solver::procedure_of(params)?;
    if !lib.eq_ignore_ascii_case(STEADY.0) {
        return None;
    }
    if entry.eq_ignore_ascii_case(HARMONIC.1) {
        Some(true)
    } else if entry.eq_ignore_ascii_case(STEADY.1) {
        Some(false)
    } else {
        None
    }
}

/// The `Element` of `WhitneyAVSolver_Init0` / `WhitneyAVHarmonicSolver_Init0`
/// for the lowest-order (non-Piola) Whitney basis when the section gives
/// none: `n:1 e:1`, or `n:0 e:1` in steady state without `Static
/// Conductivity` (explicit, else set by a velocity body force or an
/// `Electric Current Density` BC). With `Use Piola Transform` alone (steady
/// state, no Lagrange gauge) the first-order Piola one of
/// [`crate::piola::init_element`]. Sections with the options that pick
/// other bases are left alone (the solvers refuse them).
pub fn whitney_init_section(params: &mut ParamMap, body_forces: &[ParamMap], bcs: &[ParamMap]) {
    let Some(harmonic) = whitney_kind(params) else { return };
    if params.get("element").is_some() {
        return;
    }
    let prefix = |s: &ParamMap, p: &str| s.0.keys().any(|k| k.starts_with(p));
    let static_cond = params.logical("static conductivity").unwrap_or_else(|| {
        body_forces.iter().any(|b| prefix(b, "angular velocity") || prefix(b, "lorentz velocity"))
            || bcs.iter().any(|b| prefix(b, "electric current density"))
    });
    let on = |k: &str| params.logical(k) == Some(true);
    if on("use tree gauge") || on("quadratic approximation") || on("second kind basis") {
        return;
    }
    if on("use piola transform") {
        if !harmonic && !on("use lagrange gauge") {
            params.insert("Element", Value::Str(crate::piola::init_element(static_cond).into()));
        }
        return;
    }
    let e = if harmonic || static_cond { "n:1 e:1" } else { "n:0 e:1" };
    params.insert("Element", Value::Str(e.into()));
}

pub fn steady_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(WhitneyAv { harmonic: false }))
}

pub fn harmonic_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(WhitneyAv { harmonic: true }))
}

impl From<MagnetoError> for SolverError {
    fn from(e: MagnetoError) -> Self {
        match e {
            MagnetoError::Unsupported(m) => SolverError::Unsupported(m),
            MagnetoError::Input(m) => SolverError::fatal("WhitneyAVSolver", m),
            MagnetoError::Solver(e) => e,
            MagnetoError::Linalg(e) => SolverError::Linalg(e),
        }
    }
}

/// Layout slot -> Elmer slot: each active element's `AvLayout` DOFs (nodal
/// V, then its edges) against `GetElementDOFs` of the framework layout.
fn elmer_order(case: &Case, id: usize, harmonic: bool, d: &dyn elmer_solver::SolverDofs) -> Result<Vec<u32>, SolverError> {
    let sp = case.solver(id).ok_or_else(|| SolverError::fatal("WhitneyAVSolver", "no solver section"))?;
    let with_v = harmonic || crate::steady::static_conductivity(case, sp);
    let lay = crate::layout::AvLayout::build(case, id, with_v)?;
    let refuse = |what: &str| SolverError::Unsupported(format!("WhitneyAVSolver: Elmer's DOF layout {what}"));
    if d.perm().n_slots != lay.n() {
        return Err(refuse("has another number of DOFs"));
    }
    let mut order = vec![u32::MAX; lay.n()];
    let mut slots = Vec::new();
    for (k, &e) in lay.active.iter().enumerate() {
        d.element_slots(ElemId(e as u32), &mut slots);
        let own = &lay.elem_dofs[k];
        if slots.len() != own.len() {
            return Err(refuse("has another element DOF count"));
        }
        for (&a, &f) in own.iter().zip(&slots) {
            if f == NO_DOF {
                return Err(refuse("lacks an element DOF"));
            }
            let o = &mut order[a as usize];
            if *o == u32::MAX {
                *o = f;
            } else if *o != f {
                return Err(refuse("numbers the element DOFs differently"));
            }
        }
    }
    if order.contains(&u32::MAX) {
        return Err(refuse("leaves a DOF unnumbered"));
    }
    Ok(order)
}

pub struct WhitneyAv {
    harmonic: bool,
}

impl Case {
    /// The case a solver context describes: bulk elements carry their body
    /// id and boundary elements their boundary condition id, mapped from
    /// the mesh tags through `Target Bodies` / `Target Boundaries` as
    /// `elmer-solver` does.
    pub fn from_view(mesh: &dyn MeshView, model: &elmer_solver::Model) -> Self {
        let n_bulk = mesh.n_bulk();
        let bulk = (0..n_bulk)
            .map(|e| {
                let r = mesh.element(ElemId(e as u32));
                elmer_mesh::Elem {
                    code: r.etype.0,
                    body: model.body_of_tag(r.tag).map_or(0, |b| b as i32),
                    nodes: r.nodes.to_vec(),
                    tag: r.tag,
                    edges: Vec::new(),
                }
            })
            .collect();
        let boundary = (n_bulk..n_bulk + mesh.n_boundary())
            .map(|e| {
                let r = mesh.element(ElemId(e as u32));
                elmer_mesh::BoundElem {
                    code: r.etype.0,
                    constraint: model.bc_of_tag(r.tag).map_or(0, |b| b as i32),
                    parents: r.parents.map(|p| p.map(|p| p.0 as usize).filter(|&p| p < n_bulk)),
                    nodes: r.nodes.to_vec(),
                    tag: r.tag,
                    edges: Vec::new(),
                }
            })
            .collect();
        let m = elmer_mesh::Mesh {
            dim: mesh.dim(),
            max_dim: mesh.dim(),
            nodes: (0..mesh.n_nodes() as u32).map(|i| mesh.coord(i)).collect(),
            node_tags: (1..=mesh.n_nodes() as i32).collect(),
            bulk,
            boundary,
            edges: None,
        };
        let some = |v: &[ParamMap]| v.iter().cloned().map(Some).collect();
        Self {
            mesh: m,
            simulation: model.simulation.clone(),
            constants: model.constants.clone(),
            bodies: some(&model.bodies),
            materials: some(&model.materials),
            body_forces: some(&model.body_forces),
            equations: some(&model.equations),
            bcs: some(&model.boundary_conditions),
            solvers: some(&model.solvers),
            components: Vec::new(),
        }
    }
}

impl Solver for WhitneyAv {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(SolverError::Unsupported("transient WhitneyAVSolver".into()));
        }
        if !self.harmonic && crate::piola::is_piola(ctx.params()) {
            return crate::piola::solve(ctx);
        }
        let case = Case::from_view(ctx.mesh, &ctx.sim.model);
        let id = ctx.solver_id();
        // With `Element = n:1 e:1` (or `n:0 e:1`) from elmer-run, the
        // framework holds Elmer's Perm of the solver (InitialPermutation +
        // OptimizeBandwidth over nodes and mesh edges): the system is solved
        // and stored in that slot order.
        let order = match ctx.solver_dofs() {
            Some(d) => Some(elmer_order(&case, id, self.harmonic, d)?),
            None => None,
        };
        let (values, node_slot, n_nodal, n_edges, dofs) = if self.harmonic {
            let s = crate::harmonic::solve_harmonic_ordered(&case, id, order.as_deref())?;
            let v = s.values.iter().flat_map(|z| [z.re, z.im]).collect();
            (v, s.node_slot, s.n_nodal, s.n_edges, 2)
        } else {
            let s = crate::steady::solve_steady_ordered(&case, id, order.as_deref())?;
            (s.values, s.node_slot, s.n_nodal, s.n_edges, 1)
        };
        let n_nodes = node_slot.len();
        let (values, perm): (Vec<f64>, Option<Arc<[u32]>>) = match &order {
            Some(order) => {
                let mut v = vec![0.0; values.len()];
                for (i, &o) in order.iter().enumerate() {
                    for k in 0..dofs {
                        v[dofs * o as usize + k] = values[dofs * i + k];
                    }
                }
                (v, None)
            }
            None => (
                values,
                Some(
                    node_slot
                        .iter()
                        .copied()
                        .chain((0..n_edges as u32).map(|e| n_nodal as u32 + e))
                        .collect(),
                ),
            ),
        };
        let ix = ctx
            .sim
            .slots
            .iter()
            .position(|s| s.id == id)
            .ok_or_else(|| SolverError::fatal("WhitneyAVSolver", "solver slot not found"))?;
        let var_id = ctx.sim.slots[ix]
            .variable
            .ok_or_else(|| SolverError::fatal("WhitneyAVSolver", "solver has no variable"))?;
        let var = ctx.sim.vars.by_id_mut(var_id);
        if var.dofs != dofs {
            return Err(SolverError::fatal(
                "WhitneyAVSolver",
                format!("variable '{}' has {} components, the solver needs {dofs}", var.name, var.dofs),
            ));
        }
        if let Some(p) = &perm {
            var.perm = Some(p.clone());
        }
        var.values = values;
        var.values_changed = true;
        if !self.harmonic && n_nodal == 0 {
            // WhitneyAVSolver_Init0: `Variable Output = False` without the
            // scalar potential in steady state.
            var.output = false;
        }
        let slot = &mut ctx.sim.slots[ix];
        if let Some(p) = perm {
            slot.perm = Some(p);
        }
        if !self.harmonic {
            slot.params.insert("Norm Permutation", Value::Integer(vec![n_nodes as i32 + 1]));
        }
        debug_assert!(node_slot.iter().all(|&s| s == NO_DOF || (s as usize) < n_nodal));
        ctx.nonlinear_change()?;
        Ok(())
    }
}

/// `MagnetoDynamicsCalcFields` through the simulation: the potential is
/// the target solver's variable (perm over nodes then edges, as
/// [`WhitneyAv`] leaves it); the projected fields become exported nodal
/// variables (`Magnetic Flux Density` with components `... 1..3`, ...)
/// and the solver variable `hr_dummy` (`-nooutput`) holds the last
/// component solved, whose norm `ComputeChange` then reports.
pub struct CalcFields;

impl Solver for CalcFields {
    /// `MagnetoDynamicsCalcFields` (`calcfields2d`): the target solver as
    /// `_Init0` finds it, nodal 2D or edge-element 3D.
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        let t = crate::calcfields2d::target(&ctx.sim.model.solvers, ctx.params())?;
        crate::calcfields2d::solve(ctx, t, transient)
    }
}
