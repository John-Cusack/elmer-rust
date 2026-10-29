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

use elmer_core::{ElemId, FieldKind, MeshView, ParamMap, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable, VariableStore};

use crate::case::Case;
use crate::error::MagnetoError;

/// The steady solver's `(library, entry)`.
pub const STEADY: (&str, &str) = ("MagnetoDynamics", "WhitneyAVSolver");
/// The time-harmonic solver's `(library, entry)`.
pub const HARMONIC: (&str, &str) = ("MagnetoDynamics", "WhitneyAVHarmonicSolver");

/// `MagnetoDynamicsCalcFields`' `(library, entry)`.
pub const CALC_FIELDS: (&str, &str) = ("MagnetoDynamics", "MagnetoDynamicsCalcFields");

/// Registers the three procedures.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(STEADY.0, STEADY.1, steady_factory);
    registry.register(HARMONIC.0, HARMONIC.1, harmonic_factory);
    registry.register(CALC_FIELDS.0, CALC_FIELDS.1, calc_fields_factory);
}

pub fn calc_fields_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(CalcFields))
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
        let case = Case::from_view(ctx.mesh, &ctx.sim.model);
        let id = ctx.solver_id();
        let (values, node_slot, n_nodal, n_edges, dofs) = if self.harmonic {
            let s = crate::solve_harmonic(&case, id)?;
            let v = s.values.iter().flat_map(|z| [z.re, z.im]).collect();
            (v, s.node_slot, s.n_nodal, s.n_edges, 2)
        } else {
            let s = crate::solve_steady(&case, id)?;
            (s.values, s.node_slot, s.n_nodal, s.n_edges, 1)
        };
        let n_nodes = node_slot.len();
        let perm: Arc<[u32]> = node_slot
            .iter()
            .copied()
            .chain((0..n_edges as u32).map(|e| n_nodal as u32 + e))
            .collect();
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
        var.perm = Some(perm.clone());
        var.values = values;
        var.values_changed = true;
        if !self.harmonic && n_nodal == 0 {
            // WhitneyAVSolver_Init0: `Variable Output = False` without the
            // scalar potential in steady state.
            var.output = false;
        }
        let slot = &mut ctx.sim.slots[ix];
        slot.perm = Some(perm);
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

/// Adds or overwrites a nodal variable.
fn upsert(
    vars: &mut VariableStore,
    name: &str,
    dofs: usize,
    perm: &Arc<[u32]>,
    values: Vec<f64>,
    output: bool,
) -> Result<elmer_solver::VarId, SolverError> {
    if let Some(v) = vars.get_mut(name) {
        if v.dofs != dofs {
            return Err(SolverError::fatal(
                "MagnetoDynamicsCalcFields",
                format!("variable '{name}' exists with {} components", v.dofs),
            ));
        }
        v.perm = Some(perm.clone());
        v.values = values;
        v.values_changed = true;
        return Ok(vars.id(name).expect("just found"));
    }
    let mut v = Variable::new(name, FieldKind::Nodal, dofs, Some(perm.clone()), values);
    v.output = output;
    vars.add(v)
}

impl Solver for CalcFields {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(SolverError::Unsupported("transient MagnetoDynamicsCalcFields".into()));
        }
        let case = Case::from_view(ctx.mesh, &ctx.sim.model);
        let id = ctx.solver_id();
        let sp = case
            .solver(id)
            .ok_or_else(|| SolverError::fatal("MagnetoDynamicsCalcFields", "no solver section"))?;
        let (pot_id, harmonic) = crate::calcfields::target_solver(&case, sp)?;
        if harmonic {
            return Err(SolverError::Unsupported("CalcFields of a harmonic potential".into()));
        }
        let psp = case.solver(pot_id).expect("target solver exists");
        let layout = crate::layout::AvLayout::build(&case, pot_id, crate::steady::static_conductivity(&case, psp))?;
        let pot_var = ctx
            .sim
            .slots
            .iter()
            .find(|s| s.id == pot_id)
            .and_then(|s| s.variable)
            .ok_or_else(|| SolverError::fatal("MagnetoDynamicsCalcFields", "potential variable not found"))?;
        let res = {
            let pv = ctx.sim.vars.by_id(pot_var);
            if pv.values.len() != layout.n() || pv.dofs != 1 {
                return Err(SolverError::fatal(
                    "MagnetoDynamicsCalcFields",
                    "the potential has not been solved on the edge layout yet",
                ));
            }
            let pot = crate::Potential {
                edges: &layout.edges,
                n_nodal: layout.n_nodal,
                dofs: 1,
                values: &pv.values,
            };
            crate::calc_fields(&case, id, &pot)?
        };
        let perm: Arc<[u32]> = res.node_slot.iter().copied().collect();
        // `hr_dummy` keeps the last component solved, or the one `Show Norm
        // Index` (1-based over all components) names.
        let total = 3 * res.fields.len();
        let pick = match sp.0.get("show norm index") {
            Some(Value::Integer(i)) => i.first().copied().filter(|&k| k > 0).map(|k| k as usize),
            _ => None,
        }
        .unwrap_or(total);
        let mut last = vec![0.0; res.n];
        for (fi, f) in res.fields.iter().enumerate() {
            upsert(&mut ctx.sim.vars, f.name, 3, &perm, f.values.clone(), true)?;
            if pick > 3 * fi && pick <= 3 * fi + 3 {
                let c = pick - 3 * fi - 1;
                for (i, l) in last.iter_mut().enumerate() {
                    *l = f.values[3 * i + c];
                }
            }
        }
        let dummy = upsert(&mut ctx.sim.vars, "hr_dummy", 1, &perm, last, false)?;
        let slot = ctx
            .sim
            .slots
            .iter_mut()
            .find(|s| s.id == id)
            .ok_or_else(|| SolverError::fatal("MagnetoDynamicsCalcFields", "solver slot not found"))?;
        slot.variable = Some(dummy);
        slot.perm = Some(perm);
        ctx.nonlinear_change()?;
        Ok(())
    }
}
