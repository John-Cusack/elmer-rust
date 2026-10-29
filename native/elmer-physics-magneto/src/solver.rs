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

use elmer_core::{ElemId, MeshView, ParamMap, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::case::Case;
use crate::error::MagnetoError;

/// The steady solver's `(library, entry)`.
pub const STEADY: (&str, &str) = ("MagnetoDynamics", "WhitneyAVSolver");
/// The time-harmonic solver's `(library, entry)`.
pub const HARMONIC: (&str, &str) = ("MagnetoDynamics", "WhitneyAVHarmonicSolver");

/// Registers both procedures.
pub fn register(registry: &mut SolverRegistry) {
    registry.register(STEADY.0, STEADY.1, steady_factory);
    registry.register(HARMONIC.0, HARMONIC.1, harmonic_factory);
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
