//! [`LiveModel`]: the part of Elmer's `Model_t` a `Procedure` function
//! reads on every call (`Model % Nodes`, `Model % NumberOfNodes`, `Time`,
//! `VariableGet` + `Perm`), as a shared object that follows the run.
//!
//! A user function is an [`elmer_expr::NodeFn`]: it gets a node and its
//! dependency values, not the simulation, so one that reads other fields
//! keeps an `Arc` to this view. The [`crate::Simulation`] publishes its
//! fields into it at the points where a Fortran function would see new
//! values: setup (before each solver's initial conditions), every solver
//! visit, every `DefaultInitialize` and after every `DefaultSolve`
//! ([`crate::SolverCtx::sync_model`] does it on demand). Nothing is
//! published until [`crate::Simulation::live_model`] creates the view.

use std::sync::{PoisonError, RwLock};

use elmer_core::MeshView;
use elmer_expr::FieldLookup;

use crate::variable::VariableStore;

/// The published state: node coordinates and the current values of every
/// variable (histories are not copied).
pub struct LiveModel {
    coords: Vec<[f64; 3]>,
    vars: RwLock<VariableStore>,
}

impl LiveModel {
    pub(crate) fn new(mesh: &dyn MeshView, vars: &VariableStore) -> Self {
        Self {
            coords: (0..mesh.n_nodes() as u32).map(|n| mesh.coord(n)).collect(),
            vars: RwLock::new(vars.values_snapshot()),
        }
    }

    pub(crate) fn publish(&self, vars: &VariableStore) {
        *self.vars.write().unwrap_or_else(PoisonError::into_inner) = vars.values_snapshot();
    }

    /// `Model % NumberOfNodes`.
    pub fn n_nodes(&self) -> usize {
        self.coords.len()
    }
}

impl FieldLookup for LiveModel {
    fn value(&self, var: &str, node: u32) -> Option<f64> {
        crate::expr::lookup(&self.vars.read().unwrap_or_else(PoisonError::into_inner), var, node)
    }

    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }

    fn time(&self) -> f64 {
        self.vars
            .read()
            .unwrap_or_else(PoisonError::into_inner)
            .scalar("time")
            .unwrap_or(0.0)
    }
}
