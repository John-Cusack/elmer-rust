//! Keyword helpers shared by the mappers (Lists.F90 semantics).

use elmer_core::{ParamMap, Params, Value};
use elmer_solver::expr::eval_at_nodes_with;
use elmer_solver::{SolverCtx, SolverError};

use crate::motion;

/// `ListGetReal(section, key, 1, [node])`: the keyword at one node, `None`
/// when the section lacks it.
pub fn node_real(
    ctx: &SolverCtx<'_>,
    section: &ParamMap,
    key: &str,
    node: u32,
) -> Result<Option<f64>, SolverError> {
    let Some(v) = section.get(key) else {
        return Ok(None);
    };
    let b = &ctx.sim.backends;
    let xs = eval_at_nodes_with(v, &[node], &ctx.fields(), b.matc.as_deref(), &b.udfs)?;
    Ok(xs.first().copied())
}

/// The first column of a constant real array (`ListGetConstRealArray(...)(:,1)`).
pub fn const_array(section: &ParamMap, key: &str) -> Option<Vec<f64>> {
    match section.get(key)? {
        Value::Real { rows, cols, data } => Some((0..*rows).map(|i| data[i * cols]).collect()),
        _ => None,
    }
}

/// `ListCheckPrefix`: any key of the section starts with `prefix`
/// (normalised, lower case).
pub fn has_prefix(section: &ParamMap, prefix: &str) -> bool {
    section.0.keys().any(|k| k.starts_with(prefix))
}

/// Moves the nodes through the installed [`motion::MeshMotion`]; a no-op when
/// nothing changes.
pub fn apply_coords(caller: &str, old: &[[f64; 3]], new: &[[f64; 3]]) -> Result<(), SolverError> {
    if old == new {
        return Ok(());
    }
    match motion::current() {
        Some(m) => {
            m.set_coords(new);
            Ok(())
        }
        None => Err(SolverError::NotWired(format!(
            "{caller}: moving mesh nodes (no elmer_physics_meshtools::motion::install on this thread)"
        ))),
    }
}

pub fn coords_of(mesh: &dyn elmer_core::MeshView) -> Vec<[f64; 3]> {
    (0..mesh.n_nodes() as u32).map(|n| mesh.coord(n)).collect()
}

pub fn slot_index(ctx: &SolverCtx<'_>) -> usize {
    let id = ctx.solver_id();
    ctx.sim
        .slots
        .iter()
        .position(|s| s.id == id)
        .expect("running solver has a slot")
}
