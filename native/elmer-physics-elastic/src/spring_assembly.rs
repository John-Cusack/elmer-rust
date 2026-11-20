//! `SpringAssembler` (`fem/src/modules/SpringAssembly.F90`): point springs
//! at `Target Nodes`, added to another solver's matrix. It runs as one of
//! that solver's `Assembly Solvers`, at the end of its `DefaultInitialize`.
//!
//! For each boundary condition with `Target Nodes`, `Spring i` is the
//! stiffness of displacement component `i` at each node (a one-node element
//! with a diagonal matrix, glued by `DefaultUpdateEquations`). `Mass i` is
//! only read in transient and harmonic runs, which are refused here.

use elmer_core::{normalize_key, ParamMap, Params, Value};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "SpringAssembler";

pub(crate) fn factory(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(SpringAssembler))
}

struct SpringAssembler;

/// A keyword given as one constant real.
fn constant(v: &Value, key: &str) -> Result<f64, SolverError> {
    match v {
        Value::Real { data, .. } if data.len() == 1 => Ok(data[0]),
        _ => Err(SolverError::Unsupported(format!("{CALLER}: '{key}' that is not a constant"))),
    }
}

impl Solver for SpringAssembler {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        let uname = normalize_key(ctx.params().string("displacement variable name").unwrap_or("displacement"));
        let master = ctx
            .sim
            .slots
            .iter()
            .position(|s| s.variable.is_some_and(|id| normalize_key(&ctx.sim.vars.by_id(id).name) == uname))
            .ok_or_else(|| {
                SolverError::fatal(CALLER, format!("Solver associated with Displacement Variable Name > {uname} < not found!"))
            })?;
        let mp = &ctx.sim.slots[master].params;
        let harmonic = ctx.sim.slots[master].nof_eigen > 0
            || mp.logical("harmonic mode") == Some(true)
            || mp.logical("harmonic analysis") == Some(true);
        let id = ctx.sim.slots[master].variable.expect("found by its variable");
        let (dofs, perm) = {
            let v = ctx.sim.vars.by_id(id);
            (v.dofs, v.perm.clone())
        };
        let perm = perm.ok_or_else(|| SolverError::fatal(CALLER, "displacement has no permutation"))?;

        // (dof, spring) in the order DefaultUpdateEquations adds them.
        let mut adds: Vec<(usize, f64)> = Vec::new();
        for bc in &ctx.sim.model.boundary_conditions {
            let Some(Value::Integer(nodes)) = bc.get("target nodes") else { continue };
            if bc.get("spring").is_some() {
                ctx.warn(CALLER, "Define a spring by components (Spring i = ...)");
                ctx.warn(CALLER, "Skipping a definition (Spring = ...)");
            }
            let mut springs: Vec<Option<f64>> = Vec::with_capacity(dofs);
            for i in 1..=dofs {
                let key = format!("spring {i}");
                springs.push(bc.get(&key).map(|v| constant(v, &key)).transpose()?);
            }
            let any_mass = (1..=dofs).any(|i| bc.get(&format!("mass {i}")).is_some());
            if any_mass && (transient || harmonic) {
                return Err(SolverError::Unsupported(format!("{CALLER}: point masses in a transient or harmonic run")));
            }
            if springs.iter().all(Option::is_none) {
                continue;
            }
            for &node in nodes {
                let Some(&s) = usize::try_from(node - 1).ok().and_then(|k| perm.get(k)) else { continue };
                if s == elmer_core::NO_DOF {
                    continue;
                }
                for (j, k) in springs.iter().enumerate() {
                    // K(j,j) = Spring j where given, 0 otherwise.
                    adds.push((dofs * s as usize + j, k.unwrap_or(0.0)));
                }
            }
        }
        let sys = ctx.sim.slots[master]
            .system
            .as_mut()
            .ok_or_else(|| SolverError::fatal(CALLER, "the displacement solver has no matrix yet"))?;
        for (d, k) in adds {
            let m = &mut sys.matrix;
            let (lo, hi) = (m.row_ptr[d], m.row_ptr[d + 1]);
            let pos = (lo..hi)
                .find(|&p| m.col_idx[p] as usize == d)
                .ok_or_else(|| SolverError::fatal(CALLER, format!("no diagonal entry for dof {d}")))?;
            m.values[pos] += k;
        }
        Ok(())
    }
}
