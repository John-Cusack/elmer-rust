//! Which solver's variable a `Linear System Solver = Block` solve starts from
//! and writes into. This reproduces how BlockSolve.F90 shares its block
//! matrix between solvers.
//!
//! `TotMatrix` is one module-level pointer (`TYPE(BlockMatrix_t), POINTER,
//! SAVE :: TotMatrix`, BlockSolve.F90:57). `BlockSolveInt` re-points it, via
//! `BlockInitMatrix`, only when the solver has no block matrix yet
//! (`OldMatrix = ASSOCIATED(Solver % BlockMatrix)`, :5049; the call at
//! :5156-5157; :312-329). With the default split of a multi-dof variable,
//! the block vectors `TotMatrix % SubVector(i) % Var` are the components of
//! the variable of the solver that set the pointer (`ComponentName`, :480).
//! `BlockKrylovIter` starts from those vectors (:4527) and writes its
//! solution back into them (:4640-4643); this path has no other back-copy
//! (`BlockBackCopyVar` only runs for `DoMyOwnVars`, :5378).
//!
//! So when two solvers both block-solve, the first one to come back to its
//! block solve (`OldMatrix`) finds `TotMatrix` still pointing at the other
//! solver's block matrix. It starts from the other solver's variable, its
//! solution goes into that variable, and its own variable keeps what it
//! had. BlockLinElast2bMany depends on this in its second steady iteration.
//!
//! Only that verified path is reproduced: serial, the block preconditioning
//! strategy (`BlockKrylovIter`) on the plain component split, with both
//! variables of the same size. Any other block strategy keeps its own
//! variable.

use std::cell::RefCell;

use elmer_core::{normalize_key, ParamMap, Params};

use crate::sim::Simulation;
use crate::variable::VarId;

/// The solver `TotMatrix` points at.
#[derive(Clone, Copy, Debug)]
pub(crate) struct BlockOwner {
    pub slot: usize,
    /// Its block vectors are the components of its variable.
    pub components: bool,
}

/// A block solve that is about to run.
pub(crate) struct BlockPlan {
    slot: usize,
    components: bool,
    /// The owner's variable and a copy of its values, when this solve goes
    /// through another solver's block matrix.
    alias: Option<(VarId, RefCell<Vec<f64>>)>,
}

impl BlockPlan {
    /// Not the plain solve of the solver's matrix after all: the solve still
    /// sets up a block matrix, but neither aliases nor gives `TotMatrix`
    /// component vectors.
    pub fn demote(&mut self) {
        self.components = false;
        self.alias = None;
    }

    /// The owner's values for [`crate::solver::BlockEnv::alias`].
    pub fn alias(&self) -> Option<&RefCell<Vec<f64>>> {
        self.alias.as_ref().map(|(_, v)| v)
    }
}

pub(crate) fn is_block(p: &ParamMap) -> bool {
    p.string("linear system solver")
        .is_some_and(|s| normalize_key(s) == "block")
}

/// `BlockSolveInt` uses the components of the solver's variable as block
/// vectors and runs `BlockKrylovIter`: no `DoMyOwnVars` system (:5083-5092),
/// no `Block Solvers` (:5112), no `Block Structure` or `Variable i`
/// names (:473-481), not nested or monolithic, block preconditioning
/// on (:5057, :5345-5350).
fn component_krylov(p: &ParamMap, dofs: usize) -> bool {
    const OWN_VARS: [&str; 11] = [
        "block domain system",
        "block hdiv system",
        "block re-im system",
        "block hcurl system",
        "block quadratic hcurl system",
        "block a-v system",
        "block hor-ver system",
        "block cartesian system",
        "block nodal system",
        "block nested system",
        "block monolithic",
    ];
    dofs > 1
        && OWN_VARS.iter().all(|k| p.logical(k) != Some(true))
        && p.get("block solvers").is_none()
        && p.get("block structure").is_none()
        && (1..=dofs).all(|i| p.get(&format!("variable {i}")).is_none())
        && p.logical("block preconditioner") != Some(false)
}

/// The plan for a block solve of `slot`, or `None` when it is not a block
/// solve. `main` says whether it is the plain `SolveLinearSystem` of the
/// solver's own matrix (not EliminateDirichlet's reduced system, not a
/// constraint collection), the only one that may alias.
pub(crate) fn plan(sim: &Simulation, slot: usize, main: bool, parallel: bool) -> Option<BlockPlan> {
    let s = &sim.slots[slot];
    if !is_block(&s.params) {
        return None;
    }
    let id = s.variable?;
    let var = sim.vars.by_id(id);
    let components = main && component_krylov(&s.params, var.dofs);
    let alias = match sim.block_tot {
        Some(o) if s.block_matrix && o.slot != slot && o.components && components && !parallel => {
            sim.slots[o.slot].variable.and_then(|oid| {
                let ov = sim.vars.by_id(oid);
                (ov.dofs == var.dofs && ov.values.len() == var.values.len())
                    .then(|| (oid, RefCell::new(ov.values.clone())))
            })
        }
        _ => None,
    };
    Some(BlockPlan {
        slot,
        components,
        alias,
    })
}

impl Simulation {
    /// After the block solve of `plan` ran: a first block solve points
    /// `TotMatrix` at its solver; an aliased one hands its solution to the
    /// owner's variable.
    pub(crate) fn finish_block_solve(&mut self, plan: BlockPlan) {
        if !self.slots[plan.slot].block_matrix {
            self.slots[plan.slot].block_matrix = true;
            self.block_tot = Some(BlockOwner {
                slot: plan.slot,
                components: plan.components,
            });
        }
        if let Some((oid, vals)) = plan.alias {
            let ov = self.vars.by_id_mut(oid);
            ov.values = vals.into_inner();
            ov.values_changed = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Model;
    use crate::sim::SolverSlot;
    use crate::variable::Variable;
    use elmer_core::FieldKind;
    use elmer_core::Value;

    fn block_sim() -> Simulation {
        let mut sim = Simulation::new(Model::default());
        for (k, name) in ["dispa", "dispb"].into_iter().enumerate() {
            let mut p = ParamMap::new();
            p.insert("Linear System Solver", Value::Str("block".into()));
            let mut slot = SolverSlot::new(k + 1, p).unwrap();
            let v = Variable::new(name, FieldKind::Nodal, 2, None, vec![k as f64 + 1.0; 4]);
            slot.variable = Some(sim.vars.add(v).unwrap());
            sim.slots.push(slot);
        }
        sim
    }

    /// One block solve of `slot` whose solution is `x`, as `solve_system`
    /// runs it: the solve starts from the alias when there is one and hands
    /// the solution to it, else it updates the solver's own variable.
    fn solve(sim: &mut Simulation, slot: usize, x: f64) -> Option<Vec<f64>> {
        let plan = plan(sim, slot, true, false).unwrap();
        let start = plan.alias().map(|a| a.borrow().clone());
        match plan.alias() {
            Some(a) => a.borrow_mut().fill(x),
            None => {
                let id = sim.slots[slot].variable.unwrap();
                sim.vars.by_id_mut(id).values.fill(x);
            }
        }
        sim.finish_block_solve(plan);
        start
    }

    fn values(sim: &Simulation, name: &str) -> Vec<f64> {
        sim.vars.get(name).unwrap().values.clone()
    }

    #[test]
    fn second_visit_goes_through_the_last_initialised_block_matrix() {
        let mut sim = block_sim();
        // Steady iteration 1: each solver sets up its own block matrix.
        assert_eq!(solve(&mut sim, 0, 10.0), None);
        assert_eq!(solve(&mut sim, 1, 20.0), None);
        // Iteration 2: TotMatrix still points at solver 2's block matrix.
        // Solver 1 starts from DispB, leaves its solution there and keeps
        // its own values; solver 2 then starts from that solution.
        assert_eq!(solve(&mut sim, 0, 11.0), Some(vec![20.0; 4]));
        assert_eq!(
            (values(&sim, "dispa"), values(&sim, "dispb")),
            (vec![10.0; 4], vec![11.0; 4])
        );
        assert_eq!(solve(&mut sim, 1, 21.0), None);
        assert_eq!(values(&sim, "dispb"), vec![21.0; 4]);
    }

    #[test]
    fn other_strategies_and_sizes_do_not_alias() {
        let mut sim = block_sim();
        sim.slots[0]
            .params
            .insert("Block Monolithic", Value::Logical(true));
        solve(&mut sim, 0, 10.0);
        solve(&mut sim, 1, 20.0);
        assert_eq!(solve(&mut sim, 0, 11.0), None);
        assert_eq!(values(&sim, "dispa"), vec![11.0; 4]);

        let mut sim = block_sim();
        solve(&mut sim, 0, 10.0);
        solve(&mut sim, 1, 20.0);
        sim.vars.get_mut("dispb").unwrap().values.push(0.0);
        assert_eq!(solve(&mut sim, 0, 11.0), None);

        // A solver that is not a block solve has no plan.
        let mut sim = block_sim();
        sim.slots[0]
            .params
            .insert("Linear System Solver", Value::Str("direct".into()));
        assert!(plan(&sim, 0, true, false).is_none());
    }
}
