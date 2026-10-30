//! `SlaveSolverSlots/CountSolver.F90`: `SUBROUTINE CountSolver`, which counts
//! its visits (`INTEGER, SAVE :: VisitedTimes`) and stores the count as its
//! variable's value and norm. The cases call it only as a slave solver
//! (`Pre Solvers`, `Post Solvers` and their nonlinear forms).

use elmer_core::ParamMap;
use elmer_solver::{Solver, SolverCtx, SolverError};

#[derive(Default)]
pub struct CountSolver {
    visited: u32,
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`CountSolver`].
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(CountSolver::default()))
}

impl Solver for CountSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        self.visited += 1;
        ctx.info("CountSolver", &format!("Solver visited times: {}", self.visited), 4);
        let count = f64::from(self.visited);
        let sid = ctx.solver_id();
        let id = ctx
            .sim
            .slots
            .iter()
            .find(|s| s.id == sid)
            .and_then(|s| s.variable)
            .ok_or_else(|| SolverError::fatal("CountSolver", "solver has no variable"))?;
        let v = ctx.sim.vars.by_id_mut(id);
        v.values.iter_mut().for_each(|x| *x = count);
        v.conv.norm = count;
        Ok(())
    }
}
