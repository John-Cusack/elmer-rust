//! Which provider serves which seam of `elmer-solver`. Assembly, time and the
//! linear solvers are `elmer-solver`'s default backends; this crate adds the
//! physics procedures.

use elmer_solver::SolverRegistry;

/// The procedures this build can run.
pub fn solver_registry() -> SolverRegistry {
    let mut reg = SolverRegistry::new();
    elmer_physics_heat::register(&mut reg);
    reg
}
