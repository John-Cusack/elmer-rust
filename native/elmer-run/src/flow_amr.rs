//! FlowSolve with `Adaptive Mesh Refinement` (FlowSolve.F90:1611-1616):
//! after the flow solve, `RefineMesh` with FlowSolve's residual estimators
//! ([`elmer_adaptive::flow`]). A wrapper here, as `elmer-physics-flow`
//! does not depend on `elmer-adaptive`.

use elmer_adaptive::flow::FlowEstimator;
use elmer_adaptive::RefineState;
use elmer_core::{ParamMap, Params};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

struct FlowAmr {
    inner: Box<dyn Solver>,
    refine: RefineState,
}

impl Solver for FlowAmr {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        self.inner.init(ctx)
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        self.inner.solve(ctx, dt, transient)?;
        let p = ctx.params();
        if p.logical("adaptive mesh refinement") == Some(true) && p.logical("library adaptivity") != Some(true) {
            FlowEstimator::check(ctx)?;
            elmer_adaptive::refine_mesh(ctx, &mut self.refine, &mut FlowEstimator)?;
        }
        Ok(())
    }
}

/// `FlowSolver`, wrapped when the section refines the mesh.
fn factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    let inner = elmer_physics_flow::factory(params)?;
    if params.logical("adaptive mesh refinement") == Some(true) {
        return Ok(Box::new(FlowAmr { inner, refine: RefineState::default() }));
    }
    Ok(inner)
}

/// Registers the wrapping factory over `FlowSolve` `FlowSolver`.
pub fn register(reg: &mut SolverRegistry) {
    let (lib, entry) = elmer_physics_flow::REGISTRY_ENTRY;
    reg.register(lib, entry, factory);
}
