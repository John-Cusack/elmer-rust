//! `Block Solver = True` with a `Procedure` (`BlockSolver`,
//! MainUtils.F90:4323): each coupled iteration runs the solver with its
//! linear solve disabled, scales the assembled monolithic system, solves it
//! by blocks through `elmer-block` (`BlockPickMatrix` + `BlockKrylovIter` /
//! `BlockStandardIter` under `outer:`), scales back, and stops once the
//! block solve's relative change is under `Nonlinear System Convergence
//! Tolerance`. The block-assembled form (no `Procedure`) is `block_system`.
//!
//! `Coupled Solver` / `Assembly Solver` with a `Procedure` (`CoupledSolver`
//! inheriting the matrix, MainUtils.F90:3605) is [`coupled_with_procedure`].

use elmer_core::{Communicator, Logger, MeshView, Params, Value};

use crate::error::SolverError;
use crate::sim::Simulation;
use crate::solver::{Solver, SolverCtx};

const CALLER: &str = "BlockSolver";

/// The monolithic scaling `ScaleLinearSystem(Solver, A, b, x)` applies:
/// the section's method, the right-hand side scaled too.
fn scale_options(p: &elmer_core::ParamMap) -> Result<elmer_sparse::ScaleOptions, SolverError> {
    let method = if p.logical("linear system row equilibration") == Some(true) {
        elmer_sparse::ScalingMethod::RowEquilibration
    } else {
        elmer_sparse::ScalingMethod::from_keyword(p.string("linear system scaling method").unwrap_or("diagonal"))?
    };
    Ok(elmer_sparse::ScaleOptions {
        method,
        apply: true,
        rhs_scaling: true,
    })
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn block_solver(
    sim: &mut Simulation,
    ix: usize,
    inst: &mut Box<dyn Solver>,
    dt: f64,
    transient: bool,
    mesh: &dyn MeshView,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<(), SolverError> {
    log.info(CALLER, "---------------------------------------", 5);
    log.info(CALLER, "Solving system of equations utilizing block strategies", 5);
    log.info(CALLER, "---------------------------------------", 5);
    let p = &sim.slots[ix].params;
    if p.logical("nonlinear system linesearch") == Some(true) {
        return Err(SolverError::Unsupported("Block Solver with Nonlinear System Linesearch".into()));
    }
    let id = sim.slots[ix]
        .variable
        .ok_or_else(|| SolverError::fatal(CALLER, "Block Solver needs the solver's variable"))?;
    let dofs = sim.vars.by_id(id).dofs;
    let iters = p.integer("nonlinear system max iterations").unwrap_or(1);
    let tol = sim.creal(p, "nonlinear system convergence tolerance", mesh)?.unwrap_or(0.0);
    let scale = p.logical("linear system scaling").unwrap_or(true);
    let schur = crate::solver::schur_matrix(sim, &sim.slots[ix].params);

    log.info(CALLER, "-------------------------------------------------", 6);
    for iter in 1..=iters.max(1) {
        // WRITE(Message,'(A,T35,I0)'): the count at column 35.
        log.info(CALLER, &format!("{:<34}{iter}", "Coupled iteration: "), 6);

        // SingleSolver with one nonlinear iteration and no linear solve.
        let params = &mut sim.slots[ix].params;
        params.insert("Nonlinear System Max Iterations", Value::Integer(vec![1]));
        params.insert("Linear System Solver Disabled", Value::Logical(true));
        let r = inst.solve(&mut SolverCtx::new(sim, ix, mesh, log, comm), dt, transient);
        let params = &mut sim.slots[ix].params;
        params.insert("Nonlinear System Max Iterations", Value::Integer(vec![iters]));
        params.0.remove("linear system solver disabled");
        r?;

        let slot = &mut sim.slots[ix];
        let p = &slot.params;
        let sys = slot
            .system
            .as_mut()
            .ok_or_else(|| SolverError::fatal(CALLER, "the solver assembled no system"))?;
        let n = sys.matrix.nrows;
        let var = sim.vars.by_id_mut(id);
        let x = &mut var.values[..n];
        let scaling = if scale {
            log.info(CALLER, "Applying scaling", 8);
            Some(elmer_sparse::scale_system(&mut sys.matrix, Some(&mut sys.rhs[..n]), Some(&mut *x), &scale_options(p)?))
        } else {
            None
        };
        // BlockPickMatrix (CRS_BlockMatrixPick) needs two blocks or more.
        let solved = if dofs < 2 {
            Err(SolverError::fatal("CRS_BlockMatrixPick", "No applicable to just one block!"))
        } else {
            let opts = elmer_block::BlockOptions {
                additive_namespaces: sim.model.simulation.logical("additive namespaces") == Some(true),
                outer_scaling: scaling.as_ref().map(|rs| elmer_block::OuterScaling {
                    diag: rs.diag(),
                    rhs_scaling: rs.rhs_scaling(),
                }),
                global_abort: sim.model.simulation.logical("global abort not converged").unwrap_or(true),
                schur: schur.as_ref(),
            };
            elmer_block::block_solve(p, dofs, &sys.matrix, &sys.rhs[..n], x, &opts, log, comm)
                .map_err(crate::solver::block_error)
        };
        if let Some(rs) = scaling {
            rs.unscale(&mut sys.matrix, Some(&mut sys.rhs[..n]), Some(x));
        }
        let report = solved?;
        var.values_changed = true;
        sim.publish_live();
        if report.max_change < tol {
            break;
        }
    }
    log.info(CALLER, "All done", 8);
    log.info(CALLER, "-------------------------------------------------", 6);
    Ok(())
}

/// `CoupledSolver` for a `Coupled Solver` or `Assembly Solver` section with a
/// `Procedure` (MainUtils.F90:3605, 3820-3912): each coupled iteration runs
/// the procedure with one nonlinear iteration, its linear solve disabled and
/// its change skipped, then solves the assembled system itself
/// (`DefaultSolve`) until the nonlinear change is below `Nonlinear System
/// Convergence Tolerance`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn coupled_with_procedure(
    sim: &mut Simulation,
    ix: usize,
    inst: &mut Box<dyn Solver>,
    dt: f64,
    transient: bool,
    mesh: &dyn MeshView,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<(), SolverError> {
    const CALLER: &str = "CoupledSolver";
    const RULE: &str = "---------------------------------------";
    const LONG_RULE: &str = "-------------------------------------------------";
    log.info(CALLER, RULE, 5);
    log.info(CALLER, "Inheriting matrix from procedure", 5);
    log.info(CALLER, RULE, 5);
    log.info(CALLER, "Using existing variable and matrix", 8);
    if sim.slots[ix].variable.is_none() {
        return Err(SolverError::fatal(CALLER, "In legacy solver mode the Variable should exist!"));
    }
    let p = &sim.slots[ix].params;
    if p.logical("nonlinear system linesearch") == Some(true) {
        return Err(SolverError::Unsupported("Nonlinear System Linesearch in a coupled solver".into()));
    }
    let iters = p.integer("nonlinear system max iterations").unwrap_or(1);
    let tol = sim.creal(p, "nonlinear system convergence tolerance", mesh)?.unwrap_or(0.0);

    log.info(CALLER, LONG_RULE, 5);
    for iter in 1..=iters.max(1) {
        // WRITE(Message,'(A,T35,I5)') 'Coupled iteration:',iter
        log.info(CALLER, &format!("{:<34}{iter:>5}", "Coupled iteration:"), 5);
        let params = &mut sim.slots[ix].params;
        params.insert("Nonlinear System Max Iterations", Value::Integer(vec![1]));
        params.insert("Linear System Solver Disabled", Value::Logical(true));
        params.insert("Skip Compute Nonlinear Change", Value::Logical(true));
        let r = inst.solve(&mut SolverCtx::new(sim, ix, mesh, log, comm), dt, transient);
        let params = &mut sim.slots[ix].params;
        params.insert("Nonlinear System Max Iterations", Value::Integer(vec![iters]));
        params.0.remove("linear system solver disabled");
        params.insert("Skip Compute Nonlinear Change", Value::Logical(false));
        r?;

        let mut ctx = SolverCtx::new(sim, ix, mesh, log, comm);
        ctx.default_solve()?;
        log.info(CALLER, LONG_RULE, 5);
        if ctx.convergence().map_or(0.0, |c| c.nonlin_change) < tol {
            break;
        }
    }
    log.info(CALLER, "All done", 8);
    log.info(CALLER, LONG_RULE, 5);
    Ok(())
}
