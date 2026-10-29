//! Solvers that Elmer runs without a `Procedure` of their own, chosen by the
//! solver mode keywords (MainUtils.F90:1388-1395): `Assembly Solver = True`
//! runs `CoupledSolver` (MainUtils.F90:3565) in its one-equation mode, which
//! assembles the system element by element with the routine named by
//! `Bulk Assembly Procedure` (`BlockSystemAssembly`, MainUtils.F90:4885).
//!
//! `Coupled Solver` (several variables and constraints in one matrix),
//! `Block Solver`, and an assembly solver that also names a `Procedure`
//! (`CoupledSolver`'s "inheriting matrix from procedure" mode) are refused.

use elmer_core::{ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};

use crate::error::SolverError;
use crate::solver::{procedure_value, ElementAssembly, Solver, SolverCtx, SolverRegistry};

const CALLER: &str = "CoupledSolver";
const RULE: &str = "---------------------------------------";
const LONG_RULE: &str = "-------------------------------------------------";

/// The built-in solver a section's mode keywords select, if any. Errors on
/// the modes that are not ported and on an assembly routine the registry
/// does not have.
pub fn solver_for_mode(params: &ParamMap, registry: &SolverRegistry) -> Result<Option<Box<dyn Solver>>, SolverError> {
    // The last of these that is set wins, as in the Fortran.
    let mode = ["coupled solver", "block solver", "assembly solver"]
        .into_iter()
        .rfind(|k| params.logical(k) == Some(true));
    match mode {
        None => Ok(None),
        Some("assembly solver") if params.get("procedure").is_none() => {
            Ok(Some(Box::new(AssemblySolver::from_section(params, registry)?)))
        }
        Some("assembly solver") => Err(SolverError::Unsupported(
            "Assembly Solver with a Procedure (CoupledSolver inheriting the matrix)".into(),
        )),
        Some(m) => Err(SolverError::Unsupported(format!("{m} = True"))),
    }
}

/// `CoupledSolver` for `Assembly Solver = True`: the solver's own variable
/// and matrix, one assembly routine for the (1,1) block.
pub struct AssemblySolver {
    bulk: Box<dyn ElementAssembly>,
    /// `ProcName` as Elmer prints it: library and entry.
    proc_name: String,
}

impl AssemblySolver {
    /// Resolves `Bulk Assembly Procedure 11`, else `Bulk Assembly
    /// Procedure` (BlockSystemAssembly). `Boundary Assembly Procedure` is
    /// not read: `BlockSystemAssembly` switches to boundary mode after the
    /// bulk loop but its `GOTO 100` is commented out (MainUtils.F90:5156),
    /// so Elmer never calls the boundary routine.
    pub fn from_section(params: &ParamMap, registry: &SolverRegistry) -> Result<Self, SolverError> {
        let (lib, entry) = procedure_value(params, "bulk assembly procedure 11")
            .or_else(|| procedure_value(params, "bulk assembly procedure"))
            .ok_or_else(|| SolverError::fatal("BlockSystemAssembly", "Bulk Assembly Procedure not given!"))?;
        let proc_name = format!("{lib} {entry}");
        let f = registry
            .assembly(&lib, &entry)
            .ok_or_else(|| SolverError::fatal("BlockSystemAssembly", format!("Assembly routine not found: {proc_name}")))?;
        Ok(Self { bulk: f(params)?, proc_name })
    }

    /// `BlockSystemAssembly(Solver, dt, Transient, 1, 1)`: every bulk
    /// element whose nodes all have a dof, through the routine, the
    /// solver's time discretisation and the global glue.
    fn assemble(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.info("BlockSystemAssembly", "Starting block system assembly", 8);
        ctx.info("BlockSystemAssembly", &format!("Using assembly routine: {}", self.proc_name), 8);
        // "These may be fetched within the assembly routine, if needed."
        let p = &mut ctx.sim.slots[ctx.slot].params;
        p.insert("Block Matrix Row", Value::Integer(vec![1]));
        p.insert("Block Matrix Column", Value::Integer(vec![1]));

        let dofs = ctx.primary().map_or(1, |v| v.dofs);
        let time_order = ctx.sim.slots[ctx.slot].time_order;
        let mut slots = Vec::new();
        for e in (0..ctx.mesh.n_bulk() as u32).map(ElemId) {
            ctx.element_slots(e, &mut slots);
            if slots.contains(&NO_DOF) {
                continue;
            }
            let n = slots.len() * dofs;
            let mut local = LocalSystem::zeros(n);
            local.mass = Some(vec![0.0; n * n]);
            local.damp = Some(vec![0.0; n * n]);
            self.bulk.assemble(ctx, e, dt, transient, &mut local)?;
            if transient {
                match time_order {
                    1 => ctx.default_first_order_time(&mut local, e, dt)?,
                    2 => ctx.default_second_order_time(&mut local, e, dt)?,
                    _ => {}
                }
            }
            ctx.default_update(&local, e)?;
        }
        Ok(())
    }
}

impl Solver for AssemblySolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, RULE, 5);
        ctx.info(CALLER, "Solving one equation", 5);
        ctx.info(CALLER, RULE, 5);
        ctx.info(CALLER, "Using existing variable and matrix", 8);
        if ctx.primary().is_none() {
            return Err(SolverError::fatal(CALLER, "In legacy solver mode the Variable should exist!"));
        }
        let p = ctx.params();
        if p.logical("eigen analysis") == Some(true) || p.logical("harmonic analysis") == Some(true) {
            return Err(SolverError::Unsupported("eigen/harmonic analysis in an assembly solver".into()));
        }
        if p.logical("nonlinear system linesearch") == Some(true) {
            return Err(SolverError::Unsupported("Nonlinear System Linesearch in an assembly solver".into()));
        }
        let iterations = p.integer("nonlinear system max iterations").unwrap_or(1);
        let tol = ctx.creal(ctx.params(), "nonlinear system convergence tolerance")?.unwrap_or(0.0);

        ctx.info(CALLER, LONG_RULE, 5);
        for iter in 1..=iterations {
            // WRITE(Message,'(A,T35,I5)') 'Coupled iteration:',iter
            ctx.info(CALLER, &format!("{:<34}{:>5}", "Coupled iteration:", iter), 5);
            ctx.default_initialize()?;
            self.assemble(ctx, dt, transient)?;
            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            ctx.info(CALLER, LONG_RULE, 5);
            let change = ctx.convergence().map_or(0.0, |c| c.nonlin_change);
            if change < tol {
                break;
            }
        }
        ctx.info(CALLER, "All done", 8);
        ctx.info(CALLER, LONG_RULE, 5);
        Ok(())
    }
}
