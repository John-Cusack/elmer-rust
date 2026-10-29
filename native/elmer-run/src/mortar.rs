//! `Apply Mortar BCs`: the solver's Galerkin mortar constraints.
//!
//! Elmer's `DefaultSolve` builds the boundary projectors
//! (`GenerateProjectors`) and the constraint matrix
//! (`GenerateConstraintMatrix`) and solves the system restricted by it.
//! `elmer-mortar` builds the two; `elmer-solver` solves with a solver's own
//! constraint rows (`SolverCtx::set_constraint_matrix`). [`MortarSolver`]
//! wraps the slot's solver and hands the rows over before its first solve.

use std::path::PathBuf;

use elmer_core::{CsrMatrix, NullLogger, ParamMap, Params};
use elmer_mortar::{generate_constraint_matrix, generate_mortar_projectors, ConstraintTarget, ProjectorEnv};
use elmer_solver::{ConstraintMatrix, Solver, SolverCtx, SolverError};

/// Whether the section asks for mortar constraints.
pub(crate) fn applies(params: &ParamMap) -> bool {
    params.logical("apply mortar bcs") == Some(true)
}

/// `inner` with its mortar constraint matrix installed before the first
/// solve. `inner` was built from the registry after setup, so its `init`
/// runs then too.
pub(crate) struct MortarSolver {
    pub inner: Box<dyn Solver>,
    pub coord_dim: usize,
    /// Where `Save Projector` / `Save Constraint Matrix` files go.
    pub case_dir: PathBuf,
    pub ready: bool,
}

impl MortarSolver {
    fn constraint_matrix(&self, ctx: &SolverCtx<'_>) -> Result<Option<ConstraintMatrix>, SolverError> {
        let model = &ctx.sim.model;
        let bcs: Vec<&dyn Params> = model.boundary_conditions.iter().map(|b| b as &dyn Params).collect();
        let bc_of = |tag: i32| model.bc_of_tag(tag);
        let env = ProjectorEnv {
            mesh: ctx.mesh,
            bcs: &bcs,
            bc_of: &bc_of,
            simulation: &model.simulation,
            solver: Some(ctx.params()),
            coord_dim: self.coord_dim,
            log: &NullLogger,
        };
        let Some((projectors, stop)) = generate_mortar_projectors(&env, &[], false, &self.case_dir, 0)? else {
            return Ok(None);
        };
        if stop {
            return Err(SolverError::Unsupported("Save Projector And Stop".into()));
        }
        let name = ctx.variable_name().ok_or_else(|| SolverError::fatal("GenerateConstraintMatrix", "no variable"))?;
        let var = ctx.variable(name).ok_or_else(|| SolverError::fatal("GenerateConstraintMatrix", "no variable"))?;
        let perm = var.perm.ok_or_else(|| SolverError::Unsupported("mortar BCs on a variable without Perm".into()))?;
        let target = ConstraintTarget {
            solver: ctx.params(),
            name,
            dofs: var.dofs,
            perm,
            n_rows: var.values.len(),
        };
        let Some(c) = generate_constraint_matrix(&bcs, &projectors, &target, &NullLogger)? else {
            return Ok(None);
        };
        let n = c.n_rows();
        let ncols = target.n_rows + n;
        Ok(Some(ConstraintMatrix {
            matrix: CsrMatrix { nrows: n, ncols, row_ptr: c.rows, col_idx: c.cols, values: c.values, diag: None },
            rhs: c.rhs,
            inv_perm: c.inv_perm.into_iter().map(|i| i as u32).collect(),
        }))
    }
}

impl Solver for MortarSolver {
    fn init(&mut self, _ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        if !self.ready {
            self.inner.init(ctx)?;
            let c = self.constraint_matrix(ctx)?;
            ctx.set_constraint_matrix(c);
            self.ready = true;
        }
        self.inner.solve(ctx, dt, transient)
    }
}
