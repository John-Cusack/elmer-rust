//! `Apply Mortar BCs`: the solver's Galerkin mortar constraints.
//!
//! Elmer's `DefaultSolve` builds the boundary projectors
//! (`GenerateProjectors`) and the constraint matrix
//! (`GenerateConstraintMatrix`) and solves the system restricted by it.
//! `elmer-mortar` builds the two; `elmer-solver` solves with a solver's own
//! constraint rows (`SolverCtx::set_constraint_matrix`). [`MortarSolver`]
//! wraps the slot's solver and hands the rows over before its first solve.

use std::path::PathBuf;
use std::sync::Arc;

use elmer_core::{CsrMatrix, NullLogger, ParamMap, Params};
use elmer_mortar::{
    generate_constraint_matrix, generate_mortar_projectors, periodic_active, periodic_projector, periodic_scale,
    ConstraintTarget, Projector, ProjectorEnv,
};
use elmer_solver::{ConstraintMatrix, Solver, SolverCtx, SolverError};

/// Whether the section asks for mortar constraints.
pub(crate) fn applies(params: &ParamMap) -> bool {
    params.logical("apply mortar bcs") == Some(true)
}

/// Whether a boundary condition constrains the section's variable with
/// `Periodic BC Use Lagrange Coefficient`, which `SetPeriodicBoundariesPass1`
/// turns into a mortar constraint of the solver.
pub(crate) fn lagrange_periodic(params: &ParamMap, bcs: &[ParamMap]) -> bool {
    let Some(spec) = params.string("variable") else { return false };
    let Ok(v) = elmer_solver::sim::parse_variable_spec(spec, params.integer("variable dofs")) else {
        return false;
    };
    bcs.iter().any(|bc| is_lagrange(bc, &v.name))
}

fn is_lagrange(bc: &dyn Params, name: &str) -> bool {
    bc.logical("periodic bc use lagrange coefficient") == Some(true) && periodic_active(bc, name)
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
    /// `Mesh % DiscontPerm` of a split discontinuous mesh.
    pub discont_perm: Option<Arc<Vec<u32>>>,
    /// A mortar jump was set: the constraints are rebuilt before each solve
    /// (`MortarBCsChanged`).
    pub jump: std::cell::Cell<bool>,
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
            discont_perm: self.discont_perm.as_deref().map(Vec::as_slice),
        };
        let name = ctx.variable_name().ok_or_else(|| SolverError::fatal("GenerateConstraintMatrix", "no variable"))?;
        let mut projectors: Vec<Option<Projector>> = match generate_mortar_projectors(&env, &[], false, &self.case_dir, 0)? {
            Some((p, false)) => p,
            Some((_, true)) => return Err(SolverError::Unsupported("Save Projector And Stop".into())),
            None => vec![None; bcs.len()],
        };
        // SetPeriodicBoundariesPass1: a Lagrange periodic boundary becomes a
        // constraint, SlaveScale = -Scale and MasterScale = -1.
        for (i, bc) in bcs.iter().enumerate() {
            if !is_lagrange(*bc, name) {
                continue;
            }
            if bc.get(&format!("periodic bc coefficient {name}")).is_some() {
                return Err(SolverError::Unsupported(
                    "Periodic BC Coefficient with Periodic BC Use Lagrange Coefficient (MortarBC Diag)".into(),
                ));
            }
            if projectors[i].is_some() {
                return Err(SolverError::Unsupported(format!(
                    "Boundary {}: both a mortar and a Lagrange periodic projector",
                    i + 1
                )));
            }
            let target = bc.integer("periodic bc").unwrap_or(0).max(0) as usize;
            let Some(mut p) = periodic_projector(&env, i + 1, target, None)? else { continue };
            let scale = periodic_scale(*bc, name).unwrap_or(-1.0);
            p.slave_scale = -scale;
            p.master_scale = -1.0;
            p.mortar_perm = true;
            projectors[i] = Some(p);
        }
        // SetWeightedProjectorJump: `Mortar BC Offset` at the nodes of the
        // boundary's elements, in MortarBC % Rhs by projector row. The
        // Fortran sets it with the Dirichlet conditions, so it is redone
        // before every solve.
        let mut jump = false;
        let key = format!("mortar bc offset {name}");
        for (i, p) in projectors.iter_mut().enumerate() {
            let Some(p) = p else { continue };
            let bc = &model.boundary_conditions[i];
            if bc.get(&key).is_none() {
                continue;
            }
            jump = true;
            let var = ctx.variable(name).ok_or_else(|| SolverError::fatal("SetWeightedProjectorJump", "no variable"))?;
            let perm = var.perm.ok_or_else(|| SolverError::Unsupported("mortar BCs on a variable without Perm".into()))?;
            let mut row_of = std::collections::HashMap::new();
            for (r, &k) in p.inv_perm.iter().enumerate() {
                if k != elmer_core::NO_DOF {
                    row_of.insert(k, r);
                }
            }
            let target_bc = bc.integer("mortar bc").unwrap_or(0);
            let mut rhs = vec![0.0; p.n_rows()];
            let mut done = vec![false; p.n_rows()];
            let mesh = ctx.mesh;
            for t in 0..mesh.n_boundary() {
                let e = elmer_core::ElemId((mesh.n_bulk() + t) as u32);
                let el = mesh.element(e);
                if model.bc_of_tag(el.tag) != Some(i + 1) {
                    continue;
                }
                let [left, right] = el.parents;
                if if target_bc > 0 { left.is_none() } else { left.is_none() || right.is_none() } {
                    continue;
                }
                let Some(voff) = ctx.real_at_nodes(bc, &key, e)? else { continue };
                for (u, &node) in el.nodes.iter().enumerate() {
                    if perm[node as usize] == elmer_core::NO_DOF {
                        continue;
                    }
                    let Some(&r) = row_of.get(&node) else { continue };
                    if done[r] {
                        continue;
                    }
                    done[r] = true;
                    rhs[r] = voff[u];
                }
            }
            p.mortar_perm = true;
            p.rhs = Some(rhs);
        }
        self.jump.set(jump);
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
        }
        if !self.ready || self.jump.get() {
            let c = self.constraint_matrix(ctx)?;
            ctx.set_constraint_matrix(c);
            self.ready = true;
        }
        self.inner.solve(ctx, dt, transient)
    }
}
