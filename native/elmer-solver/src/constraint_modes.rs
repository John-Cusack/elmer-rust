//! `Constraint Modes Analysis` on boundary conditions: every dof of a
//! boundary with `Constraint Modes <var> = True` (or `Constraint Modes`) is
//! a mode (`SetConstraintModesBoundaries`, BoundaryConditionUtils.F90:2765):
//! all of them are held at zero, and the linear system is solved once per
//! mode with that dof at one (`ConstraintModesDriver` around
//! `SolveLinearSystem`, SolveCore.F90:3215-3290, 3995-4290), each solve with
//! its own nonlinear change. `Constraint Mode <var> = Integer k` numbers a
//! boundary (0 is the ground, held at zero).
//!
//! Refused: lumped modes, the `rhs`, `EM Wave` and `Coils` variants,
//! `Constraint Modes Single` and mode directions, complex systems, modes over
//! body forces or components, the external loops (`Nonlinear System`,
//! `Steady State` and `Run Control Constraint Modes`), the lumped fluxes and
//! linkage, symmetric Dirichlet elimination, p elements and parallel runs.

use std::collections::BTreeMap;

use elmer_core::{Params, Value};

use crate::error::SolverError;
use crate::fortran_fmt::i2s;
use crate::solver::SolverCtx;

const CALLER: &str = "SetConstraintModesBoundaries";

/// The mode of each matrix row: `k > 0` the `k`-th mode, `< 0` the ground,
/// `0` none.
#[derive(Clone, Debug, Default)]
pub(crate) struct Modes {
    pub index: Vec<i32>,
    pub count: usize,
}

fn refuse(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    let p = ctx.params();
    for k in [
        "constraint modes lumped",
        "constraint modes rhs",
        "constraint modes em wave",
        "constraint modes coils",
        "constraint modes single",
        "nonlinear system constraint modes",
        "steady state constraint modes",
        "run control constraint modes",
        "constraint modes fluxes",
        "constraint modes linkage",
        "constraint modes fluxes norm",
        "constraint modes matrix norm",
        "constraint modes constant matrix",
        "linear system complex",
        "linear system symmetric",
        "eigen analysis",
        "solve adjoint equation",
        "apply source control",
        "monolithic slave",
    ] {
        if p.logical(k) == Some(true) {
            return Err(SolverError::Unsupported(format!("Constraint Modes Analysis with {k}")));
        }
    }
    let m = &ctx.sim.model;
    if m.body_forces.iter().any(|b| b.0.keys().any(|k| k.starts_with("constraint mode"))) {
        return Err(SolverError::Unsupported("constraint modes over body forces".into()));
    }
    if m.boundary_conditions.iter().any(|b| b.get("constraint mode direction").is_some()) {
        return Err(SolverError::Unsupported("Constraint Mode Direction".into()));
    }
    if p.get("nonlinear post solvers").is_some() {
        return Err(SolverError::Unsupported("Constraint Modes Analysis with Nonlinear Post Solvers".into()));
    }
    if ctx.comm.size() > 1 || ctx.solver_dofs().is_some() {
        return Err(SolverError::Unsupported("constraint modes in parallel or on p elements".into()));
    }
    Ok(())
}

/// `SetConstraintModesBoundaries` (from `DefaultDirichletBCs`,
/// DefUtils.F90:7125): numbers the mode dofs and holds them (and the ground)
/// at zero in `fixed`, over any other condition.
pub(crate) fn set_boundaries(ctx: &mut SolverCtx<'_>, fixed: &mut BTreeMap<u32, f64>) -> Result<(), SolverError> {
    if ctx.params().logical("constraint modes analysis") != Some(true) {
        return Ok(());
    }
    let model = &ctx.sim.model;
    if !model.boundary_conditions.iter().any(|b| b.0.keys().any(|k| k.starts_with("constraint mode"))) {
        // NumberOfConstraintModes = 0: nothing to do.
        ctx.sim.slots[ctx.slot].cmodes = None;
        return Ok(());
    }
    refuse(ctx)?;
    let id = ctx.sim.slots[ctx.slot]
        .variable
        .ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    let (name, dofs, perm, nrows) = {
        let v = ctx.sim.vars.by_id(id);
        (v.name.clone(), v.dofs, v.perm.clone(), v.values.len())
    };
    let perm = perm.ok_or_else(|| SolverError::fatal(CALLER, "solver variable has no permutation"))?;
    ctx.info(CALLER, &format!("Setting constraint modes boundaries for variable: {name}"), 7);

    // BCPerm: `Constraint Mode <var>` numbers a boundary (0 the ground),
    // `Constraint Modes <var>` the next one.
    let model = &ctx.sim.model;
    let mut bc_perm = vec![0i32; model.boundary_conditions.len()];
    let mut j = 0;
    for (i, bc) in model.boundary_conditions.iter().enumerate() {
        let k = bc.integer("constraint mode").or_else(|| bc.integer(&format!("constraint mode {name}")));
        if let Some(k) = k {
            bc_perm[i] = if k == 0 { -1 } else { k };
        } else {
            let on = bc
                .logical("constraint modes")
                .or_else(|| bc.logical(&format!("constraint modes {name}")));
            if on == Some(true) {
                j += 1;
                bc_perm[i] = j;
            }
        }
    }
    let active = bc_perm.iter().copied().max().unwrap_or(0);
    ctx.info(CALLER, &format!("Number of active constraint modes boundaries: {}", i2s(i64::from(active))), 7);
    if active == 0 {
        return Err(SolverError::fatal(CALLER, "Constraint Modes Analysis requested but no constrained BCs given!"));
    }

    let mut index = vec![0i32; nrows];
    let (nb, nt) = (ctx.mesh.n_bulk(), ctx.mesh.n_bulk() + ctx.mesh.n_boundary());
    for e in nb..nt {
        let e = elmer_core::ElemId(e as u32);
        let Some(bc) = ctx.bc_id(e) else { continue };
        let k = bc_perm[bc - 1];
        if k == 0 {
            continue;
        }
        for &node in ctx.mesh.element(e).nodes {
            let s = perm.get(node as usize).copied().unwrap_or(elmer_core::NO_DOF);
            if s == elmer_core::NO_DOF {
                continue;
            }
            for c in 0..dofs {
                index[dofs * s as usize + c] = dofs as i32 * (k - 1) + c as i32 + 1;
            }
        }
    }
    // Not lumped: every mode dof is a mode of its own, in row order.
    let mut count = 0;
    for v in index.iter_mut().filter(|v| **v > 0) {
        count += 1;
        *v = count;
    }
    ctx.info(CALLER, &format!("Number of active constraint modes: {}", i2s(i64::from(count))), 7);
    for (r, &k) in index.iter().enumerate() {
        if k != 0 {
            fixed.insert(r as u32, 0.0);
        }
    }
    ctx.sim.slots[ctx.slot].cmodes = Some(Modes { index, count: count as usize });
    Ok(())
}

/// The `ConstraintModesDriver` loop: mode `m` moves the one from mode
/// `m - 1`'s dof to its own (Dirichlet row: `b = diag * 1`), solves, and
/// takes it back. Returns the last solve's norm.
pub(crate) fn solve_modes(ctx: &mut SolverCtx<'_>, modes: Modes) -> Result<f64, SolverError> {
    const DRIVER: &str = "ConstraintModesDriver";
    ctx.info(DRIVER, &format!("Number of constraint modes is: {}", i2s(modes.count as i64)), 8);
    let set_rhs = |ctx: &mut SolverCtx<'_>, m: i32, on: bool| -> Result<(), SolverError> {
        let sys = ctx.sim.slots[ctx.slot]
            .system
            .as_mut()
            .ok_or_else(|| SolverError::fatal(DRIVER, "no assembled system"))?;
        for (r, &k) in modes.index.iter().enumerate() {
            if k != m {
                continue;
            }
            let a = &sys.matrix;
            let d = (a.row_ptr[r]..a.row_ptr[r + 1])
                .find(|&p| a.col_idx[p] as usize == r)
                .map_or(0.0, |p| a.values[p]);
            sys.rhs[r] = if on { d } else { 0.0 };
        }
        Ok(())
    };
    let params = &mut ctx.sim.slots[ctx.slot].params;
    let skip = params.get("skip zero rhs test").cloned();
    params.insert("Skip Zero Rhs Test", Value::Logical(true));
    let mut norm = 0.0;
    let mut r = Ok(());
    for m in 1..=modes.count as i32 {
        ctx.info(DRIVER, &format!("Setting up constrained mode: {}", i2s(i64::from(m))), 6);
        r = set_rhs(ctx, m, true).and_then(|_| ctx.solve_system().map(|n| norm = n));
        if r.is_err() {
            break;
        }
        r = set_rhs(ctx, m, false);
        if r.is_err() {
            break;
        }
    }
    let params = &mut ctx.sim.slots[ctx.slot].params;
    match skip {
        Some(v) => params.insert("Skip Zero Rhs Test", v),
        None => params.0.remove("skip zero rhs test"),
    };
    r.map(|_| norm)
}
