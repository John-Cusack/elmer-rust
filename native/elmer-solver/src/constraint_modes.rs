//! `Constraint Modes Analysis` on boundary conditions: every dof of a
//! boundary with `Constraint Modes <var> = True` (or `Constraint Modes`) is
//! a mode (`SetConstraintModesBoundaries`, BoundaryConditionUtils.F90:2765):
//! all of them are held at zero, and the linear system is solved once per
//! mode with that dof at one (with `Constraint Modes Lumped`, a mode per
//! boundary and component, all its dofs at one) (`ConstraintModesDriver` around
//! `SolveLinearSystem`, SolveCore.F90:3215-3290, 3995-4290), each solve with
//! its own nonlinear change. `Constraint Mode <var> = Integer k` numbers a
//! boundary (0 is the ground, held at zero).
//!
//! With `Constraint Modes Fluxes` each mode's fluxes `A_bulk x` summed over
//! the mode dofs make a row of the lumped matrix (`ConstraintModesFluxes`,
//! `StoreLumpedFluxes`), printed after the last mode (`FinalizeLumpedMatrix`,
//! SolveCore.F90:3526-3840), optionally symmetrised, its norm optionally
//! stored as `res: <var> lumped matrix norm` and its entries as `res: CMF i
//! j`.
//!
//! Refused: modes set by `mode <n>: <component>` conditions, the `rhs`,
//! `EM Wave` and `Coils` variants,
//! `Constraint Modes Single` and mode directions, complex systems, modes over
//! body forces or components, the external loops (`Nonlinear System`,
//! `Steady State` and `Run Control Constraint Modes`), linkage, writing the
//! lumped matrix to a file, symmetric Dirichlet elimination, p elements and
//! parallel runs.

use std::collections::BTreeMap;

use elmer_core::{ParamMap, Params, Value};

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

/// A solver's "compute a lumped matrix" flag, which its `_init` expands
/// into the library's constraint-mode lumping ([`lumping_init`]).
pub struct Lumping {
    /// The solver flag, e.g. `Calculate Conductance Matrix`.
    pub flag: &'static str,
    /// The default `Constraint Modes Matrix Filename`.
    pub filename: &'static str,
    /// A solver keyword renamed `Constraint Modes Matrix Filename` when given.
    pub filename_key: Option<&'static str>,
    /// The BC (and, with `body_forces`, Body Force) keyword renamed
    /// `Constraint Mode <var>`: (from, to).
    pub rename: (&'static str, &'static str),
    pub body_forces: bool,
    /// `Optimize Bandwidth = False` is set (`ListAddLogical`, overriding).
    pub no_bandwidth: bool,
}

/// `HeatSolver_init` (HeatSolve.F90:201).
pub const CONDUCTANCE: Lumping = Lumping {
    flag: "calculate conductance matrix",
    filename: "ThermalConductanceMatrix.dat",
    filename_key: None,
    rename: ("conductivity body", "Constraint Mode Temperature"),
    body_forces: false,
    no_bandwidth: false,
};

/// `StatCurrentSolver_init` (StatCurrentSolve.F90:186).
pub const CONDUCTIVITY: Lumping = Lumping {
    flag: "calculate conductivity matrix",
    filename: "ConductivityMatrix.dat",
    filename_key: None,
    rename: ("conductivity body", "Constraint Mode Potential"),
    body_forces: false,
    no_bandwidth: false,
};

/// `StatElecSolver_init` (StatElecSolve.F90:211).
pub const CAPACITANCE: Lumping = Lumping {
    flag: "calculate capacitance matrix",
    filename: "CapacitanceMatrix.dat",
    filename_key: Some("capacitance matrix filename"),
    rename: ("capacitance body", "Constraint Mode Potential"),
    body_forces: true,
    no_bandwidth: true,
};

/// The `_init` expansion of `l.flag`, on the solver section `p` and the
/// model's BCs and body forces; nothing when the flag is off. `Constraint
/// Modes Analysis`, `... Lumped`, `... Fluxes` and `... Matrix Symmetric`
/// are added unless given (`ListAddNewLogical`), the matrix filename is
/// renamed from `l.filename_key` or defaults to `l.filename`, and the
/// BC / body force keyword is renamed (`ListRenameAllBC`,
/// `ListRenameAllBodyForce`). Elmer runs this before the solver's variable
/// and permutation exist (`Optimize Bandwidth`), so a caller building the
/// model calls it before [`Simulation::setup`](crate::Simulation::setup).
/// Returns whether the flag is on; a second call changes nothing.
pub fn lumping_init(p: &mut ParamMap, bcs: &mut [ParamMap], body_forces: &mut [ParamMap], l: &Lumping) -> bool {
    if p.logical(l.flag) != Some(true) {
        return false;
    }
    for k in [
        "Constraint Modes Analysis",
        "Constraint Modes Lumped",
        "Constraint Modes Fluxes",
        "Constraint Modes Matrix Symmetric",
    ] {
        if p.get(k).is_none() {
            p.insert(k, Value::Logical(true));
        }
    }
    if let Some(v) = l.filename_key.and_then(|k| p.0.remove(k)) {
        p.insert("Constraint Modes Matrix Filename", v);
    } else if p.get("constraint modes matrix filename").is_none() {
        p.insert("Constraint Modes Matrix Filename", Value::Str(l.filename.to_string()));
    }
    let (from, to) = l.rename;
    let n = if l.body_forces { body_forces.len() } else { 0 };
    for sec in bcs.iter_mut().chain(body_forces[..n].iter_mut()) {
        if let Some(v) = sec.0.remove(from) {
            sec.insert(to, v);
        }
    }
    if l.no_bandwidth {
        p.insert("Optimize Bandwidth", Value::Logical(false));
    }
    true
}

/// [`lumping_init`] on the running solver, for a model built without it.
/// `Optimize Bandwidth` then comes too late for a variable created by
/// setup: its dofs keep the optimized order (the same lumped matrix up to
/// rounding).
pub fn constraint_modes_lumping(ctx: &mut SolverCtx<'_>, l: &Lumping) {
    let sim = &mut *ctx.sim;
    let model = &mut sim.model;
    lumping_init(&mut sim.slots[ctx.slot].params, &mut model.boundary_conditions, &mut model.body_forces, l);
}

fn refuse(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
    let p = ctx.params();
    for k in [
        "constraint modes rhs",
        "constraint modes em wave",
        "constraint modes coils",
        "constraint modes single",
        "nonlinear system constraint modes",
        "steady state constraint modes",
        "run control constraint modes",
        "constraint modes linkage",
        "consider p fluxes",
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
    // Lumped modes set by `mode <n>: <component>` boundary conditions.
    if m.boundary_conditions.iter().any(|b| b.0.keys().any(|k| k.starts_with("mode "))) {
        return Err(SolverError::Unsupported("constraint modes set by mode <n>: boundary conditions".into()));
    }
    if m.boundary_conditions.iter().chain(&m.body_forces).any(|b| b.get("constraint mode direction").is_some()) {
        return Err(SolverError::Unsupported("Constraint Mode Direction".into()));
    }
    if p.get("nonlinear post solvers").is_some() {
        return Err(SolverError::Unsupported("Constraint Modes Analysis with Nonlinear Post Solvers".into()));
    }
    if ctx.comm.size() > 1 {
        return Err(SolverError::Unsupported("constraint modes in parallel".into()));
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
    // BCMode / BFMode: the modes are numbered over the boundary conditions
    // when any has a `Constraint Mode...` keyword (ListCheckPrefixAnyBC),
    // else over the body forces.
    let has = |s: &[elmer_core::ParamMap]| s.iter().any(|b| b.0.keys().any(|k| k.starts_with("constraint mode")));
    let bc_mode = has(&model.boundary_conditions);
    if !bc_mode && !has(&model.body_forces) {
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
    let entities = if bc_mode { &model.boundary_conditions } else { &model.body_forces };
    let mut bc_perm = vec![0i32; entities.len()];
    let mut j = 0;
    for (i, bc) in entities.iter().enumerate() {
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

    // `HaveP` is whether the first bulk element is an active p element;
    // `Ignore Constraint Modes p` then leaves the p dofs out (`nb = n`).
    let have_p = ctx.mesh.n_bulk() > 0 && ctx.solver_dofs().is_some_and(|d| d.is_p_element(elmer_core::ElemId(0)));
    let ignore_p = have_p && ctx.params().logical("ignore constraint modes p") == Some(true);
    let mut index = vec![0i32; nrows];
    let (nb, nt) = (ctx.mesh.n_bulk(), ctx.mesh.n_bulk() + ctx.mesh.n_boundary());
    // BCMode: the boundary elements; BFMode: every bulk element, active
    // or not, through its body's `Body Force`.
    let range = if bc_mode { nb..nt } else { 0..nb };
    let mut slots = Vec::new();
    for e in range {
        let e = elmer_core::ElemId(e as u32);
        let ent = if bc_mode {
            ctx.bc_id(e)
        } else {
            let model = &ctx.sim.model;
            model
                .body_of_tag(ctx.mesh.element(e).tag)
                .and_then(|b| model.bodies.get(b - 1))
                .and_then(|b| b.integer("body force"))
                .filter(|&f| f > 0)
                .map(|f| f as usize)
        };
        let Some(ent) = ent else { continue };
        let Some(&k) = bc_perm.get(ent - 1) else { continue };
        if k == 0 {
            continue;
        }
        // `Perm(Indexes(1:nb))` of `mGetElementDOFs`: the nodes, then any
        // p dofs of the element.
        let nodes = ctx.mesh.element(e).nodes;
        match ctx.solver_dofs() {
            Some(d) => {
                d.element_slots(e, &mut slots);
                if ignore_p {
                    slots.truncate(nodes.len());
                }
            }
            None => {
                slots.clear();
                slots.extend(nodes.iter().map(|&n| perm.get(n as usize).copied().unwrap_or(elmer_core::NO_DOF)));
            }
        }
        for &s in &slots {
            if s == elmer_core::NO_DOF {
                continue;
            }
            for c in 0..dofs {
                index[dofs * s as usize + c] = dofs as i32 * (k - 1) + c as i32 + 1;
            }
        }
    }
    // The p dofs (perm entries past the nodes) are shifted by
    // `2*(NoModes+1)` below zero: held at zero in every mode, apart from
    // the nodal mode dofs (BoundaryConditionUtils.F90:3047).
    if have_p && !ignore_p {
        let poffset = 2 * (dofs as i32 * active + 1);
        for &s in perm.iter().skip(ctx.mesh.n_nodes()) {
            if s == elmer_core::NO_DOF {
                continue;
            }
            for c in 0..dofs {
                let v = &mut index[dofs * s as usize + c];
                if *v != 0 {
                    *v -= poffset;
                }
            }
        }
    }
    // Lumped: a mode per boundary and component; otherwise every mode dof
    // is a mode of its own, in row order.
    let count = if ctx.params().logical("constraint modes lumped") == Some(true) {
        dofs as i32 * active
    } else {
        let mut count = 0;
        for v in index.iter_mut().filter(|v| **v > 0) {
            count += 1;
            *v = count;
        }
        ctx.info(CALLER, &format!("Number of active constraint modes: {}", i2s(i64::from(count))), 7);
        count
    };
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
    let fluxes = ctx.params().logical("constraint modes fluxes") == Some(true);
    let mut lumped: Vec<Vec<f64>> = Vec::new();
    let mut norm = 0.0;
    let mut r = Ok(());
    for m in 1..=modes.count as i32 {
        ctx.info(DRIVER, &format!("Setting up constrained mode: {}", i2s(i64::from(m))), 6);
        r = set_rhs(ctx, m, true).and_then(|_| ctx.solve_system().map(|n| norm = n));
        if r.is_err() {
            break;
        }
        // SolveCore.F90:4250: `Var % ConstraintModes(NMode,:) = x`.
        if let Some(id) = ctx.sim.slots[ctx.slot].variable {
            let v = ctx.sim.vars.by_id_mut(id);
            let mode = v.values.clone();
            let k = m as usize;
            v.constraint_modes.resize(modes.count.max(k), Vec::new());
            v.constraint_modes[k - 1] = mode;
        }
        ctx.info(DRIVER, &format!("Mode {} computed, doing some postprocessing", i2s(i64::from(m))), 10);
        r = set_rhs(ctx, m, false);
        if r.is_err() {
            break;
        }
        if fluxes {
            ctx.info(DRIVER, "Computing lumped fluxes", 8);
            match flux_row(ctx, &modes, m) {
                Ok(row) => {
                    ctx.info("StoreLumpedFluxes", "Storing lumped fluxes", 10);
                    if lumped.is_empty() {
                        ctx.info(
                            "StoreLumpedFluxes",
                            &format!("Allocating lumped model of size: {}", i2s(modes.count as i64)),
                            5,
                        );
                    }
                    lumped.push(row);
                }
                Err(e) => {
                    r = Err(e);
                    break;
                }
            }
        }
    }
    if r.is_ok() && fluxes {
        finalize(ctx, lumped);
    }
    let params = &mut ctx.sim.slots[ctx.slot].params;
    match skip {
        Some(v) => params.insert("Skip Zero Rhs Test", v),
        None => params.0.remove("skip zero rhs test"),
    };
    r.map(|_| norm)
}

/// `ConstraintModesFluxes` (real): `A_bulk x` over the mode dofs; mode `n`'s
/// row gets the fluxes of every mode dof, each other mode's entry minus its
/// own.
fn flux_row(ctx: &SolverCtx<'_>, modes: &Modes, n: i32) -> Result<Vec<f64>, SolverError> {
    let slot = &ctx.sim.slots[ctx.slot];
    let sys = slot.system.as_ref().ok_or_else(|| SolverError::fatal("ConstraintModesDriver", "no assembled system"))?;
    let (values, _) = slot
        .asm
        .bulk
        .as_ref()
        .ok_or_else(|| SolverError::fatal("ConstraintModesDriver", "BulkValues not associated!"))?;
    let a = &sys.matrix;
    let id = slot.variable.expect("set up with a variable");
    let x = &ctx.sim.vars.by_id(id).values;
    let mut row = vec![0.0; modes.count];
    for (j, &k) in modes.index.iter().enumerate() {
        if k <= 0 {
            continue;
        }
        let flux: f64 = (a.row_ptr[j]..a.row_ptr[j + 1]).map(|p| values[p] * x[a.col_idx[p] as usize]).sum();
        if k != n {
            row[k as usize - 1] -= flux;
        }
        row[n as usize - 1] += flux;
    }
    Ok(row)
}

/// `FinalizeLumpedMatrix` (real).
fn finalize(ctx: &mut SolverCtx<'_>, mut c: Vec<Vec<f64>>) {
    const FIN: &str = "FinalizeLumpedMatrix";
    ctx.info(FIN, "Finalizing lumped matrix", 8);
    let n = c.len();
    let p = ctx.params();
    let symmetric = p
        .logical("constraint modes fluxes symmetric")
        .or_else(|| p.logical("constraint modes matrix symmetric"))
        == Some(true);
    if symmetric {
        ctx.info(FIN, "Enforcing symmetry of reduced system!", 8);
        for i in 0..n {
            for j in i + 1..n {
                ctx.info(FIN, &format!("{:>3}{:>3}{}", i + 1, j + 1, crate::fortran_fmt::es(c[i][j] - c[j][i], 17, 9)), 10);
            }
        }
        let t = c.clone();
        for i in 0..n {
            for j in 0..n {
                c[i][j] = 0.5 * (t[i][j] + t[j][i]);
            }
        }
    }
    ctx.info(FIN, "Lumped Matrix", 5);
    for (i, row) in c.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            if symmetric && j < i {
                continue;
            }
            ctx.info(FIN, &format!("{:>3}{:>3}{}", i + 1, j + 1, crate::fortran_fmt::es(*v, 17, 9)), 5);
        }
    }
    // `Constraint Modes Fluxes Filename` (or `... Matrix Filename`): the
    // matrix one row per record, list-directed (`WRITE(10,*)
    // FluxesMatrix(i,:)`), in the working directory.
    let file = p
        .string("constraint modes fluxes filename")
        .or_else(|| p.string("constraint modes matrix filename"))
        .map(str::to_string);
    if let Some(file) = file {
        let text: String = c
            .iter()
            .map(|row| row.iter().map(|&v| elmer_saveio::fortran::list_real(v)).collect::<String>() + "\n")
            .collect();
        if let Err(e) = std::fs::write(&file, text) {
            ctx.warn(FIN, &format!("could not write {file}: {e}"));
        }
    }
    let nrm: f64 = c.iter().flatten().map(|v| v.abs()).sum();
    let name = ctx
        .sim
        .slots[ctx.slot]
        .variable
        .map(|id| ctx.sim.vars.by_id(id).name.clone())
        .unwrap_or_default();
    ctx.info(FIN, &format!("{name} lumped matrix norm: {}", crate::fortran_fmt::es(nrm, 12, 5)), 6);
    let p = ctx.params();
    let store_norm = p
        .logical("constraint modes fluxes norm")
        .or_else(|| p.logical("constraint modes matrix norm"))
        == Some(true);
    let results = p.logical("constraint modes matrix results") == Some(true);
    if results {
        ctx.info(FIN, "Adding Constraint Modes Fluxes with \"res:\" to list", 5);
    }
    let sim = &mut ctx.sim.model.simulation;
    if store_norm {
        sim.insert(&format!("res: {name} lumped matrix norm"), Value::real(nrm));
    }
    if results {
        for (i, row) in c.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                sim.insert(&format!("res: CMF {} {}", i + 1, j + 1), Value::real(*v));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `Calculate Capacitance Matrix` (StatElecSolve.F90:211): the
    /// constraint-mode keywords are added unless given, `Capacitance Matrix
    /// Filename` becomes the matrix filename, `Capacitance Body` of BCs and
    /// body forces becomes `Constraint Mode Potential`, bandwidth
    /// optimization is off; a second call changes nothing.
    #[test]
    fn capacitance_flag_expands_to_constraint_mode_lumping() {
        let mut p = ParamMap::new();
        p.insert("Calculate Capacitance Matrix", Value::Logical(true));
        p.insert("Constraint Modes Lumped", Value::Logical(false));
        p.insert("Capacitance Matrix Filename", Value::Str("c.dat".into()));
        let mut bc = ParamMap::new();
        bc.insert("Capacitance Body", Value::Integer(vec![2]));
        let mut bf = ParamMap::new();
        bf.insert("Capacitance Body", Value::Integer(vec![1]));
        let (mut bcs, mut bfs) = (vec![bc], vec![bf]);
        for _ in 0..2 {
            assert!(lumping_init(&mut p, &mut bcs, &mut bfs, &CAPACITANCE));
            assert_eq!(p.logical("constraint modes analysis"), Some(true));
            assert_eq!(p.logical("constraint modes lumped"), Some(false));
            assert_eq!(p.string("constraint modes matrix filename"), Some("c.dat"));
            assert!(p.get("capacitance matrix filename").is_none());
            assert_eq!(p.logical("optimize bandwidth"), Some(false));
            assert_eq!(bcs[0].integer("constraint mode potential"), Some(2));
            assert_eq!(bfs[0].integer("constraint mode potential"), Some(1));
            assert!(bcs[0].get("capacitance body").is_none());
        }
        // HeatSolve: BCs only, the default filename, bandwidth untouched.
        let mut h = ParamMap::new();
        h.insert("Calculate Conductance Matrix", Value::Logical(true));
        let mut bfs = vec![ParamMap::new()];
        bfs[0].insert("Conductivity Body", Value::Integer(vec![1]));
        assert!(lumping_init(&mut h, &mut [], &mut bfs, &CONDUCTANCE));
        assert_eq!(h.string("constraint modes matrix filename"), Some("ThermalConductanceMatrix.dat"));
        assert!(h.get("optimize bandwidth").is_none());
        assert_eq!(bfs[0].integer("conductivity body"), Some(1));
        assert!(!lumping_init(&mut ParamMap::new(), &mut [], &mut [], &CONDUCTANCE));
    }
}
