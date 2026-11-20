//! `ReloadSolution` of the `ReloadData` module (fem/src/modules/
//! ReloadData.F90) and the `LoadRestartFile` it calls
//! (ModelDescription.F90:4325).
//!
//! Each visit reads one stored step of `Reload Solution File` (under the
//! mesh directory, `Mesh % Name`) into the existing variables. The step
//! cycles over `Reload Range Minimum`..`Maximum`, starting from `Reload
//! Starting Position`, every `Reload Reading Intervals` visits. The loaded
//! step also sets the `Time` and `Timestep` globals, unless the Simulation
//! has `Restart Time Ignore`, and the main loop carries on from them.
//!
//! The file may be a version 2 or 3 `.result`, ASCII or binary. Values go
//! in as `LoadRestartFile` puts them: the value of node `j` goes to the
//! variable's own slot `Perm(j)`; the file's slot is only checked against
//! the field's size.
//!
//! Refused: variables the file holds but the case lacks (`Restart Create
//! Variables`), `Restart Variable i`, `Restart File Cycle`, `Transient
//! Restart`, and runs on more than one partition.

use std::path::PathBuf;

use elmer_core::{ParamMap, Params, NO_DOF};
use elmer_saveio::result::{read_result, ResultFile};
use elmer_solver::variable::VarHandle;
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "ReloadSolution";
const LOAD: &str = "LoadRestartFile";

pub(crate) struct ReloadSolution {
    file: Option<String>,
    first: i32,
    last: i32,
    between: i32,
    round: i32,
    visit: i32,
    loaded: Option<(PathBuf, ResultFile)>,
}

impl ReloadSolution {
    pub(crate) fn from_params(p: &ParamMap) -> Result<Self, SolverError> {
        let first = p.integer("reload range minimum").unwrap_or(1);
        let last = p.integer("reload range maximum").unwrap_or(100_000);
        let mut round = 0;
        if let Some(start) = p.integer("reload starting position") {
            round += start - first;
        }
        let between = p.integer("reload reading intervals").unwrap_or(1);
        if between < 1 {
            return Err(SolverError::fatal(CALLER, "Reload Reading Intervals must be positive"));
        }
        if last < first {
            return Err(SolverError::fatal(CALLER, "Reload Range Maximum is below the minimum"));
        }
        Ok(ReloadSolution {
            file: p.string("reload solution file").map(str::to_string),
            first,
            last,
            between,
            round,
            visit: 0,
            loaded: None,
        })
    }
}

impl Solver for ReloadSolution {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let visit = self.visit;
        self.visit += 1;
        if visit % self.between != 0 {
            return Ok(());
        }
        let Some(file) = self.file.clone() else { return Ok(()) };
        if ctx.comm.size() > 1 {
            return Err(SolverError::Unsupported("ReloadSolution: parallel runs".into()));
        }
        if visit == 0 {
            ctx.info(CALLER, &format!(" Reading old solution from file {file}"), 4);
        }
        self.round += 1;
        let span = 1 + self.last - self.first;
        let k = match self.round % span {
            0 => self.last,
            r => r + self.first - 1,
        };
        let path = match ctx.sim.model.header.string("mesh name") {
            Some(m) if !m.trim().is_empty() => format!("{}/{file}", m.trim()),
            _ => file,
        };
        ctx.info(CALLER, &format!(" Opening file: {path}"), 4);
        ctx.info(CALLER, &format!(" Loading Timestep: {k:>12}"), 4);
        let path = PathBuf::from(path);
        if self.loaded.as_ref().is_none_or(|(p, _)| *p != path) {
            let rf = read_result(&path).map_err(|e| SolverError::fatal(LOAD, format!("{}: {e}", path.display())))?;
            self.loaded = Some((path, rf));
        }
        let (_, rf) = self.loaded.as_ref().expect("just loaded");
        load_restart(ctx, rf, k)
    }
}

/// `LoadRestartFile(…, TimeCount = k, Mesh)` of an already read file: the
/// `k`-th stored step (the last one for `k = 0` or past the end).
fn load_restart(ctx: &mut SolverCtx<'_>, rf: &ResultFile, k: i32) -> Result<(), SolverError> {
    let simp = ctx.simulation_params();
    for key in ["restart variable 1", "restart file cycle"] {
        if simp.get(key).is_some() {
            return Err(SolverError::Unsupported(format!("LoadRestartFile: {key}")));
        }
    }
    let ignore_time = simp.logical("restart time ignore") == Some(true);
    if ctx.sim.slots.iter().any(|s| s.params.logical("transient restart") == Some(true)) {
        return Err(SolverError::Unsupported("LoadRestartFile: Transient Restart".into()));
    }
    // The header pass: every variable must exist with the file's dofs.
    for h in &rf.header.vars {
        let Some(handle) = ctx.sim.vars.resolve(&h.name) else {
            return Err(SolverError::Unsupported(format!(
                "LoadRestartFile: creating variable {} (absent from the case)",
                h.name
            )));
        };
        let (dofs, size, perm_size) = match handle {
            VarHandle::Primary(id) => {
                let v = ctx.sim.vars.by_id(id);
                (v.dofs, v.values.len(), v.perm.as_ref().map(|p| p.len()))
            }
            VarHandle::Component { parent, .. } => {
                let v = ctx.sim.vars.by_id(parent);
                (1, v.values.len() / v.dofs.max(1), v.perm.as_ref().map(|p| p.len()))
            }
        };
        if dofs != h.dofs {
            return Err(SolverError::fatal(
                LOAD,
                format!(
                    "Fields have different number of components ({} vs. {dofs}): {}",
                    h.dofs, h.name
                ),
            ));
        }
        if h.field_size != size {
            ctx.warn(LOAD, &format!("Fields are of different size ({} vs. {size}): {}", h.field_size, h.name));
        }
        match perm_size {
            Some(n) if n != h.perm_size => ctx.warn(
                LOAD,
                &format!("Permutations are of different size ({} vs. {n}): {}", h.perm_size, h.name),
            ),
            None if h.perm_size > 0 => ctx.warn(
                LOAD,
                &format!(
                    "Existing variable defined without perm: {} but size in restart file is: {}",
                    h.name, h.perm_size
                ),
            ),
            _ => {}
        }
    }
    let n = rf.steps.len();
    if n == 0 {
        if k != 0 {
            ctx.warn(LOAD, "Did not find data at the requested point");
        }
        return Ok(());
    }
    let ix = if k <= 0 {
        n - 1
    } else if k as usize > n {
        ctx.warn(LOAD, "Did not find data at the requested point");
        ctx.warn(LOAD, &format!("Reading from the last existing point: {n:>6}"));
        n - 1
    } else {
        k as usize - 1
    };
    let step = &rf.steps[ix];
    if !ignore_time {
        ctx.sim.vars.set_scalar("time", step.time);
        ctx.sim.vars.set_scalar("timestep", f64::from(step.timestep));
    }
    ctx.info(LOAD, &format!("Reading time sequence: {}", elmer_saveio::fortran::fortran_es(step.time, 12, 3)), 4);
    ctx.info(LOAD, &format!("Reading variables on timestep: {}", step.timestep), 4);
    for rec in &step.vars {
        let f = rf
            .restart_field(ix, &rec.name)
            .map_err(|e| SolverError::fatal(LOAD, e.to_string()))?
            .expect("the step holds the variable");
        let Some(handle) = ctx.sim.vars.resolve(&rec.name) else {
            return Err(SolverError::fatal(LOAD, format!("Variable is not present for reading: {}", rec.name)));
        };
        let (id, offset, stride) = match handle {
            VarHandle::Primary(id) => (id, 0, 1),
            VarHandle::Component { parent, component } => (parent, component, ctx.sim.vars.by_id(parent).dofs),
        };
        let has_table = !matches!(rec.perm, elmer_saveio::result::PermRecord::Null);
        let var = ctx.sim.vars.by_id_mut(id);
        if has_table != var.perm.is_some() {
            return Err(SolverError::fatal(LOAD, "Permutation should either exist or not!"));
        }
        if var.kind == elmer_core::FieldKind::Ip {
            return Err(SolverError::Unsupported("LoadRestartFile: integration-point variables".into()));
        }
        let size = var.values.len() / stride;
        let perm = var.perm.clone();
        if perm.as_ref().is_some_and(|p| p.iter().all(|&s| s == NO_DOF)) {
            return Err(SolverError::Unsupported("LoadRestartFile: a variable without any dof (new Perm)".into()));
        }
        for ((&node, &slot), &val) in f.node.iter().zip(&f.slot).zip(&f.values) {
            let (node, slot) = (node as usize, slot as usize);
            if slot >= size {
                continue;
            }
            let at = match &perm {
                None => node,
                Some(p) => match p.get(node) {
                    Some(&s) if s != NO_DOF => s as usize,
                    _ => continue,
                },
            };
            if at < size {
                var.values[offset + stride * at] = val;
            }
        }
    }
    ctx.info(LOAD, "All done", 4);
    Ok(())
}
