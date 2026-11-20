//! `LoadRestartFile` (ModelDescription.F90:4325) on an already parsed
//! `.result` file, and the global `Restart` of `SetInitialConditions`
//! (ElmerSolver.F90:2519) that calls it for the Simulation's `Restart
//! File` after the initial conditions.
//!
//! `LoadRestartFile` reads one stored step (`Restart Position`, the last for
//! 0) into the existing variables, all of them or the ones `Restart Variable
//! i` names (a name matches a file variable with the same start, the rest
//! of the file's name being blanks and digits), stored under `Target
//! Variable i` when given. The value of node `j` goes to the variable's own
//! slot `Perm(j)`; the file's slot is only checked against the field's
//! size. The step's time and timestep go to `Time` and `Timestep` unless
//! `Restart Time Ignore`; `Restart Time` then overrides `Time`.
//!
//! With `Mesh Levels` the file belongs to the finest mesh, the only one
//! here (`Skipping restart for child mesh` for the coarser ones).
//!
//! Refused: variables the file holds but the case lacks (`Restart Create
//! Variables`), `Restart File Cycle`, `Transient Restart`, integration-point
//! variables, a solver's own `Restart File`, `Restart File Number`, `Meshes
//! To Restart` and `Restart Before Initial Conditions`.

use std::path::Path;

use elmer_core::{normalize_key, FieldKind, Logger, ParamMap, Params, NO_DOF};
use elmer_saveio::result::{read_result, PermRecord, ResultFile};

use crate::sim::Simulation;
use crate::variable::VarHandle;
use crate::SolverError;

const LOAD: &str = "LoadRestartFile";

/// The Simulation's `Restart File`, read from `dir` (`Mesh % Name`), after
/// the initial conditions. Nothing without one.
pub fn restart(sim: &mut Simulation, dir: &Path, log: &dyn Logger) -> Result<(), SolverError> {
    if let Some(i) = sim.model.solvers.iter().position(|p| p.get("restart file").is_some()) {
        return Err(SolverError::Unsupported(format!("Solver {}: Restart File (a solver-specific restart)", i + 1)));
    }
    let simp = sim.model.simulation.clone();
    let Some(file) = simp.string("restart file") else { return Ok(()) };
    for key in ["restart file number", "meshes to restart", "restart before initial conditions"] {
        if simp.get(key).is_some() {
            return Err(SolverError::Unsupported(format!("Simulation: {key}")));
        }
    }
    let k = simp.integer("restart position").unwrap_or(0);
    if k < 0 {
        return Err(SolverError::fatal("Restart", "Restart Position must not be negative"));
    }
    let path = dir.join(file.trim());
    let rf = read_result(&path).map_err(|e| SolverError::fatal(LOAD, format!("{}: {e}", path.display())))?;
    log.info(LOAD, &format!("Reading data from file: {}", path.display()), 4);
    load_restart(sim, &rf, k, &simp, log)?;
    if let Some(t) = simp.real("restart time") {
        sim.vars.set_scalar("time", t);
    }
    Ok(())
}

/// `LoadRestartFile(…, TimeCount = k, Mesh)` of a parsed file, with the
/// restart keywords of `res` (the Simulation section): the `k`-th stored
/// step, the last one for `k = 0` or, with a warning, past the end.
pub fn load_restart(sim: &mut Simulation, rf: &ResultFile, k: i32, res: &ParamMap, log: &dyn Logger) -> Result<(), SolverError> {
    if res.get("restart file cycle").is_some() {
        return Err(SolverError::Unsupported("LoadRestartFile: Restart File Cycle".into()));
    }
    if sim.slots.iter().any(|s| s.params.logical("transient restart") == Some(true)) {
        return Err(SolverError::Unsupported("LoadRestartFile: Transient Restart".into()));
    }
    let ignore_time = res.logical("restart time ignore") == Some(true);
    // `Restart Variable 1, 2, ...` up to the first missing one.
    let mut listed: Vec<(String, Option<String>)> = Vec::new();
    while let Some(v) = res.string(&format!("restart variable {}", listed.len() + 1)) {
        let target = res.string(&format!("target variable {}", listed.len() + 1)).map(normalize_key);
        listed.push((normalize_key(v), target));
    }
    let mut found = vec![false; listed.len()];
    // The variable a file name loads into, if any (`LoadThis`, `NewName`).
    let mut target_of = |name: &str| -> Option<String> {
        if listed.is_empty() {
            return Some(name.to_string());
        }
        for (j, (want, target)) in listed.iter().enumerate() {
            let Some(rest) = name.strip_prefix(want.as_str()) else { continue };
            if rest.chars().all(|c| c == ' ' || c.is_ascii_digit()) {
                found[j] = true;
                return Some(target.clone().unwrap_or_else(|| name.to_string()));
            }
        }
        None
    };
    // The header pass: every variable loaded must exist with the file's dofs.
    let mut scalar_target: std::collections::HashMap<String, Option<String>> = std::collections::HashMap::new();
    for h in &rf.header.vars {
        let base = h.name.split('[').next().unwrap_or(&h.name).trim().to_string();
        let Some(new) = target_of(&base) else {
            if h.dofs == 1 {
                scalar_target.insert(base, None);
            }
            continue;
        };
        if h.dofs == 1 {
            scalar_target.insert(base.clone(), Some(new.clone()));
        }
        let Some(handle) = sim.vars.resolve(&new) else {
            return Err(SolverError::Unsupported(format!(
                "LoadRestartFile: creating variable {new} (absent from the case)"
            )));
        };
        let (dofs, size, perm_size) = match handle {
            VarHandle::Primary(id) => {
                let v = sim.vars.by_id(id);
                (v.dofs, v.values.len(), v.perm.as_ref().map(|p| p.len()))
            }
            VarHandle::Component { parent, .. } => {
                let v = sim.vars.by_id(parent);
                (1, v.values.len() / v.dofs.max(1), v.perm.as_ref().map(|p| p.len()))
            }
        };
        if dofs != h.dofs {
            return Err(SolverError::fatal(
                LOAD,
                format!("Fields have different number of components ({} vs. {dofs}): {base}", h.dofs),
            ));
        }
        if h.field_size != size {
            log.warn(LOAD, &format!("Fields are of different size ({} vs. {size}): {base}", h.field_size));
        }
        match perm_size {
            Some(n) if n != h.perm_size => {
                log.warn(LOAD, &format!("Permutations are of different size ({} vs. {n}): {base}", h.perm_size))
            }
            None if h.perm_size > 0 => log.warn(
                LOAD,
                &format!("Existing variable defined without perm: {base} but size in restart file is: {}", h.perm_size),
            ),
            _ => {}
        }
    }
    let n = rf.steps.len();
    if n > 0 {
        let ix = if k <= 0 {
            n - 1
        } else if k as usize > n {
            log.warn(LOAD, "Did not find data at the requested point");
            log.warn(LOAD, &format!("Reading from the last existing point: {n:>6}"));
            n - 1
        } else {
            k as usize - 1
        };
        let step = &rf.steps[ix];
        if !ignore_time {
            sim.vars.set_scalar("time", step.time);
            sim.vars.set_scalar("timestep", f64::from(step.timestep));
        }
        log.info(LOAD, &format!("Reading time sequence: {}", elmer_saveio::fortran::fortran_es(step.time, 12, 3)), 4);
        log.info(LOAD, &format!("Reading variables on timestep: {}", step.timestep), 4);
        for rec in &step.vars {
            let Some(Some(new)) = scalar_target.get(&rec.name) else { continue };
            load_one(sim, rf, ix, &rec.name, new, &rec.perm)?;
        }
    } else if k != 0 {
        log.warn(LOAD, "Did not find data at the requested point");
    }
    for (j, f) in found.iter().enumerate() {
        if !f {
            log.warn(LOAD, &format!("Could not find restart variable: {}", j + 1));
        }
    }
    log.info(LOAD, "All done", 4);
    Ok(())
}

/// The values of file variable `name` of step `ix` into case variable `new`.
fn load_one(sim: &mut Simulation, rf: &ResultFile, ix: usize, name: &str, new: &str, perm: &PermRecord) -> Result<(), SolverError> {
    let f = rf
        .restart_field(ix, name)
        .map_err(|e| SolverError::fatal(LOAD, e.to_string()))?
        .ok_or_else(|| SolverError::fatal(LOAD, format!("{name} is not in the step")))?;
    let Some(handle) = sim.vars.resolve(new) else {
        return Err(SolverError::fatal(LOAD, format!("Variable is not present for reading: {new}")));
    };
    let (id, offset, stride) = match handle {
        VarHandle::Primary(id) => (id, 0, 1),
        VarHandle::Component { parent, component } => (parent, component, sim.vars.by_id(parent).dofs),
    };
    let var = sim.vars.by_id_mut(id);
    if !matches!(perm, PermRecord::Null) != var.perm.is_some() {
        return Err(SolverError::fatal(LOAD, "Permutation should either exist or not!"));
    }
    if var.kind == FieldKind::Ip {
        return Err(SolverError::Unsupported("LoadRestartFile: integration-point variables".into()));
    }
    let size = var.values.len() / stride;
    let vperm = var.perm.clone();
    if vperm.as_ref().is_some_and(|p| p.iter().all(|&s| s == NO_DOF)) {
        return Err(SolverError::Unsupported("LoadRestartFile: a variable without any dof (new Perm)".into()));
    }
    for ((&node, &slot), &val) in f.node.iter().zip(&f.slot).zip(&f.values) {
        let (node, slot) = (node as usize, slot as usize);
        if slot >= size {
            continue;
        }
        let at = match &vperm {
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
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Model;
    use crate::variable::Variable;
    use elmer_core::{NullLogger, Value};
    use elmer_saveio::result::{ResultFormat, ResultHeader, StepRecord, VarHeader, VarRecord};
    use std::sync::Arc;

    fn header(name: &str, perm_size: usize) -> VarHeader {
        VarHeader { name: name.into(), field_size: 2, perm_size, dofs: 1, equation: "e".into() }
    }

    /// `MeshLevelsRestart`'s `Restart Variable 1 = Temperature` over a file
    /// with coordinates and globals the case lacks: only the listed
    /// variable loads (through the case's own perm), into `Target Variable
    /// 1` when given, and the step's time and timestep carry over.
    #[test]
    fn listed_variables_load_through_the_case_perm() {
        let table = PermRecord::Table { size: 2, entries: vec![(1, 1), (2, 2)] };
        let rf = ResultFile {
            header: ResultHeader {
                format: ResultFormat::Ascii,
                started_at: String::new(),
                vars: vec![header("coordinate 1", 0), header("temperature", 2), header("temperature 2", 2)],
                n_nodes: 2,
            },
            steps: vec![StepRecord {
                saved_count: 1,
                timestep: 10,
                time: 20.0,
                vars: vec![
                    VarRecord { name: "coordinate 1".into(), perm: PermRecord::Null, values: vec![0.0, 1.0] },
                    VarRecord { name: "temperature".into(), perm: table, values: vec![3.0, 4.0] },
                    VarRecord { name: "temperature 2".into(), perm: PermRecord::UsePrevious, values: vec![5.0, 6.0] },
                ],
            }],
        };
        let run = |target: Option<&str>| {
            let mut res = ParamMap::default();
            res.insert("Restart Variable 1", Value::Str("Temperature".into()));
            if let Some(t) = target {
                res.insert("Target Variable 1", Value::Str(t.into()));
            }
            let mut sim = Simulation::new(Model::default());
            let perm: Arc<[u32]> = Arc::from(vec![1u32, 0]);
            for n in ["temperature", "temperature 2", "other"] {
                sim.vars.add(Variable::new(n, FieldKind::Nodal, 1, Some(perm.clone()), vec![0.0, 0.0])).unwrap();
            }
            load_restart(&mut sim, &rf, 0, &res, &NullLogger).unwrap();
            sim
        };
        let sim = run(None);
        // "temperature 2" matches "temperature" too: the rest is a blank and digits.
        assert_eq!(sim.vars.get("temperature").unwrap().values, vec![4.0, 3.0]);
        assert_eq!(sim.vars.get("temperature 2").unwrap().values, vec![6.0, 5.0]);
        assert_eq!((sim.vars.scalar("time"), sim.vars.scalar("timestep")), (Some(20.0), Some(10.0)));
        let sim = run(Some("Other"));
        assert_eq!(sim.vars.get("other").unwrap().values, vec![6.0, 5.0]);
        assert_eq!(sim.vars.get("temperature").unwrap().values, vec![0.0, 0.0]);
    }
}
