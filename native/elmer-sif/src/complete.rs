//! Model completion replicated explicitly in `elmer-sif` (report §3.3):
//! legacy-solver injection (`AddEquationBasics`, `MainUtils.F90`), the
//! `Post File` / `Scalars File` output-solver hacks (`ElmerSolver.F90`),
//! and the small `LoadModel` additions (`Solver Input File`, obsolete
//! `Transmittivity` check from `CompleteModelKeywords`).
//!
//! Mortar/contact detection and other mesh-dependent completion belong to
//! later lanes and are not done here.

use crate::model::{Entry, Model, Sif, Value};

/// Apply all completion passes to every run block. Errors are fatal in
/// Elmer (`Fatal`), so they are returned.
pub fn complete(sif: &mut Sif, sif_path: &str) -> Result<(), String> {
    for run in &mut sif.runs {
        complete_run(run, sif_path)?;
    }
    Ok(())
}

/// Complete one run block (one `RUN` segment). Used by the corpus proof to
/// attribute completion failures to the right file.
pub fn complete_run(run: &mut Model, sif_path: &str) -> Result<(), String> {
    legacy_solver_injection(run);
    add_vtu_output_solver(run);
    add_save_scalars_solver(run);
    if run.simulation.is_none() {
        run.simulation = Some(crate::model::Section::default());
    }
    if let Some(sim) = run.simulation.as_mut() {
        sim.entries.push(Entry {
            name: "solver input file".to_string(),
            value: Value::Text { text: sif_path.to_string(), is_file: true },
            partag: None,
            distribute: false,
        });
    }
    check_obsolete_transmittivity(run)?;
    Ok(())
}



fn find_solver_by_proc(solvers: &[crate::model::Section], proc: &str) -> Option<usize> {
    solvers.iter().position(|s| {
        matches!(
            s.get("procedure"),
            Some(Entry { value: Value::Text { text, .. }, .. })
            if text.split_whitespace().last().unwrap_or("") == proc
                || text.split_whitespace().next().unwrap_or("") == proc
        )
    })
}

fn append_solver(run: &mut Model) -> &mut crate::model::Section {
    let id = run.solvers.iter().map(|s| s.id).max().unwrap_or(0) + 1;
    run.solvers.push(crate::model::Section { id, entries: Vec::new() });
    run.solvers.last_mut().unwrap()
}

/// Legacy solvers (`MainUtils.F90:AddEquationBasics`): an `Equation`
/// section whose `Name` is a historical solver name gains the default
/// `Variable` and `Procedure` keywords when absent.
fn legacy_solver_injection(run: &mut Model) {
    // (equation name, variable, dofs-expr, procedure)
    for eq in &mut run.equations {
        let Some(name_entry) = eq.get("name") else { continue };
        let Value::Text { text: eq_name, .. } = &name_entry.value else { continue };
        let (var, proc): (String, &str) = match eq_name.trim().to_lowercase().as_str() {
            "heat equation" => ("temperature".to_string(), "heatsolve heatsolver"),
            "navier-stokes" => ("flow solution[velocity:2 pressure:1]".to_string(), "flowsolve flowsolver"),
            "magnetic induction" => ("magnetic field".to_string(), "magneticsolve magneticsolver"),
            "stress analysis" => ("displacement".to_string(), "stresssolve stresssolver"),
            "mesh update" => ("mesh update".to_string(), "meshsolve meshsolver"),
            _ => continue,
        };
        // NOTE: Elmer keys legacy behaviour off the *solver's* Equation
        // string, not the Equation section Name; the section-level Name
        // form here covers the corpus usage (`Equation = "Heat Equation"`
        // style inputs route through the same defaults when a solver
        // references them). Solver-level strings are handled by the
        // solver lane; see `elmer-solver`.
        if eq.get("variable").is_none() {
            eq.set(Entry {
                name: "variable".to_string(),
                value: Value::Text { text: var, is_file: false },
                partag: None,
                distribute: false,
            });
        }
        if eq.get("procedure").is_none() {
            eq.set(Entry {
                name: "procedure".to_string(),
                value: Value::Text { text: proc.to_string(), is_file: false },
                partag: None,
                distribute: false,
            });
        }
    }
}

/// `Post File = case.vtu` gains an internal ResultOutputSolver
/// (`AddVtuOutputSolverHack`), and `Post File` is removed from Simulation.
fn add_vtu_output_solver(run: &mut Model) {
    let sim = match run.simulation.as_mut() {
        Some(s) => s,
        None => return,
    };
    let post = match sim.get("post file") {
        Some(Entry { value: Value::Text { text, .. }, .. }) => text.clone(),
        _ => return,
    };
    let Some(k) = post.find(".vtu") else { return };
    if find_solver_by_proc(&run.solvers, "resultoutputsolver").is_some() {
        return;
    }
    sim.remove("post file");
    let stem = post[..k].to_string();
    let prefix_hits: Vec<(String, Value)> = sim
        .entries
        .iter()
        .filter(|e| e.name.starts_with("vtu:"))
        .map(|e| {
            let mut n = e.name["vtu:".len()..].to_string();
            if n.starts_with(' ') {
                n.remove(0);
            }
            (n, e.value.clone())
        })
        .collect();
    let ascii = sim.get("ascii output").map(|e| e.value.clone());
    let binary = sim.get("binary output").map(|e| e.value.clone());
    let s = append_solver(run);
    set_entry(s, "procedure", Value::Text { text: "resultoutputsolve resultoutputsolver".to_string(), is_file: false });
    set_entry(s, "equation", Value::Text { text: "internalvtuoutputsolver".to_string(), is_file: false });
    set_entry(s, "output format", Value::Text { text: "vtu".to_string(), is_file: false });
    set_entry(s, "output file name", Value::Text { text: stem, is_file: false });
    set_entry(s, "no matrix", Value::Logical(true));
    set_entry(s, "variable", Value::Text { text: "-global vtu_internal_dummy".to_string(), is_file: false });
    set_entry(s, "exec solver", Value::Text { text: "after saving".to_string(), is_file: false });
    set_entry(s, "save geometry ids", Value::Logical(true));
    for (n, v) in prefix_hits {
        set_entry(s, &n, v);
    }
    if let Some(v) = ascii {
        if s.get("ascii output").is_none() {
            set_entry(s, "ascii output", v);
        }
    }
    if let Some(v) = binary {
        if s.get("binary output").is_none() {
            set_entry(s, "binary output", v);
        }
    }
}

fn set_entry(s: &mut crate::model::Section, name: &str, value: Value) {
    s.set(Entry { name: name.to_string(), value, partag: None, distribute: false });
}

/// `Scalars File = f.dat` gains an internal SaveScalars solver
/// (`AddSaveScalarsHack`).
fn add_save_scalars_solver(run: &mut Model) {
    let sim = match run.simulation.as_ref() {
        Some(s) => s,
        None => return,
    };
    let scalars = match sim.get("scalars file") {
        Some(Entry { value: Value::Text { text, .. }, .. }) => text.clone(),
        _ => return,
    };
    if find_solver_by_proc(&run.solvers, "savescalars").is_some() {
        return;
    }
    let prefix_hits: Vec<(String, Value)> = sim
        .entries
        .iter()
        .filter(|e| e.name.starts_with("scalars:"))
        .map(|e| {
            let mut n = e.name["scalars:".len()..].to_string();
            if n.starts_with(' ') {
                n.remove(0);
            }
            (n, e.value.clone())
        })
        .collect();
    let s = append_solver(run);
    set_entry(s, "procedure", Value::Text { text: "savedata savescalars".to_string(), is_file: false });
    set_entry(s, "equation", Value::Text { text: "internalsavescalars".to_string(), is_file: false });
    set_entry(s, "filename", Value::Text { text: scalars, is_file: false });
    set_entry(s, "exec solver", Value::Text { text: "after saving".to_string(), is_file: false });
    for (n, v) in prefix_hits {
        set_entry(s, &n, v);
    }
}

fn check_obsolete_transmittivity(run: &Model) -> Result<(), String> {
    for bc in run.bcs.iter().chain(run.materials.iter()) {
        if bc.get("transmittivity").is_some() {
            return Err(
                "obsolete keyword 'transmittivity': use 'transmissivity' instead".to_string(),
            );
        }
    }
    Ok(())
}
