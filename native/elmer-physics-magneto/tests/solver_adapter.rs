//! `WhitneyAVSolver` and `WhitneyAVHarmonicSolver` through `elmer-solver`'s
//! `Simulation`, registered under the `.sif` procedure names, on the raw
//! mesh database (tags mapped by the model, as `elmer-run` passes them).

use std::path::PathBuf;
use std::sync::Mutex;

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_solver::{Model, Simulation, SolverRegistry};

static CWD: Mutex<()> = Mutex::new(());

fn numbered(v: &[elmer_sif::model::Section]) -> Vec<ParamMap> {
    let n = v.iter().map(|s| s.id as usize).max().unwrap_or(0);
    let mut out = vec![ParamMap::new(); n];
    for s in v.iter().filter(|s| s.id >= 1) {
        out[s.id as usize - 1] = ParamMap::from(s);
    }
    out
}

/// Runs the first `keep` solvers of `case/sif` (the others dropped) and
/// returns the last kept solver's variable norm, that variable's value
/// count and the reference tally.
fn run_n(case: &str, sif: &str, keep: usize, tweak: impl Fn(&mut Model)) -> (f64, usize, usize) {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(case);
    let parsed = {
        let _guard = CWD.lock().unwrap_or_else(|e| e.into_inner());
        std::env::set_current_dir(&dir).unwrap();
        let path = dir.join(sif);
        let mut parsed = elmer_sif::parse_file(path.to_str().unwrap()).unwrap();
        elmer_sif::complete(&mut parsed.sif, path.to_str().unwrap()).unwrap();
        parsed
    };
    let m = &parsed.sif.runs[0];
    let mut solvers = numbered(&m.solvers);
    solvers.truncate(keep);
    let mut model = Model {
        header: ParamMap::new(),
        simulation: m.simulation.as_ref().map(ParamMap::from).unwrap_or_default(),
        constants: m.constants.as_ref().map(ParamMap::from).unwrap_or_default(),
        solvers,
        equations: numbered(&m.equations),
        bodies: numbered(&m.bodies),
        materials: numbered(&m.materials),
        body_forces: numbered(&m.body_forces),
        initial_conditions: numbered(&m.ics),
        boundary_conditions: numbered(&m.bcs),
    };
    tweak(&mut model);
    let (db, name) = m.header.mesh_db.clone().unwrap();
    let mesh = elmer_mesh::load(dir.join(db).join(name).to_str().unwrap()).unwrap();
    let mut registry = SolverRegistry::new();
    elmer_physics_magneto::register(&mut registry);
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &registry, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();
    let id = sim.slot(keep).unwrap().variable.unwrap();
    let var = sim.vars.by_id(id);
    (var.conv.norm, var.values.len(), sim.reference.pass_count)
}

fn run(case: &str, sif: &str, tweak: impl Fn(&mut Model)) -> (f64, usize, usize) {
    run_n(case, sif, 1, tweak)
}

#[test]
fn steady_solver_through_the_simulation() {
    let (norm, n, passed) = run("mgdyn_airgap", "airgap.sif", |_| {});
    // 265 nodes + 1202 edges; norm over the edges only (Norm Permutation).
    assert_eq!(n, 265 + 1202);
    let reference = 7.08538965e-05;
    assert!((norm - reference).abs() / reference < 1e-5, "norm {norm}");
    assert_eq!(passed, 1, "CompareToReferenceSolution for Solver 1");
}

#[test]
fn harmonic_solver_through_the_simulation() {
    // elmer-solver's variable spec parser takes one-word component names
    // only, so `P[P re:1 P im:1]` becomes `P[Pre:1 Pim:1]` with the
    // condition keywords renamed to match.
    let (norm, n, _) = run("mgdyn_airgap_harmonic", "airgap.sif", |m| {
        let s = &mut m.solvers[0];
        s.insert("Linear System Convergence Tolerance", Value::real(1e-10));
        s.insert("Variable", Value::Str("P[Pre:1 Pim:1]".into()));
        for sec in m.boundary_conditions.iter_mut().chain(m.initial_conditions.iter_mut()) {
            let renamed = sec
                .0
                .iter()
                .map(|(k, v)| {
                    let k = match (k.strip_prefix("p re"), k.strip_prefix("p im")) {
                        (Some(rest), _) => format!("pre{rest}"),
                        (_, Some(rest)) => format!("pim{rest}"),
                        _ => k.clone(),
                    };
                    (k, v.clone())
                })
                .collect();
            sec.0 = renamed;
        }
    });
    assert_eq!(n, 2 * (265 + 1202));
    // The pristine run at the same tolerance (tests/data/README.md).
    let reference = 2.42424125e-03;
    assert!((norm - reference).abs() / reference < 1e-5, "norm {norm}");
}

#[test]
fn calc_fields_through_the_simulation() {
    let (norm, n, passed) = run_n("mgdyn_airgap", "airgap.sif", 2, |_| {});
    // hr_dummy over the 265 nodes holds Magnetic Field Strength 3.
    assert_eq!(n, 265);
    let reference = 3.06334400e3;
    assert!((norm - reference).abs() / reference < 1e-5, "norm {norm}");
    assert_eq!(passed, 2, "CompareToReferenceSolution for Solvers 1 and 2");
}

#[test]
fn mgdyn2d_harmonic_induction_heating3() {
    // fem/tests/InductionHeating3, Solver 1 (MagnetoDynamics2DHarmonic,
    // axisymmetric). The section's `Variable = Potential[Potential Re:1
    // Potential Im:1]` is the solver's default, which `init` creates; it is
    // dropped here because elmer-solver's variable-spec parser takes
    // one-word component names only.
    let (norm, n, passed) = run("InductionHeating3", "crucible.sif", |m| {
        m.solvers[0].0.remove("variable");
    });
    assert_eq!(n, 2 * 5194);
    // Pristine: Solver 1 Norm = 8.88220240E-06 (RefNorm 0.88822024E-05).
    let reference = 8.88220240e-06;
    eprintln!("InductionHeating3: norm {norm:.9e}, rel err {:.3e}", (norm - reference).abs() / reference);
    assert!((norm - reference).abs() / reference < 1e-5, "norm {norm}");
    assert_eq!(passed, 1, "CompareToReferenceSolution for Solver 1");
}

#[test]
fn mgdyn2d_harmonic_induction_heating4_surface_current() {
    // fem/tests/InductionHeating4: InductionHeating3 with the coil current
    // given as `Surface Current` on the coil boundaries (LocalMatrixBC)
    // instead of a body current density. `Variable` dropped as above.
    let (norm, n, passed) = run("InductionHeating4", "crucible.sif", |m| {
        m.solvers[0].0.remove("variable");
    });
    assert!(n > 0);
    // Pristine: Solver 1 Norm = 8.85766857E-06 (RefNorm 8.85766857E-06).
    let reference = 8.85766857e-06;
    eprintln!("InductionHeating4: norm {norm:.9e}, rel err {:.3e}", (norm - reference).abs() / reference);
    assert!((norm - reference).abs() / reference < 1e-5, "norm {norm}");
    assert_eq!(passed, 1, "CompareToReferenceSolution for Solver 1");
}
