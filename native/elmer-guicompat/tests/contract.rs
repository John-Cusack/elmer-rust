//! Fixture, typing, loader and checker tests that need no real solver.

use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::time::Duration;

use elmer_guicompat::check::{run_all, Config, Status};
use elmer_guicompat::fixtures::{build_all, Class, EdfSet, EDF_FILES};
use elmer_guicompat::logparse::parse_log;
use elmer_guicompat::sif::scan;
use elmer_guicompat::typecheck::check_fixture;
use elmer_guicompat::vtu::{load_ep, load_vtu};

/// `ELMER_GUICOMPAT_EDF_DIR` overrides the default fork checkout.
fn edf_dir() -> String {
    std::env::var("ELMER_GUICOMPAT_EDF_DIR")
        .unwrap_or_else(|_| "/home/john/runs/elmer-fork/ElmerGUI/Application/edf".into())
}

fn edfs() -> Option<EdfSet> {
    let dir = edf_dir();
    if !Path::new(&dir).is_dir() {
        eprintln!("skipping: {dir} not present");
        return None;
    }
    Some(EdfSet::load(Path::new(&dir)).expect("EDF files parse"))
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("elmer-guicompat-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn seven_fixtures_one_per_edf() {
    let Some(set) = edfs() else { return };
    let all = build_all(&set);
    assert_eq!(all.len(), 7);
    assert_eq!(EDF_FILES.len(), 7);
    let classes: Vec<_> = all.iter().map(|f| (f.id.as_str(), f.class)).collect();
    assert!(classes.contains(&("heatequation", Class::Heat)));
    assert!(classes.contains(&("linearelasticity", Class::Elasticity)));
    assert!(classes.contains(&("resultoutput", Class::Output)));
    assert_eq!(
        all.iter().filter(|f| f.class == Class::Unsupported).count(),
        4
    );
    for f in &all {
        let names: Vec<&str> = f.files.iter().map(|(n, _)| n.as_str()).collect();
        assert!(names.contains(&"case.sif") && names.contains(&"ELMERSOLVER_STARTINFO"));
        assert!(names.contains(&"mesh.header"), "{}", f.id);
        assert_eq!(f.files[1].1, "case.sif\n1\n");
    }
}

#[test]
fn procedure_and_variable_come_from_the_edf() {
    let Some(set) = edfs() else { return };
    let get = |id: &str| build_all(&set).into_iter().find(|f| f.id == id).unwrap();
    let heat = get("heatequation");
    assert!(heat
        .sif()
        .contains("Procedure = \"HeatSolve\" \"HeatSolver\""));
    assert!(heat.sif().contains("Variable = Temperature"));
    assert!(get("navier-stokes")
        .sif()
        .contains("Variable = Flow Solution[Velocity:2 Pressure:1]"));
    assert!(get("meshdeform")
        .sif()
        .contains("Variable = -dofs 2 Mesh Update"));
    assert!(get("helmholtz")
        .sif()
        .contains("Variable = -dofs 2 Pressure Wave"));
    assert!(get("resultoutput")
        .sif()
        .contains("Procedure = \"ResultOutputSolve\" \"ResultOutputSolver\""));
}

#[test]
fn every_fixture_keyword_resolves_through_the_catalog() {
    let Some(set) = edfs() else { return };
    for f in build_all(&set) {
        let edf = set.get(if f.id == "resultoutput" {
            "heatequation"
        } else {
            &f.id
        });
        let errors: Vec<_> = check_fixture(&f, edf)
            .into_iter()
            .filter(|x| !x.note)
            .collect();
        assert!(errors.is_empty(), "{}: {errors:?}", f.id);
        assert!(!f.edf_keys.is_empty(), "{} uses no EDF keywords", f.id);
    }
}

/// The checkbox rule of `SifGenerator::handleCheckBox`: a box left at its
/// EDF default is not written, so the SIF must not carry such keywords.
#[test]
fn gui_default_checkboxes_are_not_written() {
    let Some(set) = edfs() else { return };
    let out = build_all(&set)
        .into_iter()
        .find(|f| f.id == "resultoutput")
        .unwrap();
    assert!(!out
        .sif()
        .to_ascii_lowercase()
        .contains("no matrix equation"));
    let keys: Vec<_> = scan(out.sif()).into_iter().map(|a| a.key).collect();
    assert!(keys.contains(&"Output Format".to_string()));
}

#[test]
fn untyped_keyword_outside_the_catalog_is_caught() {
    let Some(set) = edfs() else { return };
    let mut f = build_all(&set).remove(0);
    f.files[0].1 = f.files[0].1.replace(
        "Heat Conductivity = 1.0",
        "Heat Conductivity = 1.0\n  Bogus Keyword = 2",
    );
    let edf = set.get("heatequation");
    let errs: Vec<_> = check_fixture(&f, edf)
        .into_iter()
        .filter(|x| !x.note)
        .collect();
    assert_eq!(errs.len(), 1);
    assert!(errs[0].problem.contains("unlisted"));
    f.files[0].1 = f.files[0]
        .1
        .replace("Bogus Keyword = 2", "Heat Capacity = abc");
    let errs: Vec<_> = check_fixture(&f, edf)
        .into_iter()
        .filter(|x| !x.note)
        .collect();
    assert!(errs.iter().any(|e| e.problem.contains("does not fit")));
}

/// Files written by the pristine Fortran ElmerSolver for the heat and
/// result-output fixtures (v26.2, serial, 2x2 square).
#[test]
fn loads_real_solver_output() {
    let vtu = load_vtu(include_bytes!("data/elmer_heat_t0001.vtu")).unwrap();
    assert_eq!((vtu.points, vtu.cells), (9, 16));
    let t = vtu.point_array("temperature").unwrap();
    assert_eq!(t.tuples(), 9);
    assert!(t.values.iter().all(|v| v.is_finite()));
    // Dirichlet rows: y=0 is 0, y=1 is 1.
    for (i, xyz) in vtu.coordinates.chunks(3).enumerate() {
        if xyz[1] == 0.0 {
            assert!(t.values[i].abs() < 1e-9);
        }
        if xyz[1] == 1.0 {
            assert!((t.values[i] - 1.0).abs() < 1e-9);
        }
    }
    let ep = load_ep(include_str!("data/elmer_heat.ep")).unwrap();
    assert_eq!(
        (ep.nodes, ep.elements, ep.components, ep.timesteps),
        (9, 16, 1, 1)
    );
    assert_eq!(
        ep.fields,
        [("scalar".to_string(), "Temperature".to_string())]
    );
}

#[test]
fn parses_a_real_solver_log() {
    let log = parse_log(
        " Time: 1/3:   1.000E-01\n\
         ComputeChange: NS (ITER=1) (NRM,RELC): ( 0.68083894      2.0000000     ) :: heat equation\n\
         ComputeChange: SS (ITER=1) (NRM,RELC): ( 0.68083894      2.0000000     ) :: heat equation\n\
         MAIN: *** Elmer Solver: ALL DONE ***\n",
    );
    assert_eq!(log.times, [0.1]);
    assert_eq!(log.changes.len(), 2);
    assert!(log.finished);
}

fn fake_solver(dir: &Path, body: &str) -> PathBuf {
    let p = dir.join("ElmerSolver");
    std::fs::write(&p, format!("#!/bin/sh\n{body}\n")).unwrap();
    std::fs::set_permissions(&p, std::fs::Permissions::from_mode(0o755)).unwrap();
    p
}

fn config(dir: &Path, solver: PathBuf, only: &str) -> Config {
    Config {
        solver,
        solver_mpi: None,
        work: dir.join("work"),
        edf_dir: PathBuf::from(edf_dir()),
        timeout: Duration::from_secs(10),
        only: vec![only.to_string()],
        native_only: true,
    }
}

/// Test scaffolding for the checker itself: scripts that refuse, crash, or
/// claim a solve. They are never offered as solvers.
#[test]
fn checker_classifies_refusal_failure_and_fake_solves() {
    if edfs().is_none() {
        return;
    }
    let d = scratch("classify");
    let refuse = fake_solver(
        &d,
        "[ \"$1\" = -v ] && { echo ELMER SOLVER; exit 0; }\necho refused >&2; exit 78",
    );
    let rep = run_all(&config(&d, refuse, "navier-stokes"));
    let run = rep
        .steps
        .iter()
        .find(|s| s.name.starts_with("run:"))
        .unwrap();
    assert_eq!(run.status, Status::Unsupported);
    assert!(rep.ok(false), "{}", rep.to_text());
    assert!(
        rep.ok(true),
        "refusing unsupported physics is not a native gap"
    );

    let heat = run_all(&config(
        &d,
        fake_solver(&d, "[ \"$1\" = -v ] && exit 0\nexit 78"),
        "heatequation",
    ));
    assert!(heat.ok(false));
    assert!(
        !heat.ok(true),
        "refusing a native-target fixture fails --strict"
    );

    let crash = fake_solver(&d, "[ \"$1\" = -v ] && exit 0\nexit 3");
    let rep = run_all(&config(&d, crash, "heatequation"));
    assert!(!rep.ok(false));

    // Exit 0 without any solver output for physics nobody implements.
    let liar = fake_solver(&d, "[ \"$1\" = -v ] && exit 0\nexit 0");
    let rep = run_all(&config(&d, liar, "navier-stokes"));
    let run = rep
        .steps
        .iter()
        .find(|s| s.name.starts_with("run:"))
        .unwrap();
    assert_eq!(run.status, Status::Fail);
    assert!(run.detail.contains("no native implementation"));
}

#[test]
fn solver_that_cannot_start_fails_the_probe() {
    if edfs().is_none() {
        return;
    }
    let d = scratch("nostart");
    let rep = run_all(&config(&d, d.join("missing"), "heatequation"));
    assert_eq!(rep.steps[0].status, Status::Fail);
    assert!(!rep.ok(false));
}

/// Opt-in: `ELMER_GUICOMPAT_SOLVER=/abs/ElmerSolver` (plus `ELMER_HOME` and
/// `ELMER_LIB` for a Fortran build) runs the whole contract against it.
#[test]
fn real_solver_when_requested() {
    let Ok(solver) = std::env::var("ELMER_GUICOMPAT_SOLVER") else {
        eprintln!("skipping: ELMER_GUICOMPAT_SOLVER not set");
        return;
    };
    if edfs().is_none() {
        return;
    }
    let d = scratch("real");
    let cfg = Config {
        solver: PathBuf::from(solver),
        solver_mpi: std::env::var("ELMER_GUICOMPAT_SOLVER_MPI")
            .ok()
            .map(PathBuf::from),
        work: d.join("work"),
        edf_dir: PathBuf::from(edf_dir()),
        timeout: Duration::from_secs(120),
        only: Vec::new(),
        native_only: false,
    };
    let rep = run_all(&cfg);
    assert!(rep.ok(false), "{}", rep.to_text());
}
