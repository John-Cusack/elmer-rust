//! Periodic cases end to end: the pristine `.sif` through `elmer-sif`, the
//! mesh through `elmer-mesh`, the projectors from this crate and a
//! `Simulation` with [`PeriodicAssembly`] as its assembly backend, compared
//! with the pristine solver norm and `Output File` (the bandwidth-optimised
//! `Perm`, which includes the projector couplings, and the nodal values).
//!
//! `elmer-run`'s `run_case` refuses `Periodic BC`, so the steps it takes
//! for these cases are repeated here.

use std::path::{Path, PathBuf};
use std::rc::Rc;

use elmer_core::{MeshView, NullLogger, Params, SerialComm, Value};
use elmer_mesh::adapters::MeshViewAdapter;
use elmer_mortar::{generate_periodic_projectors, PeriodicAssembly, ProjectorEnv};
use elmer_run::{coordinate_system, install_fe, solver_registry, SifMatc};
use elmer_solver::Simulation;

fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for e in std::fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let dst = to.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_dir(&e.path(), &dst);
        } else {
            std::fs::copy(e.path(), dst).unwrap();
        }
    }
}

/// A scratch copy of `tests/fixtures/<case>` without the expected files.
fn scratch(case: &str) -> PathBuf {
    let src = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(case);
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("elmer-mortar-e2e")
        .join(case);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::copy(src.join("case.sif"), dir.join("case.sif")).unwrap();
    copy_dir(&src.join("mesh"), &dir.join("mesh"));
    dir
}

/// Solver 1's variable after the run.
struct Solution {
    norm: f64,
    /// Node -> 0-based slot.
    perm: Vec<u32>,
    /// Nodal values in node order.
    values: Vec<f64>,
}

/// The single-field, single-timestep pristine `Output File`: the 1-based
/// `Perm` pairs, then one value per node.
fn read_result(path: &Path) -> (Vec<u32>, Vec<f64>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines().skip_while(|l| !l.starts_with("Perm:"));
    let n: usize = lines
        .next()
        .unwrap()
        .split_whitespace()
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    let mut perm = vec![u32::MAX; n];
    for l in lines.by_ref().take(n) {
        let mut f = l.split_whitespace().map(|x| x.parse::<usize>().unwrap());
        let (node, slot) = (f.next().unwrap(), f.next().unwrap());
        perm[node - 1] = slot as u32 - 1;
    }
    let values = lines.take(n).map(|l| l.trim().parse().unwrap()).collect();
    (perm, values)
}

/// Runs `tests/fixtures/<case>`.
fn run(case: &str) -> Solution {
    let dir = scratch(case);
    let sif = dir.join("case.sif");
    let path = sif.to_str().unwrap();
    let mut parsed = elmer_sif::parse_file(path).unwrap();
    elmer_sif::complete(&mut parsed.sif, path).unwrap();
    let model = elmer_run::sif::to_model(&parsed.sif.runs[0]).unwrap();
    let mesh = elmer_mesh::load(dir.join("mesh").to_str().unwrap()).unwrap();
    let view = MeshViewAdapter::new(&mesh);
    let coords = coordinate_system(&model, &view).unwrap();
    install_fe(coords.dim).unwrap();

    let bc_of = |tag: i32| {
        model
            .boundary_conditions
            .iter()
            .position(|bc| matches!(bc.get("target boundaries"), Some(Value::Integer(v)) if v.contains(&tag)))
            .map(|i| i + 1)
    };
    let bcs: Vec<&dyn Params> = model
        .boundary_conditions
        .iter()
        .map(|b| b as &dyn Params)
        .collect();
    let env = ProjectorEnv {
        mesh: &view,
        bcs: &bcs,
        bc_of: &bc_of,
        simulation: &model.simulation,
        solver: model.solvers.first().map(|s| s as &dyn Params),
        coord_dim: coords.dim,
        log: &NullLogger,
    };
    let (projectors, stop) = generate_periodic_projectors(&env, &dir, 0).unwrap();
    assert!(!stop);
    let backend = PeriodicAssembly::new(
        model.boundary_conditions.clone(),
        projectors,
        view.n_nodes(),
    );

    let mut sim = Simulation::new(model.clone());
    sim.set_matc(Rc::new(SifMatc));
    sim.backends.assembly = Some(Box::new(backend));
    sim.setup(&view, &solver_registry(), &NullLogger, &SerialComm)
        .unwrap();
    sim.run(&view, &NullLogger, &SerialComm).unwrap();
    let slot = sim.slots.iter().find(|s| s.id == 1).unwrap();
    let v = sim.vars.by_id(slot.variable.unwrap());
    let perm = v.perm.as_deref().unwrap().to_vec();
    let values = perm.iter().map(|&p| v.values[p as usize]).collect();
    Solution {
        norm: v.conv.norm,
        perm,
        values,
    }
}

/// Solver 1's norm against the 8 significant digits the pristine run
/// prints, its permutation exactly and its nodal values to the linear
/// solver's tolerance, relative to the largest value.
fn check(case: &str, pristine: f64) {
    let s = run(case);
    let rel = (s.norm - pristine).abs() / pristine;
    println!(
        "{case}: native norm {:.17e}, pristine {pristine:e}, rel {rel:e}",
        s.norm
    );
    assert!(
        rel <= 1e-8,
        "{case}: norm {}, pristine {pristine}, rel {rel:e}",
        s.norm
    );

    let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(case);
    let (perm, values) = read_result(&dir.join("expected/case.result"));
    assert_eq!(s.perm, perm, "{case}: Perm");
    let scale = values.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    let err = s
        .values
        .iter()
        .zip(&values)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
    println!("{case}: max nodal difference {err:e} (max |u| {scale:e})");
    assert!(
        err <= 1e-8 * scale,
        "{case}: max nodal difference {err:e}, max |u| {scale:e}"
    );
}

#[test]
fn periodic1_two_pairs_norm() {
    check("periodic1", 1.33502443e-2);
}

#[test]
fn periodic2_default_overlay_norm() {
    check("periodic2", 2.43762993);
}

#[test]
fn periodic_offset_norm() {
    check("periodic_offset", 2.48409204);
}

#[test]
fn periodic_nonconforming_norm() {
    check("periodic_nonconforming", 3.75473557e-2);
}
