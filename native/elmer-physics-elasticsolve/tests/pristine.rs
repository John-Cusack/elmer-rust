//! Quick tests that run `ElasticSolve`, from their pristine `.sif` and mesh,
//! against the nodal fields the pristine ElmerSolver wrote for the same case
//! (`pristine.tsv` beside each `.sif`: variable, L2 norm, max abs, node count of
//! the last `.result`, 17 digits). Meshes come from the pristine ElmerGrid.
//!
//! The per-case checks are the L2 norm and the largest absolute value of every
//! field this solver produces, which are what `Reference Norm` and
//! `SaveScalars ... max abs` read in the quick tests.

mod support;

use std::collections::BTreeMap;

use support::{data, Run};

/// `name -> (l2 norm, max abs)`.
fn pristine_file(case: &str, file: &str) -> BTreeMap<String, (f64, f64)> {
    std::fs::read_to_string(data(case).join(file))
        .unwrap()
        .lines()
        .map(|l| {
            let t: Vec<&str> = l.split('\t').collect();
            (t[0].to_string(), (t[1].parse().unwrap(), t[2].parse().unwrap()))
        })
        .collect()
}

fn norms(run: &Run, name: &str) -> Option<(f64, f64)> {
    let s = run.sim.vars.strided(name)?;
    let n = s.len() as f64;
    let l2 = (s.iter().map(|v| v * v).sum::<f64>() / n).sqrt();
    let mx = s.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    Some((l2, mx))
}

fn close(a: f64, b: f64, rel: f64, abs: f64) -> bool {
    (a - b).abs() <= rel * b.abs() + abs
}

/// Checks every variable of the pristine file that this run has; `must` have to be there.
fn check(case: &str, run: &Run, must: &[&str], rel: f64) {
    check_file(case, "pristine.tsv", run, must, rel);
}

fn check_file(case: &str, file: &str, run: &Run, must: &[&str], rel: f64) {
    let reference = pristine_file(case, file);
    for m in must {
        assert!(reference.contains_key(*m), "{m} not in pristine.tsv");
    }
    for (name, (l2, mx)) in &reference {
        let Some((a2, amx)) = norms(run, name) else {
            assert!(!must.contains(&name.as_str()), "{case}: no variable {name}");
            continue;
        };
        println!("{case} {name}: l2 {a2:.15e} (pristine {l2:.15e}) max {amx:.15e} (pristine {mx:.15e})");
        let scale = l2.abs().max(mx.abs());
        assert!(close(a2, *l2, rel, 1e-12 * scale), "{case} {name}: l2 {a2:e} vs {l2:e}");
        assert!(close(amx, *mx, rel, 1e-12 * scale), "{case} {name}: max {amx:e} vs {mx:e}");
    }
}

#[test]
fn elastic_cylindric() {
    let run = support::run_sif("elasticcylindric", "case.sif", |_| {}).unwrap();
    check("elasticcylindric", &run, &["displacement 1", "displacement 2", "stress_xx", "strain_zz", "vonmises"], 1e-8);
    assert_eq!(run.passed, Some(true));
}

#[test]
fn elastic_stress_load_plane() {
    let run = support::run_sif("elasticstressloadplane", "case.sif", |_| {}).unwrap();
    check("elasticstressloadplane", &run, &["displacement 1", "stress_xy", "stress_xx", "strain_zz"], 1e-8);
}

#[test]
fn elastic_anisotropic_axi() {
    let run = support::run_sif("elasticanisotropicaxi", "case.sif", |_| {}).unwrap();
    check("elasticanisotropicaxi", &run, &["displacement 1", "displacement 2", "stress_yy", "stress_zz"], 1e-8);
}

#[test]
fn elastic_thermal_bimetal() {
    // `Exported Variable 1 = Temperature` with `Initial Condition 1: Temperature = 400`.
    let run = support::run_sif("elasticthermalbimetal", "case.sif", |_| {}).unwrap();
    check("elasticthermalbimetal", &run, &["displacement 1", "displacement 2", "temperature"], 1e-8);
    assert_eq!(run.passed, Some(true));
}

#[test]
fn neohookean2d_fields() {
    let run = support::run_sif("neohookean2d", "cooksmembrane.sif", |_| {}).unwrap();
    check("neohookean2d", &run, &["displacement 1", "displacement 2"], 1e-9);
}
