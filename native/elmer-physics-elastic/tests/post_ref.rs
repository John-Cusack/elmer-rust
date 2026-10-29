//! Stress and strain postprocessing against `SaveScalars` of pristine
//! ElmerSolver (`tests/data/<case>/post_ref.txt`, from `<case>_post.sif`):
//! the sum over the nodes and the largest-magnitude value of every projected
//! component. Assembly and the direct solve are test-local stand-ins.

mod support;

use elmer_solver::Model;
use support::cases::{a2d, c3d, e2d, run};
use support::sim::Mesh;

/// `(variable, component)` of a `SaveScalars` name such as `stress_xy`.
fn locate(name: &str) -> (&'static str, usize) {
    let comp = |s: &str| ["xx", "yy", "zz", "xy", "yz", "xz"].iter().position(|c| *c == s).unwrap();
    match name {
        "vonmises" => ("vonmises", 0),
        "tresca" => ("tresca", 0),
        "principal stress 1" => ("principal stress", 0),
        "principal stress 3" => ("principal stress", 2),
        "principal strain 1" => ("principal strain", 0),
        n if n.starts_with("stress_") => ("stress", comp(&n[7..])),
        n if n.starts_with("strain_") => ("strain", comp(&n[7..])),
        n => panic!("unknown {n}"),
    }
}

fn check(case: &str, model: Model, mesh: Mesh, tol: f64) {
    let sim = run(model, &mesh);
    let reference = std::fs::read_to_string(
        [env!("CARGO_MANIFEST_DIR"), "tests", "data", case, "post_ref.txt"].iter().collect::<std::path::PathBuf>(),
    )
    .unwrap();
    let mut checked = 0;
    let mut bad = Vec::new();
    for line in reference.lines().filter(|l| !l.starts_with('#')) {
        let (key, val) = line.split_once(" = ").unwrap();
        let (op, name) = key.split_once(": ").unwrap();
        let want: f64 = val.trim().parse().unwrap();
        let (var, comp) = locate(name);
        let v = sim.vars.get(var).unwrap_or_else(|| panic!("{case}: no variable {var}"));
        let xs: Vec<f64> = v.values.iter().skip(comp).step_by(v.dofs).copied().collect();
        let got = match op {
            "sum" => xs.iter().sum::<f64>(),
            "max abs" => xs.iter().copied().max_by(|a, b| a.abs().total_cmp(&b.abs())).unwrap(),
            o => panic!("op {o}"),
        };
        // Cancellation in a sum leaves an absolute floor set by the field's size.
        let scale = xs.iter().fold(0.0_f64, |m, x| m.max(x.abs())).max(want.abs() * 1e-3);
        if (got - want).abs() > tol * (want.abs() + scale) {
            bad.push(format!("{case} {key}: got {got:e}, pristine {want:e}"));
        }
        checked += 1;
    }
    assert!(bad.is_empty(), "{}", bad.join("\n"));
    assert!(checked >= 20, "{case}: only {checked} reference values");
}

#[test]
fn plane_strain_2d_matches_pristine() {
    let (m, mesh) = a2d();
    check("a2d", m, mesh, 1e-8);
}

#[test]
fn axisymmetric_matches_pristine() {
    let (m, mesh) = e2d();
    check("e2d", m, mesh, 1e-8);
}

#[test]
fn anisotropic_3d_matches_pristine() {
    let (m, mesh) = c3d();
    check("c3d", m, mesh, 1e-8);
}
