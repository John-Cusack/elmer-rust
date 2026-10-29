//! `Linear System Complex`: the interleaved real system HelmholtzSolve
//! assembles, solved in complex arithmetic. Checked against pristine
//! ElmerSolver (26.2-devel, Rev 9f6af2f85) runs of HelmholtzSolve on the
//! 1dtests mesh (k = 2 with damping 0.3, complex Dirichlet values at both
//! ends; `tests/data/complex/`): UMFPACK, and complex BiCGStab with the
//! diagonal preconditioner (112 iterations to 1e-10 in Elmer). HelmholtzSolve
//! is the test-local stand-in; its `_init` sets `Linear System Complex`,
//! which the test adds to the section.

mod support;

use std::path::Path;

use elmer_core::{SerialComm, Value};
use elmer_solver::{Simulation, SolverError, SolverRegistry};
use support::poisson::helmholtz;
use support::{sif_model, split_change, Capture};

fn run(sif: &str, extra: &[(&str, Value)]) -> Result<(Simulation, Capture), SolverError> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join("1d"), 1);
    let mut model = sif_model(&dir.join("complex").join(sif));
    for (k, v) in extra {
        model.solvers[0].insert(k, v.clone());
    }
    let mut reg = SolverRegistry::new();
    reg.register("HelmholtzSolve", "HelmholtzSolver", helmholtz);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm)?;
    sim.run(&mesh, &log, &SerialComm)?;
    Ok((sim, log))
}

fn check(sif: &str, stdout: &str) -> Vec<f64> {
    let (sim, log) = run(sif, &[("Linear System Complex", Value::Logical(true))]).unwrap();
    let got = log.lines(&["ComputeChange"]);
    let got: Vec<&String> = got.iter().filter(|l| l.contains("(NRM,RELC)")).collect();
    let want: Vec<&str> = stdout.lines().collect();
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        assert!((gr - wr).abs() <= 1e-9 * wr, "{g}\n{w}");
    }
    sim.vars.get("pressure").unwrap().values.clone()
}

#[test]
fn helmholtz_direct_complex() {
    let x = check("direct.sif", include_str!("data/complex/direct_stdout.txt"));
    // The same system solved as the real interleaved one.
    let (real, _) = run("direct.sif", &[]).unwrap();
    for (a, b) in x.iter().zip(&real.vars.get("pressure").unwrap().values) {
        assert!((a - b).abs() <= 1e-12, "{a} vs {b}");
    }
}

#[test]
fn helmholtz_complex_bicgstab_diagonal() {
    let x = check("iterative.sif", include_str!("data/complex/iterative_stdout.txt"));
    let (direct, _) = run("direct.sif", &[("Linear System Complex", Value::Logical(true))]).unwrap();
    for (a, b) in x.iter().zip(&direct.vars.get("pressure").unwrap().values) {
        assert!((a - b).abs() <= 1e-8, "{a} vs {b}");
    }
}

#[test]
fn complex_dispatch_refuses_what_is_not_ported() {
    // BiCGStab(2) runs on the real system but IterSolve remaps it to
    // BiCGStab(l) for a complex one, which is refused here.
    let m = [("Linear System Iterative Method", Value::Str("bicgstab2".into()))];
    assert!(run("iterative.sif", &m).is_ok());
    let mut c = m.to_vec();
    c.push(("Linear System Complex", Value::Logical(true)));
    let err = run("iterative.sif", &c).err().expect("refused");
    assert!(err.to_string().contains("bicgstab2"), "{err}");
    let err = run("direct.sif", &[
        ("Linear System Complex", Value::Logical(true)),
        ("Linear System Direct Method", Value::Str("cholmod".into())),
    ])
    .err()
    .expect("refused");
    assert!(err.to_string().contains("Cholesky"), "{err}");
}
