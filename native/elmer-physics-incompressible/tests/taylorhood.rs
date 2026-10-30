//! The Taylor-Hood pair (`Element = "p:2"`, no `Pressure Stabilization`)
//! against the pristine ElmerSolver: `StokesStabilized/quad_taylorhood.sif`
//! and `tri_taylorhood.sif` with a direct solve (`tests/data/README.md`),
//! the p dofs from `elmer-pdofs`. Serendipity quads and p:2 triangles have
//! no bubbles: the velocity is quadratic, the edge pressures are pinned to
//! zero. Every slot's row, velocity and pressure of the `case.result`, and
//! the norm the run printed (velocity only).

mod support;

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_pdofs::PDofsBackend;
use elmer_physics_incompressible::{register, solver_defaults};
use elmer_solver::{Simulation, SolverRegistry};
use support::cases::{data, ints, model, ns_solver, profile, real, sec, st, steady};
use support::mesh::TestMesh;

/// The whitespace-split lines of a fixture, `#` comments skipped.
fn lines(path: &std::path::Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .filter(|l: &Vec<String>| !l.is_empty())
        .collect()
}

fn bc(tags: &[i32], entries: &[(&str, Value)]) -> ParamMap {
    let mut v = vec![("Target Boundaries", ints(tags))];
    v.extend(entries.iter().cloned());
    sec(&v)
}

fn run(case: &str, simulation: ParamMap) -> Result<Simulation, elmer_solver::SolverError> {
    let dir = data(case);
    let mesh = TestMesh::read(&dir, 2);
    let mut solver = ns_solver(&[
        ("Element", st("p:2")),
        ("Stokes Flow", Value::Logical(true)),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Nonlinear System Convergence Tolerance", real(1.0e-3)),
    ]);
    solver_defaults(&mut solver, 2, false);
    let bcs = vec![
        bc(&[1], &[("Velocity 1", profile(&mesh, 1, 2, |y| (2.0 - y) * (y - 1.0))), ("Velocity 2", real(0.0))]),
        bc(&[2], &[("Velocity 2", real(0.0))]),
        bc(&[3], &[("Velocity 1", real(0.0)), ("Velocity 2", real(0.0))]),
    ];
    let material = sec(&[("Density", real(1.0)), ("Viscosity", real(1.0e3))]);
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(model(simulation, solver, material, None, bcs));
    sim.backends.dofs = Some(Box::new(PDofsBackend::new()));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(&mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

fn run_case(case: &str, n_slots: usize, ref_norm: f64) {
    let dir = data(case);
    let sim = run(case, steady("Cartesian 2D")).unwrap();

    let u = sim.vars.get("flow solution").unwrap();
    assert_eq!(u.values.len(), 3 * n_slots, "{case}: dofs");
    let norm = u.conv.norm;

    // perm.dat: `slot row`; solution.dat: `slot u1 u2 p`, x(Perm(slot)).
    let perm = u.perm.as_deref().unwrap();
    let rows = lines(&dir.join("perm.dat"));
    assert_eq!(rows.len(), n_slots, "{case}: Perm");
    for l in &rows {
        let (slot, row): (usize, u32) = (l[0].parse().unwrap(), l[1].parse().unwrap());
        assert_eq!(perm[slot - 1] + 1, row, "{case}: Perm({slot})");
    }
    let want: Vec<Vec<f64>> = lines(&dir.join("solution.dat"))
        .iter()
        .map(|l| l[1..].iter().map(|x| x.parse().unwrap()).collect())
        .collect();
    assert_eq!(want.len(), n_slots);
    let mut scale = [0.0f64; 3];
    for w in &want {
        for c in 0..3 {
            scale[c] = scale[c].max(w[c].abs());
        }
    }
    let mut err = [0.0f64; 3];
    for (slot, w) in want.iter().enumerate() {
        let s = perm[slot] as usize * 3;
        for c in 0..3 {
            err[c] = err[c].max((u.values[s + c] - w[c]).abs());
        }
    }
    eprintln!("{case}: norm {norm:.17e}, max error / scale per component {:?}", [0, 1, 2].map(|c| err[c] / scale[c]));
    assert!((norm - ref_norm).abs() / ref_norm < 1e-10, "{case}: norm {norm:.17e} vs {ref_norm:.17e}");
    for c in 0..3 {
        assert!(err[c] <= 1e-10 * scale[c], "{case}: component {c} max error {} (scale {})", err[c], scale[c]);
    }
}

/// 531 serendipity quads: 595 nodes + 1125 edges.
#[test]
fn taylorhood_quad_matches_elmer() {
    run_case("taylorhood_quad", 1720, 0.04399609386373855);
}

/// 1062 triangles: 595 nodes + 1656 edges.
#[test]
fn taylorhood_tri_matches_elmer() {
    run_case("taylorhood_tri", 2251, 0.03851002030822979);
}

/// Full p:2 quads have a bubble each, which `LocalBulkMatrix` would
/// condense (`LCondensate`): refused.
#[test]
fn taylorhood_with_bubbles_is_refused() {
    let mut simulation = steady("Cartesian 2D");
    simulation.insert("Serendipity P Elements", Value::Logical(false));
    let err = run("taylorhood_quad", simulation).err().expect("refused");
    assert!(err.to_string().contains("bubbles"), "{err}");
}
