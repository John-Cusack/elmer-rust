//! StatElecSolve on p elements through `Simulation` with `elmer-pdofs`'
//! `PDofsBackend`, against pristine runs of `StatElecSolveVec` with
//! `Element = p:2` on CapacitanceMatrix3's `plates.grd` (coarser, see
//! `tests/data/statelec_p2/*.sif`): Dirichlet potentials on three
//! boundaries and an electric flux on the fourth, with the bubbles in the
//! global system (`global.sif`, 780 dofs) and condensed locally
//! (`local.sif`, 620 dofs, `CondensateP`).

mod support;

use elmer_core::Value;
use elmer_solver::Model;
use support::cases::*;
use support::meshio::read_elmer_mesh;
use support::simmesh::SimMesh;

fn model(global_bubbles: bool) -> Model {
    let mut solver = vec![
        ("Equation", st("Stat Elec Solver")),
        ("Variable", st("Potential")),
        ("Procedure", st("\"StatElecSolveVec\" \"StatElecSolver\"")),
        ("Linear System Solver", st("Direct")),
        ("Linear System Direct Method", st("UMFPack")),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Element", st("p:2")),
    ];
    if !global_bubbles {
        solver.push(("Bubbles in Global System", Value::Logical(false)));
    }
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian")),
            ("Simulation Type", st("Steady State")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Serendipity P Elements", Value::Logical(false)),
        ]),
        constants: sec(&[("Permittivity Of Vacuum", real(1.0))]),
        solvers: vec![sec(&solver)],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[("Equation", ints(&[1])), ("Material", ints(&[1]))])],
        materials: vec![sec(&[("Relative Permittivity", real(1.0))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", ints(&[1])), ("Potential", real(0.0))]),
            sec(&[("Target Boundaries", ints(&[2])), ("Potential", real(1.0))]),
            sec(&[("Target Boundaries", ints(&[3])), ("Potential", real(2.0))]),
            sec(&[
                ("Target Boundaries", ints(&[4])),
                ("Electric Flux", real(0.5)),
            ]),
        ],
        ..Model::default()
    }
}

/// Rows of a fixture file, `#` headers skipped.
fn rows(name: &str) -> Vec<Vec<f64>> {
    let path = data("statelec_p2").join(name);
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|t| t.parse().unwrap()).collect())
        .collect()
}

fn check(global_bubbles: bool, case: &str, n_dofs: usize, norm: f64) {
    let mesh = SimMesh::new(read_elmer_mesh(&data("statelec_p2"), 2));
    let (mut sim, reg) = simulation(model(global_bubbles));
    sim.backends.dofs = Some(Box::new(elmer_pdofs::PDofsBackend::new()));
    sim.setup(
        &mesh,
        &reg,
        &elmer_core::NullLogger,
        &elmer_core::SerialComm,
    )
    .unwrap();
    sim.run(&mesh, &elmer_core::NullLogger, &elmer_core::SerialComm)
        .unwrap();
    let v = sim.vars.get("potential").unwrap();
    assert_eq!(v.values.len(), n_dofs, "{case}");
    let perm = v.perm.as_deref().unwrap();
    let want_perm = rows(&format!("{case}_perm.dat"));
    let want = rows(&format!("{case}_solution.dat"));
    let scale = want.iter().fold(0.0f64, |m, r| m.max(r[1].abs()));
    for (pr, sr) in want_perm.iter().zip(&want) {
        let slot = pr[0] as usize;
        assert_eq!(
            f64::from(perm[slot - 1] + 1),
            pr[1],
            "{case}: Perm of slot {slot}"
        );
        let ours = v.values[perm[slot - 1] as usize];
        assert!(
            (ours - sr[1]).abs() <= 1e-10 * scale,
            "{case}: slot {slot}: {ours:e} vs {:e}",
            sr[1]
        );
    }
    assert!(
        (v.conv.norm - norm).abs() < 1e-12,
        "{case}: norm {}",
        v.conv.norm
    );
}

#[test]
fn stat_elec_p2_global_bubbles_matches_pristine() {
    check(true, "global", 780, 0.471068247747862);
}

#[test]
fn stat_elec_p2_local_bubbles_matches_pristine() {
    check(false, "local", 620, 0.5283447349145732);
}
