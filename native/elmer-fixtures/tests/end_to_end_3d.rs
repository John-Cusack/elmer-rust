//! `Poisson` `PoissonSolver` (the `ExtrudedMesh/Poisson.f90` source) on an
//! extruded mesh, end to end against the pristine ElmerSolver (26.2-devel,
//! Rev 9f6af2f85) run of `fem/tests/ExtrudedMesh` (`ElmerGrid 1 2 square`,
//! `Extruded Mesh Layers = 3`).

mod case;

#[test]
fn extruded_mesh_matches_pristine() {
    let dir = case::stage("ExtrudedMesh", None);
    let out = case::run("ExtrudedMesh", &dir, "ptest.sif").unwrap_or_else(|e| panic!("ExtrudedMesh: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    case::assert_norm_prints("ExtrudedMesh", out.norm(1), "3.78170795E-02");
}

/// `ExtrudedMeshBaseline1` runs `PoissonSolver` on the baseline, boundary
/// elements that `Body ID = 2` makes a body. `elmer-solver` has no such
/// active elements, so the solver refuses instead of solving an empty
/// system (the pristine norm is 3.35220719E+01).
#[test]
fn extruded_baseline_poisson_is_refused_by_name() {
    let dir = case::stage("ExtrudedMeshBaseline1", None);
    let err = case::run("ExtrudedMeshBaseline1", &dir, "case_baseline.sif").err().expect("refused");
    assert!(matches!(err, elmer_solver::SolverError::Unsupported(_)), "{err}");
    assert!(err.to_string().contains("Boundary Condition 5 Body ID"), "{err}");
    assert!(!dir.join("TEST.PASSED").exists(), "refused before any output");
}

/// `ElementUnitTest_h1`: twelve `RUN` blocks, one element per visit
/// (808, 706, 605, 504, then the boundary 404 and 303), first the element
/// measure `SUM(MASS)`, then `SUM(ABS(STIFF))`, each against the `.sif`'s
/// reference value (the pristine run with `SingleElements.so` passes all 12).
#[test]
fn element_unit_test_h1_passes_all_twelve_runs() {
    let dir = case::stage("ElementUnitTest_h1", None);
    let out = case::run("ElementUnitTest_h1", &dir, "case.sif").unwrap_or_else(|e| panic!("ElementUnitTest_h1: {e}"));
    assert!(out.log.contains("PASSED all 12 tests"), "{}", out.log);
    assert_eq!(out.passed, Some(true));
}
