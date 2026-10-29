//! `Poisson` `PoissonSolver` (the `ExtrudedMesh/Poisson.f90` source) on an
//! extruded mesh, end to end against the pristine ElmerSolver (26.2-devel,
//! Rev 9f6af2f85) run of `fem/tests/ExtrudedMesh` (`ElmerGrid 1 2 square`,
//! `Extruded Mesh Layers = 3`).

mod case;

use std::path::Path;

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

/// `ElementUnitTest_edge`, `_edge_piola` and `_edge_quadratic`: the
/// Whitney, lowest-order Piola and second-order Piola edge bases on one
/// element per visit, twelve `RUN` blocks each against the `.sif`'s
/// reference values (the pristine runs with `SingleEdgeElements.so` pass
/// all 12; `elmer-fe-hp`'s `edge_unit_norms` test checks the 17-digit
/// values those runs print).
#[test]
fn element_unit_test_edge_cases_pass_all_twelve_runs() {
    for name in ["ElementUnitTest_edge", "ElementUnitTest_edge_piola", "ElementUnitTest_edge_quadratic"] {
        let dir = case::stage(name, None);
        let mesh = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/meshes/unitelems_edge");
        case::copy_dir(&mesh, &dir.join("unitelems"));
        let out = case::run(name, &dir, "case.sif").unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(out.log.contains("PASSED all 12 tests"), "{name}\n{}", out.log);
        assert_eq!(out.passed, Some(true), "{name}");
    }
}

/// `ElementUnitTest_p1` and `SD_ElementUnitTest_p1`: `SingleElements` on
/// `Element = "p:1"`, with each element's p basis and Gauss rule from
/// `elmer-pdofs`. Its p layout does not cover tetrahedra yet (the unit mesh
/// has a 504), so the solver refuses by name before the run.
#[test]
fn element_unit_test_p1_cases_are_refused_until_p_tetrahedra() {
    for name in ["ElementUnitTest_p1", "SD_ElementUnitTest_p1"] {
        let dir = case::stage(name, None);
        let err = case::run(name, &dir, "case.sif").err().unwrap_or_else(|| panic!("{name}: ran"));
        assert!(matches!(err, elmer_solver::SolverError::Unsupported(_)), "{name}: {err}");
        assert!(err.to_string().contains("tetrahedral p elements"), "{name}: {err}");
        assert!(!dir.join("TEST.PASSED").exists(), "{name}: refused before any output");
    }
}

/// The p path of `SingleElements` where `elmer-pdofs` has a layout: the
/// brick of ElementUnitTest_p1's unit mesh (`tests/data/ElementUnitTest_p1_brick`),
/// its mass and stiffness sums checked against the pristine run's 17-digit
/// values to 1e-15.
#[test]
fn element_unit_test_p1_brick_matches_pristine_digits() {
    let dir = case::stage("ElementUnitTest_p1_brick", None);
    let out = case::run("ElementUnitTest_p1", &dir, "case.sif").unwrap_or_else(|e| panic!("p1 brick: {e}"));
    assert!(out.log.contains("PASSED all 2 tests"), "{}", out.log);
}

/// `fem/tests/L2ndDerivatives` (`L2ndDerivativesTest.so`): no `Body` or
/// `Equation` sections, so `LoadModel`'s defaults activate the solver on all
/// 20 elements, each passing the basis and second-derivative checks; the
/// solver's nameless null variable has the reference norm 0, as in the
/// pristine run.
#[test]
fn l2nd_derivatives_matches_pristine() {
    let dir = case::stage("L2ndDerivatives", None);
    let out = case::run("L2ndDerivatives", &dir, "test.sif").unwrap_or_else(|e| panic!("L2ndDerivatives: {e}"));
    assert_eq!(out.passed, Some(true), "{}", out.log);
}
