//! The `PoissonAssembly` cases end to end against the pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) runs of `fem/tests/CoupledPoisson1`,
//! `3`..`9` with the compiled `PoissonAssembly.so`: the cases' own
//! `case.sif`, meshes from `ElmerGrid 1 2 angle|squares`. Each must pass
//! its reference check and print the norm the pristine run printed. Also
//! `fem/tests/tresca` (StressSolve, then `Tresca.so`'s `TrescaSolver`).

mod case;

/// `(case, mesh, pristine "Norm =" of solver 1)`.
const CASES: [(&str, &str, &str); 8] = [
    ("CoupledPoisson1", "angle", "7.78240309E-02"),
    ("CoupledPoisson3", "angle", "6.96052509E-02"),
    ("CoupledPoisson4", "angle", "5.10670175E+01"),
    ("CoupledPoisson5", "angle", "5.85178077E+01"),
    ("CoupledPoisson6", "angle", "7.56321899E-01"),
    ("CoupledPoisson7", "angle", "7.84091700E-02"),
    ("CoupledPoisson8", "squares", "5.25188761E-01"),
    ("CoupledPoisson9", "squares", "2.40813064E+00"),
];

#[test]
fn poisson_assembly_cases_match_pristine() {
    for (name, mesh, pristine) in CASES {
        let dir = case::stage(name, Some(mesh));
        let out = case::run(name, &dir, "case.sif").unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(out.passed, Some(true), "{name}: reference check\n{}", out.log);
        case::assert_norm_prints(name, out.norm(1), pristine);
        // Once per run, as in the pristine log: blocks naming one routine share its instance.
        assert_eq!(out.log.matches("PoissonBulkAssembly: 1st time").count(), 1, "{name}");
        assert_eq!(std::fs::read_to_string(dir.join("TEST.PASSED")).unwrap().trim(), "1", "{name}");
    }
}

#[test]
fn tresca_matches_pristine() {
    let dir = case::stage("tresca", None);
    let out = case::run("tresca", &dir, "Tresca.sif").unwrap_or_else(|e| panic!("tresca: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    case::assert_norm_prints("tresca", out.norm(1), "2.19639953E+01");
}

/// `fem/tests/ModelPDEhandle` (`ModelPDEhandle.so`, `ElmerGrid 1 2 square`):
/// convection-diffusion with Dirichlet, Robin and flux boundaries.
#[test]
fn model_pde_handle_matches_pristine() {
    let dir = case::stage("ModelPDEhandle", None);
    let out = case::run("ModelPDEhandle", &dir, "case.sif").unwrap_or_else(|e| panic!("ModelPDEhandle: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    case::assert_norm_prints("ModelPDEhandle", out.norm(1), "1.79722562E+00");
}

/// `fem/tests/StokesProj` (`Pressure.so`, `Velocity.so`): the projection
/// scheme converges in 45 coupled iterations, as the pristine run does.
#[test]
fn stokes_proj_matches_pristine() {
    let dir = case::stage("StokesProj", None);
    let out = case::run("StokesProj", &dir, "ns.sif").unwrap_or_else(|e| panic!("StokesProj: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    case::assert_norm_prints("StokesProj", out.norm(2), "3.79706752E-01");
    assert_eq!(out.log.matches("SS (ITER=").count(), 2 * 45, "45 coupled iterations of both solvers");
}

/// `fem/tests/bentonite` (`BentoniteSolver.so`, axisymmetric 1D mesh):
/// 29 nonlinear iterations, as the pristine run.
#[test]
fn bentonite_matches_pristine() {
    let dir = case::stage("bentonite", None);
    let out = case::run("bentonite", &dir, "Bentonite.sif").unwrap_or_else(|e| panic!("bentonite: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    case::assert_norm_prints("bentonite", out.norm(1), "1.95292467E+02");
    assert_eq!(out.log.matches("Bentonite iteration").count(), 29);
}
