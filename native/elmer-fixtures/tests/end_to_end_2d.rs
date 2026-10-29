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

/// `fem/tests/ModelPDEsplitted` and `SD_ModelPDEsplitted`
/// (`ModelPDEsplitted.so`, `ElmerGrid 1 2 square`): ModelPDE as `_bulk` and
/// `_boundary` steps on `Element = "p:2"` (full and serendipity p
/// elements), with the dofs from `elmer-pdofs`.
#[test]
fn model_pde_splitted_cases_match_pristine() {
    for (name, pristine) in [("ModelPDEsplitted", "9.23486039E-01"), ("SD_ModelPDEsplitted", "1.05785926E+00")] {
        let dir = case::stage(name, None);
        let out = case::run(name, &dir, "case.sif").unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(out.passed, Some(true), "{name}: reference check\n{}", out.log);
        case::assert_norm_prints(name, out.norm(1), pristine);
    }
}

/// `fem/tests/PoissonML` (`PoissonML.so`, `ElmerGrid 1 2 square`): the
/// matrix-free CG on `Element = "p:9"` quads, with the dofs from
/// `elmer-pdofs`.
#[test]
fn poisson_ml_matches_pristine() {
    let dir = case::stage("PoissonML", None);
    let out = case::run("PoissonML", &dir, "ptest.sif").unwrap_or_else(|e| panic!("PoissonML: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    case::assert_norm_prints("PoissonML", out.norm(1), "1.71580608E-01");
}

/// `fem/tests/SD_p-FEM_two_solvers` (`BenchmarkShell.so`, `ElmerGrid 1 2
/// eighth`): the shell solver on `p:4` serendipity elements next to
/// ModelPDE on `p:2`, both norms as pristine. The pristine run prints
/// `Relative energy error = 3.4142444709098718E-002`; the port gives
/// `3.4142444709132004E-002`.
#[test]
fn sd_p_fem_two_solvers_matches_pristine() {
    let name = "SD_p-FEM_two_solvers";
    let dir = case::stage(name, None);
    let out = case::run(name, &dir, "two_solvers.sif").unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    case::assert_norm_prints(name, out.norm(1), "1.72402774E-02");
    case::assert_norm_prints(name, out.norm(2), "2.10847920E-01");
}

/// `p-FEM_two_solvers` gives its p degrees body-wise (`Solver 1: Element`
/// in each Body), which the solver's dof layout does not see: the shell
/// solver refuses instead of running on nodal elements (pristine norm
/// 9.50681481E-03).
#[test]
fn p_fem_two_solvers_body_wise_elements_are_refused() {
    let dir = case::stage("p-FEM_two_solvers", None);
    let err = case::run("p-FEM_two_solvers", &dir, "two_solvers.sif").err().expect("refused");
    assert!(err.to_string().contains("body-wise element definitions"), "{err}");
}

/// `fem/tests/StokesPFEM` (`Stokes1.so`): Stokes on `p:2`. The case as
/// given is refused (`Set Dirichlet BCs By BC Numbering`); without that
/// keyword the pristine run prints the same norm, 5.41277436E+00, from a
/// BiCGStabL solve stopped at 1e-6. With a direct solver instead
/// (`ptest_direct.sif`) the pristine norm is 5.41277409E+00, which the port
/// must print; with the iterative solve it must stay within that solve's
/// tolerance.
#[test]
fn stokes_pfem_matches_pristine_without_bc_numbering() {
    let dir = case::stage("StokesPFEM", None);
    let out = case::run("StokesPFEM", &dir, "ptest_direct.sif").unwrap_or_else(|e| panic!("StokesPFEM: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    case::assert_norm_prints("StokesPFEM", out.norm(1), "5.41277409E+00");

    let dir = case::stage("StokesPFEM", None);
    let out = case::run("StokesPFEM", &dir, "ptest_without_bc_numbering.sif").unwrap_or_else(|e| panic!("StokesPFEM: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    let rel = (out.norm(1) - 5.41277436).abs() / 5.41277436;
    assert!(rel < 1e-6, "iterative norm {} vs pristine 5.41277436E+00", out.norm(1));
}

/// The case as given is refused by name before anything runs.
#[test]
fn stokes_pfem_as_given_is_refused_for_bc_numbering() {
    let dir = case::stage("StokesPFEM", None);
    let err = case::run("StokesPFEM", &dir, "ptest.sif").err().expect("refused");
    assert!(err.to_string().contains("BC Numbering"), "{err}");
}

/// `HeatControlExplicit` has its `AlternatingSource` registered, and is
/// refused by name before anything runs: its argument `cpar` comes from
/// `Apply Explicit Control`, which is not implemented.
#[test]
fn heat_control_explicit_is_refused_for_explicit_control() {
    let dir = case::stage("HeatControlExplicit", None);
    let err = case::run("HeatControlExplicit", &dir, "case.sif").err().expect("refused");
    assert!(err.to_string().contains("Apply Explicit Control"), "{err}");
}

/// `fem/tests/p-FEM_with_varying_p` (`BenchmarkShell.so`): the shell solver
/// with `Element = "p:%bodywisep"`, p = 4 in body 1 and 8 in body 2 from
/// the MATC function `bodywisep(cx)` (`GetMatcReal('bodywisep(cx)', 4,
/// [body, x, y, z], 'cx')`); norm as pristine (Reference Norm
/// 9.50681481E-03).
#[test]
fn p_fem_with_varying_p_matches_pristine() {
    let name = "p-FEM_with_varying_p";
    let dir = case::stage(name, None);
    let out = case::run(name, &dir, "cylindrical_shell.sif").unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(out.passed, Some(true), "reference check\n{}", out.log);
    case::assert_norm_prints(name, out.norm(1), "9.50681481E-03");
}
