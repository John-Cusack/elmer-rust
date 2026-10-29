//! `Apply Integral BCs`, checked against the pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) runs of `fem/tests/IntegralConstraintPoisson`
//! (an integral condition on a boundary) and `IntegralConstraintPoissonBF`
//! (on a body force) from their `case.sif` files on their ElmerGrid meshes
//! (`tests/data/intbc/`), run at `Max Output Level = 10`: the projector,
//! constraint matrix and restricted-solve messages, the `ComputeChange` lines
//! (the nonlinear norm includes the Lagrange multiplier) and the
//! reference-norm check. HeatSolve is the test-local stand-in.

mod support;

use std::path::Path;

use elmer_core::SerialComm;
use elmer_solver::{FeBackend, Simulation, SolverRegistry};
use support::poisson::heat;
use support::{sif_model, split_change, Capture};

/// The messages up to `Max Output Level = 10`.
struct Level10(Capture);

impl elmer_core::Logger for Level10 {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        if level <= 10 {
            self.0.info(caller, msg, level);
        }
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.0.warn(caller, msg);
    }
}

const CALLERS: &[&str] = &[
    "IntegralProjector",
    "GenerateConstraintMatrix",
    "SolveWithLinearRestriction",
    "ComputeChange",
    "CompareTo",
];

fn check(case: &str) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/intbc").join(case);
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let model = sif_model(&dir.join("case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat);
    let log = Level10(Capture::default());
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap_or_else(|e| panic!("{case}: {e}"));

    let stdout = std::fs::read_to_string(dir.join("stdout.txt")).unwrap();
    let keep = |l: &&str| !l.contains("Relative Error") && !l.contains("all 1 tests");
    let got = log.0.lines(CALLERS);
    let got: Vec<&str> = got.iter().map(String::as_str).filter(keep).collect();
    let want: Vec<&str> = stdout.lines().filter(keep).collect();
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        if w.contains("(NRM,RELC)") {
            let (gh, gn, gr, gt) = split_change(g);
            let (wh, wn, wr, wt) = split_change(w);
            assert_eq!((gh, gn, gt), (wh, wn, wt), "{case} line {i}");
            assert!((gr - wr).abs() <= 1e-6 * wr, "{case} line {i}\n{g}\n{w}");
        } else {
            assert_eq!(g, w, "{case} line {i}");
        }
    }
    assert_eq!(got.len(), want.len(), "{case}: {got:#?}");
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1), "{case}");
}

#[test]
fn integral_condition_on_a_boundary() {
    check("poisson");
}

#[test]
fn integral_condition_on_a_body_force() {
    check("poissonbf");
}

/// The `IntegralConstraintPoisson` case (without `Apply Integral BCs` unless
/// `integral`) with the one row `x[col] = value`, standing for dof `rep`,
/// handed over as the solver's own constraint matrix.
fn own_constraints(row: impl Fn(&[u32]) -> (u32, u32, f64), integral: bool) -> Result<Simulation, elmer_solver::SolverError> {
    use elmer_core::CsrMatrix;
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/intbc/poisson");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let mut model = sif_model(&dir.join("case.sif"));
    if !integral {
        model.solvers[0].0.remove("apply integral bcs");
    }
    model.solvers[0].0.remove("reference norm");
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    let v = sim.vars.get("temperature").unwrap();
    let n = v.values.len();
    let (col, rep, value) = row(v.perm.as_deref().unwrap());
    let c = elmer_solver::ConstraintMatrix {
        matrix: CsrMatrix {
            nrows: 1,
            ncols: n,
            row_ptr: vec![0, 1],
            col_idx: vec![col],
            values: vec![1.0],
            diag: None,
        },
        rhs: vec![value],
        inv_perm: vec![rep],
    };
    elmer_solver::SolverCtx::new(&mut sim, 0, &mesh, &log, &SerialComm).set_constraint_matrix(Some(c));
    sim.run(&mesh, &log, &SerialComm)?;
    Ok(sim)
}

#[test]
fn a_solvers_own_constraint_is_enforced() {
    // Node 127, the integral condition's representative, has no Dirichlet condition.
    let sim = own_constraints(|perm| (perm[126], perm[126], 0.25), false).unwrap();
    let v = sim.vars.get("temperature").unwrap();
    let x = v.values[v.perm.as_deref().unwrap()[126] as usize];
    assert!((x - 0.25).abs() < 1e-12, "{x}");
    assert!(v.values.iter().any(|&t| t != 0.0 && t != 0.25));
}

#[test]
fn a_constraint_on_a_dirichlet_dof_is_dropped() {
    // Node 1 is on boundary 2, Temperature = 0: the row becomes lambda = 0,
    // which leaves the unforced problem with its zero solution.
    let sim = own_constraints(|perm| (perm[0], perm[0], 0.25), false).unwrap();
    let v = sim.vars.get("temperature").unwrap();
    assert!(v.values.iter().all(|&t| t.abs() < 1e-14), "{:?}", v.values);
}

#[test]
fn own_constraints_with_integral_bcs_are_refused() {
    let e = own_constraints(|perm| (perm[126], perm[126], 0.25), true).err().unwrap();
    assert!(matches!(e, elmer_solver::SolverError::Unsupported(_)), "{e}");
}
