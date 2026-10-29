//! Done criteria for elmer-core (swarm plan, lane L00).

use elmer_core::*;
use num_complex::Complex64;

#[test]
fn normalize_key_lowercases_and_collapses_whitespace() {
    assert_eq!(
        normalize_key("Linear System Solver"),
        "linear system solver"
    );
    assert_eq!(
        normalize_key("  Linear   System\tSolver  "),
        "linear system solver"
    );
    assert_eq!(normalize_key("MAX OUTPUT LEVEL"), "max output level");
    assert_eq!(
        normalize_key("Steady State Max Iterations"),
        normalize_key("steady  state max\t\titerations")
    );
    assert_eq!(normalize_key(""), "");
}

#[test]
fn param_map_lookups_are_key_normalised() {
    let mut p = ParamMap::new();
    p.insert("Linear System  Max Iterations", Value::Integer(vec![500]));
    p.insert("Linear System Convergence Tolerance", Value::real(1e-8));
    p.insert("Stabilize", Value::Logical(true));
    p.insert("Equation", Value::Str("Heat Equation".into()));
    p.insert("Mesh DB", Value::File("mesh".into()));
    p.insert(
        "Heat Conductivity",
        Value::Dependent(Dependent {
            variables: vec!["temperature".into()],
            kind: DependentKind::Table {
                x: vec![0.0, 1.0],
                y: vec![1.0, 2.0],
                cubic: false,
            },
            rows: 1,
            cols: 1,
        }),
    );

    assert_eq!(p.integer("linear system max iterations"), Some(500));
    assert_eq!(p.integer("LINEAR SYSTEM\tMAX ITERATIONS"), Some(500));
    assert_eq!(p.real("linear system convergence tolerance"), Some(1e-8));
    assert_eq!(p.logical("stabilize"), Some(true));
    assert_eq!(p.string("equation"), Some("Heat Equation"));
    assert_eq!(p.string("mesh db"), Some("mesh"));
    // A dependent value is not a constant.
    assert_eq!(p.real("heat conductivity"), None);
    assert!(p.get("heat conductivity").is_some());
    // Wrong type or absent.
    assert_eq!(p.real("stabilize"), None);
    assert_eq!(p.integer("nonexistent"), None);
    // Re-insert under a differently spelled key replaces.
    assert!(p.insert("STABILIZE", Value::Logical(false)).is_some());
    assert_eq!(p.logical("Stabilize"), Some(false));
    assert_eq!(p.0.len(), 6);
}

#[test]
fn serial_comm_sums_and_maxes_as_identity() {
    let c = SerialComm;
    assert_eq!((c.rank(), c.size()), (0, 1));
    let mut buf = [1.5, -2.0, 3.0];
    c.sum_f64(&mut buf);
    assert_eq!(buf, [1.5, -2.0, 3.0]);
    let s: f64 = buf.iter().sum();
    let mut total = [s];
    c.sum_f64(&mut total);
    assert_eq!(total[0], 2.5);
    c.max_f64(&mut buf);
    assert_eq!(buf, [1.5, -2.0, 3.0]);
    c.barrier();
    let _: &dyn Communicator = &c;
}

/// Every `Code` / `Nodes` pair in fem/src/elements.def (pin 9f6af2f85).
const ELEMENTS_DEF: &[(u16, usize)] = &[
    (101, 1),
    (102, 2),
    (202, 2),
    (203, 3),
    (204, 4),
    (303, 3),
    (306, 6),
    (310, 10),
    (404, 4),
    (408, 8),
    (409, 9),
    (412, 12),
    (416, 16),
    (504, 4),
    (510, 10),
    (605, 5),
    (613, 13),
    (614, 14),
    (706, 6),
    (715, 15),
    (718, 18),
    (808, 8),
    (820, 20),
    (827, 27),
];

#[test]
fn element_type_n_nodes_matches_elements_def() {
    for &(code, nodes) in ELEMENTS_DEF {
        let t = ElementType(code);
        assert_eq!(t.n_nodes(), nodes, "code {code}");
        assert_eq!(t.family(), code / 100, "code {code}");
    }
    assert_eq!(ElementType(827).family(), 8);
}

/// [[4,1,0],[1,4,1],[0,1,4]] with an extra column 3 entry in row 0.
fn sample() -> CsrMatrix {
    CsrMatrix {
        nrows: 3,
        ncols: 4,
        row_ptr: vec![0, 3, 6, 8],
        col_idx: vec![0, 1, 3, 0, 1, 2, 1, 2],
        values: vec![4.0, 1.0, 2.0, 1.0, 4.0, 1.0, 1.0, 4.0],
        diag: Some(vec![0, 4, 7]),
    }
}

#[test]
fn csr_apply_and_transpose() {
    let a = sample();
    assert_eq!(a.nnz(), 8);
    assert_eq!((a.nrows(), a.ncols()), (3, 4));
    let mut y = [0.0; 3];
    a.apply(&[1.0, 2.0, 3.0, 4.0], &mut y);
    assert_eq!(y, [4.0 + 2.0 + 8.0, 1.0 + 8.0 + 3.0, 2.0 + 12.0]);
    let mut yt = [9.0; 4];
    a.apply_transpose(&[1.0, 2.0, 3.0], &mut yt).unwrap();
    assert_eq!(yt, [4.0 + 2.0, 1.0 + 8.0 + 3.0, 2.0 + 12.0, 2.0]);
    let mut short = [0.0; 2];
    assert!(matches!(
        a.apply_transpose(&[1.0, 2.0, 3.0], &mut short),
        Err(LinalgError::Dim(_))
    ));
    let _: &dyn LinearOperator = &a;
}

#[test]
fn complex_csr_and_scalar() {
    let i = Complex64::new(0.0, 1.0);
    let a = CsrMatrix::<Complex64> {
        nrows: 1,
        ncols: 1,
        row_ptr: vec![0, 1],
        col_idx: vec![0],
        values: vec![i],
        diag: None,
    };
    let mut y = [Complex64::zero()];
    a.apply(&[Complex64::one()], &mut y);
    assert_eq!(y[0], i);
    assert_eq!(Scalar::conj(i), Complex64::new(0.0, -1.0));
    assert_eq!(Scalar::abs(Complex64::new(3.0, 4.0)), 5.0);
    assert_eq!(<f64 as Scalar>::abs(-2.0), 2.0);
    assert_eq!(
        <Complex64 as Scalar>::from_f64(2.0),
        Complex64::new(2.0, 0.0)
    );
}

#[test]
fn identity_precond_copies() {
    let mut z = [0.0; 3];
    Preconditioner::<f64>::apply(&IdentityPrecond, &[1.0, 2.0, 3.0], &mut z);
    assert_eq!(z, [1.0, 2.0, 3.0]);
    let mut zt = [0.0; 3];
    Preconditioner::<f64>::apply_transpose(&IdentityPrecond, &[3.0, 2.0, 1.0], &mut zt).unwrap();
    assert_eq!(zt, [3.0, 2.0, 1.0]);
}

/// The shared traits must stay object safe: lanes pass them as `dyn`.
struct Jacobi;
impl IterativeSolver for Jacobi {
    fn solve(
        &mut self,
        a: &dyn LinearOperator,
        m: &dyn Preconditioner,
        b: &[f64],
        x: &mut [f64],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError> {
        let n = a.nrows();
        let (mut ax, mut r, mut z) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        for it in 1..=200 {
            a.apply(x, &mut ax);
            for k in 0..n {
                r[k] = b[k] - ax[k];
            }
            let mut nrm = [r.iter().map(|v| v * v).sum::<f64>()];
            comm.sum_f64(&mut nrm);
            if nrm[0].sqrt() < 1e-12 {
                return Ok(SolveStats {
                    iterations: it,
                    residual: nrm[0].sqrt(),
                    status: SolveStatus::Converged,
                });
            }
            m.apply(&r, &mut z);
            for k in 0..n {
                x[k] += z[k];
            }
        }
        Ok(SolveStats {
            iterations: 200,
            residual: f64::NAN,
            status: SolveStatus::MaxIterations,
        })
    }
}

struct DiagPrecond(Vec<f64>);
impl Preconditioner for DiagPrecond {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for (k, zk) in z.iter_mut().enumerate() {
            *zk = r[k] / self.0[k];
        }
    }
}

#[test]
fn traits_compose_through_dyn() {
    let mut a = sample();
    a.ncols = 3;
    a.row_ptr = vec![0, 2, 5, 7];
    a.col_idx = vec![0, 1, 0, 1, 2, 1, 2];
    a.values = vec![4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0];
    let m = DiagPrecond(vec![4.0; 3]);
    let mut x = [0.0; 3];
    let b = [5.0, 6.0, 5.0];
    let st = Jacobi.solve(&a, &m, &b, &mut x, &SerialComm).unwrap();
    assert_eq!(st.status, SolveStatus::Converged);
    for v in x {
        assert!((v - 1.0).abs() < 1e-10);
    }
}

#[test]
fn local_system_zeros_and_logger() {
    let l = LocalSystem::zeros(4);
    assert_eq!((l.n, l.stiff.len(), l.force.len()), (4, 16, 4));
    assert!(l.mass.is_none() && l.damp.is_none());
    let log: &dyn Logger = &NullLogger;
    log.info("test", "msg", 3);
    log.warn("test", "msg");
    assert_eq!(NO_DOF, u32::MAX);
}
