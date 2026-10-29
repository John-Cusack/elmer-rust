//! The default backends end to end: perm, pattern, glue and Dirichlet from
//! `elmer-assembly`, a table-valued BC from `elmer-expr`, and a test-local
//! dense direct solver standing in for `elmer-direct`.
//!
//! These are local-correctness checks (exact discrete answers), not Elmer
//! parity captures.

use elmer_core::{
    Dependent, DependentKind, DirectSolver, ElemId, ElementRef, ElementType, LinalgError, LocalSystem, MeshView,
    NullLogger, ParamMap, SerialComm, Value,
};
use elmer_solver::backend::DirectFactory;
use elmer_solver::{Model, Simulation, Solver, SolverCtx, SolverError, SolverRegistry};

/// `n` line elements on [0, 1], then a point boundary element at each end
/// (tags 1 and 2).
struct Line {
    n: usize,
    conn: Vec<[u32; 2]>,
    ends: [u32; 2],
}

impl Line {
    fn new(n: usize) -> Self {
        Self {
            n,
            conn: (0..n as u32).map(|i| [i, i + 1]).collect(),
            ends: [0, n as u32],
        }
    }
}

impl MeshView for Line {
    fn dim(&self) -> usize {
        1
    }
    fn n_nodes(&self) -> usize {
        self.n + 1
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        [node as f64 / self.n as f64, 0.0, 0.0]
    }
    fn n_bulk(&self) -> usize {
        self.n
    }
    fn n_boundary(&self) -> usize {
        2
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let i = e.0 as usize;
        if i < self.n {
            ElementRef {
                etype: ElementType(202),
                nodes: &self.conn[i],
                tag: 1,
                parents: [None, None],
            }
        } else {
            let b = i - self.n;
            ElementRef {
                etype: ElementType(101),
                nodes: std::slice::from_ref(&self.ends[b]),
                tag: b as i32 + 1,
                parents: [Some(ElemId(if b == 0 { 0 } else { self.n as u32 - 1 })), None],
            }
        }
    }
}

/// Dense LU with partial pivoting.
#[derive(Default)]
struct DenseLu {
    n: usize,
    lu: Vec<f64>,
    piv: Vec<usize>,
}

impl DirectSolver<f64> for DenseLu {
    fn factor(&mut self, a: &elmer_core::CsrMatrix<f64>) -> Result<(), LinalgError> {
        let n = a.nrows;
        let mut m = vec![0.0; n * n];
        for i in 0..n {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                m[i * n + a.col_idx[k] as usize] += a.values[k];
            }
        }
        let mut piv: Vec<usize> = (0..n).collect();
        for c in 0..n {
            let p = (c..n).max_by(|&x, &y| m[x * n + c].abs().total_cmp(&m[y * n + c].abs())).unwrap();
            if m[p * n + c] == 0.0 {
                return Err(LinalgError::Singular(c));
            }
            if p != c {
                for j in 0..n {
                    m.swap(p * n + j, c * n + j);
                }
                piv.swap(p, c);
            }
            for r in c + 1..n {
                let f = m[r * n + c] / m[c * n + c];
                m[r * n + c] = f;
                for j in c + 1..n {
                    m[r * n + j] -= f * m[c * n + j];
                }
            }
        }
        *self = Self { n, lu: m, piv };
        Ok(())
    }

    fn solve_in_place(&self, b: &mut [f64]) -> Result<(), LinalgError> {
        let n = self.n;
        let mut x: Vec<f64> = self.piv.iter().map(|&p| b[p]).collect();
        for i in 0..n {
            for j in 0..i {
                x[i] -= self.lu[i * n + j] * x[j];
            }
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                x[i] -= self.lu[i * n + j] * x[j];
            }
            x[i] /= self.lu[i * n + i];
        }
        b.copy_from_slice(&x);
        Ok(())
    }
}

fn dense(_: &ParamMap) -> Result<Box<dyn DirectSolver<f64> + Send>, SolverError> {
    Ok(Box::new(DenseLu::default()))
}

/// -u'' = 1 with linear elements and the framework's Default* calls.
struct Poisson;

impl Solver for Poisson {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        for e in active {
            let nodes = ctx.mesh.element(e).nodes;
            let h = ctx.mesh.coord(nodes[1])[0] - ctx.mesh.coord(nodes[0])[0];
            let mut l = LocalSystem::zeros(2);
            l.stiff = vec![1.0 / h, -1.0 / h, -1.0 / h, 1.0 / h];
            l.force = vec![h / 2.0, h / 2.0];
            ctx.default_update(&l, e)?;
        }
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}

fn poisson(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Poisson))
}

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

/// -u'' = f on [0, 1] with u(0) = 0 and u(1) = a table on `Coordinate 1`
/// giving 2. `linsys` are the solver section's `Linear System *` entries.
/// Returns the simulation after `run` and the max nodal error against
/// x (1 - x) f / 2 + 2 x.
fn run_poisson(n: usize, linsys: &[(&str, Value)], extra: Option<(&str, DirectFactory)>) -> (Simulation, f64) {
    // u(1) = Variable Coordinate 1, Real 0 0 / 1 2: a table evaluated by elmer-expr.
    let table = Value::Dependent(Dependent {
        variables: vec!["Coordinate 1".into()],
        kind: DependentKind::Table {
            x: vec![0.0, 1.0],
            y: vec![0.0, 2.0],
            cubic: false,
        },
        rows: 1,
        cols: 1,
    });
    let mut solver = vec![
        ("Equation", Value::Str("poisson".into())),
        ("Variable", Value::Str("u".into())),
        ("Procedure", Value::Str("\"Test\" \"Poisson\"".into())),
        ("Steady State Convergence Tolerance", Value::real(1e-8)),
    ];
    solver.extend(linsys.iter().cloned());
    let model = Model {
        simulation: sec(&[
            ("Simulation Type", Value::Str("steady state".into())),
            ("Steady State Max Iterations", Value::Integer(vec![10])),
        ]),
        solvers: vec![sec(&solver)],
        equations: vec![sec(&[("Active Solvers", Value::Integer(vec![1]))])],
        bodies: vec![sec(&[("Equation", Value::Integer(vec![1]))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", Value::Integer(vec![1])), ("u", Value::real(0.0))]),
            sec(&[("Target Boundaries", Value::Integer(vec![2])), ("u", table)]),
        ],
        ..Model::default()
    };
    let mesh = Line::new(n);
    let mut reg = SolverRegistry::new();
    reg.register("Test", "Poisson", poisson);
    let mut sim = Simulation::new(model);
    if let Some((name, f)) = extra {
        sim.backends.linear.register_direct(name, f);
    }
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    let report = sim.run(&mesh, &NullLogger, &SerialComm).unwrap();
    // Linear elements are nodally exact in 1D; the second steady pass
    // reproduces the first, so the norm change is (near) zero.
    assert_eq!(report.steps, 2);
    let u = sim.vars.get("u").unwrap();
    let perm = u.perm.as_deref().unwrap();
    let mut err: f64 = 0.0;
    for node in 0..=n as u32 {
        let x = mesh.coord(node)[0];
        let exact = x * (1.0 - x) / 2.0 + 2.0 * x;
        err = err.max((u.values[perm[node as usize] as usize] - exact).abs());
    }
    assert_eq!(sim.slots[0].system.as_ref().unwrap().matrix.nrows, n + 1);
    (sim, err)
}

fn st(s: &str) -> Value {
    Value::Str(s.into())
}

#[test]
fn poisson_test_local_direct() {
    let ls = [("Linear System Solver", st("direct")), ("Linear System Direct Method", st("test dense"))];
    let (_, err) = run_poisson(12, &ls, Some(("test dense", dense)));
    assert!(err < 1e-12, "{err}");
}

/// No `Linear System Direct Method`: Elmer's default `banded` (elmer-direct band LU).
#[test]
fn poisson_default_banded_direct() {
    let (_, err) = run_poisson(12, &[("Linear System Solver", st("direct"))], None);
    assert!(err < 1e-12, "{err}");
}

#[test]
fn poisson_sparse_lu() {
    let ls = [("Linear System Solver", st("direct")), ("Linear System Direct Method", st("umfpack"))];
    let (_, err) = run_poisson(12, &ls, None);
    assert!(err < 1e-12, "{err}");
}

/// HUTI CG with ILU0 on the diagonally scaled system (Linear System Scaling defaults on).
#[test]
fn poisson_cg_ilu0() {
    let ls = [
        ("Linear System Solver", st("iterative")),
        ("Linear System Iterative Method", st("cg")),
        ("Linear System Preconditioning", st("ilu0")),
        ("Linear System Max Iterations", Value::Integer(vec![200])),
        ("Linear System Convergence Tolerance", Value::real(1e-13)),
        ("Linear System Residual Output", Value::Integer(vec![0])),
    ];
    let (_, err) = run_poisson(20, &ls, None);
    assert!(err < 1e-10, "{err}");
}

#[test]
fn poisson_bicgstab_diagonal_unscaled() {
    let ls = [
        ("Linear System Solver", st("iterative")),
        ("Linear System Iterative Method", st("bicgstab")),
        ("Linear System Preconditioning", st("diagonal")),
        ("Linear System Scaling", Value::Logical(false)),
        ("Linear System Max Iterations", Value::Integer(vec![200])),
        ("Linear System Convergence Tolerance", Value::real(1e-13)),
        ("Linear System Residual Output", Value::Integer(vec![0])),
    ];
    let (_, err) = run_poisson(20, &ls, None);
    assert!(err < 1e-10, "{err}");
}

#[test]
fn poisson_gcr_and_idrs_ilu1() {
    for m in ["gcr", "idrs", "bicgstabl"] {
        let ls = [
            ("Linear System Solver", st("iterative")),
            ("Linear System Iterative Method", st(m)),
            ("Linear System Preconditioning", st("ilu1")),
            ("Linear System Max Iterations", Value::Integer(vec![200])),
            ("Linear System Convergence Tolerance", Value::real(1e-13)),
            ("Linear System Residual Output", Value::Integer(vec![0])),
        ];
        let (_, err) = run_poisson(20, &ls, None);
        assert!(err < 1e-10, "{m}: {err}");
    }
}

#[test]
#[should_panic(expected = "Linear System Solver")]
fn missing_linear_system_solver_is_fatal() {
    run_poisson(4, &[], None);
}

/// A second solver sees the first one's update in the same coupled pass.
struct Writer;
impl Solver for Writer {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _t: bool) -> Result<(), SolverError> {
        let it = ctx.scalar("coupled iter").unwrap();
        ctx.primary_mut().unwrap().values.iter_mut().for_each(|v| *v = 10.0 + it);
        ctx.nonlinear_change()?;
        Ok(())
    }
}
struct Reader;
impl Solver for Reader {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _t: bool) -> Result<(), SolverError> {
        let a = ctx.variable("a").unwrap().values[0];
        ctx.primary_mut().unwrap().values.iter_mut().for_each(|v| *v = 2.0 * a);
        ctx.nonlinear_change()?;
        Ok(())
    }
}

#[test]
fn second_solver_sees_first_update() {
    let solver = |eq: &str, var: &str, entry: &str| {
        sec(&[
            ("Equation", Value::Str(eq.into())),
            ("Variable", Value::Str(var.into())),
            ("Procedure", Value::Str(format!("\"Test\" \"{entry}\""))),
            // Without a tolerance SteadyConverged stays -1, which SolveCoupled
            // counts as done after one pass.
            ("Steady State Convergence Tolerance", Value::real(1e-12)),
        ])
    };
    let model = Model {
        simulation: sec(&[
            ("Simulation Type", Value::Str("steady state".into())),
            ("Steady State Max Iterations", Value::Integer(vec![3])),
        ]),
        solvers: vec![solver("w", "a", "W"), solver("r", "b", "R")],
        equations: vec![sec(&[("Active Solvers", Value::Integer(vec![1, 2]))])],
        bodies: vec![sec(&[("Equation", Value::Integer(vec![1]))])],
        ..Model::default()
    };
    let mesh = Line::new(2);
    let mut reg = SolverRegistry::new();
    reg.register("Test", "W", |_| Ok(Box::new(Writer)));
    reg.register("Test", "R", |_| Ok(Box::new(Reader)));
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    let report = sim.run(&mesh, &NullLogger, &SerialComm).unwrap();
    assert_eq!(report.steps, 3);
    assert_eq!(sim.vars.get("a").unwrap().values, vec![13.0; 3]);
    assert_eq!(sim.vars.get("b").unwrap().values, vec![26.0; 3]);
}
