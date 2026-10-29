//! The `Backends::dofs` seam for non-nodal (p-element) solvers, with a
//! test-local `DofsBackend`: one bubble dof per line element after the
//! nodes. Linear elements plus bubbles solve `-u'' = 1` exactly, so the
//! nodal values and the bubble coefficients (`h^2 / 8`) are checked to
//! rounding, and the Perm, variable size, pattern, element dofs and the
//! extra Dirichlet rows come from the backend.
//!
//! Local-correctness checks, not Elmer parity captures.

use std::sync::Mutex;

use elmer_assembly::{initialize_matrix, ListMatrix};
use elmer_core::{
    CsrMatrix, Dependent, DependentKind, ElemId, ElementBasis, ElementRef, ElementType, LocalSystem, MeshView,
    NullLogger, ParamMap, QuadPoint, SerialComm, Value, NO_DOF,
};
use elmer_solver::{DofsBackend, Model, Perm, Simulation, Solver, SolverCtx, SolverDofs, SolverError, SolverRegistry};

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

/// Nodes in slots `0..=n`, the bubble of element `j` in slot `n + 1 + j`;
/// two unused slots at the end, as Elmer's `Ndeg` can exceed the used
/// count. With `pin_last_bubble`, the right end's boundary element reports
/// the last bubble as an extra Dirichlet slot.
struct Bubbles {
    n: usize,
    pin_last_bubble: bool,
}

impl SolverDofs for Bubbles {
    fn perm(&self) -> Perm {
        let used = 2 * self.n + 1;
        let mut perm: Vec<u32> = (0..used as u32).collect();
        perm.extend([NO_DOF, NO_DOF]);
        Perm { perm, n_slots: used }
    }
    fn create_matrix(&self, dofs: usize) -> CsrMatrix<f64> {
        let mut list = ListMatrix::new(2 * self.n + 1);
        let mut out = Vec::new();
        for j in 0..self.n as u32 {
            self.element_slots(ElemId(j), &mut out);
            list.add_clique(&mut out);
        }
        initialize_matrix(&list, None, dofs)
    }
    fn element_slots(&self, e: ElemId, out: &mut Vec<u32>) {
        out.clear();
        let j = e.0;
        if (j as usize) < self.n {
            out.extend([j, j + 1, self.n as u32 + 1 + j]);
        } else {
            out.push(if j as usize == self.n { 0 } else { self.n as u32 });
        }
    }
    fn boundary_extra_slots(&self, e: ElemId, out: &mut Vec<u32>) {
        out.clear();
        if self.pin_last_bubble && e.0 as usize == self.n + 1 {
            out.push(2 * self.n as u32);
        }
    }
    fn element_basis(&self, _e: ElemId) -> Result<Box<dyn ElementBasis + Send + Sync>, String> {
        Err("not used".into())
    }
    fn element_rule(&self, _e: ElemId) -> Result<Vec<QuadPoint>, String> {
        Err("not used".into())
    }
    fn element_nodes(&self, e: ElemId) -> Vec<u32> {
        let j = e.0;
        vec![j, j + 1]
    }
    fn local_bubbles(&self, _e: ElemId) -> usize {
        0
    }
}

struct BubbleBackend {
    pin_last_bubble: bool,
}

impl DofsBackend for BubbleBackend {
    fn for_solver(
        &self,
        mesh: &dyn MeshView,
        _model: &Model,
        solver: &ParamMap,
        active: &[ElemId],
        dofs: usize,
    ) -> Option<Result<Box<dyn SolverDofs>, SolverError>> {
        use elmer_core::Params;
        solver.string("element")?;
        assert_eq!((active.len(), dofs), (mesh.n_bulk(), 1));
        Some(Ok(Box::new(Bubbles {
            n: mesh.n_bulk(),
            pin_last_bubble: self.pin_last_bubble,
        })))
    }
}

/// The last element's dofs and the Dirichlet values, as the solver saw them.
type Seen = (Vec<u32>, Vec<(u32, f64)>);

/// What each run's solver saw, for the assertions.
static SEEN: Mutex<Vec<Seen>> = Mutex::new(Vec::new());

/// `-u'' = 1` with linear nodal functions and the bubble `1 - xi^2`
/// (orthogonal to them in energy), through the framework's `Default*` calls.
struct Poisson;

impl Solver for Poisson {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.default_initialize()?;
        assert!(ctx.solver_dofs().is_some());
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut idx = Vec::new();
        for &e in &active {
            ctx.element_dofs(e, &mut idx);
            assert_eq!(idx.len(), 3);
            let nodes = ctx.solver_dofs().unwrap().element_nodes(e);
            let h = ctx.mesh.coord(nodes[1])[0] - ctx.mesh.coord(nodes[0])[0];
            let mut l = LocalSystem::zeros(3);
            l.stiff = vec![
                1.0 / h,
                -1.0 / h,
                0.0,
                -1.0 / h,
                1.0 / h,
                0.0,
                0.0,
                0.0,
                16.0 / (3.0 * h),
            ];
            l.force = vec![h / 2.0, h / 2.0, 2.0 * h / 3.0];
            ctx.default_update(&l, e)?;
        }
        ctx.element_dofs(active[active.len() - 1], &mut idx);
        SEEN.lock().unwrap().push((idx.clone(), ctx.dirichlet_values()?));
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

fn run(n: usize, pin_last_bubble: bool, right: Value) -> Result<Simulation, SolverError> {
    let model = Model {
        simulation: sec(&[
            ("Simulation Type", Value::Str("steady state".into())),
            ("Steady State Max Iterations", Value::Integer(vec![1])),
            // Keeps a procedure BC value from being evaluated before the solve.
            ("Initialize Dirichlet Conditions", Value::Logical(false)),
        ]),
        solvers: vec![sec(&[
            ("Equation", Value::Str("poisson".into())),
            ("Variable", Value::Str("u".into())),
            ("Procedure", Value::Str("\"Test\" \"Poisson\"".into())),
            ("Element", Value::Str("p:2".into())),
            ("Linear System Solver", Value::Str("direct".into())),
        ])],
        equations: vec![sec(&[("Active Solvers", Value::Integer(vec![1]))])],
        bodies: vec![sec(&[("Equation", Value::Integer(vec![1]))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", Value::Integer(vec![1])), ("u", Value::real(0.0))]),
            sec(&[("Target Boundaries", Value::Integer(vec![2])), ("u", right)]),
        ],
        ..Model::default()
    };
    let mesh = Line::new(n);
    let mut reg = SolverRegistry::new();
    reg.register("Test", "Poisson", poisson);
    let mut sim = Simulation::new(model);
    sim.backends.dofs = Some(Box::new(BubbleBackend { pin_last_bubble }));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(&mesh, &NullLogger, &SerialComm)?;
    Ok(sim)
}

#[test]
fn bubble_dofs_solve_exactly() {
    let n = 4;
    let sim = run(n, false, Value::real(0.0)).unwrap();
    let u = sim.vars.get("u").unwrap();
    let perm = u.perm.as_deref().unwrap();
    // Perm and values come from the backend: Ndeg entries, 2n + 1 used.
    assert_eq!(perm.len(), 2 * n + 3);
    assert_eq!(u.values.len(), 2 * n + 1);
    assert_eq!(sim.slots[0].system.as_ref().unwrap().matrix.nrows, 2 * n + 1);
    let h = 1.0 / n as f64;
    for node in 0..=n {
        let x = node as f64 * h;
        let got = u.values[perm[node] as usize];
        assert!((got - x * (1.0 - x) / 2.0).abs() < 1e-14, "node {node}: {got}");
    }
    for j in 0..n {
        let got = u.values[perm[n + 1 + j] as usize];
        assert!((got - h * h / 8.0).abs() < 1e-14, "bubble {j}: {got}");
    }
    let seen = SEEN.lock().unwrap();
    let (idx, fixed) = seen
        .iter()
        .find(|(_, f)| f.len() == 2)
        .expect("run without the pinned bubble");
    assert_eq!(idx, &[n as u32 - 1, n as u32, 2 * n as u32]);
    assert_eq!(fixed, &[(0, 0.0), (n as u32, 0.0)]);
}

#[test]
fn extra_boundary_slots_are_fixed_to_zero() {
    let n = 3;
    let sim = run(n, true, Value::real(0.0)).unwrap();
    let u = sim.vars.get("u").unwrap();
    let perm = u.perm.as_deref().unwrap();
    // The last bubble is pinned, the others keep h^2 / 8.
    assert_eq!(u.values[perm[2 * n] as usize], 0.0);
    let h = 1.0 / n as f64;
    assert!((u.values[perm[n + 1] as usize] - h * h / 8.0).abs() < 1e-14);
    let seen = SEEN.lock().unwrap();
    assert!(seen
        .iter()
        .any(|(_, f)| f == &[(0, 0.0), (n as u32, 0.0), (2 * n as u32, 0.0)]));
}

#[test]
fn element_rule_np_defaults_to_a_refusal() {
    let b = Bubbles {
        n: 2,
        pin_last_bubble: false,
    };
    let err = b.element_rule_np(ElemId(0), 9).unwrap_err();
    assert!(err.contains("9 points"), "{err}");
}

#[test]
fn procedure_values_on_extra_slots_are_refused() {
    let proc_value = Value::Dependent(Dependent {
        variables: vec!["Coordinate 1".into()],
        kind: DependentKind::Procedure {
            library: "lib".into(),
            entry: "f".into(),
        },
        rows: 1,
        cols: 1,
    });
    let err = run(2, true, proc_value).err().expect("refused");
    assert!(err.to_string().contains("LocalBcBDOFs"), "{err}");
}
