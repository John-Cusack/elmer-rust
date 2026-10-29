//! Checks against the pristine ElmerSolver run of `fem/tests/DirichletNeumann`
//! (fixtures and provenance in `tests/data/`).

use std::sync::{Arc, Mutex};

use elmer_core::{ElemId, ElementRef, ElementType, FieldKind, Logger, MeshView, ParamMap, SerialComm, Value};
use elmer_solver::norm::{compute_norm, NormSpec};
use elmer_solver::{
    compare_to_reference, finalize_reference, Model, Simulation, Solver, SolverCtx, SolverError, SolverRegistry,
    SolverSlot, Variable,
};

const RESULT: &str = include_str!("data/dirichlet_neumann.result");
const STDOUT: &str = include_str!("data/dirichlet_neumann_stdout.txt");

/// Collects lines as Elmer prints them with the default output prefixes.
#[derive(Default)]
struct Capture(Mutex<Vec<String>>);

impl Logger for Capture {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        if level <= 5 {
            self.0.lock().unwrap().push(format!("{caller}: {msg}"));
        }
    }
    fn warn(&self, caller: &str, msg: &str) {
        self.0.lock().unwrap().push(format!("WARNING:: {caller}: {msg}"));
    }
}

impl Capture {
    fn lines(&self, prefix: &str) -> Vec<String> {
        self.0.lock().unwrap().iter().filter(|l| l.starts_with(prefix)).cloned().collect()
    }
}

/// A test-local `.result` reader: (name, perm as 0-based slots, values in
/// storage order) for each variable.
fn read_result(text: &str) -> Vec<(String, Vec<u32>, Vec<f64>)> {
    let lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    let mut out = Vec::new();
    let mut prev: Option<(usize, Vec<(usize, usize)>)> = None;
    let mut i = 0;
    while i < lines.len() {
        let name = lines[i].trim().to_string();
        i += 1;
        let p: Vec<&str> = lines[i].split_whitespace().collect();
        i += 1;
        let (size, pairs) = if p[1] == "use" {
            prev.clone().unwrap()
        } else {
            let size: usize = p[1].parse().unwrap();
            let nz: usize = p[2].parse().unwrap();
            let pairs: Vec<(usize, usize)> = lines[i..i + nz]
                .iter()
                .map(|l| {
                    let v: Vec<usize> = l.split_whitespace().map(|t| t.parse().unwrap()).collect();
                    (v[0], v[1])
                })
                .collect();
            i += nz;
            prev = Some((size, pairs.clone()));
            (size, pairs)
        };
        let mut perm = vec![u32::MAX; size];
        let mut vals = vec![0.0; pairs.len()];
        for &(node, k) in &pairs {
            perm[node - 1] = (k - 1) as u32;
            vals[k - 1] = lines[i].trim().parse().unwrap();
            i += 1;
        }
        out.push((name, perm, vals));
    }
    out
}

fn field<'a>(f: &'a [(String, Vec<u32>, Vec<f64>)], name: &str) -> &'a (String, Vec<u32>, Vec<f64>) {
    f.iter().find(|(n, _, _)| n == name).unwrap()
}

/// A mesh with `n` nodes and one point element per node, tag = body.
struct PointMesh {
    nodes: Vec<u32>,
    tags: Vec<i32>,
}

impl MeshView for PointMesh {
    fn dim(&self) -> usize {
        2
    }
    fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
    fn coord(&self, _node: u32) -> [f64; 3] {
        [0.0; 3]
    }
    fn n_bulk(&self) -> usize {
        self.nodes.len()
    }
    fn n_boundary(&self) -> usize {
        0
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let i = e.0 as usize;
        ElementRef {
            etype: ElementType(101),
            nodes: std::slice::from_ref(&self.nodes[i]),
            tag: self.tags[i],
            parents: [None, None],
        }
    }
}

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

/// `Var % Norm` of the final fields reproduces the bits behind Fortran's
/// `CompareToReferenceSolution` lines: the absolute error to a reference
/// norm 1e-11 away is printed to 7 digits, which pins the norm far below
/// one ulp. The whole message lines must match byte for byte.
#[test]
fn reference_norm_lines_are_byte_identical() {
    let fields = read_result(RESULT);
    let mut sim = Simulation::new(Model::default());
    let refs = [("tempb", 1, 1.54633630735), ("tempa", 2, 0.6228119325)];
    for (name, id, refnorm) in refs {
        let (_, perm, vals) = field(&fields, name);
        let spec = NormSpec::l2(1);
        let mut v = Variable::new(name, FieldKind::Nodal, 1, Some(Arc::from(perm.as_slice())), vals.clone());
        v.conv.norm = compute_norm(&spec, vals, vals.len(), &SerialComm);
        let mut slot = SolverSlot::new(
            id,
            sec(&[("Reference Norm", Value::real(refnorm)), ("Reference Norm Absolute", Value::Logical(true))]),
        )
        .unwrap();
        slot.variable = Some(sim.vars.add(v).unwrap());
        sim.slots.push(slot);
    }
    let log = Capture::default();
    let mesh = PointMesh {
        nodes: (0..240).collect(),
        tags: vec![1; 240],
    };
    compare_to_reference(&mut sim, &mesh, &log, &SerialComm).unwrap();
    let dir = std::env::temp_dir().join(format!("elmer-solver-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    assert_eq!(finalize_reference(&sim.reference, &dir, &log, &SerialComm).unwrap(), Some(true));
    assert_eq!(std::fs::read_to_string(dir.join("TEST.PASSED")).unwrap(), "1\n");
    std::fs::remove_dir_all(&dir).unwrap();

    let want: Vec<String> = RESULT
        .lines()
        .filter_map(|l| l.strip_prefix("#> "))
        .map(str::to_string)
        .collect();
    assert_eq!(log.lines("CompareToReferenceSolution"), want);
}

/// The unmodified test's relative-error line (reference norm from the sif,
/// printed relative error 1.069559E-15).
#[test]
#[allow(clippy::excessive_precision)] // the literal exactly as case.sif gives it
fn relative_error_line_matches() {
    let fields = read_result(RESULT);
    let (_, perm, vals) = field(&fields, "tempa");
    let mut sim = Simulation::new(Model::default());
    let mut v = Variable::new("tempa", FieldKind::Nodal, 1, Some(Arc::from(perm.as_slice())), vals.clone());
    v.conv.norm = compute_norm(&NormSpec::l2(1), vals, vals.len(), &SerialComm);
    let mut slot = SolverSlot::new(2, sec(&[("Reference Norm", Value::real(0.62281193251726608))])).unwrap();
    slot.variable = Some(sim.vars.add(v).unwrap());
    sim.slots.push(slot);
    let log = Capture::default();
    let mesh = PointMesh {
        nodes: (0..240).collect(),
        tags: vec![1; 240],
    };
    compare_to_reference(&mut sim, &mesh, &log, &SerialComm).unwrap();
    let want: Vec<&str> = STDOUT.lines().filter(|l| l.starts_with("CompareToReferenceSolution:")).take(2).collect();
    assert_eq!(log.lines("CompareToReferenceSolution"), want);
}

/// Replays the solver norms printed by the Fortran run.
struct Scripted {
    norms: Vec<f64>,
    call: usize,
}

impl Solver for Scripted {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let x = self.norms[self.call];
        self.call += 1;
        ctx.primary_mut().unwrap().values[0] = x;
        ctx.nonlinear_change()?;
        Ok(())
    }
}

fn scripted_norms(name: &str) -> Vec<f64> {
    STDOUT
        .lines()
        .filter(|l| l.starts_with("ComputeChange: NS") && l.ends_with(&format!(":: {name}")))
        .map(|l| {
            let inner = l.split("(NRM,RELC): (").nth(1).unwrap();
            inner.split_whitespace().next().unwrap().parse().unwrap()
        })
        .collect()
}

thread_local! {
    static NORMS: std::cell::RefCell<Vec<Vec<f64>>> = const { std::cell::RefCell::new(Vec::new()) };
}

fn make_b(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Scripted {
        norms: NORMS.with(|n| n.borrow()[0].clone()),
        call: 0,
    }))
}

fn make_a(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Scripted {
        norms: NORMS.with(|n| n.borrow()[1].clone()),
        call: 0,
    }))
}

/// Two solvers on two bodies in a steady coupled loop, the DirichletNeumann
/// shape: the scheduler must run 14 coupled iterations, solver 1 before
/// solver 2, and print the same NS/SS lines (NRM exactly; RELC is computed
/// from 8-digit printed norms, so only its first digits can agree).
#[test]
fn coupled_loop_iteration_count() {
    NORMS.with(|n| *n.borrow_mut() = vec![scripted_norms("heatb"), scripted_norms("heata")]);
    let solver = |eq: &str, var: &str, proc_: &str| {
        sec(&[
            ("Equation", Value::Str(eq.into())),
            ("Variable", Value::Str(var.into())),
            ("Procedure", Value::Str(format!("\"Stub\" \"{proc_}\""))),
            ("Steady State Convergence Tolerance", Value::real(1.0e-5)),
            ("Nonlinear System Convergence Tolerance", Value::real(1.0e-5)),
            ("Nonlinear System Max Iterations", Value::Integer(vec![1])),
        ])
    };
    let model = Model {
        simulation: sec(&[
            ("Simulation Type", Value::Str("Steady State".into())),
            ("Steady State Max Iterations", Value::Integer(vec![100])),
            ("Output Intervals", Value::Integer(vec![0])),
        ]),
        // The SIF reader lowercases string values (stdout prints `heatb`).
        solvers: vec![solver("heatb", "tempb", "B"), solver("heata", "tempa", "A")],
        equations: vec![
            sec(&[("Active Solvers", Value::Integer(vec![2]))]),
            sec(&[("Active Solvers", Value::Integer(vec![1]))]),
        ],
        bodies: vec![
            sec(&[("Equation", Value::Integer(vec![1]))]),
            sec(&[("Equation", Value::Integer(vec![2]))]),
        ],
        ..Model::default()
    };
    let mesh = PointMesh {
        nodes: vec![0, 1],
        tags: vec![1, 2],
    };
    let mut reg = SolverRegistry::new();
    reg.register("Stub", "B", make_b);
    reg.register("Stub", "A", make_a);
    let mut sim = Simulation::new(model);
    let log = Capture::default();
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    let report = sim.run(&mesh, &log, &SerialComm).unwrap();
    assert_eq!(report.steps, 14);
    assert!(report.steady_state_reached);

    let got = log.lines("ComputeChange");
    let want: Vec<&str> = STDOUT.lines().filter(|l| l.starts_with("ComputeChange")).collect();
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(&want) {
        let split = |l: &str| {
            let (head, rest) = l.split_once("(NRM,RELC): (").unwrap();
            let nrm = rest[..15].to_string();
            let relc: f64 = rest[15..30].trim().parse().unwrap();
            let tail = rest[30..].to_string();
            (head.to_string(), nrm, relc, tail)
        };
        let (gh, gn, gr, gt) = split(g);
        let (wh, wn, wr, wt) = split(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        assert!((gr - wr).abs() <= 0.05 * wr, "{g}\n{w}");
    }
}
