//! `Block Solver` without a `Procedure`, beyond the diagonal blocks:
//!
//! - Off-diagonal blocks: CoupledPoisson3 (`tests/data/coupled/cp3.sif`, two
//!   fields on the angle mesh) with a mass coupling `-0.5 M` in blocks (1,2)
//!   and (2,1), solved by blocks (GCR over UMFPACK blocks) and as the
//!   monolithic `Coupled Solver` system (UMFPACK): the same fields.
//! - Two-component block variables: the same with `Variable n Dofs = 2`
//!   and component-wise kernels.
//! - p-element block variables: two fields on line elements with a bubble
//!   dof each (the test-local `DofsBackend` of `solver_dofs.rs`), each block
//!   `-u'' = 1`, exact at the nodes and in the bubble coefficients; a bubble
//!   the boundary reports as an extra Dirichlet dof is fixed to zero in each
//!   block variable.
//!
//! Local-correctness checks, not Elmer parity captures (the pristine block
//! tests with these features need the BlockLinElast plugin).

mod support;

use std::path::Path;
use std::sync::Arc;

use elmer_assembly::{initialize_matrix, ListMatrix};
use elmer_core::{
    CsrMatrix, ElemId, ElementBasis, ElementRef, ElementType, LocalSystem, MeshView, NullLogger, ParamMap, Params,
    QuadPoint, Quadrature, SerialComm, Value, NO_DOF,
};
use elmer_solver::{
    DofsBackend, ElementAssembly, FeBackend, Model, Perm, Simulation, SolverCtx, SolverDofs, SolverError,
    SolverRegistry,
};
use support::fixtures::Bulk;
use support::sif_model;

/// `-0.5 M` on the element's nodes, component by component.
struct MassCoupling;

/// `grad u . grad v` and a unit source, component by component.
struct VecLaplace;

/// One Gauss point: weight, basis values and gradients.
type Point = (f64, Vec<f64>, Vec<[f64; 3]>);

/// The element's basis at its Gauss points.
fn at_points(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Vec<Point>, SolverError> {
    let el = ctx.mesh.element(e);
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let mut ev = elmer_core::BasisEval::default();
    let mut out = Vec::new();
    for qp in support::quad::TestQuadrature.gauss_points(el.etype, None) {
        support::basis::TestBasis
            .eval(el.etype, &xyz, qp.uvw, &mut ev)
            .map_err(|m| SolverError::fatal("test kernel", m))?;
        out.push((qp.weight * ev.det_j, ev.basis.clone(), ev.dbasis_dx.clone()));
    }
    Ok(out)
}

impl ElementAssembly for VecLaplace {
    fn assemble(
        &mut self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        _dt: f64,
        _transient: bool,
        local: &mut LocalSystem,
    ) -> Result<(), SolverError> {
        let n = ctx.mesh.element(e).nodes.len();
        let dofs = local.force.len() / n;
        let m = n * dofs;
        for (w, b, g) in at_points(ctx, e)? {
            for p in 0..n {
                for q in 0..n {
                    let k = w * (g[p][0] * g[q][0] + g[p][1] * g[q][1] + g[p][2] * g[q][2]);
                    for c in 0..dofs {
                        local.stiff[(dofs * p + c) * m + dofs * q + c] += k;
                    }
                }
                for c in 0..dofs {
                    local.force[dofs * p + c] += w * b[p] * (c + 1) as f64;
                }
            }
        }
        Ok(())
    }
}

impl ElementAssembly for MassCoupling {
    fn assemble(
        &mut self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        _dt: f64,
        _transient: bool,
        local: &mut LocalSystem,
    ) -> Result<(), SolverError> {
        let n = ctx.mesh.element(e).nodes.len();
        let dofs = local.force.len() / n;
        let m = n * dofs;
        for (w, b, _) in at_points(ctx, e)? {
            for p in 0..n {
                for q in 0..n {
                    for c in 0..dofs {
                        local.stiff[(dofs * p + c) * m + dofs * q + c] -= 0.5 * w * b[p] * b[q];
                    }
                }
            }
        }
        Ok(())
    }
}

fn cp3(block: bool) -> Simulation {
    cp3_with(block, 1)
}

/// CoupledPoisson3 with the coupling blocks; with `dofs` 2 the two fields
/// have two components each, assembled by the component-wise kernels.
fn cp3_with(block: bool, dofs: i32) -> Simulation {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join("angle"), 2);
    let mut model = sif_model(&dir.join("coupled/cp3.sif"));
    model.solvers.truncate(1);
    let p = &mut model.solvers[0];
    for k in ["Bulk Assembly Procedure 12", "Bulk Assembly Procedure 21"] {
        p.insert(k, Value::File("\"./Coupling\" \"Mass\"".into()));
    }
    p.0.remove("reference norm");
    if dofs > 1 {
        for i in 1..=2 {
            p.insert(&format!("Variable {i} Dofs"), Value::Integer(vec![dofs]));
            p.insert(&format!("Bulk Assembly Procedure {i}{i}"), Value::File("\"./Coupling\" \"Laplace\"".into()));
        }
        for bc in &mut model.boundary_conditions {
            for name in ["t1", "t2"] {
                if let Some(v) = bc.0.remove(name) {
                    for c in 1..=dofs {
                        bc.insert(&format!("{name} {c}"), v.clone());
                    }
                }
            }
        }
    }
    let p = &mut model.solvers[0];
    if block {
        p.0.remove("coupled solver");
        let s = |v: &str| Value::Str(v.into());
        for (k, v) in [
            ("Block Solver", Value::Logical(true)),
            ("Nonlinear System Max Iterations", Value::Integer(vec![1])),
            ("Block Preconditioner", Value::Logical(true)),
            ("Outer: Linear System Iterative Method", s("gcr")),
            ("Outer: Linear System Max Iterations", Value::Integer(vec![200])),
            ("Outer: Linear System Convergence Tolerance", Value::real(1e-13)),
            ("block 11: Linear System Solver", s("direct")),
            ("block 11: Linear System Direct Method", s("umfpack")),
            ("block 22: Linear System Solver", s("direct")),
            ("block 22: Linear System Direct Method", s("umfpack")),
        ] {
            p.insert(k, v);
        }
    } else {
        p.insert("Linear System Solver", Value::Str("direct".into()));
        p.insert("Linear System Direct Method", Value::Str("umfpack".into()));
    }
    let mut reg = SolverRegistry::new();
    reg.register_assembly("PoissonAssembly", "BulkAssembly", Arc::new(|_| Ok(Box::new(Bulk::default()))));
    reg.register_assembly("Coupling", "Mass", Arc::new(|_| Ok(Box::new(MassCoupling))));
    reg.register_assembly("Coupling", "Laplace", Arc::new(|_| Ok(Box::new(VecLaplace))));
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();
    sim
}

#[test]
fn off_diagonal_blocks_give_the_coupled_solution() {
    let (blocks, coupled) = (cp3(true), cp3(false));
    for name in ["t1", "t2"] {
        let (a, b) = (&blocks.vars.get(name).unwrap().values, &coupled.vars.get(name).unwrap().values);
        assert_eq!(a.len(), b.len());
        let scale = b.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(scale > 0.0);
        for (x, y) in a.iter().zip(b) {
            assert!((x - y).abs() <= 1e-9 * scale, "{name}: {x} vs {y}");
        }
    }
}

#[test]
fn two_component_block_variables() {
    let (blocks, coupled) = (cp3_with(true, 2), cp3_with(false, 2));
    for name in ["t1", "t2"] {
        let (a, b) = (blocks.vars.get(name).unwrap(), coupled.vars.get(name).unwrap());
        assert_eq!((a.dofs, a.values.len()), (2, b.values.len()));
        let scale = b.values.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(scale > 0.0);
        for (x, y) in a.values.iter().zip(&b.values) {
            assert!((x - y).abs() <= 1e-9 * scale, "{name}: {x} vs {y}");
        }
    }
}

/// `n` line elements on [0, 1], then a point boundary element at each end
/// (tags 1 and 2).
struct Line {
    n: usize,
    conn: Vec<[u32; 2]>,
    ends: [u32; 2],
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
            ElementRef { etype: ElementType(202), nodes: &self.conn[i], tag: 1, parents: [None, None] }
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

/// Nodes in slots `0..=n`, the bubble of element `j` in slot `n + 1 + j`,
/// two unused dofs at the end; the right end reports the last bubble as an
/// extra Dirichlet dof.
struct Bubbles {
    n: usize,
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
        if e.0 as usize == self.n + 1 {
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
        vec![e.0, e.0 + 1]
    }
    fn local_bubbles(&self, _e: ElemId) -> usize {
        0
    }
}

struct BubbleBackend;

impl DofsBackend for BubbleBackend {
    fn for_solver(
        &self,
        mesh: &dyn MeshView,
        _model: &Model,
        solver: &ParamMap,
        _active: &[ElemId],
        _dofs: usize,
    ) -> Option<Result<Box<dyn SolverDofs>, SolverError>> {
        solver.string("element")?;
        Some(Ok(Box::new(Bubbles { n: mesh.n_bulk() })))
    }
}

/// `-u'' = 1` on an element's two nodes and its bubble (`1 - xi^2`).
struct BubbleLaplace;

impl ElementAssembly for BubbleLaplace {
    fn assemble(
        &mut self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        _dt: f64,
        _transient: bool,
        local: &mut LocalSystem,
    ) -> Result<(), SolverError> {
        assert_eq!(local.force.len(), 3);
        let nodes = ctx.solver_dofs().expect("p layout").element_nodes(e);
        let h = ctx.mesh.coord(nodes[1])[0] - ctx.mesh.coord(nodes[0])[0];
        local.stiff = vec![1.0 / h, -1.0 / h, 0.0, -1.0 / h, 1.0 / h, 0.0, 0.0, 0.0, 16.0 / (3.0 * h)];
        local.force = vec![h / 2.0, h / 2.0, 2.0 * h / 3.0];
        Ok(())
    }
}

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

#[test]
fn p_element_block_variables() {
    let n = 4;
    let s = |v: &str| Value::Str(v.into());
    let model = Model {
        simulation: sec(&[
            ("Simulation Type", s("steady state")),
            ("Steady State Max Iterations", Value::Integer(vec![1])),
        ]),
        solvers: vec![sec(&[
            ("Equation", s("bubbles")),
            ("Block Solver", Value::Logical(true)),
            ("Element", s("p:2")),
            ("Variable 1", s("u")),
            ("Variable 2", s("v")),
            ("Bulk Assembly Procedure 11", Value::File("\"Test\" \"BubbleLaplace\"".into())),
            ("Bulk Assembly Procedure 22", Value::File("\"Test\" \"BubbleLaplace\"".into())),
            ("Block Preconditioner", Value::Logical(true)),
            ("Outer: Linear System Iterative Method", s("gcr")),
            ("Outer: Linear System Convergence Tolerance", Value::real(1e-14)),
            ("Outer: Linear System Max Iterations", Value::Integer(vec![50])),
            ("block 11: Linear System Solver", s("direct")),
            ("block 22: Linear System Solver", s("direct")),
        ])],
        equations: vec![sec(&[
            ("Active Solvers", Value::Integer(vec![1])),
            ("Active Variables", Value::Integer(vec![1, 2])),
        ])],
        bodies: vec![sec(&[("Equation", Value::Integer(vec![1]))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", Value::Integer(vec![1])), ("u", Value::real(0.0)), ("v", Value::real(0.0))]),
            // Only `v` is given on the right: only its last bubble is pinned.
            sec(&[("Target Boundaries", Value::Integer(vec![2])), ("v", Value::real(0.0))]),
        ],
        ..Model::default()
    };
    let mesh = Line { n, conn: (0..n as u32).map(|i| [i, i + 1]).collect(), ends: [0, n as u32] };
    let mut reg = SolverRegistry::new();
    reg.register_assembly("Test", "BubbleLaplace", Arc::new(|_| Ok(Box::new(BubbleLaplace))));
    let mut sim = Simulation::new(model);
    sim.backends.dofs = Some(Box::new(BubbleBackend));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();

    let h = 1.0 / n as f64;
    let v = sim.vars.get("v").unwrap();
    let perm = v.perm.as_deref().unwrap();
    // Numbered over the layout's dofs: 2n + 3 entries, 2n + 1 used.
    assert_eq!((perm.len(), v.values.len()), (2 * n + 3, 2 * n + 1));
    for node in 0..=n {
        let x = node as f64 * h;
        assert!((v.values[perm[node] as usize] - x * (1.0 - x) / 2.0).abs() < 1e-13, "v node {node}");
    }
    for j in 0..n - 1 {
        assert!((v.values[perm[n + 1 + j] as usize] - h * h / 8.0).abs() < 1e-13, "v bubble {j}");
    }
    // The extra Dirichlet dof of the right end, in v's own numbering.
    assert_eq!(v.values[perm[2 * n] as usize], 0.0);
    // u is free on the right: -u'' = 1, u(0) = 0, u'(1) = 0.
    let u = sim.vars.get("u").unwrap();
    let up = u.perm.as_deref().unwrap();
    for node in 0..=n {
        let x = node as f64 * h;
        assert!((u.values[up[node] as usize] - x * (2.0 - x) / 2.0).abs() < 1e-12, "u node {node}");
    }
    assert!((u.values[up[2 * n] as usize] - h * h / 8.0).abs() < 1e-13, "u last bubble");
}
