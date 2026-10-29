//! `StressSolver` driven through `elmer_solver::Simulation`: a uniaxial
//! plane-stress bar has an exactly linear displacement field, which any mesh
//! of bilinear quads reproduces. Assembly and the direct solve are test-local
//! stand-ins for `elmer-assembly` and `elmer-direct`.

mod support;

use elmer_core::{
    CsrMatrix, DirectSolver, ElemId, ElementRef, ElementType, LinalgError, LocalSystem, MeshView, NullLogger,
    ParamMap, SerialComm, Value, NO_DOF,
};
use elmer_physics_elastic::{install_fe, register, FeServices};
use elmer_solver::backend::{AssemblyBackend, Perm};
use elmer_solver::{Model, Simulation, SolverError, SolverRegistry};
use support::{Gauss, P1Basis};

type Elem = (u16, i32, Vec<u32>, [Option<ElemId>; 2]);

struct Mesh {
    xyz: Vec<[f64; 3]>,
    /// (etype, tag, nodes, parents)
    elems: Vec<Elem>,
    n_bulk: usize,
}

impl MeshView for Mesh {
    fn dim(&self) -> usize {
        2
    }
    fn n_nodes(&self) -> usize {
        self.xyz.len()
    }
    fn coord(&self, n: u32) -> [f64; 3] {
        self.xyz[n as usize]
    }
    fn n_bulk(&self) -> usize {
        self.n_bulk
    }
    fn n_boundary(&self) -> usize {
        self.elems.len() - self.n_bulk
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let (t, tag, nodes, parents) = &self.elems[e.0 as usize];
        ElementRef { etype: ElementType(*t), nodes, tag: *tag, parents: *parents }
    }
}

fn bar_mesh() -> Mesh {
    let mut xyz = vec![];
    for j in 0..3 {
        for i in 0..3 {
            xyz.push([i as f64, j as f64 * 0.5, 0.0]);
        }
    }
    xyz[4] = [1.15, 0.55, 0.0]; // off-centre interior node
    let id = |i: u32, j: u32| 3 * j + i;
    let mut elems = vec![];
    for j in 0..2 {
        for i in 0..2 {
            elems.push((404, 1, vec![id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)], [None, None]));
        }
    }
    let n_bulk = elems.len();
    // left (tag 1), bottom (tag 2), right (tag 3)
    elems.push((202, 1, vec![id(0, 0), id(0, 1)], [Some(ElemId(0)), None]));
    elems.push((202, 1, vec![id(0, 1), id(0, 2)], [Some(ElemId(2)), None]));
    elems.push((202, 2, vec![id(0, 0), id(1, 0)], [Some(ElemId(0)), None]));
    elems.push((202, 2, vec![id(1, 0), id(2, 0)], [Some(ElemId(1)), None]));
    elems.push((202, 3, vec![id(2, 0), id(2, 1)], [Some(ElemId(1)), None]));
    elems.push((202, 3, vec![id(2, 1), id(2, 2)], [Some(ElemId(3)), None]));
    Mesh { xyz, elems, n_bulk }
}

struct DenseAssembly;

impl AssemblyBackend for DenseAssembly {
    fn build_perm(&self, mesh: &dyn MeshView, active: &[ElemId], _opt: bool) -> Perm {
        let mut perm = vec![NO_DOF; mesh.n_nodes()];
        let mut n = 0;
        for &e in active {
            for &nd in mesh.element(e).nodes {
                if perm[nd as usize] == NO_DOF {
                    perm[nd as usize] = n;
                    n += 1;
                }
            }
        }
        Perm { perm, n_slots: n as usize }
    }

    fn create_matrix(&self, _m: &dyn MeshView, _a: &[ElemId], perm: &Perm, dofs: usize) -> CsrMatrix<f64> {
        let n = perm.n_slots * dofs;
        CsrMatrix {
            nrows: n,
            ncols: n,
            row_ptr: (0..=n).map(|i| i * n).collect(),
            col_idx: (0..n).flat_map(|_| 0..n as u32).collect(),
            values: vec![0.0; n * n],
            diag: Some((0..n).map(|i| i * n + i).collect()),
        }
    }

    fn add_local(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], idx: &[u32], l: &LocalSystem) {
        let n = a.nrows;
        for (i, &gi) in idx.iter().enumerate() {
            if gi == NO_DOF {
                continue;
            }
            rhs[gi as usize] += l.force[i];
            for (j, &gj) in idx.iter().enumerate() {
                if gj != NO_DOF {
                    a.values[gi as usize * n + gj as usize] += l.stiff[i * l.n + j];
                }
            }
        }
    }

    fn apply_dirichlet(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], fixed: &[(u32, f64)], _p: &ParamMap) {
        let n = a.nrows;
        for &(d, v) in fixed {
            let d = d as usize;
            for j in 0..n {
                a.values[d * n + j] = 0.0;
            }
            a.values[d * n + d] = 1.0;
            rhs[d] = v;
        }
    }
}

struct Dense {
    a: Vec<f64>,
    n: usize,
}

impl DirectSolver<f64> for Dense {
    fn factor(&mut self, a: &CsrMatrix<f64>) -> Result<(), LinalgError> {
        self.n = a.nrows;
        self.a = a.values.clone();
        Ok(())
    }

    fn solve_in_place(&self, b: &mut [f64]) -> Result<(), LinalgError> {
        let n = self.n;
        let mut m = self.a.clone();
        for c in 0..n {
            let p = (c..n).max_by(|&x, &y| m[x * n + c].abs().partial_cmp(&m[y * n + c].abs()).unwrap()).unwrap();
            if m[p * n + c] == 0.0 {
                return Err(LinalgError::Singular(c));
            }
            for j in 0..n {
                m.swap(c * n + j, p * n + j);
            }
            b.swap(c, p);
            for r in c + 1..n {
                let f = m[r * n + c] / m[c * n + c];
                for j in c..n {
                    m[r * n + j] -= f * m[c * n + j];
                }
                b[r] -= f * b[c];
            }
        }
        for c in (0..n).rev() {
            let s: f64 = (c + 1..n).map(|j| m[c * n + j] * b[j]).sum();
            b[c] = (b[c] - s) / m[c * n + c];
        }
        Ok(())
    }
}

fn dense_factory(_: &ParamMap) -> Result<Box<dyn DirectSolver<f64> + Send>, SolverError> {
    Ok(Box::new(Dense { a: vec![], n: 0 }))
}

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::default();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn model(plane_stress: bool, extra_solver: &[(&str, Value)]) -> Model {
    let mut solver = vec![
        ("Equation", Value::Str("Stress Solver".into())),
        ("Procedure", Value::Str("\"StressSolve\" \"StressSolver\"".into())),
        ("Variable", Value::Str("Displacement".into())),
        ("Variable DOFs", Value::Integer(vec![2])),
        ("Linear System Solver", Value::Str("Direct".into())),
        ("Linear System Direct Method", Value::Str("umfpack".into())),
    ];
    solver.extend(extra_solver.iter().map(|(k, v)| (*k, v.clone())));
    Model {
        simulation: sec(&[
            ("Coordinate System", Value::Str("Cartesian 2D".into())),
            ("Simulation Type", Value::Str("Steady State".into())),
            ("Steady State Max Iterations", Value::Integer(vec![1])),
        ]),
        solvers: vec![sec(&solver)],
        equations: vec![sec(&[
            ("Active Solvers", Value::Integer(vec![1])),
            ("Plane Stress", Value::Logical(plane_stress)),
        ])],
        bodies: vec![sec(&[("Equation", Value::Integer(vec![1])), ("Material", Value::Integer(vec![1]))])],
        materials: vec![sec(&[
            ("Youngs Modulus", Value::real(100.0)),
            ("Poisson Ratio", Value::real(0.3)),
            ("Density", Value::real(1.0)),
        ])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", Value::Integer(vec![1])), ("Displacement 1", Value::real(0.0))]),
            sec(&[("Target Boundaries", Value::Integer(vec![2])), ("Displacement 2", Value::real(0.0))]),
            sec(&[("Target Boundaries", Value::Integer(vec![3])), ("Surface Traction 1", Value::real(5.0))]),
        ],
        ..Model::default()
    }
}

fn run(model: Model) -> Result<(Simulation, Mesh), SolverError> {
    let _ = install_fe(FeServices { basis: Box::new(P1Basis), quad: Box::new(Gauss) });
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mesh = bar_mesh();
    let mut sim = Simulation::new(model);
    sim.backends.assembly = Some(Box::new(DenseAssembly));
    sim.backends.linear.register_direct("umfpack", dense_factory);
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(&mesh, &NullLogger, &SerialComm)?;
    Ok((sim, mesh))
}

#[test]
fn uniaxial_plane_stress_bar_is_exact_and_aliases_are_views() {
    let (sim, mesh) = run(model(true, &[])).unwrap();
    let u = sim.vars.get("displacement").unwrap();
    let (e, nu, s) = (100.0, 0.3, 5.0);
    for n in 0..mesh.n_nodes() as u32 {
        let slot = u.slot(n).unwrap();
        let x = mesh.coord(n);
        let (ux, uy) = (u.values[2 * slot], u.values[2 * slot + 1]);
        assert!((ux - s * x[0] / e).abs() < 1e-12, "node {n}: ux {ux}");
        assert!((uy + nu * s * x[1] / e).abs() < 1e-12, "node {n}: uy {uy}");
    }
    // `Displacement 1` / `Displacement 2` are strided views of the same storage.
    let c2 = sim.vars.strided("displacement 2").unwrap();
    let slot = u.slot(8).unwrap();
    assert_eq!(c2.get(slot), u.values[2 * slot + 1]);
}

#[test]
fn plane_strain_differs_from_plane_stress() {
    let (a, _) = run(model(true, &[])).unwrap();
    let (b, _) = run(model(false, &[])).unwrap();
    let ux = |s: &Simulation| s.vars.get("displacement").unwrap().values[2 * 8];
    assert!((ux(&a) - ux(&b)).abs() > 1e-4);
}

#[test]
fn unimplemented_options_are_refused() {
    let r = run(model(true, &[("Large Deflection", Value::Logical(true))]));
    assert!(matches!(r, Err(SolverError::Unsupported(_))), "{:?}", r.err());
    let r = run(model(true, &[("Calculate Stresses", Value::Logical(true))]));
    assert!(matches!(r, Err(SolverError::Unsupported(_))));
}
