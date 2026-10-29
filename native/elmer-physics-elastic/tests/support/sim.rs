//! Test-local stand-ins for `elmer-assembly` and `elmer-direct` (dense), and a
//! mesh read from ElmerGrid files, for driving `StressSolver` through
//! `elmer_solver::Simulation`.
#![allow(dead_code)]

use elmer_core::{
    CsrMatrix, DirectSolver, ElemId, ElementRef, ElementType, LinalgError, LocalSystem, MeshView, ParamMap, Value,
    NO_DOF,
};
use elmer_solver::backend::{AssemblyBackend, Perm};

pub struct DenseAssembly;

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

    fn add_local(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], slots: &[u32], l: &LocalSystem) {
        let dofs = l.n / slots.len();
        let idx: Vec<u32> = slots
            .iter()
            .flat_map(|&s| (0..dofs as u32).map(move |c| if s == NO_DOF { NO_DOF } else { s * dofs as u32 + c }))
            .collect();
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
            // Elmer keeps the diagonal, so the row stays on the matrix's scale.
            let diag = a.values[d * n + d];
            for j in 0..n {
                a.values[d * n + j] = 0.0;
            }
            a.values[d * n + d] = diag;
            rhs[d] = diag * v;
        }
    }
}

pub struct Dense {
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
            let p = (c..n).max_by(|&x, &y| m[x * n + c].abs().total_cmp(&m[y * n + c].abs())).unwrap();
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

pub fn dense_factory(_: &ParamMap) -> Result<Box<dyn DirectSolver<f64> + Send>, elmer_solver::SolverError> {
    Ok(Box::new(Dense { a: vec![], n: 0 }))
}

pub fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::default();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

type Elem = (u16, i32, Vec<u32>, [Option<ElemId>; 2]);

/// A mesh in memory: `elems` lists the bulk elements first, then the boundary.
pub struct Mesh {
    pub dim: usize,
    pub xyz: Vec<[f64; 3]>,
    pub elems: Vec<Elem>,
    pub n_bulk: usize,
}

impl MeshView for Mesh {
    fn dim(&self) -> usize {
        self.dim
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

/// Reads `tests/data/<case>/mesh.{nodes,elements,boundary}`.
pub fn read_mesh(case: &str, dim: usize) -> Mesh {
    let read = |f: &str| {
        let p: std::path::PathBuf = [env!("CARGO_MANIFEST_DIR"), "tests", "data", case, f].iter().collect();
        std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
    };
    let xyz: Vec<[f64; 3]> = read("mesh.nodes")
        .lines()
        .map(|l| {
            let t: Vec<f64> = l.split_whitespace().skip(2).map(|x| x.parse().unwrap()).collect();
            [t[0], t[1], t[2]]
        })
        .collect();
    let mut elems: Vec<Elem> = read("mesh.elements")
        .lines()
        .map(|l| {
            let t: Vec<usize> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
            (t[2] as u16, t[1] as i32, t[3..].iter().map(|n| (n - 1) as u32).collect(), [None, None])
        })
        .collect();
    let n_bulk = elems.len();
    for l in read("mesh.boundary").lines() {
        let t: Vec<usize> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
        // id tag left right type nodes...
        let left = (t[2] != 0).then(|| ElemId((t[2] - 1) as u32));
        let right = (t[3] != 0).then(|| ElemId((t[3] - 1) as u32));
        elems.push((t[4] as u16, t[1] as i32, t[5..].iter().map(|n| (n - 1) as u32).collect(), [left, right]));
    }
    Mesh { dim, xyz, elems, n_bulk }
}
