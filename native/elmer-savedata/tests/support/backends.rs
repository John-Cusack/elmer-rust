//! Test-local stand-ins for `elmer-assembly` (matrix structure, scatter,
//! Dirichlet), `elmer-timeint` (backward Euler) and a dense direct solver,
//! plus a `MeshView` over the test mesh. Replace with the real crates when
//! they merge.

use elmer_core::{
    CsrMatrix, DirectSolver, ElemId, ElementRef, ElementType, LinalgError, LocalSystem, MeshView,
    ParamMap, NO_DOF,
};
use elmer_solver::backend::{AssemblyBackend, Perm, TimeBackend, TimeTerms};
use elmer_solver::SolverError;

use super::meshio::Mesh;

pub struct SimMesh {
    pub mesh: Mesh,
    nodes: Vec<Vec<u32>>,
}

impl SimMesh {
    pub fn new(mesh: Mesh) -> Self {
        let nodes = mesh
            .bulk
            .iter()
            .chain(&mesh.boundary)
            .map(|e| e.nodes.iter().map(|&n| n as u32).collect())
            .collect();
        Self { mesh, nodes }
    }
}

impl MeshView for SimMesh {
    fn dim(&self) -> usize {
        self.mesh.dim
    }
    fn n_nodes(&self) -> usize {
        self.mesh.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.mesh.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.mesh.bulk.len()
    }
    fn n_boundary(&self) -> usize {
        self.mesh.boundary.len()
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let i = e.0 as usize;
        let el = if i < self.mesh.bulk.len() {
            &self.mesh.bulk[i]
        } else {
            &self.mesh.boundary[i - self.mesh.bulk.len()]
        };
        ElementRef {
            etype: ElementType(el.etype.0),
            nodes: &self.nodes[i],
            tag: el.tag,
            parents: [None, None],
        }
    }
}

pub struct TestAssembly;

impl AssemblyBackend for TestAssembly {
    fn build_perm(&self, mesh: &dyn MeshView, active: &[ElemId], _opt: bool) -> Perm {
        let mut perm = vec![NO_DOF; mesh.n_nodes()];
        for &e in active {
            for &n in mesh.element(e).nodes {
                perm[n as usize] = 0;
            }
        }
        let mut next = 0;
        for p in perm.iter_mut().filter(|p| **p == 0) {
            *p = next;
            next += 1;
        }
        Perm {
            perm,
            n_slots: next as usize,
        }
    }

    fn create_matrix(
        &self,
        mesh: &dyn MeshView,
        active: &[ElemId],
        perm: &Perm,
        dofs: usize,
    ) -> CsrMatrix<f64> {
        assert_eq!(dofs, 1);
        let n = perm.n_slots;
        let mut rows: Vec<std::collections::BTreeSet<u32>> = vec![Default::default(); n];
        for &e in active {
            let slots: Vec<u32> = mesh
                .element(e)
                .nodes
                .iter()
                .map(|&n| perm.perm[n as usize])
                .collect();
            for &i in &slots {
                rows[i as usize].extend(slots.iter().copied());
            }
        }
        let mut row_ptr = vec![0];
        let mut col_idx = Vec::new();
        let mut diag = Vec::new();
        for (i, r) in rows.iter().enumerate() {
            for &c in r {
                if c as usize == i {
                    diag.push(col_idx.len());
                }
                col_idx.push(c);
            }
            row_ptr.push(col_idx.len());
        }
        CsrMatrix {
            nrows: n,
            ncols: n,
            row_ptr,
            values: vec![0.0; col_idx.len()],
            col_idx,
            diag: Some(diag),
        }
    }

    fn add_local(&self, a: &mut CsrMatrix<f64>, rhs: &mut [f64], dof: &[u32], local: &LocalSystem) {
        let n = local.n;
        for p in 0..n {
            let r = dof[p] as usize;
            for q in 0..n {
                let (lo, hi) = (a.row_ptr[r], a.row_ptr[r + 1]);
                let k = lo + a.col_idx[lo..hi].binary_search(&dof[q]).unwrap();
                a.values[k] += local.stiff[p * n + q];
            }
            rhs[r] += local.force[p];
        }
    }

    /// Row zeroed except the diagonal, `b = diag * value` (non-symmetric Elmer).
    fn apply_dirichlet(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        fixed: &[(u32, f64)],
        _params: &ParamMap,
    ) {
        for &(dof, v) in fixed {
            let r = dof as usize;
            let d = a.diag.as_ref().unwrap()[r];
            let diag = a.values[d];
            for k in a.row_ptr[r]..a.row_ptr[r + 1] {
                if k != d {
                    a.values[k] = 0.0;
                }
            }
            rhs[r] = diag * v;
        }
    }
}

/// Backward Euler: `K += M/dt`, `f += M x_prev / dt`.
pub struct BackwardEuler;

impl TimeBackend for BackwardEuler {
    fn add_time_terms(
        &self,
        local: &mut LocalSystem,
        t: &TimeTerms<'_>,
    ) -> Result<(), SolverError> {
        let n = local.n;
        let mass = local.mass.as_ref().expect("transient element has mass");
        let prev = &t.prev[0];
        for p in 0..n {
            for q in 0..n {
                local.stiff[p * n + q] += mass[p * n + q] / t.dt;
                local.force[p] += mass[p * n + q] * prev[q] / t.dt;
            }
        }
        Ok(())
    }
}

#[derive(Default)]
pub struct DenseLu {
    n: usize,
    lu: Vec<f64>,
    piv: Vec<usize>,
}

impl DirectSolver<f64> for DenseLu {
    fn factor(&mut self, a: &CsrMatrix<f64>) -> Result<(), LinalgError> {
        let n = a.nrows;
        let mut lu = vec![0.0; n * n];
        for i in 0..n {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                lu[i * n + a.col_idx[k] as usize] = a.values[k];
            }
        }
        let mut piv = vec![0; n];
        for c in 0..n {
            let p = (c..n)
                .max_by(|&x, &y| lu[x * n + c].abs().total_cmp(&lu[y * n + c].abs()))
                .unwrap();
            if lu[p * n + c] == 0.0 {
                return Err(LinalgError::Singular(c));
            }
            piv[c] = p;
            if p != c {
                for j in 0..n {
                    lu.swap(c * n + j, p * n + j);
                }
            }
            for i in c + 1..n {
                let f = lu[i * n + c] / lu[c * n + c];
                lu[i * n + c] = f;
                for j in c + 1..n {
                    lu[i * n + j] -= f * lu[c * n + j];
                }
            }
        }
        *self = Self { n, lu, piv };
        Ok(())
    }

    fn solve_in_place(&self, b: &mut [f64]) -> Result<(), LinalgError> {
        let n = self.n;
        for c in 0..n {
            b.swap(c, self.piv[c]);
            for i in c + 1..n {
                b[i] -= self.lu[i * n + c] * b[c];
            }
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                b[i] -= self.lu[i * n + j] * b[j];
            }
            b[i] /= self.lu[i * n + i];
        }
        Ok(())
    }
}
