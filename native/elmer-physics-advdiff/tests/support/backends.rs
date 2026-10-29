//! A `MeshView` over the test mesh (until `elmer-mesh` merges) and a dense
//! direct solver standing in for UMFPACK.

use elmer_core::{CsrMatrix, DirectSolver, ElemId, ElementRef, ElementType, LinalgError, MeshView};

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

/// Dense LU with partial pivoting, standing in for UMFPACK.
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
