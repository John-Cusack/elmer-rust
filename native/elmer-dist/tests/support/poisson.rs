//! Test-local stand-in for elmer-fe/elmer-assembly (S1 lane): `-div grad u = 1`
//! on a partitioned 2D mesh, `u = 0` where x or y is 0 or 1. P1 triangles
//! (303) and bilinear quads (404, 2x2 Gauss). Each rank assembles only the
//! bulk elements it owns, so shared rows hold partial sums.
//!
//! Dirichlet rows follow `EnforceDirichletConditions` in parallel: every rank
//! holding the dof zeroes its row and column and puts 1 on the diagonal (the
//! summed row has the number of holders there), rhs 0.
#![allow(dead_code, clippy::needless_range_loop)]

use std::path::Path;

use elmer_core::{CsrMatrix, MeshView};
use elmer_dist::{DofLayout, Exchange};
use elmer_mpi::PartitionMesh;

pub struct LocalProblem {
    pub mesh: PartitionMesh,
    pub layout: DofLayout,
    pub a: CsrMatrix<f64>,
    pub b: Vec<f64>,
}

pub fn data(sub: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(sub)
}

/// Element stiffness and load (`f = 1`) of one 2D element.
fn element(xy: &[[f64; 3]], code: u16) -> (Vec<f64>, Vec<f64>) {
    let n = xy.len();
    let mut k = vec![0.0; n * n];
    let mut f = vec![0.0; n];
    // Reference points, weights, basis and reference gradients.
    let pts: Vec<([f64; 2], f64)> = match code {
        303 => vec![([1.0 / 3.0, 1.0 / 3.0], 0.5)],
        404 => {
            let g = 1.0 / 3f64.sqrt();
            vec![
                ([-g, -g], 1.0),
                ([g, -g], 1.0),
                ([g, g], 1.0),
                ([-g, g], 1.0),
            ]
        }
        _ => panic!("element type {code} not in the Poisson stand-in"),
    };
    for (uv, w) in pts {
        let (u, v) = (uv[0], uv[1]);
        let (phi, dphi): (Vec<f64>, Vec<[f64; 2]>) = if code == 303 {
            (
                vec![1.0 - u - v, u, v],
                vec![[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]],
            )
        } else {
            let s = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
            (
                s.iter()
                    .map(|&(a, b)| 0.25 * (1.0 + a * u) * (1.0 + b * v))
                    .collect(),
                s.iter()
                    .map(|&(a, b)| [0.25 * a * (1.0 + b * v), 0.25 * b * (1.0 + a * u)])
                    .collect(),
            )
        };
        let mut j = [[0.0; 2]; 2];
        for (p, d) in xy.iter().zip(&dphi) {
            for r in 0..2 {
                for c in 0..2 {
                    j[r][c] += d[r] * p[c];
                }
            }
        }
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let inv = [
            [j[1][1] / det, -j[0][1] / det],
            [-j[1][0] / det, j[0][0] / det],
        ];
        let grad: Vec<[f64; 2]> = dphi
            .iter()
            .map(|d| {
                [
                    inv[0][0] * d[0] + inv[0][1] * d[1],
                    inv[1][0] * d[0] + inv[1][1] * d[1],
                ]
            })
            .collect();
        let wd = w * det.abs();
        for a in 0..n {
            f[a] += wd * phi[a];
            for b in 0..n {
                k[a * n + b] += wd * (grad[a][0] * grad[b][0] + grad[a][1] * grad[b][1]);
            }
        }
    }
    (k, f)
}

fn on_boundary(p: [f64; 3]) -> bool {
    let e = 1e-12;
    p[0].abs() < e || p[1].abs() < e || (p[0] - 1.0).abs() < e || (p[1] - 1.0).abs() < e
}

/// Read rank `rank` of `<mesh_dir>/partitioning.<size>` and assemble its part.
pub fn assemble(comm: &dyn Exchange, mesh_dir: &Path) -> LocalProblem {
    let mesh = PartitionMesh::read(mesh_dir, comm.rank(), comm.size()).expect("read partition");
    let layout = DofLayout::nodal(comm, &mesh, 1).expect("layout");
    let n = mesh.n_nodes();
    let fixed: Vec<bool> = (0..n as u32).map(|i| on_boundary(mesh.coord(i))).collect();

    let mut t: Vec<(usize, usize, f64)> = Vec::new();
    let mut b = vec![0.0; n];
    for e in mesh
        .bulk_elements()
        .iter()
        .filter(|e| e.owner == comm.rank())
    {
        let xy: Vec<[f64; 3]> = e.nodes.iter().map(|&i| mesh.coord(i)).collect();
        let (k, f) = element(&xy, e.etype.0);
        let m = e.nodes.len();
        for a in 0..m {
            let i = e.nodes[a] as usize;
            if fixed[i] {
                continue;
            }
            b[i] += f[a];
            for c in 0..m {
                let j = e.nodes[c] as usize;
                if !fixed[j] {
                    t.push((i, j, k[a * m + c]));
                }
            }
        }
    }
    for (i, _) in fixed.iter().enumerate().filter(|(_, &f)| f) {
        t.push((i, i, 1.0));
    }
    // Every node gets a diagonal slot, even one no owned element touches.
    t.extend((0..n).map(|i| (i, i, 0.0)));
    let a = elmer_sparse::from_triplets(n, n, t.into_iter().map(|(i, j, v)| (i, j as u32, v)));
    LocalProblem { mesh, layout, a, b }
}
