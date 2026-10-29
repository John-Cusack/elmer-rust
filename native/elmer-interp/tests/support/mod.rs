//! Test-local stand-ins for `elmer-mesh` and `elmer-fe` (not merged yet).
#![allow(dead_code)]

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementRef, ElementType, MeshView};
use elmer_interp::shape::Shape;

pub struct TestMesh {
    pub dim: usize,
    pub coords: Vec<[f64; 3]>,
    /// (type code, tag, nodes); bulk first then boundary.
    pub elems: Vec<(u16, i32, Vec<u32>)>,
    pub n_bulk: usize,
}

impl TestMesh {
    pub fn new(dim: usize, coords: Vec<[f64; 3]>, bulk: Vec<(u16, i32, Vec<u32>)>) -> Self {
        let n_bulk = bulk.len();
        TestMesh {
            dim,
            coords,
            elems: bulk,
            n_bulk,
        }
    }
}

impl MeshView for TestMesh {
    fn dim(&self) -> usize {
        self.dim
    }
    fn n_nodes(&self) -> usize {
        self.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.n_bulk
    }
    fn n_boundary(&self) -> usize {
        self.elems.len() - self.n_bulk
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let (code, tag, nodes) = &self.elems[e.0 as usize];
        ElementRef {
            etype: ElementType(*code),
            nodes,
            tag: *tag,
            parents: [None, None],
        }
    }
}

/// Structured triangle mesh of `[x0,x1] x [y0,y1]`, `nx` by `ny` squares each
/// cut along the same diagonal into two 303s.
pub fn tri_grid(x0: f64, x1: f64, y0: f64, y1: f64, nx: usize, ny: usize) -> TestMesh {
    let mut coords = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            coords.push([
                x0 + (x1 - x0) * i as f64 / nx as f64,
                y0 + (y1 - y0) * j as f64 / ny as f64,
                0.0,
            ]);
        }
    }
    let id = |i: usize, j: usize| (j * (nx + 1) + i) as u32;
    let mut el = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            el.push((303, 1, vec![id(i, j), id(i + 1, j), id(i + 1, j + 1)]));
            el.push((303, 1, vec![id(i, j), id(i + 1, j + 1), id(i, j + 1)]));
        }
    }
    TestMesh::new(2, coords, el)
}

pub fn quad_grid(x0: f64, x1: f64, y0: f64, y1: f64, nx: usize, ny: usize) -> TestMesh {
    let t = tri_grid(x0, x1, y0, y1, nx, ny);
    let id = |i: usize, j: usize| (j * (nx + 1) + i) as u32;
    let mut el = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            el.push((
                404,
                1,
                vec![id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)],
            ));
        }
    }
    TestMesh::new(2, t.coords, el)
}

/// Brick grid of 808s on the box.
pub fn hex_grid(lo: [f64; 3], hi: [f64; 3], n: [usize; 3]) -> TestMesh {
    let mut coords = Vec::new();
    for k in 0..=n[2] {
        for j in 0..=n[1] {
            for i in 0..=n[0] {
                coords.push([
                    lo[0] + (hi[0] - lo[0]) * i as f64 / n[0] as f64,
                    lo[1] + (hi[1] - lo[1]) * j as f64 / n[1] as f64,
                    lo[2] + (hi[2] - lo[2]) * k as f64 / n[2] as f64,
                ]);
            }
        }
    }
    let id =
        |i: usize, j: usize, k: usize| (k * (n[1] + 1) * (n[0] + 1) + j * (n[0] + 1) + i) as u32;
    let mut el = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                el.push((
                    808,
                    1,
                    vec![
                        id(i, j, k),
                        id(i + 1, j, k),
                        id(i + 1, j + 1, k),
                        id(i, j + 1, k),
                        id(i, j, k + 1),
                        id(i + 1, j, k + 1),
                        id(i + 1, j + 1, k + 1),
                        id(i, j + 1, k + 1),
                    ],
                ));
            }
        }
    }
    TestMesh::new(3, coords, el)
}

/// The hex grid with every brick cut into 6 tets around the main diagonal.
pub fn tet_grid(lo: [f64; 3], hi: [f64; 3], n: [usize; 3]) -> TestMesh {
    let h = hex_grid(lo, hi, n);
    let mut el = Vec::new();
    for (_, tag, v) in &h.elems {
        for t in [
            [0, 1, 2, 6],
            [0, 2, 3, 6],
            [0, 3, 7, 6],
            [0, 7, 4, 6],
            [0, 4, 5, 6],
            [0, 5, 1, 6],
        ] {
            el.push((504, *tag, t.iter().map(|&i| v[i]).collect()));
        }
    }
    TestMesh::new(3, h.coords, el)
}

/// Thin `ElementBasis` over the crate's nodal tables, with physical
/// gradients from the Jacobian.
pub struct NodalBasis;

impl ElementBasis for NodalBasis {
    fn eval(
        &self,
        etype: ElementType,
        xyz: &[[f64; 3]],
        uvw: [f64; 3],
        out: &mut BasisEval,
    ) -> Result<(), String> {
        let sh = Shape::get(etype.0).ok_or_else(|| format!("no table for {}", etype.0))?;
        let n = sh.n_nodes();
        out.basis.resize(n, 0.0);
        out.dbasis_dx.resize(n, [0.0; 3]);
        sh.basis(uvw, &mut out.basis);
        let mut g = vec![[0.0; 3]; n];
        sh.dbasis(uvw, &mut g);
        let d = sh.dim();
        let mut j = [[0.0; 3]; 3];
        for r in 0..3 {
            for c in 0..d {
                j[r][c] = (0..n).map(|k| g[k][c] * xyz[k][r]).sum();
            }
        }
        if d == 3 {
            let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
                - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
                + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
            out.det_j = det;
        } else if d == 2 {
            out.det_j = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        } else {
            out.det_j = j[0][0];
        }
        Ok(())
    }
}

pub fn assert_close(a: f64, b: f64, tol: f64, what: &str) {
    assert!((a - b).abs() <= tol * (1.0 + b.abs()), "{what}: {a} vs {b}");
}

/// Read an Elmer mesh directory (`mesh.nodes`, `mesh.elements`); boundary
/// elements are ignored. A thin stand-in for `elmer-mesh`.
pub fn read_elmer_mesh(dir: &std::path::Path, dim: usize) -> TestMesh {
    let nodes = std::fs::read_to_string(dir.join("mesh.nodes")).unwrap();
    let coords: Vec<[f64; 3]> = nodes
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let t: Vec<&str> = l.split_whitespace().collect();
            [
                t[2].parse().unwrap(),
                t[3].parse().unwrap(),
                t[4].parse().unwrap(),
            ]
        })
        .collect();
    let elems = std::fs::read_to_string(dir.join("mesh.elements")).unwrap();
    let bulk = elems
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let t: Vec<&str> = l.split_whitespace().collect();
            let nodes = t[3..]
                .iter()
                .map(|n| n.parse::<u32>().unwrap() - 1)
                .collect();
            (t[2].parse().unwrap(), t[1].parse().unwrap(), nodes)
        })
        .collect();
    TestMesh::new(dim, coords, bulk)
}
