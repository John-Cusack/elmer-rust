//! Elmer mesh DB reader (`mesh.nodes`, `mesh.elements`, `mesh.boundary`) as a
//! `MeshView` with the boundary parents resolved, standing in for S1's
//! `elmer-mesh`.
#![allow(dead_code)]

use std::path::Path;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

type Elem = (u16, i32, Vec<u32>, [Option<ElemId>; 2]);

pub struct FileMesh {
    dim: usize,
    xyz: Vec<[f64; 3]>,
    elems: Vec<Elem>,
    n_bulk: usize,
}

fn rows(path: &Path) -> Vec<Vec<usize>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|x| x.parse().unwrap()).collect())
        .collect()
}

impl FileMesh {
    pub fn read(dir: &Path, dim: usize) -> Self {
        let xyz: Vec<[f64; 3]> = std::fs::read_to_string(dir.join("mesh.nodes"))
            .unwrap()
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| {
                let t: Vec<f64> = l.split_whitespace().skip(2).map(|x| x.parse().unwrap()).collect();
                [t[0], t[1], t[2]]
            })
            .collect();
        let mut elems: Vec<Elem> = rows(&dir.join("mesh.elements"))
            .iter()
            .map(|t| (t[2] as u16, t[1] as i32, t[3..].iter().map(|n| (n - 1) as u32).collect(), [None, None]))
            .collect();
        let n_bulk = elems.len();
        for t in rows(&dir.join("mesh.boundary")) {
            // id tag left right type nodes...
            let left = (t[2] != 0).then(|| ElemId((t[2] - 1) as u32));
            let right = (t[3] != 0).then(|| ElemId((t[3] - 1) as u32));
            elems.push((t[4] as u16, t[1] as i32, t[5..].iter().map(|n| (n - 1) as u32).collect(), [left, right]));
        }
        FileMesh { dim, xyz, elems, n_bulk }
    }
}

impl MeshView for FileMesh {
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
