//! Test-only readers for Elmer mesh directories and `Linear System Save`
//! dumps. Not a mesh crate: just enough to drive `elmer-assembly` against
//! pristine captures.
#![allow(dead_code)]

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

pub fn data_dir(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(case)
}

/// Non-comment lines of a fixture file.
pub fn lines(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(str::to_owned)
        .collect()
}

struct Elem {
    etype: ElementType,
    nodes: Vec<u32>,
    tag: i32,
    parents: [Option<ElemId>; 2],
}

/// An Elmer mesh directory (`mesh.nodes`, `mesh.elements`, `mesh.boundary`),
/// nodes and elements in file order as ElmerSolver loads them.
pub struct TestMesh {
    dim: usize,
    coords: Vec<[f64; 3]>,
    bulk: Vec<Elem>,
    boundary: Vec<Elem>,
}

impl TestMesh {
    pub fn read(dir: &Path, dim: usize) -> Self {
        let mut node_index = HashMap::new();
        let mut coords = Vec::new();
        for l in lines(&dir.join("mesh.nodes")) {
            let f: Vec<&str> = l.split_whitespace().collect();
            node_index.insert(f[0].parse::<u64>().unwrap(), coords.len() as u32);
            coords.push([f[2], f[3], f[4]].map(|s| s.parse::<f64>().unwrap()));
        }
        let node = |s: &str| node_index[&s.parse::<u64>().unwrap()];

        let mut elem_index = HashMap::new();
        let mut bulk = Vec::new();
        for l in lines(&dir.join("mesh.elements")) {
            let f: Vec<&str> = l.split_whitespace().collect();
            elem_index.insert(f[0].parse::<u64>().unwrap(), bulk.len() as u32);
            let etype = ElementType(f[2].parse().unwrap());
            let nodes: Vec<u32> = f[3..].iter().map(|s| node(s)).collect();
            assert_eq!(nodes.len(), etype.n_nodes());
            bulk.push(Elem {
                etype,
                nodes,
                tag: f[1].parse().unwrap(),
                parents: [None, None],
            });
        }
        let parent = |s: &str| match s.parse::<u64>().unwrap() {
            0 => None,
            id => Some(ElemId(elem_index[&id])),
        };
        let mut boundary = Vec::new();
        for l in lines(&dir.join("mesh.boundary")) {
            let f: Vec<&str> = l.split_whitespace().collect();
            let etype = ElementType(f[4].parse().unwrap());
            boundary.push(Elem {
                etype,
                nodes: f[5..].iter().map(|s| node(s)).collect(),
                tag: f[1].parse().unwrap(),
                parents: [parent(f[2]), parent(f[3])],
            });
        }
        Self {
            dim,
            coords,
            bulk,
            boundary,
        }
    }

    pub fn is_boundary(&self, e: ElemId) -> bool {
        e.0 as usize >= self.bulk.len()
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
        self.bulk.len()
    }
    fn n_boundary(&self) -> usize {
        self.boundary.len()
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let i = e.0 as usize;
        let el = if i < self.bulk.len() {
            &self.bulk[i]
        } else {
            &self.boundary[i - self.bulk.len()]
        };
        ElementRef {
            etype: el.etype,
            nodes: &el.nodes,
            tag: el.tag,
            parents: el.parents,
        }
    }
}

/// `Linear System Save` matrix dump: `(row, col, value)`, 0-based, in file
/// (= CSR) order.
pub fn read_matrix(path: &Path) -> Vec<(usize, usize, f64)> {
    lines(path)
        .iter()
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            (
                f[0].parse::<usize>().unwrap() - 1,
                f[1].parse::<usize>().unwrap() - 1,
                f[2].parse::<f64>().unwrap(),
            )
        })
        .collect()
}

/// `PrintRHS` dump. `None` marks Elmer's literal `0.0`, written for any
/// value with `|v| <= tiny` (sign and subnormals not recorded).
pub fn read_rhs(path: &Path) -> Vec<Option<f64>> {
    lines(path)
        .iter()
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            match f[1] {
                "0.0" => None,
                v => Some(v.parse::<f64>().unwrap()),
            }
        })
        .collect()
}

/// `Linear System Save Perm` dump: node -> 1-based slot, 0 = absent.
pub fn read_perm(path: &Path) -> Vec<u32> {
    lines(path)
        .iter()
        .map(|l| l.split_whitespace().nth(1).unwrap().parse().unwrap())
        .collect()
}
