//! Elmer mesh files (mesh.header / .nodes / .elements / .boundary) as a
//! `MeshView`, keeping the boundary elements' parents.

use std::path::Path;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

pub struct Elem {
    pub etype: ElementType,
    /// 0-based node indices, Elmer local order.
    pub nodes: Vec<u32>,
    /// Body id (bulk) or boundary tag.
    pub tag: i32,
    pub parents: [Option<ElemId>; 2],
}

pub struct TestMesh {
    pub dim: usize,
    pub coords: Vec<[f64; 3]>,
    pub bulk: Vec<Elem>,
    pub boundary: Vec<Elem>,
}

fn lines(path: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .filter(|l: &Vec<String>| !l.is_empty())
        .collect()
}

impl TestMesh {
    /// Reads `dir/mesh.*`; element ids are assumed to be 1..n in file order.
    pub fn read(dir: &Path, dim: usize) -> Self {
        let coords = lines(&dir.join("mesh.nodes"))
            .iter()
            .map(|l| [l[2].parse().unwrap(), l[3].parse().unwrap(), l[4].parse().unwrap()])
            .collect();
        let nodes = |l: &[String]| l.iter().map(|s| s.parse::<u32>().unwrap() - 1).collect();
        let bulk = lines(&dir.join("mesh.elements"))
            .iter()
            .map(|l| Elem {
                tag: l[1].parse().unwrap(),
                etype: ElementType(l[2].parse().unwrap()),
                nodes: nodes(&l[3..]),
                parents: [None, None],
            })
            .collect();
        let parent = |s: &str| match s.parse::<u32>().unwrap() {
            0 => None,
            p => Some(ElemId(p - 1)),
        };
        let boundary = lines(&dir.join("mesh.boundary"))
            .iter()
            .map(|l| Elem {
                tag: l[1].parse().unwrap(),
                parents: [parent(&l[2]), parent(&l[3])],
                etype: ElementType(l[4].parse().unwrap()),
                nodes: nodes(&l[5..]),
            })
            .collect();
        Self { dim, coords, bulk, boundary }
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
