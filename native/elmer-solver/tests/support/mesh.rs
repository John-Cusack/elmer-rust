//! Elmer mesh DB reader (`mesh.nodes`, `mesh.elements`, `mesh.boundary`) as
//! a `MeshView`, with the boundary elements' parents.

use std::path::Path;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

struct Elem {
    etype: ElementType,
    nodes: Vec<u32>,
    tag: i32,
    parents: [Option<ElemId>; 2],
}

pub struct FileMesh {
    dim: usize,
    coords: Vec<[f64; 3]>,
    /// Bulk elements, then boundary elements.
    elems: Vec<Elem>,
    n_bulk: usize,
}

fn rows(path: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .map(|l| l.split_whitespace().map(str::to_string).collect::<Vec<_>>())
        .filter(|l| !l.is_empty())
        .collect()
}

/// `first` is the column of the element type; node ids follow it.
fn elem(row: &[String], first: usize) -> Elem {
    Elem {
        tag: row[1].parse().unwrap(),
        etype: ElementType(row[first].parse().unwrap()),
        nodes: row[first + 1..].iter().map(|s| s.parse::<u32>().unwrap() - 1).collect(),
        parents: [None, None],
    }
}

impl FileMesh {
    /// Reads `dir/mesh.*`; `dim` is the coordinate system dimension.
    pub fn read(dir: &Path, dim: usize) -> Self {
        let coords = rows(&dir.join("mesh.nodes"))
            .iter()
            .map(|l| [l[2].parse().unwrap(), l[3].parse().unwrap(), l[4].parse().unwrap()])
            .collect();
        let bulk = rows(&dir.join("mesh.elements"));
        let index: std::collections::HashMap<i64, u32> =
            bulk.iter().enumerate().map(|(i, r)| (r[0].parse().unwrap(), i as u32)).collect();
        let mut elems: Vec<Elem> = bulk.iter().map(|r| elem(r, 2)).collect();
        let n_bulk = elems.len();
        elems.extend(rows(&dir.join("mesh.boundary")).iter().map(|r| {
            let mut e = elem(r, 4);
            let parent = |c: &String| c.parse::<i64>().ok().and_then(|id| index.get(&id)).map(|&i| ElemId(i));
            e.parents = [parent(&r[2]), parent(&r[3])];
            e
        }));
        FileMesh { dim, coords, elems, n_bulk }
    }
}

impl MeshView for FileMesh {
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
        let el = &self.elems[e.0 as usize];
        ElementRef { etype: el.etype, nodes: &el.nodes, tag: el.tag, parents: el.parents }
    }
}
