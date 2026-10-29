//! Elmer mesh DB reader (`mesh.nodes`, `mesh.elements`, `mesh.boundary`) as
//! a `MeshView`. Boundary parents are not resolved.

use std::path::Path;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView, ParamMap, Params, Value};

struct Elem {
    etype: ElementType,
    nodes: Vec<u32>,
    tag: i32,
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
    }
}

impl FileMesh {
    /// Reads `dir/mesh.*`; `dim` is the coordinate system dimension.
    pub fn read(dir: &Path, dim: usize) -> Self {
        let coords = rows(&dir.join("mesh.nodes"))
            .iter()
            .map(|l| [l[2].parse().unwrap(), l[3].parse().unwrap(), l[4].parse().unwrap()])
            .collect();
        let mut elems: Vec<Elem> = rows(&dir.join("mesh.elements")).iter().map(|r| elem(r, 2)).collect();
        let n_bulk = elems.len();
        elems.extend(rows(&dir.join("mesh.boundary")).iter().map(|r| elem(r, 4)));
        FileMesh { dim, coords, elems, n_bulk }
    }

    /// `Coordinate Scaling` applied at mesh load (MeshUtils `CoordinateScaling`
    /// for a scalar factor), standing in for `elmer-mesh`.
    pub fn scaled(mut self, factor: f64) -> Self {
        for c in &mut self.coords {
            for x in c.iter_mut() {
                *x *= factor;
            }
        }
        self
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
        ElementRef { etype: el.etype, nodes: &el.nodes, tag: el.tag, parents: [None, None] }
    }
}

/// [`elmer_run::MeshLoader`]: the mesh DB with its dimension taken from the
/// bulk elements and a scalar `Coordinate Scaling` applied.
pub fn load(dir: &Path, simulation: &ParamMap) -> Result<Box<dyn MeshView>, String> {
    let mesh = FileMesh::read(dir, 3);
    let dim = match mesh.elems[..mesh.n_bulk].iter().map(|e| e.etype.family()).max() {
        Some(2) => 1,
        Some(3 | 4) => 2,
        _ => 3,
    };
    let scale = match simulation.get("coordinate scaling") {
        None => 1.0,
        Some(Value::Real { data, .. }) if data.len() == 1 => data[0],
        Some(_) => return Err("Coordinate Scaling other than one real".into()),
    };
    Ok(Box::new(FileMesh { dim, ..mesh }.scaled(scale)))
}
