//! Elmer mesh DB reader (`mesh.nodes`, `mesh.elements`, `mesh.boundary`) as
//! a `MeshView`, standing in for S1's `elmer-mesh`, with boundary parents
//! resolved from `mesh.boundary`. Unlike the copy in `elmer-run`'s tests, a malformed file is an
//! error, not a panic.

use std::path::Path;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView, ParamMap, Params};

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

fn rows(path: &Path) -> Result<Vec<Vec<String>>, String> {
    Ok(std::fs::read_to_string(path)
        .map_err(|e| format!("{}: {e}", path.display()))?
        .lines()
        .map(|l| l.split_whitespace().map(str::to_string).collect::<Vec<_>>())
        .filter(|l| !l.is_empty())
        .collect())
}

fn num<T: std::str::FromStr>(row: &[String], i: usize, what: &str) -> Result<T, String> {
    row.get(i).and_then(|s| s.parse().ok()).ok_or_else(|| format!("mesh: bad {what} in `{}`", row.join(" ")))
}

/// `first` is the column of the element type; node ids follow it.
fn elem(row: &[String], first: usize) -> Result<Elem, String> {
    let nodes = row
        .get(first + 1..)
        .unwrap_or(&[])
        .iter()
        .map(|s| s.parse::<u32>().ok().and_then(|n| n.checked_sub(1)).ok_or_else(|| format!("mesh: bad node id `{s}`")))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Elem { tag: num(row, 1, "tag")?, etype: ElementType(num(row, first, "element type")?), nodes, parents: [None, None] })
}

impl FileMesh {
    /// Reads `dir/mesh.*`; `dim` is the mesh dimension and
    /// `scale` the `Coordinate Scaling` factor (MeshUtils applies a scalar
    /// factor at load).
    pub fn read(dir: &Path, dim: usize, scale: f64) -> Result<Self, String> {
        let coords = rows(&dir.join("mesh.nodes"))?
            .iter()
            .map(|l| Ok([scale * num::<f64>(l, 2, "x")?, scale * num::<f64>(l, 3, "y")?, scale * num::<f64>(l, 4, "z")?]))
            .collect::<Result<Vec<_>, String>>()?;
        let bulk_rows = rows(&dir.join("mesh.elements"))?;
        // Parent columns hold mesh.elements ids, which need not be 1..n.
        let mut index = std::collections::HashMap::new();
        for (i, r) in bulk_rows.iter().enumerate() {
            index.insert(num::<u64>(r, 0, "element id")?, i as u32);
        }
        let mut elems = bulk_rows.iter().map(|r| elem(r, 2)).collect::<Result<Vec<_>, _>>()?;
        let n_bulk = elems.len();
        for r in rows(&dir.join("mesh.boundary"))? {
            let mut e = elem(&r, 4)?;
            for k in 0..2 {
                e.parents[k] = index.get(&num::<u64>(&r, 2 + k, "parent")?).map(|&i| ElemId(i));
            }
            elems.push(e);
        }
        if elems.iter().any(|e| e.nodes.iter().any(|&n| n as usize >= coords.len())) {
            return Err("mesh: element node id beyond mesh.nodes".into());
        }
        Ok(FileMesh { dim, coords, elems, n_bulk })
    }

    /// Element types in the mesh, bulk and boundary.
    pub fn element_types(&self) -> std::collections::BTreeSet<u16> {
        self.elems.iter().map(|e| e.etype.0).collect()
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

/// `elmer-run`'s `MeshLoader`: the mesh DB in `dir` with `Coordinate
/// Scaling` applied, refusing element types `elmer-fe` does not define.
pub fn load(dir: &Path, simulation: &ParamMap) -> Result<Box<dyn MeshView>, String> {
    let scale = match simulation.get("coordinate scaling") {
        None => 1.0,
        Some(_) => simulation.real("coordinate scaling").ok_or("Coordinate Scaling is not a constant")?,
    };
    let mut mesh = FileMesh::read(dir, 3, scale)?;
    // The mesh dimension (Elmer's `Mesh % MaxDim`, inherited by a plain
    // "cartesian" coordinate system) from the bulk elements.
    mesh.dim = match mesh.elems[..mesh.n_bulk].iter().map(|e| e.etype.family()).max() {
        Some(1 | 2) => 1,
        Some(3 | 4) => 2,
        _ => 3,
    };
    if let Some(t) = mesh.element_types().into_iter().find(|&t| elmer_fe::element_type(t).is_err()) {
        return Err(format!("element type {t} is not an elmer-fe Lagrange type"));
    }
    Ok(Box::new(mesh))
}
