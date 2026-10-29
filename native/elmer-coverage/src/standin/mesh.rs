//! Elmer mesh DB reader (`mesh.nodes`, `mesh.elements`, `mesh.boundary`) as
//! a `MeshView`, standing in for S1's `elmer-mesh`. Boundary parents are not
//! resolved. Unlike the copy in `elmer-run`'s tests, a malformed file is an
//! error, not a panic.

use std::path::Path;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView, ParamMap, Params};

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
    Ok(Elem { tag: num(row, 1, "tag")?, etype: ElementType(num(row, first, "element type")?), nodes })
}

impl FileMesh {
    /// Reads `dir/mesh.*`; `dim` is the coordinate system dimension and
    /// `scale` the `Coordinate Scaling` factor (MeshUtils applies a scalar
    /// factor at load).
    pub fn read(dir: &Path, dim: usize, scale: f64) -> Result<Self, String> {
        let coords = rows(&dir.join("mesh.nodes"))?
            .iter()
            .map(|l| Ok([scale * num::<f64>(l, 2, "x")?, scale * num::<f64>(l, 3, "y")?, scale * num::<f64>(l, 4, "z")?]))
            .collect::<Result<Vec<_>, String>>()?;
        let mut elems = rows(&dir.join("mesh.elements"))?.iter().map(|r| elem(r, 2)).collect::<Result<Vec<_>, _>>()?;
        let n_bulk = elems.len();
        elems.extend(rows(&dir.join("mesh.boundary"))?.iter().map(|r| elem(r, 4)).collect::<Result<Vec<_>, _>>()?);
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
        ElementRef { etype: el.etype, nodes: &el.nodes, tag: el.tag, parents: [None, None] }
    }
}

/// Dimension of the `Coordinate System` (`Cartesian 2D`, `Axi Symmetric`, ...).
fn dimension(simulation: &ParamMap) -> Result<usize, String> {
    let cs = simulation.string("coordinate system").unwrap_or("cartesian").to_ascii_lowercase();
    match cs.as_str() {
        "cartesian 1d" => Ok(1),
        "cartesian 2d" | "axi symmetric" | "cylindric symmetric" | "cylindric" => Ok(2),
        "cartesian" | "cartesian 3d" => Ok(3),
        other => Err(format!("stand-in: Coordinate System = {other}")),
    }
}

/// `elmer-run`'s `MeshLoader`: the mesh DB in `dir` with `Coordinate
/// Scaling` applied, refusing element types `TestBasis` cannot evaluate.
pub fn load(dir: &Path, simulation: &ParamMap) -> Result<Box<dyn MeshView>, String> {
    let scale = match simulation.get("coordinate scaling") {
        None => 1.0,
        Some(_) => simulation.real("coordinate scaling").ok_or("Coordinate Scaling is not a constant")?,
    };
    let mesh = FileMesh::read(dir, dimension(simulation)?, scale)?;
    if let Some(t) = mesh.element_types().into_iter().find(|t| !super::SUPPORTED_ELEMENTS.contains(t)) {
        return Err(format!("stand-in: element type {t} (elmer-fe not merged)"));
    }
    Ok(Box::new(mesh))
}
