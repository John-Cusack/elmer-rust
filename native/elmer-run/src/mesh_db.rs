//! Elmer mesh DB reader (`mesh.nodes`, `mesh.elements`, `mesh.boundary`) as
//! a `MeshView`, standing in for S1's `elmer-mesh`, with boundary parents
//! resolved from `mesh.boundary`. The `elmer-run` binary and
//! `elmer-coverage` read meshes with it. Unlike the copy in `elmer-run`'s
//! tests, a malformed file is an error, not a panic.
#![allow(clippy::needless_range_loop)]

use std::path::Path;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView, ParamMap, Params, Value};

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
            // A parent id <= 0 (`-1` on point elements) is no parent.
            for k in 0..2 {
                let id = num::<i64>(&r, 2 + k, "parent")?;
                e.parents[k] = u64::try_from(id).ok().and_then(|id| index.get(&id)).map(|&i| ElemId(i));
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

/// `MapCoordinates` (MeshUtils): `Coordinate Mapping(3) = i j k` makes
/// the coordinate vectors `(x, y, z)` the mesh file's `i`th, `j`th and `k`th.
fn map_coordinates(coords: &mut [[f64; 3]], simulation: &ParamMap) -> Result<(), String> {
    let Some(v) = simulation.get("coordinate mapping") else { return Ok(()) };
    let map = match v {
        Value::Integer(m) if m.len() == 3 && [1, 2, 3].iter().all(|c| m.contains(c)) => [m[0], m[1], m[2]],
        v => return Err(format!("inconsistent Coordinate Mapping: {v:?}")),
    };
    for c in coords {
        let old = *c;
        for k in 0..3 {
            c[k] = old[(map[k] - 1) as usize];
        }
    }
    Ok(())
}

/// The Simulation's `Coordinate Transformation` (ElmerSolver.F90:462,
/// `CoordinateTransformation(..., Irreversible = .TRUE.)`, MeshTransform.F90:
/// 1114): the nodes are mapped once, after loading, with `Coordinate
/// Transformation Radius` `r0`: "cartesian to polar|cylindrical" gives
/// `(r0 + sqrt(x^2 + y^2), c * atan2(y, x), z)`, "polar|cylindrical to
/// cartesian" `(cos(y / c) (x + r0), sin(y / c) (x + r0), z)`, with `c =
/// 180 / pi` under `Coordinate Transformation Use Degrees`. The callers
/// refuse it together with `Mesh Levels` and extrusion, which Elmer applies
/// before it.
fn coordinate_transformation(coords: &mut [[f64; 3]], simulation: &ParamMap) -> Result<(), String> {
    let Some(name) = simulation.string("coordinate transformation") else {
        return match simulation.get("coordinate transformation") {
            None => Ok(()),
            Some(v) => Err(format!("Coordinate Transformation {v:?} is not a string")),
        };
    };
    let mode = match elmer_core::normalize_key(name).as_str() {
        "cartesian to polar" | "cartesian to cylindrical" => 1,
        "polar to cartesian" | "cylindrical to cartesian" => -1,
        t => return Err(format!("CoordinateTransformation: Unknown transformation: {t}")),
    };
    let coeff = if simulation.logical("coordinate transformation use degrees") == Some(true) {
        180.0 / std::f64::consts::PI
    } else {
        1.0
    };
    let rad0 = simulation.real("coordinate transformation radius").unwrap_or(0.0);
    for c in coords {
        let r0 = *c;
        *c = if mode == 1 {
            [rad0 + (r0[0].powi(2) + r0[1].powi(2)).sqrt(), coeff * r0[1].atan2(r0[0]), r0[2]]
        } else {
            [(r0[1] / coeff).cos() * (r0[0] + rad0), (r0[1] / coeff).sin() * (r0[0] + rad0), r0[2]]
        };
    }
    Ok(())
}

/// `elmer-run`'s `MeshLoader`: the mesh DB in `dir` with `Coordinate
/// Scaling`, `Coordinate Mapping` and the Simulation's `Coordinate
/// Transformation` applied, refusing element types `elmer-fe` does not
/// define.
pub fn load(dir: &Path, simulation: &ParamMap) -> Result<Box<dyn MeshView>, String> {
    let scale = match simulation.get("coordinate scaling") {
        None => 1.0,
        Some(_) => simulation.real("coordinate scaling").ok_or("Coordinate Scaling is not a constant")?,
    };
    let mut mesh = FileMesh::read(dir, 3, scale)?;
    map_coordinates(&mut mesh.coords, simulation)?;
    coordinate_transformation(&mut mesh.coords, simulation)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn sim(entries: &[(&str, Value)]) -> ParamMap {
        let mut p = ParamMap::new();
        for (k, v) in entries {
            p.insert(k, v.clone());
        }
        p
    }

    #[test]
    fn mapping_permutes_the_coordinate_vectors() {
        let mut c = [[1.0, 2.0, 3.0]];
        map_coordinates(&mut c, &sim(&[("Coordinate Mapping", Value::Integer(vec![1, 3, 2]))])).unwrap();
        assert_eq!(c, [[1.0, 3.0, 2.0]]);
        assert!(map_coordinates(&mut c, &sim(&[("Coordinate Mapping", Value::Integer(vec![1, 1, 2]))])).is_err());
    }

    #[test]
    fn cylindrical_and_cartesian_transformations() {
        let half_pi = std::f64::consts::FRAC_PI_2;
        let r = |v: f64| Value::Real { rows: 1, cols: 1, data: vec![v] };
        let mut c = [[1.0, half_pi, 5.0]];
        let to_cart = sim(&[
            ("Coordinate Transformation", Value::Str("Cylindrical to Cartesian".into())),
            ("Coordinate Transformation Radius", r(1.0)),
        ]);
        coordinate_transformation(&mut c, &to_cart).unwrap();
        assert!((c[0][0] - 0.0).abs() < 1e-15 && c[0][1] == 2.0 && c[0][2] == 5.0);
        let mut c = [[0.0, 2.0, 5.0]];
        let to_cyl = sim(&[
            ("Coordinate Transformation", Value::Str("cartesian to cylindrical".into())),
            ("Coordinate Transformation Use Degrees", Value::Logical(true)),
        ]);
        coordinate_transformation(&mut c, &to_cyl).unwrap();
        assert_eq!(c, [[2.0, 90.0, 5.0]]);
        let bad = sim(&[("Coordinate Transformation", Value::Str("spherical".into()))]);
        assert!(coordinate_transformation(&mut c, &bad).is_err());
    }
}
