//! Test-local stand-ins until `elmer-mesh` merges: an Elmer mesh DB reader
//! and the `.sif` boundary conditions as parameter maps.

#![allow(dead_code)]

use std::path::Path;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView, ParamMap, Params};

struct Elem {
    etype: ElementType,
    nodes: Vec<u32>,
    tag: i32,
}

/// `mesh.nodes`, `mesh.elements`, `mesh.boundary` with node ids `1..n` in
/// file order (as ElmerGrid writes them). Boundary parents are not resolved.
pub struct FileMesh {
    dim: usize,
    coords: Vec<[f64; 3]>,
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
        nodes: row[first + 1..]
            .iter()
            .map(|s| s.parse::<u32>().unwrap() - 1)
            .collect(),
    }
}

impl FileMesh {
    pub fn read(dir: &Path, dim: usize) -> Self {
        let coords = rows(&dir.join("mesh.nodes"))
            .iter()
            .map(|l| {
                [
                    l[2].parse().unwrap(),
                    l[3].parse().unwrap(),
                    l[4].parse().unwrap(),
                ]
            })
            .collect();
        let mut elems: Vec<Elem> = rows(&dir.join("mesh.elements"))
            .iter()
            .map(|r| elem(r, 2))
            .collect();
        let n_bulk = elems.len();
        elems.extend(rows(&dir.join("mesh.boundary")).iter().map(|r| elem(r, 4)));
        FileMesh {
            dim,
            coords,
            elems,
            n_bulk,
        }
    }
}

impl FileMesh {
    /// The dimension of the bulk elements (Elmer `MeshDim`).
    pub fn bulk_dim(&self) -> usize {
        self.elems[..self.n_bulk]
            .iter()
            .map(|e| match e.etype.0 / 100 {
                1 => 0,
                2 => 1,
                3 | 4 => 2,
                _ => 3,
            })
            .max()
            .unwrap_or(0)
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
        ElementRef {
            etype: el.etype,
            nodes: &el.nodes,
            tag: el.tag,
            parents: [None, None],
        }
    }
}

/// The parsed sections a test needs.
pub struct Case {
    pub simulation: ParamMap,
    pub solvers: Vec<ParamMap>,
    /// `Boundary Condition 1..n` at `0..n`.
    pub bcs: Vec<ParamMap>,
}

impl Case {
    pub fn read(sif: &Path) -> Self {
        let parsed = elmer_sif::parse_file(sif.to_str().unwrap()).unwrap();
        let m = parsed.sif.merged();
        let numbered = |s: &[elmer_sif::Section]| {
            let n = s.iter().map(|s| s.id as usize).max().unwrap_or(0);
            let mut out = vec![ParamMap::new(); n];
            for s in s {
                out[s.id as usize - 1] = elmer_sif::section_to_param_map(s);
            }
            out
        };
        Case {
            simulation: m
                .simulation
                .as_ref()
                .map(elmer_sif::section_to_param_map)
                .unwrap_or_default(),
            solvers: numbered(&m.solvers),
            bcs: numbered(&m.bcs),
        }
    }

    /// `BoundaryInfo % Constraint`: the first boundary condition whose
    /// `Target Boundaries` lists the tag.
    pub fn bc_of(&self, tag: i32) -> Option<usize> {
        self.bcs
            .iter()
            .position(|bc| match bc.get("target boundaries") {
                Some(elmer_core::Value::Integer(v)) => v.contains(&tag),
                _ => false,
            })
            .map(|i| i + 1)
    }

    /// `CoordinateSystemDimension()` for the explicit `Cartesian 1D/2D/3D`
    /// and axisymmetric systems; plain `Cartesian` takes the mesh dimension.
    pub fn coord_dim(&self, mesh_dim: usize) -> usize {
        let s = self
            .simulation
            .string("coordinate system")
            .unwrap_or("")
            .to_ascii_lowercase();
        if s.contains("2d") || s.contains("axi") {
            2
        } else if s.contains("1d") {
            1
        } else if s.contains("3d") {
            3
        } else {
            mesh_dim
        }
    }
}
