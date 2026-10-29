//! Elmer mesh-file reader (mesh.header / .nodes / .elements / .boundary).

use std::path::Path;

use elmer_core::ElementType;

#[derive(Clone, Debug)]
pub struct Elem {
    pub etype: ElementType,
    /// 0-based node indices, Elmer local order.
    pub nodes: Vec<usize>,
    /// Body id (bulk) or boundary tag.
    pub tag: i32,
    /// Bulk elements on either side of a boundary element (0-based), as
    /// `mesh.boundary` lists them; empty for bulk elements.
    pub parents: [Option<usize>; 2],
}

#[derive(Clone, Debug, Default)]
pub struct Mesh {
    pub dim: usize,
    pub coords: Vec<[f64; 3]>,
    pub bulk: Vec<Elem>,
    pub boundary: Vec<Elem>,
}

impl Mesh {
    pub fn xyz(&self, e: &Elem) -> Vec<[f64; 3]> {
        e.nodes.iter().map(|&i| self.coords[i]).collect()
    }
}

fn lines(path: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .filter(|l: &Vec<String>| !l.is_empty())
        .collect()
}

/// Reads `dir/mesh.*`. `dim` is the coordinate-system dimension.
pub fn read_elmer_mesh(dir: &Path, dim: usize) -> Mesh {
    let nodes = lines(&dir.join("mesh.nodes"));
    let coords = nodes
        .iter()
        .map(|l| {
            [
                l[2].parse().unwrap(),
                l[3].parse().unwrap(),
                l[4].parse().unwrap(),
            ]
        })
        .collect();
    let bulk = lines(&dir.join("mesh.elements"))
        .iter()
        .map(|l| Elem {
            tag: l[1].parse().unwrap(),
            etype: ElementType(l[2].parse().unwrap()),
            nodes: l[3..]
                .iter()
                .map(|s| s.parse::<usize>().unwrap() - 1)
                .collect(),
            parents: [None, None],
        })
        .collect();
    let boundary = lines(&dir.join("mesh.boundary"))
        .iter()
        .map(|l| Elem {
            tag: l[1].parse().unwrap(),
            parents: [2, 3].map(|k| l[k].parse::<usize>().unwrap().checked_sub(1)),
            etype: ElementType(l[4].parse().unwrap()),
            nodes: l[5..]
                .iter()
                .map(|s| s.parse::<usize>().unwrap() - 1)
                .collect(),
        })
        .collect();
    Mesh {
        dim,
        coords,
        bulk,
        boundary,
    }
}

