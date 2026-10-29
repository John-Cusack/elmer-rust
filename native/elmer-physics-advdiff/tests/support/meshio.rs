//! Elmer mesh-file reader (mesh.header / .nodes / .elements / .boundary) and
//! structured unit-square meshes for the convergence tests.

use std::collections::HashMap;
use std::path::Path;

use elmer_core::ElementType;

#[derive(Clone, Debug)]
pub struct Elem {
    pub etype: ElementType,
    /// 0-based node indices, Elmer local order.
    pub nodes: Vec<usize>,
    /// Body id (bulk) or boundary tag.
    pub tag: i32,
}

#[derive(Clone, Debug, Default)]
pub struct Mesh {
    pub dim: usize,
    pub coords: Vec<[f64; 3]>,
    pub bulk: Vec<Elem>,
    pub boundary: Vec<Elem>,
    /// Parent bulk elements (0-based) of each boundary element, from
    /// `mesh.boundary`; empty for the generated meshes.
    pub parents: Vec<[Option<u32>; 2]>,
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
        })
        .collect();
    let parent = |s: &str| match s.parse::<u32>().unwrap() {
        0 => None,
        id => Some(id - 1),
    };
    let parents = lines(&dir.join("mesh.boundary"))
        .iter()
        .map(|l| [parent(&l[2]), parent(&l[3])])
        .collect();
    let boundary = lines(&dir.join("mesh.boundary"))
        .iter()
        .map(|l| Elem {
            tag: l[1].parse().unwrap(),
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
        parents,
    }
}

/// `n x n` cells on the unit square. `quad`: 404 cells, else two 303 per cell.
/// `second_order`: add mid-edge nodes (306, 408). Boundary lines are tagged
/// 1 bottom, 2 right, 3 top, 4 left (202 / 203).
pub fn unit_square(n: usize, quad: bool, second_order: bool) -> Mesh {
    let idx = |i: usize, j: usize| j * (n + 1) + i;
    let mut m = Mesh {
        dim: 2,
        ..Default::default()
    };
    for j in 0..=n {
        for i in 0..=n {
            m.coords
                .push([i as f64 / n as f64, j as f64 / n as f64, 0.0]);
        }
    }
    let mut bulk = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1));
            if quad {
                bulk.push((404, vec![a, b, c, d]));
            } else {
                bulk.push((303, vec![a, b, c]));
                bulk.push((303, vec![a, c, d]));
            }
        }
    }
    let mut bnd = Vec::new();
    for i in 0..n {
        bnd.push((1, vec![idx(i, 0), idx(i + 1, 0)]));
        bnd.push((2, vec![idx(n, i), idx(n, i + 1)]));
        bnd.push((3, vec![idx(i + 1, n), idx(i, n)]));
        bnd.push((4, vec![idx(0, i + 1), idx(0, i)]));
    }
    let mut mids: HashMap<(usize, usize), usize> = HashMap::new();
    let mut mid = |m: &mut Mesh, a: usize, b: usize| -> usize {
        let key = (a.min(b), a.max(b));
        *mids.entry(key).or_insert_with(|| {
            let (p, q) = (m.coords[a], m.coords[b]);
            m.coords.push([
                (p[0] + q[0]) / 2.0,
                (p[1] + q[1]) / 2.0,
                (p[2] + q[2]) / 2.0,
            ]);
            m.coords.len() - 1
        })
    };
    for (t, nodes) in bulk {
        let mut nodes = nodes;
        let mut et = t;
        if second_order {
            let k = nodes.len();
            for e in 0..k {
                let x = mid(&mut m, nodes[e], nodes[(e + 1) % k]);
                nodes.push(x);
            }
            et = if t == 303 { 306 } else { 408 };
        }
        m.bulk.push(Elem {
            etype: ElementType(et),
            nodes,
            tag: 1,
        });
    }
    for (tag, nodes) in bnd {
        let mut nodes = nodes;
        let mut et = 202;
        if second_order {
            let x = mid(&mut m, nodes[0], nodes[1]);
            nodes.push(x);
            et = 203;
        }
        m.boundary.push(Elem {
            etype: ElementType(et),
            nodes,
            tag,
        });
    }
    m
}
