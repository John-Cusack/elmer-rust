//! The meshes of the refinement: an owned triangle mesh (`Mesh_t` as
//! SplitOneLevel builds it) and the 2-D edge table of `FindMeshEdges2D`
//! (MeshBasics.F90).

use std::collections::HashMap;

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

/// One element of an owned mesh: type, nodes (0-based), body or boundary
/// tag, and the bulk parents of a boundary element.
#[derive(Clone, Debug)]
pub struct Elem {
    pub etype: ElementType,
    pub nodes: Vec<u32>,
    pub tag: i32,
    pub parents: [Option<ElemId>; 2],
}

/// A mesh the refinement made. Bulk elements come first.
#[derive(Clone, Debug)]
pub struct OwnedMesh {
    pub dim: usize,
    pub coords: Vec<[f64; 3]>,
    pub elems: Vec<Elem>,
    pub n_bulk: usize,
}

impl OwnedMesh {
    /// A copy of any mesh.
    pub fn from_view(m: &dyn MeshView) -> OwnedMesh {
        let elems = (0..(m.n_bulk() + m.n_boundary()) as u32)
            .map(|e| {
                let r = m.element(ElemId(e));
                Elem {
                    etype: r.etype,
                    nodes: r.nodes.to_vec(),
                    tag: r.tag,
                    parents: r.parents,
                }
            })
            .collect();
        OwnedMesh {
            dim: m.dim(),
            coords: (0..m.n_nodes() as u32).map(|k| m.coord(k)).collect(),
            elems,
            n_bulk: m.n_bulk(),
        }
    }
}

impl MeshView for OwnedMesh {
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
            parents: el.parents,
        }
    }
}

/// One edge of `FindMeshEdges2D`: end nodes in the order of the first
/// element with it, the middle node of a quadratic edge (`NodeIndexes(3)`,
/// type 203), and its bulk parents (`BoundaryInfo % Left`, `Right`).
#[derive(Clone, Debug)]
pub struct Edge {
    pub nodes: [u32; 2],
    pub mid: Option<u32>,
    pub left: Option<usize>,
    pub right: Option<usize>,
}

impl Edge {
    /// `Edge % TYPE`: 202, or 203 with a middle node.
    pub fn etype(&self) -> ElementType {
        ElementType(if self.mid.is_some() { 203 } else { 202 })
    }

    /// `Edge % NodeIndexes`.
    pub fn element_nodes(&self) -> Vec<u32> {
        let mut v = self.nodes.to_vec();
        v.extend(self.mid);
        v
    }
}

/// `Mesh % Edges` and every bulk element's `EdgeIndexes`.
pub struct Edges {
    pub edges: Vec<Edge>,
    pub of_element: Vec<Vec<usize>>,
}

impl Edges {
    /// `FindMeshEdges2D` over the bulk elements (linear or quadratic
    /// triangles and quadrilaterals) and then the boundary elements, which
    /// only meet edges the bulk already made.
    pub fn find(mesh: &dyn MeshView) -> Edges {
        let nb = mesh.n_bulk();
        let mut edges: Vec<Edge> = Vec::new();
        let mut hash: HashMap<(u32, u32), usize> = HashMap::new();
        let mut of_element = Vec::with_capacity(nb);
        for t in 0..nb + mesh.n_boundary() {
            let el = mesh.element(ElemId(t as u32));
            let n = match el.etype.0 / 100 {
                2 => 1,
                3 => 3,
                4 => 4,
                _ => 0,
            };
            let mut idx = Vec::with_capacity(n);
            for k in 0..n {
                let a = el.nodes[k];
                let b = if n == 1 {
                    el.nodes[1]
                } else {
                    el.nodes[(k + 1) % n]
                };
                let key = (a.min(b), a.max(b));
                // `NodeIndexes(k+n+j-2)` for `j = 2..Degree`: the middle node of
                // a quadratic element's edge.
                let mid = (n > 1 && matches!(el.etype.0, 306 | 408 | 409)).then(|| el.nodes[k + n]);
                let e = *hash.entry(key).or_insert_with(|| {
                    edges.push(Edge {
                        nodes: [a, b],
                        mid,
                        left: None,
                        right: None,
                    });
                    edges.len() - 1
                });
                if t < nb {
                    if edges[e].left.is_some() {
                        edges[e].right = Some(t);
                    } else {
                        edges[e].left = Some(t);
                    }
                }
                idx.push(e);
            }
            if t < nb {
                of_element.push(idx);
            }
        }
        Edges { edges, of_element }
    }
}

/// Fortran `Ew.d` as gfortran writes it: `0.ddd...E+xx`, right-aligned.
pub fn e_format(x: f64, w: usize, d: usize) -> String {
    let body = if x == 0.0 {
        format!("0.{}E+00", "0".repeat(d))
    } else {
        // d significant digits, then shift the point left by one.
        let s = format!("{:.*e}", d - 1, x.abs());
        let (mant, exp) = s.split_once('e').expect("exponent");
        let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
        let e: i32 = exp.parse::<i32>().expect("exponent") + 1;
        let sign = if e < 0 { '-' } else { '+' };
        let body = format!("0.{digits}E{sign}{:02}", e.abs());
        if x < 0.0 {
            format!("-{body}")
        } else {
            body
        }
    };
    format!("{body:>w$}")
}

/// `WriteMeshToDisk` (MeshIO.F90): `mesh.header`, `mesh.nodes`,
/// `mesh.elements`, `mesh.boundary` of `mesh` in `dir`.
pub fn write_mesh_to_disk(
    mesh: &dyn MeshView,
    dir: &std::path::Path,
    bc_tag: &dyn Fn(i32) -> i32,
) -> std::io::Result<()> {
    use std::fmt::Write as _;
    std::fs::create_dir_all(dir)?;
    let nb = mesh.n_bulk();
    let nbd = mesh.n_boundary();
    let el = |i: usize| mesh.element(ElemId(i as u32));
    let biggest = |range: std::ops::Range<usize>| {
        let (mut code, mut max) = (0u16, 0usize);
        for i in range {
            let e = el(i);
            if e.nodes.len() > max {
                code = e.etype.0;
                max = e.nodes.len();
            }
        }
        code
    };
    let mut h = String::new();
    let _ = writeln!(h, "{} {} {}", mesh.n_nodes(), nb, nbd);
    let _ = writeln!(h, "2");
    let _ = writeln!(h, "{} {}", biggest(nb..nb + nbd), nbd);
    let _ = writeln!(h, "{} {}", biggest(0..nb), nb);
    std::fs::write(dir.join("mesh.header"), h)?;
    let mut s = String::new();
    for i in 0..mesh.n_nodes() {
        let c = mesh.coord(i as u32);
        let _ = writeln!(
            s,
            "{} -1 {}{}{} ",
            i + 1,
            e_format(c[0], 23, 15),
            e_format(c[1], 23, 15),
            e_format(c[2], 23, 15)
        );
    }
    std::fs::write(dir.join("mesh.nodes"), s)?;
    let mut s = String::new();
    for i in 0..nb {
        let e = el(i);
        let _ = write!(s, "{} {} {} ", i + 1, e.tag, e.etype.0);
        for n in e.nodes {
            let _ = write!(s, "{} ", n + 1);
        }
        let _ = writeln!(s, " ");
    }
    std::fs::write(dir.join("mesh.elements"), s)?;
    let mut s = String::new();
    for i in 0..nbd {
        let e = el(nb + i);
        let p = |x: Option<ElemId>| x.map_or(0, |p| p.0 as usize + 1);
        let _ = write!(
            s,
            "{} {} {} {} {} ",
            i + 1,
            bc_tag(e.tag),
            p(e.parents[0]),
            p(e.parents[1]),
            e.etype.0
        );
        for n in e.nodes {
            let _ = write!(s, "{} ", n + 1);
        }
        let _ = writeln!(s, " ");
    }
    std::fs::write(dir.join("mesh.boundary"), s)
}
