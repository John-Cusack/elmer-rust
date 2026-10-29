//! The two boundary meshes of a projector (`CreateInterfaceMeshes` and
//! `MapInterfaceCoordinate`, MortarUtils.F90).

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

use crate::MortarError;

#[derive(Clone, Debug)]
struct IfaceElem {
    etype: ElementType,
    /// Local node indices.
    nodes: Vec<u32>,
    /// The boundary condition the parent boundary element belongs to.
    tag: i32,
}

/// One side of an interface: the boundary elements of one boundary
/// condition, renumbered to their own nodes. Its elements are the "bulk"
/// elements of the view, as in Elmer's `BMesh`. `dim()` is 0 like an
/// `AllocateMesh` result, so point location uses the coordinate system
/// dimension.
#[derive(Clone, Debug)]
pub struct InterfaceMesh {
    /// Local node coordinates; the overlay transforms these in place.
    pub coords: Vec<[f64; 3]>,
    /// Local node -> parent mesh node (0-based), ascending (`InvPerm`).
    pub inv_perm: Vec<u32>,
    elems: Vec<IfaceElem>,
}

impl InterfaceMesh {
    pub fn n_elements(&self) -> usize {
        self.elems.len()
    }
}

impl MeshView for InterfaceMesh {
    fn dim(&self) -> usize {
        0
    }
    fn n_nodes(&self) -> usize {
        self.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.elems.len()
    }
    fn n_boundary(&self) -> usize {
        0
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

/// `Mortar BCs Split Quadratic`: corner/midside sub-elements of a quadratic
/// face, in Elmer's local node numbering (1-based as in the Fortran).
fn split_map(code: u16) -> Option<&'static [&'static [usize]]> {
    match code {
        408 => Some(&[
            &[5, 6, 7, 8],
            &[1, 5, 8],
            &[2, 6, 5],
            &[3, 7, 6],
            &[4, 8, 7],
        ]),
        409 => Some(&[&[1, 5, 9, 8], &[2, 6, 9, 5], &[3, 7, 9, 6], &[4, 8, 9, 7]]),
        306 => Some(&[&[1, 4, 6], &[2, 5, 4], &[3, 6, 5], &[4, 5, 6]]),
        _ => None,
    }
}

/// `CreateInterfaceMeshes` (serial, nodal): boundary elements whose boundary
/// condition is `this` go to the first mesh, those of `target` to the
/// second. `bc_of` maps a boundary element's tag to its 1-based boundary
/// condition (`BoundaryInfo % Constraint`). `split_quadratic` is
/// `Mortar BCs Split Quadratic`, already cleared by the caller when the mesh
/// has edges or faces. Returns `None` when either side is empty (`Success =
/// .FALSE.`).
pub fn create_interface_meshes(
    mesh: &dyn MeshView,
    bc_of: &dyn Fn(i32) -> Option<usize>,
    this: usize,
    target: usize,
    split_quadratic: bool,
) -> Result<Option<(InterfaceMesh, InterfaceMesh)>, MortarError> {
    if this == 0 || target == 0 {
        return Err(MortarError::Fatal(
            "CreateInterfaceMeshes: Invalid target boundaries".into(),
        ));
    }
    let nb = mesh.n_bulk();
    let n_nodes = mesh.n_nodes();
    let mut sides: [Vec<IfaceElem>; 2] = [Vec::new(), Vec::new()];
    let mut perms = [vec![0u32; n_nodes], vec![0u32; n_nodes]];
    for i in 0..mesh.n_boundary() {
        let el = mesh.element(ElemId((nb + i) as u32));
        let constraint = bc_of(el.tag).unwrap_or(0);
        let side = if constraint == this {
            0
        } else if constraint == target {
            1
        } else {
            continue;
        };
        let split = if split_quadratic {
            split_map(el.etype.0)
        } else {
            None
        };
        match split {
            Some(parts) => {
                for part in parts {
                    sides[side].push(IfaceElem {
                        etype: ElementType(101 * part.len() as u16),
                        nodes: part.iter().map(|&k| el.nodes[k - 1]).collect(),
                        tag: constraint as i32,
                    });
                }
            }
            None => sides[side].push(IfaceElem {
                etype: el.etype,
                nodes: el.nodes.to_vec(),
                tag: constraint as i32,
            }),
        }
        for &n in el.nodes {
            perms[side][n as usize] = 1;
        }
    }
    if sides[0].is_empty() || sides[1].is_empty() {
        return Ok(None);
    }
    let [elems1, elems2] = sides;
    let [perm1, perm2] = perms;
    let build = |mut perm: Vec<u32>, mut elems: Vec<IfaceElem>| {
        let mut inv_perm = Vec::new();
        let mut coords = Vec::new();
        for (i, p) in perm.iter_mut().enumerate() {
            if *p > 0 {
                inv_perm.push(i as u32);
                coords.push(mesh.coord(i as u32));
                *p = inv_perm.len() as u32;
            }
        }
        for e in &mut elems {
            for n in &mut e.nodes {
                *n = perm[*n as usize] - 1;
            }
        }
        InterfaceMesh {
            coords,
            inv_perm,
            elems,
        }
    };
    Ok(Some((build(perm1, elems1), build(perm2, elems2))))
}

/// `MapInterfaceCoordinate`: with `Projector Coordinate Mapping = a b c`
/// (a permutation of 1 2 3) the interface x, y, z become the parent's
/// coordinates a, b, c on both meshes.
pub fn map_interface_coordinate(
    b1: &mut InterfaceMesh,
    b2: &mut InterfaceMesh,
    map: Option<&[i32]>,
) -> Result<(), MortarError> {
    let Some(map) = map else {
        return Ok(());
    };
    if map.len() != 3 || (1..=3).any(|c| !map.contains(&c)) {
        return Err(MortarError::Fatal(format!(
            "MapInterfaceCoordinate: Inconsistent Coordinate Mapping: {map:?}"
        )));
    }
    let pick = |c: i32| (c - 1) as usize;
    for b in [b1, b2] {
        for x in &mut b.coords {
            *x = [x[pick(map[0])], x[pick(map[1])], x[pick(map[2])]];
        }
    }
    Ok(())
}
