//! `IncreaseElementOrder` (MeshBasics.F90:3564): every linear element gets
//! one node at the middle of each of its edges and becomes the serendipity
//! element with that many nodes (202 to 203, 303 to 306, 404 to 408, 504 to
//! 510, 706 to 715, 808 to 820).
//!
//! The new nodes follow the old ones in the order of the global edges of
//! `FindMeshEdges`, and every element lists its mid-edge nodes in its local
//! edge order after its corner nodes. A brick is the exception: its 820
//! node list has the four mid-edge nodes of the top face before the four
//! vertical ones, while its edge map has them the other way round, so those
//! two groups of four swap places.

use elmer_core::{ElementType, MeshView};
use elmer_mesh_levels::topology::find_mesh_edges;
use elmer_mesh_levels::Mesh;

use crate::MeshOpsError;

/// The quadratic mesh of the linear `mesh`, as `PrepareMesh` makes it for
/// `Increase Element Order = True`, without the `FollowCurvedBoundary`
/// passes it brackets the call with (callers refuse `Follow * Boundary`).
/// Element order, tags and boundary parents are unchanged.
pub fn increase_element_order(mesh: &dyn MeshView) -> Result<Mesh, MeshOpsError> {
    let mut mesh = Mesh::from_view(mesh);
    if let Some(e) = mesh.elements.iter().find(|e| e.etype.family() == 6) {
        return Err(MeshOpsError::Unsupported(format!(
            "pyramid element type {} in Increase Element Order",
            e.etype.0
        )));
    }
    if let Some(e) = mesh
        .elements
        .iter()
        .find(|e| e.etype.n_nodes() != nodes_of_linear(e.etype))
    {
        return Err(MeshOpsError::Invalid(format!(
            "Cannot increase element order for: {}",
            e.etype.0
        )));
    }
    let topo = find_mesh_edges(&mesh)?;
    let n0 = mesh.coords.len() as u32;
    for &[a, b] in &topo.edges {
        let (a, b) = (mesh.coords[a as usize], mesh.coords[b as usize]);
        mesh.coords.push([
            0.5 * (a[0] + b[0]),
            0.5 * (a[1] + b[1]),
            0.5 * (a[2] + b[2]),
        ]);
    }
    for (el, edges) in mesh.elements.iter_mut().zip(&topo.element_edges) {
        if el.etype.0 == 101 {
            continue;
        }
        el.etype = ElementType(el.etype.0 + edges.len() as u16);
        el.nodes.extend(edges.iter().map(|&e| n0 + e));
        if el.etype.0 == 820 {
            // Tinds = N(17:20); N(17:20) = N(13:16); N(13:16) = Tinds.
            let (head, tail) = el.nodes.split_at_mut(16);
            head[12..16].swap_with_slice(tail);
        }
    }
    Ok(mesh)
}

/// The node count of the linear element of the family of `t`; a type with
/// another count is already of higher order.
fn nodes_of_linear(t: ElementType) -> usize {
    match t.family() {
        1 => 1,
        2 => 2,
        3 => 3,
        4 => 4,
        5 => 4,
        7 => 6,
        8 => 8,
        _ => 0,
    }
}
