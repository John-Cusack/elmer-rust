//! `CreateDiscontMesh` (MeshBasics.F90:377): duplicate the nodes on the BCs
//! marked `Discontinuous Boundary = True` so the bodies on the two sides no
//! longer share them.
//!
//! The bulk elements of the `Discontinuous Target Bodies` (`... Body`; else
//! the one body that is a parent of every discontinuous boundary element)
//! move to the new nodes, which are numbered after the old ones in the order
//! of the nodes they copy. A boundary element on the interface keeps the old
//! nodes and puts the moving parent on the right; one touching the interface
//! follows its parents, or the majority of its nodes.
//!
//! With `Discontinuous BC`, `Periodic BC`, `Mortar BC` or `Contact BC` on
//! the BC, the boundary is doubled instead: every interface element gets a
//! twin on the new nodes, tagged as the BC it names, and each of the two is
//! left with the single parent it still touches.
//!
//! Tags are those after `MapBodiesAndBCs`: a bulk tag is the body index and a
//! boundary tag the BC index, as `Element % BodyId` and
//! `BoundaryInfo % Constraint` are in the Fortran. Serial meshes only.

use elmer_core::{ElemId, MeshView, ParamMap, Params, Value};
use elmer_mesh_levels::topology::find_mesh_edges;
use elmer_mesh_levels::Mesh;

use crate::MeshOpsError;

/// A mesh with discontinuous boundaries and the bookkeeping the Fortran
/// keeps on `Mesh` (`DiscontMesh`, `DiscontPerm`, `DiscontNodes`).
#[derive(Clone, Debug)]
pub struct DiscontMesh {
    pub mesh: Mesh,
    /// `Mesh % DiscontPerm`: for each old node, the 1-based number of the
    /// node that copies it (that node is `n_old + perm - 1`), 0 if none.
    pub perm: Vec<u32>,
    /// `Mesh % DiscontNodes`, the number of nodes added; 0 when the boundary
    /// is doubled, where the Fortran keeps no record.
    pub n_discont_nodes: usize,
    /// `Mesh % DiscontMesh`: false when the boundary was doubled.
    pub discont_mesh: bool,
}

fn bc_logical(bcs: &[ParamMap], tag: i32, key: &str) -> bool {
    tag >= 1
        && bcs
            .get(tag as usize - 1)
            .is_some_and(|bc| bc.logical(key) == Some(true))
}

fn int_array(bc: &ParamMap, key: &str) -> Option<Vec<i32>> {
    match bc.get(key) {
        Some(Value::Integer(v)) => Some(v.clone()),
        _ => None,
    }
}

fn parent_body(mesh: &Mesh, p: Option<ElemId>) -> Option<i32> {
    p.map(|p| mesh.elements[p.0 as usize].tag)
}

/// `CreateDiscontMesh(Model, Mesh)` for a serial mesh. `None` when no BC is
/// discontinuous or no node needs duplicating, where the Fortran leaves the
/// mesh as it is.
pub fn create_discont_mesh(
    mesh: &dyn MeshView,
    simulation: &ParamMap,
    bcs: &[ParamMap],
) -> Result<Option<DiscontMesh>, MeshOpsError> {
    let mut double_bc = false;
    let mut active_bcs = 0;
    let mut bc_list = None;
    for bc in bcs {
        if bc.logical("discontinuous boundary") == Some(true) {
            let any = ["discontinuous bc", "periodic bc", "mortar bc", "contact bc"]
                .iter()
                .map(|k| bc.integer(k).unwrap_or(0))
                .sum::<i32>();
            double_bc = any > 0;
            active_bcs += 1;
            bc_list = Some(bc);
        }
    }
    let Some(bc_list) = bc_list.filter(|_| active_bcs > 0) else {
        return Ok(None);
    };

    let mut mesh = Mesh::from_view(mesh);
    let no_nodes = mesh.coords.len();
    let n_bulk = mesh.n_bulk;
    let n_bound = mesh.elements.len() - n_bulk;

    // Mark the nodes that should be discontinuous.
    let mut node = vec![false; no_nodes];
    let mut elem = vec![false; n_bound];
    let mut parent_used = vec![false; n_bulk];
    for t in 0..n_bound {
        let el = &mesh.elements[n_bulk + t];
        if !bc_logical(bcs, el.tag, "discontinuous boundary") {
            continue;
        }
        for &j in &el.nodes {
            node[j as usize] = true;
        }
        elem[t] = true;
        for p in el.parents.into_iter().flatten() {
            parent_used[p.0 as usize] = true;
        }
    }

    non_greedy(&mesh, simulation, &mut node, &elem, &parent_used)?;
    if node.iter().all(|&b| !b) {
        return Ok(None);
    }

    // For historical reasons there is both a single 'body' and multiple
    // 'bodies' that tell which side the new nodes are on.
    let mut target_bodies = int_array(bc_list, "discontinuous target bodies")
        .or_else(|| int_array(bc_list, "discontinuous target body"));
    if target_bodies.is_none() {
        // If either parent is consistently one of the bodies then that body
        // moves.
        let (mut cand_a, mut cand_b) = (-1, -1);
        for t in (0..n_bound).filter(|&t| elem[t]) {
            let [left, right] = mesh.elements[n_bulk + t].parents;
            let (Some(left), Some(right)) = (parent_body(&mesh, left), parent_body(&mesh, right))
            else {
                return Err(MeshOpsError::Invalid(
                    "Alternative strategy requires all parent elements!".into(),
                ));
            };
            if cand_a == -1 {
                cand_a = left;
            } else if cand_a != 0 && cand_a != left && cand_a != right {
                cand_a = 0;
            }
            if cand_b == -1 {
                cand_b = right;
            } else if cand_b != 0 && cand_b != left && cand_b != right {
                cand_b = 0;
            }
        }
        // The bigger one honors the old convention and rules out zero.
        let body = cand_a.max(cand_b);
        if body <= 0 {
            return Err(MeshOpsError::Invalid(
                "No simple rules available for determining discontinuous body".into(),
            ));
        }
        target_bodies = Some(vec![body]);
    }
    let target_bodies = target_bodies.unwrap_or_default();
    let is_target = |body: i32| target_bodies.contains(&body);

    // Nodes related to the moving nodes by the moving elements.
    let mut moving_node = vec![false; no_nodes];
    let mut staying_node = vec![false; no_nodes];
    for el in &mesh.elements[..n_bulk] {
        if !el.nodes.iter().any(|&j| node[j as usize]) {
            continue;
        }
        let marks = if is_target(el.tag) {
            &mut moving_node
        } else {
            &mut staying_node
        };
        for &j in &el.nodes {
            marks[j as usize] = true;
        }
    }

    // A node is doubled only if a moving element has it.
    for (n, m) in node.iter_mut().zip(&moving_node) {
        *n &= *m;
    }
    let mut perm = vec![0u32; no_nodes];
    let mut n_discont = 0;
    for (p, &n) in perm.iter_mut().zip(&node) {
        if n {
            n_discont += 1;
            *p = n_discont;
        }
    }
    let new_index = |perm: &[u32], j: u32| match perm[j as usize] {
        0 => j,
        p => no_nodes as u32 + p - 1,
    };

    // The moving bulk elements take the new nodes.
    for el in &mut mesh.elements[..n_bulk] {
        if !el.nodes.iter().any(|&j| node[j as usize]) || !is_target(el.tag) {
            continue;
        }
        for j in &mut el.nodes {
            *j = new_index(&perm, *j);
        }
    }

    // The boundary elements follow the ownership of their parents, or the
    // majority if that conflicts.
    for t in 0..n_bound {
        let [left, right] = mesh.elements[n_bulk + t].parents;
        if mesh.elements[n_bulk + t].tag == 0 {
            continue;
        }
        if elem[t] {
            let moving = match left {
                Some(l) => is_target(mesh.elements[l.0 as usize].tag),
                None => {
                    let r = right.ok_or_else(|| {
                        MeshOpsError::Invalid("Boundary BC has no parents!".into())
                    })?;
                    !is_target(mesh.elements[r.0 as usize].tag)
                }
            };
            if moving {
                mesh.elements[n_bulk + t].parents = [right, left];
            }
            continue;
        }
        if !mesh.elements[n_bulk + t]
            .nodes
            .iter()
            .any(|&j| node[j as usize])
        {
            continue;
        }
        let mut set = true;
        let moving;
        match (left, right) {
            (Some(l), Some(r)) => {
                let m1 = is_target(mesh.elements[l.0 as usize].tag);
                let m2 = is_target(mesh.elements[r.0 as usize].tag);
                moving = m1;
                if m1 != m2 {
                    // Conflicting moving information.
                    set = false;
                } else if m1 {
                    mesh.elements[n_bulk + t].parents = [right, left];
                }
            }
            (Some(p), None) | (None, Some(p)) => {
                moving = mesh.elements[p.0 as usize]
                    .nodes
                    .iter()
                    .any(|&j| j as usize >= no_nodes);
            }
            (None, None) => {
                return Err(MeshOpsError::Invalid("Boundary BC has no parents!".into()))
            }
        }
        let mut moving = moving;
        if !set {
            let nodes = &mesh.elements[n_bulk + t].nodes;
            let n_moving = nodes.iter().filter(|&&j| moving_node[j as usize]).count();
            let n_staying = nodes.iter().filter(|&&j| staying_node[j as usize]).count();
            if n_staying != n_moving {
                moving = n_moving > n_staying;
                set = true;
            }
        }
        if set && moving {
            for j in &mut mesh.elements[n_bulk + t].nodes {
                *j = new_index(&perm, *j);
            }
        }
    }

    // The new nodes sit where the nodes they copy do.
    mesh.coords.resize(no_nodes + n_discont as usize, [0.0; 3]);
    for (i, &p) in perm.iter().enumerate() {
        if p > 0 {
            mesh.coords[no_nodes + p as usize - 1] = mesh.coords[i];
        }
    }

    if double_bc {
        // The boundary is duplicated and no information of it is saved.
        for t in (0..n_bound).filter(|&t| elem[t]) {
            let el = mesh.elements[n_bulk + t].clone();
            let target = ["discontinuous bc", "mortar bc", "periodic bc", "contact bc"]
                .iter()
                .find_map(|k| bcs.get(el.tag as usize - 1).and_then(|bc| bc.integer(k)));
            let Some(target @ 1..) = target else {
                return Err(MeshOpsError::Invalid(
                    "Nonzero target boundary must be given for all, if any bc!".into(),
                ));
            };
            let mut other = el.clone();
            other.parents = [el.parents[1], None];
            other.nodes = el.nodes.iter().map(|&j| new_index(&perm, j)).collect();
            other.tag = target;
            mesh.elements[n_bulk + t].parents[1] = None;
            mesh.elements.push(other);
        }
        drop_false_parents(&mut mesh);
    }

    Ok(Some(DiscontMesh {
        mesh,
        perm: if double_bc { vec![0; no_nodes] } else { perm },
        n_discont_nodes: if double_bc { 0 } else { n_discont as usize },
        discont_mesh: !double_bc,
    }))
}

/// A parent whose nodes are not all among the element's nodes is no parent
/// any more after the nodes moved; a lone right parent becomes the left one.
fn drop_false_parents(mesh: &mut Mesh) {
    for t in mesh.n_bulk..mesh.elements.len() {
        let mut parents = mesh.elements[t].parents;
        for right in 0..2 {
            let Some(p) = parents[right] else {
                continue;
            };
            let all_in = {
                let parent = &mesh.elements[p.0 as usize].nodes;
                mesh.elements[t].nodes.iter().all(|j| parent.contains(j))
            };
            if all_in {
                if right == 1 && parents[0].is_none() {
                    parents[0] = parents[1];
                    parents[1] = None;
                }
            } else {
                parents[right] = None;
            }
        }
        mesh.elements[t].parents = parents;
    }
}

/// `NonGreedyDiscontinuity`: by default all nodes of the elements at the
/// discontinuous boundary are doubled; the user may release the nodes that
/// other elements (`Discontinuous Bulk Greedy`, `Discontinuous Boundary
/// Greedy`) or the angle sum (`Discontinuous Boundary Full Angle`) need.
fn non_greedy(
    mesh: &Mesh,
    simulation: &ParamMap,
    node: &mut [bool],
    elem: &[bool],
    parent_used: &[bool],
) -> Result<(), MeshOpsError> {
    if node.iter().all(|&b| !b) {
        return Ok(());
    }
    let n_bulk = mesh.n_bulk;
    let greedy_bulk = simulation
        .logical("discontinuous bulk greedy")
        .unwrap_or(true);
    let greedy_bc = simulation
        .logical("discontinuous boundary greedy")
        .unwrap_or(true);
    if !(greedy_bc && greedy_bulk) {
        for (t, el) in mesh.elements.iter().enumerate() {
            if t < n_bulk {
                if greedy_bulk || parent_used[t] {
                    continue;
                }
            } else {
                // An external BC is left alone.
                if greedy_bc || elem[t - n_bulk] || el.parents.iter().any(|p| p.is_none()) {
                    continue;
                }
            }
            if el.nodes.iter().any(|&j| node[j as usize]) {
                for &j in &el.nodes {
                    node[j as usize] = false;
                }
            }
        }
    }

    if simulation.logical("discontinuous boundary full angle") == Some(true) {
        full_angle(mesh, node, elem)?;
    }
    Ok(())
}

/// `Discontinuous Boundary Full Angle`: a node stays single unless the
/// angles of the discontinuous elements around it add up to a full circle
/// (360 degrees inside the body, 90 on its boundary, where the corner nodes
/// are released).
fn full_angle(mesh: &Mesh, node: &mut [bool], elem: &[bool]) -> Result<(), MeshOpsError> {
    if mesh.dim != 3 {
        return Err(MeshOpsError::Unsupported(
            "Discontinuous Boundary Full Angle in a mesh that is not 3D".into(),
        ));
    }
    let n_bulk = mesh.n_bulk;
    let mut node_phi = vec![0.0f64; mesh.coords.len()];
    for t in (0..mesh.elements.len() - n_bulk).filter(|&t| elem[t]) {
        let el = &mesh.elements[n_bulk + t];
        let n = (el.etype.0 / 100) as usize;
        // Coords(1:n,1:3) = (y, x, z)
        let c: Vec<[f64; 3]> = el.nodes[..n]
            .iter()
            .map(|&j| {
                let x = mesh.coords[j as usize];
                [x[1], x[0], x[2]]
            })
            .collect();
        for i in 0..n {
            let i1 = (i + 1) % n;
            let i2 = (n + i + n - 1) % n;
            let unit = |k: usize| {
                let e = [c[k][0] - c[i][0], c[k][1] - c[i][1], c[k][2] - c[i][2]];
                let norm = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
                [e[0] / norm, e[1] / norm, e[2] / norm]
            };
            let (e1, e2) = (unit(i1), unit(i2));
            let phi = (e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2]).acos();
            node_phi[el.nodes[i] as usize] += phi;
        }
    }
    for p in &mut node_phi {
        *p = 180.0 * *p / std::f64::consts::PI;
    }

    // A node is on the boundary if it is on a face with only one bulk element.
    let topo = find_mesh_edges(mesh)?;
    let mut face_elems = vec![0u32; topo.faces.len()];
    for faces in &topo.element_faces[..n_bulk] {
        for &f in faces {
            face_elems[f as usize] += 1;
        }
    }
    let mut boundary_node = vec![false; mesh.coords.len()];
    for (f, _) in face_elems.iter().enumerate().filter(|(_, &c)| c == 1) {
        for &j in &topo.faces[f] {
            boundary_node[j as usize] = true;
        }
    }

    for i in 0..node.len() {
        if !node[i] {
            continue;
        }
        // On the boundary the discontinuity is released when the angle is
        // about 90 degrees, i.e. on corner nodes; elsewhere when it is below
        // 360 degrees.
        let limit = if boundary_node[i] { 100.0 } else { 350.0 };
        if node_phi[i] < limit {
            node[i] = false;
        }
    }
    Ok(())
}
