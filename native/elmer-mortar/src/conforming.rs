//! The conforming part of `GeneratePeriodicProjectors` (MortarUtils.F90):
//! for every `Conforming BC = k` boundary condition, `PeriodicPermutation`
//! pairs each of its nodes with the nearest node of boundary condition `k`
//! (`ConformingNodePerm`), and chains through a second map are resolved
//! once. The result is `Mesh % PeriodicPerm` and `Mesh % PeriodicFlip` as an
//! [`elmer_solver::conforming::ConformingMap`], which `elmer-solver` applies
//! to the solvers with `Apply Conforming BCs`.
//!
//! Nodes only: the edge and face permutations of edge and p elements
//! (`ConformingEdgePerm`, `ConformingFacePerm`) are not ported.

use elmer_core::{Params, Value, NO_DOF};
use elmer_solver::conforming::ConformingMap;

use crate::interface::{create_interface_meshes, map_interface_coordinate, InterfaceMesh};
use crate::overlay::overlay_interface_meshes;
use crate::projector::ProjectorEnv;
use crate::MortarError;

const CALLER: &str = "PeriodicPermutation";

/// Projector families `PeriodicPermutation` handles through interface
/// meshes that are not ported.
const REFUSED: [&str; 11] = [
    "rotational projector",
    "anti rotational projector",
    "cylindrical projector",
    "radial projector",
    "anti radial projector",
    "axial projector",
    "anti axial projector",
    "sliding projector",
    "anti sliding projector",
    "flat projector",
    "level projector",
];

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// `ListGetInteger(BC, 'Conforming BC')`.
fn conforming_target(bc: &dyn Params) -> Option<usize> {
    match bc.get("conforming bc")? {
        Value::Integer(v) => v.first().map(|&k| k.max(0) as usize),
        Value::Real { data, .. } => data.first().map(|&k| k.round().max(0.0) as usize),
        _ => None,
    }
}

/// `Mesh % PeriodicPerm` and `Mesh % PeriodicFlip` of the model, or `None`
/// when no boundary condition has `Conforming BC`.
pub fn generate_conforming_map(env: &ProjectorEnv<'_>) -> Result<Option<ConformingMap>, MortarError> {
    if !env.bcs.iter().any(|b| b.get("conforming bc").is_some()) {
        return Ok(None);
    }
    let n = env.mesh.n_nodes();
    // 1-based partner, 0 where none, as the Fortran tables.
    let mut perm = vec![0usize; n];
    let mut flip = vec![false; n];
    for this in 1..=env.bcs.len() {
        if let Some(target) = conforming_target(env.bcs[this - 1]) {
            periodic_permutation(env, this, target, &mut perm, &mut flip)?;
        }
    }
    let (mut noconf, mut nocyclic) = (0usize, 0usize);
    for i in 0..n {
        let j = perm[i];
        if j > 0 {
            noconf += 1;
            if perm[j - 1] > 0 {
                perm[i] = perm[j - 1];
                flip[i] = if flip[i] { !flip[j - 1] } else { flip[j - 1] };
                nocyclic += 1;
            }
        }
    }
    let noflip = flip.iter().filter(|&&f| f).count();
    let gen = "GeneratePeriodicProjectors";
    env.log.info(gen, &format!("Number of conforming maps: {noconf}"), 8);
    if nocyclic > 0 {
        env.log.info(gen, &format!("Number of cyclic maps: {nocyclic}"), 8);
    }
    if noflip > 0 {
        env.log.info(gen, &format!("Number of periodic flips: {noflip}"), 8);
    }
    Ok(Some(ConformingMap {
        partner: perm
            .iter()
            .map(|&p| if p > 0 { (p - 1) as u32 } else { NO_DOF })
            .collect(),
        flip,
    }))
}

/// `PeriodicPermutation(Model, Mesh, This, Trgt, PerPerm, PerFlip)` for the
/// nodes.
fn periodic_permutation(
    env: &ProjectorEnv<'_>,
    this: usize,
    target: usize,
    perm: &mut [usize],
    flip: &mut [bool],
) -> Result<(), MortarError> {
    if this == 0 || target == 0 {
        return Ok(());
    }
    env.log.info(
        CALLER,
        &format!("Creating mapping between BCs {this} and {target}"),
        8,
    );
    let bc = env.bcs[this - 1];
    let split = env
        .simulation
        .logical("mortar bcs split quadratic")
        .unwrap_or(false);
    let Some((mut b1, mut b2)) = create_interface_meshes(env.mesh, env.bc_of, this, target, split)?
    else {
        return Ok(());
    };
    let map: Option<Vec<i32>> = match bc.get("projector coordinate mapping") {
        Some(Value::Integer(v)) => Some(v.clone()),
        _ => None,
    };
    map_interface_coordinate(&mut b1, &mut b2, map.as_deref())?;

    if let Some(k) = REFUSED.iter().find(|k| flag(bc, k)) {
        return Err(MortarError::Unsupported(format!(
            "Conforming BC with {k} (only the plane and the overlaid generic interface are ported)"
        )));
    }
    let anti_plane = flag(bc, "anti plane projector");
    let plane = anti_plane || flag(bc, "plane projector");
    let anti = bc.logical("antisymmetric bc").unwrap_or(anti_plane);

    // `PlaneInterfaceMeshes` would turn both interface meshes into the same
    // orthonormal frame (tangents and plane normal). The nearest-node match
    // below is invariant under that, so it is not applied. Otherwise the
    // overlay defaults to on (not a sliding projector).
    let enforce = bc.logical("mortar bc enforce overlay");
    if enforce.unwrap_or(!plane) {
        overlay_interface_meshes(&b1, &mut b2, bc, env.log)?;
    }
    if !flag(bc, "projector skip nodes") {
        conforming_node_perm(env, &b1, &b2, perm, flip, anti)?;
    }
    env.log.info(CALLER, "Periodic permutation created, now exiting...", 8);
    Ok(())
}

/// `ConformingNodePerm`: each node of `b1` not yet mapped takes the nearest
/// node of `b2` (the first at the smallest squared distance, stopping at one
/// closer than machine epsilon). A node on both sides keeps its own dof.
fn conforming_node_perm(
    env: &ProjectorEnv<'_>,
    b1: &InterfaceMesh,
    b2: &InterfaceMesh,
    perm: &mut [usize],
    flip: &mut [bool],
    anti: bool,
) -> Result<(), MortarError> {
    let caller = "ConformingNodePerm";
    env.log.info(caller, "Creating permutations for conforming nodes", 8);
    if perm.is_empty() || b1.coords.is_empty() || b2.coords.is_empty() {
        return Ok(());
    }
    let mut used = vec![false; b2.coords.len()];
    let (mut same, mut double, mut hits) = (0usize, 0usize, 0usize);
    let mut max_min = 0.0f64;
    for (i1, x1) in b1.coords.iter().enumerate() {
        let j1 = b1.inv_perm[i1] as usize;
        if perm[j1] > 0 {
            continue;
        }
        let mut min = f64::MAX;
        let mut mini = 0;
        for (i2, x2) in b2.coords.iter().enumerate() {
            let (dx, dy, dz) = (x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]);
            let ss = dx * dx + dy * dy + dz * dz;
            if ss < min {
                min = ss;
                mini = i2;
            }
            if min < f64::EPSILON {
                break;
            }
        }
        let j2 = b2.inv_perm[mini] as usize;
        if j1 == j2 {
            same += 1;
            continue;
        }
        if used[mini] {
            double += 1;
        } else {
            used[mini] = true;
            hits += 1;
        }
        perm[j1] = j2 + 1;
        max_min = max_min.max(min);
        if anti {
            flip[j1] = true;
        }
    }
    if same > 0 {
        env.log.info(caller, &format!("Number of nodes are the same: {same}"), 8);
    }
    env.log.info(caller, &format!("Number of conforming nodes found: {hits}"), 8);
    env.log.info(
        caller,
        &format!("Maximum minimum deviation in node coords: {:12.4E}", max_min.sqrt()),
        10,
    );
    if double > 0 {
        return Err(MortarError::Fatal(format!(
            "{caller}: This is not conforming! Number of nodes used twice: {double} (out of {})",
            b1.coords.len()
        )));
    }
    Ok(())
}
