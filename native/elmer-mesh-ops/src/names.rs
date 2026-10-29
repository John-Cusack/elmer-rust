//! `ReadTargetNames` (MeshIO.F90:959): `Use Mesh Names` gives the sections
//! `Target Bodies` and `Target Boundaries` from the `$ name = tags` lines of
//! `mesh.names`, matching the lowercase name against the sections' `Name`.

use std::path::Path;

use elmer_core::{ParamMap, Params, Value};

use crate::MeshOpsError;

/// Reads `file` (the mesh directory's `mesh.names`) when the Simulation asks
/// for it (`Use Mesh Names`, refined by `Use Mesh Body Names` and `Use Mesh
/// Boundary Names`) and adds the targets to the sections that do not have
/// them. `bodies` and `bcs` are the Body and Boundary Condition sections by
/// position.
pub fn read_target_names(
    file: &Path,
    simulation: &ParamMap,
    bodies: &mut [ParamMap],
    bcs: &mut [ParamMap],
) -> Result<(), MeshOpsError> {
    let all = simulation.logical("use mesh names") == Some(true);
    let do_bodies = simulation.logical("use mesh body names").unwrap_or(all);
    let do_bcs = simulation.logical("use mesh boundary names").unwrap_or(all);
    if !(do_bodies || do_bcs) {
        return Ok(());
    }
    let text = std::fs::read_to_string(file).map_err(|_| {
        MeshOpsError::Invalid(format!(
            "Requested the use of entity names but this file does not exits: {}",
            file.display()
        ))
    })?;
    for line in text.lines() {
        let (Some(i), Some(j)) = (line.find('$'), line.find('=')) else {
            continue;
        };
        let name = line[i + 1..j].trim_start().to_lowercase();
        let name = name.trim_end();
        let tags: Vec<i32> = line[j + 1..]
            .split_whitespace()
            .map_while(|w| w.parse().ok())
            .collect();
        if tags.is_empty() {
            return Err(MeshOpsError::Invalid(format!(
                "Could not find arguments for: {name}"
            )));
        }
        let mut already_set = false;
        for (on, sections, key) in [
            (do_bcs, &mut *bcs, "Target Boundaries"),
            (do_bodies, &mut *bodies, "Target Bodies"),
        ] {
            if !on {
                continue;
            }
            for section in sections.iter_mut() {
                if section.string("name").map(str::to_lowercase).as_deref() != Some(name) {
                    continue;
                }
                if already_set {
                    return Err(MeshOpsError::Invalid(format!(
                        "Mapping of name is not unique: {name}"
                    )));
                }
                if section.get(key).is_none() {
                    section.insert(key, Value::Integer(tags.clone()));
                    already_set = true;
                }
            }
        }
    }
    Ok(())
}
