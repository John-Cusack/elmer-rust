//! Transform proof driver: parse a `.sif`, load its mesh, apply the
//! S1 load-time transforms, and write the transformed mesh DB back out,
//! plus a B-side `.sif` whose triggers are already applied.
//!
//! Usage: `meshtransform -- <sif> <meshdir> <outmeshdir> <outb.sif>`
//!
//! The B sif is the merged model with model-dependent triggers removed
//! (Target Bodies/Boundaries/Nodes/Coordinates, Default Internal/
//! External BC, Intersection BC/Body, Increase Element Order, Coordinate
//! Mapping/Scaling) and explicit `Boundary N` sections pinning each
//! nonzero transformed constraint to its BC. Running the pristine solver
//! on `<meshdir>` + `<sif>` vs `<outmeshdir>` + `<outb.sif>` and comparing
//! norms proves the native transforms match Elmer's.

use std::collections::BTreeMap;
use std::io::Write;

fn strip_triggers(model: &mut elmer_sif::Model) {
    // Drop only non-idempotent triggers and intersection markers.
    // Target lists below are rewritten as the inverse image of the
    // transformed mesh, so re-running Elmer's mapping is a fixed point.
    for s in model
        .bodies
        .iter_mut()
        .chain(model.bcs.iter_mut())
        .chain(model.ics.iter_mut())
        .chain(model.materials.iter_mut())
        .chain(model.body_forces.iter_mut())
        .chain(model.equations.iter_mut())
        .chain(model.solvers.iter_mut())
        .chain(model.components.iter_mut())
    {
        s.entries.retain(|e| {
            !matches!(
                e.name.as_str(),
                "intersection bc" | "intersection body"
            )
        });
    }
    if let Some(sim) = model.simulation.as_mut() {
        sim.entries.retain(|e| {
            !matches!(
                e.name.as_str(),
                "coordinate mapping"
                    | "coordinate scaling"
                    | "increase element order"
                    | "split mesh quads"
                    | "split mesh prisms"
            )
        });
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (sif_path, mesh_dir, out_dir) = (&args[1], &args[2], &args[3]);
    let out_b: Option<&String> = args.get(4);

    let parsed = elmer_sif::parse_file(sif_path).unwrap_or_else(|e| {
        eprintln!("sif parse failed: {e}");
        std::process::exit(2);
    });
    let mut sif = parsed.sif;
    elmer_sif::complete(&mut sif, sif_path).unwrap_or_else(|e| {
        eprintln!("complete failed: {e}");
        std::process::exit(2);
    });
    let model = sif.merged();

    let mut mesh = elmer_mesh::load(mesh_dir).unwrap_or_else(|e| {
        eprintln!("mesh load failed: {e}");
        std::process::exit(2);
    });
    elmer_mesh::transform(&mut mesh, &model).unwrap_or_else(|e| {
        eprintln!("transform failed: {e}");
        std::process::exit(2);
    });
    elmer_mesh::mesh::write_db(&mesh, out_dir).unwrap_or_else(|e| {
        eprintln!("mesh write failed: {e}");
        std::process::exit(2);
    });

    // B sif: identity Targets so Elmer's re-mapping is a fixed point.
    // (Transformed constraints/bodies are already BC/body indices;
    // Target Bodies/Boundaries of [id] map id->id. Elmer's Boundary
    // sections do NOT substitute for Targets, so none are emitted.)
    let mut bmodel = model.clone();
    strip_triggers(&mut bmodel);
    for (i, b) in bmodel.bodies.iter_mut().enumerate() {
        let id = i as i32 + 1;
        b.entries.retain(|e| e.name != "target bodies");
        b.entries.push(elmer_sif::Entry {
            name: "target bodies".to_string(),
            value: elmer_sif::Value::IntArray { values: vec![id] },
            partag: None,
            distribute: false,
        });
    }
    bmodel.boundaries.clear();
    for (i, bc) in bmodel.bcs.iter_mut().enumerate() {
        let id = i as i32 + 1;
        bc.entries.retain(|e| e.name != "target boundaries");
        bc.entries.push(elmer_sif::Entry {
            name: "target boundaries".to_string(),
            value: elmer_sif::Value::IntArray { values: vec![id] },
            partag: None,
            distribute: false,
        });
    }
    let mut bsif = elmer_sif::Sif { runs: vec![bmodel] };
    let _ = elmer_sif::complete(&mut bsif, sif_path);
    let out_mesh_name = std::path::Path::new(out_dir.as_str())
        .file_name()
        .unwrap()
        .to_string_lossy()
        .into_owned();
    let text = elmer_sif::write_sif(&bsif, Some((".", out_mesh_name.as_str())));
    if let Some(out_b) = out_b {
        std::fs::write(out_b, &text).unwrap_or_else(|e| {
            eprintln!("b sif write failed: {e}");
            std::process::exit(2);
        });
    }

    // Summary for the proof log.
    let mut f = std::fs::File::create(format!("{out_dir}.summary")).unwrap();
    writeln!(
        f,
        "nodes={} bulk={} boundary={} dim={}",
        mesh.nodes.len(),
        mesh.bulk.len(),
        mesh.boundary.len(),
        mesh.dim
    )
    .unwrap();
    let mut codes: BTreeMap<u16, usize> = Default::default();
    for e in &mesh.bulk {
        *codes.entry(e.code).or_insert(0) += 1;
    }
    for b in &mesh.boundary {
        *codes.entry(b.code).or_insert(0) += 1;
    }
    writeln!(f, "codes={codes:?}").unwrap();
    println!(
        "nodes={} bulk={} boundary={} dim={} codes={codes:?}",
        mesh.nodes.len(),
        mesh.bulk.len(),
        mesh.boundary.len(),
        mesh.dim
    );
}
