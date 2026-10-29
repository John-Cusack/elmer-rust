//! Output for the Fortran field values must equal what the pristine
//! `ElmerSolver` wrote. The `.vtu`/`.pvd`/`.pvtu` fixtures in `tests/fixtures/`
//! are the output of `ElmerSolver` (Post File = case.vtu) on the meshes in
//! `tests/fixtures/meshes/`; the Rust writer is fed the mesh files and the field
//! arrays parsed back out of the fixture. ASCII and binary files are compared
//! byte for byte (the arrays decode equal as a consequence).

mod support;

use elmer_core::{ParamMap, Value};
use elmer_vtu::{
    pvd_time, render_pvd, render_pvtu, render_vtu, GeometryIds, PvdEntry, VtuOptions,
};
use support::*;

fn opts(ascii: bool, scalars: &[&str], vectors: &[&str], bc_tags: &[i32]) -> VtuOptions {
    let mut p = ParamMap::new();
    p.insert("Ascii Output", Value::Logical(ascii));
    p.insert("Skip Metainfo", Value::Logical(true));
    p.insert("Save Geometry Ids", Value::Logical(true));
    for (i, s) in scalars.iter().enumerate() {
        p.insert(&format!("Scalar Field {}", i + 1), Value::Str((*s).into()));
    }
    for (i, s) in vectors.iter().enumerate() {
        p.insert(&format!("Vector Field {}", i + 1), Value::Str((*s).into()));
    }
    VtuOptions::from_params(&p).with_geometry(GeometryIds {
        n_bodies: 1,
        bc_tags: bc_tags.to_vec(),
        ..Default::default()
    })
}

fn render(mesh_dir: &str, fixture_name: &str, o: &VtuOptions) -> (Vec<u8>, Vec<u8>) {
    let mesh = TestMesh::read(mesh_dir);
    let want = fixture_bytes(fixture_name);
    let fields = TestFields::from_vtu(&parse_vtu(&want));
    let got = render_vtu(&mesh, &fields, o, None).unwrap().unwrap();
    (got, want)
}

#[test]
fn triangle_ascii_and_binary() {
    for (name, ascii) in [("mini_asc.vtu", true), ("mini_bin.vtu", false)] {
        let o = opts(ascii, &["temperature"], &[], &[2, 3]);
        let (got, want) = render("triangle", name, &o);
        assert_bytes_eq(&got, &want, name);
        assert_eq!(parse_vtu(&got), parse_vtu(&want));
    }
}

#[test]
fn mixed_linear_elements_ascii_and_binary() {
    // 808, 706, 605, 504 bulk; 404, 303 boundary; scalar and 3-vector.
    for (name, ascii) in [("unit_asc.vtu", true), ("unit_bin.vtu", false)] {
        let o = opts(ascii, &["temperature"], &["displacement"], &[1, 2]);
        let (got, want) = render("unitelems", name, &o);
        assert_bytes_eq(&got, &want, name);
    }
}

#[test]
fn single_precision_binary_from_double_fields() {
    // unit_bin holds the double-precision solution; unit_bin32 was written from
    // the same solution with `Single Precision = True`.
    let mesh = TestMesh::read("unitelems");
    let fields = TestFields::from_vtu(&parse_vtu(&fixture_bytes("unit_bin.vtu")));
    let mut o = opts(false, &["temperature"], &["displacement"], &[1, 2]);
    o.single_precision = true;
    let got = render_vtu(&mesh, &fields, &o, None).unwrap().unwrap();
    assert_bytes_eq(&got, &fixture_bytes("unit_bin32.vtu"), "unit_bin32");
}

#[test]
fn boundaries_only_compacts_nodes_and_offsets_ids() {
    // Same fields as the full file; only the 404 and 303 boundary elements are
    // saved, so nodes 1-4 are renumbered 0-3 and BC ids start from 1, not 101.
    let mesh = TestMesh::read("unitelems");
    let full = parse_vtu(&fixture_bytes("unit_asc.vtu"));
    let fields = TestFields::from_vtu(&full);
    let mut o = opts(true, &["temperature"], &["displacement"], &[1, 2]);
    o.save_boundaries_only = true;
    let got = render_vtu(&mesh, &fields, &o, None).unwrap().unwrap();
    // The fixture's four saved nodes carry the first four values of each field.
    assert_bytes_eq(&got, &fixture_bytes("unit_bnd.vtu"), "unit_bnd");
}

#[test]
fn quadratic_elements_reordered() {
    // 510 tetra and 820 hexahedron, 306 and 408 boundary. The 820 connectivity
    // is reordered (Order820) and 306/408 map to VTK 22/23.
    for (name, ascii) in [("quad_asc.vtu", true), ("quad_bin.vtu", false)] {
        let o = opts(ascii, &["temperature"], &["displacement"], &[1, 2]);
        let (got, want) = render("quadmesh", name, &o);
        assert_bytes_eq(&got, &want, name);
    }
}

#[test]
fn save_linear_elements() {
    // Only the 4 + 8 corner nodes are saved and renumbered; the fields are the
    // full 30-node solution of the quadratic run.
    let mesh = TestMesh::read("quadmesh");
    let fields = TestFields::from_vtu(&parse_vtu(&fixture_bytes("quad_asc.vtu")));
    let mut o = opts(true, &["temperature"], &["displacement"], &[1, 2]);
    o.save_linear = true;
    let got = render_vtu(&mesh, &fields, &o, None).unwrap().unwrap();
    assert_bytes_eq(&got, &fixture_bytes("quad_lin.vtu"), "quad_lin");
}

#[test]
fn transient_collection_files() {
    let mut o = opts(true, &["temperature"], &[], &[2, 3]);
    o.time_collection = true;
    let mesh = TestMesh::read("triangle");

    let mut entries = Vec::new();
    for (n, t) in [(1, 0.25), (2, 0.5), (3, 0.75)] {
        let vtu = fixture_bytes(&format!("mini_tr_t000{n}.vtu"));
        let fields = TestFields::from_vtu(&parse_vtu(&vtu));
        let got = render_vtu(&mesh, &fields, &o, None).unwrap().unwrap();
        assert_bytes_eq(&got, &vtu, &format!("mini_tr_t000{n}.vtu"));

        // `Enforce Parallel format` was on in this run: a .pvtu wrapper.
        let pvtu = render_pvtu(&fields, &o, &[format!("case_t000{n}.vtu")]).unwrap();
        assert_bytes_eq(pvtu.as_bytes(), &fixture_bytes(&format!("mini_tr_t000{n}.pvtu")), "pvtu");

        entries.push(PvdEntry {
            timestep: pvd_time(t, 0.25, &o),
            group: 0,
            file: format!("case_t000{n}.pvtu"),
        });
    }
    let pvd = render_pvd(&entries).unwrap();
    assert_bytes_eq(pvd.as_bytes(), &fixture_bytes("mini_tr.pvd"), "pvd");
}

#[test]
fn g0_time_strings_match_fortran() {
    let table = String::from_utf8(fixture_bytes("g0_times.tsv")).unwrap();
    for line in table.lines() {
        let (dt, want) = line.split_once('\t').unwrap();
        let entry = PvdEntry { timestep: dt.parse().unwrap(), group: 0, file: "f.vtu".into() };
        let pvd = render_pvd(&[entry]).unwrap();
        assert!(
            pvd.contains(&format!("timestep=\"{want}\"")),
            "dt {dt}: wanted {want} in {pvd}"
        );
    }
}
