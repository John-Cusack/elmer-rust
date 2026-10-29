//! Interpolation between small meshes: exact for linear fields on linear
//! elements and quadratic fields on quadratic ones.

mod support;

use elmer_core::{ElemId, MeshView, NullLogger, NO_DOF};
use elmer_interp::{
    build_projector, interpolate, DgIndexes, Extras, FieldData, InterpError, InterpOptions,
    Location, Locator,
};
use support::*;

fn lin2(p: [f64; 3]) -> f64 {
    1.0 + 2.0 * p[0] - 3.0 * p[1]
}

fn lin3(p: [f64; 3]) -> f64 {
    0.5 + 2.0 * p[0] - 3.0 * p[1] + 1.5 * p[2]
}

fn field_of(mesh: &TestMesh, f: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
    mesh.coords.iter().map(|&c| f(c)).collect()
}

fn run(
    old: &TestMesh,
    new: &TestMesh,
    opts: InterpOptions,
    old_f: Vec<f64>,
    sentinel: f64,
) -> (Vec<f64>, elmer_interp::InterpReport) {
    let loc = Locator::new(old, opts).unwrap();
    let of = vec![FieldData::nodal("temp", old_f)];
    let mut nf = vec![FieldData::nodal("temp", vec![sentinel; new.n_nodes()])];
    let rep = interpolate(&loc, new, &of, &mut nf, &Extras::default(), &NullLogger).unwrap();
    (nf.remove(0).values, rep)
}

#[test]
fn tri_to_quad_linear_exact_and_outside_untouched() {
    let old = tri_grid(0.0, 1.0, 0.0, 1.0, 5, 4);
    // New mesh pokes out of the old one on the right and top.
    let new = quad_grid(0.1, 1.3, 0.05, 1.1, 7, 6);
    let (vals, rep) = run(
        &old,
        &new,
        InterpOptions::default(),
        field_of(&old, lin2),
        -99.0,
    );
    let mut inside = 0;
    for (i, c) in new.coords.iter().enumerate() {
        let want_in = c[0] <= 1.0 + 1e-12 && c[1] <= 1.0 + 1e-12;
        assert_eq!(rep.found[i], want_in, "node {i} at {c:?}");
        if want_in {
            inside += 1;
            assert_close(vals[i], lin2(*c), 1e-12, "value");
        } else {
            assert_eq!(vals[i], -99.0);
        }
    }
    assert_eq!(rep.n_not_found, new.n_nodes() - inside);
}

#[test]
fn quad_to_tri_linear_exact() {
    let old = quad_grid(-1.0, 2.0, -1.0, 2.0, 6, 6);
    let new = tri_grid(0.0, 1.0, 0.0, 1.0, 9, 9);
    let (vals, rep) = run(
        &old,
        &new,
        InterpOptions::default(),
        field_of(&old, lin2),
        0.0,
    );
    assert_eq!(rep.n_not_found, 0);
    for (v, c) in vals.iter().zip(&new.coords) {
        assert_close(*v, lin2(*c), 1e-12, "value");
    }
}

#[test]
fn tet_to_hex_and_hex_to_tet_linear_exact() {
    let tets = tet_grid([0.0; 3], [1.0; 3], [3, 3, 3]);
    let hexes = hex_grid([0.1; 3], [0.9; 3], [4, 3, 5]);
    let (v, rep) = run(
        &tets,
        &hexes,
        InterpOptions::default(),
        field_of(&tets, lin3),
        0.0,
    );
    assert_eq!(rep.n_not_found, 0);
    for (v, c) in v.iter().zip(&hexes.coords) {
        assert_close(*v, lin3(*c), 1e-11, "tet->hex");
    }
    let hexes = hex_grid([0.0; 3], [1.0; 3], [3, 4, 2]);
    let tets = tet_grid([0.2; 3], [0.7; 3], [2, 2, 2]);
    let (v, rep) = run(
        &hexes,
        &tets,
        InterpOptions::default(),
        field_of(&hexes, lin3),
        0.0,
    );
    assert_eq!(rep.n_not_found, 0);
    for (v, c) in v.iter().zip(&tets.coords) {
        assert_close(*v, lin3(*c), 1e-11, "hex->tet");
    }
}

#[test]
fn located_local_coordinates_reproduce_the_point() {
    let old = tet_grid([0.0; 3], [1.0; 3], [3, 3, 3]);
    let loc = Locator::new(&old, InterpOptions::default()).unwrap();
    let new = hex_grid([0.03; 3], [0.97; 3], [6, 6, 6]);
    for i in 0..new.n_nodes() {
        let p = new.coord(i as u32);
        let hit = loc.locate(loc.search_point(p), None).expect("inside");
        let el = old.element(hit.element);
        let sh = elmer_interp::shape::Shape::get(el.etype.0).unwrap();
        let mut b = vec![0.0; sh.n_nodes()];
        sh.basis(hit.uvw, &mut b);
        for (d, &pd) in p.iter().enumerate() {
            let x: f64 = b
                .iter()
                .zip(el.nodes)
                .map(|(b, &n)| b * old.coord(n)[d])
                .sum();
            assert_close(x, pd, 1e-12, "reconstructed coordinate");
        }
    }
}

#[test]
fn tree_and_exhaustive_search_agree_on_a_large_mesh() {
    // 30 x 30 x 2 triangles: well past the 8-element leaf limit.
    let old = tri_grid(0.0, 1.0, 0.0, 1.0, 30, 30);
    let new = quad_grid(0.0, 1.0, 0.0, 1.0, 37, 41);
    let f = field_of(&old, |p| lin2(p) + 0.25 * p[0] * 0.0);
    let (a, ra) = run(&old, &new, InterpOptions::default(), f.clone(), 0.0);
    let opts = InterpOptions {
        use_quadrant_tree: false,
        ..Default::default()
    };
    let (b, rb) = run(&old, &new, opts, f, 0.0);
    assert_eq!(ra.n_not_found, 0);
    assert_eq!(rb.n_not_found, 0);
    assert_eq!(ra.qtree_fails, 0);
    for (x, y) in a.iter().zip(&b) {
        assert_close(*x, *y, 1e-12, "tree vs scan");
    }
}

#[test]
fn passive_coordinate_flattens_the_search() {
    // Old 2-D mesh in the xy plane; the new 3-D nodes sit at arbitrary z.
    let old = tri_grid(0.0, 1.0, 0.0, 1.0, 4, 4);
    let new = hex_grid([0.1, 0.1, 5.0], [0.9, 0.9, 8.0], [3, 3, 2]);
    let opts = InterpOptions {
        passive_coordinate: 3,
        coordinate_system_dim: Some(3),
        ..Default::default()
    };
    let (v, rep) = run(&old, &new, opts, field_of(&old, lin2), 0.0);
    assert_eq!(rep.n_not_found, 0);
    for (v, c) in v.iter().zip(&new.coords) {
        assert_close(*v, lin2(*c), 1e-12, "passive z");
    }
    // Without the passive axis the points are far outside the old mesh.
    let (_, rep) = run(
        &old,
        &new,
        InterpOptions::default(),
        field_of(&old, lin2),
        0.0,
    );
    assert_eq!(rep.n_not_found, new.n_nodes());
}

#[test]
fn cylindric_search_uses_radius_and_height() {
    // Old mesh: (r, z) half-plane. New mesh: 3-D points; value depends on r, z.
    let old = tri_grid(0.0, 2.0, 0.0, 1.0, 8, 4);
    let mut new = hex_grid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2, 2, 2]);
    // Make the old mesh's second coordinate be z: remap its y -> stays; the
    // search point is (r, z, 0) so the old mesh y axis plays the role of z.
    new.coords.iter_mut().for_each(|c| *c = [c[0], c[1], c[2]]);
    let f = |p: [f64; 3]| 1.0 + p[0] + 2.0 * p[1];
    let opts = InterpOptions {
        cylindric: true,
        ..Default::default()
    };
    let (v, rep) = run(&old, &new, opts, field_of(&old, f), 0.0);
    assert_eq!(rep.n_not_found, 0);
    for (v, c) in v.iter().zip(&new.coords) {
        let r = (c[0] * c[0] + c[1] * c[1]).sqrt();
        assert_close(*v, 1.0 + r + 2.0 * c[2], 1e-12, "cylindric");
    }
}

#[test]
fn perm_history_columns_and_partial_hit() {
    let old = tri_grid(0.0, 1.0, 0.0, 1.0, 3, 3);
    let new = tri_grid(0.25, 0.75, 0.25, 0.75, 2, 2);
    // Old variable lives on all nodes except node 5; slots are reversed.
    let n = old.n_nodes();
    let mut perm = vec![NO_DOF; n];
    let mut slot = 0;
    for i in (0..n).rev() {
        if i != 5 {
            perm[i] = slot;
            slot += 1;
        }
    }
    let mut vals = vec![0.0; slot as usize];
    let mut prev = vec![0.0; slot as usize];
    for i in 0..n {
        if perm[i] != NO_DOF {
            vals[perm[i] as usize] = lin2(old.coords[i]);
            prev[perm[i] as usize] = 10.0 * lin2(old.coords[i]);
        }
    }
    let of = vec![FieldData::nodal("t", vals)
        .with_perm(perm)
        .with_prev(vec![prev])];
    let mk_new = || {
        vec![FieldData::nodal("t", vec![-1.0; new.n_nodes()])
            .with_prev(vec![vec![-1.0; new.n_nodes()]])]
    };

    // Strict: elements touching node 5 are skipped, the others interpolate.
    let loc = Locator::new(&old, InterpOptions::default()).unwrap();
    let mut nf = mk_new();
    interpolate(&loc, &new, &of, &mut nf, &Extras::default(), &NullLogger).unwrap();
    let mut touched_untouched = 0;
    for (i, c) in new.coords.iter().enumerate() {
        let hit = loc.locate(*c, None).unwrap();
        let el_nodes = old.element(hit.element).nodes.to_vec();
        if el_nodes.contains(&5) {
            assert_eq!(nf[0].values[i], -1.0);
            touched_untouched += 1;
        } else {
            assert_close(nf[0].values[i], lin2(*c), 1e-12, "value");
            assert_close(nf[0].prev_values[0][i], 10.0 * lin2(*c), 1e-12, "history");
        }
    }
    assert!(touched_untouched > 0);

    // Partial hit fills them too (missing node counts as zero).
    let opts = InterpOptions {
        partial_hit: true,
        ..Default::default()
    };
    let loc = Locator::new(&old, opts).unwrap();
    let mut nf = mk_new();
    interpolate(&loc, &new, &of, &mut nf, &Extras::default(), &NullLogger).unwrap();
    assert!(nf[0].values.iter().all(|&v| v != -1.0));
}

#[test]
fn body_mask_restricts_tree_search() {
    let mut old = tri_grid(0.0, 1.0, 0.0, 1.0, 4, 4);
    // Left half is body 1, right half body 2 (by first node x).
    for e in old.elems.iter_mut() {
        let x = old.coords[e.2[0] as usize][0];
        e.1 = if x < 0.5 { 1 } else { 2 };
    }
    let new = tri_grid(0.1, 0.9, 0.1, 0.9, 4, 4);
    let loc = Locator::new(
        &old,
        InterpOptions {
            try_linear_search: false,
            ..Default::default()
        },
    )
    .unwrap();
    let of = vec![FieldData::nodal("t", field_of(&old, lin2))];
    let mut nf = vec![FieldData::nodal("t", vec![-1.0; new.n_nodes()])];
    let mask = |body: i32| body == 2;
    let extras = Extras {
        body_mask: Some(&mask),
        ..Default::default()
    };
    let rep = interpolate(&loc, &new, &of, &mut nf, &extras, &NullLogger).unwrap();
    for (i, c) in new.coords.iter().enumerate() {
        if c[0] >= 0.5 {
            assert!(rep.found[i], "{c:?}");
            assert_close(nf[0].values[i], lin2(*c), 1e-12, "masked");
        } else if c[0] < 0.45 {
            assert!(!rep.found[i], "{c:?}");
        }
    }
}

#[test]
fn abort_not_found_is_an_error() {
    let old = tri_grid(0.0, 1.0, 0.0, 1.0, 2, 2);
    let new = tri_grid(0.5, 1.5, 0.0, 1.0, 2, 2);
    let opts = InterpOptions {
        abort_not_found: true,
        ..Default::default()
    };
    let loc = Locator::new(&old, opts).unwrap();
    let of = vec![FieldData::nodal("t", field_of(&old, lin2))];
    let mut nf = vec![FieldData::nodal("t", vec![0.0; new.n_nodes()])];
    let err = interpolate(&loc, &new, &of, &mut nf, &Extras::default(), &NullLogger).unwrap_err();
    assert!(matches!(err, InterpError::NotFound(n) if n > 0));
}

#[test]
fn new_mask_perm_skips_nodes() {
    let old = tri_grid(0.0, 1.0, 0.0, 1.0, 3, 3);
    let new = tri_grid(0.1, 0.9, 0.1, 0.9, 2, 2);
    let mut mask = vec![0u32; new.n_nodes()];
    mask[4] = NO_DOF;
    let loc = Locator::new(&old, InterpOptions::default()).unwrap();
    let of = vec![FieldData::nodal("t", field_of(&old, lin2))];
    let mut nf = vec![FieldData::nodal("t", vec![-1.0; new.n_nodes()])];
    let extras = Extras {
        new_mask_perm: Some(&mask),
        ..Default::default()
    };
    let rep = interpolate(&loc, &new, &of, &mut nf, &extras, &NullLogger).unwrap();
    assert!(!rep.found[4]);
    assert_eq!(nf[0].values[4], -1.0);
    assert_close(nf[0].values[0], lin2(new.coords[0]), 1e-12, "other nodes");
}

#[test]
fn projector_matches_direct_interpolation() {
    let old = tet_grid([0.0; 3], [1.0; 3], [3, 3, 3]);
    let new = hex_grid([0.05; 3], [0.95; 3], [4, 4, 4]);
    let loc = Locator::new(&old, InterpOptions::default()).unwrap();
    let p = build_projector(&loc, &new, &NodalBasis, &Extras::default(), &NullLogger).unwrap();
    assert_eq!(p.n_rows(), new.n_nodes());
    for i in 0..p.n_rows() {
        let s: f64 = p.values[p.rows[i]..p.rows[i + 1]].iter().sum();
        assert_close(s, 1.0, 1e-13, "row sum");
    }
    // Transpose has the same entries.
    assert_eq!(p.t_cols.len(), p.cols.len());
    let of = vec![FieldData::nodal("t", field_of(&old, lin3))];
    let mut via_p = vec![FieldData::nodal("t", vec![0.0; new.n_nodes()])];
    p.apply(&of, &mut via_p);
    let mut direct = vec![FieldData::nodal("t", vec![0.0; new.n_nodes()])];
    interpolate(
        &loc,
        &new,
        &of,
        &mut direct,
        &Extras::default(),
        &NullLogger,
    )
    .unwrap();
    for i in 0..new.n_nodes() {
        assert_close(
            via_p[0].values[i],
            direct[0].values[i],
            1e-12,
            "projector vs direct",
        );
        assert_close(
            via_p[0].values[i],
            lin3(new.coords[i]),
            1e-11,
            "linear exact",
        );
    }
}

#[test]
fn projector_unfound_rows_empty_or_compacted() {
    let old = tri_grid(0.0, 1.0, 0.0, 1.0, 2, 2);
    let new = tri_grid(0.5, 1.5, 0.0, 1.0, 2, 2);
    let loc = Locator::new(&old, InterpOptions::default()).unwrap();
    let p = build_projector(&loc, &new, &NodalBasis, &Extras::default(), &NullLogger).unwrap();
    assert_eq!(p.n_rows(), new.n_nodes());
    for i in 0..p.n_rows() {
        assert_eq!(p.rows[i + 1] > p.rows[i], p.found[i]);
    }
    let opts = InterpOptions {
        keep_unfound_nodes: false,
        ..Default::default()
    };
    let loc = Locator::new(&old, opts).unwrap();
    let q = build_projector(&loc, &new, &NodalBasis, &Extras::default(), &NullLogger).unwrap();
    assert_eq!(q.n_rows(), p.found.iter().filter(|f| **f).count());
}

#[test]
fn dg_variable_uses_representative_index() {
    // Two triangles sharing an edge; DG slot = 3*element + local node.
    let old = tri_grid(0.0, 1.0, 0.0, 1.0, 1, 1);
    let new = tri_grid(0.2, 0.8, 0.2, 0.8, 1, 1);
    let dg = |m: &TestMesh| {
        DgIndexes::new(
            (0..m.n_bulk())
                .map(|e| (0..3).map(|j| (3 * e + j) as u32).collect())
                .collect(),
        )
    };
    let (odg, ndg) = (dg(&old), dg(&new));
    let mut ov = vec![0.0; 3 * old.n_bulk()];
    for e in 0..old.n_bulk() {
        for (j, &n) in old.element(ElemId(e as u32)).nodes.iter().enumerate() {
            ov[3 * e + j] = lin2(old.coords[n as usize]);
        }
    }
    let of = vec![FieldData::nodal("t", ov).with_location(Location::NodesOnElements)];
    let mut nf = vec![FieldData::nodal("t", vec![-1.0; 3 * new.n_bulk()])
        .with_location(Location::NodesOnElements)];
    let loc = Locator::new(&old, InterpOptions::default()).unwrap();
    let extras = Extras {
        old_dg: Some(&odg),
        new_dg: Some(&ndg),
        ..Default::default()
    };
    interpolate(&loc, &new, &of, &mut nf, &extras, &NullLogger).unwrap();
    // First element's DG slots hold the interpolated nodal values.
    for (j, &n) in new.element(ElemId(0)).nodes.iter().enumerate() {
        assert_close(nf[0].values[j], lin2(new.coords[n as usize]), 1e-12, "dg");
    }
    // Missing tables are an error, not a panic.
    let err = interpolate(&loc, &new, &of, &mut nf, &Extras::default(), &NullLogger).unwrap_err();
    assert!(matches!(err, InterpError::MissingDgIndexes("old")));
}

#[test]
fn options_come_from_params() {
    use elmer_core::{ParamMap, Value};
    let mut sim = ParamMap::new();
    sim.insert("Interpolation Global Epsilon", Value::real(1e-6));
    sim.insert("Interpolation Max Iterations", Value::Integer(vec![3]));
    sim.insert("Try Linear Search If Qtree Fails", Value::Logical(false));
    let mut sol = ParamMap::new();
    sol.insert("Interpolation Passive Coordinate", Value::Integer(vec![3]));
    let o = InterpOptions::from_params(&sim, Some(&sol));
    assert_eq!(o.global_epsilon, 1e-6);
    assert_eq!(o.max_iterations, 3);
    assert!(!o.try_linear_search);
    assert_eq!(o.passive_coordinate, 3);
    assert_eq!(o.local_epsilon, 1e-10);
}
