//! The axisymmetric view factors of `fem/tests/radiation` against the
//! `ViewFactors.dat` the pristine `ViewFactors` program writes for it
//! (`ElmerGrid 1 2 radiation`, then `ViewFactors` in the mesh directory).

use std::fs;
use std::path::PathBuf;

use elmer_viewfactors::{
    line_direction, normalize_factors, saved_rows, view_factors_axis, ParentInfo,
};

fn data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/radiation")
}

fn lines(name: &str) -> Vec<Vec<f64>> {
    fs::read_to_string(data().join(name))
        .unwrap()
        .lines()
        .map(|l| l.split_whitespace().map(|t| t.parse().unwrap()).collect())
        .collect()
}

/// (surf, coord, areas) of the boundary elements tagged 1 or 2, the two
/// `Radiation = Diffuse Gray` conditions with `Radiation Target Body = -1`.
fn surfaces() -> (Vec<usize>, Vec<f64>, Vec<f64>) {
    let nodes: Vec<[f64; 3]> = lines("mesh.nodes")
        .iter()
        .map(|l| [l[2], l[3], l[4]])
        .collect();
    let bulk = lines("mesh.elements");
    let parent_info = |id: usize, nds: &[usize]| -> Option<ParentInfo> {
        if id == 0 {
            return None;
        }
        let e = &bulk[id - 1];
        let mut last = None;
        for &n in &e[3..] {
            if !nds.contains(&(n as usize)) {
                last = Some(nodes[n as usize - 1]);
            }
        }
        Some(ParentInfo {
            body: e[1] as i32,
            node: last.unwrap(),
        })
    };
    let coord: Vec<f64> = nodes.iter().flat_map(|n| [n[0], n[1]]).collect();

    let mut surf = Vec::new();
    let mut areas = Vec::new();
    for l in lines("mesh.boundary") {
        if l[1] != 1.0 && l[1] != 2.0 {
            continue;
        }
        let (n1, n2) = (l[5] as usize, l[6] as usize);
        let (a, b) = (nodes[n1 - 1], nodes[n2 - 1]);
        let left = parent_info(l[2] as usize, &[n1, n2]);
        let right = parent_info(l[3] as usize, &[n1, n2]);
        let dir = line_direction(a, b, left, right, Some(-1), (false, false)).unwrap();
        if dir > 0.0 {
            surf.extend([n2 - 1, n1 - 1]);
        } else {
            surf.extend([n1 - 1, n2 - 1]);
        }
        let len = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
        areas.push(len * (a[0] + b[0]).abs() / 2.0);
    }
    (surf, coord, areas)
}

#[test]
fn matches_the_viewfactors_program() {
    let (surf, coord, areas) = surfaces();
    let n = areas.len();
    assert_eq!(n, 39);
    let mut vf = view_factors_axis(&surf, &coord, 1, true);

    // Every ring of a closed cavity sees all the others.
    normalize_factors(&mut vf, &areas, false);
    let rows = saved_rows(&vf, n, n, 1e-20, false);

    // ViewFactors.dat: the count of a row, then `row column factor` lines.
    let dat = lines("ViewFactors.dat");
    let mut at = 0;
    let mut worst: f64 = 0.0;
    for (i, row) in rows.iter().enumerate() {
        let k = dat[at][0] as usize;
        at += 1;
        assert_eq!(k, row.len(), "row {i}");
        for (j, v) in row {
            let want = &dat[at];
            at += 1;
            assert_eq!((want[0] as usize, want[1] as usize), (i + 1, j + 1));
            worst = worst.max((v - want[2]).abs() / want[2].abs().max(1e-12));
        }
    }
    assert_eq!(at, dat.len());
    assert!(worst < 1e-9, "worst relative difference {worst:e}");
}
