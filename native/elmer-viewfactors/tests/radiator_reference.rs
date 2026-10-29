//! Radiator factors of `fem/tests/radiator3d` (two point radiators inside a
//! cube, `Radiation Target Body = 1`) against the `RadiatorFactors.dat` the
//! pristine `ViewFactors -radiators` writes for it.

use std::fs;
use std::path::PathBuf;

use elmer_viewfactors::{
    radiator_factors, saved_rows, surface_direction, ParentInfo, RadiatorParams, RadiatorSurface,
};

fn data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/radiator3d")
}

fn lines(name: &str) -> Vec<Vec<f64>> {
    fs::read_to_string(data().join(name))
        .unwrap()
        .lines()
        .map(|l| l.split_whitespace().map(|t| t.parse().unwrap()).collect())
        .collect()
}

/// The oriented surfaces of the boundary elements tagged 1, 2 or 3.
fn surfaces() -> Vec<RadiatorSurface> {
    let nodes: Vec<[f64; 3]> = lines("mesh.nodes")
        .iter()
        .map(|l| [l[2], l[3], l[4]])
        .collect();
    let bulk = lines("mesh.elements");
    let parent = |id: usize, own: &[usize]| -> Option<ParentInfo> {
        if id == 0 {
            return None;
        }
        let e = &bulk[id - 1];
        let last = e[3..]
            .iter()
            .rev()
            .map(|&n| n as usize)
            .find(|n| !own.contains(n))?;
        Some(ParentInfo {
            body: e[1] as i32,
            node: nodes[last - 1],
        })
    };
    let mut out = Vec::new();
    for l in lines("mesh.boundary") {
        if ![1.0, 2.0, 3.0].contains(&l[1]) {
            continue;
        }
        let own: Vec<usize> = l[5..].iter().map(|&n| n as usize).collect();
        let p: Vec<[f64; 3]> = own.iter().map(|&n| nodes[n - 1]).collect();
        // NormalVector at the element centre (u = v = 0 on [-1, 1]^2).
        let (xu, xv) = (
            (0..3)
                .map(|i| ((p[1][i] + p[2][i]) - (p[0][i] + p[3][i])) / 4.0)
                .collect::<Vec<_>>(),
            (0..3)
                .map(|i| ((p[2][i] + p[3][i]) - (p[0][i] + p[1][i])) / 4.0)
                .collect::<Vec<_>>(),
        );
        let c = [
            xu[1] * xv[2] - xu[2] * xv[1],
            xu[2] * xv[0] - xu[0] * xv[2],
            xu[0] * xv[1] - xu[1] * xv[0],
        ];
        let m = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        let nrm = [c[0] / m, c[1] / m, c[2] / m];
        let centre = [0, 1, 2].map(|i| p.iter().map(|q| q[i]).sum::<f64>() / 4.0);
        let dir = surface_direction(
            nrm,
            centre,
            parent(l[2] as usize, &own),
            parent(l[3] as usize, &own),
            Some(1),
            (false, false),
        )
        .unwrap();
        let normal = if dir > 0.0 { nrm } else { nrm.map(|c| -c) };
        out.push(RadiatorSurface { nodes: p, normal });
    }
    out
}

#[test]
fn matches_the_viewfactors_program() {
    let s = surfaces();
    assert_eq!(s.len(), 600);
    // Viewfactor area tolerance 0.001, factor tolerance 0.001, 8 rays.
    let params = RadiatorParams::from_keywords(Some(0.001), Some(0.001), Some(8), None, None, None);
    let radiators = [[0.15, 0.15, 0.15], [0.85, 0.85, 0.85]];
    let mut f = radiator_factors(&s, &radiators, &params).unwrap();
    let n = s.len();
    for row in f.chunks_mut(n) {
        let sum: f64 = row.iter().sum();
        row.iter_mut().for_each(|x| *x /= sum);
    }
    let rows = saved_rows(&f, n, 2, 1e-20, false);

    let dat = lines("RadiatorFactors.dat");
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
            worst = worst.max((v - want[2]).abs() / want[2].abs());
        }
    }
    assert!(worst < 1e-9, "worst relative difference {worst:e}");
}
