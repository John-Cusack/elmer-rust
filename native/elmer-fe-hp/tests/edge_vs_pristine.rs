//! Per-function comparison of the edge (H(curl)) bases (values, curls, detJ)
//! against what the pristine ElmerSolver's `ElementInfo(..., EdgeBasis,
//! RotBasis)` returns, captured with gdb at 17 digits (see
//! `tests/fixtures/README.md`), for the lowest-order non-Piola, lowest-order
//! Piola and second-order configurations of `ElementUnitTest_edge*`, with
//! identity and shuffled global node numbering. The Piola paths are bitwise
//! equal; the non-Piola 706/808 paths differ by at most one rounding.

use elmer_core::mesh::ElementType;
use elmer_fe_hp::edge::{edge_element_info, get_edge_basis};

struct Rec {
    code: u16,
    mode: u8,
    ids: Vec<i64>,
    xyz: Vec<[f64; 3]>,
    uvw: [f64; 3],
    detj: f64,
    edge: Vec<[f64; 3]>,
    curl: Vec<[f64; 3]>,
}

fn parse(text: &str) -> Vec<Rec> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let t: Vec<f64> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
            let (code, nn, nf, mode) = (t[0] as u16, t[1] as usize, t[2] as usize, t[3] as u8);
            let mut i = 4;
            let (mut ids, mut xyz) = (vec![], vec![]);
            for _ in 0..nn {
                ids.push(t[i] as i64);
                xyz.push([t[i + 1], t[i + 2], t[i + 3]]);
                i += 4;
            }
            let uvw = [t[i], t[i + 1], t[i + 2]];
            let detj = t[i + 3];
            i += 4;
            let v3 = |s: usize| -> Vec<[f64; 3]> { (0..nf).map(|k| [t[s + 3 * k], t[s + 3 * k + 1], t[s + 3 * k + 2]]).collect() };
            let edge = v3(i);
            let curl = v3(i + 3 * nf);
            Rec { code, mode, ids, xyz, uvw, detj, edge, curl }
        })
        .collect()
}


fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * (1.0 + b.abs())
}

#[test]
fn edge_matches_pristine_per_function() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/pristine");
    let mut files: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.file_name().unwrap().to_str().unwrap().starts_with("e_"))
        .collect();
    files.sort();
    assert_eq!(files.len(), 6, "fixtures missing");
    let mut checked = 0usize;
    for file in files {
        for r in parse(&std::fs::read_to_string(&file).unwrap()) {
            let et = ElementType(r.code);
            let ev = match r.mode {
                0 => get_edge_basis(et, &r.xyz, &r.ids, r.uvw, 3),
                m => edge_element_info(et, &r.xyz, &r.ids, r.uvw, m, 3),
            }
            .unwrap();
            let ctx = format!("{:?} element {} mode {} uvw={:?}", file.file_name().unwrap(), r.code, r.mode, r.uvw);
            assert_eq!(ev.edge.len(), r.edge.len(), "{ctx}: function count");
            // Everything is bitwise except the lowest-order non-Piola wedge/brick.
            let exact = !(r.mode == 0 && (r.code == 706 || r.code == 808));
            let eq = |a: f64, b: f64| if exact { a == b } else { close(a, b, 1e-15) };
            assert!(eq(ev.det, r.detj), "{ctx}: detJ {} vs {}", ev.det, r.detj);
            for k in 0..r.edge.len() {
                for c in 0..3 {
                    assert!(eq(ev.edge[k][c], r.edge[k][c]), "{ctx}: fn {k} edge[{c}] {} vs {}", ev.edge[k][c], r.edge[k][c]);
                    assert!(eq(ev.curl[k][c], r.curl[k][c]), "{ctx}: fn {k} curl[{c}] {} vs {}", ev.curl[k][c], r.curl[k][c]);
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 6_000, "only {checked} components checked");
}
