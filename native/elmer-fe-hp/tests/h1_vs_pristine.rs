//! Per-function comparison of the p-hierarchic H1 basis (values, global
//! gradients, detJ) against what the pristine ElmerSolver's scalar
//! `ElementInfo` returns for p elements, captured with gdb at 17 digits (see
//! `tests/fixtures/README.md`). Covers p=2..5, 2D and 3D elements, identity and
//! shuffled global node numbering (so edge/face directions and tetra types are
//! exercised) and the serendipity scheme.
//!
//! Not bitwise: `ElementInfo` evaluates `PElementBase.F90`, this crate ports
//! `H1Basis.F90` (what `ElementInfoVec` calls). The two agree to rounding
//! except pyramids, which are bit-for-bit identical, and the wedge square-face
//! ordering, handled by `PElementSpec::scalar_face_order`.

use elmer_core::fe::{BasisEval, ElementBasis};
use elmer_core::mesh::ElementType;
use elmer_fe_hp::maps::{edge_direction, face_direction, tetra_edge_direction, tetra_face_direction};
use elmer_fe_hp::{PElementBasis, PElementSpec};

struct Rec {
    code: u16,
    ttype: i32,
    p: u8,
    sd: bool,
    ids: Vec<i64>,
    xyz: Vec<[f64; 3]>,
    uvw: [f64; 3],
    detj: f64,
    basis: Vec<f64>,
    grad: Vec<[f64; 3]>,
}

fn parse(text: &str) -> Vec<Rec> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let t: Vec<f64> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
            let (code, nn, n, tt, p, sd) = (t[0] as u16, t[1] as usize, t[2] as usize, t[3] as i32, t[4] as u8, t[5] != 0.0);
            let mut i = 6;
            let (mut ids, mut xyz) = (vec![], vec![]);
            for _ in 0..nn {
                ids.push(t[i] as i64);
                xyz.push([t[i + 1], t[i + 2], t[i + 3]]);
                i += 4;
            }
            let uvw = [t[i], t[i + 1], t[i + 2]];
            let detj = t[i + 3];
            i += 4;
            let basis = t[i..i + n].to_vec();
            i += n;
            let grad = (0..n).map(|k| [t[i + 3 * k], t[i + 3 * k + 1], t[i + 3 * k + 2]]).collect();
            Rec { code, ttype: tt, p, sd, ids, xyz, uvw, detj, basis, grad }
        })
        .collect()
}

fn spec_for(r: &Rec) -> PElementSpec {
    let ne = match r.code { 808 => 12, 706 => 9, 605 => 8, 504 => 6, 404 => 4, 303 => 3, _ => 0 };
    let nf = match r.code { 808 => 6, 706 | 605 => 5, 504 => 4, _ => 0 };
    let nf_dir = if nf == 0 { 0 } else { nf };
    let (ed, fd): (Vec<[usize; 2]>, Vec<[usize; 4]>) = if r.code == 504 {
        (tetra_edge_direction(r.ttype).to_vec(), tetra_face_direction(r.ttype).to_vec())
    } else {
        (
            edge_direction(r.code, ne, &r.ids).unwrap().to_vec(),
            if nf_dir == 0 { vec![[0; 4]; 6] } else { face_direction(r.code, nf, &r.ids).unwrap().to_vec() },
        )
    };
    let mut s = PElementSpec::uniform(r.code, r.p, Some(&ed), Some(&fd)).unwrap();
    s.serendipity = r.sd;
    s
}


fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * (1.0 + b.abs())
}

#[test]
fn h1_matches_pristine_per_function() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/pristine");
    let mut files: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.file_name().unwrap().to_str().unwrap().starts_with('h'))
        .collect();
    files.sort();
    assert!(files.len() >= 12, "fixtures missing");
    let mut checked = 0usize;
    for file in files {
        let recs = parse(&std::fs::read_to_string(&file).unwrap());
        for r in &recs {
            let b = PElementBasis::new(spec_for(r));
            let mut out = BasisEval::default();
            b.eval(ElementType(r.code), &r.xyz, r.uvw, &mut out).unwrap();
            let ctx = format!("{:?} element {} p={} sd={} uvw={:?}", file.file_name().unwrap(), r.code, r.p, r.sd, r.uvw);
            assert_eq!(out.basis.len(), r.basis.len(), "{ctx}: function count");
            assert!(close(out.det_j, r.detj, 1e-14), "{ctx}: detJ {} vs {}", out.det_j, r.detj);
            for k in 0..r.basis.len() {
                if r.code == 605 {
                    assert_eq!(out.basis[k], r.basis[k], "{ctx}: pyramid fn {k} not bitwise");
                } else {
                    assert!(close(out.basis[k], r.basis[k], 5e-13), "{ctx}: fn {k} value {} vs {}", out.basis[k], r.basis[k]);
                }
                for c in 0..3 {
                    assert!(close(out.dbasis_dx[k][c], r.grad[k][c], 5e-13), "{ctx}: fn {k} d/dx{c} {} vs {}", out.dbasis_dx[k][c], r.grad[k][c]);
                }
                checked += 1;
            }
        }
    }
    assert!(checked > 10_000, "only {checked} function evaluations checked");
}
