use elmer_core::fe::{BasisEval, ElementBasis};
use elmer_core::mesh::ElementType;
use elmer_fe_hp::maps::p_total_dofs;
use elmer_fe_hp::{h1_p_basis, PElementBasis, PElementSpec};

const CODES: [u16; 7] = [202, 303, 404, 504, 605, 706, 808];

fn interior(ecode: u16) -> [f64; 3] {
    match ecode {
        202 => [0.3, 0.0, 0.0],
        303 => [0.1, 0.5, 0.0],
        404 => [0.3, -0.2, 0.0],
        504 => [-0.1, 0.4, 0.2],
        605 => [0.1, -0.1, 0.3],
        706 => [0.1, 0.5, 0.2],
        _ => [0.2, -0.3, 0.1],
    }
}

#[test]
fn function_counts_match_dof_formulas() {
    for &code in &CODES {
        for p in 1..=6u8 {
            let spec = PElementSpec::uniform(code, p, None, None).unwrap();
            let n = spec.n_basis().unwrap();
            let want = p_total_dofs(code, p, false).unwrap();
            assert_eq!(n, want, "element {code} p={p}");
        }
    }
}

#[test]
fn serendipity_counts_match_dof_formulas() {
    for &code in &[404u16, 706, 808] {
        for p in 2..=6u8 {
            let mut spec = PElementSpec::uniform(code, p, None, None).unwrap();
            spec.serendipity = true;
            let n = spec.n_basis().unwrap();
            let want = p_total_dofs(code, p, true).unwrap();
            assert_eq!(n, want, "serendipity element {code} p={p}");
        }
    }
}

#[test]
fn h1_p_basis_front_end_agrees_with_spec() {
    let (mut v, mut g) = (vec![0.0; 600], vec![[0.0; 3]; 600]);
    let n = h1_p_basis(ElementType(808), 3, [0.2, 0.1, -0.4], &elmer_fe_hp::maps::edge_map(808).unwrap(), &elmer_fe_hp::maps::face_map(808).unwrap(), &mut v, &mut g).unwrap();
    assert_eq!(n, p_total_dofs(808, 3, false).unwrap());
    // Too small an output is an error, not a panic.
    let mut small = vec![0.0; 10];
    let mut sg = vec![[0.0; 3]; 10];
    assert!(h1_p_basis(ElementType(808), 3, [0.0; 3], &elmer_fe_hp::maps::edge_map(808).unwrap(), &elmer_fe_hp::maps::face_map(808).unwrap(), &mut small, &mut sg).is_err());
}

#[test]
fn brick_element_basis_global_gradients() {
    // Brick [0,2] x [0,4] x [0,6]: dx/du = 1, dy/dv = 2, dz/dw = 3.
    let xyz: Vec<[f64; 3]> = [
        [0., 0., 0.],
        [2., 0., 0.],
        [2., 4., 0.],
        [0., 4., 0.],
        [0., 0., 6.],
        [2., 0., 6.],
        [2., 4., 6.],
        [0., 4., 6.],
    ]
    .to_vec();
    let b = PElementBasis::new(PElementSpec::uniform(808, 3, None, None).unwrap());
    let mut out = BasisEval::default();
    b.eval(ElementType(808), &xyz, [0.2, -0.3, 0.1], &mut out).unwrap();
    assert!((out.det_j - 1.0 * 2.0 * 3.0).abs() < 1e-13);
    // Nodal functions interpolate x: sum N_k x_k = x, gradient (1,0,0).
    let mut x = 0.0;
    let mut gx = [0.0; 3];
    for k in 0..8 {
        x += out.basis[k] * xyz[k][0];
        for c in 0..3 {
            gx[c] += out.dbasis_dx[k][c] * xyz[k][0];
        }
    }
    assert!((x - (1.0 + 0.2)).abs() < 1e-13);
    assert!((gx[0] - 1.0).abs() < 1e-13 && gx[1].abs() < 1e-13 && gx[2].abs() < 1e-13);
    assert_eq!(out.basis.len(), p_total_dofs(808, 3, false).unwrap());
}

#[test]
fn every_element_evaluates_with_an_identity_map() {
    for &code in &CODES {
        let spec = PElementSpec::uniform(code, 3, None, None).unwrap();
        let (mut v, mut g) = (vec![0.0; 700], vec![[0.0; 3]; 700]);
        let n = elmer_fe_hp::h1_p_basis_spec(&spec, interior(code), &mut v, &mut g).unwrap();
        assert!(n > 0);
        assert!(v[..n].iter().all(|x| x.is_finite()));
    }
}
