//! The element kernels on single elements, with hand-computed matrices.

use elmer_core::{ElementType, LocalSystem};
use elmer_fe::{ElmerQuadrature, LagrangeBasis};
use elmer_physics_advdiff::{
    advdiff_boundary_local, advdiff_bulk_local, bubble_gauss_points, condensate, eval_with_bubbles,
    model_pde_boundary_local, model_pde_bulk_local, AdvDiffConfig, AdvDiffNodal, Fe,
    ModelPdeBcNodal, ModelPdeNodal, Stabilization,
};

const TRI: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];

fn eye(k: f64) -> [[f64; 3]; 3] {
    let mut m = [[0.0; 3]; 3];
    for (i, r) in m.iter_mut().enumerate() {
        r[i] = k;
    }
    m
}

fn nodal(n: usize, ct: f64, c1: f64, k: f64, v: [f64; 3], load: f64) -> AdvDiffNodal {
    AdvDiffNodal {
        ct: vec![ct; n],
        c1: vec![c1; n],
        c2: vec![eye(k); n],
        velocity: vec![v; n],
        load: vec![load; n],
        ..Default::default()
    }
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-14
}

#[test]
fn diffusion_mass_and_source_on_the_unit_triangle() {
    let basis = LagrangeBasis::new(2);
    let fe = Fe::new(&basis, &ElmerQuadrature);
    let mut out = LocalSystem::zeros(0);
    advdiff_bulk_local(
        &fe,
        ElementType(303),
        &TRI,
        &nodal(3, 2.0, 0.0, 3.0, [0.0; 3], 6.0),
        &AdvDiffConfig::new(2),
        &mut out,
    )
    .unwrap();
    // k * area * grad N_p . grad N_q with area 1/2.
    let k = [[2.0, -1.0, -1.0], [-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]];
    let mass = out.mass.as_ref().unwrap();
    for p in 0..3 {
        for q in 0..3 {
            assert!(close(out.stiff[p * 3 + q], 1.5 * k[p][q]));
            let m = if p == q { 2.0 } else { 1.0 } * 2.0 * 0.5 / 12.0;
            assert!(close(mass[p * 3 + q], m));
        }
        assert!(close(out.force[p], 6.0 * 0.5 / 3.0));
    }
}

#[test]
fn convection_row_sums_and_orientation() {
    // (C1 u . grad N_q, N_p): each row p sums to C1 (u . grad(sum N)) * ... = 0,
    // and each column q sums to C1 (u . grad N_q) * area.
    let basis = LagrangeBasis::new(2);
    let fe = Fe::new(&basis, &ElmerQuadrature);
    let mut out = LocalSystem::zeros(0);
    advdiff_bulk_local(
        &fe,
        ElementType(303),
        &TRI,
        &nodal(3, 0.0, 2.0, 0.0, [1.0, 0.0, 0.0], 0.0),
        &AdvDiffConfig::new(2),
        &mut out,
    )
    .unwrap();
    // u = (1, 0): grad N = (-1,-1), (1,0), (0,1) so u.grad N = -1, 1, 0.
    let du = [-1.0, 1.0, 0.0];
    for p in 0..3 {
        for q in 0..3 {
            assert!(close(out.stiff[p * 3 + q], 2.0 * du[q] * 0.5 / 3.0));
        }
    }
}

#[test]
fn bubble_basis_vanishes_on_the_boundary_and_condenses_a_pure_source() {
    let basis = LagrangeBasis::new(2);
    let mut be = elmer_core::BasisEval::default();
    // On the edge v = 0 the bubble functions are zero.
    eval_with_bubbles(&basis, ElementType(303), &TRI, [0.3, 0.0, 0.0], 2, &mut be).unwrap();
    assert_eq!(be.basis.len(), 6);
    assert!(be.basis[3..].iter().all(|b| b.abs() < 1e-15));
    // Nodal functions sum to one.
    assert!(close(be.basis[..3].iter().sum(), 1.0));

    assert_eq!(bubble_gauss_points(ElementType(303)), Some(11));
    assert_eq!(bubble_gauss_points(ElementType(404)), Some(16));

    // Condensation: K = [[2, 1], [1, 4]], n = nb = 1, f = [1, 2]
    // => K' = 2 - 1/4 = 1.75 and f' = 1 - 2/4 = 0.5.
    let mut k = vec![2.0, 1.0, 1.0, 4.0];
    let mut f = vec![1.0, 2.0];
    condensate(1, 1, &mut k, &mut f, None).unwrap();
    assert!(close(k[0], 1.75) && close(f[0], 0.5));
}

#[test]
fn bubbles_apply_only_with_convection() {
    let basis = LagrangeBasis::new(2);
    let fe = Fe::new(&basis, &ElmerQuadrature);
    let mut cfg = AdvDiffConfig::new(2);
    cfg.stabilization = Stabilization::Bubbles;
    let mut out = LocalSystem::zeros(0);
    advdiff_bulk_local(
        &fe,
        ElementType(303),
        &TRI,
        &nodal(3, 1.0, 0.0, 1.0, [0.0; 3], 1.0),
        &cfg,
        &mut out,
    )
    .unwrap();
    assert_eq!(out.n, 3);
    advdiff_bulk_local(
        &fe,
        ElementType(303),
        &TRI,
        &nodal(3, 1.0, 1.0, 1.0, [1.0, 0.0, 0.0], 1.0),
        &cfg,
        &mut out,
    )
    .unwrap();
    assert_eq!(out.n, 6);
    assert_eq!(out.stiff.len(), 36);
}

#[test]
fn supg_without_hessian_is_unsupported() {
    let basis = LagrangeBasis::new(2);
    let fe = Fe::new(&basis, &ElmerQuadrature);
    let mut cfg = AdvDiffConfig::new(2);
    cfg.stabilization = Stabilization::Supg {
        h_k: 1.0,
        m_k: 1.0 / 3.0,
    };
    let mut out = LocalSystem::zeros(0);
    let e = advdiff_bulk_local(
        &fe,
        ElementType(303),
        &TRI,
        &nodal(3, 1.0, 1.0, 1.0, [1.0, 0.0, 0.0], 1.0),
        &cfg,
        &mut out,
    )
    .unwrap_err();
    assert!(e.to_string().contains("second derivatives"), "{e}");
}

#[test]
fn boundary_mass_transfer_on_a_unit_line() {
    let basis = LagrangeBasis::new(2);
    let fe = Fe::new(&basis, &ElmerQuadrature);
    let xyz = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let mut out = LocalSystem::zeros(0);
    advdiff_boundary_local(
        &fe,
        ElementType(202),
        &xyz,
        &[3.0, 3.0],
        &[1.0, 1.0],
        &mut out,
    )
    .unwrap();
    // Length 2: int alpha N_p N_q = alpha L/6 [[2,1],[1,2]], int load N_q = L/2.
    assert!(close(out.stiff[0], 3.0 * 2.0 / 6.0 * 2.0));
    assert!(close(out.stiff[1], 3.0 * 2.0 / 6.0));
    assert!(close(out.force[0], 1.0) && close(out.force[1], 1.0));
}

#[test]
fn model_pde_terms_on_the_unit_triangle() {
    let basis = LagrangeBasis::new(2);
    let fe = Fe::new(&basis, &ElmerQuadrature);
    let n = 3;
    let nod = ModelPdeNodal {
        diffusion: vec![1.0; n],
        convection: vec![2.0; n],
        reaction: vec![4.0; n],
        time: vec![6.0; n],
        velocity: vec![[1.0, 0.0, 0.0]; n],
        source: vec![12.0; n],
    };
    let mut out = LocalSystem::zeros(0);
    model_pde_bulk_local(&fe, ElementType(303), &TRI, &nod, 2, None, &mut out).unwrap();
    let k = [[2.0, -1.0, -1.0], [-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]];
    let du = [-1.0, 1.0, 0.0];
    let mass = out.mass.as_ref().unwrap();
    for p in 0..3 {
        for q in 0..3 {
            let m = (if p == q { 2.0 } else { 1.0 }) / 24.0;
            let want = 0.5 * k[p][q] + 2.0 * du[q] * 0.5 / 3.0 + 4.0 * m;
            assert!(close(out.stiff[p * 3 + q], want), "{p}{q}");
            assert!(close(mass[p * 3 + q], 6.0 * m));
        }
        assert!(close(out.force[p], 12.0 * 0.5 / 3.0));
    }

    let line = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let bc = ModelPdeBcNodal {
        flux: vec![10.0; 2],
        robin: vec![10.0; 2],
        external: vec![5.0; 2],
    };
    model_pde_boundary_local(&fe, ElementType(202), &line, &bc, &mut out).unwrap();
    assert!(close(out.stiff[0], 10.0 / 3.0) && close(out.stiff[1], 10.0 / 6.0));
    // (flux + c u_ext) L / 2 = 60 / 2.
    assert!(close(out.force[0], 30.0) && close(out.force[1], 30.0));
}
