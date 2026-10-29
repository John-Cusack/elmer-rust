//! NeoHookean2D from its pristine `.sif`: Cook's membrane in plane strain with
//! the compressible neo-Hookean material, `Pseudo-Traction`, Newton iteration.
//! The pristine ElmerSolver reaches `Reference Norm = 2.0045949` in five
//! iterations; the displacement norms it printed per iteration are
//! 2.1536113, 2.0601079, 2.0049654, 2.0045957 and 2.0045949.

mod support;

#[test]
fn neohookean2d_matches_the_reference_norm() {
    let run = support::run_sif("neohookean2d", "cooksmembrane.sif", |_| {}).unwrap();
    let norm = run.norm();
    println!("NeoHookean2D: norm {norm:.17e}");
    assert!((norm - 2.0045949).abs() / 2.0045949 < 1.0e-5, "norm {norm}");
    assert_eq!(run.passed, Some(true));
}
