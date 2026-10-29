//! The three assembled systems under `tests/fixtures` come from the pristine
//! Fortran `ElmerSolver` (UMFPACK path); see the README there.
//!
//! Elmer writes its solution vector with `ES15.6`, so agreement with it is
//! checked to 1e-5 (six digits). Tight (~1e-12) checks are between the
//! independent backends here and by residual.
mod common;
use common::*;

use elmer_core::DirectSolver;
use elmer_direct::{set_single_thread, BandLu, FaerCholesky, FaerLu};

fn solve(s: &mut dyn DirectSolver<f64>, d: &Dump) -> Vec<f64> {
    s.factor(&d.a).unwrap();
    let mut x = d.b.clone();
    s.solve_in_place(&mut x).unwrap();
    x
}

fn max_rel(x: &[f64], y: &[f64]) -> f64 {
    let scale = y.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    x.iter().zip(y).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max) / scale
}

#[test]
fn backends_agree_with_fortran_umfpack_solution() {
    set_single_thread();
    for name in ["bodydir", "bodyload", "pfem"] {
        let d = fixture(name);
        assert_eq!(d.a.nrows, 441, "{name}");
        let sol = d.sol.as_ref().unwrap();

        let lu = solve(&mut FaerLu::new(), &d);
        let band = solve(&mut BandLu::new(), &d);

        let r_lu = rel_residual(&d.a, &lu, &d.b);
        let r_band = rel_residual(&d.a, &band, &d.b);
        let vs_fortran = max_rel(&lu, sol);
        let lu_vs_band = max_rel(&lu, &band);
        println!("{name}: n=441 residual lu {r_lu:.2e} band {r_band:.2e}; max rel diff lu vs fortran {vs_fortran:.2e}; lu vs band {lu_vs_band:.2e}");

        assert!(r_lu < 1e-12 && r_band < 1e-12, "{name}: residuals {r_lu:e} {r_band:e}");
        assert!(vs_fortran < 1e-5, "{name}: vs Fortran {vs_fortran:e}");
        assert!(lu_vs_band < 1e-8, "{name}: lu vs band {lu_vs_band:e}");

        if is_symmetric(&d.a) {
            let ch = solve(&mut FaerCholesky::new(), &d);
            println!("{name}: symmetric, cholesky vs lu {:.2e}", max_rel(&ch, &lu));
            assert!(max_rel(&ch, &lu) < 1e-8, "{name}");
        }
    }
}
