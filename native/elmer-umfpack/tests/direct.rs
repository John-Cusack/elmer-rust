//! `UmfpackLu` (the elmer-solver backend) gives the same bits as the
//! validated `umf4sym`/`umf4num`/`umf4sol` path on Elmer's CRS arrays, and
//! solves `A x = b` for the CRS matrix `A`.

#[path = "../examples/common/harness_io.rs"]
mod harness_io;

mod support;

use elmer_core::{CsrMatrix, DirectSolver};
use elmer_umfpack::UmfpackLu;

fn csr(m: &harness_io::Matrix) -> CsrMatrix<f64> {
    CsrMatrix {
        nrows: m.n as usize,
        ncols: m.n as usize,
        row_ptr: m.ap.iter().map(|&p| p as usize).collect(),
        col_idx: m.ai.iter().map(|&c| c as u32).collect(),
        values: m.ax.clone(),
        diag: None,
    }
}

fn faer_unused() -> Box<dyn DirectSolver<f64> + Send> {
    panic!("valid Elmer CRS matrices never take the fallback")
}

#[test]
fn umfpack_lu_matches_umf4sol_on_elmer_systems() {
    for name in ["bentonite", "ElasticAnisotropic2D", "ElasticStressLoad"] {
        let m = harness_io::parse(&support::data(&format!("elmer_{name}.mat")));
        let want = harness_io::run("sol", &m);
        let want_x = want
            .lines()
            .find(|l| l.starts_with("X "))
            .unwrap()
            .to_string();
        let a = csr(&m);
        let mut lu = UmfpackLu::new(faer_unused);
        lu.factor(&a).unwrap();
        let mut x = m.b.clone().unwrap();
        lu.solve_in_place(&mut x).unwrap();
        let mut got = harness_io::Out::default();
        got.dbls("X", &x);
        assert_eq!(got.0.trim_end(), want_x, "{name}");
        // A x = b for the CRS matrix
        let b = m.b.as_ref().unwrap();
        let mut rmax: f64 = 0.0;
        for i in 0..a.nrows {
            let mut s = 0.0;
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                s += a.values[k] * x[a.col_idx[k] as usize];
            }
            rmax = rmax.max((s - b[i]).abs());
        }
        let bmax = b.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
        assert!(rmax <= 1e-8 * bmax.max(1.0), "{name}: residual {rmax}");
    }
}

#[test]
fn unsorted_columns_use_the_fallback() {
    let a = CsrMatrix {
        nrows: 2,
        ncols: 2,
        row_ptr: vec![0, 2, 4],
        col_idx: vec![1, 0, 0, 1],
        values: vec![1.0, 2.0, 3.0, 4.0],
        diag: None,
    };
    let mut lu = UmfpackLu::new(|| Box::new(Dense2));
    lu.factor(&a).unwrap();
    let mut b = vec![1.0, 2.0];
    lu.solve_in_place(&mut b).unwrap();
    assert_eq!(b, vec![7.0, 7.0]);
}

/// Stand-in fallback: records that it was used by returning a fixed answer.
#[derive(Default)]
struct Dense2;
impl DirectSolver<f64> for Dense2 {
    fn factor(&mut self, _a: &CsrMatrix<f64>) -> Result<(), elmer_core::LinalgError> {
        Ok(())
    }
    fn solve_in_place(&self, b: &mut [f64]) -> Result<(), elmer_core::LinalgError> {
        b.fill(7.0);
        Ok(())
    }
}
