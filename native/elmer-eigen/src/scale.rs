//! Post-processing of real-matrix eigenpairs: `ScaleEigenVectors`
//! (EigenSolve.F90:622, the `.NOT. Matrix % COMPLEX` branch) and
//! `CheckResiduals` (EigenSolve.F90:775).

use elmer_core::CsrMatrix;
use num_complex::Complex64;

use crate::solve::{apply_mass, Mass};

/// Rescale each vector in place: to unit `M`-norm, or, with
/// `normalize_to_unity`, so that its largest component magnitude is 1. A
/// vector already within `EPSILON` of the target and a zero vector are left
/// alone, as in Elmer (which warns about the latter). Returns the indices of
/// zero-amplitude vectors.
pub fn scale_eigen_vectors(
    mass: Mass<'_>,
    vectors: &mut [Vec<Complex64>],
    normalize_to_unity: bool,
) -> Vec<usize> {
    let mut zero = Vec::new();
    for (i, x) in vectors.iter_mut().enumerate() {
        let r = if normalize_to_unity {
            x.iter().fold(0.0_f64, |m, z| m.max(z.norm()))
        } else {
            // r = sum conj(x_j) M_jl x_l, kept real as the Fortran REAL does.
            let r: f64 = match mass {
                Mass::Identity => x.iter().map(|z| z.norm_sqr()).sum(),
                Mass::Lumped(d) => x.iter().zip(d).map(|(z, di)| z.norm_sqr() * di).sum(),
                Mass::Consistent(m) => {
                    let mut s = 0.0;
                    for (j, xj) in x.iter().enumerate() {
                        for l in m.row_ptr[j]..m.row_ptr[j + 1] {
                            s += (xj.conj() * m.values[l] * x[m.col_idx[l] as usize]).re;
                        }
                    }
                    s
                }
            };
            r.sqrt()
        };
        if (r - 1.0).abs() < f64::EPSILON {
            continue;
        }
        if r > 0.0 {
            let s = 1.0 / r;
            x.iter_mut().for_each(|z| *z *= s);
        } else {
            zero.push(i);
        }
    }
    zero
}

/// `||(K - Re(lambda_i) M) Re(x_i)||_2` per pair.
pub fn check_residuals(
    k: &CsrMatrix<f64>,
    mass: Mass<'_>,
    values: &[Complex64],
    vectors: &[Vec<Complex64>],
) -> Vec<f64> {
    use elmer_core::LinearOperator;
    let n = k.nrows;
    let (mut kx, mut mx) = (vec![0.0; n], vec![0.0; n]);
    values
        .iter()
        .zip(vectors)
        .map(|(lam, x)| {
            let xr: Vec<f64> = x.iter().map(|z| z.re).collect();
            k.apply(&xr, &mut kx);
            apply_mass(mass, &xr, &mut mx);
            kx.iter()
                .zip(&mx)
                .map(|(a, b)| (a - lam.re * b).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(x: f64) -> Complex64 {
        Complex64::new(x, 0.0)
    }

    #[test]
    fn unity_and_mass_scaling() {
        let mut v = vec![vec![c(1.0), c(-4.0), c(2.0)]];
        scale_eigen_vectors(Mass::Identity, &mut v, true);
        assert_eq!(v[0], vec![c(0.25), c(-1.0), c(0.5)]);

        let d = [2.0, 1.0, 0.5];
        let mut v = vec![vec![c(1.0), c(1.0), c(2.0)]];
        scale_eigen_vectors(Mass::Lumped(&d), &mut v, false);
        let r: f64 = v[0].iter().zip(&d).map(|(z, di)| z.norm_sqr() * di).sum();
        assert!((r - 1.0).abs() < 1e-15);

        let mut v = vec![vec![c(0.0); 3]];
        assert_eq!(scale_eigen_vectors(Mass::Identity, &mut v, false), vec![0]);
    }

    #[test]
    fn consistent_mass_norm() {
        // M = [[2, 1], [1, 2]]
        let m = CsrMatrix {
            nrows: 2,
            ncols: 2,
            row_ptr: vec![0, 2, 4],
            col_idx: vec![0, 1, 0, 1],
            values: vec![2.0, 1.0, 1.0, 2.0],
            diag: Some(vec![0, 3]),
        };
        let mut v = vec![vec![c(1.0), c(1.0)]];
        scale_eigen_vectors(Mass::Consistent(&m), &mut v, false);
        assert!((v[0][0].re - 1.0 / 6f64.sqrt()).abs() < 1e-15);
    }
}
