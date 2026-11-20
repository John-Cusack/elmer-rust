//! What PARPACK (ARPACK's MPI version, `p`-prefixed routines) adds to the
//! serial routines: each rank holds its rows of the Arnoldi vectors, and the
//! inner products, norms and `V^T x` products are reduced over the ranks
//! (`MPI_ALLREDUCE`, one call per reduction, as the Fortran makes them).
//!
//! The serial routines take a [`Par`] and run unchanged when it is `None`.
//!
//! Under PARPACK a rank may hold no rows (`n = 0`): it takes part in every
//! reduction with empty vectors, where Elmer leaves such a partition out of
//! the solver's communicator (`comm_active`); the argument checks accept it.

use elmer_core::Communicator;

use super::blas::{ddot, dnrm2};
use super::lapack::EPS;

/// The communicator of a PARPACK call, or `None` for serial ARPACK.
pub type Par<'a> = Option<&'a dyn Communicator>;

/// `ddot`, then `MPI_ALLREDUCE(SUM)` in PARPACK.
pub(crate) fn pddot(par: Par<'_>, n: usize, x: &[f64], y: &[f64]) -> f64 {
    let d = ddot(n, x, y);
    match par {
        None => d,
        Some(c) => {
            let mut b = [d];
            c.sum_f64(&mut b);
            b[0]
        }
    }
}

/// `dnrm2`, or PARPACK's `pdnorm2`: the local norms' maximum `m`, then
/// `m * sqrt(|sum over ranks of (local/m)**2|)`.
pub(crate) fn pdnorm2(par: Par<'_>, n: usize, x: &[f64]) -> f64 {
    let local = dnrm2(n, x);
    let Some(c) = par else { return local };
    let mut max = [local];
    c.max_f64(&mut max);
    if max[0] == 0.0 {
        return 0.0;
    }
    // `(pdnorm2/max)**2.0`, which gfortran evaluates as a product.
    let r = local / max[0];
    let mut s = [r * r];
    c.sum_f64(&mut s);
    max[0] * s[0].abs().sqrt()
}

/// `MPI_ALLREDUCE(SUM)` of `buf` in place under PARPACK; nothing serially.
pub(crate) fn psum(par: Par<'_>, buf: &mut [f64]) {
    if let Some(c) = par {
        c.sum_f64(buf);
    }
}

/// `eps23` of `_saup2`, `_naup2` and `_seupd`: `eps**(2/3)`, the exponent
/// a double-precision literal in ARPACK and a default-REAL one, `2.0/3.0`,
/// in PARPACK.
pub(crate) fn eps23(par: Par<'_>) -> f64 {
    if par.is_some() {
        EPS.powf(f64::from(2.0f32 / 3.0f32))
    } else {
        EPS.powf(2.0 / 3.0)
    }
}

#[cfg(test)]
mod tests {
    use super::super::naupd::{dnaupd, pdnaupd, NWhich};
    use super::super::neupd::{dneupd, pdneupd};
    use super::super::saupd::{dsaupd, pdsaupd, Arpack, Which};
    use super::super::seupd::{dseupd, pdseupd};
    use super::*;

    /// One rank: every reduction is the identity.
    struct Solo;

    impl Communicator for Solo {
        fn rank(&self) -> usize {
            0
        }
        fn size(&self) -> usize {
            1
        }
        fn sum_f64(&self, _buf: &mut [f64]) {}
        fn max_f64(&self, _buf: &mut [f64]) {}
        fn barrier(&self) {}
    }

    /// `y = A x`, `A` the 1D Laplacian `tridiag(-1, 2, -1)`, plus `shift`
    /// on the diagonal of row `i` times `i` (to break the symmetry of the
    /// spectrum, and of `A` when `skew` is set).
    fn apply(n: usize, skew: f64, x: &[f64], y: &mut [f64]) {
        for i in 0..n {
            let mut s = 2.0 * x[i] + 0.01 * i as f64 * x[i];
            if i > 0 {
                s -= (1.0 + skew) * x[i - 1];
            }
            if i + 1 < n {
                s -= (1.0 - skew) * x[i + 1];
            }
            y[i] = s;
        }
    }

    #[test]
    fn one_rank_parpack_finds_arpacks_eigenvalues() {
        let (n, nev, ncv) = (60usize, 4usize, 12usize);
        let mut iparam = [0i32; 11];
        iparam[0] = 1;
        iparam[2] = 300;
        iparam[6] = 1;
        let run = |par: bool| {
            let mut op = |ido: i32, ip: [usize; 3], wd: &mut [f64]| -> Result<(), String> {
                assert!(ido == -1 || ido == 1);
                let x = wd[ip[0]..ip[0] + n].to_vec();
                apply(n, 0.0, &x, &mut wd[ip[1]..ip[1] + n]);
                Ok(())
            };
            let (mut resid, mut v) = (vec![0.0; n], vec![0.0; n * ncv]);
            let (mut workd, mut workl) = (vec![0.0; 3 * n], vec![0.0; ncv * ncv + 8 * ncv]);
            let mut tol = 1e-12;
            let mut arp = Arpack::default();
            let au = if par {
                pdsaupd(&Solo, &mut op, &mut arp, false, n, Which::LM, nev, &mut tol, &mut resid, ncv, &mut v, n, iparam, &mut workd, &mut workl, 0)
            } else {
                dsaupd(&mut op, &mut arp, false, n, Which::LM, nev, &mut tol, &mut resid, ncv, &mut v, n, iparam, &mut workd, &mut workl, 0)
            }
            .unwrap();
            assert_eq!(au.info, 0);
            let (d, ierr) = if par {
                pdseupd(&Solo, 0.0, false, n, Which::LM, nev, tol, &resid, ncv, &mut v, n, &au, &mut workd, &mut workl)
            } else {
                dseupd(0.0, false, n, Which::LM, nev, tol, &resid, ncv, &mut v, n, &au, &mut workd, &mut workl)
            };
            assert_eq!(ierr, 0);
            d[..nev].to_vec()
        };
        let (s, p) = (run(false), run(true));
        for (a, b) in s.iter().zip(&p) {
            assert!((a - b).abs() <= 1e-10 * a.abs(), "{s:?} {p:?}");
        }
    }

    #[test]
    fn one_rank_parpack_finds_arpacks_nonsymmetric_eigenvalues() {
        let (n, nev, ncv) = (60usize, 4usize, 14usize);
        let mut iparam = [0i32; 11];
        iparam[0] = 1;
        iparam[2] = 300;
        iparam[6] = 1;
        let run = |par: bool| {
            let mut op = |ido: i32, ip: [usize; 3], wd: &mut [f64]| -> Result<(), String> {
                assert!(ido == -1 || ido == 1);
                let x = wd[ip[0]..ip[0] + n].to_vec();
                apply(n, 0.3, &x, &mut wd[ip[1]..ip[1] + n]);
                Ok(())
            };
            let (mut resid, mut v) = (vec![0.0; n], vec![0.0; n * ncv]);
            let (mut workd, mut workl) = (vec![0.0; 3 * n], vec![0.0; 3 * ncv * ncv + 6 * ncv]);
            let mut tol = 1e-12;
            let mut arp = Arpack::default();
            let au = if par {
                pdnaupd(&Solo, &mut op, &mut arp, false, n, NWhich::LM, nev, &mut tol, &mut resid, ncv, &mut v, n, iparam, &mut workd, &mut workl, 0)
            } else {
                dnaupd(&mut op, &mut arp, false, n, NWhich::LM, nev, &mut tol, &mut resid, ncv, &mut v, n, iparam, &mut workd, &mut workl, 0)
            }
            .unwrap();
            assert_eq!(au.info, 0);
            let r = if par {
                pdneupd(&Solo, 0.0, 0.0, false, n, NWhich::LM, nev, tol, &resid, ncv, &mut v, n, &au, &mut workd, &mut workl)
            } else {
                dneupd(0.0, 0.0, false, n, NWhich::LM, nev, tol, &resid, ncv, &mut v, n, &au, &mut workd, &mut workl)
            };
            assert_eq!(r.info, 0);
            let mut d: Vec<(f64, f64)> = (0..nev).map(|i| (r.dr[i], r.di[i])).collect();
            d.sort_by(|a, b| a.partial_cmp(b).unwrap());
            d
        };
        let (s, p) = (run(false), run(true));
        for (a, b) in s.iter().zip(&p) {
            assert!((a.0 - b.0).abs() <= 1e-10 * a.0.abs() && (a.1 - b.1).abs() <= 1e-10 * a.0.abs(), "{s:?} {p:?}");
        }
    }

    #[test]
    fn pdnorm2_of_one_rank_is_dnrm2() {
        let x = [3.0, -4.0, 12.0];
        assert_eq!(pdnorm2(Some(&Solo), 3, &x), 13.0);
        assert_eq!(pdnorm2(None, 3, &x), dnrm2(3, &x));
    }
}
