//! The port of ARPACK's complex `ZNEUPD` (after `ZNAUPD`) against the
//! pristine build's `libarpack.so`, bit for bit.
//!
//! Each line of `tests/data/zneupd_pristine.txt` is `znaupd n nev ncv which
//! tol mxiter mode OP B neupd SIGMA<TAB>info ierr NC D Z WORKL`: the problem
//! of `tests/znaupd_pristine.rs`, then what the library's `zneupd_(.true.,
//! 'A', ...)` returned with `Z` the same array as `V`, as Elmer's
//! `ArpackEigenSolveComplex` calls it. Recorded
//! by driving the library's `znaupd_` through ctypes on a random dense
//! problem (`mode` 1: `OP = A`; `mode` 3: `OP = A^-1`, `B` Hermitian
//! positive definite, as Elmer's shift-invert calls hand it), complex
//! entries as pairs of f64 bit patterns. The problems run in one process in
//! file order, so `zgetv0`'s `SAVE`d seed and the `SAVE`d `smlnum` carry
//! over from one to the next as they do in the library, and here in one
//! [`ZArpack`].
//!
//! `WORKL` is compared by value: the library can leave some exact zeros of
//! its scratch (the triangular eigenvector matrix, the Ritz estimates) with
//! the other sign. `D` and the Ritz vectors are compared bit for bit.

use elmer_eigen::arpack::naupd::NWhich;
use elmer_eigen::arpack::znaupd::{znaupd, ZArpack};
use elmer_eigen::arpack::zneupd::zneupd;
use num_complex::Complex64 as C;

fn f(t: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(t, 16).unwrap())
}

fn zv(t: &str) -> Vec<C> {
    let v: Vec<f64> = t.split(',').map(f).collect();
    v.chunks(2).map(|p| C::new(p[0], p[1])).collect()
}

fn hz(v: &[C]) -> String {
    if v.is_empty() {
        return "-".into();
    }
    v.iter().map(|z| format!("{:016x},{:016x}", z.re.to_bits(), z.im.to_bits())).collect::<Vec<_>>().join(",")
}

#[test]
fn zneupd_matches_the_pristine_libarpack_bit_for_bit() {
    let data = include_str!("data/zneupd_pristine.txt");
    let mut zarp = ZArpack::default();
    let mut cases = 0;
    for line in data.lines() {
        let (input, want) = line.split_once('\t').unwrap();
        let t: Vec<&str> = input.split_whitespace().collect();
        assert_eq!(t[0], "znaupd");
        let u = |k: usize| t[k].parse::<usize>().unwrap();
        let (n, nev, ncv) = (u(1), u(2), u(3));
        let which = match t[4] {
            "LM" => NWhich::LM,
            "SM" => NWhich::SM,
            "LR" => NWhich::LR,
            "SR" => NWhich::SR,
            "LI" => NWhich::LI,
            _ => NWhich::SI,
        };
        let (mut tol, mxiter, mode) = (f(t[5]), t[6].parse::<i32>().unwrap(), t[7].parse::<i32>().unwrap());
        let (opm, b) = (zv(t[8]), zv(t[9]));
        let bmat_g = mode != 1;
        let mv = |m: &[C], x: &[C], y: &mut [C]| {
            for r in 0..n {
                let mut s = C::new(0.0, 0.0);
                for c in 0..n {
                    s += m[r + c * n] * x[c];
                }
                y[r] = s;
            }
        };
        let mut op = |ido: i32, ip: [usize; 3], wd: &mut [C]| -> Result<(), String> {
            let mut y = vec![C::new(0.0, 0.0); n];
            match (ido, bmat_g) {
                (-1, true) => {
                    let mut bx = vec![C::new(0.0, 0.0); n];
                    mv(&b, &wd[ip[0]..ip[0] + n], &mut bx);
                    mv(&opm, &bx, &mut y);
                }
                (1, true) => mv(&opm, &wd[ip[2]..ip[2] + n], &mut y),
                (-1 | 1, false) => mv(&opm, &wd[ip[0]..ip[0] + n], &mut y),
                (2, _) => mv(&b, &wd[ip[0]..ip[0] + n], &mut y),
                _ => return Err(format!("ido {ido}")),
            }
            wd[ip[1]..ip[1] + n].copy_from_slice(&y);
            Ok(())
        };
        let mut resid = vec![C::new(0.0, 0.0); n];
        let mut v = vec![C::new(0.0, 0.0); n * ncv];
        let mut workd = vec![C::new(0.0, 0.0); 3 * n];
        let mut workl = vec![C::new(0.0, 0.0); 3 * ncv * ncv + 5 * ncv];
        let mut iparam = [0i32; 11];
        iparam[0] = 1;
        iparam[2] = mxiter;
        iparam[6] = mode;
        let r = znaupd(&mut op, &mut zarp, bmat_g, n, which, nev, &mut tol, &mut resid, ncv, &mut v, n, iparam, &mut workd, &mut workl, 0)
            .unwrap();
        let w: Vec<&str> = want.split(' ').collect();
        if w[1] == "99" {
            assert_eq!(r.info.to_string(), w[0]);
            cases += 1;
            continue;
        }
        let sigma = zv(t[11])[0];
        let e = zneupd(true, sigma, bmat_g, n, which, nev, tol, &resid, ncv, &mut v, n, &r, &mut workd, &mut workl);
        assert_eq!(format!("{} {}", r.info, e.info), w[..2].join(" "), "{}", &input[..30]);
        assert_eq!(hz(&e.d[..e.nconv]), w[3], "D, {}", &input[..30]);
        assert_eq!(hz(&v[..n * e.nconv]), w[4], "Z, {}", &input[..30]);
        let wl = zv(w[5]);
        assert!(wl.len() == workl.len() && wl.iter().zip(&workl).all(|(a, b)| a == b), "WORKL, {}", &input[..30]);
        cases += 1;
    }
    assert!(cases >= 8);
}
