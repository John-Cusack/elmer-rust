//! The port of ARPACK's `DNEUPD` (and the LAPACK `DTRSEN` it reorders the
//! Schur form with) against the pristine build's `libarpack.so` and the
//! system reference LAPACK 3.12 the pristine Elmer links, bit for bit.
//!
//! `tests/data/neupd_pristine.txt` holds pairs of lines. The first is the
//! problem of `tests/naupd_pristine.rs` followed by `sigmar` (0 in mode
//! 1). The second is what the pristine library returned for `dnaupd_`
//! followed by `dneupd_(.true., 'A', ...)` with `Z` the same array as `V`
//! and `sigmai = 0`, as Elmer's `ArpackEigenSolve` calls it: `info ierr NC
//! dr di Z(:,1:nconv)` (`info 99` when `dnaupd_` did not converge). Four of
//! the ten problems go through `dneupd`'s `DTRSEN` reordering.
//!
//! `tests/data/trsen_pristine.txt` holds `DTRSEN('N', 'V')` cases: `dtrsen n
//! T Q select` on a Schur form from the library's `DLAHQR`, then `M INFO
//! WR WI T Q` as the library returned them.

use elmer_eigen::arpack::lapack_ns::dtrsen_nv;
use elmer_eigen::arpack::naupd::{dnaupd, NWhich};
use elmer_eigen::arpack::neupd::dneupd;
use elmer_eigen::arpack::saupd::Arpack;

fn f(t: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(t, 16).unwrap())
}

fn v(t: &str) -> Vec<f64> {
    let t = t.trim_start_matches('[').trim_end_matches(']');
    if t.is_empty() {
        return vec![];
    }
    t.split(',').map(f).collect()
}

fn hv(x: &[f64]) -> String {
    format!("[{}]", x.iter().map(|y| format!("{:016x}", y.to_bits())).collect::<Vec<_>>().join(","))
}

fn solve(line: &str) -> String {
    let t: Vec<&str> = line.split_whitespace().collect();
    assert_eq!(t[0], "dnaupd");
    let (n, nev, ncv): (usize, usize, usize) = (t[1].parse().unwrap(), t[2].parse().unwrap(), t[3].parse().unwrap());
    let which = match t[4] {
        "LM" => NWhich::LM,
        "SM" => NWhich::SM,
        "LR" => NWhich::LR,
        "SR" => NWhich::SR,
        "LI" => NWhich::LI,
        _ => NWhich::SI,
    };
    let (mut tol, mxiter, mode) = (f(t[5]), t[6].parse::<i32>().unwrap(), t[7].parse::<i32>().unwrap());
    let (opm, b, sigmar) = (v(t[8]), v(t[9]), f(t[10]));
    let bmat_g = mode != 1;
    let mv = |m: &[f64], x: &[f64], y: &mut [f64]| {
        for r in 0..n {
            let mut s = 0.0;
            for c in 0..n {
                s += m[r * n + c] * x[c];
            }
            y[r] = s;
        }
    };
    let mut op = |ido: i32, ip: [usize; 3], wd: &mut [f64]| -> Result<(), String> {
        let mut y = vec![0.0; n];
        match (ido, bmat_g) {
            (-1, true) => {
                let mut bx = vec![0.0; n];
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
    let mut arp = Arpack::default();
    let mut resid = vec![0.0; n];
    let mut vv = vec![0.0; n * ncv];
    let mut workd = vec![0.0; 3 * n];
    let mut workl = vec![0.0; 3 * ncv * ncv + 6 * ncv];
    let mut iparam = [0i32; 11];
    iparam[0] = 1;
    iparam[2] = mxiter;
    iparam[6] = mode;
    let au = dnaupd(&mut op, &mut arp, bmat_g, n, which, nev, &mut tol, &mut resid, ncv, &mut vv, n, iparam, &mut workd, &mut workl, 0)
        .unwrap();
    if au.info != 0 {
        return format!("{} 99", au.info);
    }
    let r = dneupd(sigmar, 0.0, bmat_g, n, which, nev, tol, &resid, ncv, &mut vv, n, &au, &mut workd, &mut workl);
    let nc = r.nconv;
    assert_eq!(nc, au.iparam[4] as usize, "DTRSEN's M is IPARAM(5) here");
    format!("{} {} NC {} {} {}", au.info, r.info, hv(&r.dr[..nc]), hv(&r.di[..nc]), hv(&vv[..n * nc]))
}

#[test]
fn dneupd_matches_the_pristine_libarpack_bit_for_bit() {
    let data = include_str!("data/neupd_pristine.txt");
    let lines: Vec<&str> = data.lines().collect();
    assert!(lines.len() >= 16 && lines.len().is_multiple_of(2));
    for pair in lines.chunks(2) {
        assert_eq!(solve(pair[0]), pair[1], "problem {}", &pair[0][..60]);
    }
}

#[test]
fn dtrsen_matches_the_system_lapack_bit_for_bit() {
    let data = include_str!("data/trsen_pristine.txt");
    let lines: Vec<&str> = data.lines().collect();
    assert!(lines.len() >= 20 && lines.len().is_multiple_of(2));
    for pair in lines.chunks(2) {
        let t: Vec<&str> = pair[0].split_whitespace().collect();
        assert_eq!(t[0], "dtrsen");
        let n: usize = t[1].parse().unwrap();
        let (mut tm, mut q) = (v(t[2]), v(t[3]));
        let sel: Vec<bool> = t[4].chars().map(|c| c == '1').collect();
        let (mut wr, mut wi) = (vec![0.0; n], vec![0.0; n]);
        let (m, info) = dtrsen_nv(&sel, n, &mut tm, n, &mut q, n, &mut wr, &mut wi);
        let got = format!("{} {} {} {} {} {}", m, info, hv(&wr), hv(&wi), hv(&tm), hv(&q));
        assert_eq!(got, pair[1], "{}", &pair[0][..60]);
    }
}
