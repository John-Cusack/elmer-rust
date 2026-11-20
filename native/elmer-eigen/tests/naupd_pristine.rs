//! The port of ARPACK's `DNAUPD` against the pristine build's `libarpack.so`
//! (ARPACK 2.4 from Elmer's `mathlibs/src/arpack`, with the reference
//! BLAS/LAPACK 3.12 the pristine Elmer links), bit for bit.
//!
//! `tests/data/naupd_pristine.txt` holds pairs of lines. The first is a
//! problem: `n nev ncv which tol mxiter mode OP B`, dense row-major, floats
//! as `f64::to_bits` in hex. In mode 1 (`bmat = 'I'`) `OP` is applied
//! directly; in mode 3 (`bmat = 'G'`, shift-invert) `ido = -1` gives
//! `OP (B x)`, `ido = 1` `OP` of the `B x` at `ipntr(3)`, `ido = 2` `B x`.
//! The second is what `dnaupd_` of the pristine library returned when
//! driven through reverse communication with the same products (one
//! accumulator, columns in order), recorded once on the super server:
//! `info IPARAM(3) IPARAM(5) IPARAM(9) IPARAM(10) IPARAM(11)`, then the
//! real and imaginary Ritz values, the Ritz estimates and `H` from `workl`,
//! `V` (column-major) and `RESID`. Each problem starts from a fresh process
//! state (`dgetv0`'s seed 1, 3, 5, 7 and the unset `SAVE`d `smlnum`s).

use elmer_eigen::arpack::naupd::{dnaupd, NWhich};
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
    let (opm, b) = (v(t[8]), v(t[9]));
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
    let w = |a: usize, len: usize| hv(&workl[a - 1..a - 1 + len]);
    format!(
        "{} {} {} {} {} {} {} {} {} {} {} {}",
        au.info,
        au.iparam[2],
        au.iparam[4],
        au.iparam[8],
        au.iparam[9],
        au.iparam[10],
        w(au.ipntr[5], ncv),
        w(au.ipntr[6], ncv),
        w(au.ipntr[7], ncv),
        w(au.ipntr[4], ncv * ncv),
        hv(&vv),
        hv(&resid)
    )
}

#[test]
fn dnaupd_matches_the_pristine_libarpack_bit_for_bit() {
    let data = include_str!("data/naupd_pristine.txt");
    let lines: Vec<&str> = data.lines().collect();
    assert!(lines.len() >= 12 && lines.len().is_multiple_of(2));
    for pair in lines.chunks(2) {
        assert_eq!(solve(pair[0]), pair[1], "problem {}", &pair[0][..60]);
    }
}
