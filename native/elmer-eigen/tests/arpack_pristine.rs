//! The ARPACK port against the pristine build's `libarpack.so` (ARPACK 2.4
//! from Elmer's `mathlibs/src/arpack`, with the reference BLAS/LAPACK 3.12
//! the pristine Elmer links), bit for bit.
//!
//! `tests/data/arpack_pristine.txt` holds pairs of lines. The first is a
//! generalized symmetric problem `K x = lambda M x` in shift-invert mode 3,
//! as Elmer's `ArpackEigenSolve` drives it: `n nev ncv which sigma tol
//! mxiter OP B`, with `OP = (K - sigma M)^-1` and `B = M` dense and
//! row-major, floats as `f64::to_bits` in hex. The second is what
//! `dsaupd_`/`dseupd_` of the pristine library returned when driven through
//! reverse communication with the same `OP` and `B` products (one
//! accumulator, columns in order), recorded once on the super server:
//! `info IPARAM(3) IPARAM(5) IPARAM(9) IPARAM(10) IPARAM(11) ierr d v`, the
//! Ritz values and the Ritz vectors (column-major) in hex. Each problem
//! starts from a fresh `dgetv0` seed (1, 3, 5, 7), as one process does.

use elmer_eigen::arpack::saupd::{dsaupd, Arpack, Which};
use elmer_eigen::arpack::seupd::dseupd;

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
    assert_eq!(t[0], "arpack");
    let (n, nev, ncv): (usize, usize, usize) = (t[1].parse().unwrap(), t[2].parse().unwrap(), t[3].parse().unwrap());
    let which = match t[4] {
        "LM" => Which::LM,
        "SM" => Which::SM,
        "LA" => Which::LA,
        "SA" => Which::SA,
        _ => Which::BE,
    };
    let (sigma, mut tol, mxiter) = (f(t[5]), f(t[6]), t[7].parse::<i32>().unwrap());
    let (opm, b) = (v(t[8]), v(t[9]));
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
        match ido {
            -1 => {
                let mut bx = vec![0.0; n];
                mv(&b, &wd[ip[0]..ip[0] + n], &mut bx);
                mv(&opm, &bx, &mut y);
            }
            1 => mv(&opm, &wd[ip[2]..ip[2] + n], &mut y),
            2 => mv(&b, &wd[ip[0]..ip[0] + n], &mut y),
            _ => return Err(format!("ido {ido}")),
        }
        wd[ip[1]..ip[1] + n].copy_from_slice(&y);
        Ok(())
    };
    let mut arp = Arpack::default();
    let mut resid = vec![0.0; n];
    let mut vv = vec![0.0; n * ncv];
    let mut workd = vec![0.0; 3 * n];
    let mut workl = vec![0.0; ncv * ncv + 8 * ncv];
    let mut iparam = [0i32; 11];
    iparam[0] = 1;
    iparam[2] = mxiter;
    iparam[6] = 3;
    let au = dsaupd(&mut op, &mut arp, true, n, which, nev, &mut tol, &mut resid, ncv, &mut vv, n, iparam, &mut workd, &mut workl, 0)
        .unwrap();
    // The recording calls dseupd only after a clean dsaupd (99 otherwise).
    let (d, ierr) = if au.info == 0 {
        dseupd(sigma, true, n, which, nev, tol, &resid, ncv, &mut vv, n, &au, &mut workd, &mut workl)
    } else {
        (vec![], 99)
    };
    let nconv = au.iparam[4].max(0) as usize;
    format!(
        "{} {} {} {} {} {} {} {} {}",
        au.info,
        au.iparam[2],
        au.iparam[4],
        au.iparam[8],
        au.iparam[9],
        au.iparam[10],
        ierr,
        hv(&d[..nconv]),
        hv(&vv[..n * nconv])
    )
}

#[test]
fn dsaupd_and_dseupd_match_pristine_arpack_bitwise() {
    let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/arpack_pristine.txt")).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len() % 2, 0);
    for pair in lines.chunks(2) {
        let got = solve(pair[0]);
        let head: Vec<&str> = pair[0].split_whitespace().take(8).collect();
        assert_eq!(got, pair[1], "problem {head:?}");
    }
}
