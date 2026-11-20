//! The complex BLAS and LAPACK 3.12 ports under ARPACK's complex driver
//! (`zblas`, `zlapack`) against the system `libblas.so.3` and
//! `liblapack.so.3` the pristine Elmer build links, bit for bit.
//!
//! Each line of `tests/data/zlapack_lapack.txt` is `op args<TAB>result`,
//! recorded from the library through ctypes on random inputs (complex
//! entries as pairs of f64 bit patterns): `ZGEMV`, `ZGERU`, `ZTRMM`,
//! `ZDOTC`, `DZNRM2`, `IZAMAX`, `DZASUM`, `ZAXPY`, `ZSCAL`, `ZDSCAL`, `ZLAHQR`
//! (with and without the Schur form and vectors, and on a sub-block),
//! `ZLARTG`, `ZLARFG`, `ZLARNV`, `ZLANHS`, `ZLASCL`, `ZLADIV`, gfortran's
//! complex `SQRT`, `ZTREVC('R','B')`, `ZLATRS` (both its paths), `ZGEQR2`,
//! `ZUNM2R`, `ZLARF` and `ZTRSEN('N','V')`.
//!
//! One known difference: on a sub-block (`ILO > 1`) the library's `ZLAHQR`
//! leaves the sign of some zero entries of `Z` outside the block different;
//! ARPACK only calls it on the whole matrix, so those cases compare `Z` by
//! value.

use elmer_eigen::arpack::zblas::*;
use elmer_eigen::arpack::zlapack::*;
use num_complex::Complex64 as C;

fn f(t: &str) -> f64 { f64::from_bits(u64::from_str_radix(t, 16).unwrap()) }
fn zv(t: &str) -> Vec<C> {
    if t == "-" { return vec![]; }
    let v: Vec<f64> = t.split(',').map(f).collect();
    v.chunks(2).map(|p| C::new(p[0], p[1])).collect()
}
fn hz(v: &[C]) -> String { if v.is_empty() { return "-".into(); } v.iter().map(|z| format!("{:016x},{:016x}", z.re.to_bits(), z.im.to_bits())).collect::<Vec<_>>().join(",") }
fn hd(x: f64) -> String { format!("{:016x}", x.to_bits()) }


fn run(line: &str) -> String {
    let t: Vec<&str> = line.split_whitespace().collect();
    let u = |k: usize| t[k].parse::<usize>().unwrap();
    match t[0] {
            "zgemv_n" | "zgemv_c" => {
                let (m, n) = (u(1), u(2));
                let (al, be) = (zv(t[3])[0], zv(t[4])[0]);
                let (a, x, mut y) = (zv(t[5]), zv(t[6]), zv(t[7]));
                if t[0] == "zgemv_n" { zgemv_n(m, n, al, &a, m, &x, be, &mut y) } else { zgemv_c(m, n, al, &a, m, &x, be, &mut y) }
                hz(&y)
            }
            "zgeru" => {
                let (m, n) = (u(1), u(2));
                let al = zv(t[3])[0];
                let (x, y, mut a) = (zv(t[4]), zv(t[5]), zv(t[6]));
                zgeru(m, n, al, &x, &y, &mut a, m);
                hz(&a)
            }
            "ztrmm" => {
                let (m, n) = (u(1), u(2));
                let al = zv(t[3])[0];
                let (a, mut b) = (zv(t[4]), zv(t[5]));
                ztrmm_runn(m, n, al, &a, n, &mut b, m);
                hz(&b)
            }
            "zdotc" => { let (x, y) = (zv(t[2]), zv(t[3])); hz(&[zdotc(u(1), &x, &y)]) }
            "dznrm2" => hd(dznrm2(u(1), &zv(t[2]))),
            "izamax" => izamax(u(1), &zv(t[2])).to_string(),
            "dzasum" => hd(dzasum(u(1), &zv(t[2]))),
            "zaxpy" => { let al = zv(t[2])[0]; let (x, mut y) = (zv(t[3]), zv(t[4])); zaxpy(u(1), al, &x, &mut y); hz(&y) }
            "zdscal" => { let mut x = zv(t[3]); zdscal(u(1), f(t[2]), &mut x); hz(&x) }
            "zscal" => { let al = zv(t[2])[0]; let mut x = zv(t[3]); zscal(u(1), al, &mut x); hz(&x) }
            "cdiv" => { let (a, b) = (zv(t[1])[0], zv(t[2])[0]); hz(&[cdiv(a, b)]) }
            "cabs" => hd(cabs(zv(t[1])[0])),
            "csqrt" => hz(&[csqrt(zv(t[1])[0])]),
            "zladiv" => hz(&[zladiv(zv(t[1])[0], zv(t[2])[0])]),
            "zlartg" => { let (c, s, r) = zlartg(zv(t[1])[0], zv(t[2])[0]); format!("{} {} {}", hd(c), hz(&[s]), hz(&[r])) }
            "zlarfg" => { let n = u(1); let mut a = zv(t[2])[0]; let mut x = zv(t[3]); let tau = zlarfg(n, &mut a, &mut x); format!("{} {} {}", hz(&[a]), hz(&x), hz(&[tau])) }
            "zlarnv" => { let mut s4 = [0i32; 4]; for k in 0..4 { s4[k] = t[2 + k].parse().unwrap(); } let mut x = vec![C::new(0.0, 0.0); u(1)]; zlarnv2(&mut s4, u(1), &mut x); format!("{} {:?}", hz(&x), s4) }
            "zlanhs" => hd(zlanhs_one(u(1), &zv(t[2]), u(1))),
            "zlascl" => { let (m, n) = (u(1), u(2)); let mut a = zv(t[5]); zlascl_g(f(t[3]), f(t[4]), m, n, &mut a, m); hz(&a) }
            "zlahqr" => {
                let (wantt, wantz, n, ilo, ihi) = (t[1] == "1", t[2] == "1", u(3), u(4), u(5));
                let (mut h, mut z) = (zv(t[6]), zv(t[7]));
                let mut w = vec![C::new(0.0, 0.0); n];
                let info = zlahqr(wantt, wantz, n, ilo, ihi, &mut h, n, &mut w, 1, n, &mut z, n);
                format!("{} {} {} {}", info, hz(&h), hz(&w), hz(&z))
            }
            "ztrevc" => { let n = u(1); let (mut tt, mut vr) = (zv(t[2]), zv(t[3])); ztrevc_rb(n, &mut tt, n, &mut vr, n); format!("{} {}", hz(&tt), hz(&vr)) }
            "zlatrs" => { let n = u(1); let a = zv(t[2]); let mut x = zv(t[3]); let mut cn: Vec<f64> = t[4].split(',').map(f).collect(); let sc = zlatrs_unn(n, &a, n, &mut x, &mut cn); format!("{} {} {}", hz(&x), hd(sc), cn.iter().map(|v| hd(*v)).collect::<Vec<_>>().join(",")) }
            "zgeqr2" => { let (m, n) = (u(1), u(2)); let mut a = zv(t[3]); let tau = zgeqr2(m, n, &mut a, m); format!("{} {}", hz(&a), hz(&tau)) }
            "zunm2r" => { let (left, conj, m, n, k) = (t[1] == "L", t[2] == "C", u(3), u(4), u(5)); let nq = if left { m } else { n }; let mut a = zv(t[6]); let tau = zv(t[7]); let mut c = zv(t[8]); zunm2r(left, conj, m, n, k, &mut a, nq, &tau, &mut c, m); hz(&c) }
            "zlarf" => { let (left, m, n) = (t[1] == "L", u(2), u(3)); let v = zv(t[4]); let tau = zv(t[5])[0]; let mut c = zv(t[6]); zlarf(left, m, n, &v, tau, &mut c, m); hz(&c) }
            "ztrsen" => { let n = u(1); let sel: Vec<bool> = t[2].chars().map(|c| c == '1').collect(); let (mut tt, mut q) = (zv(t[3]), zv(t[4])); let mut w = vec![C::new(0.0, 0.0); n]; let m = ztrsen_nv(&sel, n, &mut tt, n, &mut q, n, &mut w); format!("{} {} {} {}", m, hz(&tt), hz(&q), hz(&w)) }
            op => panic!("op {op}"),
        }
}

fn same_value(a: &str, b: &str) -> bool {
    let (va, vb): (Vec<&str>, Vec<&str>) = (a.split(',').collect(), b.split(',').collect());
    va.len() == vb.len() && va.iter().zip(&vb).all(|(x, y)| f(x) == f(y))
}

#[test]
fn complex_blas_and_lapack_match_the_system_libraries_bit_for_bit() {
    let data = include_str!("data/zlapack_lapack.txt");
    let mut n = 0;
    for line in data.lines() {
        let (input, want) = line.split_once('\t').unwrap();
        let got = run(input);
        let t: Vec<&str> = input.split_whitespace().collect();
        if t[0] == "zlahqr" && t[4] != "1" {
            // Sub-block: INFO, H and W bit for bit, Z by value.
            let (w, g): (Vec<&str>, Vec<&str>) = (want.split(' ').collect(), got.split(' ').collect());
            assert_eq!(&w[..3], &g[..3], "{}", &input[..40]);
            assert!(same_value(w[3], g[3]), "{}", &input[..40]);
        } else {
            assert_eq!(got, want, "{}", &input[..input.len().min(60)]);
        }
        n += 1;
    }
    assert!(n >= 150);
}
