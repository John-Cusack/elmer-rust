//! The LAPACK 3.12 routines of ARPACK's non-symmetric path
//! (`elmer_eigen::arpack::lapack_ns`) against the system reference
//! `liblapack.so.3` / `libblas.so.3` the pristine Elmer build links, bit for
//! bit.
//!
//! `tests/data/lapack_ns_pristine.txt` holds pairs of lines, recorded once on
//! the super server by calling the system libraries: a routine name with its
//! inputs, then what the library returned. Floats are `f64::to_bits` in hex,
//! vectors and column-major matrices `[h,h,...]`.
//!
//! * `dladiv a b c d` -> `p q`
//! * `dlanv2 a b c d` -> `a b c d rt1r rt1i rt2r rt2i cs sn`
//! * `dlanhs n A` -> `DLANHS('1')`
//! * `dtrmm m n alpha A B` -> `B` (`'Right','Upper','N','Non-unit'`)
//! * `dlaln2 na nw smin ca A d1 d2 B wr wi` -> `X scale xnorm info`
//!   (`LTRANS = .FALSE.`, 2x2 `A`, `B`, `X`; the parts of `X` the routine
//!   does not set are zero)
//! * `dlahqr wantt wantz n H Z` -> `info wr wi H Z` (`ILO = ILOZ = 1`,
//!   `IHI = IHIZ = N`)
//! * `dtrevc n T sel` -> `m VR(:,1:m) sel` (`SIDE = 'R'`; `sel` is `A` for
//!   `HOWMNY = 'A'`, else the `SELECT` flags for 'S', before and after)

use elmer_eigen::arpack::lapack_ns as l;

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

fn hf(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

fn hv(x: &[f64]) -> String {
    format!("[{}]", x.iter().map(|&y| hf(y)).collect::<Vec<_>>().join(","))
}

fn i(t: &str) -> usize {
    t.parse().unwrap()
}

fn eval(line: &str) -> String {
    let t: Vec<&str> = line.split_whitespace().collect();
    match t[0] {
        "dladiv" => {
            let (p, q) = l::dladiv(f(t[1]), f(t[2]), f(t[3]), f(t[4]));
            format!("{} {}", hf(p), hf(q))
        }
        "dlanv2" => {
            let r = l::dlanv2(f(t[1]), f(t[2]), f(t[3]), f(t[4]));
            [r.0, r.1, r.2, r.3, r.4, r.5, r.6, r.7, r.8, r.9].iter().map(|&x| hf(x)).collect::<Vec<_>>().join(" ")
        }
        "dlanhs" => {
            let (n, a) = (i(t[1]), v(t[2]));
            hf(l::dlanhs_one(n, &a, n))
        }
        "dtrmm" => {
            let (m, n, alpha, a, mut b) = (i(t[1]), i(t[2]), f(t[3]), v(t[4]), v(t[5]));
            l::dtrmm_run(m, n, alpha, &a, n, &mut b, m);
            hv(&b)
        }
        "dlaln2" => {
            let b = v(t[8]);
            let (x, scale, xnorm, info) =
                l::dlaln2(false, i(t[1]), i(t[2]), f(t[3]), f(t[4]), &v(t[5]), 2, f(t[6]), f(t[7]), [b[0], b[1], b[2], b[3]], f(t[9]), f(t[10]));
            format!("{} {} {} {}", hv(&x), hf(scale), hf(xnorm), info)
        }
        "dlahqr" => {
            let (wantt, wantz, n, mut h, mut z) = (t[1] == "1", t[2] == "1", i(t[3]), v(t[4]), v(t[5]));
            let (mut wr, mut wi) = (vec![0.0; n], vec![0.0; n]);
            let info = l::dlahqr(wantt, wantz, n, 1, n, &mut h, n, &mut wr, &mut wi, 1, n, &mut z, n);
            format!("{} {} {} {} {}", info, hv(&wr), hv(&wi), hv(&h), hv(&z))
        }
        "dtrevc" => {
            let (n, tm, some) = (i(t[1]), v(t[2]), t[3] != "A");
            let mut sel: Vec<bool> = t[3].chars().map(|c| c == '1').collect();
            let mut vr = vec![0.0; n * n];
            let mut work = vec![0.0; 3 * n];
            let m = l::dtrevc_right(some.then_some(&mut sel[..]), n, &tm, n, &mut vr, n, &mut work);
            let sel: String = if some { sel.iter().map(|&b| if b { '1' } else { '0' }).collect() } else { "A".into() };
            format!("{} {} {}", m, hv(&vr[..n * m]), sel)
        }
        other => panic!("unknown routine {other}"),
    }
}

#[test]
fn nonsymmetric_lapack_routines_match_the_system_library_bit_for_bit() {
    let data = include_str!("data/lapack_ns_pristine.txt");
    let lines: Vec<&str> = data.lines().collect();
    assert_eq!(lines.len() % 2, 0);
    let mut checked = std::collections::BTreeMap::<&str, usize>::new();
    for pair in lines.chunks(2) {
        let got = eval(pair[0]);
        assert_eq!(got, pair[1], "{}", &pair[0][..pair[0].len().min(200)]);
        *checked.entry(pair[0].split_whitespace().next().unwrap()).or_default() += 1;
    }
    for name in ["dladiv", "dlanv2", "dlanhs", "dtrmm", "dlaln2", "dlahqr", "dtrevc"] {
        assert!(checked.get(name).copied().unwrap_or(0) >= 10, "{name}: {checked:?}");
    }
}
