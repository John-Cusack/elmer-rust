//! The arithmetic of gfortran's `MATMUL` as the pristine Elmer build runs
//! it, where a sum that Elmer writes as `MATMUL` feeds a result the port
//! must reproduce to round-off.
//!
//! The pristine build is a Debug (`-g`, unoptimized) one, so `MATMUL` is
//! not inlined: it calls libgfortran's `_gfortran_matmul_r8`, which
//! dispatches by CPU. On the AMD Zen machine the golden outputs come from,
//! that is the `avx128_fma3` variant, and for a result with at most three
//! rows (every remainder row of its 4x4-blocked kernel) each entry is the
//! compiler-vectorized reduction of its `f11 += a * b` loop:
//!
//! * a main loop of 16 terms per iteration into eight 2-lane accumulators
//!   (fused multiply-add);
//! * then fresh accumulators for one chunk each of 8, 4 and 2 terms, as
//!   many as fit;
//! * the accumulators folded lane by lane in that order, the two lanes
//!   added;
//! * the terms left over, one by one with a fused multiply-add onto the sum.
//!
//! Checked bit for bit against the library on 400 random `2 x k` by `k x 2`
//! products, `k` from 1 to 27.

/// `SUM over l of a(l) * b(l)` as one entry of libgfortran's `MATMUL`
/// computes it (see the module notes); `term(l)` gives `(a(l), b(l))`.
pub fn matmul_dot(k: usize, term: impl Fn(usize) -> (f64, f64)) -> f64 {
    const MAIN: usize = 16;
    let mut lanes: Vec<[f64; 2]> = Vec::with_capacity(MAIN / 2 + 7);
    let mut l = 0;
    if k >= MAIN {
        let mut main = [[0.0f64; 2]; MAIN / 2];
        while l + MAIN <= k {
            for x in 0..MAIN {
                let (a, b) = term(l + x);
                let v = &mut main[x / 2][x % 2];
                *v = a.mul_add(b, *v);
            }
            l += MAIN;
        }
        lanes.extend_from_slice(&main);
    }
    let mut size = MAIN / 2;
    while size >= 2 {
        if l + size <= k {
            let mut chunk = vec![[0.0f64; 2]; size / 2];
            for x in 0..size {
                let (a, b) = term(l + x);
                let v = &mut chunk[x / 2][x % 2];
                *v = a.mul_add(b, *v);
            }
            lanes.extend_from_slice(&chunk);
            l += size;
        }
        size /= 2;
    }
    let mut s = match lanes.split_first() {
        Some((first, rest)) => {
            let mut acc = *first;
            for v in rest {
                acc[0] += v[0];
                acc[1] += v[1];
            }
            acc[0] + acc[1]
        }
        None => 0.0,
    };
    for x in l..k {
        let (a, b) = term(x);
        s = a.mul_add(b, s);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::matmul_dot;

    /// Values recorded from `_gfortran_matmul_r8` on the pristine build's
    /// machine: the displacement gradient of the first integration point of
    /// the buckling case's first element (2 x 9 times 9 x 2), which a plain
    /// left-to-right sum misses by up to 7 ulps.
    #[test]
    fn reproduces_the_library_on_a_recorded_gradient() {
        let u: [[f64; 2]; 9] = [
            [0.0, 0.0],
            [-1.1422239198041311e-05, -4.8261681219575974e-06],
            [-1.082559979129621e-05, -1.6487887814422957e-06],
            [0.0, 0.0],
            [-6.2146081095587278e-06, -4.007919469577654e-06],
            [-1.1180472665968105e-05, -3.3597397592789929e-06],
            [-5.0482163142985282e-06, -1.1160283648824584e-06],
            [0.0, 0.0],
            [-5.2692317302622973e-06, -2.3830604425882542e-06],
        ];
        let db: [[f64; 9]; 2] = [
            [
                -1.9177465533747342,
                -8.9016133230340593,
                70.082253446625316,
                15.098386676965937,
                10.819359876408804,
                40.787093415727433,
                -85.180640123591232,
                8.7870934157274743,
                -49.574186831454959,
            ],
            [
                -2.8766198300456525,
                22.647580015319406,
                105.12338016933688,
                -13.352419984474741,
                13.180640123515843,
                -127.77096018465623,
                61.180640123241304,
                16.229039814520409,
                -74.361280246757246,
            ],
        ];
        let want = [
            [-0.00048903496095962259, 3.2891208024753985e-05],
            [-3.9784999190351014e-05, 0.00020275104789825598],
        ];
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(matmul_dot(9, |k| (u[k][i], db[j][k])), want[i][j], "Grad({},{})", i + 1, j + 1);
            }
        }
    }

    /// Short sums: one term is one rounded product; two terms one 2-lane
    /// vector, its lanes added.
    #[test]
    fn short_sums() {
        assert_eq!(matmul_dot(0, |_| (1.0, 1.0)), 0.0);
        assert_eq!(matmul_dot(1, |_| (0.1, 0.3)), 0.1 * 0.3);
        let t = [(0.1, 0.7), (0.2, 0.9)];
        assert_eq!(matmul_dot(2, |l| t[l]), 0.1 * 0.7 + 0.2 * 0.9);
    }
}
