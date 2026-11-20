//! `AMD_1`: build the pattern of `A + A'` in the workspace and call `AMD_2`.
//!
//! AMD Version 1.1, Copyright (c) 2004 by Timothy A. Davis, Patrick R.
//! Amestoy, and Iain S. Duff.  Translated to Rust.

use super::*;

/// `s` has `slen` entries (`slen >= nzaat + nzaat/5 + 7n`, Pinv excluded);
/// `len` holds the column counts from [`amd_aat`] and is overwritten.
#[allow(clippy::too_many_arguments)]
pub fn amd_1(
    n: Int,
    ap: &[Int],
    ai: &[Int],
    p_out: &mut [Int],
    pinv: &mut [Int],
    len: &mut [Int],
    slen: Int,
    s: &mut [Int],
    control: Option<&[f64]>,
    info: Option<&mut [f64]>,
) {
    let nu = n as usize;
    let iwlen = slen - 6 * n;
    let (pe, rest) = s.split_at_mut(nu);
    let (nv, rest) = rest.split_at_mut(nu);
    let (head, rest) = rest.split_at_mut(nu);
    let (elen, rest) = rest.split_at_mut(nu);
    let (degree, rest) = rest.split_at_mut(nu);
    let (w, rest) = rest.split_at_mut(nu);
    let iw = &mut rest[..iwlen as usize];

    // Sp = Nv, Tp = W
    let mut pfree: Int = 0;
    for j in 0..nu {
        pe[j] = pfree;
        nv[j] = pfree;
        pfree += len[j];
    }
    {
        let sp = &mut *nv;
        let tp = &mut *w;
        for k in 0..n {
            let ku = k as usize;
            let p1 = ap[ku];
            let p2 = ap[ku + 1];
            let mut p = p1;
            while p < p2 {
                let j = ai[p as usize];
                if j < k {
                    let ju = j as usize;
                    iw[sp[ju] as usize] = k;
                    sp[ju] += 1;
                    iw[sp[ku] as usize] = j;
                    sp[ku] += 1;
                    p += 1;
                } else if j == k {
                    p += 1;
                    break;
                } else {
                    break;
                }
                let ju = j as usize;
                let pj2 = ap[ju + 1];
                let mut pj = tp[ju];
                while pj < pj2 {
                    let i = ai[pj as usize];
                    if i < k {
                        let iu = i as usize;
                        iw[sp[iu] as usize] = j;
                        sp[iu] += 1;
                        iw[sp[ju] as usize] = i;
                        sp[ju] += 1;
                        pj += 1;
                    } else if i == k {
                        pj += 1;
                        break;
                    } else {
                        break;
                    }
                }
                tp[ju] = pj;
            }
            tp[ku] = p;
        }
        for j in 0..nu {
            for pj in tp[j]..ap[j + 1] {
                let i = ai[pj as usize];
                let iu = i as usize;
                iw[sp[iu] as usize] = j as Int;
                sp[iu] += 1;
                iw[sp[j] as usize] = i;
                sp[j] += 1;
            }
        }
    }
    amd_2(
        n, pe, iw, len, iwlen, pfree, nv, pinv, p_out, head, elen, degree, w, control, info,
    );
}
