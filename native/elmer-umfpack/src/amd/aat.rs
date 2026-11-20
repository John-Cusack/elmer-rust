//! `AMD_aat`: column counts of `A + A'` (diagonal excluded) and the pattern
//! symmetry statistics.
//!
//! AMD Version 1.1, Copyright (c) 2004 by Timothy A. Davis, Patrick R.
//! Amestoy, and Iain S. Duff.  Translated to Rust.

use super::*;

/// Returns `nzaat`; `len[k]` gets the number of off-diagonal entries in column
/// `k` of `A + A'`; `tp` is workspace of size `n`.
pub fn amd_aat(
    n: Int,
    ap: &[Int],
    ai: &[Int],
    len: &mut [Int],
    tp: &mut [Int],
    info: Option<&mut [f64]>,
) -> Int {
    let mut info = info;
    if let Some(info) = info.as_deref_mut() {
        for v in info.iter_mut().take(AMD_INFO) {
            *v = EMPTY as f64;
        }
        info[AMD_STATUS] = AMD_OK as f64;
    }
    let nu = n as usize;
    for k in 0..nu {
        len[k] = 0;
    }
    let mut nzdiag: Int = 0;
    let mut nzboth: Int = 0;
    let nz = ap[nu];
    for k in 0..n {
        let ku = k as usize;
        let p1 = ap[ku];
        let p2 = ap[ku + 1];
        let mut p = p1;
        while p < p2 {
            let j = ai[p as usize];
            if j < k {
                len[j as usize] += 1;
                len[ku] += 1;
                p += 1;
            } else if j == k {
                p += 1;
                nzdiag += 1;
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
                    len[i as usize] += 1;
                    len[ju] += 1;
                    pj += 1;
                } else if i == k {
                    pj += 1;
                    nzboth += 1;
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
            let i = ai[pj as usize] as usize;
            len[i] += 1;
            len[j] += 1;
        }
    }
    let sym: f64 = if nz == nzdiag {
        1.0
    } else {
        (2 * nzboth) as f64 / (nz - nzdiag) as f64
    };
    let mut nzaat: Int = 0;
    for k in 0..nu {
        nzaat += len[k];
    }
    if let Some(info) = info {
        info[AMD_STATUS] = AMD_OK as f64;
        info[AMD_N] = n as f64;
        info[AMD_NZ] = nz as f64;
        info[AMD_SYMMETRY] = sym;
        info[AMD_NZDIAG] = nzdiag as f64;
        info[AMD_NZ_A_PLUS_AT] = nzaat as f64;
    }
    nzaat
}
