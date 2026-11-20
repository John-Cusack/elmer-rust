//! AMD 1.1 approximate minimum degree ordering (`amd_*`, `int` version).
//!
//! AMD Version 1.1 (Jan. 21, 2004), Copyright (c) 2004 by Timothy A. Davis,
//! Patrick R. Amestoy, and Iain S. Duff.  See the crate documentation for
//! the licence.  Translated to Rust; see the modification notice there.

mod aat;
mod amd1;
mod amd2;
mod postorder;

pub use aat::amd_aat;
pub use amd1::amd_1;
pub use amd2::amd_2;
pub use postorder::{amd_post_tree, amd_postorder};

/// The C `Int` of the `di` build.
pub type Int = i32;

pub const EMPTY: Int = -1;

#[inline]
pub(crate) const fn flip(i: Int) -> Int {
    -i - 2
}

pub const AMD_CONTROL: usize = 5;
pub const AMD_INFO: usize = 20;

pub const AMD_DENSE: usize = 0;
pub const AMD_AGGRESSIVE: usize = 1;

pub const AMD_DEFAULT_DENSE: f64 = 10.0;
pub const AMD_DEFAULT_AGGRESSIVE: f64 = 1.0;

pub const AMD_STATUS: usize = 0;
pub const AMD_N: usize = 1;
pub const AMD_NZ: usize = 2;
pub const AMD_SYMMETRY: usize = 3;
pub const AMD_NZDIAG: usize = 4;
pub const AMD_NZ_A_PLUS_AT: usize = 5;
pub const AMD_NDENSE: usize = 6;
pub const AMD_MEMORY: usize = 7;
pub const AMD_NCMPA: usize = 8;
pub const AMD_LNZ: usize = 9;
pub const AMD_NDIV: usize = 10;
pub const AMD_NMULTSUBS_LDL: usize = 11;
pub const AMD_NMULTSUBS_LU: usize = 12;
pub const AMD_DMAX: usize = 13;

pub const AMD_OK: Int = 0;
pub const AMD_OUT_OF_MEMORY: Int = -1;
pub const AMD_INVALID: Int = -2;

/// `AMD_defaults`.
pub fn amd_defaults(control: &mut [f64]) {
    for c in control.iter_mut().take(AMD_CONTROL) {
        *c = 0.0;
    }
    control[AMD_DENSE] = AMD_DEFAULT_DENSE;
    control[AMD_AGGRESSIVE] = AMD_DEFAULT_AGGRESSIVE;
}

/// `AMD_valid`: columns sorted, no duplicates, row indices in range.
pub fn amd_valid(n_row: Int, n_col: Int, ap: &[Int], ai: &[Int]) -> bool {
    if n_row < 0 || n_col < 0 {
        return false;
    }
    let nz = ap[n_col as usize];
    if ap[0] != 0 || nz < 0 {
        return false;
    }
    for j in 0..n_col as usize {
        let p1 = ap[j];
        let p2 = ap[j + 1];
        if p1 > p2 {
            return false;
        }
        let mut ilast = EMPTY;
        for p in p1..p2 {
            let i = ai[p as usize];
            if i <= ilast || i >= n_row {
                return false;
            }
            ilast = i;
        }
    }
    true
}

/// `AMD_preprocess_valid`.
pub fn amd_preprocess_valid(n: Int, ap: &[Int], ai: &[Int]) -> bool {
    if n < 0 {
        return false;
    }
    let nz = ap[n as usize];
    if ap[0] != 0 || nz < 0 {
        return false;
    }
    for j in 0..n as usize {
        if ap[j] > ap[j + 1] {
            return false;
        }
    }
    for p in 0..nz as usize {
        let i = ai[p];
        if i < 0 || i >= n {
            return false;
        }
    }
    true
}

/// `AMD_wpreprocess`: `R = pattern of A'`, duplicates removed, sorted.
pub fn amd_wpreprocess(
    n: Int,
    ap: &[Int],
    ai: &[Int],
    rp: &mut [Int],
    ri: &mut [Int],
    w: &mut [Int],
    flag: &mut [Int],
) {
    let n = n as usize;
    for i in 0..n {
        w[i] = 0;
        flag[i] = EMPTY;
    }
    for j in 0..n {
        let p2 = ap[j + 1];
        for p in ap[j]..p2 {
            let i = ai[p as usize] as usize;
            if flag[i] != j as Int {
                w[i] += 1;
                flag[i] = j as Int;
            }
        }
    }
    rp[0] = 0;
    for i in 0..n {
        rp[i + 1] = rp[i] + w[i];
    }
    for i in 0..n {
        w[i] = rp[i];
        flag[i] = EMPTY;
    }
    for j in 0..n {
        let p2 = ap[j + 1];
        for p in ap[j]..p2 {
            let i = ai[p as usize] as usize;
            if flag[i] != j as Int {
                ri[w[i] as usize] = j as Int;
                w[i] += 1;
                flag[i] = j as Int;
            }
        }
    }
}

/// `AMD_preprocess`.
pub fn amd_preprocess(n: Int, ap: &[Int], ai: &[Int], rp: &mut [Int], ri: &mut [Int]) -> Int {
    if !amd_preprocess_valid(n, ap, ai) {
        return AMD_INVALID;
    }
    let m = n.max(1) as usize;
    let mut w = vec![0; m];
    let mut flag = vec![0; m];
    amd_wpreprocess(n, ap, ai, rp, ri, &mut w, &mut flag);
    AMD_OK
}

/// `AMD_order`: the user-callable ordering. `p` receives the permutation;
/// `info`, when given, must hold `AMD_INFO` entries.
pub fn amd_order(
    n: Int,
    ap: &[Int],
    ai: &[Int],
    p: &mut [Int],
    control: Option<&[f64]>,
    mut info: Option<&mut [f64]>,
) -> Int {
    if let Some(info) = info.as_deref_mut() {
        for v in info.iter_mut().take(AMD_INFO) {
            *v = EMPTY as f64;
        }
        info[AMD_N] = n as f64;
        info[AMD_STATUS] = AMD_OK as f64;
    }
    if n < 0 {
        if let Some(info) = info.as_deref_mut() {
            info[AMD_STATUS] = AMD_INVALID as f64;
        }
        return AMD_INVALID;
    }
    if n == 0 {
        return AMD_OK;
    }
    let nz = ap[n as usize];
    if let Some(info) = info.as_deref_mut() {
        info[AMD_NZ] = nz as f64;
    }
    if nz < 0 {
        if let Some(info) = info.as_deref_mut() {
            info[AMD_STATUS] = AMD_INVALID as f64;
        }
        return AMD_INVALID;
    }
    if (2.4 * nz as f64 + 8.0 * n as f64) > Int::MAX as f64 / core::mem::size_of::<Int>() as f64 {
        if let Some(info) = info.as_deref_mut() {
            info[AMD_STATUS] = AMD_OUT_OF_MEMORY as f64;
        }
        return AMD_OUT_OF_MEMORY;
    }
    if !amd_valid(n, n, ap, ai) {
        if let Some(info) = info.as_deref_mut() {
            info[AMD_STATUS] = AMD_INVALID as f64;
        }
        return AMD_INVALID;
    }
    let nu = n as usize;
    let mut len = vec![0 as Int; nu];
    let nzaat = amd_aat(n, ap, ai, &mut len, p, info.as_deref_mut());
    let mut slen = (nzaat + nzaat / 5 + n) + 7 * n;
    if let Some(info) = info.as_deref_mut() {
        info[AMD_MEMORY] = (slen as f64 + n as f64) * core::mem::size_of::<Int>() as f64;
    }
    let mut s = vec![0 as Int; slen as usize];
    slen -= n;
    let (s_work, pinv) = s.split_at_mut(slen as usize);
    amd_1(n, ap, ai, p, pinv, &mut len, slen, s_work, control, info);
    AMD_OK
}
