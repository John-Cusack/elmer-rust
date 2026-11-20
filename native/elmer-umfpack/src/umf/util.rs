//! Small UMFPACK helpers: the `umf_internal.h` size macros, `GET_CONTROL`,
//! `UMF_transpose`, `UMF_fsize`, `UMF_apply_order`, `UMF_is_permutation`,
//! `UMF_symbolic_usage`, `UMF_set_stats`.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::symbolic::Symbolic;
use crate::amd::{Int, EMPTY};
use crate::ptr::P;

// sizeof of the C types of the x86_64 `di` build (checked with the pristine
// headers): Unit and Align are 8 bytes, Element 7 Ints, Tuple 2 Ints.
pub const SIZEOF_INT: usize = 4;
pub const SIZEOF_UNIT: usize = 8;
pub const SIZEOF_ENTRY: usize = 8;
pub const SIZEOF_PTR: usize = 8;
pub const SIZEOF_ELEMENT: usize = 28;
pub const SIZEOF_TUPLE: usize = 8;
pub const SIZEOF_SYMBOLIC: usize = 264;
pub const SIZEOF_NUMERIC: usize = 312;

/// `UNITS (type, n)`: `CEILING (sizeof (type) * n, sizeof (Unit))` in size_t.
#[inline]
pub fn units(size: usize, n: Int) -> Int {
    let b = (size as u64).wrapping_mul(n as i64 as u64);
    (b.wrapping_add(SIZEOF_UNIT as u64 - 1) / SIZEOF_UNIT as u64) as Int
}

/// `DUNITS (type, n)`: `ceil (sizeof (type) * (double) n / sizeof (Unit))`.
#[inline]
pub fn dunits(size: usize, n: f64) -> f64 {
    (size as f64 * n / SIZEOF_UNIT as f64).ceil()
}

/// C `MAX (a, b)`: `a > b ? a : b`.
#[inline]
pub fn cmax(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

/// C `MIN (a, b)`: `a < b ? a : b`.
#[inline]
pub fn cmin(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}

/// `GET_CONTROL (i, default)`.
#[inline]
pub fn get_control(control: Option<&[f64]>, i: usize, default: f64) -> f64 {
    match control {
        Some(c) if !c[i].is_nan() => c[i],
        _ => default,
    }
}

/// `INT_OVERFLOW (x)`.
#[inline]
pub fn int_overflow(x: f64) -> bool {
    !(x * (1.0 + 1e-8) <= Int::MAX as f64) || x.is_nan()
}

/// `UMFPACK_DENSE_DEGREE_THRESHOLD (alpha, n)`.
#[inline]
pub fn dense_degree_threshold(alpha: f64, n: Int) -> Int {
    cmax(16.0, alpha * 16.0 * (n as f64).sqrt()) as Int
}

/// `UMF_transpose` without the argument checks (`check` = FALSE in every
/// caller on the paths ported here). `ax`/`rx` null skip the values.
pub fn umf_transpose(
    n_row: Int,
    n_col: Int,
    ap: P<Int>,
    ai: P<Int>,
    ax: P<f64>,
    p: P<Int>,
    q: P<Int>,
    nq: Int,
    mut rp: P<Int>,
    mut ri: P<Int>,
    mut rx: P<f64>,
    mut w: P<Int>,
) -> Int {
    for i in 0..n_row {
        w[i] = 0;
        rp[i] = 0;
    }
    if !q.is_null() {
        for newj in 0..nq {
            let j = q[newj];
            for pp in ap[j]..ap[j + 1] {
                w[ai[pp]] += 1;
            }
        }
    } else {
        for j in 0..n_col {
            for pp in ap[j]..ap[j + 1] {
                w[ai[pp]] += 1;
            }
        }
    }
    if !p.is_null() {
        rp[0] = 0;
        for k in 0..n_row {
            let i = p[k];
            rp[k + 1] = rp[k] + w[i];
        }
        for k in 0..n_row {
            let i = p[k];
            w[i] = rp[k];
        }
    } else {
        rp[0] = 0;
        for i in 0..n_row {
            rp[i + 1] = rp[i] + w[i];
        }
        for i in 0..n_row {
            w[i] = rp[i];
        }
    }
    let do_values = !ax.is_null() && !rx.is_null();
    if !q.is_null() {
        for newj in 0..nq {
            let j = q[newj];
            for pp in ap[j]..ap[j + 1] {
                let bp = w[ai[pp]];
                w[ai[pp]] += 1;
                ri[bp] = newj;
                if do_values {
                    rx[bp] = ax[pp];
                }
            }
        }
    } else {
        for j in 0..n_col {
            for pp in ap[j]..ap[j + 1] {
                let bp = w[ai[pp]];
                w[ai[pp]] += 1;
                ri[bp] = j;
                if do_values {
                    rx[bp] = ax[pp];
                }
            }
        }
    }
    UMFPACK_OK
}

/// `UMF_fsize`.
pub fn umf_fsize(
    nn: Int,
    mut fsize: P<Int>,
    fnrows: P<Int>,
    fncols: P<Int>,
    parent: P<Int>,
    npiv: P<Int>,
) {
    for j in 0..nn {
        fsize[j] = EMPTY;
    }
    for j in 0..nn {
        if npiv[j] > 0 {
            let par = parent[j];
            let r = fnrows[j];
            let c = fncols[j];
            let mut frsize = r.wrapping_mul(c);
            if int_overflow(r as f64 * c as f64) {
                frsize = Int::MAX;
            }
            fsize[j] = fsize[j].max(frsize);
            if par != EMPTY {
                fsize[par] = fsize[par].max(fsize[j]);
            }
        }
    }
}

/// `UMF_apply_order`.
pub fn umf_apply_order(mut front: P<Int>, order: P<Int>, mut temp: P<Int>, nn: Int, nfr: Int) {
    for i in 0..nn {
        let k = order[i];
        if k != EMPTY {
            temp[k] = front[i];
        }
    }
    for k in 0..nfr {
        front[k] = temp[k];
    }
}

/// `UMF_is_permutation`.
pub fn umf_is_permutation(p: P<Int>, mut w: P<Int>, n: Int, r: Int) -> bool {
    if p.is_null() {
        return true;
    }
    for i in 0..n {
        w[i] = 0;
    }
    for k in 0..r {
        let i = p[k];
        if i < 0 || i >= n {
            return false;
        }
        if w[i] != 0 {
            return false;
        }
        w[i] = 1;
    }
    true
}

/// `UMF_symbolic_usage`.
pub fn umf_symbolic_usage(
    n_row: Int,
    n_col: Int,
    nchains: Int,
    nfr: Int,
    esize: Int,
    prefer_diagonal: bool,
) -> f64 {
    let mut u = dunits(SIZEOF_SYMBOLIC, 1.0)
        + 2.0 * dunits(SIZEOF_INT, (n_col + 1) as f64)
        + 2.0 * dunits(SIZEOF_INT, (n_row + 1) as f64)
        + 3.0 * dunits(SIZEOF_INT, (nchains + 1) as f64)
        + 4.0 * dunits(SIZEOF_INT, (nfr + 1) as f64);
    u += dunits(SIZEOF_INT, esize as f64);
    if prefer_diagonal {
        u += dunits(SIZEOF_INT, (n_col + 1) as f64);
    }
    u
}

/// `UMF_set_stats`.
#[allow(clippy::too_many_arguments)]
pub fn umf_set_stats(
    info: &mut [f64],
    sym: &Symbolic,
    max_usage: f64,
    num_mem_size: f64,
    flops: f64,
    lnz: f64,
    unz: f64,
    maxfrsize: f64,
    ulen: f64,
    npiv: f64,
    maxnrows: f64,
    maxncols: f64,
    scale: bool,
    prefer_diagonal: bool,
    what: isize,
) {
    let at = |base: usize| (base as isize + what) as usize;
    let n_col = sym.n_col as f64;
    let n_row = sym.n_row as f64;
    let n1 = sym.n1 as f64;
    let nn = cmax(n_row, n_col);
    let n_inner = cmin(n_row, n_col);
    let sym_maxncols = cmin((sym.maxncols + sym.nb) as f64, n_col);
    let sym_maxnrows = cmin((sym.maxnrows + sym.nb) as f64, n_row);
    let elen = (n_col - n1) + (n_row - n1) + cmin(n_col - n1, n_row - n1) + 1.0;
    let sym_size = umf_symbolic_usage(
        sym.n_row,
        sym.n_col,
        sym.nchains,
        sym.nfr,
        sym.esize,
        prefer_diagonal,
    );
    let scale_part = if scale {
        dunits(SIZEOF_ENTRY, n_row)
    } else {
        0.0
    };
    let num_on_size1 = dunits(SIZEOF_NUMERIC, 1.0)
        + dunits(SIZEOF_ENTRY, n_inner + 1.0)
        + 4.0 * dunits(SIZEOF_INT, n_row + 1.0)
        + 4.0 * dunits(SIZEOF_INT, n_col + 1.0)
        + scale_part;
    let num_on_size2 = dunits(SIZEOF_NUMERIC, 1.0)
        + dunits(SIZEOF_ENTRY, n_inner + 1.0)
        + dunits(SIZEOF_INT, n_row + 1.0)
        + dunits(SIZEOF_INT, n_col + 1.0)
        + 6.0 * dunits(SIZEOF_INT, npiv + 1.0)
        + scale_part;
    info[at(UMFPACK_VARIABLE_PEAK)] = max_usage;
    info[at(UMFPACK_VARIABLE_FINAL)] = num_mem_size;
    info[at(UMFPACK_NUMERIC_SIZE)] = num_on_size2 + num_mem_size + dunits(SIZEOF_INT, ulen + 1.0);
    info[at(UMFPACK_MAX_FRONT_SIZE)] = maxfrsize;
    info[at(UMFPACK_MAX_FRONT_NROWS)] = maxnrows;
    info[at(UMFPACK_MAX_FRONT_NCOLS)] = maxncols;
    let work_usage = 2.0 * dunits(SIZEOF_ENTRY, sym_maxnrows + 1.0)
        + 2.0 * dunits(SIZEOF_INT, n_row + 1.0)
        + 2.0 * dunits(SIZEOF_INT, n_col + 1.0)
        + dunits(SIZEOF_INT, nn + 1.0)
        + dunits(SIZEOF_INT, cmax(n_col, sym_maxnrows) + 1.0)
        + 2.0 * dunits(SIZEOF_INT, sym_maxnrows + 1.0)
        + 3.0 * dunits(SIZEOF_INT, sym_maxncols + 1.0)
        + dunits(SIZEOF_INT, cmax(sym_maxnrows, sym_maxncols) + 1.0)
        + dunits(SIZEOF_INT, elen)
        + dunits(SIZEOF_INT, (sym.nfr + 1) as f64)
        + if n_row == n_col {
            2.0 * dunits(SIZEOF_INT, nn)
        } else {
            0.0
        };
    let num_usage = sym_size + num_on_size1 + work_usage + max_usage;
    info[at(UMFPACK_PEAK_MEMORY)] = cmax(sym.peak_sym_usage, num_usage);
    info[at(UMFPACK_FLOPS)] = flops;
    info[at(UMFPACK_LNZ)] = lnz;
    info[at(UMFPACK_UNZ)] = unz;
}
