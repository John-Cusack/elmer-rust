//! `UMF_analyze`: symbolic factorisation of `A'A` (column elimination tree
//! and frontal matrices) for a given column ordering.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::util::{umf_apply_order, umf_fsize};
use crate::amd::{amd_postorder, Int, EMPTY};
use crate::ptr::P;

/// Returns `false` if the workspace is too small (`UMFPACK_ERROR_internal_error`).
#[allow(clippy::too_many_arguments)]
pub fn umf_analyze(
    n_row: Int,
    n_col: Int,
    mut ai: P<Int>,
    ap: P<Int>,
    mut up: P<Int>,
    fix_q: bool,
    mut w: P<Int>,
    mut link: P<Int>,
    mut front_ncols: P<Int>,
    mut front_nrows: P<Int>,
    mut front_npivcol: P<Int>,
    mut front_parent: P<Int>,
    nfr_out: &mut Int,
    p_ncompactions: &mut Int,
) -> bool {
    let mut nfr: Int = 0;
    for j in 0..n_col {
        link[j] = EMPTY;
        w[j] = EMPTY;
        up[j] = EMPTY;
        front_npivcol[j] = 0;
        front_nrows[j] = 0;
        front_ncols[j] = 0;
        front_parent[j] = EMPTY;
    }
    let mut krow: Int = 0;
    let mut pfirst = ap[0];
    let mut jlast = EMPTY;
    let mut jnext;
    let mut wflag: Int = 0;
    let mut pdest: Int = 0;
    let mut ncompactions = 0;
    let mut j: Int = 0;
    while j < n_col {
        if pdest + (n_col - j) > pfirst {
            pdest = 0;
            ncompactions += 1;
            for row in 0..j {
                if up[row] != EMPTY {
                    let mut p = up[row];
                    let p2 = p + (front_ncols[row] - front_npivcol[row]);
                    up[row] = pdest;
                    while p < p2 {
                        ai[pdest] = ai[p];
                        pdest += 1;
                        p += 1;
                    }
                }
            }
        }
        if pdest + (n_col - j) > pfirst {
            return false;
        }
        let mut parent;
        let mut thickness;
        if jlast != EMPTY && link[j] == jlast {
            up[j] = up[jlast];
            up[jlast] = EMPTY;
            parent = n_col;
            let mut p = up[j];
            while p < pdest {
                let j3 = ai[p];
                if j == j3 {
                    pdest -= 1;
                    ai[p] = ai[pdest];
                } else {
                    if j3 < parent {
                        parent = j3;
                    }
                    p += 1;
                }
            }
            link[j] = link[jlast];
            thickness = front_nrows[jlast] - front_npivcol[jlast];
        } else {
            up[j] = pdest;
            parent = n_col;
            thickness = 0;
            wflag = j;
        }
        w[j] = wflag;
        jnext = n_col;
        let mut knext = krow;
        while knext < n_row {
            jnext = ai[ap[knext]];
            if jnext != j {
                break;
            }
            knext += 1;
        }
        if knext == n_row {
            jnext = n_col;
        }
        for k in krow..knext {
            let p2 = ap[k + 1];
            for p in ap[k]..p2 {
                let col = ai[p];
                if w[col] != wflag {
                    ai[pdest] = col;
                    pdest += 1;
                    w[col] = wflag;
                    if col < parent {
                        parent = col;
                    }
                }
            }
            thickness += 1;
        }
        krow = knext;
        pfirst = ap[knext];
        let mut k = link[j];
        while k != EMPTY {
            let mut p = up[k];
            let p2 = p + (front_ncols[k] - front_npivcol[k]);
            while p < p2 {
                let col = ai[p];
                if w[col] != wflag {
                    ai[pdest] = col;
                    pdest += 1;
                    w[col] = wflag;
                    if col < parent {
                        parent = col;
                    }
                }
                p += 1;
            }
            up[k] = EMPTY;
            thickness += front_nrows[k] - front_npivcol[k];
            k = link[k];
        }
        let mut j2 = j + 1;
        while j2 < jnext {
            if w[j2] != wflag || link[j2] != EMPTY {
                break;
            }
            j2 += 1;
        }
        jnext = j2;
        j2 -= 1;
        let npivots = j2 - j + 1;
        if j2 > j {
            parent = n_col;
            let p2 = pdest;
            pdest = up[j];
            for p in up[j]..p2 {
                let col = ai[p];
                if col > j2 {
                    ai[pdest] = col;
                    pdest += 1;
                    if col < parent {
                        parent = col;
                    }
                }
            }
        }
        if parent == n_col {
            parent = EMPTY;
        }
        let npivcol = npivots;
        let fallrows = thickness;
        let fallcols = npivots + pdest - up[j];
        let fpiv = npivcol.min(fallrows);
        let frows = fallrows - fpiv;
        let fcols = fallcols - fpiv;
        if frows == 0 || fcols == 0 {
            up[j] = EMPTY;
            parent = EMPTY;
        }
        front_npivcol[j] = npivots;
        front_nrows[j] = fallrows;
        front_ncols[j] = fallcols;
        front_parent[j] = parent;
        nfr += 1;
        if parent != EMPTY {
            link[j] = link[parent];
            link[parent] = j;
        }
        jlast = j;
        j = jnext;
    }
    *nfr_out = nfr;

    let mut front_order = w;
    if fix_q {
        let mut k = 0;
        for j in 0..n_col {
            if front_npivcol[j] > 0 {
                front_order[j] = k;
                k += 1;
            } else {
                front_order[j] = EMPTY;
            }
        }
    } else {
        let front_child = ap;
        let front_sibling = link;
        let front_stack = ai;
        let front_size = front_stack + n_col;
        umf_fsize(
            n_col,
            front_size,
            front_nrows,
            front_ncols,
            front_parent,
            front_npivcol,
        );
        {
            let n = n_col as usize;
            let sl = |q: P<Int>| -> &mut [Int] {
                if n == 0 {
                    return &mut [];
                }
                // SAFETY: disjoint n_col-long pieces, as in the C code.
                unsafe { core::slice::from_raw_parts_mut(q.r(0) as *mut Int, n) }
            };
            amd_postorder(
                n_col,
                sl(front_parent),
                sl(front_npivcol),
                sl(front_size),
                sl(front_order),
                sl(front_child),
                sl(front_sibling),
                sl(front_stack),
            );
        }
        let mut winv = ai;
        for k in 0..nfr {
            winv[k] = EMPTY;
        }
        for i in 0..n_col {
            let k = front_order[i];
            if k != EMPTY {
                winv[k] = i;
            }
        }
        let mut kk = 0;
        for k in 0..nfr {
            let i = winv[k];
            for npiv in 0..front_npivcol[i] {
                up[kk] = i + npiv;
                kk += 1;
            }
        }
    }
    umf_apply_order(front_npivcol, front_order, ai, n_col, nfr);
    umf_apply_order(front_nrows, front_order, ai, n_col, nfr);
    umf_apply_order(front_ncols, front_order, ai, n_col, nfr);
    umf_apply_order(front_parent, front_order, ai, n_col, nfr);
    for i in 0..nfr {
        let parent = front_parent[i];
        if parent != EMPTY {
            front_parent[i] = front_order[parent];
        }
    }
    *p_ncompactions = ncompactions;
    true
}
