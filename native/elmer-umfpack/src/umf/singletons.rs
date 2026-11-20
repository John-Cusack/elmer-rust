//! `UMF_singletons`: find column and row singletons (no user ordering,
//! the only case Elmer reaches).
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use crate::amd::{flip, Int, EMPTY};
use crate::ptr::P;

fn create_row_form(
    n_row: Int,
    n_col: Int,
    ap: P<Int>,
    ai: P<Int>,
    rdeg: P<Int>,
    mut rp: P<Int>,
    mut ri: P<Int>,
    mut w: P<Int>,
) {
    rp[0] = 0;
    w[0] = 0;
    for row in 0..n_row {
        rp[row + 1] = rp[row] + rdeg[row];
        w[row] = rp[row];
    }
    for col in 0..n_col {
        let p2 = ap[col + 1];
        for p in ap[col]..p2 {
            let r = ai[p];
            ri[w[r]] = col;
            w[r] += 1;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn order_singletons(
    mut k: Int,
    mut head: Int,
    mut tail: Int,
    mut next: P<Int>,
    mut xdeg: P<Int>,
    mut xperm: P<Int>,
    xp: P<Int>,
    xi: P<Int>,
    mut ydeg: P<Int>,
    mut yperm: P<Int>,
    yp: P<Int>,
    yi: P<Int>,
) -> Int {
    while head != EMPTY {
        let xpivot = head;
        head = next[xpivot];
        if head == EMPTY {
            tail = EMPTY;
        }
        if xdeg[xpivot] != 1 {
            continue;
        }
        let mut ypivot = EMPTY;
        let p2 = xp[xpivot + 1];
        for p in xp[xpivot]..p2 {
            let y = xi[p];
            if ydeg[y] >= 0 {
                ypivot = y;
                break;
            }
        }
        let p2 = yp[ypivot + 1];
        for p in yp[ypivot]..p2 {
            let x = yi[p];
            if xdeg[x] < 0 {
                continue;
            }
            if x == xpivot {
                continue;
            }
            xdeg[x] -= 1;
            let deg = xdeg[x];
            if deg == 1 {
                next[x] = EMPTY;
                if head == EMPTY {
                    head = x;
                } else {
                    next[tail] = x;
                }
                tail = x;
            }
        }
        xdeg[xpivot] = flip(1);
        ydeg[ypivot] = flip(ydeg[ypivot]);
        xperm[k] = xpivot;
        yperm[k] = ypivot;
        k += 1;
    }
    k
}

#[allow(clippy::too_many_arguments)]
fn find_any_singletons(
    n_row: Int,
    n_col: Int,
    ap: P<Int>,
    ai: P<Int>,
    cdeg: P<Int>,
    rdeg: P<Int>,
    cperm: P<Int>,
    rperm: P<Int>,
    p_n1r: &mut Int,
    p_n1c: &mut Int,
    rp: P<Int>,
    ri: P<Int>,
    w: P<Int>,
    mut next: P<Int>,
) -> Int {
    let mut n1 = 0;
    let mut n1r = 0;
    let mut n1c = 0;
    let mut row_form = false;
    let mut head = EMPTY;
    let mut tail = EMPTY;
    for col in (0..n_col).rev() {
        if cdeg[col] == 1 {
            if head == EMPTY {
                tail = col;
            }
            next[col] = head;
            head = col;
        }
    }
    if head != EMPTY {
        create_row_form(n_row, n_col, ap, ai, rdeg, rp, ri, w);
        row_form = true;
        n1 = order_singletons(
            0, head, tail, next, cdeg, cperm, ap, ai, rdeg, rperm, rp, ri,
        );
        n1c = n1;
    }
    head = EMPTY;
    tail = EMPTY;
    for row in (0..n_row).rev() {
        if rdeg[row] == 1 {
            if head == EMPTY {
                tail = row;
            }
            next[row] = head;
            head = row;
        }
    }
    if head != EMPTY {
        if !row_form {
            create_row_form(n_row, n_col, ap, ai, rdeg, rp, ri, w);
        }
        n1 = order_singletons(
            n1, head, tail, next, rdeg, rperm, rp, ri, cdeg, cperm, ap, ai,
        );
        n1r = n1 - n1c;
    }
    *p_n1r = n1r;
    *p_n1c = n1c;
    n1
}

fn finish_permutation(
    n1: Int,
    nx: Int,
    mut xdeg: P<Int>,
    mut xperm: P<Int>,
    p_max_deg: &mut Int,
) -> Int {
    let mut nempty = 0;
    let mut s = n1;
    let mut max_deg = 0;
    for k in 0..nx {
        let x = k;
        let deg = xdeg[x];
        if deg == 0 {
            nempty += 1;
            xperm[nx - nempty] = x;
        } else if deg > 0 {
            xperm[s] = x;
            s += 1;
            max_deg = max_deg.max(deg);
        } else {
            xdeg[x] = flip(deg);
        }
    }
    *p_max_deg = max_deg;
    nempty
}

/// Outputs of `UMF_singletons`.
#[derive(Default, Debug)]
pub struct SingletonInfo {
    pub n1: Int,
    pub n1c: Int,
    pub n1r: Int,
    pub nempty_col: Int,
    pub nempty_row: Int,
    pub is_sym: bool,
    pub max_rdeg: Int,
}

/// `UMF_singletons` with `Quser == NULL`.
#[allow(clippy::too_many_arguments)]
pub fn umf_singletons(
    n_row: Int,
    n_col: Int,
    ap: P<Int>,
    ai: P<Int>,
    mut cdeg: P<Int>,
    cperm: P<Int>,
    mut rdeg: P<Int>,
    rperm: P<Int>,
    mut inv_rperm: P<Int>,
    out: &mut SingletonInfo,
    rp: P<Int>,
    ri: P<Int>,
    w: P<Int>,
    next: P<Int>,
) -> Int {
    if ap[0] != 0 || ap[n_col] < 0 {
        return UMFPACK_ERROR_invalid_matrix;
    }
    for row in 0..n_row {
        rdeg[row] = 0;
    }
    for col in 0..n_col {
        let p1 = ap[col];
        let p2 = ap[col + 1];
        let cd = p2 - p1;
        if cd < 0 {
            return UMFPACK_ERROR_invalid_matrix;
        }
        let mut last_row = EMPTY;
        for p in p1..p2 {
            let row = ai[p];
            if row <= last_row || row >= n_row {
                return UMFPACK_ERROR_invalid_matrix;
            }
            rdeg[row] += 1;
            last_row = row;
        }
        cdeg[col] = cd;
    }
    let (mut n1r, mut n1c) = (0, 0);
    let n1 = find_any_singletons(
        n_row, n_col, ap, ai, cdeg, rdeg, cperm, rperm, &mut n1r, &mut n1c, rp, ri, w, next,
    );
    let mut max_cdeg = 0;
    let nempty_col = finish_permutation(n1, n_col, cdeg, cperm, &mut max_cdeg);
    let mut max_rdeg = 0;
    let nempty_row = finish_permutation(n1, n_row, rdeg, rperm, &mut max_rdeg);
    for k in 0..n_row {
        inv_rperm[rperm[k]] = k;
    }
    let mut is_sym;
    if n_row == n_col && nempty_row == nempty_col {
        is_sym = true;
        for s in n1..n_col - nempty_col {
            if cperm[s] != rperm[s] {
                is_sym = false;
                break;
            }
        }
    } else {
        is_sym = false;
    }
    *out = SingletonInfo {
        n1,
        n1c,
        n1r,
        nempty_col,
        nempty_row,
        is_sym,
        max_rdeg,
    };
    UMFPACK_OK
}
