//! `UMF_2by2`: row permutation that puts large entries on the diagonal by
//! swapping pairs of rows (the 2-by-2 strategy).
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::util::cmax;
use crate::amd::{Int, EMPTY};
use crate::ptr::P;

const IS_WEAK: Int = -2;

#[inline]
fn weak(x: f64, ctol: f64) -> bool {
    (x - ctol) < 0.0
}

#[inline]
fn scalar_abs(x: f64) -> f64 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}

#[allow(clippy::too_many_arguments)]
fn two_by_two(
    n2: Int,
    cp: P<Int>,
    ci: P<Int>,
    degree: P<Int>,
    mut next: P<Int>,
    mut ri: P<Int>,
    mut p: P<Int>,
    mut rp: P<Int>,
    mut head: P<Int>,
) -> Int {
    for deg in 0..n2 {
        head[deg] = EMPTY;
    }
    let mut maxdeg = 0;
    let mut mindeg = Int::MAX;
    for newcol in (0..n2).rev() {
        if next[newcol] == IS_WEAK {
            let deg = degree[newcol];
            next[newcol] = head[deg];
            head[deg] = newcol;
            maxdeg = maxdeg.max(deg);
            mindeg = mindeg.min(deg);
        }
    }
    rp[0] = 0;
    for row in 0..n2 {
        rp[row + 1] = rp[row] + ri[row];
        p[row] = rp[row];
    }
    for col in 0..n2 {
        let p2 = cp[col + 1];
        for pp in cp[col]..p2 {
            let r = ci[pp];
            ri[p[r]] = col;
            p[r] += 1;
        }
    }
    for row in 0..n2 {
        p[row] = EMPTY;
    }
    let mut unmatched = 0;
    let mut jdiff: Int;
    let mut jdeg: Int = EMPTY;
    let mut deg = mindeg;
    while deg <= maxdeg {
        let mut k = head[deg];
        while k != EMPTY {
            if p[k] == EMPTY {
                let mut j_best = EMPTY;
                let mut jdiff_best = Int::MAX;
                let mut jdeg_best = Int::MAX;
                let cp1 = cp[k];
                let cp2 = cp[k + 1];
                let mut c = cp1;
                let rp1 = rp[k];
                let rp2 = rp[k + 1];
                let mut r = rp1;
                jdiff = EMPTY;
                loop {
                    if c >= cp2 {
                        break;
                    }
                    let j = ci[c];
                    if r >= rp2 {
                        break;
                    }
                    let j2 = ri[r];
                    if j < j2 {
                        c += 1;
                        continue;
                    }
                    if j2 < j {
                        r += 1;
                        continue;
                    }
                    let mut best = false;
                    if p[j] == EMPTY {
                        jdeg = degree[j];
                        jdiff = (k - j).abs();
                        if j_best == EMPTY {
                            best = true;
                        } else if jdeg < jdeg_best {
                            best = true;
                        } else if jdeg == jdeg_best {
                            if jdiff < jdiff_best {
                                best = true;
                            } else if jdiff == jdiff_best {
                                best = j < j_best;
                            }
                        } else {
                            best = false;
                        }
                    }
                    if best {
                        j_best = j;
                        jdeg_best = jdeg;
                        jdiff_best = jdiff;
                    }
                    c += 1;
                    r += 1;
                }
                let _ = jdiff;
                if j_best != EMPTY {
                    let j = j_best;
                    p[k] = j;
                    p[j] = k;
                } else {
                    unmatched += 1;
                }
            }
            k = next[k];
        }
        deg += 1;
    }
    for k in 0..n2 {
        if p[k] == EMPTY {
            p[k] = k;
        }
    }
    unmatched
}

/// `UMF_2by2`. Returns `(nweak, unmatched)`.
#[allow(clippy::too_many_arguments)]
pub fn umf_2by2(
    n: Int,
    ap: P<Int>,
    ai: P<Int>,
    ax: P<f64>,
    tol: f64,
    scale: Int,
    cperm1: P<Int>,
    inv_rperm1: P<Int>,
    n1: Int,
    nempty: Int,
    degree: P<Int>,
    mut p: P<Int>,
    mut ri: P<Int>,
    rp: P<Int>,
    mut rs: P<f64>,
    head: P<Int>,
    mut next: P<Int>,
    mut ci: P<Int>,
    mut cp: P<Int>,
) -> (Int, Int) {
    let mut do_recip = false;
    let do_values = tol > 0.0 && !ax.is_null();
    let (do_sum, do_max) = if do_values && !rs.is_null() {
        (scale == UMFPACK_SCALE_SUM, scale == UMFPACK_SCALE_MAX)
    } else {
        (false, false)
    };
    let do_scale = do_max || do_sum;
    if do_scale {
        for row in 0..n {
            rs[row] = 0.0;
        }
        for col in 0..n {
            let p2 = ap[col + 1];
            for pp in ap[col]..p2 {
                let row = ai[pp];
                let value = scalar_abs(ax[pp]);
                let r = rs[row];
                if !r.is_nan() {
                    if value.is_nan() {
                        rs[row] = value;
                    } else if do_max {
                        rs[row] = cmax(r, value);
                    } else {
                        rs[row] += value;
                    }
                }
            }
        }
        let mut rsmin = rs[0];
        if rsmin == 0.0 || rsmin.is_nan() {
            rsmin = 1.0;
        }
        for row in 0..n {
            let r = rs[row];
            if r == 0.0 || r.is_nan() {
                rs[row] = 1.0;
            }
            rsmin = if rsmin < rs[row] { rsmin } else { rs[row] };
        }
        do_recip = rsmin >= RECIPROCAL_TOLERANCE;
        if do_recip {
            for row in 0..n {
                rs[row] = 1.0 / rs[row];
            }
        }
    }
    let mut nweak = 0;
    let n2 = n - n1 - nempty;
    for row in 0..n2 {
        ri[row] = 0;
    }
    let mut pp_out = 0;
    let mut ctol = 0.0;
    let mut dvalue = 1.0;
    // value of one entry, scaled as the column scan does
    let scaled = |pp: Int, oldrow: Int| -> f64 {
        let v = scalar_abs(ax[pp]);
        if !do_scale {
            v
        } else if do_recip {
            v * rs[oldrow]
        } else {
            v / rs[oldrow]
        }
    };
    for k in n1..n - nempty {
        let oldcol = cperm1[k];
        let newcol = k - n1;
        next[newcol] = EMPTY;
        cp[newcol] = pp_out;
        let mut dfound = false;
        let p1 = ap[oldcol];
        let p2 = ap[oldcol + 1];
        if do_values {
            let mut cmx = 0.0;
            dvalue = 0.0;
            for pp in p1..p2 {
                let oldrow = ai[pp];
                let newrow = inv_rperm1[oldrow] - n1;
                if newrow < 0 {
                    continue;
                }
                let value = scaled(pp, oldrow);
                if !f64::is_nan(cmx) {
                    if value.is_nan() {
                        cmx = value;
                    } else {
                        cmx = cmax(cmx, value);
                    }
                }
                if oldrow == oldcol {
                    dvalue = value;
                    dfound = true;
                }
            }
            ctol = tol * cmx;
        } else {
            for pp in p1..p2 {
                let oldrow = ai[pp];
                let newrow = inv_rperm1[oldrow] - n1;
                if newrow < 0 {
                    continue;
                }
                ci[pp_out] = newrow;
                pp_out += 1;
                if oldrow == oldcol {
                    dfound = true;
                }
                ri[newrow] += 1;
            }
        }
        let is_weak = if !dfound {
            true
        } else if do_values {
            weak(dvalue, ctol)
        } else {
            false
        };
        if is_weak {
            next[newcol] = IS_WEAK;
            nweak += 1;
        }
        if do_values {
            for pp in p1..p2 {
                let oldrow = ai[pp];
                let newrow = inv_rperm1[oldrow] - n1;
                if newrow < 0 {
                    continue;
                }
                let value = scaled(pp, oldrow);
                if !weak(value, ctol) {
                    ci[pp_out] = newrow;
                    pp_out += 1;
                    ri[newrow] += 1;
                }
            }
        }
    }
    cp[n2] = pp_out;
    if nweak == 0 {
        for k in 0..n {
            p[k] = k;
        }
        return (0, 0);
    }
    let unmatched = two_by_two(n2, cp, ci, degree, next, ri, p, rp, head);
    (nweak, unmatched)
}
