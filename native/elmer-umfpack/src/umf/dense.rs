//! Dense kernels on the frontal matrix (`NBLAS` build): `UMF_scale`,
//! `UMF_blas3_update`, `UMF_scale_column`.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::numtypes::*;
use crate::amd::{Int, EMPTY};
use crate::ptr::P;

/// `SCALAR_ABS`.
#[inline(always)]
pub fn sabs(x: f64) -> f64 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}

/// `UMF_scale`: `X [0..n-1] /= pivot` (by the reciprocal unless tiny).
pub fn umf_scale(n: Int, mut pivot: Entry, mut x: P<Entry>) {
    let s = sabs(pivot);
    if s < RECIPROCAL_TOLERANCE || pivot.is_nan() {
        for i in 0..n {
            let xi = x[i];
            if xi != 0.0 {
                x[i] = xi / pivot;
            }
        }
    } else {
        pivot = 1.0 / pivot;
        for i in 0..n {
            let xi = x[i];
            x[i] = xi * pivot;
        }
    }
}

/// `UMF_blas3_update` (C loops, no BLAS).
pub fn umf_blas3_update(work: &Work) {
    let k = work.fnpiv;
    if k == 0 {
        return;
    }
    let m = work.fnrows;
    let n = work.fncols;
    let d = work.fnr_curr;
    let dc = work.fnc_curr;
    let nb = work.nb;
    let c = work.fcblock;
    let l = work.flblock;
    let u = work.fublock;
    let lu = work.flublock;
    if k == 1 {
        for j in 0..n {
            let u_j = u[j];
            if u_j != 0.0 {
                let mut c_ij = c + j * d;
                let mut l_is = l;
                for _i in 0..m {
                    *c_ij.r(0) -= l_is[0] * u_j;
                    c_ij = c_ij + 1;
                    l_is = l_is + 1;
                }
            }
        }
    } else {
        for s in 0..k {
            for i in s + 1..k {
                let l_is = lu[i + s * nb];
                if l_is != 0.0 {
                    let mut u_ij = u + i * dc;
                    let mut u_sj = u + s * dc;
                    for _j in 0..n {
                        *u_ij.r(0) -= l_is * u_sj[0];
                        u_ij = u_ij + 1;
                        u_sj = u_sj + 1;
                    }
                }
            }
        }
        for s in 0..k {
            for j in 0..n {
                let u_sj = u[j + s * dc];
                if u_sj != 0.0 {
                    let mut c_ij = c + j * d;
                    let mut l_is = l + s * d;
                    for _i in 0..m {
                        *c_ij.r(0) -= l_is[0] * u_sj;
                        c_ij = c_ij + 1;
                        l_is = l_is + 1;
                    }
                }
            }
        }
    }
}

/// `UMF_scale_column`: move the pivot row and column out of the
/// contribution block, scale the pivot column, record the pivot.
pub fn umf_scale_column(numeric: &mut Numeric, work: &mut Work) {
    let mut fnrows = work.fnrows;
    let mut fncols;
    let fnpiv = work.fnpiv;
    let mut rperm = numeric.rperm_p();
    let mut cperm = numeric.cperm_p();
    let flublock = work.flublock;
    let flblock = work.flblock;
    let fublock = work.fublock;
    let fcblock = work.fcblock;
    let fnr_curr = work.fnr_curr;
    let fnc_curr = work.fnc_curr;
    let mut frpos = work.frpos;
    let mut fcpos = work.fcpos;
    let mut frows = work.frows;
    let mut fcols = work.fcols;
    let pivrow = work.pivrow;
    let pivcol = work.pivcol;
    let mut row_tuples = numeric.uip_p();
    let mut col_tuples = numeric.lip_p();
    let nb = work.nb;

    let mut fspos = fcpos[pivcol];
    work.fncols -= 1;
    fncols = work.fncols;
    if fspos != fncols * fnr_curr {
        let fs = fspos / fnr_curr;
        {
            let mut f_s = fcblock + fspos;
            let f_e = fcblock + fncols * fnr_curr;
            for i in 0..fnrows {
                f_s[i] = f_e[i];
            }
        }
        {
            let mut f_s = fublock + fs;
            let f_e = fublock + fncols;
            for i in 0..fnpiv {
                f_s[i * fnc_curr] = f_e[i * fnc_curr];
            }
        }
        let col2 = fcols[fncols];
        fcols[fs] = col2;
        fcpos[col2] = fspos;
    }
    fcpos[pivcol] = EMPTY;

    fspos = frpos[pivrow];
    work.fnrows -= 1;
    fnrows = work.fnrows;
    let _ = &mut fncols;
    if fspos == fnrows {
        {
            let mut fd = fublock + fnpiv * fnc_curr;
            let fsb = fcblock + fspos;
            for j in 0..fncols {
                fd[j] = fsb[j * fnr_curr];
            }
        }
        if work.pivrow_in_front {
            let mut fd = flublock + fnpiv;
            let fsb = flblock + fspos;
            for j in 0..=fnpiv {
                fd[j * nb] = fsb[j * fnr_curr];
            }
        } else {
            let mut fd = flublock + fnpiv;
            let fsb = flblock + fspos;
            for j in 0..fnpiv {
                fd[j * nb] = 0.0;
            }
            fd[fnpiv * nb] = fsb[fnpiv * fnr_curr];
        }
    } else {
        {
            // shift_pivot_row
            let mut fd = fublock + fnpiv * fnc_curr;
            let mut fsb = fcblock + fspos;
            let fe = fcblock + fnrows;
            for j in 0..fncols {
                fd[j] = fsb[j * fnr_curr];
                fsb[j * fnr_curr] = fe[j * fnr_curr];
            }
        }
        if work.pivrow_in_front {
            let mut fd = flublock + fnpiv;
            let mut fsb = flblock + fspos;
            let fe = flblock + fnrows;
            for j in 0..=fnpiv {
                fd[j * nb] = fsb[j * fnr_curr];
                fsb[j * fnr_curr] = fe[j * fnr_curr];
            }
        } else {
            let mut fd = flublock + fnpiv;
            let mut fsb = flblock + fspos;
            let fe = flblock + fnrows;
            for j in 0..fnpiv {
                fd[j * nb] = 0.0;
                fsb[j * fnr_curr] = fe[j * fnr_curr];
            }
            fd[fnpiv * nb] = fsb[fnpiv * fnr_curr];
            fsb[fnpiv * fnr_curr] = fe[fnpiv * fnr_curr];
        }
        let row2 = frows[fnrows];
        frows[fspos] = row2;
        frpos[row2] = fspos;
    }
    frpos[pivrow] = EMPTY;

    let fcol = flblock + fnpiv * fnr_curr;
    let pivot_value = flublock[fnpiv + fnpiv * nb];
    let k = work.npiv + fnpiv;
    umf_scale(fnrows, pivot_value, fcol);
    let (rt, ct) = (row_tuples[pivrow], col_tuples[pivcol]);
    umf_mem_free_tail_block(numeric, rt);
    umf_mem_free_tail_block(numeric, ct);
    row_tuples[pivrow] = 0;
    col_tuples[pivcol] = 0;
    let k1 = -k - 1;
    rperm[pivrow] = k1;
    cperm[pivcol] = k1;
    work.pivrow_arr[fnpiv as usize] = pivrow;
    work.pivcol_arr[fnpiv as usize] = pivcol;
    work.fnpiv += 1;
}
