//! Frontal matrix management: `UMF_start_front`, `UMF_grow_front`,
//! `UMF_init_front`, `UMF_extend_front`, `UMF_create_element`.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::gc::umf_get_memory;
use super::numtypes::*;
use super::symbolic::Symbolic;
use super::util::*;
use crate::amd::{flip, Int, EMPTY};
use crate::ptr::P;

/// `UMF_start_front`.
pub fn umf_start_front(
    chain: Int,
    numeric: &mut Numeric,
    work: &mut Work,
    symbolic: &Symbolic,
) -> bool {
    let nb = symbolic.nb;
    let fnrows_max = symbolic.chain_maxrows[chain as usize];
    let fncols_max = symbolic.chain_maxcols[chain as usize];
    work.fnrows_max = fnrows_max;
    work.fncols_max = fncols_max;
    work.any_skip = false;
    let maxbytes = SIZEOF_ENTRY as f64 * (fnrows_max + nb) as f64 * (fncols_max + nb) as f64;
    let fcurr_size = work.fcurr_size;
    let mut cdeg: Int;
    if symbolic.prefer_diagonal {
        let e_arr = work.e;
        let m = numeric.mem();
        let col_tuples = numeric.lip_p();
        let col_tlen = numeric.lilen_p();
        let col = work.nextcand;
        let tpi = col_tuples[col];
        let tp0: P<Tuple> = m.cast::<Tuple>() + tpi;
        cdeg = 0;
        for t in 0..col_tlen[col] {
            let tp = *tp0.r(t);
            let e = tp.e;
            if e_arr[e] == 0 {
                continue;
            }
            let f = tp.f;
            let p = m + e_arr[e];
            let ep: P<Element> = p.cast();
            let cols: P<Int> = (p + UNITS_ELEMENT).cast();
            if cols[f] == EMPTY {
                continue;
            }
            cdeg += ep.r(0).nrowsleft;
        }
        if symbolic.amd_dmax > 0 {
            cdeg = cdeg.min(symbolic.amd_dmax);
        }
        cdeg += 2;
        cdeg = cdeg.min(fnrows_max);
    } else {
        cdeg = 0;
    }
    let overflow = int_overflow(maxbytes);
    let maxfrsize = if overflow {
        (Int::MAX as usize / SIZEOF_ENTRY) as Int
    } else {
        (fnrows_max + nb) * (fncols_max + nb)
    };
    let mut fsize: Int;
    if numeric.front_alloc_init < 0.0 {
        fsize = (-numeric.front_alloc_init) as Int;
        fsize = fsize.max(1);
    } else {
        if int_overflow(numeric.front_alloc_init * maxbytes) {
            fsize = (Int::MAX as usize / SIZEOF_ENTRY) as Int;
        } else {
            fsize = (numeric.front_alloc_init * maxfrsize as f64) as Int;
        }
        if cdeg > 0 {
            cdeg += nb;
            let fsize2 = if int_overflow((cdeg as f64 * cdeg as f64) * SIZEOF_ENTRY as f64) {
                (Int::MAX as usize / SIZEOF_ENTRY) as Int
            } else {
                (cdeg * cdeg).max(fcurr_size)
            };
            fsize = fsize.min(fsize2);
        }
    }
    fsize = fsize.max(2 * nb * nb);
    work.fnrows_new = 0;
    work.fncols_new = 0;
    let mut fnr2;
    let mut fnc2;
    if fsize >= maxfrsize && !overflow {
        fnr2 = fnrows_max + nb;
        fnc2 = fncols_max + nb;
        fsize = maxfrsize;
    } else if fnrows_max <= fncols_max {
        fnr2 = (fsize as f64).sqrt() as Int;
        fnr2 = fnr2.max(1);
        if fnr2 % 2 == 0 {
            fnr2 += 1;
        }
        fnr2 = fnr2.min(fnrows_max + nb);
        fnc2 = fsize / fnr2;
    } else {
        fnc2 = (fsize as f64).sqrt() as Int;
        fnc2 = fnc2.min(fncols_max + nb);
        fnr2 = fsize / fnc2;
        fnr2 = fnr2.max(1);
        if fnr2 % 2 == 0 {
            fnr2 += 1;
            fnc2 = fsize / fnr2;
        }
    }
    fnr2 = fnr2.min(fnrows_max + nb);
    fnc2 = fnc2.min(fncols_max + nb);
    fnr2 -= nb;
    fnc2 -= nb;
    if fsize > fcurr_size {
        work.do_grow = true;
        if !umf_grow_front(numeric, fnr2, fnc2, work, -1) {
            return false;
        }
    } else {
        work.fnr_curr = fnr2;
        work.fnc_curr = fnc2;
        work.flblock = work.flublock + nb * nb;
        work.fublock = work.flblock + nb * fnr2;
        work.fcblock = work.fublock + nb * fnc2;
    }
    true
}

/// C `MAX (Int, double)` / `MIN` assigned back to an `Int`.
#[inline]
fn imax_f(a: Int, b: f64) -> Int {
    (if a as f64 > b { a as f64 } else { b }) as Int
}
#[inline]
fn imin_f(a: Int, b: f64) -> Int {
    (if (a as f64) < b { a as f64 } else { b }) as Int
}

/// `UMF_grow_front`.
pub fn umf_grow_front(
    numeric: &mut Numeric,
    mut fnr2: Int,
    mut fnc2: Int,
    work: &mut Work,
    do_what: Int,
) -> bool {
    let fcols = work.fcols;
    let mut fcpos = work.fcpos;
    let mut e_arr = work.e;
    let nb = work.nb;
    let fnrows_max = work.fnrows_max + nb;
    let fncols_max = work.fncols_max + nb;
    let mut fnrows_new = work.fnrows_new + 1;
    let mut fncols_new = work.fncols_new + 1;
    if fnrows_new % 2 == 0 {
        fnrows_new += 1;
    }
    fnrows_new += nb;
    fncols_new += nb;
    let fnr_min = fnrows_new.min(fnrows_max);
    let fnc_min = fncols_new.min(fncols_max);
    let minsize = fnr_min * fnc_min;
    if int_overflow(fnr_min as f64 * fnc_min as f64 * SIZEOF_ENTRY as f64) {
        return false;
    }
    fnr2 += nb;
    fnc2 += nb;
    if fnr2 % 2 == 0 {
        fnr2 += 1;
    }
    fnr2 = fnr2.max(fnr_min);
    fnc2 = fnc2.max(fnc_min);
    fnr2 = fnr2.min(fnrows_max);
    fnc2 = fnc2.min(fncols_max);
    let s = fnr2 as f64 * fnc2 as f64;
    let mut newsize;
    if int_overflow(s * SIZEOF_ENTRY as f64) {
        let a = 0.9 * (((Int::MAX as usize / SIZEOF_ENTRY) as f64) / s).sqrt();
        fnr2 = imax_f(fnr_min, a * fnr2 as f64);
        fnc2 = imax_f(fnc_min, a * fnc2 as f64);
        newsize = fnr2 * fnc2;
        if fnr2 % 2 == 0 {
            fnr2 += 1;
        }
        fnc2 = newsize / fnr2;
    }
    fnr2 = fnr2.max(fnr_min);
    fnc2 = fnc2.max(fnc_min);
    newsize = fnr2 * fnc2;
    if e_arr[0] != 0 && do_what != 1 {
        let e0 = e_arr[0];
        super::numtypes::umf_mem_free_tail_block(numeric, e0);
        e_arr[0] = 0;
        work.flublock = P::null();
        work.flblock = P::null();
        work.fublock = P::null();
        work.fcblock = P::null();
    }
    let mut eloc = umf_mem_alloc_tail_block(numeric, units(SIZEOF_ENTRY, newsize));
    if eloc == 0 {
        let (fr, fc) = (work.fnrows, work.fncols);
        if !umf_get_memory(
            numeric,
            work,
            1 + units(SIZEOF_ENTRY, newsize),
            fr,
            fc,
            false,
        ) {
            return false;
        }
        eloc = umf_mem_alloc_tail_block(numeric, units(SIZEOF_ENTRY, newsize));
    }
    while (fnr2 != fnr_min || fnc2 != fnc_min) && eloc == 0 {
        fnr2 = imin_f(fnr2 - 2, fnr2 as f64 * UMF_REALLOC_REDUCTION);
        fnc2 = imin_f(fnc2 - 2, fnc2 as f64 * UMF_REALLOC_REDUCTION);
        fnr2 = fnr_min.max(fnr2);
        fnc2 = fnc_min.max(fnc2);
        if fnr2 % 2 == 0 {
            fnr2 += 1;
        }
        newsize = fnr2 * fnc2;
        eloc = umf_mem_alloc_tail_block(numeric, units(SIZEOF_ENTRY, newsize));
    }
    if eloc == 0 {
        fnr2 = fnr_min;
        fnc2 = fnc_min;
        newsize = minsize;
        eloc = umf_mem_alloc_tail_block(numeric, units(SIZEOF_ENTRY, newsize));
    }
    if eloc == 0 {
        return false;
    }
    let fnr_curr = work.fnr_curr;
    let fnrows = work.fnrows;
    let fncols = work.fncols;
    let mut fcold = work.fcblock;
    fnr2 -= nb;
    fnc2 -= nb;
    work.flublock = (numeric.mem() + eloc).cast();
    work.flblock = work.flublock + nb * nb;
    work.fublock = work.flblock + nb * fnr2;
    work.fcblock = work.fublock + nb * fnc2;
    let mut fcnew = work.fcblock;
    if e_arr[0] != 0 {
        for j in 0..fncols {
            let col = fcols[j];
            for i in 0..fnrows {
                fcnew[i] = fcold[i];
            }
            fcnew = fcnew + fnr2;
            fcold = fcold + fnr_curr;
            fcpos[col] = j * fnr2;
        }
    } else if do_what == 2 {
        for j in 0..fncols {
            let col = fcols[j];
            fcpos[col] = j * fnr2;
        }
    }
    let e0 = e_arr[0];
    super::numtypes::umf_mem_free_tail_block(numeric, e0);
    e_arr[0] = eloc;
    work.fnr_curr = fnr2;
    work.fnc_curr = fnc2;
    work.fcurr_size = newsize;
    work.do_grow = false;
    true
}

/// `UMF_init_front`.
pub fn umf_init_front(numeric: &mut Numeric, work: &mut Work) -> bool {
    if work.do_grow {
        let fnr2 = (UMF_FRONTAL_GROWTH * work.fnrows_new as f64 + 2.0) as Int;
        let fnc2 = (UMF_FRONTAL_GROWTH * work.fncols_new as f64 + 2.0) as Int;
        let what = if work.pivrow_in_front { 2 } else { 0 };
        if !umf_grow_front(numeric, fnr2, fnc2, work, what) {
            return false;
        }
    }
    let fnr_curr = work.fnr_curr;
    let mut frows = work.frows;
    let mut fcols = work.fcols;
    let mut frpos = work.frpos;
    let mut fcpos = work.fcpos;
    work.fnzeros = 0;
    let ccdeg = work.ccdeg;
    let rrdeg = work.rrdeg;
    let mut fnrows = work.fnrows;
    let mut fncols = work.fncols;
    let mut fl = work.flblock;
    if work.pivcol_in_front {
        work.fscan_row = fnrows;
        work.new_rows = work.wrp;
        let wy = work.wy;
        for i in 0..fnrows {
            fl[i] = wy[i];
        }
        let fnrows_extended = fnrows + ccdeg;
        for i in fnrows..fnrows_extended {
            fl[i] = wy[i];
            let row = frows[i];
            work.new_rows.set(i, flip(row));
        }
        fnrows = fnrows_extended;
    } else {
        work.fscan_row = 0;
        work.new_rows = frows;
        let wm = work.wm;
        let wx = work.wx;
        for i in 0..ccdeg {
            fl[i] = wx[i];
            let row = wm[i];
            frows[i] = row;
            frpos[row] = i;
        }
        fnrows = ccdeg;
    }
    work.fnrows = fnrows;
    let wrow = work.wrow;
    if work.pivrow_in_front {
        work.fscan_col = fncols;
        work.new_cols = work.wp;
        if wrow == fcols {
            for j in fncols..rrdeg {
                let col = wrow[j];
                work.new_cols.set(j, flip(col));
                fcpos[col] = j * fnr_curr;
            }
        } else {
            for j in fncols..rrdeg {
                let col = wrow[j];
                fcols[j] = col;
                work.new_cols.set(j, flip(col));
                fcpos[col] = j * fnr_curr;
            }
        }
    } else {
        work.fscan_col = 0;
        work.new_cols = fcols;
        for j in 0..rrdeg {
            let col = wrow[j];
            fcols[j] = col;
            fcpos[col] = j * fnr_curr;
        }
    }
    fncols = rrdeg;
    work.fncols = fncols;
    // zero_init_front (fncols, fnrows, Fcblock, fnr_curr)
    let mut fj = work.fcblock;
    for _j in 0..fncols {
        let mut f = fj;
        fj = fj + fnr_curr;
        for _i in 0..fnrows {
            f.set(0, 0.0);
            f = f + 1;
        }
    }
    let _ = &mut fl;
    true
}

/// `UMF_extend_front`.
pub fn umf_extend_front(numeric: &mut Numeric, work: &mut Work) -> bool {
    let fnpiv = work.fnpiv;
    if work.do_grow {
        let fnr2 = (UMF_FRONTAL_GROWTH * work.fnrows_new as f64 + 2.0) as Int;
        let fnc2 = (UMF_FRONTAL_GROWTH * work.fncols_new as f64 + 2.0) as Int;
        if !umf_grow_front(numeric, fnr2, fnc2, work, 1) {
            return false;
        }
    }
    let fnr_curr = work.fnr_curr;
    let fnc_curr = work.fnc_curr;
    let mut frows = work.frows;
    let mut frpos = work.frpos;
    let mut fcols = work.fcols;
    let mut fcpos = work.fcpos;
    let fnrows = work.fnrows;
    let fncols = work.fncols;
    let rrdeg = work.rrdeg;
    let ccdeg = work.ccdeg;
    work.fscan_col = fncols;
    work.new_cols = fcols;
    work.fscan_row = fnrows;
    work.new_rows = frows;
    let mut fnrows_extended = fnrows;
    let mut fncols_extended = fncols;
    let mut fl = work.flblock + fnpiv * fnr_curr;
    if work.pivcol_in_front {
        fnrows_extended += ccdeg;
        let wy = work.wy;
        for i in 0..fnrows_extended {
            fl[i] = wy[i];
        }
    } else {
        let fu = work.flublock + fnpiv * work.nb;
        let wm = work.wm;
        let wx = work.wx;
        let mut f = fu;
        for _ in 0..fnpiv {
            f.set(0, 0.0);
            f = f + 1;
        }
        let mut f = fl;
        for _ in 0..fnrows {
            f.set(0, 0.0);
            f = f + 1;
        }
        for i in 0..ccdeg {
            let row = wm[i];
            let mut pos = frpos[row];
            if pos < 0 {
                pos = fnrows_extended;
                fnrows_extended += 1;
                frows[pos] = row;
                frpos[row] = pos;
            }
            fl[pos] = wx[i];
        }
    }
    if work.pivrow_in_front {
        if work.pivcol_in_front {
            for j in fncols..rrdeg {
                fcpos[fcols[j]] = j * fnr_curr;
            }
        } else {
            let wrow = work.wrow;
            if wrow == fcols {
                for j in fncols..rrdeg {
                    let col = wrow[j];
                    fcpos[col] = j * fnr_curr;
                }
            } else {
                for j in fncols..rrdeg {
                    let col = wrow[j];
                    fcols[j] = col;
                    fcpos[col] = j * fnr_curr;
                }
            }
        }
        fncols_extended = rrdeg;
    } else {
        let wrow = work.wrow;
        for j in 0..rrdeg {
            let col = wrow[j];
            if fcpos[col] < 0 {
                fcols[fncols_extended] = col;
                fcpos[col] = fncols_extended * fnr_curr;
                fncols_extended += 1;
            }
        }
    }
    // zero_front
    {
        let mut fj = work.fcblock + fnrows;
        for _j in 0..fncols {
            let mut f = fj;
            fj = fj + fnr_curr;
            for _i in fnrows..fnrows_extended {
                f.set(0, 0.0);
                f = f + 1;
            }
        }
        fj = fj - fnrows;
        for _j in fncols..fncols_extended {
            let mut f = fj;
            fj = fj + fnr_curr;
            for _i in 0..fnrows_extended {
                f.set(0, 0.0);
                f = f + 1;
            }
        }
        let mut fj = work.flblock + fnrows;
        for _j in 0..fnpiv {
            let mut f = fj;
            fj = fj + fnr_curr;
            for _i in fnrows..fnrows_extended {
                f.set(0, 0.0);
                f = f + 1;
            }
        }
        let mut fi = work.fublock + fncols;
        for _i in 0..fnpiv {
            let mut f = fi;
            fi = fi + fnc_curr;
            for _j in fncols..fncols_extended {
                f.set(0, 0.0);
                f = f + 1;
            }
        }
    }
    let _ = &mut fl;
    work.fnrows = fnrows_extended;
    work.fncols = fncols_extended;
    true
}

/// `UMF_create_element`.
pub fn umf_create_element(numeric: &mut Numeric, work: &mut Work, symbolic: &Symbolic) -> bool {
    let mut row_degree = numeric.rperm_p();
    let mut row_tuples = numeric.uip_p();
    let mut row_tlen = numeric.uilen_p();
    let mut col_degree = numeric.cperm_p();
    let mut col_tuples = numeric.lip_p();
    let mut col_tlen = numeric.lilen_p();
    let n_row = work.n_row;
    let n_col = work.n_col;
    let nn = n_row.max(n_col);
    let fcols = work.fcols;
    let frows = work.frows;
    let mut fcpos = work.fcpos;
    let mut frpos = work.frpos;
    let mut m = numeric.mem();
    let fncols = work.fncols;
    let fnrows = work.fnrows;
    if !symbolic.fix_q {
        for j in 0..fncols {
            col_degree[fcols[j]] += fnrows;
        }
    }
    for i in 0..fnrows {
        row_degree[frows[i]] += fncols;
    }
    let mut e_arr = work.e;
    let max_mark = Int::MAX - (2 * nn + 1);
    if !work.pivcol_in_front {
        work.cdeg0 += nn + 1;
        if work.cdeg0 >= max_mark {
            work.cdeg0 = 1;
            for e in 1..=work.nel {
                if e_arr[e] != 0 {
                    (m + e_arr[e]).cast::<Element>().r(0).cdeg = 0;
                }
            }
        }
    }
    if !work.pivrow_in_front {
        work.rdeg0 += nn + 1;
        if work.rdeg0 >= max_mark {
            work.rdeg0 = 1;
            for e in 1..=work.nel {
                if e_arr[e] != 0 {
                    (m + e_arr[e]).cast::<Element>().r(0).rdeg = 0;
                }
            }
        }
    }
    if !work.pivrow_in_front {
        for j in 0..fncols {
            fcpos[fcols[j]] = EMPTY;
        }
    }
    if !work.pivcol_in_front {
        for i in 0..fnrows {
            frpos[frows[i]] = EMPTY;
        }
    }
    if fncols <= 0 || fnrows <= 0 {
        work.prior_element = EMPTY;
        return true;
    }
    let mut needunits: Int = 0;
    let mut got_memory = false;
    let (mut eloc, mut view) = umf_mem_alloc_element(numeric, fnrows, fncols, &mut needunits);
    let (r2, c2, do_fcpos) = if work.do_grow {
        (fnrows, fncols, false)
    } else {
        (
            fnrows.max(work.fnrows_new + 1),
            fncols.max(work.fncols_new + 1),
            work.pivrow_in_front,
        )
    };
    if eloc == 0 {
        if !umf_get_memory(numeric, work, needunits, r2, c2, do_fcpos) {
            return false;
        }
        got_memory = true;
        m = numeric.mem();
        let r = umf_mem_alloc_element(numeric, fnrows, fncols, &mut needunits);
        eloc = r.0;
        view = r.1;
        if eloc == 0 {
            return false;
        }
        // Rperm & co. are not reallocated by UMF_get_memory
        row_tuples = numeric.uip_p();
        row_tlen = numeric.uilen_p();
        col_tuples = numeric.lip_p();
        col_tlen = numeric.lilen_p();
        row_degree = numeric.rperm_p();
        col_degree = numeric.cperm_p();
    }
    let _ = (&mut row_degree, &mut col_degree);
    let v = view.unwrap();
    work.nel += 1;
    let e = work.nel;
    work.prior_element = e;
    e_arr[e] = eloc;
    if work.pivcol_in_front {
        v.ep.r(0).cdeg = work.cdeg0;
    }
    if work.pivrow_in_front {
        v.ep.r(0).rdeg = work.rdeg0;
    }
    let mut rows = v.rows;
    let mut cols = v.cols;
    for i in 0..fnrows {
        rows[i] = frows[i];
    }
    for i in 0..fncols {
        cols[i] = fcols[i];
    }
    let mut fcol = work.fcblock;
    let fnr_curr = work.fnr_curr;
    let mut c = v.c;
    for _j in 0..fncols {
        for i in 0..fnrows {
            c[i] = fcol[i];
        }
        fcol = fcol + fnr_curr;
        c = c + fnrows;
    }
    let tuple_e = e;
    if got_memory {
        for f in 0..fncols {
            let col = fcols[f];
            let tp: P<Tuple> = (m + col_tuples[col]).cast();
            let k = col_tlen[col];
            col_tlen[col] += 1;
            *tp.r(k) = Tuple { e: tuple_e, f };
        }
        for f in 0..fnrows {
            let row = frows[f];
            let tp: P<Tuple> = (m + row_tuples[row]).cast();
            let k = row_tlen[row];
            row_tlen[row] += 1;
            *tp.r(k) = Tuple { e: tuple_e, f };
        }
    } else {
        for which in 0..2 {
            let (count, list, mut tuples_p, mut tlen) = if which == 0 {
                (fncols, fcols, col_tuples, col_tlen)
            } else {
                (fnrows, frows, row_tuples, row_tlen)
            };
            let mut tp2: P<Tuple> = P::null();
            for f in 0..count {
                let x = list[f];
                let t1 = tuples_p[x];
                let mut size: Int = 0;
                let mut len: Int = 0;
                let mut tp: P<Tuple> = P::null();
                if t1 != 0 {
                    let p = m + t1;
                    tp = p.cast();
                    size = m.r(t1 - 1).size;
                    len = tlen[x];
                    tp2 = tp + len;
                }
                let mut needunits = units(SIZEOF_TUPLE, len + 1);
                let mut tp1: P<Tuple> = P::null();
                if needunits > size && t1 != 0 {
                    tp1 = tp;
                    tp2 = tp;
                    for t in 0..len {
                        let tt = *tp.r(t);
                        let ee = tt.e;
                        if e_arr[ee] == 0 {
                            continue;
                        }
                        let ff = tt.f;
                        let p = m + e_arr[ee];
                        let ep: P<Element> = p.cast();
                        let ecols: P<Int> = (p + UNITS_ELEMENT).cast();
                        let idx = if which == 0 {
                            ecols
                        } else {
                            ecols + ep.r(0).ncols
                        };
                        if idx[ff] == EMPTY {
                            continue;
                        }
                        *tp2.r(0) = tt;
                        tp2 = tp2 + 1;
                    }
                    len = tp2.diff(tp1);
                    tlen[x] = len;
                    needunits = units(SIZEOF_TUPLE, len + 1);
                }
                if needunits > size {
                    needunits = (2 * needunits).min(units(SIZEOF_TUPLE, nn));
                    let t2 = umf_mem_alloc_tail_block(numeric, needunits);
                    if t2 == 0 {
                        return umf_get_memory(numeric, work, 0, r2, c2, do_fcpos);
                    }
                    tuples_p[x] = t2;
                    tp2 = (m + t2).cast();
                    if t1 != 0 {
                        for i in 0..len {
                            *tp2.r(0) = *tp1.r(i);
                            tp2 = tp2 + 1;
                        }
                        umf_mem_free_tail_block(numeric, t1);
                    }
                }
                tlen[x] += 1;
                *tp2.r(0) = Tuple { e: tuple_e, f };
            }
        }
    }
    true
}
