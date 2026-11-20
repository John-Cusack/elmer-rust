//! `UMF_get_memory` and `UMF_garbage_collection`.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::numtypes::*;
use super::util::*;
use crate::amd::{Int, EMPTY};
use crate::ptr::P;

/// `UMF_get_memory`: grow `Memory`, compact it, rebuild the tuples.
pub fn umf_get_memory(
    numeric: &mut Numeric,
    work: &mut Work,
    mut needunits: Int,
    r2: Int,
    c2: Int,
    do_fcpos: bool,
) -> bool {
    let n_row = work.n_row;
    let n_col = work.n_col;
    let row_degree = numeric.rperm_p();
    let col_degree = numeric.cperm_p();
    let mut row_tlen = numeric.uilen_p();
    let mut col_tlen = numeric.lilen_p();
    for row in 0..n_row {
        if row_degree[row] >= 0 {
            row_tlen[row] = 0;
        }
    }
    for col in 0..n_col {
        if col_degree[col] >= 0 {
            col_tlen[col] = 0;
        }
    }
    let mut nsize = needunits as f64 + 2.0;
    let mut tsize = 0.0;
    needunits += umf_tuple_lengths(numeric, work, &mut tsize);
    nsize += tsize;
    needunits += 2;
    let minsize = numeric.size + needunits;
    nsize += numeric.size as f64;
    let bsize = Int::MAX as f64 / SIZEOF_UNIT as f64 - 1.0;
    let mut newsize = (UMF_REALLOC_INCREASE * minsize as f64) as Int;
    nsize *= UMF_REALLOC_INCREASE;
    nsize += 1.0;
    if newsize < 0 || nsize > bsize {
        newsize = bsize as Int;
    } else {
        newsize = newsize.max(minsize);
    }
    newsize = newsize.max(numeric.size);
    numeric.ibig = EMPTY;

    // UMF_realloc (Memory, newsize): always succeeds here
    let old_ptr = numeric.memory.as_ptr();
    numeric
        .memory
        .resize(newsize.max(1) as usize, Unit::default());
    let costly = numeric.memory.as_ptr() != old_ptr;

    if work.e[0] != 0 {
        let nb = work.nb;
        let dr = work.fnr_curr;
        let dc = work.fnc_curr;
        work.flublock = (numeric.mem() + work.e[0]).cast();
        work.flblock = work.flublock + nb * nb;
        work.fublock = work.flblock + dr * nb;
        work.fcblock = work.fublock + nb * dc;
    }
    let newmem = newsize - numeric.size;
    if newmem >= 2 {
        let m = numeric.mem();
        let mut p = numeric.size - 2;
        m.r(p).size = newmem - 1;
        let i = numeric.size - 1;
        p += newmem;
        m.r(p).prevsize = newmem - 1;
        m.r(p).size = 1;
        numeric.size = newsize;
        umf_mem_free_tail_block(numeric, i);
        numeric.nrealloc += 1;
        if costly {
            numeric.ncostly += 1;
        }
    }
    umf_garbage_collection(numeric, work, r2, c2, do_fcpos);
    umf_build_tuples(numeric, work)
}

/// `UMF_garbage_collection`.
pub fn umf_garbage_collection(
    numeric: &mut Numeric,
    work: &mut Work,
    mut drnew: Int,
    dcnew: Int,
    do_fcpos: bool,
) {
    let col_degree = numeric.cperm_p();
    let row_degree = numeric.rperm_p();
    let mut row_tuples = numeric.uip_p();
    let mut col_tuples = numeric.lip_p();
    let mut e_arr = work.e;
    let n_row = work.n_row;
    let n_col = work.n_col;
    let m = numeric.mem();
    numeric.ngarbage += 1;
    for row in 0..n_row {
        if row_degree[row] >= 0 && row_tuples[row] != 0 {
            let p = m.r(row_tuples[row] - 1);
            p.size = -p.size;
            row_tuples[row] = 0;
        }
    }
    for col in 0..n_col {
        if col_degree[col] >= 0 && col_tuples[col] != 0 {
            let p = m.r(col_tuples[col] - 1);
            p.size = -p.size;
            col_tuples[col] = 0;
        }
    }
    let nel = work.nel;
    let mut e2: Int = 0;
    for e in 0..=nel {
        if e_arr[e] != 0 {
            let psrc = e_arr[e] - 1;
            if e > 0 {
                e2 += 1;
            }
            m.r(psrc).size = e2;
            e_arr[e] = 0;
            if e == work.prior_element {
                work.prior_element = e2;
            }
        }
    }
    work.nel = e2;
    let nel = work.nel;

    let mut psrc = numeric.size - 2;
    let mut pdest = psrc;
    let mut prevsize = m.r(psrc).prevsize;
    while prevsize > 0 {
        let size = prevsize;
        psrc -= size + 1;
        let e = m.r(psrc).size;
        prevsize = m.r(psrc).prevsize;
        if e == 0 {
            // the current frontal matrix
            let f1: P<Entry> = (m + (psrc + 1)).cast();
            let k = work.fnpiv;
            let dr = work.fnr_curr;
            let dc = work.fnc_curr;
            let r = work.fnrows;
            let c = work.fncols;
            let nb = work.nb;
            if drnew % 2 == 0 {
                drnew += 1;
            }
            drnew = dr.min(drnew);
            let pnext = pdest;
            let mut fsrc = work.flblock;
            let mut fdst = work.flblock;
            let mut gap1 = dr - r;
            let mut gap2 = drnew - r;
            for _j in 0..k {
                for _i in 0..r {
                    fdst.set(0, fsrc[0]);
                    fdst = fdst + 1;
                    fsrc = fsrc + 1;
                }
                fsrc = fsrc + gap1;
                fdst = fdst + gap2;
            }
            fdst = fdst + drnew * (nb - k);
            fsrc = work.fublock;
            gap1 = dc - c;
            gap2 = dcnew - c;
            for _i in 0..k {
                for _j in 0..c {
                    fdst.set(0, fsrc[0]);
                    fdst = fdst + 1;
                    fsrc = fsrc + 1;
                }
                fsrc = fsrc + gap1;
                fdst = fdst + gap2;
            }
            fdst = fdst + dcnew * (nb - k);
            fsrc = work.fcblock;
            gap1 = dr - r;
            gap2 = drnew - r;
            for _j in 0..c {
                for _i in 0..r {
                    fdst.set(0, fsrc[0]);
                    fdst = fdst + 1;
                    fsrc = fsrc + 1;
                }
                fsrc = fsrc + gap1;
                fdst = fdst + gap2;
            }
            if do_fcpos {
                let fcols = work.fcols;
                let mut fcpos = work.fcpos;
                for j in 0..c {
                    let col = fcols[j];
                    fcpos[col] = j * drnew;
                }
            }
            work.fnr_curr = drnew;
            work.fnc_curr = dcnew;
            work.fcurr_size = (drnew + nb) * (dcnew + nb);
            let mut size = units(SIZEOF_ENTRY, work.fcurr_size);
            size = size.max(1);
            m.r(pnext).prevsize = size;
            pdest -= size + 1;
            let f2: P<Entry> = (m + (pdest + 1)).cast();
            let off = nb * nb + drnew * nb + nb * dcnew + drnew * c;
            let mut fsrc = f1 + off;
            let mut fdst = f2 + off;
            let mut gap = drnew - r;
            for _j in (0..c).rev() {
                fsrc = fsrc - gap;
                fdst = fdst - gap;
                for _i in (0..r).rev() {
                    fdst = fdst - 1;
                    fsrc = fsrc - 1;
                    fdst.set(0, fsrc[0]);
                }
            }
            fsrc = fsrc - dcnew * (nb - k);
            fdst = fdst - dcnew * (nb - k);
            gap = dcnew - c;
            for _i in (0..k).rev() {
                fsrc = fsrc - gap;
                fdst = fdst - gap;
                for _j in (0..c).rev() {
                    fdst = fdst - 1;
                    fsrc = fsrc - 1;
                    fdst.set(0, fsrc[0]);
                }
            }
            fsrc = fsrc - drnew * (nb - k);
            fdst = fdst - drnew * (nb - k);
            gap = drnew - r;
            for _j in (0..k).rev() {
                fsrc = fsrc - gap;
                fdst = fdst - gap;
                for _i in (0..r).rev() {
                    fdst = fdst - 1;
                    fsrc = fsrc - 1;
                    fdst.set(0, fsrc[0]);
                }
            }
            fsrc = fsrc - nb * (nb - k);
            fdst = fdst - nb * (nb - k);
            gap = nb - k;
            for _j in (0..k).rev() {
                fsrc = fsrc - gap;
                fdst = fdst - gap;
                for _i in (0..k).rev() {
                    fdst = fdst - 1;
                    fsrc = fsrc - 1;
                    fdst.set(0, fsrc[0]);
                }
            }
            e_arr[0] = pdest + 1;
            work.flublock = (m + e_arr[0]).cast();
            work.flblock = work.flublock + nb * nb;
            work.fublock = work.flblock + drnew * nb;
            work.fcblock = work.fublock + nb * dcnew;
            m.r(pdest).prevsize = 0;
            m.r(pdest).size = size;
        } else if e > 0 {
            let v = get_element(m + (psrc + 1));
            let (cols, rows, ncols, nrows) = (v.cols, v.rows, v.ncols, v.nrows);
            let mut c = v.c;
            let es = *v.ep.r(0);
            let nrowsleft = es.nrowsleft;
            let ncolsleft = es.ncolsleft;
            let cdeg = es.cdeg;
            let rdeg = es.rdeg;
            let csize = nrowsleft * ncolsleft;
            let size2 = UNITS_ELEMENT
                + units(SIZEOF_INT, nrowsleft + ncolsleft)
                + units(SIZEOF_ENTRY, csize);
            let pnext = pdest;
            m.r(pnext).prevsize = size2;
            pdest -= size2 + 1;
            let p = m + (pdest + 1);
            let q = p + UNITS_ELEMENT;
            let mut cols2: P<Int> = q.cast();
            let mut rows2 = cols2 + ncolsleft;
            let q = q + units(SIZEOF_INT, nrowsleft + ncolsleft);
            let mut c2: P<Entry> = q.cast();
            if nrowsleft < nrows || ncolsleft < ncols {
                let mut c1 = c;
                let mut c3 = c;
                for j in 0..ncols {
                    if cols[j] >= 0 {
                        for i in 0..nrows {
                            if rows[i] >= 0 {
                                c3.set(0, c1[i]);
                                c3 = c3 + 1;
                            }
                        }
                    }
                    c1 = c1 + nrows;
                }
            }
            c = c + csize;
            c2 = c2 + csize;
            for _ in 0..csize {
                c2 = c2 - 1;
                c = c - 1;
                c2.set(0, c[0]);
            }
            let mut i2 = nrowsleft;
            for i in (0..nrows).rev() {
                if rows[i] >= 0 {
                    i2 -= 1;
                    rows2[i2] = rows[i];
                }
            }
            let mut j2 = ncolsleft;
            for j in (0..ncols).rev() {
                if cols[j] >= 0 {
                    j2 -= 1;
                    cols2[j2] = cols[j];
                }
            }
            e_arr[e] = pdest + 1;
            let epdest = (m + (pdest + 1)).cast::<Element>().r(0);
            epdest.next = EMPTY;
            epdest.ncols = ncolsleft;
            epdest.nrows = nrowsleft;
            epdest.ncolsleft = ncolsleft;
            epdest.nrowsleft = nrowsleft;
            epdest.rdeg = rdeg;
            epdest.cdeg = cdeg;
            m.r(pdest).prevsize = 0;
            m.r(pdest).size = size2;
        }
    }
    numeric.itail = pdest;
    m.r(pdest).prevsize = 0;
    numeric.ibig = EMPTY;
    numeric.tail_usage = numeric.size - numeric.itail;
    for e in nel + 1..work.elen {
        e_arr[e] = 0;
    }
}
