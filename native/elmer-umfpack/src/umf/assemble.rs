//! `UMF_assemble` and `UMF_assemble_fixq` (the same source compiled with
//! and without `FIXQ`; here a flag): degree updates and assembly of the
//! prior elements into the current frontal matrix.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::numtypes::*;
use crate::amd::{flip, Int, EMPTY};
use crate::ptr::P;

fn row_assemble(row: Int, numeric: &mut Numeric, work: &Work, fixq: bool) {
    let mut col_degree = numeric.cperm_p();
    let row_tuples = numeric.uip_p();
    let tpi = row_tuples[row];
    if tpi == 0 {
        return;
    }
    let m = numeric.mem();
    let e_arr = work.e;
    let fcpos = work.fcpos;
    let frpos = work.frpos;
    let mut row_degree = numeric.rperm_p();
    let mut row_tlen = numeric.uilen_p();
    let rdeg0 = work.rdeg0;
    let fcblock = work.fcblock;
    let tp1: P<Tuple> = (m + tpi).cast();
    let mut tp2 = tp1;
    let ntp = row_tlen[row];
    for t in 0..ntp {
        let tt = *tp1.r(t);
        let e = tt.e;
        if e_arr[e] == 0 {
            continue;
        }
        let f = tt.f;
        let p = m + e_arr[e];
        let ep = p.cast::<Element>().r(0);
        let cols: P<Int> = (p + UNITS_ELEMENT).cast();
        let mut rows = cols + ep.ncols;
        if rows[f] == EMPTY {
            continue;
        }
        if ep.rdeg == rdeg0 {
            rows[f] = EMPTY;
            let nrows = ep.nrows;
            let ncols = ep.ncols;
            let mut s: P<Entry> =
                (p + UNITS_ELEMENT + super::util::units(4, ncols + nrows)).cast::<Entry>() + f;
            let ncolsleft = ep.ncolsleft;
            let mut frow = fcblock + frpos[row];
            row_degree[row] -= ncolsleft;
            if ncols == ncolsleft {
                for j in 0..ncols {
                    let col = cols[j];
                    if !fixq {
                        col_degree[col] -= 1;
                    }
                    frow[fcpos[col]] += s[0];
                    s = s + nrows;
                }
            } else {
                for j in 0..ncols {
                    let col = cols[j];
                    if col >= 0 {
                        if !fixq {
                            col_degree[col] -= 1;
                        }
                        frow[fcpos[col]] += s[0];
                    }
                    s = s + nrows;
                }
            }
            ep.nrowsleft -= 1;
        } else {
            *tp2.r(0) = tt;
            tp2 = tp2 + 1;
        }
    }
    row_tlen[row] = tp2.diff(tp1);
}

fn col_assemble(col: Int, numeric: &mut Numeric, work: &Work, fixq: bool) {
    let mut col_degree = numeric.cperm_p();
    let col_tuples = numeric.lip_p();
    let tpi = col_tuples[col];
    if tpi == 0 {
        return;
    }
    let m = numeric.mem();
    let e_arr = work.e;
    let fcpos = work.fcpos;
    let frpos = work.frpos;
    let mut row_degree = numeric.rperm_p();
    let mut col_tlen = numeric.lilen_p();
    let cdeg0 = work.cdeg0;
    let fcblock = work.fcblock;
    let tp1: P<Tuple> = (m + tpi).cast();
    let mut tp2 = tp1;
    let ntp = col_tlen[col];
    for t in 0..ntp {
        let tt = *tp1.r(t);
        let e = tt.e;
        if e_arr[e] == 0 {
            continue;
        }
        let f = tt.f;
        let p = m + e_arr[e];
        let ep = p.cast::<Element>().r(0);
        let mut cols: P<Int> = (p + UNITS_ELEMENT).cast();
        if cols[f] == EMPTY {
            continue;
        }
        if ep.cdeg == cdeg0 {
            cols[f] = EMPTY;
            let nrows = ep.nrows;
            let ncols = ep.ncols;
            let rows = cols + ncols;
            let s: P<Entry> = (p + UNITS_ELEMENT + super::util::units(4, ncols + nrows))
                .cast::<Entry>()
                + f * nrows;
            let nrowsleft = ep.nrowsleft;
            let mut fcol = fcblock + fcpos[col];
            if !fixq {
                col_degree[col] -= nrowsleft;
            }
            if nrows == nrowsleft {
                for i in 0..nrows {
                    let row = rows[i];
                    row_degree[row] -= 1;
                    fcol[frpos[row]] += s[i];
                }
            } else {
                for i in 0..nrows {
                    let row = rows[i];
                    if row >= 0 {
                        row_degree[row] -= 1;
                        fcol[frpos[row]] += s[i];
                    }
                }
            }
            ep.ncolsleft -= 1;
        } else {
            *tp2.r(0) = tt;
            tp2 = tp2 + 1;
        }
    }
    col_tlen[col] = tp2.diff(tp1);
}

/// `UMF_assemble` (`fixq == false`) / `UMF_assemble_fixq` (`fixq == true`).
pub fn umf_assemble(numeric: &mut Numeric, work: &mut Work, fixq: bool) {
    let mut col_degree = numeric.cperm_p();
    let fncols = work.fncols;
    let fnrows = work.fnrows;
    let fcpos = work.fcpos;
    let frpos = work.frpos;
    let mut row_degree = numeric.rperm_p();
    let row_tuples = numeric.uip_p();
    let mut row_tlen = numeric.uilen_p();
    let col_tuples = numeric.lip_p();
    let mut col_tlen = numeric.lilen_p();
    let mut e_arr = work.e;
    let m = numeric.mem();
    let mut wm = work.wm;
    let mut woo = work.woo;
    let rdeg0 = work.rdeg0;
    let cdeg0 = work.cdeg0;
    let fcblock = work.fcblock;
    numeric.maxnrows = numeric.maxnrows.max(fnrows);
    numeric.maxncols = numeric.maxncols.max(fncols);
    numeric.maxfrsize = numeric.maxfrsize.max(fnrows * fncols);

    let mut son_list: Int = 0;
    if !work.do_extend {
        let prior_uson = work.pivcol_in_front && !work.pivrow_in_front;
        let prior_lson = !work.pivcol_in_front && work.pivrow_in_front;
        if prior_uson || prior_lson {
            let e = work.prior_element;
            if e != EMPTY {
                (m + e_arr[e]).cast::<Element>().r(0).next = son_list;
                son_list = e;
            }
        }
    }
    work.prior_element = EMPTY;

    // scan the new rows
    for i2 in work.fscan_row..fnrows {
        let mut row = work.new_rows[i2];
        if row < 0 {
            row = flip(row);
        }
        let tpi = row_tuples[row];
        if tpi == 0 {
            continue;
        }
        let tp1: P<Tuple> = (m + tpi).cast();
        let mut tp2 = tp1;
        let ntp = row_tlen[row];
        for t in 0..ntp {
            let tt = *tp1.r(t);
            let e = tt.e;
            if e_arr[e] == 0 {
                continue;
            }
            let f = tt.f;
            let p = m + e_arr[e];
            let ep = p.cast::<Element>().r(0);
            let rows: P<Int> = (p + UNITS_ELEMENT).cast::<Int>() + ep.ncols;
            if rows[f] == EMPTY {
                continue;
            }
            if ep.cdeg < cdeg0 {
                ep.cdeg = ep.nrowsleft + cdeg0;
            }
            ep.cdeg -= 1;
            if ep.cdeg == cdeg0 && ep.next == EMPTY {
                ep.next = son_list;
                son_list = e;
            }
            *tp2.r(0) = tt;
            tp2 = tp2 + 1;
        }
        row_tlen[row] = tp2.diff(tp1);
    }

    // scan the new columns
    for j2 in work.fscan_col..fncols {
        let mut col = work.new_cols[j2];
        if col < 0 {
            col = flip(col);
        }
        let tpi = col_tuples[col];
        if tpi == 0 {
            continue;
        }
        let tp1: P<Tuple> = (m + tpi).cast();
        let mut tp2 = tp1;
        let ntp = col_tlen[col];
        for t in 0..ntp {
            let tt = *tp1.r(t);
            let e = tt.e;
            if e_arr[e] == 0 {
                continue;
            }
            let f = tt.f;
            let p = m + e_arr[e];
            let ep = p.cast::<Element>().r(0);
            let cols: P<Int> = (p + UNITS_ELEMENT).cast();
            if cols[f] == EMPTY {
                continue;
            }
            if ep.rdeg < rdeg0 {
                ep.rdeg = ep.ncolsleft + rdeg0;
            }
            ep.rdeg -= 1;
            if ep.rdeg == rdeg0 && ep.next == EMPTY {
                ep.next = son_list;
                son_list = e;
            }
            *tp2.r(0) = tt;
            tp2 = tp2 + 1;
        }
        col_tlen[col] = tp2.diff(tp1);
    }

    // assemble the sons
    let mut next: Int;
    let mut e = son_list;
    while e > 0 {
        let v = get_element(m + e_arr[e]);
        let (mut cols, mut rows, ncols, nrows) = (v.cols, v.rows, v.ncols, v.nrows);
        let mut s = v.c;
        let ep = v.ep.r(0);
        let nrowsleft = ep.nrowsleft;
        let ncolsleft = ep.ncolsleft;
        next = ep.next;
        ep.next = EMPTY;
        let extrdeg = if ep.rdeg < rdeg0 {
            ncolsleft
        } else {
            ep.rdeg - rdeg0
        };
        let extcdeg = if ep.cdeg < cdeg0 {
            nrowsleft
        } else {
            ep.cdeg - cdeg0
        };
        let ncols_to_assemble = ncolsleft - extrdeg;
        let nrows_to_assemble = nrowsleft - extcdeg;

        if extrdeg == 0 && extcdeg == 0 {
            // the whole element
            if nrows == nrowsleft {
                for i in 0..nrows {
                    let row = rows[i];
                    row_degree[row] -= ncolsleft;
                    wm[i] = frpos[row];
                }
                for j in 0..ncols {
                    let col = cols[j];
                    if ncols == ncolsleft || col >= 0 {
                        if !fixq {
                            col_degree[col] -= nrowsleft;
                        }
                        let mut fcol = fcblock + fcpos[col];
                        for i in 0..nrows {
                            fcol[wm[i]] += s[i];
                        }
                    }
                    s = s + nrows;
                }
            } else {
                let mut ngetrows = 0;
                for i in 0..nrows {
                    let row = rows[i];
                    if row >= 0 {
                        row_degree[row] -= ncolsleft;
                        woo[ngetrows] = i;
                        wm[ngetrows] = frpos[row];
                        ngetrows += 1;
                    }
                }
                for j in 0..ncols {
                    let col = cols[j];
                    if ncols == ncolsleft || col >= 0 {
                        if !fixq {
                            col_degree[col] -= nrowsleft;
                        }
                        let mut fcol = fcblock + fcpos[col];
                        for i in 0..nrowsleft {
                            fcol[wm[i]] += s[woo[i]];
                        }
                    }
                    s = s + nrows;
                }
            }
            let ee = e_arr[e];
            umf_mem_free_tail_block(numeric, ee);
            e_arr[e] = 0;
        } else if extcdeg == 0 {
            // Lson: assemble the columns in the front
            if ncols_to_assemble > 0 {
                let mut skip = false;
                if ncols_to_assemble * 16 < ncols && nrows == 1 {
                    skip = true;
                    work.any_skip = true;
                }
                if !skip {
                    if nrows == nrowsleft {
                        for i in 0..nrows {
                            let row = rows[i];
                            row_degree[row] -= ncols_to_assemble;
                            wm[i] = frpos[row];
                        }
                        for j in 0..ncols {
                            let col = cols[j];
                            if col >= 0 && fcpos[col] >= 0 {
                                if !fixq {
                                    col_degree[col] -= nrowsleft;
                                }
                                let mut fcol = fcblock + fcpos[col];
                                for i in 0..nrows {
                                    fcol[wm[i]] += s[i];
                                }
                                cols[j] = EMPTY;
                            }
                            s = s + nrows;
                        }
                    } else {
                        let mut ngetrows = 0;
                        for i in 0..nrows {
                            let row = rows[i];
                            if row >= 0 {
                                row_degree[row] -= ncols_to_assemble;
                                woo[ngetrows] = i;
                                wm[ngetrows] = frpos[row];
                                ngetrows += 1;
                            }
                        }
                        for j in 0..ncols {
                            let col = cols[j];
                            if col >= 0 && fcpos[col] >= 0 {
                                if !fixq {
                                    col_degree[col] -= nrowsleft;
                                }
                                let mut fcol = fcblock + fcpos[col];
                                for i in 0..nrowsleft {
                                    fcol[wm[i]] += s[woo[i]];
                                }
                                cols[j] = EMPTY;
                            }
                            s = s + nrows;
                        }
                    }
                    ep.ncolsleft = extrdeg;
                }
            }
        } else {
            // Uson: assemble the rows in the front
            if nrows_to_assemble > 0 {
                let mut skip = false;
                if nrows_to_assemble * 16 < nrows && ncols == 1 {
                    work.any_skip = true;
                    skip = true;
                }
                if !skip {
                    let mut ngetrows = 0;
                    for i in 0..nrows {
                        let row = rows[i];
                        if row >= 0 && frpos[row] >= 0 {
                            row_degree[row] -= ncolsleft;
                            woo[ngetrows] = i;
                            wm[ngetrows] = frpos[row];
                            ngetrows += 1;
                            rows[i] = EMPTY;
                        }
                    }
                    for j in 0..ncols {
                        let col = cols[j];
                        if ncols == ncolsleft || col >= 0 {
                            if !fixq {
                                col_degree[col] -= nrows_to_assemble;
                            }
                            let mut fcol = fcblock + fcpos[col];
                            for i in 0..nrows_to_assemble {
                                fcol[wm[i]] += s[woo[i]];
                            }
                        }
                        s = s + nrows;
                    }
                    ep.nrowsleft = extcdeg;
                }
            }
        }
        e = next;
    }

    if work.any_skip {
        row_assemble(work.pivrow, numeric, work, fixq);
    }
    if work.do_scan2row {
        for i2 in work.fscan_row..fnrows {
            let mut row = work.new_rows[i2];
            if row < 0 {
                row = flip(row);
            }
            if !(row == work.pivrow && work.any_skip) {
                row_assemble(row, numeric, work, fixq);
            }
        }
    }
    if work.any_skip {
        col_assemble(work.pivcol, numeric, work, fixq);
    }
    if work.do_scan2col {
        for j2 in work.fscan_col..fncols {
            let mut col = work.new_cols[j2];
            if col < 0 {
                col = flip(col);
            }
            if !(col == work.pivcol && work.any_skip) {
                col_assemble(col, numeric, work, fixq);
            }
        }
    }
}
