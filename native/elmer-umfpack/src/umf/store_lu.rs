//! `UMF_store_lu` (without `DROP`; Elmer's default `droptol` is 0): copy the
//! pivot rows and columns of the front into the LU factors, as L and U
//! chains in the head of `Memory`.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::gc::umf_get_memory;
use super::numtypes::*;
use super::util::*;
use crate::amd::{Int, EMPTY};
use crate::ptr::P;

/// `UMF_store_lu`.
pub fn umf_store_lu(numeric: &mut Numeric, work: &mut Work) -> bool {
    let fnrows = work.fnrows;
    let fncols = work.fncols;
    let fnpiv = work.fnpiv;
    let mut lpos = numeric.lpos_p();
    let mut upos = numeric.upos_p();
    let mut lilen = numeric.lilen_p();
    let mut uilen = numeric.uilen_p();
    let mut lip = numeric.lip_p();
    let mut uip = numeric.uip_p();
    let mut d = numeric.d_p();
    let mut flublock = work.flublock;
    let mut flblock = work.flblock;
    let mut fublock = work.fublock;
    let mut fnr_curr = work.fnr_curr;
    let mut fnc_curr = work.fnc_curr;
    let frows = work.frows;
    let fcols = work.fcols;
    let mut lpattern = work.lpattern;
    let mut llen = work.llen;
    let mut upattern = work.upattern;
    let mut ulen = work.ulen;
    let nb = work.nb;
    let pivrow_arr = work.pivrow_arr;
    let pivcol_arr = work.pivcol_arr;
    let mut lnz_arr = [0 as Int; MAXNB as usize];

    for kk in 0..fnpiv {
        let k = work.npiv + kk;
        let pivrow = pivrow_arr[kk as usize];
        let pivcol = pivcol_arr[kk as usize];
        let mut pivrow_position = lpos[pivrow];
        if pivrow_position != EMPTY {
            llen -= 1;
            let row = lpattern[llen];
            lpattern[pivrow_position] = row;
            lpos[row] = pivrow_position;
            lpos[pivrow] = EMPTY;
        }
        let mut fl1 = flublock + kk * nb;
        let mut fl2 = flblock + kk * fnr_curr;
        let pivot_value = fl1[kk];
        d[k] = pivot_value;
        let zero_pivot = pivot_value == 0.0;
        let mut lnz: Int = 0;
        let mut lnz2i: Int = 0;
        let mut lnz2x = llen;
        for i in kk + 1..fnpiv {
            if fl1[i] == 0.0 {
                continue;
            }
            lnz += 1;
            if lpos[pivrow_arr[i as usize]] == EMPTY {
                lnz2i += 1;
            }
        }
        for i in 0..fnrows {
            if fl2[i] == 0.0 {
                continue;
            }
            lnz += 1;
            if lpos[frows[i]] == EMPTY {
                lnz2i += 1;
            }
        }
        lnz2x += lnz2i;
        let new_lchain = if llen == 0 || zero_pivot {
            true
        } else {
            units(SIZEOF_ENTRY, lnz) + units(SIZEOF_INT, lnz)
                <= units(SIZEOF_ENTRY, lnz2x) + units(SIZEOF_INT, lnz2i)
        };
        let (lnzi, lnzx);
        if new_lchain {
            pivrow_position = EMPTY;
            for i in 0..llen {
                let row = lpattern[i];
                lpos[row] = EMPTY;
            }
            llen = 0;
            lnzi = lnz;
            lnzx = lnz;
        } else {
            lnzi = lnz2i;
            lnzx = lnz2x;
        }
        let size = units(SIZEOF_INT, lnzi) + units(SIZEOF_ENTRY, lnzx);
        let mut p = umf_mem_alloc_head_block(numeric, size);
        if p == 0 {
            let (r2, c2) = if work.do_grow {
                (fnrows, fncols)
            } else {
                (
                    fnrows.max(work.fnrows_new + 1),
                    fncols.max(work.fncols_new + 1),
                )
            };
            if !umf_get_memory(numeric, work, size, r2, c2, true) {
                return false;
            }
            p = umf_mem_alloc_head_block(numeric, size);
            if p == 0 {
                return false;
            }
            fnc_curr = work.fnc_curr;
            fnr_curr = work.fnr_curr;
            flublock = work.flublock;
            flblock = work.flblock;
            fublock = work.fublock;
            fl1 = flublock + kk * nb;
            fl2 = flblock + kk * fnr_curr;
        }
        let m = numeric.mem();
        let mut lipv = p;
        let mut li: P<Int> = (m + p).cast();
        p += units(SIZEOF_INT, lnzi);
        let mut lval: P<Entry> = (m + p).cast();
        for i in 0..lnzx {
            lval[i] = 0.0;
        }
        if new_lchain {
            lipv = -lipv;
            for i in kk + 1..fnpiv {
                let x = fl1[i];
                if x == 0.0 {
                    continue;
                }
                let row2 = pivrow_arr[i as usize];
                let pos = llen;
                llen += 1;
                lpattern[pos] = row2;
                lpos[row2] = pos;
                li[pos] = row2;
                lval[pos] = x;
            }
            for i in 0..fnrows {
                let x = fl2[i];
                if x == 0.0 {
                    continue;
                }
                let row2 = frows[i];
                let pos = llen;
                llen += 1;
                lpattern[pos] = row2;
                lpos[row2] = pos;
                li[pos] = row2;
                lval[pos] = x;
            }
        } else {
            for i in kk + 1..fnpiv {
                let x = fl1[i];
                if x == 0.0 {
                    continue;
                }
                let row2 = pivrow_arr[i as usize];
                let mut pos = lpos[row2];
                if pos == EMPTY {
                    pos = llen;
                    llen += 1;
                    lpattern[pos] = row2;
                    lpos[row2] = pos;
                    li[0] = row2;
                    li = li + 1;
                }
                lval[pos] = x;
            }
            for i in 0..fnrows {
                let x = fl2[i];
                if x == 0.0 {
                    continue;
                }
                let row2 = frows[i];
                let mut pos = lpos[row2];
                if pos == EMPTY {
                    pos = llen;
                    llen += 1;
                    lpattern[pos] = row2;
                    lpos[row2] = pos;
                    li[0] = row2;
                    li = li + 1;
                }
                lval[pos] = x;
            }
        }
        numeric.lnz += lnz;
        numeric.all_lnz += lnz;
        lnz_arr[kk as usize] = lnz;
        numeric.n_lentries += lnzx;
        work.llen = llen;
        numeric.isize += lnzi;
        lpos[pivrow] = pivrow_position;
        lip[pivcol] = lipv;
        lilen[pivcol] = lnzi;
    }

    for kk in 0..fnpiv {
        let k = work.npiv + kk;
        let pivrow = pivrow_arr[kk as usize];
        let pivcol = pivcol_arr[kk as usize];
        let zero_pivot = d[k] == 0.0;
        let mut fu1 = flublock + kk;
        let mut fu2 = fublock + kk * fnc_curr;
        let mut unz: Int = 0;
        let mut unz2i: Int = 0;
        let mut unz2x = ulen;
        let mut pivcol_position = upos[pivcol];
        if pivcol_position != EMPTY {
            unz2x -= 1;
        }
        for i in kk + 1..fnpiv {
            if fu1[i * nb] == 0.0 {
                continue;
            }
            unz += 1;
            if upos[pivcol_arr[i as usize]] == EMPTY {
                unz2i += 1;
            }
        }
        for i in 0..fncols {
            if fu2[i] == 0.0 {
                continue;
            }
            unz += 1;
            if upos[fcols[i]] == EMPTY {
                unz2i += 1;
            }
        }
        unz2x += unz2i;
        let new_uchain = if ulen == 0 || zero_pivot {
            true
        } else {
            units(SIZEOF_ENTRY, unz) + units(SIZEOF_INT, unz)
                <= units(SIZEOF_ENTRY, unz2x) + units(SIZEOF_INT, unz2i)
        };
        let mut size: Int = 0;
        let unzx;
        if new_uchain {
            size += units(SIZEOF_INT, ulen);
            unzx = unz;
        } else {
            unzx = unz2x;
        }
        size += units(SIZEOF_ENTRY, unzx);
        let mut p = umf_mem_alloc_head_block(numeric, size);
        if p == 0 {
            let (r2, c2) = if work.do_grow {
                (fnrows, fncols)
            } else {
                (
                    fnrows.max(work.fnrows_new + 1),
                    fncols.max(work.fncols_new + 1),
                )
            };
            if !umf_get_memory(numeric, work, size, r2, c2, true) {
                return false;
            }
            p = umf_mem_alloc_head_block(numeric, size);
            if p == 0 {
                return false;
            }
            fnc_curr = work.fnc_curr;
            fnr_curr = work.fnr_curr;
            flublock = work.flublock;
            flblock = work.flblock;
            fublock = work.fublock;
            fu1 = flublock + kk;
            fu2 = fublock + kk * fnc_curr;
        }
        let _ = (fnr_curr, flblock);
        let m = numeric.mem();
        let mut uipv = p;
        let uilenv;
        if new_uchain {
            uipv = -uipv;
            pivcol_position = EMPTY;
            uilenv = ulen;
            let mut ui: P<Int> = (m + p).cast();
            numeric.isize += ulen;
            p += units(SIZEOF_INT, ulen);
            for i in 0..ulen {
                let col = upattern[i];
                upos[col] = EMPTY;
                ui[i] = col;
            }
            ulen = 0;
        } else {
            if pivcol_position != EMPTY {
                ulen -= 1;
                let col = upattern[ulen];
                upattern[pivcol_position] = col;
                upos[col] = pivcol_position;
                upos[pivcol] = EMPTY;
            }
            uilenv = unz2i;
        }
        let mut uval: P<Entry> = (m + p).cast();
        for i in 0..unzx {
            uval[i] = 0.0;
        }
        if new_uchain {
            for i in kk + 1..fnpiv {
                let x = fu1[i * nb];
                if x == 0.0 {
                    continue;
                }
                let col2 = pivcol_arr[i as usize];
                let pos = ulen;
                ulen += 1;
                upattern[pos] = col2;
                upos[col2] = pos;
                uval[pos] = x;
            }
            for i in 0..fncols {
                let x = fu2[i];
                if x == 0.0 {
                    continue;
                }
                let col2 = fcols[i];
                let pos = ulen;
                ulen += 1;
                upattern[pos] = col2;
                upos[col2] = pos;
                uval[pos] = x;
            }
        } else {
            for i in kk + 1..fnpiv {
                let x = fu1[i * nb];
                if x == 0.0 {
                    continue;
                }
                let col2 = pivcol_arr[i as usize];
                let mut pos = upos[col2];
                if pos == EMPTY {
                    pos = ulen;
                    ulen += 1;
                    upattern[pos] = col2;
                    upos[col2] = pos;
                }
                uval[pos] = x;
            }
            for i in 0..fncols {
                let x = fu2[i];
                if x == 0.0 {
                    continue;
                }
                let col2 = fcols[i];
                let mut pos = upos[col2];
                if pos == EMPTY {
                    pos = ulen;
                    ulen += 1;
                    upattern[pos] = col2;
                    upos[col2] = pos;
                }
                uval[pos] = x;
            }
        }
        numeric.unz += unz;
        numeric.all_unz += unz;
        let lk = lnz_arr[kk as usize];
        numeric.flops += DIV_FLOPS * lk as f64 + MULTSUB_FLOPS * (lk.wrapping_mul(unz)) as f64;
        numeric.n_uentries += unzx;
        work.ulen = ulen;
        upos[pivcol] = pivcol_position;
        uip[pivrow] = uipv;
        uilen[pivrow] = uilenv;
    }
    work.npiv += fnpiv;
    work.fnpiv = 0;
    work.fnzeros = 0;
    true
}
