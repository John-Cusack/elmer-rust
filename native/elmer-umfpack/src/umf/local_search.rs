//! `UMF_local_search`: choose the next pivot (row and column) among the
//! candidate columns, in or out of the current front.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::numtypes::*;
use super::row_search::{umf_row_search, IN, OUT};
use super::symbolic::Symbolic;
use crate::amd::{flip, Int, EMPTY};
use crate::ptr::P;

const RELAX1: f64 = 0.25;
const SYM_RELAX1: f64 = 0.0;
const RELAX2: f64 = 0.1;
const RELAX3: f64 = 0.125;

const IN_IN: Int = 0;
const IN_OUT: Int = 1;
const OUT_IN: Int = 2;
const OUT_OUT: Int = 3;

fn remove_candidate(jj: Int, work: &mut Work, symbolic: &Symbolic) {
    if symbolic.fix_q {
        if work.ncand > 1 {
            work.candidates[0] = work.nextcand;
            work.nextcand += 1;
        } else {
            work.n_candidates = 0;
        }
    } else if work.ncand > MAX_CANDIDATES as Int {
        work.candidates[jj as usize] = work.nextcand;
        work.nextcand += 1;
    } else {
        work.candidates[jj as usize] = work.candidates[(work.ncand - 1) as usize];
        work.candidates[(work.ncand - 1) as usize] = EMPTY;
        work.n_candidates -= 1;
    }
    work.ncand -= 1;
}

#[inline]
fn unflip(i: Int) -> Int {
    if i < EMPTY {
        flip(i)
    } else {
        i
    }
}

/// Gather a candidate column from its element tuples into `(rows, pos, w)`.
#[allow(clippy::too_many_arguments)]
fn scan_col(
    numeric: &mut Numeric,
    e_arr: P<Int>,
    col: Int,
    mut pos: P<Int>,
    mut rows_out: P<Int>,
    mut w: P<Entry>,
    deg: &mut Int,
    max_cdeg: Int,
) -> bool {
    let m = numeric.mem();
    let col_tuples = numeric.lip_p();
    let mut col_tlen = numeric.lilen_p();
    let tpi = col_tuples[col];
    if tpi != 0 {
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
            let rows = cols + ep.ncols;
            let nrows = ep.nrows;
            let c: P<Entry> = (p + UNITS_ELEMENT + super::util::units(4, ep.ncols + nrows))
                .cast::<Entry>()
                + f * nrows;
            for i in 0..nrows {
                let row = rows[i];
                if row >= 0 {
                    let ps = pos[row];
                    if ps < 0 {
                        if *deg >= max_cdeg {
                            return false;
                        }
                        pos[row] = *deg;
                        rows_out[*deg] = row;
                        w[*deg] = c[i];
                        *deg += 1;
                    } else {
                        w[ps] += c[i];
                    }
                }
            }
            *tp2.r(0) = tt;
            tp2 = tp2 + 1;
        }
        col_tlen[col] = tp2.diff(tp1);
    }
    true
}

/// `UMF_local_search`.
pub fn umf_local_search(numeric: &mut Numeric, work: &mut Work, symbolic: &Symbolic) -> Int {
    let e_arr = work.e;
    let col_degree = numeric.cperm_p();
    let mut col_tuples = numeric.lip_p();
    let mut col_tlen = numeric.lilen_p();
    let mut wx = work.wx;
    let mut wy = work.wy;
    let mut wp = work.wp;
    let mut wm = work.wm;
    let mut woi = work.woi;
    let wio = work.wio;
    let mut woo = work.woo;
    let fcpos = work.fcpos;
    let mut frpos = work.frpos;
    let frows = work.frows;
    let fcols = work.fcols;
    let fnrows = work.fnrows;
    let fncols = work.fncols;
    let nb = work.nb;
    let fnr_curr = work.fnr_curr;
    let fnc_curr = work.fnc_curr;
    let fnpiv = work.fnpiv;
    let nothing = [EMPTY, EMPTY];
    let relax1 = if symbolic.prefer_diagonal {
        SYM_RELAX1
    } else {
        RELAX1
    };
    let mut fnzeros = work.fnzeros;
    let mut new_fnzeros = fnzeros;
    let mut jj = EMPTY;
    let fcblock = work.fcblock;
    let flblock = work.flblock;
    let fublock = work.fublock;
    let flublock = work.flublock;
    let max_cdeg = work.fnrows_max;

    if fnrows == 0 && fncols == 0 {
        work.firstsuper = work.ksuper;
    }
    let mut pivcol = [EMPTY, EMPTY];
    let mut cdeg_in = Int::MAX;
    let mut cdeg_out = Int::MAX;
    let mut pivrow = [[EMPTY, EMPTY], [EMPTY, EMPTY]];
    let mut rdeg = [[Int::MAX, Int::MAX], [Int::MAX, Int::MAX]];
    let mut freebie = [false, false];
    work.pivot_case = EMPTY;
    let mut bestcost = EMPTY;
    let mut nr_out = EMPTY;
    let mut nr_in = EMPTY;
    let mut jcand = [EMPTY, EMPTY];
    let mut fnrows_new = [[EMPTY, EMPTY], [EMPTY, EMPTY]];
    let mut fncols_new = [[EMPTY, EMPTY], [EMPTY, EMPTY]];

    // find the best candidate columns in and out of the front
    let col = work.candidates[0];
    let deg = if symbolic.fix_q {
        EMPTY
    } else {
        col_degree[col]
    };
    if fcpos[col] >= 0 {
        pivcol[IN] = col;
        cdeg_in = deg;
        jcand[IN] = 0;
    } else {
        pivcol[OUT] = col;
        cdeg_out = deg;
        jcand[OUT] = 0;
    }
    for j in 1..work.n_candidates {
        let col = work.candidates[j as usize];
        let deg = col_degree[col];
        if fcpos[col] >= 0 {
            if deg < cdeg_in || (deg == cdeg_in && col < pivcol[IN]) {
                pivcol[IN] = col;
                cdeg_in = deg;
                jcand[IN] = j;
            }
        } else if deg < cdeg_out || (deg == cdeg_out && col < pivcol[OUT]) {
            pivcol[OUT] = col;
            cdeg_out = deg;
            jcand[OUT] = j;
        }
    }
    cdeg_in = EMPTY;
    cdeg_out = EMPTY;

    // ---- pivot column in the front
    if pivcol[IN] != EMPTY {
        let fs = fcblock + fcpos[pivcol[IN]];
        let fu = fublock + fcpos[pivcol[IN]] / fnr_curr;
        let mut flu = flublock + fnpiv * nb;
        for i in 0..fnpiv {
            flu[i] = fu[i * fnc_curr];
        }
        if fnpiv > 1 {
            let mut flub = flublock;
            for j in 0..fnpiv {
                let fuj = flu[j];
                for i in j + 1..fnpiv {
                    flu[i] -= flub[i] * fuj;
                }
                flub = flub + nb;
            }
        }
        for i in 0..fnrows {
            wy[i] = fs[i];
        }
        for j in 0..fnpiv {
            let flub = flblock + j * fnr_curr;
            let fuj = flu[j];
            if fuj != 0.0 {
                for i in 0..fnrows {
                    wy[i] -= flub[i] * fuj;
                }
            }
        }
        cdeg_in = fnrows;
        if !scan_col(
            numeric,
            e_arr,
            pivcol[IN],
            frpos,
            frows,
            wy,
            &mut cdeg_in,
            max_cdeg,
        ) {
            return UMFPACK_ERROR_different_pattern;
        }
        nr_in = cdeg_in - fnrows;
        let status = umf_row_search(
            numeric,
            work,
            symbolic,
            fnrows,
            cdeg_in,
            frows,
            frpos,
            &mut pivrow[IN],
            &mut rdeg[IN],
            fcols,
            wio,
            &nothing,
            wy,
            pivcol[IN],
            &mut freebie,
        );
        if status == UMFPACK_ERROR_different_pattern {
            return UMFPACK_ERROR_different_pattern;
        }
        if pivrow[IN][IN] != EMPTY {
            let nc = rdeg[IN][IN] - fncols;
            let thiscost = (nr_in * (fncols - 1)) + (nc * (cdeg_in - 1));
            fnrows_new[IN][IN] = (fnrows - 1) + nr_in;
            fncols_new[IN][IN] = (fncols - 1) + nc;
            let do_extend = true;
            fnzeros = work.fnzeros + fnpiv * (nr_in + nc);
            let new_lusize = (fnpiv + 1) * (fnrows + nr_in + fncols + nc);
            let do_update = fnpiv > 0 && (fnzeros as f64 / new_lusize as f64) > RELAX2;
            work.pivot_case = IN_IN;
            bestcost = thiscost;
            work.do_extend = do_extend;
            work.do_update = do_update;
            new_fnzeros = fnzeros;
        }
        if pivrow[IN][OUT] != EMPTY {
            let mut nc = 0;
            for i in 0..rdeg[IN][OUT] {
                let col = wio[i];
                if fcpos[col] < 0 {
                    nc += 1;
                }
            }
            let mut thiscost = (nc * fnrows) + ((nr_in - 1) * (rdeg[IN][OUT] - 1));
            let extra_cols = ((fncols - 1) + nc) - (rdeg[IN][OUT] - 1);
            let extra_zeros = (nr_in - 1) * extra_cols;
            fnrows_new[IN][OUT] = fnrows + (nr_in - 1);
            fncols_new[IN][OUT] = (fncols - 1) + nc;
            let relaxed_front = fnrows_new[IN][OUT] * fncols_new[IN][OUT];
            let do_extend = if extra_zeros == 0 {
                true
            } else {
                (extra_zeros as f64) < relax1 * relaxed_front as f64
            };
            let do_update;
            if do_extend {
                thiscost += extra_zeros;
                fnzeros = work.fnzeros + fnpiv * (nr_in + nc);
                let new_lusize = (fnpiv + 1) * (fnrows + nr_in + fncols + nc);
                do_update = fnpiv > 0 && (fnzeros as f64 / new_lusize as f64) > RELAX3;
            } else {
                do_update = fnpiv > 0;
                fnzeros = 0;
                fnrows_new[IN][OUT] = cdeg_in - 1;
                fncols_new[IN][OUT] = rdeg[IN][OUT] - 1;
            }
            if bestcost == EMPTY || thiscost < bestcost {
                work.pivot_case = IN_OUT;
                bestcost = thiscost;
                work.do_extend = do_extend;
                work.do_update = do_update;
                new_fnzeros = fnzeros;
            }
        }
    }

    // ---- pivot column out of the front
    let mut search_pivcol_out = bestcost != 0 && pivcol[OUT] != EMPTY;
    if symbolic.prefer_diagonal {
        search_pivcol_out = search_pivcol_out && pivrow[IN][IN] == EMPTY;
    }
    if search_pivcol_out {
        cdeg_out = 0;
        if !scan_col(
            numeric,
            e_arr,
            pivcol[OUT],
            wp,
            wm,
            wx,
            &mut cdeg_out,
            max_cdeg,
        ) {
            return UMFPACK_ERROR_different_pattern;
        }
        let prior = pivrow[IN];
        let status = umf_row_search(
            numeric,
            work,
            symbolic,
            0,
            cdeg_out,
            wm,
            wp,
            &mut pivrow[OUT],
            &mut rdeg[OUT],
            woi,
            woo,
            &prior,
            wx,
            pivcol[OUT],
            &mut freebie,
        );
        if status == UMFPACK_ERROR_different_pattern {
            return UMFPACK_ERROR_different_pattern;
        }
        for i in 0..cdeg_out {
            wp[wm[i]] = EMPTY;
        }
        if status == UMFPACK_WARNING_singular_matrix {
            remove_candidate(jcand[OUT], work, symbolic);
            work.ndiscard += 1;
            col_tlen[pivcol[OUT]] = 0;
            let t = col_tuples[pivcol[OUT]];
            umf_mem_free_tail_block(numeric, t);
            col_tuples[pivcol[OUT]] = 0;
            return UMFPACK_WARNING_singular_matrix;
        }
        if freebie[IN] {
            woi = fcols;
            rdeg[OUT][IN] = rdeg[IN][IN];
        }
        if freebie[OUT] {
            woo = wio;
            rdeg[OUT][OUT] = rdeg[IN][OUT];
        }
        if pivrow[OUT][IN] != EMPTY {
            let did_rowmerge = cdeg_out == 0;
            if did_rowmerge {
                wm[0] = pivrow[OUT][IN];
                wx[0] = 0.0;
                cdeg_out = 1;
            }
            let nc = rdeg[OUT][IN] - fncols;
            nr_out = 0;
            for i in 0..cdeg_out {
                let row = wm[i];
                if frpos[row] < 0 || frpos[row] >= fnrows {
                    nr_out += 1;
                }
            }
            let mut thiscost = (nr_out * fncols) + ((nc - 1) * (cdeg_out - 1));
            let extra_rows = ((fnrows - 1) + nr_out) - (cdeg_out - 1);
            let extra_zeros = (nc - 1) * extra_rows;
            fnrows_new[OUT][IN] = (fnrows - 1) + nr_out;
            fncols_new[OUT][IN] = fncols + (nc - 1);
            let relaxed_front = fnrows_new[OUT][IN] * fncols_new[OUT][IN];
            let do_extend = if did_rowmerge {
                false
            } else if extra_zeros == 0 {
                true
            } else {
                (extra_zeros as f64) < relax1 * relaxed_front as f64
            };
            let do_update;
            if do_extend {
                thiscost += extra_zeros;
                fnzeros = work.fnzeros + fnpiv * (nr_out + nc);
                let new_lusize = (fnpiv + 1) * (fnrows + nr_out + fncols + nc);
                do_update = fnpiv > 0 && (fnzeros as f64 / new_lusize as f64) > RELAX3;
            } else {
                do_update = fnpiv > 0;
                fnzeros = 0;
                fnrows_new[OUT][IN] = cdeg_out - 1;
                fncols_new[OUT][IN] = rdeg[OUT][IN] - 1;
            }
            if bestcost == EMPTY || thiscost < bestcost {
                work.pivot_case = OUT_IN;
                bestcost = thiscost;
                work.do_extend = do_extend;
                work.do_update = do_update;
                new_fnzeros = fnzeros;
            }
        }
        if pivrow[OUT][OUT] != EMPTY {
            let did_rowmerge = cdeg_out == 0;
            if did_rowmerge {
                wm[0] = pivrow[OUT][OUT];
                wx[0] = 0.0;
                cdeg_out = 1;
                nr_out = 1;
            }
            let nc;
            let extra_zeros;
            let mut thiscost;
            let relaxed_front;
            if fnrows == 0 && fncols == 0 {
                nc = rdeg[OUT][OUT];
                nr_out = cdeg_out;
                extra_zeros = 0;
                thiscost = (nc - 1) * (cdeg_out - 1);
                fnrows_new[OUT][OUT] = nr_out - 1;
                fncols_new[OUT][OUT] = nc - 1;
                relaxed_front = fnrows_new[OUT][OUT] * fncols_new[OUT][OUT];
            } else {
                if nr_out == EMPTY {
                    nr_out = 0;
                    for i in 0..cdeg_out {
                        let row = wm[i];
                        if frpos[row] < 0 || frpos[row] >= fnrows {
                            nr_out += 1;
                        }
                    }
                }
                let mut ncc = 0;
                for i in 0..rdeg[OUT][OUT] {
                    let col = woo[i];
                    if fcpos[col] < 0 {
                        ncc += 1;
                    }
                }
                nc = ncc;
                let extra_cols = (fncols + (nc - 1)) - (rdeg[OUT][OUT] - 1);
                let extra_rows = (fnrows + (nr_out - 1)) - (cdeg_out - 1);
                extra_zeros = ((nc - 1) * extra_rows) + ((nr_out - 1) * extra_cols);
                thiscost = ((nc - 1) * (cdeg_out - 1)) + ((nr_out - 1) * (fncols - extra_cols));
                fnrows_new[OUT][OUT] = fnrows + (nr_out - 1);
                fncols_new[OUT][OUT] = fncols + (nc - 1);
                relaxed_front = fnrows_new[OUT][OUT] * fncols_new[OUT][OUT];
            }
            let do_extend = if did_rowmerge {
                false
            } else if extra_zeros == 0 {
                true
            } else {
                (extra_zeros as f64) < relax1 * relaxed_front as f64
            };
            let do_update;
            if do_extend {
                thiscost += extra_zeros;
                fnzeros = work.fnzeros + fnpiv * (nr_out + nc);
                let new_lusize = (fnpiv + 1) * (fnrows + nr_out + fncols + nc);
                do_update = fnpiv > 0 && (fnzeros as f64 / new_lusize as f64) > RELAX3;
            } else {
                do_update = fnpiv > 0;
                fnzeros = 0;
                fnrows_new[OUT][OUT] = cdeg_out - 1;
                fncols_new[OUT][OUT] = rdeg[OUT][OUT] - 1;
            }
            if bestcost == EMPTY || thiscost < bestcost {
                work.pivot_case = OUT_OUT;
                bestcost = thiscost;
                work.do_extend = do_extend;
                work.do_update = do_update;
                new_fnzeros = fnzeros;
            }
        }
    }
    let _ = bestcost;

    work.fnzeros = new_fnzeros;
    match work.pivot_case {
        IN_IN => {
            work.pivcol_in_front = true;
            work.pivrow_in_front = true;
            work.pivcol = pivcol[IN];
            work.pivrow = pivrow[IN][IN];
            work.ccdeg = nr_in;
            work.wrow = fcols;
            work.rrdeg = rdeg[IN][IN];
            jj = jcand[IN];
            work.fnrows_new = fnrows_new[IN][IN];
            work.fncols_new = fncols_new[IN][IN];
        }
        IN_OUT => {
            work.pivcol_in_front = true;
            work.pivrow_in_front = false;
            work.pivcol = pivcol[IN];
            work.pivrow = pivrow[IN][OUT];
            work.ccdeg = nr_in;
            work.wrow = wio;
            work.rrdeg = rdeg[IN][OUT];
            jj = jcand[IN];
            work.fnrows_new = fnrows_new[IN][OUT];
            work.fncols_new = fncols_new[IN][OUT];
        }
        OUT_IN => {
            work.pivcol_in_front = false;
            work.pivrow_in_front = true;
            work.pivcol = pivcol[OUT];
            work.pivrow = pivrow[OUT][IN];
            work.ccdeg = cdeg_out;
            work.wrow = woi;
            work.rrdeg = rdeg[OUT][IN];
            jj = jcand[OUT];
            work.fnrows_new = fnrows_new[OUT][IN];
            work.fncols_new = fncols_new[OUT][IN];
        }
        OUT_OUT => {
            work.pivcol_in_front = false;
            work.pivrow_in_front = false;
            work.pivcol = pivcol[OUT];
            work.pivrow = pivrow[OUT][OUT];
            work.ccdeg = cdeg_out;
            work.wrow = woo;
            work.rrdeg = rdeg[OUT][OUT];
            jj = jcand[OUT];
            work.fnrows_new = fnrows_new[OUT][OUT];
            work.fncols_new = fncols_new[OUT][OUT];
        }
        _ => {}
    }

    if !work.pivcol_in_front && pivcol[IN] != EMPTY {
        for i in fnrows..cdeg_in {
            frpos[frows[i]] = EMPTY;
        }
    }
    remove_candidate(jj, work, symbolic);

    work.do_grow = work.fnrows_new + 1 > fnr_curr || work.fncols_new + 1 > fnc_curr;
    if work.do_grow && !work.do_update && fnpiv > 0 {
        work.nforced += 1;
        work.do_update = true;
    }
    if work.do_extend {
        work.do_scan2row = fncols > 0;
        work.do_scan2col = fnrows > 0;
    } else {
        work.do_scan2row = fncols > 0 && work.pivrow_in_front;
        work.do_scan2col = fnrows > 0 && work.pivcol_in_front;
    }

    if symbolic.prefer_diagonal && work.pivcol < work.n_col - symbolic.nempty_col {
        let mut diagonal_map = work.diagonal_map;
        let mut diagonal_imap = work.diagonal_imap;
        let mut row2 = diagonal_map[work.pivcol];
        let col2 = diagonal_imap[work.pivrow];
        if row2 < 0 {
            work.noff_diagonal += 1;
            row2 = unflip(row2);
        }
        if row2 != work.pivrow {
            diagonal_map[work.pivcol] = flip(work.pivrow);
            diagonal_imap[work.pivrow] = work.pivcol;
            diagonal_map[col2] = flip(row2);
            diagonal_imap[row2] = col2;
        }
    }
    UMFPACK_OK
}
