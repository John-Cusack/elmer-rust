//! `UMF_row_search`: choose the pivot row(s) for a candidate pivot column
//! and build their patterns.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::dense::sabs;
use super::numtypes::*;
use super::symbolic::Symbolic;
use super::util::cmax;
use crate::amd::{flip, Int, EMPTY};
use crate::ptr::P;

pub const IN: usize = 0;
pub const OUT: usize = 1;

/// `UMF_row_search`.
#[allow(clippy::too_many_arguments)]
pub fn umf_row_search(
    numeric: &mut Numeric,
    work: &mut Work,
    symbolic: &Symbolic,
    cdeg0: Int,
    cdeg1: Int,
    pattern: P<Int>,
    pos: P<Int>,
    pivrow: &mut [Int; 2],
    rdeg: &mut [Int; 2],
    mut w_i: P<Int>,
    mut w_o: P<Int>,
    prior_pivrow: &[Int; 2],
    wxy: P<Entry>,
    pivcol: Int,
    freebie: &mut [bool; 2],
) -> Int {
    let mut pivot = [0.0f64; 2];
    let row_degree = numeric.rperm_p();
    let row_tuples = numeric.uip_p();
    let mut row_tlen = numeric.uilen_p();
    let mut wrp = work.wrp;
    let frpos = work.frpos;
    let e_arr = work.e;
    let m = numeric.mem();
    let fnrows = work.fnrows;
    let prefer_diagonal = symbolic.prefer_diagonal;
    let diagonal_map = work.diagonal_map;
    let mut diag_row;
    let was_offdiag;
    if !diagonal_map.is_null() {
        diag_row = diagonal_map[pivcol];
        was_offdiag = diag_row < 0;
        if was_offdiag {
            diag_row = flip(diag_row);
        }
    } else {
        diag_row = EMPTY;
        was_offdiag = true; // EMPTY (nonzero) in C
    }
    let max_rdeg = work.fncols_max;

    let mut maxval = 0.0f64;
    let mut nans_in_col = false;
    for i in 0..cdeg1 {
        let value = sabs(wxy[i]);
        if value.is_nan() {
            nans_in_col = true;
            maxval = value;
            break;
        }
        maxval = cmax(maxval, value);
    }
    let toler = numeric.relpt * maxval;
    let mut toler2 = numeric.relpt2 * maxval;
    toler2 = if was_offdiag { toler } else { toler2 };
    if toler.is_nan() || toler2.is_nan() {
        nans_in_col = true;
    }

    if !nans_in_col {
        let mut found = false;
        if prefer_diagonal {
            let i = pos[diag_row];
            if i >= 0 {
                let a = sabs(wxy[i]);
                if a != 0.0 && a >= toler2 {
                    found = true;
                    if frpos[diag_row] >= 0 && frpos[diag_row] < fnrows {
                        pivrow[IN] = diag_row;
                        pivrow[OUT] = EMPTY;
                    } else {
                        pivrow[IN] = EMPTY;
                        pivrow[OUT] = diag_row;
                    }
                }
            }
        }
        if !found {
            if cdeg0 > 0 {
                let mut i = 0;
                while i < cdeg0 {
                    let a = sabs(wxy[i]);
                    if a != 0.0 && a >= toler {
                        let row = pattern[i];
                        let deg = row_degree[row];
                        if deg < rdeg[IN] || (deg == rdeg[IN] && a > pivot[IN]) {
                            pivrow[IN] = row;
                            rdeg[IN] = deg;
                            pivot[IN] = a;
                        }
                    }
                    i += 1;
                }
                while i < cdeg1 {
                    let a = sabs(wxy[i]);
                    if a != 0.0 && a >= toler {
                        let row = pattern[i];
                        let deg = row_degree[row];
                        if deg < rdeg[OUT] || (deg == rdeg[OUT] && a > pivot[OUT]) {
                            pivrow[OUT] = row;
                            rdeg[OUT] = deg;
                            pivot[OUT] = a;
                        }
                    }
                    i += 1;
                }
            } else {
                for i in 0..cdeg1 {
                    let a = sabs(wxy[i]);
                    if a != 0.0 && a >= toler {
                        let row = pattern[i];
                        let deg = row_degree[row];
                        if frpos[row] >= 0 && frpos[row] < fnrows {
                            if deg < rdeg[IN] || (deg == rdeg[IN] && a > pivot[IN]) {
                                pivrow[IN] = row;
                                rdeg[IN] = deg;
                                pivot[IN] = a;
                            }
                        } else if deg < rdeg[OUT] || (deg == rdeg[OUT] && a > pivot[OUT]) {
                            pivrow[OUT] = row;
                            rdeg[OUT] = deg;
                            pivot[OUT] = a;
                        }
                    }
                }
            }
        }
    }

    if cdeg1 > 0 && pivrow[IN] == EMPTY && pivrow[OUT] == EMPTY {
        let row = pattern[0];
        let deg = row_degree[row];
        if frpos[row] >= 0 && frpos[row] < fnrows {
            pivrow[IN] = row;
            rdeg[IN] = deg;
        } else {
            pivrow[OUT] = row;
            rdeg[OUT] = deg;
        }
        for i in 1..cdeg1 {
            let row = pattern[i];
            let deg = row_degree[row];
            if frpos[row] >= 0 && frpos[row] < fnrows {
                if deg < rdeg[IN] || (deg == rdeg[IN] && row == diag_row) {
                    pivrow[IN] = row;
                    rdeg[IN] = deg;
                }
            } else if deg < rdeg[OUT] || (deg == rdeg[OUT] && row == diag_row) {
                pivrow[OUT] = row;
                rdeg[OUT] = deg;
            }
        }
    }

    if cdeg1 == 0 {
        if fnrows > 0 {
            pivrow[IN] = work.frows[0];
        } else {
            let front_leftmostdesc = &symbolic.front_leftmostdesc;
            let front_1strow = &symbolic.front_1strow;
            let mut front_new1strow = work.front_new1strow;
            let nfr = symbolic.nfr;
            let n_row = numeric.n_row;
            let frontid = work.frontid;
            let fleftmost = front_leftmostdesc[frontid as usize];
            let row1 = front_new1strow[fleftmost];
            let row2 = front_1strow[frontid as usize + 1] - 1;
            let mut row = row1;
            while row <= row2 {
                if row_degree[row] >= 0 {
                    pivrow[OUT] = row;
                    break;
                }
                row += 1;
            }
            front_new1strow[fleftmost] = row;
            if pivrow[OUT] == EMPTY {
                let row1 = front_new1strow[nfr];
                let row2 = n_row - 1;
                let mut row = row1;
                while row <= row2 {
                    if row_degree[row] >= 0 {
                        pivrow[OUT] = row;
                        break;
                    }
                    row += 1;
                }
                front_new1strow[nfr] = row;
            }
            if pivrow[OUT] == EMPTY {
                return UMFPACK_WARNING_singular_matrix;
            }
        }
    }

    if pivrow[IN] != EMPTY {
        freebie[IN] = pivrow[IN] == prior_pivrow[IN] && cdeg1 > 0;
        if !freebie[IN] {
            let fcpos = work.fcpos;
            let fncols = work.fncols;
            let wrpflag = work.wrpflag;
            rdeg[IN] = fncols;
            if fcpos[pivcol] < 0 {
                if rdeg[IN] >= max_rdeg {
                    return UMFPACK_ERROR_different_pattern;
                }
                wrp[pivcol] = wrpflag;
                w_i[rdeg[IN]] = pivcol;
                rdeg[IN] += 1;
            }
            let r = pivrow[IN];
            let tpi = row_tuples[r];
            if tpi != 0 {
                let tp1: P<Tuple> = (m + tpi).cast();
                let mut tp2 = tp1;
                let ntp = row_tlen[r];
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
                    let ncols = ep.ncols;
                    let rows = cols + ncols;
                    if rows[f] == EMPTY {
                        continue;
                    }
                    for j in 0..ncols {
                        let col = cols[j];
                        if col >= 0 && wrp[col] != wrpflag && fcpos[col] < 0 {
                            if rdeg[IN] >= max_rdeg {
                                return UMFPACK_ERROR_different_pattern;
                            }
                            wrp[col] = wrpflag;
                            w_i[rdeg[IN]] = col;
                            rdeg[IN] += 1;
                        }
                    }
                    *tp2.r(0) = tt;
                    tp2 = tp2 + 1;
                }
                row_tlen[r] = tp2.diff(tp1);
            }
            work.wrpflag += 1;
        }
    }

    if pivrow[OUT] != EMPTY {
        freebie[OUT] = pivrow[OUT] == prior_pivrow[OUT] && cdeg1 > 0;
        if !freebie[OUT] {
            let wrpflag = work.wrpflag;
            rdeg[OUT] = 0;
            if rdeg[OUT] >= max_rdeg {
                return UMFPACK_ERROR_different_pattern;
            }
            wrp[pivcol] = wrpflag;
            w_o[rdeg[OUT]] = pivcol;
            rdeg[OUT] += 1;
            let r = pivrow[OUT];
            let tpi = row_tuples[r];
            if tpi != 0 {
                let tp1: P<Tuple> = (m + tpi).cast();
                let mut tp2 = tp1;
                let ntp = row_tlen[r];
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
                    let ncols = ep.ncols;
                    let rows = cols + ncols;
                    if rows[f] == EMPTY {
                        continue;
                    }
                    for j in 0..ncols {
                        let col = cols[j];
                        if col >= 0 && wrp[col] != wrpflag {
                            if rdeg[OUT] >= max_rdeg {
                                return UMFPACK_ERROR_different_pattern;
                            }
                            wrp[col] = wrpflag;
                            w_o[rdeg[OUT]] = col;
                            rdeg[OUT] += 1;
                        }
                    }
                    *tp2.r(0) = tt;
                    tp2 = tp2 + 1;
                }
                row_tlen[r] = tp2.diff(tp1);
            }
            work.wrpflag += 1;
        }
    }
    UMFPACK_OK
}
