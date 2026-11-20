//! `UMFPACK_solve` / `UMF_solve` for `sys = UMFPACK_A`, `UMFPACK_At` and
//! `UMFPACK_Aat` without iterative refinement (Elmer's `umf4sol` sets
//! `Control [UMFPACK_IRSTEP] = 0` and solves `A.'x = b`), with
//! `UMF_lsolve`, `UMF_usolve`, `UMF_ltsolve`, `UMF_utsolve`.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::numtypes::*;
use super::util::*;
use crate::amd::{Int, EMPTY};
use crate::ptr::P;

fn lsolve(numeric: &mut Numeric, mut x: P<Entry>, mut pattern: P<Int>) -> f64 {
    if numeric.n_row != numeric.n_col {
        return 0.0;
    }
    let npiv = numeric.npiv;
    let lpos = numeric.lpos_p();
    let lilen = numeric.lilen_p();
    let lip = numeric.lip_p();
    let n1 = numeric.n1;
    let m = numeric.mem();
    for k in 0..n1 {
        let xk = x[k];
        let deg = lilen[k];
        if deg > 0 && xk != 0.0 {
            let lp = lip[k];
            let li: P<Int> = (m + lp).cast();
            let lval: P<Entry> = (m + (lp + units(SIZEOF_INT, deg))).cast();
            for j in 0..deg {
                x[li[j]] -= xk * lval[j];
            }
        }
    }
    let mut deg: Int = 0;
    for k in n1..npiv {
        let mut lp = lip[k];
        if lp < 0 {
            lp = -lp;
            deg = 0;
        }
        let pos = lpos[k];
        if pos != EMPTY {
            deg -= 1;
            pattern[pos] = pattern[deg];
        }
        let ip: P<Int> = (m + lp).cast();
        let llen = lilen[k];
        for j in 0..llen {
            pattern[deg] = ip[j];
            deg += 1;
        }
        let xk = x[k];
        if xk != 0.0 {
            let xp: P<Entry> = (m + (lp + units(SIZEOF_INT, llen))).cast();
            for j in 0..deg {
                x[pattern[j]] -= xk * xp[j];
            }
        }
    }
    MULTSUB_FLOPS * numeric.lnz as f64
}

fn usolve(numeric: &mut Numeric, mut x: P<Entry>, mut pattern: P<Int>) -> f64 {
    if numeric.n_row != numeric.n_col {
        return 0.0;
    }
    let n = numeric.n_row;
    let npiv = numeric.npiv;
    let upos = numeric.upos_p();
    let uilen = numeric.uilen_p();
    let uip = numeric.uip_p();
    let d = numeric.d_p();
    let n1 = numeric.n1;
    let m = numeric.mem();
    let mut k = n - 1;
    while k >= npiv {
        let xk = x[k];
        x[k] = xk / d[k];
        k -= 1;
    }
    let mut deg = numeric.ulen;
    for j in 0..deg {
        pattern[j] = numeric.upattern[j as usize];
    }
    let mut k = npiv - 1;
    while k >= n1 {
        let mut up = uip[k];
        let ulen = uilen[k];
        let new_uchain = up < 0;
        let xp: P<Entry> = if new_uchain {
            up = -up;
            (m + (up + units(SIZEOF_INT, ulen))).cast()
        } else {
            (m + up).cast()
        };
        let mut xk = x[k];
        for j in 0..deg {
            xk -= x[pattern[j]] * xp[j];
        }
        x[k] = xk / d[k];
        if k == n1 {
            break;
        }
        if new_uchain {
            deg = ulen;
            let ip: P<Int> = (m + up).cast();
            for j in 0..deg {
                pattern[j] = ip[j];
            }
        } else {
            deg -= ulen;
            let pos = upos[k];
            if pos != EMPTY {
                pattern[deg] = pattern[pos];
                deg += 1;
                pattern[pos] = k;
            }
        }
        k -= 1;
    }
    let mut k = n1 - 1;
    while k >= 0 {
        let deg = uilen[k];
        let mut xk = x[k];
        if deg > 0 {
            let up = uip[k];
            let ui: P<Int> = (m + up).cast();
            let uval: P<Entry> = (m + (up + units(SIZEOF_INT, deg))).cast();
            for j in 0..deg {
                xk -= x[ui[j]] * uval[j];
            }
        }
        x[k] = xk / d[k];
        k -= 1;
    }
    DIV_FLOPS * n as f64 + MULTSUB_FLOPS * numeric.unz as f64
}

fn ltsolve(numeric: &mut Numeric, mut x: P<Entry>, mut pattern: P<Int>) -> f64 {
    if numeric.n_row != numeric.n_col {
        return 0.0;
    }
    let npiv = numeric.npiv;
    let lpos = numeric.lpos_p();
    let lilen = numeric.lilen_p();
    let lip = numeric.lip_p();
    let n1 = numeric.n1;
    let m = numeric.mem();
    let mut kstart = npiv;
    let mut kend = npiv - 1;
    while kend >= n1 {
        kstart = kend;
        while kstart >= 0 && lip[kstart] > 0 {
            kstart -= 1;
        }
        let mut deg: Int = 0;
        for k in kstart..=kend {
            let pos = lpos[k];
            if pos != EMPTY {
                deg -= 1;
                pattern[pos] = pattern[deg];
            }
            let mut lp = lip[k];
            if k == kstart {
                lp = -lp;
            }
            let ip: P<Int> = (m + lp).cast();
            let llen = lilen[k];
            for j in 0..llen {
                pattern[deg] = ip[j];
                deg += 1;
            }
        }
        let mut k = kend;
        while k >= kstart {
            let mut lp = lip[k];
            if k == kstart {
                lp = -lp;
            }
            let llen = lilen[k];
            let xp: P<Entry> = (m + (lp + units(SIZEOF_INT, llen))).cast();
            let mut xk = x[k];
            for j in 0..deg {
                xk -= x[pattern[j]] * xp[j];
            }
            x[k] = xk;
            deg -= llen;
            let pos = lpos[k];
            if pos != EMPTY {
                pattern[deg] = pattern[pos];
                deg += 1;
                pattern[pos] = k;
            }
            k -= 1;
        }
        kend = kstart - 1;
    }
    let _ = kstart;
    let mut k = n1 - 1;
    while k >= 0 {
        let deg = lilen[k];
        if deg > 0 {
            let mut xk = x[k];
            let lp = lip[k];
            let li: P<Int> = (m + lp).cast();
            let lval: P<Entry> = (m + (lp + units(SIZEOF_INT, deg))).cast();
            for j in 0..deg {
                xk -= x[li[j]] * lval[j];
            }
            x[k] = xk;
        }
        k -= 1;
    }
    MULTSUB_FLOPS * numeric.lnz as f64
}

fn utsolve(numeric: &mut Numeric, mut x: P<Entry>, mut pattern: P<Int>) -> f64 {
    if numeric.n_row != numeric.n_col {
        return 0.0;
    }
    let n = numeric.n_row;
    let npiv = numeric.npiv;
    let upos = numeric.upos_p();
    let uilen = numeric.uilen_p();
    let uip = numeric.uip_p();
    let d = numeric.d_p();
    let n1 = numeric.n1;
    let m = numeric.mem();
    for k in 0..n1 {
        let xk = x[k] / d[k];
        x[k] = xk;
        let deg = uilen[k];
        if deg > 0 && xk != 0.0 {
            let up = uip[k];
            let ui: P<Int> = (m + up).cast();
            let uval: P<Entry> = (m + (up + units(SIZEOF_INT, deg))).cast();
            for j in 0..deg {
                x[ui[j]] -= xk * uval[j];
            }
        }
    }
    let mut kstart = n1;
    while kstart < npiv {
        let mut kend = kstart;
        while kend < npiv && uip[kend + 1] > 0 {
            kend += 1;
        }
        let k = kend + 1;
        let mut deg: Int;
        if k == npiv {
            deg = numeric.ulen;
            for j in 0..deg {
                pattern[j] = numeric.upattern[j as usize];
            }
        } else {
            let up = -uip[k];
            deg = uilen[k];
            let ip: P<Int> = (m + up).cast();
            for j in 0..deg {
                pattern[j] = ip[j];
            }
        }
        let mut uhead = n;
        let mut k = kend;
        while k > kstart {
            let ulen = uilen[k];
            for _ in 0..ulen {
                uhead -= 1;
                deg -= 1;
                pattern[uhead] = pattern[deg];
            }
            let pos = upos[k];
            if pos != EMPTY {
                pattern[deg] = pattern[pos];
                deg += 1;
                pattern[pos] = k;
            }
            k -= 1;
        }
        for k in kstart..=kend {
            let pos = upos[k];
            if pos != EMPTY {
                deg -= 1;
                pattern[pos] = pattern[deg];
            }
            let mut up = uip[k];
            let ulen = uilen[k];
            if k > kstart {
                for _ in 0..ulen {
                    pattern[deg] = pattern[uhead];
                    deg += 1;
                    uhead += 1;
                }
            }
            let xk = x[k] / d[k];
            x[k] = xk;
            if xk != 0.0 {
                let xp: P<Entry> = if k == kstart {
                    up = -up;
                    (m + (up + units(SIZEOF_INT, ulen))).cast()
                } else {
                    (m + up).cast()
                };
                for j in 0..deg {
                    x[pattern[j]] -= xk * xp[j];
                }
            }
        }
        kstart = kend + 1;
    }
    for k in npiv..n {
        let xk = x[k] / d[k];
        x[k] = xk;
    }
    DIV_FLOPS * n as f64 + MULTSUB_FLOPS * numeric.unz as f64
}

/// `UMFPACK_solve` with `irstep = 0` for `sys` = `UMFPACK_A`, `UMFPACK_At`
/// or `UMFPACK_Aat` (the latter two are the same for real matrices).
/// Returns the status (`UMFPACK_OK` or `UMFPACK_WARNING_singular_matrix`).
pub fn umfpack_solve(
    sys: Int,
    xx: &mut [f64],
    bx: &[f64],
    numeric: &mut Numeric,
    info: Option<&mut [f64]>,
) -> Int {
    let mut info2 = [EMPTY as f64; UMFPACK_INFO];
    let info: &mut [f64] = match info {
        Some(i) => {
            for k in UMFPACK_IR_TAKEN..=UMFPACK_SOLVE_TIME {
                i[k] = EMPTY as f64;
            }
            i
        }
        None => &mut info2,
    };
    info[UMFPACK_STATUS] = UMFPACK_OK as f64;
    info[UMFPACK_SOLVE_FLOPS] = 0.0;
    if numeric.valid != NUMERIC_VALID {
        info[UMFPACK_STATUS] = UMFPACK_ERROR_invalid_Numeric_object as f64;
        return UMFPACK_ERROR_invalid_Numeric_object;
    }
    info[UMFPACK_NROW] = numeric.n_row as f64;
    info[UMFPACK_NCOL] = numeric.n_col as f64;
    if numeric.n_row != numeric.n_col {
        info[UMFPACK_STATUS] = UMFPACK_ERROR_invalid_system as f64;
        return UMFPACK_ERROR_invalid_system;
    }
    let n = numeric.n_row;
    let mut pattern_v = vec![0 as Int; n.max(1) as usize];
    let mut w_v = vec![0.0f64; n.max(1) as usize];
    let pattern = P::new(&mut pattern_v);
    let mut w = P::new(&mut w_v);

    // UMF_solve
    let rperm = numeric.rperm_p();
    let cperm = numeric.cperm_p();
    let rs = numeric.rs_p();
    let do_scale = !rs.is_null();
    let do_recip = numeric.do_recip;
    let mut flops = 0.0;
    info[UMFPACK_SOLVE_FLOPS] = 0.0;
    info[UMFPACK_IR_TAKEN] = 0.0;
    info[UMFPACK_IR_ATTEMPTED] = 0.0;
    let status = if numeric.nnzpiv < n || numeric.rcond == 0.0 || numeric.rcond.is_nan() {
        UMFPACK_WARNING_singular_matrix
    } else {
        UMFPACK_OK
    };
    let mut x = P::new(xx);
    let b = bx;
    if sys == UMFPACK_A {
        if do_scale {
            if do_recip {
                for i in 0..n {
                    x[i] = b[i as usize];
                    x[i] *= rs[i];
                }
            } else {
                for i in 0..n {
                    x[i] = b[i as usize];
                    x[i] /= rs[i];
                }
            }
            flops += SCALE_FLOPS * n as f64;
            for i in 0..n {
                w[i] = x[rperm[i]];
            }
        } else {
            for i in 0..n {
                w[i] = b[rperm[i] as usize];
            }
        }
        flops += lsolve(numeric, w, pattern);
        flops += usolve(numeric, w, pattern);
        for i in 0..n {
            x[cperm[i]] = w[i];
        }
    } else if sys == UMFPACK_At || sys == UMFPACK_Aat {
        for i in 0..n {
            w[i] = b[cperm[i] as usize];
        }
        flops += utsolve(numeric, w, pattern);
        flops += ltsolve(numeric, w, pattern);
        for i in 0..n {
            x[rperm[i]] = w[i];
        }
        if do_scale {
            if do_recip {
                for i in 0..n {
                    x[i] *= rs[i];
                }
            } else {
                for i in 0..n {
                    x[i] /= rs[i];
                }
            }
            flops += SCALE_FLOPS * n as f64;
        }
    } else {
        info[UMFPACK_STATUS] = UMFPACK_ERROR_invalid_system as f64;
        return UMFPACK_ERROR_invalid_system;
    }
    info[UMFPACK_SOLVE_FLOPS] = flops;
    info[UMFPACK_STATUS] = status as f64;
    status
}
