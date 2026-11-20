//! `UMFPACK_symbolic` / `UMFPACK_qsymbolic` (no user column ordering): the
//! symbolic analysis, ending in the `Symbolic` object the numeric
//! factorisation reads.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::analyze::umf_analyze;
use super::colamd::*;
use super::consts::*;
use super::singletons::{umf_singletons, SingletonInfo};
use super::twobytwo::umf_2by2;
use super::util::*;
use crate::amd::{self, Int, EMPTY};
use crate::ptr::P;

/// `SymbolicType`.
#[derive(Default, Clone, Debug)]
pub struct Symbolic {
    pub num_mem_usage_est: f64,
    pub num_mem_size_est: f64,
    pub peak_sym_usage: f64,
    pub sym: f64,
    pub dnum_mem_init_usage: f64,
    pub amd_lunz: f64,
    pub lunz_bound: f64,
    pub valid: Int,
    pub max_nchains: Int,
    pub nchains: Int,
    pub chain_start: Vec<Int>,
    pub chain_maxrows: Vec<Int>,
    pub chain_maxcols: Vec<Int>,
    pub maxnrows: Int,
    pub maxncols: Int,
    pub front_npivcol: Vec<Int>,
    pub front_1strow: Vec<Int>,
    pub front_leftmostdesc: Vec<Int>,
    pub front_parent: Vec<Int>,
    pub cperm_init: Vec<Int>,
    pub rperm_init: Vec<Int>,
    pub cdeg: Vec<Int>,
    pub rdeg: Vec<Int>,
    /// `Esize` (`None` when `esize == 0`).
    pub esize_vec: Option<Vec<Int>>,
    pub dense_row_threshold: Int,
    pub n1: Int,
    pub nempty: Int,
    pub diagonal_map: Option<Vec<Int>>,
    pub esize: Int,
    pub nfr: Int,
    pub n_row: Int,
    pub n_col: Int,
    pub nz: Int,
    pub nb: Int,
    pub num_mem_init_usage: Int,
    pub nempty_row: Int,
    pub nempty_col: Int,
    pub strategy: Int,
    pub ordering: Int,
    pub fix_q: bool,
    pub prefer_diagonal: bool,
    pub nzaat: Int,
    pub nzdiag: Int,
    pub amd_dmax: Int,
}

/// `UMFPACK_symbolic`.
pub fn umfpack_symbolic(
    n_row: Int,
    n_col: Int,
    ap: &[Int],
    ai: &[Int],
    ax: Option<&[f64]>,
    control: Option<&[f64]>,
    info: Option<&mut [f64]>,
) -> Result<Symbolic, Int> {
    umfpack_qsymbolic(n_row, n_col, ap, ai, ax, control, info)
}

#[allow(clippy::too_many_arguments)]
fn do_amd(
    n: Int,
    ap: P<Int>,
    ai: P<Int>,
    q: P<Int>,
    qinv: P<Int>,
    sdeg: P<Int>,
    clen: Int,
    ci: P<Int>,
    amd_control: &[f64],
    amd_info: &mut [f64],
    symbolic: &mut Symbolic,
    info: &mut [f64],
) {
    if n == 0 {
        symbolic.amd_dmax = 0;
        symbolic.amd_lunz = 0.0;
        info[UMFPACK_SYMMETRIC_LUNZ] = 0.0;
        info[UMFPACK_SYMMETRIC_FLOPS] = 0.0;
        info[UMFPACK_SYMMETRIC_DMAX] = 0.0;
        info[UMFPACK_SYMMETRIC_NDENSE] = 0.0;
    } else {
        let nu = n as usize;
        // SAFETY: AMD_1 reads Ap [0..n], Ai [0..Ap[n]] and writes the first n
        // entries of Q, Qinv, Sdeg and Clen entries of Ci, all disjoint.
        unsafe {
            let apl = core::slice::from_raw_parts(ap.r(0) as *const Int, nu + 1);
            let nz = apl[nu] as usize;
            let ail: &[Int] = if nz == 0 {
                &[]
            } else {
                core::slice::from_raw_parts(ai.r(0) as *const Int, nz)
            };
            amd::amd_1(
                n,
                apl,
                ail,
                core::slice::from_raw_parts_mut(q.r(0) as *mut Int, nu),
                core::slice::from_raw_parts_mut(qinv.r(0) as *mut Int, nu),
                core::slice::from_raw_parts_mut(sdeg.r(0) as *mut Int, nu),
                clen,
                core::slice::from_raw_parts_mut(ci.r(0) as *mut Int, clen as usize),
                Some(amd_control),
                Some(amd_info),
            );
        }
        symbolic.amd_dmax = amd_info[amd::AMD_DMAX] as Int;
        symbolic.amd_lunz = 2.0 * amd_info[amd::AMD_LNZ] + n as f64;
        info[UMFPACK_SYMMETRIC_LUNZ] = symbolic.amd_lunz;
        info[UMFPACK_SYMMETRIC_FLOPS] =
            DIV_FLOPS * amd_info[amd::AMD_NDIV] + MULTSUB_FLOPS * amd_info[amd::AMD_NMULTSUBS_LU];
        info[UMFPACK_SYMMETRIC_DMAX] = symbolic.amd_dmax as f64;
        info[UMFPACK_SYMMETRIC_NDENSE] = amd_info[amd::AMD_NDENSE];
        info[UMFPACK_SYMBOLIC_DEFRAG] += amd_info[amd::AMD_NCMPA];
    }
}

#[allow(clippy::too_many_arguments)]
fn prune_singletons(
    n1: Int,
    n_col: Int,
    ap: P<Int>,
    ai: P<Int>,
    ax: P<f64>,
    cperm1: P<Int>,
    inv_rperm1: P<Int>,
    mut si: P<Int>,
    mut sp: P<Int>,
) -> Int {
    let mut nzdiag = 0;
    let do_nzdiag = !ax.is_null();
    let mut pp = 0;
    for k in n1..n_col {
        let oldcol = cperm1[k];
        let newcol = k - n1;
        sp[newcol] = pp;
        for p in ap[oldcol]..ap[oldcol + 1] {
            let row = ai[p];
            let newrow = inv_rperm1[row] - n1;
            if newrow >= 0 {
                si[pp] = newrow;
                pp += 1;
                if do_nzdiag && newrow == newcol && ax[p] != 0.0 {
                    nzdiag += 1;
                }
            }
        }
    }
    sp[n_col - n1] = pp;
    nzdiag
}

fn combine_ordering(
    n1: Int,
    nempty_col: Int,
    n_col: Int,
    mut cperm_init: P<Int>,
    cperm1: P<Int>,
    qinv: P<Int>,
) {
    for k in 0..n1 {
        cperm_init[k] = cperm1[k];
    }
    for k in n1..n_col - nempty_col {
        let oldcol = cperm1[k];
        let newcol = k - n1;
        let knew = qinv[newcol] + n1;
        cperm_init[knew] = oldcol;
    }
    for k in n_col - nempty_col..n_col {
        cperm_init[k] = cperm1[k];
    }
}

/// `ELEMENT_SIZE (r, c)` of qsymbolic (doubles).
fn element_size(r: f64, c: f64) -> f64 {
    dget_element_size(r, c) + 1.0 + (r + c) * units(SIZEOF_TUPLE, 1) as f64
}

/// `DGET_ELEMENT_SIZE (nr, nc)`.
pub fn dget_element_size(nr: f64, nc: f64) -> f64 {
    dunits(SIZEOF_ELEMENT, 1.0) + dunits(SIZEOF_INT, nc + nr) + dunits(SIZEOF_ENTRY, nc * nr)
}

/// `GET_ELEMENT_SIZE (nr, nc)`.
pub fn get_element_size(nr: Int, nc: Int) -> Int {
    units(SIZEOF_ELEMENT, 1) + units(SIZEOF_INT, nc + nr) + units(SIZEOF_ENTRY, nc * nr)
}

/// `TUPLES (t)`.
pub fn tuples(t: Int) -> Int {
    4.max(t + 1)
}

/// `UMFPACK_qsymbolic` with `Quser == NULL`.
pub fn umfpack_qsymbolic(
    n_row: Int,
    n_col: Int,
    ap_s: &[Int],
    ai_s: &[Int],
    ax_s: Option<&[f64]>,
    control: Option<&[f64]>,
    user_info: Option<&mut [f64]>,
) -> Result<Symbolic, Int> {
    let mut info2 = [0.0f64; UMFPACK_INFO];
    let info: &mut [f64] = match user_info {
        Some(i) => i,
        None => &mut info2,
    };

    let drow = get_control(control, UMFPACK_DENSE_ROW, UMFPACK_DEFAULT_DENSE_ROW);
    let dcol = get_control(control, UMFPACK_DENSE_COL, UMFPACK_DEFAULT_DENSE_COL);
    let mut nb = get_control(control, UMFPACK_BLOCK_SIZE, UMFPACK_DEFAULT_BLOCK_SIZE) as Int;
    let mut strategy =
        get_control(control, UMFPACK_STRATEGY, UMFPACK_DEFAULT_STRATEGY as f64) as Int;
    let mut tol = get_control(
        control,
        UMFPACK_2BY2_TOLERANCE,
        UMFPACK_DEFAULT_2BY2_TOLERANCE,
    );
    let mut scale = get_control(control, UMFPACK_SCALE, UMFPACK_DEFAULT_SCALE as f64) as Int;
    let force_fix_q = get_control(control, UMFPACK_FIXQ, UMFPACK_DEFAULT_FIXQ);
    let mut amd_control = [0.0f64; amd::AMD_CONTROL];
    amd::amd_defaults(&mut amd_control);
    amd_control[amd::AMD_DENSE] =
        get_control(control, UMFPACK_AMD_DENSE, UMFPACK_DEFAULT_AMD_DENSE);
    let aggressive = get_control(control, UMFPACK_AGGRESSIVE, UMFPACK_DEFAULT_AGGRESSIVE) != 0.0;
    amd_control[amd::AMD_AGGRESSIVE] = if aggressive { 1.0 } else { 0.0 };
    nb = nb.max(2);
    nb = nb.min(MAXNB);
    if nb % 2 == 1 {
        nb += 1;
    }
    tol = cmax(0.0, cmin(tol, 1.0));
    if scale != UMFPACK_SCALE_NONE && scale != UMFPACK_SCALE_MAX {
        scale = UMFPACK_DEFAULT_SCALE;
    }

    for v in info.iter_mut().take(UMFPACK_INFO) {
        *v = EMPTY as f64;
    }
    let nn = n_row.max(n_col);
    let n_inner = n_row.min(n_col);
    info[UMFPACK_STATUS] = UMFPACK_OK as f64;
    info[UMFPACK_NROW] = n_row as f64;
    info[UMFPACK_NCOL] = n_col as f64;
    info[UMFPACK_SIZE_OF_UNIT] = SIZEOF_UNIT as f64;
    info[UMFPACK_SIZE_OF_INT] = 4.0;
    info[UMFPACK_SIZE_OF_LONG] = 8.0;
    info[UMFPACK_SIZE_OF_POINTER] = SIZEOF_PTR as f64;
    info[UMFPACK_SIZE_OF_ENTRY] = SIZEOF_ENTRY as f64;
    info[UMFPACK_SYMBOLIC_DEFRAG] = 0.0;

    if n_row <= 0 || n_col <= 0 {
        info[UMFPACK_STATUS] = UMFPACK_ERROR_n_nonpositive as f64;
        return Err(UMFPACK_ERROR_n_nonpositive);
    }
    let nz = ap_s[n_col as usize];
    info[UMFPACK_NZ] = nz as f64;
    if nz < 0 {
        info[UMFPACK_STATUS] = UMFPACK_ERROR_invalid_matrix as f64;
        return Err(UMFPACK_ERROR_invalid_matrix);
    }
    if n_row != n_col {
        strategy = UMFPACK_STRATEGY_UNSYMMETRIC;
    }
    if !(UMFPACK_STRATEGY_AUTO..=UMFPACK_STRATEGY_SYMMETRIC).contains(&strategy) {
        strategy = UMFPACK_STRATEGY_AUTO;
    }
    let user_auto_strategy = strategy == UMFPACK_STRATEGY_AUTO;

    // sizes of the workspace
    let mut d_clen = colamd_recommended_f(nz as f64, n_row as f64, n_col as f64);
    let d_clen_analyze =
        n_col as f64 + cmax(nz as f64, n_col as f64) + 3.0 * nn as f64 + 1.0 + n_col as f64;
    d_clen = cmax(d_clen, d_clen_analyze);
    let mut d_clen_amd = 2.4 * nz as f64 + 8.0 * n_inner as f64;
    d_clen_amd += nn.max(nz) as f64;
    d_clen = cmax(d_clen, d_clen_amd);
    let sym_work_usage = |clen: f64| -> f64 {
        dunits(SIZEOF_INT, clen)
            + dunits(SIZEOF_INT, nz as f64)
            + 4.0 * dunits(SIZEOF_INT, n_row as f64)
            + 4.0 * dunits(SIZEOF_INT, n_col as f64)
            + 2.0 * dunits(SIZEOF_INT, (n_col + 1) as f64)
            + dunits(SIZEOF_ENTRY, n_row as f64)
    };
    info[UMFPACK_SYMBOLIC_PEAK_MEMORY] =
        sym_work_usage(d_clen) + umf_symbolic_usage(n_row, n_col, n_col, n_col, n_col, true);
    if int_overflow(d_clen * SIZEOF_INT as f64) {
        info[UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory as f64;
        return Err(UMFPACK_ERROR_out_of_memory);
    }
    let mut clen = colamd_recommended(nz, n_row, n_col);
    let clen_analyze = n_col + nz.max(n_col) + 3 * nn + 1 + n_col;
    clen = clen.max(clen_analyze);
    let mut clen_amd = (2.4 * nz as f64 + (8 * n_inner) as f64) as Int;
    clen_amd += nn.max(nz);
    clen = clen.max(clen_amd);

    let mut symbolic = Symbolic {
        valid: 0,
        esize: 0,
        ..Default::default()
    };
    symbolic.cperm_init = vec![0; n_col as usize + 1];
    symbolic.rperm_init = vec![0; n_row as usize + 1];
    symbolic.cdeg = vec![0; n_col as usize + 1];
    symbolic.rdeg = vec![0; n_row as usize + 1];
    let mut cperm_init = P::new(&mut symbolic.cperm_init);
    let mut rperm_init = P::new(&mut symbolic.rperm_init);
    let mut cdeg = P::new(&mut symbolic.cdeg);
    let mut rdeg = P::new(&mut symbolic.rdeg);
    symbolic.n_row = n_row;
    symbolic.n_col = n_col;
    symbolic.nz = nz;
    symbolic.nb = nb;

    // the input matrix (read only)
    let mut ap_v = ap_s.to_vec();
    let mut ai_v = if ai_s.is_empty() {
        vec![0]
    } else {
        ai_s[..nz as usize].to_vec()
    };
    let mut ax_v: Vec<f64> = match ax_s {
        Some(x) if !x.is_empty() => x[..nz as usize].to_vec(),
        Some(_) => vec![0.0],
        None => Vec::new(),
    };
    let ap = P::new(&mut ap_v);
    let ai = P::new(&mut ai_v);
    let ax: P<f64> = if ax_s.is_some() {
        P::new(&mut ax_v)
    } else {
        P::null()
    };

    // SW workspace
    let mut sw_si = vec![0 as Int; nz.max(1) as usize];
    let mut sw_sp = vec![0 as Int; n_col as usize + 1];
    let mut sw_inv_rperm1 = vec![0 as Int; n_row as usize];
    let mut sw_cperm1 = vec![0 as Int; n_col as usize];
    let mut sw_ci = vec![0 as Int; clen as usize];
    let mut sw_front_npivcol = vec![0 as Int; n_col as usize + 1];
    let mut sw_front_nrows = vec![0 as Int; n_col as usize];
    let mut sw_front_ncols = vec![0 as Int; n_col as usize];
    let mut sw_front_parent = vec![0 as Int; n_col as usize];
    let mut sw_front_cols = vec![0 as Int; n_col as usize];
    let mut sw_rperm1 = vec![0 as Int; n_row as usize];
    let mut sw_in_front = vec![0 as Int; n_row as usize];
    let mut sw_rperm_2by2: Vec<Int>;
    let mut sw_rs: Vec<f64> = Vec::new();

    let mut ci = P::new(&mut sw_ci);
    let mut fr_npivcol = P::new(&mut sw_front_npivcol);
    let fr_nrows = P::new(&mut sw_front_nrows);
    let fr_ncols = P::new(&mut sw_front_ncols);
    let fr_parent = P::new(&mut sw_front_parent);
    let fr_cols = P::new(&mut sw_front_cols);
    let cperm1 = P::new(&mut sw_cperm1);
    let rperm1 = P::new(&mut sw_rperm1);
    let si = P::new(&mut sw_si);
    let sp = P::new(&mut sw_sp);
    let inv_rperm1 = P::new(&mut sw_inv_rperm1);
    let mut rperm_2by2: P<Int> = P::null();
    let mut in_front = P::new(&mut sw_in_front);

    // singletons
    let mut sg = SingletonInfo::default();
    let status = umf_singletons(
        n_row,
        n_col,
        ap,
        ai,
        cdeg,
        cperm1,
        rdeg,
        rperm1,
        inv_rperm1,
        &mut sg,
        rperm_init,
        ci,
        ci + nz,
        ci + (nz + n_row),
    );
    if status != UMFPACK_OK {
        info[UMFPACK_STATUS] = status as f64;
        return Err(status);
    }
    let SingletonInfo {
        n1,
        n1c,
        n1r,
        nempty_col,
        mut nempty_row,
        is_sym,
        max_rdeg,
    } = sg;
    info[UMFPACK_NEMPTY_COL] = nempty_col as f64;
    info[UMFPACK_NEMPTY_ROW] = nempty_row as f64;
    info[UMFPACK_NDENSE_COL] = 0.0;
    info[UMFPACK_NDENSE_ROW] = 0.0;
    info[UMFPACK_COL_SINGLETONS] = n1c as f64;
    info[UMFPACK_ROW_SINGLETONS] = n1r as f64;
    info[UMFPACK_S_SYMMETRIC] = if is_sym { 1.0 } else { 0.0 };
    let nempty = nempty_col.min(nempty_row);
    symbolic.nempty_row = nempty_row;
    symbolic.nempty_col = nempty_col;
    symbolic.n1 = n1;
    symbolic.nempty = nempty;
    let n2 = nn - n1 - nempty;
    let dense_row_threshold = dense_degree_threshold(drow, n_col - n1 - nempty_col);
    symbolic.dense_row_threshold = dense_row_threshold;
    if !is_sym {
        strategy = UMFPACK_STRATEGY_UNSYMMETRIC;
    }

    let wq = rperm_init;
    let sdeg = cperm_init;
    let mut sym = EMPTY as f64;
    let mut nzaat = EMPTY;
    let mut nzdiag = EMPTY;
    let mut amd_info = [EMPTY as f64; amd::AMD_INFO];

    if strategy != UMFPACK_STRATEGY_UNSYMMETRIC {
        nzdiag = prune_singletons(n1, nn, ap, ai, ax, cperm1, inv_rperm1, si, sp);
        nzaat = amd_aat_p(n2, sp, si, sdeg, wq, &mut amd_info);
        sym = amd_info[amd::AMD_SYMMETRY];
        info[UMFPACK_N2] = n2 as f64;
    }
    symbolic.sym = sym;
    symbolic.nzaat = nzaat;
    symbolic.nzdiag = nzdiag;
    symbolic.amd_dmax = EMPTY;
    info[UMFPACK_PATTERN_SYMMETRY] = sym;
    info[UMFPACK_NZ_A_PLUS_AT] = nzaat as f64;
    info[UMFPACK_NZDIAG] = nzdiag as f64;

    if strategy == UMFPACK_STRATEGY_AUTO {
        if sym < 0.10 {
            strategy = UMFPACK_STRATEGY_UNSYMMETRIC;
        } else if sym >= 0.7 && nzdiag == n2 {
            strategy = UMFPACK_STRATEGY_SYMMETRIC;
        } else {
            strategy = UMFPACK_STRATEGY_2BY2;
        }
    }

    if strategy == UMFPACK_STRATEGY_2BY2 {
        sw_rperm_2by2 = vec![0; nn as usize];
        rperm_2by2 = P::new(&mut sw_rperm_2by2);
        let clen3 = clen - (nn.max(nz) + 1);
        let ri = ci + clen3;
        let rp = fr_npivcol;
        let blen = fr_ncols;
        let mut w = fr_cols;
        let rs: P<f64> = if scale != UMFPACK_SCALE_NONE {
            sw_rs = vec![0.0; nn as usize];
            P::new(&mut sw_rs)
        } else {
            P::null()
        };
        let (nweak, unmatched) = umf_2by2(
            nn, ap, ai, ax, tol, scale, cperm1, inv_rperm1, n1, nempty, sdeg, rperm_2by2, ri, rp,
            rs, blen, w, ci, wq,
        );
        info[UMFPACK_2BY2_NWEAK] = nweak as f64;
        info[UMFPACK_2BY2_UNMATCHED] = unmatched as f64;
        drop(core::mem::take(&mut sw_rs));
        umf_transpose(
            n2,
            n2,
            sp,
            si,
            P::null(),
            rperm_2by2,
            P::null(),
            0,
            rp,
            ri,
            P::null(),
            w,
        );
        for v in amd_info.iter_mut() {
            *v = EMPTY as f64;
        }
        let nz_papat = amd_aat_p(n2, rp, ri, sdeg, wq, &mut amd_info);
        let sym2 = amd_info[amd::AMD_SYMMETRY];
        let nzd2 = amd_info[amd::AMD_NZDIAG];
        info[UMFPACK_2BY2_PATTERN_SYMMETRY] = sym2;
        info[UMFPACK_2BY2_NZ_PA_PLUS_PAT] = nz_papat as f64;
        info[UMFPACK_2BY2_NZDIAG] = nzd2;
        if user_auto_strategy {
            if (sym2 > 1.1 * sym) && (nzd2 > 0.9 * n2 as f64) {
                strategy = UMFPACK_STRATEGY_2BY2;
            } else if sym2 < 0.7 * sym {
                strategy = UMFPACK_STRATEGY_UNSYMMETRIC;
            } else if sym2 < 0.25 {
                strategy = UMFPACK_STRATEGY_UNSYMMETRIC;
            } else if sym2 >= 0.51 {
                strategy = UMFPACK_STRATEGY_2BY2;
            } else if sym2 >= 0.999 * sym {
                strategy = UMFPACK_STRATEGY_2BY2;
            } else {
                strategy = UMFPACK_STRATEGY_UNSYMMETRIC;
            }
        }
        if strategy == UMFPACK_STRATEGY_2BY2 {
            let qinv = fr_npivcol;
            do_amd(
                n2,
                rp,
                ri,
                wq,
                qinv,
                sdeg,
                clen3,
                ci,
                &amd_control,
                &mut amd_info,
                &mut symbolic,
                info,
            );
            combine_ordering(n1, nempty, nn, cperm_init, cperm1, qinv);
            for k in 0..n1 {
                let oldcol = cperm1[k];
                let oldrow = rperm1[k];
                w[oldcol] = oldrow;
            }
            for k in n1..nn - nempty {
                let oldcol = cperm1[k];
                let i = rperm_2by2[k - n1] + n1;
                let oldrow = rperm1[i];
                w[oldcol] = oldrow;
            }
            for k in nn - nempty..nn {
                let oldcol = cperm1[k];
                let oldrow = rperm1[k];
                w[oldcol] = oldrow;
            }
            for k in 0..nn {
                rperm_2by2[k] = w[k];
            }
        }
    }

    let fix_q;
    let prefer_diagonal;
    if strategy == UMFPACK_STRATEGY_SYMMETRIC || strategy == UMFPACK_STRATEGY_2BY2 {
        symbolic.ordering = UMFPACK_ORDERING_AMD;
        fix_q = force_fix_q >= 0.0;
        prefer_diagonal = true;
    } else {
        symbolic.ordering = UMFPACK_ORDERING_COLAMD;
        fix_q = force_fix_q > 0.0;
        prefer_diagonal = false;
    }
    symbolic.strategy = strategy;
    symbolic.fix_q = fix_q;
    symbolic.prefer_diagonal = prefer_diagonal;
    info[UMFPACK_STRATEGY_USED] = strategy as f64;
    info[UMFPACK_ORDERING_USED] = symbolic.ordering as f64;
    info[UMFPACK_QFIXED] = if fix_q { 1.0 } else { 0.0 };
    info[UMFPACK_DIAG_PREFERRED] = if prefer_diagonal { 1.0 } else { 0.0 };

    if strategy == UMFPACK_STRATEGY_SYMMETRIC {
        let qinv = fr_npivcol;
        do_amd(
            n2,
            sp,
            si,
            wq,
            qinv,
            sdeg,
            clen,
            ci,
            &amd_control,
            &mut amd_info,
            &mut symbolic,
            info,
        );
        combine_ordering(n1, nempty, nn, cperm_init, cperm1, qinv);
    }

    let do_umf_analyze;
    let mut nfr: Int = 0;
    if strategy == UMFPACK_STRATEGY_UNSYMMETRIC {
        prune_singletons(
            n1,
            n_col,
            ap,
            ai,
            P::null(),
            cperm1,
            inv_rperm1,
            ci,
            cperm_init,
        );
        let mut knobs = [0.0f64; COLAMD_KNOBS];
        umf_colamd_set_defaults(&mut knobs);
        knobs[COLAMD_DENSE_ROW] = drow;
        knobs[COLAMD_DENSE_COL] = dcol;
        knobs[COLAMD_AGGRESSIVE] = if aggressive { 1.0 } else { 0.0 };
        let mut colamd_stats = [0 as Int; COLAMD_STATS];
        umf_colamd(
            n_row - n1 - nempty_row,
            n_col - n1 - nempty_col,
            clen,
            ci,
            cperm_init,
            &knobs,
            &mut colamd_stats,
            fr_npivcol,
            fr_nrows,
            fr_ncols,
            fr_parent,
            fr_cols,
            &mut nfr,
            in_front,
        );
        info[UMFPACK_NDENSE_ROW] = colamd_stats[COLAMD_DENSE_ROW] as f64;
        info[UMFPACK_NDENSE_COL] = colamd_stats[COLAMD_DENSE_COL] as f64;
        info[UMFPACK_SYMBOLIC_DEFRAG] = colamd_stats[COLAMD_DEFRAG_COUNT] as f64;
        do_umf_analyze = colamd_stats[COLAMD_DENSE_ROW] > 0 || colamd_stats[COLAMD_DENSE_COL] > 0;
        combine_ordering(n1, nempty_col, n_col, cperm_init, cperm1, ci);
    } else {
        do_umf_analyze = true;
    }
    cperm_init[n_col] = EMPTY;

    if do_umf_analyze {
        prune_singletons(n1, n_col, ap, ai, P::null(), cperm_init, inv_rperm1, si, sp);
        let clen0 = clen - (nn + 1 + 2 * nn + n_col);
        let mut bp = ci + clen0;
        let link = bp + (nn + 1);
        let mut w = link + nn;
        let cperm2 = w + nn;
        for row in 0..n_row - n1 {
            w[row] = 0;
        }
        let mut pperm = link;
        let mut k = 0;
        for col in 0..n_col - n1 {
            for p in sp[col]..sp[col + 1] {
                let row = si[p];
                if w[row] == 0 {
                    w[row] = 1;
                    pperm[k] = row;
                    k += 1;
                }
            }
        }
        nempty_row = n_row - n1 - k;
        if k < n_row - n1 {
            for row in 0..n_row - n1 {
                if w[row] == 0 {
                    pperm[k] = row;
                    k += 1;
                }
            }
        }
        let mut clen2 = clen0;
        let snz = sp[n_col - n1];
        let bsize = snz.max(1);
        clen2 -= bsize;
        let bi = ci + clen2;
        umf_transpose(
            n_row - n1,
            n_col - n1 - nempty_col,
            sp,
            si,
            P::null(),
            pperm,
            P::null(),
            0,
            bp,
            bi,
            P::null(),
            w,
        );
        for i in 0..=n_row - n1 {
            bp[i] += clen2;
        }
        let mut analyze_compactions = 0;
        let ok = umf_analyze(
            n_row - n1 - nempty_row,
            n_col - n1 - nempty_col,
            ci,
            bp,
            cperm2,
            fix_q,
            w,
            link,
            fr_ncols,
            fr_nrows,
            fr_npivcol,
            fr_parent,
            &mut nfr,
            &mut analyze_compactions,
        );
        if !ok {
            info[UMFPACK_STATUS] = UMFPACK_ERROR_internal_error as f64;
            return Err(UMFPACK_ERROR_internal_error);
        }
        info[UMFPACK_SYMBOLIC_DEFRAG] += analyze_compactions as f64;
        if !fix_q {
            for k in 0..n_col - n1 - nempty_col {
                w[k] = cperm_init[n1 + cperm2[k]];
            }
            for k in 0..n_col - n1 - nempty_col {
                cperm_init[n1 + k] = w[k];
            }
        }
    }
    drop(core::mem::take(&mut sw_si));
    drop(core::mem::take(&mut sw_sp));
    drop(core::mem::take(&mut sw_cperm1));

    // chains and fronts
    let mut nchains = 0;
    for i in 0..nfr {
        if fr_parent[i] != i + 1 {
            nchains += 1;
        }
    }
    symbolic.nchains = nchains;
    symbolic.nfr = nfr;
    symbolic.esize = if max_rdeg > dense_row_threshold {
        n_col - n1 - nempty_col
    } else {
        0
    };
    info[UMFPACK_SYMBOLIC_SIZE] =
        umf_symbolic_usage(n_row, n_col, nchains, nfr, symbolic.esize, prefer_diagonal);
    info[UMFPACK_SYMBOLIC_PEAK_MEMORY] = sym_work_usage(clen as f64) + info[UMFPACK_SYMBOLIC_SIZE];
    symbolic.peak_sym_usage = info[UMFPACK_SYMBOLIC_PEAK_MEMORY];

    symbolic.front_npivcol = vec![0; nfr as usize + 1];
    symbolic.front_parent = vec![0; nfr as usize + 1];
    symbolic.front_1strow = vec![0; nfr as usize + 1];
    symbolic.front_leftmostdesc = vec![0; nfr as usize + 1];
    symbolic.chain_start = vec![0; nchains as usize + 1];
    symbolic.chain_maxrows = vec![0; nchains as usize + 1];
    symbolic.chain_maxcols = vec![0; nchains as usize + 1];
    if symbolic.esize > 0 {
        symbolic.esize_vec = Some(vec![0; symbolic.esize as usize]);
    }
    let mut front_npivcol = P::new(&mut symbolic.front_npivcol);
    let mut front_parent = P::new(&mut symbolic.front_parent);
    let mut front_1strow = P::new(&mut symbolic.front_1strow);
    let mut front_leftmostdesc = P::new(&mut symbolic.front_leftmostdesc);
    let mut chain_start = P::new(&mut symbolic.chain_start);
    let mut chain_maxrows = P::new(&mut symbolic.chain_maxrows);
    let mut chain_maxcols = P::new(&mut symbolic.chain_maxcols);
    let mut esize_p: P<Int> = match symbolic.esize_vec.as_mut() {
        Some(v) => P::new(v),
        None => P::null(),
    };

    if do_umf_analyze {
        for row in 0..n_row {
            in_front[row] = nfr;
        }
        for k in 0..n1 {
            let row = rperm1[k];
            in_front[row] = EMPTY;
        }
        let mut newj = n1;
        for i in 0..nfr {
            let fpivcol = fr_npivcol[i];
            let mut f1rows = 0;
            for _ in 0..fpivcol {
                let j = cperm_init[newj];
                for p in ap[j]..ap[j + 1] {
                    let row = ai[p];
                    if in_front[row] == nfr {
                        in_front[row] = i;
                        f1rows += 1;
                    }
                }
                newj += 1;
            }
            front_1strow[i] = f1rows;
        }
    } else {
        for i in 0..=nfr {
            front_1strow[i] = 0;
        }
        let mut k = 0;
        while k < n1 {
            let row = rperm1[k];
            ci[row] = EMPTY;
            k += 1;
        }
        while k < n_row - nempty_row {
            let row = rperm1[k];
            let i = in_front[k - n1];
            if i != EMPTY {
                front_1strow[i] += 1;
            }
            ci[row] = i;
            k += 1;
        }
        while k < n_row {
            let row = rperm1[k];
            ci[row] = nfr;
            k += 1;
        }
        for row in 0..n_row {
            in_front[row] = ci[row];
        }
    }

    let mut k = n1;
    for i in 0..nfr {
        let fpivcol = fr_npivcol[i];
        k += fpivcol;
        front_npivcol[i] = fpivcol;
        front_parent[i] = fr_parent[i];
    }
    front_npivcol[nfr] = n_col - k;
    front_parent[nfr] = EMPTY;

    for k in 0..n1 {
        rperm_init[k] = rperm1[k];
    }
    let mut k = n1;
    for i in 0..nfr {
        let f1rows = front_1strow[i];
        front_1strow[i] = k;
        k += f1rows;
    }
    front_1strow[nfr] = k;
    let mut f1 = ci;
    for i in 0..=nfr {
        f1[i] = front_1strow[i];
    }
    for row in 0..n_row {
        let i = in_front[row];
        if i != EMPTY {
            let newrow = f1[i];
            f1[i] += 1;
            rperm_init[newrow] = row;
        }
    }
    rperm_init[n_row] = EMPTY;

    if prefer_diagonal {
        let mut dm = vec![0 as Int; n_col as usize + 1];
        for newrow in 0..nn {
            let oldrow = rperm_init[newrow];
            ci[oldrow] = newrow;
        }
        if strategy == UMFPACK_STRATEGY_2BY2 {
            for newcol in 0..nn {
                let oldcol = cperm_init[newcol];
                let oldrow = rperm_2by2[oldcol];
                dm[newcol as usize] = ci[oldrow];
            }
        } else {
            for newcol in 0..nn {
                let oldcol = cperm_init[newcol];
                dm[newcol as usize] = ci[oldcol];
            }
        }
        symbolic.diagonal_map = Some(dm);
    }

    for i in 0..=nfr {
        front_leftmostdesc[i] = EMPTY;
    }
    for i in 0..nfr {
        let mut j = i;
        while j != EMPTY && front_leftmostdesc[j] == EMPTY {
            front_leftmostdesc[j] = i;
            j = front_parent[j];
        }
    }

    let mut maxnrows = 1;
    let mut maxncols = 1;
    let mut dmaxfrsize = 1.0f64;
    let mut nchains = 0;
    chain_start[0] = 0;
    let mut maxrows = 1;
    let mut maxcols = 1;
    for i in 0..nfr {
        let fallrows = fr_nrows[i];
        let fallcols = fr_ncols[i];
        let parent = front_parent[i];
        maxrows = maxrows.max(fallrows);
        maxcols = maxcols.max(fallcols);
        if parent != i + 1 {
            if maxrows % 2 == 0 {
                maxrows += 1;
            }
            chain_maxrows[nchains] = maxrows;
            chain_maxcols[nchains] = maxcols;
            let s = maxrows as f64 * maxcols as f64;
            dmaxfrsize = cmax(dmaxfrsize, s);
            maxnrows = maxnrows.max(maxrows);
            maxncols = maxncols.max(maxcols);
            nchains += 1;
            chain_start[nchains] = i + 1;
            maxrows = 1;
            maxcols = 1;
        }
    }
    dmaxfrsize = dmaxfrsize.ceil();
    symbolic.maxnrows = maxnrows;
    symbolic.maxncols = maxncols;

    if max_rdeg > dense_row_threshold {
        for newrow in 0..n_row {
            let oldrow = rperm_init[newrow];
            ci[oldrow] = newrow;
        }
        for col in n1..n_col - nempty_col {
            let oldcol = cperm_init[col];
            let mut esize = cdeg[oldcol];
            for p in ap[oldcol]..ap[oldcol + 1] {
                let oldrow = ai[p];
                let newrow = ci[oldrow];
                if newrow >= n1 && rdeg[oldrow] > dense_row_threshold {
                    esize -= 1;
                }
            }
            esize_p[col - n1] = esize;
        }
    }

    // permute Cdeg and Rdeg
    for k in 0..n_col {
        ci[k] = cdeg[cperm_init[k]];
    }
    for k in 0..n_col {
        cdeg[k] = ci[k];
    }
    for k in 0..n_row {
        ci[k] = rdeg[rperm_init[k]];
    }
    for k in 0..n_row {
        rdeg[k] = ci[k];
    }

    // simulate UMF_kernel_init
    let mut dlnz = n_inner as f64;
    let mut dunz = dlnz;
    let mut head_usage: Int = 1;
    let mut dhead_usage = 1.0f64;
    let mut tail_usage: Int = 2;
    let mut dtail_usage = 2.0f64;
    tail_usage += units(SIZEOF_PTR, n_row + 1) + units(SIZEOF_PTR, n_row + 1) + 2;
    dtail_usage +=
        dunits(SIZEOF_PTR, (n_row + 1) as f64) + dunits(SIZEOF_PTR, (n_row + 1) as f64) + 2.0;
    for k in 0..n1 {
        let lnz = cdeg[k] - 1;
        let unz = rdeg[k] - 1;
        dlnz += lnz as f64;
        dunz += unz as f64;
        head_usage += units(SIZEOF_INT, lnz)
            + units(SIZEOF_ENTRY, lnz)
            + units(SIZEOF_INT, unz)
            + units(SIZEOF_ENTRY, unz);
        dhead_usage += dunits(SIZEOF_INT, lnz as f64)
            + dunits(SIZEOF_ENTRY, lnz as f64)
            + dunits(SIZEOF_INT, unz as f64)
            + dunits(SIZEOF_ENTRY, unz as f64);
    }
    let has_esize = !esize_p.is_null();
    for k in n1..n_col - nempty_col {
        let esize = if has_esize { esize_p[k - n1] } else { cdeg[k] };
        if esize > 0 {
            tail_usage += get_element_size(esize, 1) + 1;
            dtail_usage += dget_element_size(esize as f64, 1.0) + 1.0;
        }
    }
    if has_esize {
        let mut nrow_elements = 0;
        for k in n1..n_row - nempty_row {
            let rd = rdeg[k];
            if rd > dense_row_threshold {
                tail_usage += get_element_size(1, rd) + 1;
                dtail_usage += (get_element_size(1, rd) + 1) as f64;
                nrow_elements += 1;
            }
        }
        info[UMFPACK_NDENSE_ROW] = nrow_elements as f64;
    }
    if has_esize {
        for row in n1..n_row {
            let rd = rdeg[row];
            let tlen = if rd > dense_row_threshold { 1 } else { rd };
            tail_usage += 1 + units(SIZEOF_TUPLE, tuples(tlen));
            dtail_usage += 1.0 + dunits(SIZEOF_TUPLE, tuples(tlen) as f64);
        }
        let mut col = n1;
        while col < n_col - nempty_col {
            let esize = esize_p[col - n1];
            let tlen = (esize > 0) as Int + (cdeg[col] - esize);
            tail_usage += 1 + units(SIZEOF_TUPLE, tuples(tlen));
            dtail_usage += 1.0 + dunits(SIZEOF_TUPLE, tuples(tlen) as f64);
            col += 1;
        }
        while col < n_col {
            tail_usage += 1 + units(SIZEOF_TUPLE, tuples(0));
            dtail_usage += 1.0 + dunits(SIZEOF_TUPLE, tuples(0) as f64);
            col += 1;
        }
    } else {
        for row in n1..n_row {
            let tlen = rdeg[row];
            tail_usage += 1 + units(SIZEOF_TUPLE, tuples(tlen));
            dtail_usage += 1.0 + dunits(SIZEOF_TUPLE, tuples(tlen) as f64);
        }
        for _col in n1..n_col {
            tail_usage += 1 + units(SIZEOF_TUPLE, tuples(1));
            dtail_usage += 1.0 + dunits(SIZEOF_TUPLE, tuples(1) as f64);
        }
    }
    symbolic.num_mem_init_usage = head_usage + tail_usage;
    let mut dmax_usage = dhead_usage + dtail_usage;
    dmax_usage = cmax(symbolic.num_mem_init_usage as f64, dmax_usage.ceil());
    info[UMFPACK_VARIABLE_INIT_ESTIMATE] = dmax_usage;
    symbolic.dnum_mem_init_usage = dmax_usage;
    let _ = tail_usage - (units(SIZEOF_PTR, n_row + 1) + units(SIZEOF_PTR, n_row + 1));
    dtail_usage -= dunits(SIZEOF_PTR, (n_row + 1) as f64) + dunits(SIZEOF_PTR, (n_row + 1) as f64);

    // simulate UMF_kernel
    let mut link = ci;
    for i in 0..nfr {
        link[i] = EMPTY;
    }
    let mut flops = 0.0f64;
    let nbf = nb as f64;
    for chain in 0..nchains {
        let f1c = chain_start[chain];
        let f2 = chain_start[chain + 1] - 1;
        let dr = chain_maxrows[chain] as f64;
        let dc = chain_maxcols[chain] as f64;
        let fsize = nbf * nbf + dr * nbf + nbf * dc + dr * dc;
        dtail_usage += dunits(SIZEOF_ENTRY, fsize);
        dmax_usage = cmax(dmax_usage, dhead_usage + dtail_usage);
        for i in f1c..=f2 {
            let fpivcol = front_npivcol[i];
            let fallrows = fr_nrows[i];
            let fallcols = fr_ncols[i];
            let parent = front_parent[i];
            let fpiv = fpivcol.min(fallrows);
            let f = fpiv as f64;
            let r = (fallrows - fpiv) as f64;
            let c = (fallcols - fpiv) as f64;
            let mut child = link[i];
            while child != EMPTY {
                let cp = front_npivcol[child].min(fr_nrows[child]) as f64;
                let cr = fr_nrows[child] as f64 - cp;
                let cc = fr_ncols[child] as f64 - cp;
                dtail_usage -= element_size(cr, cc);
                child = link[child];
            }
            flops += DIV_FLOPS * (f * r + (f - 1.0) * f / 2.0)
                + MULTSUB_FLOPS
                    * (f * r * c
                        + (r + c) * (f - 1.0) * f / 2.0
                        + (f - 1.0) * f * (2.0 * f - 1.0) / 6.0);
            let dlf = (f * f - f) / 2.0 + f * r;
            let duf = (f * f - f) / 2.0 + f * c;
            dlnz += dlf;
            dunz += duf;
            dhead_usage += dunits(SIZEOF_ENTRY, dlf + duf) + dunits(SIZEOF_INT, r + c + f);
            if parent != EMPTY {
                dtail_usage += element_size(r, c);
                link[i] = link[parent];
                link[parent] = i;
            }
            dmax_usage = cmax(dmax_usage, dhead_usage + dtail_usage);
        }
        dtail_usage -= dunits(SIZEOF_ENTRY, fsize);
    }
    dhead_usage = dhead_usage.ceil();
    dmax_usage = dmax_usage.ceil();
    symbolic.num_mem_size_est = dhead_usage;
    symbolic.num_mem_usage_est = dmax_usage;
    symbolic.lunz_bound = dlnz + dunz - n_inner as f64;

    umf_set_stats(
        info,
        &symbolic,
        dmax_usage,
        dhead_usage,
        flops,
        dlnz,
        dunz,
        dmaxfrsize,
        n_col as f64,
        n_inner as f64,
        maxnrows as f64,
        maxncols as f64,
        true,
        prefer_diagonal,
        ESTIMATE,
    );
    symbolic.valid = SYMBOLIC_VALID;
    // the timing entries are left at EMPTY (not computed here)
    let _ = (&mut fr_npivcol, &mut ap_v, &mut ai_v, &mut ax_v);
    Ok(symbolic)
}

/// `AMD_aat` on pointer views.
fn amd_aat_p(n: Int, ap: P<Int>, ai: P<Int>, len: P<Int>, tp: P<Int>, info: &mut [f64]) -> Int {
    let nu = n as usize;
    // SAFETY: AMD_aat reads Ap [0..n], Ai [0..Ap[n]], writes Len, Tp [0..n];
    // Len and Tp are disjoint from Ap and Ai.
    unsafe {
        let apl = core::slice::from_raw_parts(ap.r(0) as *const Int, nu + 1);
        let nz = apl[nu] as usize;
        let ail: &[Int] = if nz == 0 {
            &[]
        } else {
            core::slice::from_raw_parts(ai.r(0) as *const Int, nz)
        };
        let (l, t): (&mut [Int], &mut [Int]) = if nu == 0 {
            (&mut [], &mut [])
        } else {
            (
                core::slice::from_raw_parts_mut(len.r(0) as *mut Int, nu),
                core::slice::from_raw_parts_mut(tp.r(0) as *mut Int, nu),
            )
        };
        amd::amd_aat(n, apl, ail, l, t, Some(info))
    }
}
