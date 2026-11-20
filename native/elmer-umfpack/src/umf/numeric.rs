//! `UMFPACK_numeric`: the numeric LU factorisation.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::kernel::umf_kernel;
use super::numtypes::*;
use super::symbolic::Symbolic;
use super::util::*;
use crate::amd::{Int, EMPTY};
use crate::ptr::P;

fn vi(n: Int) -> Vec<Int> {
    vec![0; n.max(1) as usize]
}
fn vd(n: Int) -> Vec<f64> {
    vec![0.0; n.max(1) as usize]
}

/// `UMFPACK_numeric`. Returns the Numeric object and the status
/// (`UMFPACK_OK` or `UMFPACK_WARNING_singular_matrix`), or the error status.
pub fn umfpack_numeric(
    ap_s: &[Int],
    ai_s: &[Int],
    ax_s: &[f64],
    symbolic: &Symbolic,
    control: Option<&[f64]>,
    user_info: Option<&mut [f64]>,
) -> Result<(Numeric, Int), Int> {
    let mut relpt = get_control(
        control,
        UMFPACK_PIVOT_TOLERANCE,
        UMFPACK_DEFAULT_PIVOT_TOLERANCE,
    );
    let mut relpt2 = get_control(
        control,
        UMFPACK_SYM_PIVOT_TOLERANCE,
        UMFPACK_DEFAULT_SYM_PIVOT_TOLERANCE,
    );
    let mut alloc_init = get_control(control, UMFPACK_ALLOC_INIT, UMFPACK_DEFAULT_ALLOC_INIT);
    let mut front_alloc_init = get_control(
        control,
        UMFPACK_FRONT_ALLOC_INIT,
        UMFPACK_DEFAULT_FRONT_ALLOC_INIT,
    );
    let mut scale = get_control(control, UMFPACK_SCALE, UMFPACK_DEFAULT_SCALE as f64) as Int;
    let mut droptol = get_control(control, UMFPACK_DROPTOL, UMFPACK_DEFAULT_DROPTOL);
    relpt = cmax(0.0, cmin(relpt, 1.0));
    relpt2 = cmax(0.0, cmin(relpt2, 1.0));
    droptol = cmax(0.0, droptol);
    front_alloc_init = cmin(1.0, front_alloc_init);
    if scale != UMFPACK_SCALE_NONE && scale != UMFPACK_SCALE_MAX {
        scale = UMFPACK_DEFAULT_SCALE;
    }
    let mut info2 = [0.0f64; UMFPACK_INFO];
    let info: &mut [f64] = match user_info {
        Some(i) => {
            for k in UMFPACK_NUMERIC_SIZE..=UMFPACK_MAX_FRONT_NCOLS {
                i[k] = EMPTY as f64;
            }
            for k in UMFPACK_NUMERIC_DEFRAG..UMFPACK_IR_TAKEN {
                i[k] = EMPTY as f64;
            }
            i
        }
        None => {
            for v in info2.iter_mut() {
                *v = EMPTY as f64;
            }
            &mut info2
        }
    };
    if symbolic.valid != SYMBOLIC_VALID {
        info[UMFPACK_STATUS] = UMFPACK_ERROR_invalid_Symbolic_object as f64;
        return Err(UMFPACK_ERROR_invalid_Symbolic_object);
    }
    if symbolic.ordering == UMFPACK_ORDERING_AMD && alloc_init >= 0.0 {
        alloc_init = (symbolic.nz as f64 + symbolic.amd_lunz) / symbolic.lunz_bound;
        alloc_init = cmin(1.0, alloc_init);
        alloc_init *= UMF_REALLOC_INCREASE;
    }
    let n_row = symbolic.n_row;
    let n_col = symbolic.n_col;
    let n_inner = n_row.min(n_col);
    if int_overflow(symbolic.dnum_mem_init_usage * SIZEOF_UNIT as f64) {
        info[UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory as f64;
        return Err(UMFPACK_ERROR_out_of_memory);
    }
    info[UMFPACK_STATUS] = UMFPACK_OK as f64;
    info[UMFPACK_NROW] = n_row as f64;
    info[UMFPACK_NCOL] = n_col as f64;
    info[UMFPACK_SIZE_OF_UNIT] = SIZEOF_UNIT as f64;
    info[UMFPACK_NZ] = ap_s[n_col as usize] as f64;
    if droptol > 0.0 {
        // UMF_store_lu_drop is not ported (Elmer never sets a drop tolerance)
        info[UMFPACK_STATUS] = UMFPACK_ERROR_internal_error as f64;
        return Err(UMFPACK_ERROR_internal_error);
    }

    // work_alloc
    let nn = n_row.max(n_col);
    let nfr = symbolic.nfr;
    let n1 = symbolic.n1;
    let maxnrows = n_row.min(symbolic.maxnrows + symbolic.nb);
    let maxncols = n_col.min(symbolic.maxncols + symbolic.nb);
    let maxnrc = maxnrows.max(maxncols);
    let elen = (n_col - n1) + (n_row - n1) + (n_col - n1).min(n_row - n1) + 1;
    let mut store = WorkStore {
        wx: vd(maxnrows + 1),
        wy: vd(maxnrows + 1),
        frpos: vi(n_row + 1),
        lpattern: vi(n_row + 1),
        fcpos: vi(n_col + 1),
        wp: vi(nn + 1),
        wrp: vi(n_col.max(maxnrows) + 1),
        frows: vi(maxnrows + 1),
        wm: vi(maxnrows + 1),
        fcols: vi(maxncols + 1),
        wio: vi(maxncols + 1),
        woi: vi(maxncols + 1),
        woo: vi(maxnrc + 1),
        e: vi(elen),
        front_new1strow: vi(nfr + 1),
        diagonal_map: if symbolic.prefer_diagonal {
            vi(nn)
        } else {
            Vec::new()
        },
        diagonal_imap: if symbolic.prefer_diagonal {
            vi(nn)
        } else {
            Vec::new()
        },
        upattern: vi(n_col + 1),
    };
    let opt = |v: &mut Vec<Int>| if v.is_empty() { P::null() } else { P::new(v) };
    let mut work = Work {
        n_row,
        n_col,
        nfr,
        nb: symbolic.nb,
        n1,
        elen,
        wx: P::new(&mut store.wx),
        wy: P::new(&mut store.wy),
        frpos: P::new(&mut store.frpos),
        lpattern: P::new(&mut store.lpattern),
        fcpos: P::new(&mut store.fcpos),
        wp: P::new(&mut store.wp),
        wrp: P::new(&mut store.wrp),
        frows: P::new(&mut store.frows),
        wm: P::new(&mut store.wm),
        fcols: P::new(&mut store.fcols),
        wio: P::new(&mut store.wio),
        woi: P::new(&mut store.woi),
        woo: P::new(&mut store.woo),
        e: P::new(&mut store.e),
        front_new1strow: P::new(&mut store.front_new1strow),
        diagonal_map: opt(&mut store.diagonal_map),
        diagonal_imap: opt(&mut store.diagonal_imap),
        upattern: P::new(&mut store.upattern),
        ..Default::default()
    };

    // numeric_alloc
    let mut numeric = Numeric {
        valid: 0,
        d: vd(n_inner + 1),
        rperm: vi(n_row + 1),
        cperm: vi(n_col + 1),
        lpos: vi(n_row + 1),
        lilen: vi(n_col + 1),
        lip: vi(n_col + 1),
        upos: vi(n_col + 1),
        uilen: vi(n_row + 1),
        uip: vi(n_row + 1),
        rs: if scale != UMFPACK_SCALE_NONE {
            vd(n_row)
        } else {
            Vec::new()
        },
        ..Default::default()
    };
    let mut nsize = if alloc_init < 0.0 {
        -alloc_init
    } else {
        alloc_init * symbolic.num_mem_usage_est + 1.0
    };
    let min_usage = symbolic.num_mem_init_usage;
    nsize = cmax(min_usage as f64, nsize);
    let bsize = Int::MAX as f64 / SIZEOF_UNIT as f64 - 1.0;
    nsize = cmin(nsize, bsize);
    numeric.size = nsize as Int;
    numeric.memory = vec![Unit::default(); numeric.size.max(1) as usize];

    numeric.relpt = relpt;
    numeric.relpt2 = relpt2;
    numeric.droptol = droptol;
    numeric.alloc_init = alloc_init;
    numeric.front_alloc_init = front_alloc_init;
    numeric.scale = scale;

    // the input matrix, read through C-style pointers
    let mut ap_v = ap_s[..=n_col as usize].to_vec();
    let nz = ap_v[n_col as usize].max(0) as usize;
    let mut ai_v = if nz == 0 {
        vec![0]
    } else {
        ai_s[..nz].to_vec()
    };
    let mut ax_v = if nz == 0 {
        vec![0.0]
    } else {
        ax_s[..nz].to_vec()
    };
    let status = umf_kernel(
        P::new(&mut ap_v),
        P::new(&mut ai_v),
        P::new(&mut ax_v),
        &mut numeric,
        &mut work,
        &mut store,
        symbolic,
    );
    info[UMFPACK_STATUS] = status as f64;
    if status < UMFPACK_OK {
        return Err(status);
    }
    info[UMFPACK_FORCED_UPDATES] = work.nforced as f64;
    info[UMFPACK_VARIABLE_INIT] = numeric.init_usage as f64;
    if symbolic.prefer_diagonal {
        info[UMFPACK_NOFF_DIAG] = work.noff_diagonal as f64;
    }
    let npiv = numeric.npiv;
    let ulen = numeric.ulen;
    drop(store);

    // reduce the size of the Numeric object (UMF_realloc keeps the prefix)
    if npiv < n_row {
        numeric.lpos.truncate(npiv.max(0) as usize + 1);
        numeric.uilen.truncate(npiv.max(0) as usize + 1);
        numeric.uip.truncate(npiv.max(0) as usize + 1);
    }
    if npiv < n_col {
        numeric.upos.truncate(npiv.max(0) as usize + 1);
        numeric.lilen.truncate(npiv.max(0) as usize + 1);
        numeric.lip.truncate(npiv.max(0) as usize + 1);
    }
    if ulen > 0 && ulen < n_col {
        numeric.upattern.truncate(ulen as usize + 1);
    }
    let newsize = numeric.ihead;
    if newsize < numeric.size {
        numeric.memory.truncate(newsize.max(1) as usize);
        numeric.size = newsize;
    }
    numeric.ihead = numeric.size;
    numeric.itail = numeric.ihead;
    numeric.tail_usage = 0;
    numeric.ibig = EMPTY;

    umf_set_stats(
        info,
        symbolic,
        numeric.max_usage as f64,
        numeric.size as f64,
        numeric.flops,
        numeric.lnz as f64 + n_inner as f64,
        numeric.unz as f64 + numeric.nnzpiv as f64,
        numeric.maxfrsize as f64,
        ulen as f64,
        npiv as f64,
        numeric.maxnrows as f64,
        numeric.maxncols as f64,
        scale != UMFPACK_SCALE_NONE,
        symbolic.prefer_diagonal,
        ACTUAL,
    );
    info[UMFPACK_ALLOC_INIT_USED] = numeric.alloc_init;
    info[UMFPACK_NUMERIC_DEFRAG] = numeric.ngarbage as f64;
    info[UMFPACK_NUMERIC_REALLOC] = numeric.nrealloc as f64;
    info[UMFPACK_NUMERIC_COSTLY_REALLOC] = numeric.ncostly as f64;
    info[UMFPACK_COMPRESSED_PATTERN] = numeric.isize as f64;
    info[UMFPACK_LU_ENTRIES] = (numeric.n_lentries + numeric.n_uentries + numeric.npiv) as f64;
    info[UMFPACK_UDIAG_NZ] = numeric.nnzpiv as f64;
    info[UMFPACK_RSMIN] = numeric.rsmin;
    info[UMFPACK_RSMAX] = numeric.rsmax;
    info[UMFPACK_WAS_SCALED] = numeric.scale as f64;
    info[UMFPACK_ALL_LNZ] = (numeric.all_lnz + n_inner) as f64;
    info[UMFPACK_ALL_UNZ] = (numeric.all_unz + numeric.nnzpiv) as f64;
    info[UMFPACK_NZDROPPED] =
        ((numeric.all_lnz - numeric.lnz) + (numeric.all_unz - numeric.unz)) as f64;
    if numeric.min_udiag == 0.0
        || numeric.max_udiag == 0.0
        || numeric.min_udiag.is_nan()
        || numeric.max_udiag.is_nan()
    {
        numeric.rcond = 0.0;
    } else {
        numeric.rcond = numeric.min_udiag / numeric.max_udiag;
    }
    info[UMFPACK_UMIN] = numeric.min_udiag;
    info[UMFPACK_UMAX] = numeric.max_udiag;
    info[UMFPACK_RCOND] = numeric.rcond;
    let mut status = UMFPACK_OK;
    if numeric.nnzpiv < n_inner || numeric.rcond == 0.0 || numeric.rcond.is_nan() {
        status = UMFPACK_WARNING_singular_matrix;
        info[UMFPACK_STATUS] = status as f64;
    }
    numeric.valid = NUMERIC_VALID;
    Ok((numeric, status))
}
