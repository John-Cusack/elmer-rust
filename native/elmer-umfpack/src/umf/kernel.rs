//! `UMF_kernel`, `UMF_kernel_init`, `UMF_kernel_wrapup`.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::assemble::umf_assemble;
use super::consts::*;
use super::dense::{sabs, umf_blas3_update, umf_scale, umf_scale_column};
use super::front::*;
use super::local_search::umf_local_search;
use super::numtypes::*;
use super::store_lu::umf_store_lu;
use super::symbolic::Symbolic;
use super::util::*;
use crate::amd::{amd_valid, Int, EMPTY};
use crate::ptr::P;

/// `UMF_kernel`. `drop` (droptol > 0) is not supported: Elmer never sets it.
pub fn umf_kernel(
    ap: P<Int>,
    ai: P<Int>,
    ax: P<f64>,
    numeric: &mut Numeric,
    work: &mut Work,
    store: &mut WorkStore,
    symbolic: &Symbolic,
) -> Int {
    if !umf_kernel_init(ap, ai, ax, numeric, work, symbolic) {
        return UMFPACK_ERROR_different_pattern;
    }
    let nchains = symbolic.nchains;
    let nb = symbolic.nb;
    let fixq = symbolic.fix_q;
    macro_rules! tryo {
        ($e:expr) => {
            if !$e {
                return UMFPACK_ERROR_out_of_memory;
            }
        };
    }
    for chain in 0..nchains {
        let f1 = symbolic.chain_start[chain as usize];
        let f2 = symbolic.chain_start[chain as usize + 1] - 1;
        tryo!(umf_start_front(chain, numeric, work, symbolic));
        work.frontid = f1;
        while work.frontid <= f2 {
            work.ncand = symbolic.front_npivcol[work.frontid as usize];
            work.lo = work.nextcand;
            work.hi = work.nextcand + work.ncand - 1;
            let mut jmax = (MAX_CANDIDATES as Int).min(work.ncand);
            if fixq {
                jmax = 1;
            }
            for j in 0..jmax {
                work.candidates[j as usize] = work.nextcand;
                work.nextcand += 1;
            }
            work.n_candidates = jmax;
            while work.ncand > 0 {
                let status = umf_local_search(numeric, work, symbolic);
                if status == UMFPACK_ERROR_different_pattern {
                    return UMFPACK_ERROR_different_pattern;
                }
                if status == UMFPACK_WARNING_singular_matrix {
                    continue;
                }
                if work.do_update {
                    umf_blas3_update(work);
                    tryo!(umf_store_lu(numeric, work));
                }
                if work.do_extend {
                    tryo!(umf_extend_front(numeric, work));
                } else {
                    tryo!(umf_create_element(numeric, work, symbolic));
                    tryo!(umf_init_front(numeric, work));
                }
                umf_assemble(numeric, work, fixq);
                umf_scale_column(numeric, work);
                let evaporate = work.fnrows == 0 || work.fncols == 0;
                if work.fnpiv >= nb || evaporate {
                    umf_blas3_update(work);
                    tryo!(umf_store_lu(numeric, work));
                }
                work.pivrow_in_front = false;
                work.pivcol_in_front = false;
                if evaporate {
                    let _ = umf_create_element(numeric, work, symbolic);
                    work.fnrows = 0;
                    work.fncols = 0;
                }
            }
            work.frontid += 1;
        }
        umf_blas3_update(work);
        tryo!(umf_store_lu(numeric, work));
        work.fnrows_new = work.fnrows;
        work.fncols_new = work.fncols;
        tryo!(umf_create_element(numeric, work, symbolic));
        work.fnrows = 0;
        work.fncols = 0;
    }
    umf_kernel_wrapup(numeric, symbolic, work, store);
    UMFPACK_OK
}

/// `packsp` of `umf_kernel_init` (no dropping).
fn packsp(mut pnew: Int, p_p: &mut Int, p_len: &mut Int, m: P<Unit>) -> Int {
    let mut p = *p_p;
    let len = *p_len;
    let mut bi: P<Int> = (m + p).cast();
    p += units(SIZEOF_INT, len);
    let mut bx: P<Entry> = (m + p).cast();
    let mut len_new = 0;
    for p in 0..len {
        let i = bi[p];
        let x = bx[p];
        if x == 0.0 {
            continue;
        }
        if len_new != p {
            bi[len_new] = i;
            bx[len_new] = x;
        }
        len_new += 1;
    }
    *p_p = pnew;
    *p_len = len_new;
    let mut bi2: P<Int> = (m + pnew).cast();
    pnew += units(SIZEOF_INT, len_new);
    let mut bx2: P<Entry> = (m + pnew).cast();
    pnew += units(SIZEOF_ENTRY, len_new);
    for p in 0..len_new {
        bi2[p] = bi[p];
    }
    for p in 0..len_new {
        bx2[p] = bx[p];
    }
    pnew
}

/// `UMF_kernel_init`.
pub fn umf_kernel_init(
    ap: P<Int>,
    ai: P<Int>,
    ax: P<f64>,
    numeric: &mut Numeric,
    work: &mut Work,
    symbolic: &Symbolic,
) -> bool {
    let n_row = symbolic.n_row;
    let n_col = symbolic.n_col;
    let nn = n_row.max(n_col);
    let n_inner = n_row.min(n_col);
    let nempty_col = symbolic.nempty_col;
    let nempty_row = symbolic.nempty_row;
    let nempty = nempty_row.min(nempty_col);
    let cperm_init = &symbolic.cperm_init;
    let rperm_init = &symbolic.rperm_init;
    let cdeg = &symbolic.cdeg;
    let rdeg = &symbolic.rdeg;
    let n1 = symbolic.n1;
    let dense_row_threshold = symbolic.dense_row_threshold;
    work.nforced = 0;
    work.ndiscard = 0;
    work.noff_diagonal = 0;
    let nz = ap[n_col];
    if nz < 0 || ap[0] != 0 || nz != symbolic.nz {
        return false;
    }
    let prefer_diagonal = symbolic.prefer_diagonal;
    let mut diagonal_map = work.diagonal_map;
    let mut diagonal_imap = work.diagonal_imap;

    umf_mem_init_memoryspace(numeric);

    work.fnpiv = 0;
    work.fncols = 0;
    work.fnrows = 0;
    work.fncols_max = 0;
    work.fnrows_max = 0;
    work.fnzeros = 0;
    work.fcurr_size = 0;
    work.fnr_curr = 0;
    work.fnc_curr = 0;
    work.nz = nz;
    work.prior_element = EMPTY;
    work.ulen = 0;
    work.llen = 0;
    work.npiv = n1;
    work.frontid = 0;
    work.nextcand = n1;

    let m = numeric.mem();
    let mut rperm = numeric.rperm_p();
    let mut cperm = numeric.cperm_p();
    let mut row_degree = rperm;
    let mut col_degree = cperm;
    let mut row_tlen = numeric.uilen_p();
    let mut col_tlen = numeric.lilen_p();
    let mut lip = numeric.lip_p();
    let mut uip = numeric.uip_p();
    let mut lilen = numeric.lilen_p();
    let mut uilen = numeric.uilen_p();
    let mut frpos = work.frpos;
    let mut fcpos = work.fcpos;
    let mut wp = work.wp;
    let mut wrp = work.wrp;
    let mut d = numeric.d_p();
    let mut upos = numeric.upos_p();
    let mut lpos = numeric.lpos_p();
    for k in 0..n_inner {
        d[k] = 0.0;
    }
    let mut rs = numeric.rs_p();
    for row in 0..=n_row {
        lpos[row] = EMPTY;
        row_tlen[row] = 0;
    }
    for col in 0..=n_col {
        upos[col] = EMPTY;
        col_tlen[col] = 0;
        fcpos[col] = EMPTY;
        wrp[col] = 0;
    }
    work.wrpflag = 1;
    for i in 0..=nn {
        wp[i] = EMPTY;
    }
    work.cdeg0 = 1;
    work.rdeg0 = 1;
    let fixq = symbolic.fix_q;
    let mut e_arr = work.e;
    numeric.n_row = n_row;
    numeric.n_col = n_col;
    numeric.npiv = 0;
    numeric.nnzpiv = 0;
    numeric.min_udiag = 0.0;
    numeric.max_udiag = 0.0;
    numeric.rcond = 0.0;
    numeric.isize = 0;
    numeric.n_lentries = 0;
    numeric.n_uentries = 0;
    numeric.lnz = 0;
    numeric.unz = 0;
    numeric.all_lnz = 0;
    numeric.all_unz = 0;
    numeric.maxfrsize = 0;
    numeric.maxnrows = 0;
    numeric.maxncols = 0;
    numeric.flops = 0.0;
    numeric.n1 = n1;
    let do_scale = numeric.scale != UMFPACK_SCALE_NONE;
    let rsmin;
    let rsmax;
    let do_recip;
    if do_scale {
        let do_max = numeric.scale == UMFPACK_SCALE_MAX;
        for row in 0..n_row {
            rs[row] = 0.0;
        }
        for col in 0..n_col {
            let mut ilast = EMPTY;
            let p1 = ap[col];
            let p2 = ap[col + 1];
            if p1 > p2 {
                return false;
            }
            for p in p1..p2 {
                let row = ai[p];
                if row <= ilast || row >= n_row {
                    return false;
                }
                let value = sabs(ax[p]);
                let r = rs[row];
                if !r.is_nan() {
                    if value.is_nan() {
                        rs[row] = value;
                    } else if do_max {
                        rs[row] = cmax(r, value);
                    } else {
                        rs[row] += value;
                    }
                }
                ilast = row;
            }
        }
        for row in 0..n_row {
            let r = rs[row];
            if r == 0.0 || r.is_nan() {
                rs[row] = 1.0;
            }
        }
        let mut mn = rs[0];
        let mut mx = rs[0];
        for row in 0..n_row {
            mn = cmin(mn, rs[row]);
            mx = cmax(mx, rs[row]);
        }
        rsmin = mn;
        rsmax = mx;
        do_recip = rsmin >= RECIPROCAL_TOLERANCE;
        if do_recip {
            for row in 0..n_row {
                rs[row] = 1.0 / rs[row];
            }
        }
    } else {
        rsmin = -1.0;
        rsmax = -1.0;
        do_recip = false;
        // SAFETY: ap/ai cover n_col+1 and nz entries
        let ok = unsafe {
            let apl = core::slice::from_raw_parts(ap.r(0) as *const Int, n_col as usize + 1);
            let ail: &[Int] = if nz == 0 {
                &[]
            } else {
                core::slice::from_raw_parts(ai.r(0) as *const Int, nz as usize)
            };
            amd_valid(n_row, n_col, apl, ail)
        };
        if !ok {
            return false;
        }
    }
    numeric.rsmin = rsmin;
    numeric.rsmax = rsmax;
    numeric.do_recip = do_recip;

    for newrow in 0..n_row {
        let oldrow = rperm_init[newrow as usize];
        frpos[oldrow] = newrow;
    }
    if prefer_diagonal {
        for i in 0..nn {
            diagonal_map[i] = EMPTY;
            diagonal_imap[i] = EMPTY;
        }
        let sdm = symbolic.diagonal_map.as_ref().expect("Diagonal_map");
        for k in n1..nn - nempty {
            let newrow = sdm[k as usize];
            diagonal_map[k] = newrow;
            diagonal_imap[newrow] = k;
        }
    }

    // row-form pointers (the C code keeps these in two tail blocks)
    let rpi_blk = umf_mem_alloc_tail_block(numeric, units(SIZEOF_PTR, n_row + 1));
    let rpx_blk = umf_mem_alloc_tail_block(numeric, units(SIZEOF_PTR, n_row + 1));
    if rpi_blk == 0 || rpx_blk == 0 {
        return false;
    }
    let mut rpi: Vec<P<Int>> = vec![P::null(); n_row as usize + 1];
    let mut rpx: Vec<P<Entry>> = vec![P::null(); n_row as usize + 1];

    for k in 0..n1 {
        let lnz = cdeg[k as usize] - 1;
        let unz = rdeg[k as usize] - 1;
        let size = units(SIZEOF_INT, lnz)
            + units(SIZEOF_ENTRY, lnz)
            + units(SIZEOF_INT, unz)
            + units(SIZEOF_ENTRY, unz);
        let mut p = umf_mem_alloc_head_block(numeric, size);
        if p == 0 {
            return false;
        }
        numeric.all_lnz += lnz;
        numeric.all_unz += unz;
        let lipv = p;
        p += units(SIZEOF_INT, lnz);
        p += units(SIZEOF_ENTRY, lnz);
        let uipv = p;
        rpi[k as usize] = (m + p).cast();
        p += units(SIZEOF_INT, unz);
        rpx[k as usize] = (m + p).cast();
        lip[k] = lipv;
        lilen[k] = lnz;
        uip[k] = uipv;
        uilen[k] = unz;
        wp[k] = unz;
        let k1 = -k - 1;
        rperm[k] = k1;
        cperm[k] = k1;
    }

    let mut e: Int = 0;
    e_arr[e] = 0;
    work.flublock = P::null();
    work.flblock = P::null();
    work.fublock = P::null();
    work.fcblock = P::null();
    let esize_v = symbolic.esize_vec.as_deref();
    let mut empty_elements = false;
    for k in n1..n_col - nempty_col {
        e = k - n1 + 1;
        let esize = match esize_v {
            Some(v) => v[(k - n1) as usize],
            None => cdeg[k as usize],
        };
        if esize > 0 {
            let mut size = 0;
            let (eloc, view) = umf_mem_alloc_element(numeric, esize, 1, &mut size);
            e_arr[e] = eloc;
            if e_arr[e] <= 0 {
                return false;
            }
            view.unwrap().cols.set(0, k);
        } else {
            e_arr[e] = 0;
            empty_elements = true;
        }
    }
    if esize_v.is_some() {
        for k in n1..n_row - nempty_row {
            let rd = rdeg[k as usize];
            if rd > dense_row_threshold {
                e += 1;
                let mut size = 0;
                let (eloc, view) = umf_mem_alloc_element(numeric, 1, rd, &mut size);
                e_arr[e] = eloc;
                if e_arr[e] <= 0 {
                    return false;
                }
                let v = view.unwrap();
                v.rows.set(0, k);
                rpi[k as usize] = v.cols;
                rpx[k as usize] = v.c;
                wp[k] = rd;
            }
        }
    }
    work.nel = e;

    // the singleton columns
    for k in 0..n1 {
        let pivcol = cperm_init[k as usize];
        let p2 = ap[pivcol + 1];
        let mut p = lip[k];
        let mut li: P<Int> = (m + p).cast();
        let lilenk = lilen[k];
        p += units(SIZEOF_INT, lilenk);
        let lval: P<Entry> = (m + p).cast();
        let mut lval_m = lval;
        let mut llen = 0;
        for pa in ap[pivcol]..p2 {
            let oldrow = ai[pa];
            let newrow = frpos[oldrow];
            let mut x = ax[pa];
            if do_scale {
                if do_recip {
                    x *= rs[oldrow];
                } else {
                    x /= rs[oldrow];
                }
            }
            if newrow == k {
                d[k] = x;
            } else if newrow < k {
                wp[newrow] -= 1;
                if wp[newrow] < 0 {
                    return false;
                }
                let r = &mut rpi[newrow as usize];
                r.set(0, k);
                *r = *r + 1;
                let r = &mut rpx[newrow as usize];
                r.set(0, x);
                *r = *r + 1;
            } else {
                if llen >= lilenk {
                    return false;
                }
                li[llen] = newrow;
                lval_m[llen] = x;
                llen += 1;
            }
        }
        if llen != lilenk {
            return false;
        }
        if llen > 0 {
            let pivot_value = d[k];
            umf_scale(llen, pivot_value, lval);
        }
    }

    // the remaining columns
    for k in n1..n_col {
        let oldcol = cperm_init[k as usize];
        let p2 = ap[oldcol + 1];
        let cd = cdeg[k as usize];
        col_degree[k] = if fixq { 0 } else { cd };
        e = k - n1 + 1;
        let esize;
        let mut rows: P<Int> = P::null();
        let mut c: P<Entry> = P::null();
        if k < n_col - nempty_col {
            esize = match esize_v {
                Some(v) => v[(k - n1) as usize],
                None => cd,
            };
            if e_arr[e] != 0 {
                let v = get_element(m + e_arr[e]);
                rows = v.rows;
                c = v.c;
            }
        } else {
            esize = 0;
        }
        let mut clen = 0;
        for pa in ap[oldcol]..p2 {
            let oldrow = ai[pa];
            let newrow = frpos[oldrow];
            let mut x = ax[pa];
            if do_scale {
                if do_recip {
                    x *= rs[oldrow];
                } else {
                    x /= rs[oldrow];
                }
            }
            let rd = rdeg[newrow as usize];
            if newrow < n1 || rd > dense_row_threshold {
                wp[newrow] -= 1;
                if wp[newrow] < 0 {
                    return false;
                }
                let r = &mut rpi[newrow as usize];
                r.set(0, k);
                *r = *r + 1;
                let r = &mut rpx[newrow as usize];
                r.set(0, x);
                *r = *r + 1;
            } else {
                if clen >= esize {
                    return false;
                }
                rows[clen] = newrow;
                c[clen] = x;
                clen += 1;
            }
        }
        if clen != esize {
            return false;
        }
    }

    umf_mem_free_tail_block(numeric, rpi_blk);
    umf_mem_free_tail_block(numeric, rpx_blk);

    if n1 > 0 {
        let mut pnew = lip[0];
        for k in 0..n1 {
            pnew = packsp(pnew, lip.r(k), lilen.r(k), m);
            numeric.lnz += lilen[k];
            pnew = packsp(pnew, uip.r(k), uilen.r(k), m);
            numeric.unz += uilen[k];
        }
        numeric.ihead = pnew;
    }
    for k in 0..n1 {
        if wp[k] != 0 {
            return false;
        }
    }
    for k in n1..n_row {
        let rd = rdeg[k as usize];
        row_degree[k] = rd;
        if rd > dense_row_threshold && wp[k] != 0 {
            return false;
        }
    }
    col_degree[n_col] = 0;

    if empty_elements {
        let mut e2 = 0;
        for e in 1..=work.nel {
            if e_arr[e] != 0 {
                e2 += 1;
                e_arr[e2] = e_arr[e];
            }
        }
        work.nel = e2;
    }
    for e in work.nel + 1..work.elen {
        e_arr[e] = 0;
    }
    for row in 0..=n_row {
        frpos[row] = EMPTY;
    }
    for i in 0..=nn {
        wp[i] = EMPTY;
    }
    let mut unused = 0.0;
    let _ = umf_tuple_lengths(numeric, work, &mut unused);
    if !umf_build_tuples(numeric, work) {
        return false;
    }
    numeric.init_usage = numeric.max_usage;
    let mut fnew = work.front_new1strow;
    for i in 0..=symbolic.nfr {
        fnew[i] = symbolic.front_1strow[i as usize];
    }
    let _ = (&mut row_degree, &mut col_degree);
    true
}

/// `UMF_kernel_wrapup`.
pub fn umf_kernel_wrapup(
    numeric: &mut Numeric,
    symbolic: &Symbolic,
    work: &mut Work,
    store: &mut WorkStore,
) {
    let n_row = symbolic.n_row;
    let n_col = symbolic.n_col;
    let n_inner = n_row.min(n_col);
    let mut rperm = numeric.rperm_p();
    let mut cperm = numeric.cperm_p();
    let mut lilen = numeric.lilen_p();
    let mut uilen = numeric.uilen_p();
    let mut upos = numeric.upos_p();
    let mut lpos = numeric.lpos_p();
    let mut lip = numeric.lip_p();
    let mut uip = numeric.uip_p();
    let mut d = numeric.d_p();
    let npiv = work.npiv;
    numeric.npiv = npiv;
    numeric.ulen = work.ulen;
    for k in 0..npiv {
        let dd = sabs(d[k]);
        let zero_pivot = dd == 0.0;
        let nan_pivot = dd.is_nan();
        if !zero_pivot {
            numeric.nnzpiv += 1;
        }
        if k == 0 {
            numeric.min_udiag = dd;
            numeric.max_udiag = dd;
        } else {
            if numeric.min_udiag != 0.0 {
                if zero_pivot || nan_pivot {
                    numeric.min_udiag = dd;
                } else if !numeric.min_udiag.is_nan() {
                    numeric.min_udiag = cmin(numeric.min_udiag, dd);
                }
            }
            if nan_pivot {
                numeric.max_udiag = dd;
            } else if !numeric.max_udiag.is_nan() {
                numeric.max_udiag = cmax(numeric.max_udiag, dd);
            }
        }
    }
    if npiv < n_row {
        let mut k = npiv;
        for row in 0..n_row {
            if rperm[row] >= 0 {
                rperm[row] = -k - 1;
                lpos[row] = EMPTY;
                uip[row] = EMPTY;
                uilen[row] = 0;
                k += 1;
            }
        }
    }
    if npiv < n_col {
        let mut k = npiv;
        for col in 0..n_col {
            if cperm[col] >= 0 {
                cperm[col] = -k - 1;
                upos[col] = EMPTY;
                lip[col] = EMPTY;
                lilen[col] = 0;
                k += 1;
            }
        }
    }
    if npiv < n_inner {
        for k in npiv..n_inner {
            d[k] = 0.0;
        }
    }
    if numeric.ulen > 0 {
        numeric.upattern = core::mem::take(&mut store.upattern);
        work.upattern = P::null();
    }
    if numeric.nnzpiv < n_inner && !numeric.min_udiag.is_nan() {
        numeric.min_udiag = 0.0;
    }
    let mut frpos = work.frpos;
    let mut fcpos = work.fcpos;
    let mut wp = work.wp;
    let mut wr = work.lpattern;
    let mut wc = work.wrp;
    for pivrow in 0..n_row {
        let k = -rperm[pivrow] - 1;
        wp[k] = pivrow;
        frpos[pivrow] = k;
    }
    for k in 0..n_row {
        rperm[k] = wp[k];
    }
    for pivcol in 0..n_col {
        let k = -cperm[pivcol] - 1;
        wp[k] = pivcol;
        fcpos[pivcol] = k;
    }
    for k in 0..n_col {
        cperm[k] = wp[k];
    }
    for k in 0..npiv {
        let pivrow = rperm[k];
        wr[k] = uilen[pivrow];
        wp[k] = uip[pivrow];
    }
    for k in 0..npiv {
        uilen[k] = wr[k];
        uip[k] = wp[k];
    }
    for k in 0..npiv {
        let pivrow = rperm[k];
        wp[k] = lpos[pivrow];
    }
    for k in 0..npiv {
        lpos[k] = wp[k];
    }
    for k in 0..npiv {
        let pivcol = cperm[k];
        wc[k] = lilen[pivcol];
        wp[k] = lip[pivcol];
    }
    for k in 0..npiv {
        lilen[k] = wc[k];
        lip[k] = wp[k];
    }
    for k in 0..npiv {
        let pivcol = cperm[k];
        wp[k] = upos[pivcol];
    }
    for k in 0..npiv {
        upos[k] = wp[k];
    }
    upos[npiv] = EMPTY;
    lpos[npiv] = EMPTY;
    uip[npiv] = EMPTY;
    lip[npiv] = EMPTY;
    uilen[npiv] = 0;
    lilen[npiv] = 0;
    let n1 = symbolic.n1;
    let m = numeric.mem();
    for k in 0..n1 {
        let ulen = uilen[k];
        if ulen > 0 {
            let mut ip: P<Int> = (m + uip[k]).cast();
            for i in 0..ulen {
                ip[i] = fcpos[ip[i]];
            }
        }
    }
    for k in n1..npiv {
        let mut up = uip[k];
        if up < 0 {
            let ulen = uilen[k];
            if ulen > 0 {
                up = -up;
                let mut ip: P<Int> = (m + up).cast();
                for i in 0..ulen {
                    ip[i] = fcpos[ip[i]];
                }
            }
        }
    }
    let ulen = numeric.ulen;
    if ulen > 0 {
        for i in 0..ulen as usize {
            let col = numeric.upattern[i];
            numeric.upattern[i] = fcpos[col];
        }
    }
    for k in 0..n1 {
        let llen = lilen[k];
        if llen > 0 {
            let mut ip: P<Int> = (m + lip[k]).cast();
            for i in 0..llen {
                ip[i] = frpos[ip[i]];
            }
        }
    }
    for k in n1..npiv {
        let llen = lilen[k];
        if llen > 0 {
            let mut lp = lip[k];
            if lp < 0 {
                lp = -lp;
            }
            let mut ip: P<Int> = (m + lp).cast();
            for i in 0..llen {
                ip[i] = frpos[ip[i]];
            }
        }
    }
    for k in 0..n_row {
        rperm[k] = symbolic.rperm_init[rperm[k] as usize];
    }
    for k in 0..n_col {
        cperm[k] = symbolic.cperm_init[cperm[k] as usize];
    }
}
