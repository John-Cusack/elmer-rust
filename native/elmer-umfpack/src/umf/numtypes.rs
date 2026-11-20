//! `NumericType`, `WorkType`, `Element`, `Tuple`, `Unit` and the memory
//! manager of the numeric factorisation (`UMF_mem_*`, `UMF_tuple_lengths`,
//! `UMF_build_tuples`).
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  Translated
//! to Rust.

use super::consts::*;
use super::symbolic::{dget_element_size, get_element_size, tuples};
use super::util::*;
use crate::amd::{Int, EMPTY};
use crate::ptr::P;

pub type Entry = f64;

/// `Unit`: the 8-byte allocation unit of `Numeric->Memory` (a block header
/// or `Align`).
#[repr(C, align(8))]
#[derive(Clone, Copy, Default, Debug)]
pub struct Unit {
    pub size: Int,
    pub prevsize: Int,
}

/// `Element`.
#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct Element {
    pub cdeg: Int,
    pub rdeg: Int,
    pub nrowsleft: Int,
    pub ncolsleft: Int,
    pub nrows: Int,
    pub ncols: Int,
    pub next: Int,
}

/// `Tuple`.
#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct Tuple {
    pub e: Int,
    pub f: Int,
}

pub const UNITS_ELEMENT: Int = 4; // UNITS (Element, 1)

/// `NumericType`.
#[derive(Default)]
pub struct Numeric {
    pub flops: f64,
    pub relpt: f64,
    pub relpt2: f64,
    pub droptol: f64,
    pub alloc_init: f64,
    pub front_alloc_init: f64,
    pub rsmin: f64,
    pub rsmax: f64,
    pub min_udiag: f64,
    pub max_udiag: f64,
    pub rcond: f64,
    pub scale: Int,
    pub valid: Int,
    pub memory: Vec<Unit>,
    pub ihead: Int,
    pub itail: Int,
    pub ibig: Int,
    pub size: Int,
    pub rperm: Vec<Int>,
    pub cperm: Vec<Int>,
    pub upos: Vec<Int>,
    pub lpos: Vec<Int>,
    pub lip: Vec<Int>,
    pub lilen: Vec<Int>,
    pub uip: Vec<Int>,
    pub uilen: Vec<Int>,
    pub upattern: Vec<Int>,
    pub ulen: Int,
    pub npiv: Int,
    pub nnzpiv: Int,
    pub d: Vec<Entry>,
    pub do_recip: bool,
    pub rs: Vec<f64>,
    pub n_row: Int,
    pub n_col: Int,
    pub n1: Int,
    pub tail_usage: Int,
    pub init_usage: Int,
    pub max_usage: Int,
    pub ngarbage: Int,
    pub nrealloc: Int,
    pub ncostly: Int,
    pub isize: Int,
    pub n_lentries: Int,
    pub n_uentries: Int,
    pub lnz: Int,
    pub all_lnz: Int,
    pub unz: Int,
    pub all_unz: Int,
    pub maxfrsize: Int,
    pub maxnrows: Int,
    pub maxncols: Int,
}

impl Numeric {
    /// `Numeric->Memory` (re-derive after anything that may reallocate it).
    #[inline]
    pub fn mem(&mut self) -> P<Unit> {
        P::new(&mut self.memory)
    }
    pub fn rperm_p(&mut self) -> P<Int> {
        P::new(&mut self.rperm)
    }
    pub fn cperm_p(&mut self) -> P<Int> {
        P::new(&mut self.cperm)
    }
    pub fn upos_p(&mut self) -> P<Int> {
        P::new(&mut self.upos)
    }
    pub fn lpos_p(&mut self) -> P<Int> {
        P::new(&mut self.lpos)
    }
    pub fn lip_p(&mut self) -> P<Int> {
        P::new(&mut self.lip)
    }
    pub fn lilen_p(&mut self) -> P<Int> {
        P::new(&mut self.lilen)
    }
    pub fn uip_p(&mut self) -> P<Int> {
        P::new(&mut self.uip)
    }
    pub fn uilen_p(&mut self) -> P<Int> {
        P::new(&mut self.uilen)
    }
    pub fn d_p(&mut self) -> P<Entry> {
        P::new(&mut self.d)
    }
    pub fn rs_p(&mut self) -> P<f64> {
        if self.rs.is_empty() {
            P::null()
        } else {
            P::new(&mut self.rs)
        }
    }
}

/// The arrays `WorkType` points to.
#[derive(Default)]
pub struct WorkStore {
    pub e: Vec<Int>,
    pub wx: Vec<Entry>,
    pub wy: Vec<Entry>,
    pub wp: Vec<Int>,
    pub wrp: Vec<Int>,
    pub wm: Vec<Int>,
    pub wio: Vec<Int>,
    pub woi: Vec<Int>,
    pub woo: Vec<Int>,
    pub lpattern: Vec<Int>,
    pub upattern: Vec<Int>,
    pub diagonal_map: Vec<Int>,
    pub diagonal_imap: Vec<Int>,
    pub front_new1strow: Vec<Int>,
    pub frows: Vec<Int>,
    pub fcols: Vec<Int>,
    pub frpos: Vec<Int>,
    pub fcpos: Vec<Int>,
}

/// `WorkType` (pointer fields are views into [`WorkStore`] or `Memory`).
pub struct Work {
    pub e: P<Int>,
    pub wx: P<Entry>,
    pub wy: P<Entry>,
    pub wp: P<Int>,
    pub wrp: P<Int>,
    pub wm: P<Int>,
    pub wio: P<Int>,
    pub woi: P<Int>,
    pub woo: P<Int>,
    pub wrow: P<Int>,
    pub new_rows: P<Int>,
    pub new_cols: P<Int>,
    pub lpattern: P<Int>,
    pub upattern: P<Int>,
    pub ulen: Int,
    pub llen: Int,
    pub diagonal_map: P<Int>,
    pub diagonal_imap: P<Int>,
    pub n_row: Int,
    pub n_col: Int,
    pub nz: Int,
    pub n1: Int,
    pub elen: Int,
    pub npiv: Int,
    pub ndiscard: Int,
    pub wrpflag: Int,
    pub nel: Int,
    pub noff_diagonal: Int,
    pub prior_element: Int,
    pub rdeg0: Int,
    pub cdeg0: Int,
    pub rrdeg: Int,
    pub ccdeg: Int,
    pub candidates: [Int; MAX_CANDIDATES],
    pub n_candidates: Int,
    pub ksuper: Int,
    pub firstsuper: Int,
    pub jsuper: Int,
    pub ncand: Int,
    pub nextcand: Int,
    pub lo: Int,
    pub hi: Int,
    pub pivrow: Int,
    pub pivcol: Int,
    pub do_extend: bool,
    pub do_update: bool,
    pub nforced: Int,
    pub any_skip: bool,
    pub do_scan2row: bool,
    pub do_scan2col: bool,
    pub do_grow: bool,
    pub pivot_case: Int,
    pub frontid: Int,
    pub nfr: Int,
    pub front_new1strow: P<Int>,
    pub pivrow_arr: [Int; MAXNB as usize],
    pub pivcol_arr: [Int; MAXNB as usize],
    pub flublock: P<Entry>,
    pub flblock: P<Entry>,
    pub fublock: P<Entry>,
    pub fcblock: P<Entry>,
    pub frows: P<Int>,
    pub fcols: P<Int>,
    pub frpos: P<Int>,
    pub fcpos: P<Int>,
    pub fnrows: Int,
    pub fncols: Int,
    pub fnr_curr: Int,
    pub fnc_curr: Int,
    pub fcurr_size: Int,
    pub fnrows_max: Int,
    pub fncols_max: Int,
    pub nb: Int,
    pub fnpiv: Int,
    pub fnzeros: Int,
    pub fscan_row: Int,
    pub fscan_col: Int,
    pub fnrows_new: Int,
    pub fncols_new: Int,
    pub pivrow_in_front: bool,
    pub pivcol_in_front: bool,
}

impl Default for Work {
    fn default() -> Self {
        Work {
            e: P::null(),
            wx: P::null(),
            wy: P::null(),
            wp: P::null(),
            wrp: P::null(),
            wm: P::null(),
            wio: P::null(),
            woi: P::null(),
            woo: P::null(),
            wrow: P::null(),
            new_rows: P::null(),
            new_cols: P::null(),
            lpattern: P::null(),
            upattern: P::null(),
            ulen: 0,
            llen: 0,
            diagonal_map: P::null(),
            diagonal_imap: P::null(),
            n_row: 0,
            n_col: 0,
            nz: 0,
            n1: 0,
            elen: 0,
            npiv: 0,
            ndiscard: 0,
            wrpflag: 0,
            nel: 0,
            noff_diagonal: 0,
            prior_element: 0,
            rdeg0: 0,
            cdeg0: 0,
            rrdeg: 0,
            ccdeg: 0,
            candidates: [0; MAX_CANDIDATES],
            n_candidates: 0,
            ksuper: 0,
            firstsuper: 0,
            jsuper: 0,
            ncand: 0,
            nextcand: 0,
            lo: 0,
            hi: 0,
            pivrow: 0,
            pivcol: 0,
            do_extend: false,
            do_update: false,
            nforced: 0,
            any_skip: false,
            do_scan2row: false,
            do_scan2col: false,
            do_grow: false,
            pivot_case: 0,
            frontid: 0,
            nfr: 0,
            front_new1strow: P::null(),
            pivrow_arr: [0; MAXNB as usize],
            pivcol_arr: [0; MAXNB as usize],
            flublock: P::null(),
            flblock: P::null(),
            fublock: P::null(),
            fcblock: P::null(),
            frows: P::null(),
            fcols: P::null(),
            frpos: P::null(),
            fcpos: P::null(),
            fnrows: 0,
            fncols: 0,
            fnr_curr: 0,
            fnc_curr: 0,
            fcurr_size: 0,
            fnrows_max: 0,
            fncols_max: 0,
            nb: 0,
            fnpiv: 0,
            fnzeros: 0,
            fscan_row: 0,
            fscan_col: 0,
            fnrows_new: 0,
            fncols_new: 0,
            pivrow_in_front: false,
            pivcol_in_front: false,
        }
    }
}

/// `GET_ELEMENT` on the element at unit `p` of Memory.
#[derive(Clone, Copy)]
pub struct ElemView {
    pub ep: P<Element>,
    pub cols: P<Int>,
    pub rows: P<Int>,
    pub c: P<Entry>,
    pub ncols: Int,
    pub nrows: Int,
}

#[inline]
pub fn get_element(p: P<Unit>) -> ElemView {
    let ep: P<Element> = p.cast();
    let q = p + UNITS_ELEMENT;
    let cols: P<Int> = q.cast();
    let ncols = ep.r(0).ncols;
    let rows = cols + ncols;
    let nrows = ep.r(0).nrows;
    let q = q + units(SIZEOF_INT, ncols + nrows);
    ElemView {
        ep,
        cols,
        rows,
        c: q.cast(),
        ncols,
        nrows,
    }
}

/// `UMF_mem_init_memoryspace`.
pub fn umf_mem_init_memoryspace(numeric: &mut Numeric) {
    numeric.ngarbage = 0;
    numeric.nrealloc = 0;
    numeric.ncostly = 0;
    numeric.ibig = EMPTY;
    numeric.ihead = 0;
    numeric.itail = numeric.size;
    numeric.itail -= 2;
    let m = numeric.mem();
    let p = m.r(numeric.itail);
    numeric.tail_usage = 2;
    p.prevsize = 0;
    p.size = 1;
    numeric.ihead += 1;
    numeric.max_usage = 3;
    numeric.init_usage = numeric.max_usage;
}

/// `UMF_mem_alloc_head_block`.
pub fn umf_mem_alloc_head_block(numeric: &mut Numeric, nunits: Int) -> Int {
    if nunits > numeric.itail - numeric.ihead {
        return 0;
    }
    let p = numeric.ihead;
    numeric.ihead += nunits;
    let usage = numeric.ihead + numeric.tail_usage;
    numeric.max_usage = numeric.max_usage.max(usage);
    p
}

/// `UMF_mem_alloc_tail_block`.
pub fn umf_mem_alloc_tail_block(numeric: &mut Numeric, nunits: Int) -> Int {
    let m = numeric.mem();
    let mut bigsize: Int = 0;
    let mut pbig: Option<Int> = None;
    if numeric.ibig != EMPTY {
        pbig = Some(numeric.ibig);
        bigsize = -m.r(numeric.ibig).size;
    }
    let p: Int;
    if let Some(pb) = pbig.filter(|_| bigsize >= nunits) {
        p = pb;
        let pnext = p + 1 + bigsize;
        bigsize -= nunits + 1;
        if bigsize < 4 {
            m.r(p).size = -m.r(p).size;
            numeric.ibig = EMPTY;
        } else {
            m.r(p).size = nunits;
            numeric.ibig += nunits + 1;
            let pb = m.r(numeric.ibig);
            pb.size = -bigsize;
            pb.prevsize = nunits;
            m.r(pnext).prevsize = bigsize;
        }
    } else {
        let pnext = numeric.itail;
        if (nunits + 1) > (numeric.itail - numeric.ihead) {
            return 0;
        }
        numeric.itail -= nunits + 1;
        p = numeric.itail;
        m.r(p).size = nunits;
        m.r(p).prevsize = 0;
        m.r(pnext).prevsize = nunits;
    }
    numeric.tail_usage += m.r(p).size + 1;
    let usage = numeric.ihead + numeric.tail_usage;
    numeric.max_usage = numeric.max_usage.max(usage);
    p + 1
}

/// `UMF_mem_free_tail_block`.
pub fn umf_mem_free_tail_block(numeric: &mut Numeric, i: Int) {
    if i == EMPTY || i == 0 {
        return;
    }
    let m = numeric.mem();
    let mut p = i - 1;
    numeric.tail_usage -= m.r(p).size + 1;
    let pnext = p + 1 + m.r(p).size;
    if m.r(pnext).size < 0 {
        m.r(p).size += -m.r(pnext).size + 1;
    }
    if p > numeric.itail {
        let pprev = p - 1 - m.r(p).prevsize;
        let sprev = m.r(pprev).size;
        if sprev < 0 {
            m.r(pprev).size = m.r(p).size + (-sprev) + 1;
            p = pprev;
        }
    }
    let pnext = p + 1 + m.r(p).size;
    if p == numeric.itail {
        numeric.itail = pnext;
        m.r(pnext).prevsize = 0;
        if numeric.ibig != EMPTY && numeric.ibig <= numeric.itail {
            numeric.ibig = EMPTY;
        }
    } else {
        if numeric.ibig == EMPTY {
            numeric.ibig = p;
        } else {
            let pbig = numeric.ibig;
            if -m.r(pbig).size < m.r(p).size {
                numeric.ibig = p;
            }
        }
        m.r(pnext).prevsize = m.r(p).size;
        m.r(p).size = -m.r(p).size;
    }
}

/// `UMF_mem_alloc_element`. Returns the unit index of the element (0 on
/// failure) and its view; `size` gets the block size including the header.
pub fn umf_mem_alloc_element(
    numeric: &mut Numeric,
    nrows: Int,
    ncols: Int,
    size: &mut Int,
) -> (Int, Option<ElemView>) {
    *size = get_element_size(nrows, ncols);
    if int_overflow(dget_element_size(nrows as f64, ncols as f64) + 1.0) {
        return (0, None);
    }
    let i = umf_mem_alloc_tail_block(numeric, *size);
    *size += 1;
    if i == 0 {
        return (0, None);
    }
    let p = numeric.mem() + i;
    let ep: P<Element> = p.cast();
    let q = p + UNITS_ELEMENT;
    let cols: P<Int> = q.cast();
    let rows = cols + ncols;
    let q = q + units(SIZEOF_INT, ncols + nrows);
    let c: P<Entry> = q.cast();
    let e = ep.r(0);
    e.nrows = nrows;
    e.ncols = ncols;
    e.nrowsleft = nrows;
    e.ncolsleft = ncols;
    e.cdeg = 0;
    e.rdeg = 0;
    e.next = EMPTY;
    (
        i,
        Some(ElemView {
            ep,
            cols,
            rows,
            c,
            ncols,
            nrows,
        }),
    )
}

/// `UMF_tuple_lengths`.
pub fn umf_tuple_lengths(numeric: &mut Numeric, work: &Work, p_dusage: &mut f64) -> Int {
    let e_arr = work.e;
    let row_degree = numeric.rperm_p();
    let col_degree = numeric.cperm_p();
    let mut row_tlen = numeric.uilen_p();
    let mut col_tlen = numeric.lilen_p();
    let n_row = work.n_row;
    let n_col = work.n_col;
    let n1 = work.n1;
    let nel = work.nel;
    let m = numeric.mem();
    for e in 1..=nel {
        if e_arr[e] != 0 {
            let v = get_element(m + e_arr[e]);
            for i in 0..v.nrows {
                let row = v.rows[i];
                if row >= n1 {
                    row_tlen[row] += 1;
                }
            }
            for i in 0..v.ncols {
                let col = v.cols[i];
                if col >= n1 {
                    col_tlen[col] += 1;
                }
            }
        }
    }
    let mut usage: Int = 0;
    let mut dusage = 0.0f64;
    for col in n1..n_col {
        if col_degree[col] >= 0 {
            usage += 1 + units(SIZEOF_TUPLE, tuples(col_tlen[col]));
            dusage += 1.0 + dunits(SIZEOF_TUPLE, tuples(col_tlen[col]) as f64);
        }
    }
    for row in n1..n_row {
        if row_degree[row] >= 0 {
            usage += 1 + units(SIZEOF_TUPLE, tuples(row_tlen[row]));
            dusage += 1.0 + dunits(SIZEOF_TUPLE, tuples(row_tlen[row]) as f64);
        }
    }
    *p_dusage = dusage;
    usage
}

/// `UMF_build_tuples`.
pub fn umf_build_tuples(numeric: &mut Numeric, work: &Work) -> bool {
    let e_arr = work.e;
    let col_degree = numeric.cperm_p();
    let row_degree = numeric.rperm_p();
    let mut row_tuples = numeric.uip_p();
    let mut row_tlen = numeric.uilen_p();
    let mut col_tuples = numeric.lip_p();
    let mut col_tlen = numeric.lilen_p();
    let n_row = work.n_row;
    let n_col = work.n_col;
    let nel = work.nel;
    let n1 = work.n1;
    for row in n1..n_row {
        if row_degree[row] >= 0 {
            row_tuples[row] =
                umf_mem_alloc_tail_block(numeric, units(SIZEOF_TUPLE, tuples(row_tlen[row])));
            if row_tuples[row] == 0 {
                return false;
            }
            row_tlen[row] = 0;
        }
    }
    let mut col = n_col - 1;
    while col >= n1 {
        if col_degree[col] >= 0 {
            col_tuples[col] =
                umf_mem_alloc_tail_block(numeric, units(SIZEOF_TUPLE, tuples(col_tlen[col])));
            if col_tuples[col] == 0 {
                return false;
            }
            col_tlen[col] = 0;
        }
        col -= 1;
    }
    let m = numeric.mem();
    for e in 1..=nel {
        let v = get_element(m + e_arr[e]);
        for f in 0..v.ncols {
            let col = v.cols[f];
            let tp: P<Tuple> = (m + col_tuples[col]).cast();
            let k = col_tlen[col];
            col_tlen[col] += 1;
            *tp.r(k) = Tuple { e, f };
        }
        for f in 0..v.nrows {
            let row = v.rows[f];
            let tp: P<Tuple> = (m + row_tuples[row]).cast();
            let k = row_tlen[row];
            row_tlen[row] += 1;
            *tp.r(k) = Tuple { e, f };
        }
    }
    true
}

/// `NON_PIVOTAL_COL` / `NON_PIVOTAL_ROW` on the degree arrays.
#[inline]
pub fn non_pivotal(deg: P<Int>, i: Int) -> bool {
    deg[i] >= 0
}
