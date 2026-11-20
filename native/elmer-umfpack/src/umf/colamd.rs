//! `UMF_colamd`: UMFPACK's built-in COLAMD (column approximate minimum
//! degree) that also builds the frontal matrix tree.
//!
//! UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.  The COLAMD
//! code is derived from COLAMD by Stefan I. Larimore and Timothy A. Davis.
//! Translated to Rust.
//!
//! `Colamd_Col` (8 `Int`s) and `Colamd_Row` (6 `Int`s) live at the end of
//! the `A` workspace exactly as in C; fields are reached by offset.

use super::util::{dense_degree_threshold, umf_apply_order, umf_fsize};
use crate::amd::{amd_postorder, Int, EMPTY};
use crate::ptr::P;

pub const COLAMD_KNOBS: usize = 20;
pub const COLAMD_STATS: usize = 20;
pub const COLAMD_DENSE_ROW: usize = 0;
pub const COLAMD_DENSE_COL: usize = 1;
pub const COLAMD_AGGRESSIVE: usize = 2;
pub const COLAMD_DEFRAG_COUNT: usize = 2;
pub const COLAMD_STATUS: usize = 3;
pub const COLAMD_INFO1: usize = 4;
pub const COLAMD_INFO2: usize = 5;
pub const COLAMD_INFO3: usize = 6;
pub const COLAMD_EMPTY_ROW: usize = 7;
pub const COLAMD_EMPTY_COL: usize = 8;
pub const COLAMD_NEWLY_EMPTY_ROW: usize = 9;
pub const COLAMD_NEWLY_EMPTY_COL: usize = 10;
pub const COLAMD_OK: Int = 0;

const COL_INTS: Int = 8;
const ROW_INTS: Int = 6;

/// `UMF_COLAMD_C (n_col)` in `Int`s.
pub fn colamd_c(n_col: Int) -> Int {
    ((n_col as i64 + 1) * (COL_INTS as i64 * 4) / 4) as Int
}
/// `UMF_COLAMD_R (n_row)` in `Int`s.
pub fn colamd_r(n_row: Int) -> Int {
    ((n_row as i64 + 1) * (ROW_INTS as i64 * 4) / 4) as Int
}

/// `UMF_COLAMD_RECOMMENDED (nnz, n_row, n_col)` with `Int` arguments.
pub fn colamd_recommended(nnz: Int, n_row: Int, n_col: Int) -> Int {
    if nnz < 0 || n_row < 0 || n_col < 0 {
        -1
    } else {
        (2 * nnz).max(4 * n_col) + colamd_c(n_col) + colamd_r(n_row) + n_col + nnz / 5
    }
}

/// `UMF_COLAMD_RECOMMENDED` with `double` arguments.
pub fn colamd_recommended_f(nnz: f64, n_row: f64, n_col: f64) -> f64 {
    if nnz < 0.0 || n_row < 0.0 || n_col < 0.0 {
        -1.0
    } else {
        let m = if 2.0 * nnz > 4.0 * n_col {
            2.0 * nnz
        } else {
            4.0 * n_col
        };
        let c = ((n_col + 1.0) * 32.0 / 4.0) as Int as f64;
        let r = ((n_row + 1.0) * 24.0 / 4.0) as Int as f64;
        m + c + r + n_col + nnz / 5.0
    }
}

pub fn umf_colamd_set_defaults(knobs: &mut [f64]) {
    for k in knobs.iter_mut().take(COLAMD_KNOBS) {
        *k = 0.0;
    }
    knobs[COLAMD_DENSE_ROW] = 0.2;
    knobs[COLAMD_DENSE_COL] = 0.2;
    knobs[COLAMD_AGGRESSIVE] = 1.0;
}

// Colamd_Col fields
const C_START: Int = 0;
const C_LENGTH: Int = 1;
const C_THICKNESS: Int = 2; // shared1: thickness / parent
const C_PARENT: Int = 2;
const C_SCORE: Int = 3; // shared2: score / order
const C_ORDER: Int = 3;
const C_HEADHASH: Int = 4; // shared3: headhash / hash / prev
const C_HASH: Int = 4;
const C_PREV: Int = 4;
const C_DEGREE_NEXT: Int = 5; // shared4: degree_next / hash_next
const C_HASH_NEXT: Int = 5;
const C_NEXTCOL: Int = 6;
const C_LASTCOL: Int = 7;
// Colamd_Row fields
const R_START: Int = 0;
const R_LENGTH: Int = 1;
const R_DEGREE: Int = 2; // shared1: degree / p
const R_P: Int = 2;
const R_MARK: Int = 3; // shared2: mark / first_column
const R_FIRST_COLUMN: Int = 3;
const R_THICKNESS: Int = 4;
const R_FRONT: Int = 5;

const ALIVE: Int = 0;
const DEAD: Int = -1;
const DEAD_PRINCIPAL: Int = -1;
const DEAD_NON_PRINCIPAL: Int = -2;

#[derive(Clone, Copy)]
struct Cols(P<Int>);
#[derive(Clone, Copy)]
struct Rows(P<Int>);

impl Cols {
    #[inline(always)]
    fn g(&self, c: Int, f: Int) -> Int {
        self.0[c * COL_INTS + f]
    }
    #[inline(always)]
    fn s(&self, c: Int, f: Int, v: Int) {
        self.0.set(c * COL_INTS + f, v)
    }
    #[inline(always)]
    fn is_dead(&self, c: Int) -> bool {
        self.g(c, C_START) < ALIVE
    }
    #[inline(always)]
    fn is_alive(&self, c: Int) -> bool {
        self.g(c, C_START) >= ALIVE
    }
    #[inline(always)]
    fn kill_principal(&self, c: Int) {
        self.s(c, C_START, DEAD_PRINCIPAL)
    }
    #[inline(always)]
    fn kill_non_principal(&self, c: Int) {
        self.s(c, C_START, DEAD_NON_PRINCIPAL)
    }
}

impl Rows {
    #[inline(always)]
    fn g(&self, r: Int, f: Int) -> Int {
        self.0[r * ROW_INTS + f]
    }
    #[inline(always)]
    fn s(&self, r: Int, f: Int, v: Int) {
        self.0.set(r * ROW_INTS + f, v)
    }
    #[inline(always)]
    fn is_dead(&self, r: Int) -> bool {
        self.g(r, R_MARK) < ALIVE
    }
    #[inline(always)]
    fn is_alive(&self, r: Int) -> bool {
        self.g(r, R_MARK) >= ALIVE
    }
    #[inline(always)]
    fn kill(&self, r: Int) {
        self.s(r, R_MARK, DEAD)
    }
}

/// `UMF_colamd`. `stats` gets `COLAMD_STATS` entries.
#[allow(clippy::too_many_arguments)]
pub fn umf_colamd(
    n_row: Int,
    n_col: Int,
    mut alen: Int,
    mut a: P<Int>,
    mut p: P<Int>,
    knobs: &[f64],
    stats: &mut [Int],
    mut front_npivcol: P<Int>,
    mut front_nrows: P<Int>,
    mut front_ncols: P<Int>,
    mut front_parent: P<Int>,
    mut front_cols: P<Int>,
    p_nfr: &mut Int,
    mut in_front: P<Int>,
) -> bool {
    for s in stats.iter_mut().take(COLAMD_STATS) {
        *s = 0;
    }
    stats[COLAMD_STATUS] = COLAMD_OK;
    stats[COLAMD_INFO1] = -1;
    stats[COLAMD_INFO2] = -1;
    let nnz = p[n_col];
    let aggressive = knobs[COLAMD_AGGRESSIVE] != 0.0;
    let col_size = colamd_c(n_col);
    let row_size = colamd_r(n_row);
    alen -= col_size + row_size;
    let col = Cols(a + alen);
    let row = Rows(a + (alen + col_size));

    init_rows_cols(n_row, n_col, row, col, a, p);
    for c in 0..n_col {
        front_npivcol[c] = 0;
        front_nrows[c] = 0;
        front_ncols[c] = 0;
        front_parent[c] = EMPTY;
        front_cols[c] = EMPTY;
    }
    let sc = init_scoring(n_row, n_col, row, col, a, p, knobs);
    let mut nfr = 0;
    let ngarbage = find_ordering(
        n_row,
        n_col,
        alen,
        row,
        col,
        a,
        p,
        sc.n_col2,
        sc.max_deg,
        2 * nnz,
        front_npivcol,
        front_nrows,
        front_ncols,
        front_parent,
        front_cols,
        &mut nfr,
        aggressive,
        in_front,
    );

    let front_child = a;
    let front_sibling = front_child + nfr;
    let front_stack = front_sibling + nfr;
    let front_order = front_stack + nfr;
    let front_size = front_order + nfr;
    umf_fsize(
        nfr,
        front_size,
        front_nrows,
        front_ncols,
        front_parent,
        front_npivcol,
    );
    {
        // AMD_postorder on P views of the same workspace
        let n = nfr as usize;
        let sl = |q: P<Int>| -> &mut [Int] {
            if n == 0 {
                &mut []
            } else {
                // SAFETY: the five Front_* arrays are disjoint n-long pieces
                // of A (and the caller's arrays), as in the C code.
                unsafe { core::slice::from_raw_parts_mut(q.r(0) as *mut Int, n) }
            }
        };
        amd_postorder(
            nfr,
            sl(front_parent),
            sl(front_npivcol),
            sl(front_size),
            sl(front_order),
            sl(front_child),
            sl(front_sibling),
            sl(front_stack),
        );
    }
    umf_apply_order(front_npivcol, front_order, a, nfr, nfr);
    umf_apply_order(front_nrows, front_order, a, nfr, nfr);
    umf_apply_order(front_ncols, front_order, a, nfr, nfr);
    umf_apply_order(front_parent, front_order, a, nfr, nfr);
    umf_apply_order(front_cols, front_order, a, nfr, nfr);
    for i in 0..nfr {
        let parent = front_parent[i];
        if parent != EMPTY {
            front_parent[i] = front_order[parent];
        }
    }
    for r in 0..n_row {
        let i = in_front[r];
        if i != EMPTY {
            in_front[r] = front_order[i];
        }
    }
    for i in 0..n_col {
        a[i] = EMPTY;
    }
    let mut k = 0;
    for i in 0..nfr {
        let mut c = front_cols[i];
        while c != EMPTY {
            p[k] = c;
            a[c] = k;
            k += 1;
            c = col.g(c, C_NEXTCOL);
        }
    }
    if sc.n_col2 < n_col {
        for c in 0..n_col {
            if a[c] == EMPTY {
                let k = col.g(c, C_ORDER);
                p[k] = c;
                a[c] = k;
            }
        }
    }
    stats[COLAMD_DENSE_ROW] = sc.ndense_row;
    stats[COLAMD_EMPTY_ROW] = sc.nempty_row;
    stats[COLAMD_NEWLY_EMPTY_ROW] = sc.nnewlyempty_row;
    stats[COLAMD_DENSE_COL] = sc.ndense_col;
    stats[COLAMD_EMPTY_COL] = sc.nempty_col;
    stats[COLAMD_NEWLY_EMPTY_COL] = sc.nnewlyempty_col;
    stats[COLAMD_DEFRAG_COUNT] = ngarbage;
    *p_nfr = nfr;
    true
}

fn init_rows_cols(n_row: Int, n_col: Int, row: Rows, col: Cols, mut a: P<Int>, p: P<Int>) {
    for c in 0..n_col {
        col.s(c, C_START, p[c]);
        col.s(c, C_LENGTH, p[c + 1] - p[c]);
        col.s(c, C_THICKNESS, 1);
        col.s(c, C_SCORE, 0);
        col.s(c, C_PREV, EMPTY);
        col.s(c, C_DEGREE_NEXT, EMPTY);
        col.s(c, C_NEXTCOL, EMPTY);
        col.s(c, C_LASTCOL, c);
    }
    for r in 0..n_row {
        row.s(r, R_LENGTH, 0);
        row.s(r, R_THICKNESS, 1);
        row.s(r, R_FRONT, EMPTY);
    }
    for c in 0..n_col {
        for q in p[c]..p[c + 1] {
            let r = a[q];
            row.s(r, R_LENGTH, row.g(r, R_LENGTH) + 1);
        }
    }
    row.s(0, R_START, p[n_col]);
    row.s(0, R_P, row.g(0, R_START));
    for r in 1..n_row {
        row.s(r, R_START, row.g(r - 1, R_START) + row.g(r - 1, R_LENGTH));
        row.s(r, R_P, row.g(r, R_START));
    }
    for c in 0..n_col {
        for q in p[c]..p[c + 1] {
            let r = a[q];
            let rp = row.g(r, R_P);
            row.s(r, R_P, rp + 1);
            a[rp] = c;
        }
    }
    for r in 0..n_row {
        row.s(r, R_MARK, 0);
        row.s(r, R_DEGREE, row.g(r, R_LENGTH));
    }
}

struct Scoring {
    n_col2: Int,
    max_deg: Int,
    ndense_row: Int,
    nempty_row: Int,
    nnewlyempty_row: Int,
    ndense_col: Int,
    nempty_col: Int,
    nnewlyempty_col: Int,
}

fn init_scoring(
    n_row: Int,
    n_col: Int,
    row: Rows,
    col: Cols,
    mut a: P<Int>,
    mut head: P<Int>,
    knobs: &[f64],
) -> Scoring {
    let mut dense_row_count = dense_degree_threshold(knobs[COLAMD_DENSE_ROW], n_col);
    let mut dense_col_count = dense_degree_threshold(knobs[COLAMD_DENSE_COL], n_row);
    dense_row_count = 0.max(dense_row_count.min(n_col));
    dense_col_count = 0.max(dense_col_count.min(n_row));
    let mut max_deg = 0;
    let mut n_col2 = n_col;
    let mut ndense_col = 0;
    let nempty_col = 0;
    let mut nnewlyempty_col = 0;
    let mut ndense_row = 0;
    let nempty_row = 0;

    for c in (0..n_col).rev() {
        let deg = col.g(c, C_LENGTH);
        if deg > dense_col_count {
            n_col2 -= 1;
            col.s(c, C_ORDER, n_col2);
            ndense_col += 1;
            let start = col.g(c, C_START);
            for q in start..start + col.g(c, C_LENGTH) {
                let r = a[q];
                row.s(r, R_DEGREE, row.g(r, R_DEGREE) - 1);
            }
            col.kill_principal(c);
        }
    }
    let mut ne = 0;
    for r in 0..n_row {
        let deg = row.g(r, R_DEGREE);
        if deg > dense_row_count {
            ndense_row += 1;
        }
        if deg == 0 {
            ne += 1;
        }
        if deg > dense_row_count || deg == 0 {
            row.kill(r);
            row.s(r, R_THICKNESS, 0);
        } else {
            max_deg = max_deg.max(deg);
        }
    }
    let nnewlyempty_row = ne - nempty_row;
    for c in (0..n_col).rev() {
        if col.is_dead(c) {
            continue;
        }
        let mut score: Int = 0;
        let start = col.g(c, C_START);
        let mut new_cp = start;
        let cp_end = start + col.g(c, C_LENGTH);
        let mut cp = start;
        while cp < cp_end {
            let r = a[cp];
            cp += 1;
            if row.is_dead(r) {
                continue;
            }
            a[new_cp] = r;
            new_cp += 1;
            score += row.g(r, R_DEGREE) - 1;
            score = score.min(n_col);
        }
        let col_length = new_cp - start;
        if col_length == 0 {
            n_col2 -= 1;
            col.s(c, C_ORDER, n_col2);
            col.kill_principal(c);
            nnewlyempty_col += 1;
        } else {
            col.s(c, C_LENGTH, col_length);
            col.s(c, C_SCORE, score);
        }
    }
    for c in 0..=n_col {
        head[c] = EMPTY;
    }
    for c in (0..n_col).rev() {
        if col.is_alive(c) {
            let score = col.g(c, C_SCORE);
            let next_col = head[score];
            col.s(c, C_PREV, EMPTY);
            col.s(c, C_DEGREE_NEXT, next_col);
            if next_col != EMPTY {
                col.s(next_col, C_PREV, c);
            }
            head[score] = c;
        }
    }
    Scoring {
        n_col2,
        max_deg,
        ndense_row,
        nempty_row,
        nnewlyempty_row,
        ndense_col,
        nempty_col,
        nnewlyempty_col,
    }
}

fn clear_mark(n_row: Int, row: Rows) -> Int {
    for r in 0..n_row {
        if row.is_alive(r) {
            row.s(r, R_MARK, 0);
        }
    }
    1
}

#[allow(clippy::too_many_arguments)]
fn find_ordering(
    n_row: Int,
    n_col: Int,
    alen: Int,
    row: Rows,
    col: Cols,
    mut a: P<Int>,
    mut head: P<Int>,
    n_col2: Int,
    mut max_deg: Int,
    mut pfree: Int,
    mut front_npivcol: P<Int>,
    mut front_nrows: P<Int>,
    mut front_ncols: P<Int>,
    mut front_parent: P<Int>,
    mut front_cols: P<Int>,
    p_nfr: &mut Int,
    aggressive: bool,
    mut in_front: P<Int>,
) -> Int {
    let mut nfr: Int = 0;
    let max_mark = Int::MAX - (2 * n_col + 1);
    let mut tag_mark = clear_mark(n_row, row);
    let mut min_score: Int = 0;
    let mut ngarbage = 0;
    for r in 0..n_row {
        in_front[r] = EMPTY;
    }
    let mut k: Int = 0;
    while k < n_col2 {
        while head[min_score] == EMPTY && min_score < n_col {
            min_score += 1;
        }
        let pivot_col = head[min_score];
        let mut next_col = col.g(pivot_col, C_DEGREE_NEXT);
        head[min_score] = next_col;
        if next_col != EMPTY {
            col.s(next_col, C_PREV, EMPTY);
        }
        let pivot_col_score = col.g(pivot_col, C_SCORE);
        col.s(pivot_col, C_ORDER, k);
        let mut pivot_col_thickness = col.g(pivot_col, C_THICKNESS);
        k += pivot_col_thickness;
        let needed_memory = pivot_col_score.min(n_col - k);
        if pfree + needed_memory >= alen {
            pfree = garbage_collection(n_row, n_col, row, col, a, pfree);
            ngarbage += 1;
            tag_mark = clear_mark(n_row, row);
        }
        let pivot_row_start = pfree;
        let mut pivot_row_degree: Int = 0;
        let mut pivot_row_thickness: Int = 0;
        col.s(pivot_col, C_THICKNESS, -pivot_col_thickness);
        let cstart = col.g(pivot_col, C_START);
        let cend = cstart + col.g(pivot_col, C_LENGTH);
        for cp in cstart..cend {
            let r = a[cp];
            if row.is_dead(r) {
                continue;
            }
            pivot_row_thickness += row.g(r, R_THICKNESS);
            let rstart = row.g(r, R_START);
            let rend = rstart + row.g(r, R_LENGTH);
            for rp in rstart..rend {
                let c = a[rp];
                let col_thickness = col.g(c, C_THICKNESS);
                if col_thickness > 0 && col.is_alive(c) {
                    col.s(c, C_THICKNESS, -col_thickness);
                    a[pfree] = c;
                    pfree += 1;
                    pivot_row_degree += col_thickness;
                }
            }
        }
        col.s(pivot_col, C_THICKNESS, pivot_col_thickness);
        max_deg = max_deg.max(pivot_row_degree);

        let cstart = col.g(pivot_col, C_START);
        let cend = cstart + col.g(pivot_col, C_LENGTH);
        for cp in cstart..cend {
            let r = a[cp];
            if row.is_alive(r) {
                if row.g(r, R_FRONT) != EMPTY {
                    let child = row.g(r, R_FRONT);
                    front_parent[child] = nfr;
                } else {
                    in_front[r] = nfr;
                }
            }
            row.kill(r);
            row.s(r, R_THICKNESS, 0);
        }

        let pivot_row_length = pfree - pivot_row_start;
        let pivot_row = if pivot_row_length > 0 {
            a[col.g(pivot_col, C_START)]
        } else {
            EMPTY
        };

        // compute set differences
        for rp in pivot_row_start..pivot_row_start + pivot_row_length {
            let c = a[rp];
            let col_thickness = -col.g(c, C_THICKNESS);
            col.s(c, C_THICKNESS, col_thickness);
            let cur_score = col.g(c, C_SCORE);
            let prev_col = col.g(c, C_PREV);
            next_col = col.g(c, C_DEGREE_NEXT);
            if prev_col == EMPTY {
                head[cur_score] = next_col;
            } else {
                col.s(prev_col, C_DEGREE_NEXT, next_col);
            }
            if next_col != EMPTY {
                col.s(next_col, C_PREV, prev_col);
            }
            let cstart = col.g(c, C_START);
            let cend = cstart + col.g(c, C_LENGTH);
            for cp in cstart..cend {
                let r = a[cp];
                let row_mark = row.g(r, R_MARK);
                if row_mark < ALIVE {
                    continue;
                }
                let mut set_difference = row_mark - tag_mark;
                if set_difference < 0 {
                    set_difference = row.g(r, R_DEGREE);
                }
                set_difference -= col_thickness;
                if set_difference == 0 && aggressive {
                    if row.g(r, R_FRONT) != EMPTY {
                        let child = row.g(r, R_FRONT);
                        front_parent[child] = nfr;
                    } else {
                        in_front[r] = nfr;
                    }
                    row.kill(r);
                    pivot_row_thickness += row.g(r, R_THICKNESS);
                    row.s(r, R_THICKNESS, 0);
                } else {
                    row.s(r, R_MARK, set_difference + tag_mark);
                }
            }
        }

        // add up set differences for each column
        for rp in pivot_row_start..pivot_row_start + pivot_row_length {
            let c = a[rp];
            let mut hash: u32 = 0;
            let mut cur_score: Int = 0;
            let cstart = col.g(c, C_START);
            let mut new_cp = cstart;
            let cend = cstart + col.g(c, C_LENGTH);
            for cp in cstart..cend {
                let r = a[cp];
                let row_mark = row.g(r, R_MARK);
                if row_mark < ALIVE {
                    continue;
                }
                a[new_cp] = r;
                new_cp += 1;
                hash = hash.wrapping_add(r as u32);
                cur_score += row_mark - tag_mark;
                cur_score = cur_score.min(n_col);
            }
            col.s(c, C_LENGTH, new_cp - cstart);
            if col.g(c, C_LENGTH) == 0 {
                col.kill_principal(c);
                pivot_row_degree -= col.g(c, C_THICKNESS);
                col.s(c, C_ORDER, k);
                k += col.g(c, C_THICKNESS);
                pivot_col_thickness += col.g(c, C_THICKNESS);
                col.s(col.g(c, C_LASTCOL), C_NEXTCOL, front_cols[nfr]);
                front_cols[nfr] = c;
            } else {
                col.s(c, C_SCORE, cur_score);
                hash %= (n_col + 1) as u32;
                let head_column = head[hash as Int];
                let first_col;
                if head_column > EMPTY {
                    first_col = col.g(head_column, C_HEADHASH);
                    col.s(head_column, C_HEADHASH, c);
                } else {
                    first_col = -(head_column + 2);
                    head[hash as Int] = -(c + 2);
                }
                col.s(c, C_HASH_NEXT, first_col);
                col.s(c, C_HASH, hash as Int);
            }
        }

        detect_super_cols(col, a, head, pivot_row_start, pivot_row_length);

        col.kill_principal(pivot_col);
        col.s(col.g(pivot_col, C_LASTCOL), C_NEXTCOL, front_cols[nfr]);
        front_cols[nfr] = pivot_col;

        tag_mark += max_deg + 1;
        if tag_mark >= max_mark {
            tag_mark = clear_mark(n_row, row);
        }

        // finalize the new pivot row, and column scores
        let mut new_rp = pivot_row_start;
        for rp in pivot_row_start..pivot_row_start + pivot_row_length {
            let c = a[rp];
            if col.is_dead(c) {
                continue;
            }
            a[new_rp] = c;
            new_rp += 1;
            let len = col.g(c, C_LENGTH);
            a[col.g(c, C_START) + len] = pivot_row;
            col.s(c, C_LENGTH, len + 1);
            let mut cur_score = col.g(c, C_SCORE) + pivot_row_degree;
            let max_score = n_col - k - col.g(c, C_THICKNESS);
            cur_score -= col.g(c, C_THICKNESS);
            cur_score = cur_score.min(max_score);
            col.s(c, C_SCORE, cur_score);
            next_col = head[cur_score];
            col.s(c, C_DEGREE_NEXT, next_col);
            col.s(c, C_PREV, EMPTY);
            if next_col != EMPTY {
                col.s(next_col, C_PREV, c);
            }
            head[cur_score] = c;
            min_score = min_score.min(cur_score);
        }

        front_npivcol[nfr] = pivot_col_thickness;
        front_nrows[nfr] = pivot_row_thickness;
        front_ncols[nfr] = pivot_col_thickness + pivot_row_degree;
        front_parent[nfr] = EMPTY;
        pivot_row_thickness -= pivot_col_thickness;
        pivot_row_thickness = 0.max(pivot_row_thickness);
        if pivot_row_degree > 0 && pivot_row_thickness > 0 {
            row.s(pivot_row, R_START, pivot_row_start);
            row.s(pivot_row, R_LENGTH, new_rp - pivot_row_start);
            row.s(pivot_row, R_DEGREE, pivot_row_degree);
            row.s(pivot_row, R_MARK, 0);
            row.s(pivot_row, R_THICKNESS, pivot_row_thickness);
            row.s(pivot_row, R_FRONT, nfr);
        }
        nfr += 1;
    }
    *p_nfr = nfr;
    ngarbage
}

fn detect_super_cols(col: Cols, a: P<Int>, mut head: P<Int>, row_start: Int, row_length: Int) {
    for rp in row_start..row_start + row_length {
        let c0 = a[rp];
        if col.is_dead(c0) {
            continue;
        }
        let hash = col.g(c0, C_HASH);
        let head_column = head[hash];
        let first_col = if head_column > EMPTY {
            col.g(head_column, C_HEADHASH)
        } else {
            -(head_column + 2)
        };
        let mut super_c = first_col;
        while super_c != EMPTY {
            let length = col.g(super_c, C_LENGTH);
            let mut prev_c = super_c;
            let mut c = col.g(super_c, C_HASH_NEXT);
            while c != EMPTY {
                if col.g(c, C_LENGTH) != length || col.g(c, C_SCORE) != col.g(super_c, C_SCORE) {
                    prev_c = c;
                    c = col.g(c, C_HASH_NEXT);
                    continue;
                }
                let cp1 = col.g(super_c, C_START);
                let cp2 = col.g(c, C_START);
                let mut i = 0;
                while i < length {
                    if a[cp1 + i] != a[cp2 + i] {
                        break;
                    }
                    i += 1;
                }
                if i != length {
                    prev_c = c;
                    c = col.g(c, C_HASH_NEXT);
                    continue;
                }
                col.s(
                    super_c,
                    C_THICKNESS,
                    col.g(super_c, C_THICKNESS) + col.g(c, C_THICKNESS),
                );
                col.s(c, C_PARENT, super_c);
                col.kill_non_principal(c);
                col.s(c, C_ORDER, EMPTY);
                col.s(prev_c, C_HASH_NEXT, col.g(c, C_HASH_NEXT));
                col.s(col.g(super_c, C_LASTCOL), C_NEXTCOL, c);
                col.s(super_c, C_LASTCOL, col.g(c, C_LASTCOL));
                c = col.g(c, C_HASH_NEXT);
            }
            super_c = col.g(super_c, C_HASH_NEXT);
        }
        if head_column > EMPTY {
            col.s(head_column, C_HEADHASH, EMPTY);
        } else {
            head[hash] = EMPTY;
        }
    }
}

fn garbage_collection(
    n_row: Int,
    n_col: Int,
    row: Rows,
    col: Cols,
    mut a: P<Int>,
    pfree: Int,
) -> Int {
    let mut pdest: Int = 0;
    for c in 0..n_col {
        if col.is_alive(c) {
            let mut psrc = col.g(c, C_START);
            col.s(c, C_START, pdest);
            let length = col.g(c, C_LENGTH);
            for _ in 0..length {
                let r = a[psrc];
                psrc += 1;
                if row.is_alive(r) {
                    a[pdest] = r;
                    pdest += 1;
                }
            }
            col.s(c, C_LENGTH, pdest - col.g(c, C_START));
        }
    }
    for r in 0..n_row {
        if row.is_alive(r) {
            if row.g(r, R_LENGTH) == 0 {
                row.kill(r);
            } else {
                let psrc = row.g(r, R_START);
                row.s(r, R_FIRST_COLUMN, a[psrc]);
                a[psrc] = -r - 1; // ONES_COMPLEMENT
            }
        }
    }
    let mut psrc = pdest;
    while psrc < pfree {
        let v = a[psrc];
        psrc += 1;
        if v < 0 {
            psrc -= 1;
            let r = -a[psrc] - 1;
            a[psrc] = row.g(r, R_FIRST_COLUMN);
            row.s(r, R_START, pdest);
            let length = row.g(r, R_LENGTH);
            for _ in 0..length {
                let c = a[psrc];
                psrc += 1;
                if col.is_alive(c) {
                    a[pdest] = c;
                    pdest += 1;
                }
            }
            row.s(r, R_LENGTH, pdest - row.g(r, R_START));
        }
    }
    pdest
}
