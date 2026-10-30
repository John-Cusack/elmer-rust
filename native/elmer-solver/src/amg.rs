//! WIP, not wired: algebraic multigrid (`MG Method = algebraic`,
//! `AMGSolve`, Multigrid.F90:1484-5019), ported routine by routine.
//!
//! Done: `AMGBonds` (the strong connections, Multigrid.F90:2223).
//! Still to port before `Linear System Preconditioning = multigrid` /
//! `Linear System Solver = Multigrid` with `MG Method = algebraic` can use
//! it: `AMGCoarse` (the C/F splitting), `InterpolateF2C` (the projector),
//! `CRS_ProjectMatrixCreate` (the coarse matrix), `AMGSolve`/`AMGSweep` (the
//! cycle, which can reuse `multigrid::Hierarchy`) and the multigrid
//! preconditioner entry (`MultiGridPrec`).
#![allow(dead_code)]

use elmer_core::{CsrMatrix, ParamMap, Params};

/// The `MG*` keywords `AMGBonds` reads.
#[derive(Clone, Debug)]
pub(crate) struct BondOptions {
    /// `MG Strong Connection Limit` (0.06).
    pub neg_lim: f64,
    /// `MG Positive Connection Allow` and `MG Positive Connection Limit` (1).
    pub allow_pos: bool,
    pub pos_lim: f64,
    /// `MG Eliminate Dirichlet` (true on the finest level) and `MG
    /// Eliminate Dirichlet Limit` (1e-8).
    pub eliminate_dirichlet: bool,
    pub dir_lim: f64,
    /// `MG Strong Connection Maximum` / `Minimum` (0: none).
    pub max_conns: usize,
    pub min_conns: usize,
}

impl BondOptions {
    pub(crate) fn from_params(p: &ParamMap, finest: bool) -> Self {
        Self {
            neg_lim: p.real("mg strong connection limit").unwrap_or(0.06),
            allow_pos: p.logical("mg positive connection allow") == Some(true),
            pos_lim: p.real("mg positive connection limit").unwrap_or(1.0),
            eliminate_dirichlet: finest && p.logical("mg eliminate dirichlet").unwrap_or(true),
            dir_lim: p.real("mg eliminate dirichlet limit").unwrap_or(1.0e-8),
            max_conns: p.integer("mg strong connection maximum").map_or(0, |k| k.max(0) as usize),
            min_conns: p.integer("mg strong connection minimum").map_or(0, |k| k.max(0) as usize),
        }
    }
}

fn diag(a: &CsrMatrix<f64>, i: usize) -> f64 {
    (a.row_ptr[i]..a.row_ptr[i + 1])
        .find(|&k| a.col_idx[k] as usize == i)
        .map_or(0.0, |k| a.values[k])
}

/// `AMGBonds`: marks the strong connections (`bonds[k]` for matrix entry
/// `k`) between candidate rows of the same component; a candidate whose
/// largest off-diagonal is below `MG Eliminate Dirichlet Limit` times its
/// diagonal is fixed (`cands = -1`, `fixed`) when Dirichlet rows are
/// eliminated. `cands[i] != 0` marks the candidates.
pub(crate) fn bonds(a: &CsrMatrix<f64>, cands: &mut [i32], fixed: &mut [bool], components: usize, o: &BondOptions) -> Vec<bool> {
    let n = a.nrows;
    let mut bonds = vec![false; a.values.len()];
    let maxconn = (0..n).map(|i| a.row_ptr[i + 1] - a.row_ptr[i]).max().unwrap_or(0);
    let mut measures = vec![0.0; maxconn];
    for ind in 0..n {
        if cands[ind] == 0 {
            continue;
        }
        let (lo, hi) = (a.row_ptr[ind], a.row_ptr[ind + 1]);
        let diagsign = if diag(a, ind) < 0.0 { -1.0 } else { 1.0 };
        let same = |cj: usize| ind % components == cj % components;
        let (mut maxbond, mut minbond, mut minmaxset) = (0.0f64, 0.0f64, false);
        for j in lo..hi {
            let cj = a.col_idx[j] as usize;
            if !same(cj) || cands[cj] == 0 || cj == ind {
                continue;
            }
            if minmaxset {
                maxbond = maxbond.max(a.values[j]);
                minbond = minbond.min(a.values[j]);
            } else {
                maxbond = a.values[j];
                minbond = a.values[j];
                minmaxset = true;
            }
        }
        if !minmaxset {
            maxbond = 0.0;
        } else if o.allow_pos {
            maxbond = maxbond.abs().max(minbond.abs());
        } else {
            if minbond * diagsign < 0.0 {
                maxbond = minbond;
            }
            maxbond = maxbond.abs();
        }
        if maxbond <= o.dir_lim * diag(a, ind).abs() {
            if o.eliminate_dirichlet {
                fixed[ind] = true;
                cands[ind] = -1;
            }
            continue;
        }
        if !minmaxset {
            continue;
        }
        let (mut posnew, mut negnew) = (0usize, 0usize);
        for j in lo..hi {
            let cj = a.col_idx[j] as usize;
            if !same(cj) || cands[cj] == 0 || cj == ind {
                continue;
            }
            let v = a.values[j];
            if diagsign * v <= 0.0 {
                let meas = v.abs() / (o.neg_lim * maxbond);
                measures[j - lo] = -meas;
                if meas > 1.0 {
                    bonds[j] = true;
                    negnew += 1;
                }
            } else if o.allow_pos {
                let meas = v.abs() / (o.pos_lim * maxbond);
                measures[j - lo] = meas;
                if meas > 1.0 {
                    bonds[j] = true;
                    posnew += 1;
                }
            } else {
                measures[j - lo] = 0.0;
            }
        }
        if o.max_conns > 0 {
            while posnew + negnew > o.max_conns {
                let mut measlim = f64::MAX;
                let mut measind = lo;
                for j in lo..hi {
                    if !bonds[j] {
                        continue;
                    }
                    let meas = measures[j - lo].abs();
                    if meas < measlim {
                        measlim = meas;
                        measind = j;
                    }
                }
                if measures[measind - lo] < 0.0 {
                    negnew -= 1;
                } else {
                    posnew -= 1;
                }
                bonds[measind] = false;
            }
        }
        if o.min_conns > 0 {
            while posnew + negnew < o.min_conns {
                let (mut measlim, mut measind) = (0.0, None);
                for j in lo..hi {
                    if bonds[j] {
                        continue;
                    }
                    let cj = a.col_idx[j] as usize;
                    if cands[cj] == 0 || cj == ind || !same(cj) {
                        continue;
                    }
                    let meas = measures[j - lo].abs();
                    if meas > measlim {
                        measlim = meas;
                        measind = Some(j);
                    }
                }
                let Some(j) = measind.filter(|_| measlim >= 1.0e-50) else { break };
                if measures[j - lo] < 0.0 {
                    negnew += 1;
                } else {
                    posnew += 1;
                }
                bonds[j] = true;
            }
        }
    }
    bonds
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The 5-point Laplacian on an `m x m` grid with an anisotropic x
    /// coupling `eps`.
    fn grid(m: usize, eps: f64) -> CsrMatrix<f64> {
        let n = m * m;
        let mut a = CsrMatrix { nrows: n, ncols: n, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
        for i in 0..n {
            let (x, y) = (i % m, i / m);
            let mut row = vec![(i, 2.0 + 2.0 * eps)];
            if y > 0 {
                row.push((i - m, -1.0));
            }
            if x > 0 {
                row.push((i - 1, -eps));
            }
            if x + 1 < m {
                row.push((i + 1, -eps));
            }
            if y + 1 < m {
                row.push((i + m, -1.0));
            }
            row.sort_by_key(|&(c, _)| c);
            for (c, v) in row {
                a.col_idx.push(c as u32);
                a.values.push(v);
            }
            a.row_ptr.push(a.values.len());
        }
        a
    }

    fn opts() -> BondOptions {
        BondOptions::from_params(&ParamMap::new(), true)
    }

    #[test]
    fn strong_connections_follow_the_dominant_coupling() {
        // With eps = 0.01 only the y couplings exceed 0.06 of the largest.
        let a = grid(4, 0.01);
        let mut cands = vec![1; 16];
        let mut fixed = vec![false; 16];
        let b = bonds(&a, &mut cands, &mut fixed, 1, &opts());
        for i in 0..16 {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k] as usize;
                let vertical = j + 4 == i || i + 4 == j;
                assert_eq!(b[k], vertical, "row {i} col {j}");
            }
        }
        assert!(!fixed.iter().any(|&f| f));
    }

    #[test]
    fn a_dirichlet_row_is_fixed_and_left_out() {
        let mut a = grid(3, 1.0);
        // Row 4 becomes a Dirichlet row: only its diagonal.
        for k in a.row_ptr[4]..a.row_ptr[5] {
            if a.col_idx[k] != 4 {
                a.values[k] = 0.0;
            }
        }
        let mut cands = vec![1; 9];
        let mut fixed = vec![false; 9];
        let b = bonds(&a, &mut cands, &mut fixed, 1, &opts());
        assert!(fixed[4] && cands[4] == -1);
        assert!((a.row_ptr[4]..a.row_ptr[5]).all(|k| !b[k]));
    }

    #[test]
    fn maximum_connections_drop_the_weakest() {
        let a = grid(4, 0.5);
        let mut cands = vec![1; 16];
        let mut fixed = vec![false; 16];
        let o = BondOptions { max_conns: 2, ..opts() };
        let b = bonds(&a, &mut cands, &mut fixed, 1, &o);
        for i in 0..16 {
            let n = (a.row_ptr[i]..a.row_ptr[i + 1]).filter(|&k| b[k]).count();
            assert!(n <= 2, "row {i}: {n}");
        }
    }
}
