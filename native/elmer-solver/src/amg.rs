//! WIP, not wired: algebraic multigrid (`MG Method = algebraic`,
//! `AMGSolve`, Multigrid.F90:1484-5019), ported routine by routine.
//!
//! Done: `AMGBonds` (the strong connections, Multigrid.F90:2223),
//! `AMGCoarse` with the default `MG Coarsening Mode = 3` (the C/F
//! splitting, Multigrid.F90:2935) and `ChooseCoarseNodes` up to the coarse
//! numbering (Multigrid.F90:2094; scalar and `MG Determining Component`
//! candidates, no complex, pseudo-geometric or compatible-relaxation
//! variants).
//! `InterpolateF2C` (the projector, Multigrid.F90:3367),
//! `CRS_ProjectMatrixCreate` (the coarse matrix, Multigrid.F90:4843) and
//! the level hierarchy and cycle of `AMGSolve`/`AMGSweep` ([`build`],
//! [`solve`] on `multigrid::Hierarchy` with `Kind::Algebraic`) are ported
//! too, unit-tested but not pristine-checked.
//! Still to port before `Linear System Preconditioning = multigrid` /
//! `Linear System Solver = Multigrid` with `MG Method = algebraic` can use
//! it: the wiring into `solve_linear` and the preconditioner entry
//! (`MultiGridPrec`), checked against amultg2.
#![allow(dead_code)]

use elmer_core::{CsrMatrix, ParamMap, Params};

use crate::error::SolverError;
use crate::multigrid::{CoarseSolve, Hierarchy, Kind, MgOptions};

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

/// `AMGCoarse` with `MG Coarsening Mode = 3` (the default): coarse nodes
/// taken in order of their count of strongly connected candidates, kept in
/// buckets (`GatLims`, `ConInd`, `RevConInd`); each pick turns its
/// undecided strong neighbours fine and raises the counts of their
/// candidate neighbours. `cf` is `> 0` for coarse nodes (in pick order) and
/// `< 0` for fine ones. Index arithmetic follows the Fortran 1-based
/// arrays.
// The 1-based index loops mirror the Fortran's arrays line for line.
#[allow(clippy::needless_range_loop)]
pub(crate) fn coarse(a: &CsrMatrix<f64>, cands: &mut [i32], bonds: &[bool], cf: &mut [i32], follow_bounds: i32) {
    let nods = a.nrows;
    // 1-based views: row r (1..=nods) spans k in rows[r]..rows[r+1].
    let rows = |r: usize| a.row_ptr[r - 1]..a.row_ptr[r];
    let col = |k: usize| a.col_idx[k] as usize + 1;
    let mut con = vec![0i64; nods + 1];
    let (mut cind, mut find) = (0i32, 0i32);
    let cf1 = |cf: &[i32], i: usize| cf[i - 1];
    let cand = |cands: &[i32], i: usize| cands[i - 1];
    for ind in 1..=nods {
        if cf1(cf, ind) == 0 {
            continue;
        }
        if cf1(cf, ind) < 0 {
            find += 1;
            cf[ind - 1] = -find;
            continue;
        }
        cind += 1;
        cf[ind - 1] = cind;
        if cand(cands, ind) == 0 {
            continue;
        }
        for k in rows(ind) {
            if !bonds[k] {
                continue;
            }
            let ci = col(k);
            if cand(cands, ci) == 0 {
                continue;
            }
            if cf1(cf, ci) == 0 {
                find += 1;
                cf[ci - 1] = -find;
            }
        }
    }
    for ind in 1..=nods {
        if cf1(cf, ind) > 0 || cand(cands, ind) <= 0 {
            continue;
        }
        for k in rows(ind) {
            if !bonds[k] {
                continue;
            }
            let ci = col(k);
            if cand(cands, ci) > 0 && cf1(cf, ci) == 0 {
                if cf1(cf, ind) == 0 {
                    con[ci] += 1;
                }
                if cf1(cf, ind) < 0 {
                    con[ci] += 2;
                }
            }
        }
    }
    if follow_bounds != 0 {
        for ind in 1..=nods {
            if cand(cands, ind) <= 0 {
                continue;
            }
            for k in rows(ind) {
                if bonds[k] && cand(cands, col(k)) == -1 {
                    con[ind] += i64::from(follow_bounds);
                }
            }
        }
        if follow_bounds < 0 {
            let m = con[1..].iter().copied().min().unwrap_or(0);
            con[1..].iter_mut().for_each(|c| *c -= m);
        }
    }
    cands.iter_mut().for_each(|c| *c = (*c).max(0));

    let mut maxcon = con[1..].iter().copied().max().unwrap_or(0);
    let gatmax = (2 * maxcon + 1).max(1) as usize;
    let mut gatlims = vec![0usize; gatmax + 1];
    let mut gatnos = vec![0usize; gatmax + 1];
    for ind in 1..=nods {
        if con[ind] != 0 {
            gatnos[con[ind] as usize] += 1;
        }
    }
    let mut anods: usize = gatnos.iter().sum();
    for i in 2..=maxcon.max(1) as usize {
        gatlims[i] = gatlims[i - 1] + gatnos[i - 1];
    }
    let mut conind = vec![0usize; nods + 1];
    let mut revconind = vec![0usize; anods + 1];
    gatnos.iter_mut().for_each(|g| *g = 0);
    for ind in 1..=nods {
        if con[ind] == 0 {
            continue;
        }
        let c = con[ind] as usize;
        gatnos[c] += 1;
        conind[ind] = gatnos[c] + gatlims[c];
    }
    for ind in 1..=nods {
        if con[ind] != 0 {
            revconind[conind[ind]] = ind;
        }
    }
    while maxcon > 0 {
        let ind = revconind[anods];
        cind += 1;
        cf[ind - 1] = cind;
        for j in rows(ind) {
            if !bonds[j] {
                continue;
            }
            let cj = col(j);
            if cf1(cf, cj) != 0 {
                continue;
            }
            find += 1;
            cf[cj - 1] = -find;
            if cand(cands, cj) == 0 {
                continue;
            }
            for i in rows(cj) {
                if !bonds[i] {
                    continue;
                }
                let ci = col(i);
                if cand(cands, ci) > 0 && cf1(cf, ci) == 0 {
                    let points = (con[ci] + 1) as usize;
                    if points >= gatlims.len() {
                        gatlims.resize(points + 1, 0);
                    }
                    let oldorder = conind[ci];
                    if gatlims[points] == 0 {
                        gatlims[points] = anods;
                    }
                    let neworder = gatlims[points];
                    con[ci] = points as i64;
                    gatlims[points] -= 1;
                    if neworder != oldorder {
                        let k = revconind[neworder];
                        conind[ci] = neworder;
                        conind[k] = oldorder;
                        revconind[neworder] = ci;
                        revconind[oldorder] = k;
                    }
                }
            }
        }
        let mut refcon = 0i64;
        while anods > 0 {
            if cf1(cf, revconind[anods]) == 0 {
                break;
            }
            maxcon = con[revconind[anods]];
            refcon = refcon.max(maxcon);
            anods -= 1;
        }
        if anods > 0 {
            maxcon = con[revconind[anods]];
            refcon = refcon.max(maxcon);
            let (lo, hi) = (maxcon as usize + 1, refcon as usize + 1);
            if hi >= gatlims.len() {
                gatlims.resize(hi + 1, 0);
            }
            for g in &mut gatlims[lo..=hi] {
                *g = 0;
            }
        } else {
            maxcon = 0;
        }
    }
    for ind in 1..=nods {
        if cand(cands, ind) > 0 && cf1(cf, ind) == 0 {
            cind += 1;
            cf[ind - 1] = cind;
        }
    }
}

/// `ChooseCoarseNodes` up to the coarse numbering: candidates (every row,
/// or every `components`-th from `MG Determining Component`), strong
/// connections, the C/F split, then `cf[i]` the 1-based coarse index of a
/// coarse row and 0 for a fine one. Returns `cf`, the coarse count and the
/// rows fixed by Dirichlet elimination.
pub(crate) fn choose_coarse_nodes(
    a: &CsrMatrix<f64>,
    p: &ParamMap,
    components: usize,
    finest: bool,
) -> Result<(Vec<i32>, usize, Vec<bool>), SolverError> {
    for k in ["mg complex matrix", "mg pseudo geometric", "mg compatible relax init", "mg compatible relax merit", "mg positive connection eliminate"] {
        if p.logical(k) == Some(true) {
            return Err(SolverError::Unsupported(format!("algebraic multigrid with {k}")));
        }
    }
    let mode = p.integer("mg coarsening mode").unwrap_or(3);
    if mode != 3 {
        return Err(SolverError::Unsupported(format!("MG Coarsening Mode = {mode}")));
    }
    let nods = a.nrows;
    let first = if components > 1 {
        p.integer("mg determining component").map_or(1, |c| c.clamp(1, components as i32) as usize)
    } else {
        1
    };
    let mut cands = vec![0i32; nods];
    for i in (first - 1..nods).step_by(components.max(1)) {
        cands[i] = 1;
    }
    let mut fixed = vec![false; nods];
    let b = bonds(a, &mut cands, &mut fixed, components, &BondOptions::from_params(p, finest));
    let mut cf = vec![0i32; nods];
    coarse(a, &mut cands, &b, &mut cf, p.integer("mg boundary priority").unwrap_or(1));
    let mut cnods = 0;
    for c in cf.iter_mut() {
        if *c < 0 {
            *c = 0;
        } else if *c > 0 {
            cnods += 1;
            *c = cnods;
        }
    }
    Ok((cf, cnods as usize, fixed))
}

/// The `MG*` keywords `InterpolateF2C` reads.
#[derive(Clone, Debug)]
pub(crate) struct InterpOptions {
    /// `MG Projection Limit` (0.6).
    pub proj_lim: f64,
    /// `MG Projection Lumping`.
    pub lumping: bool,
    /// `MG Direct Interpolate` and `MG Direct Interpolate Limit`.
    pub direct: bool,
    pub direct_limit: usize,
    /// `MG Negative Connection Favor Factor` (1).
    pub favor_neg: f64,
    /// `MG Determining Component` (1-based).
    pub component1: usize,
}

impl InterpOptions {
    pub(crate) fn from_params(p: &ParamMap, components: usize) -> Self {
        let direct = p.logical("mg direct interpolate") == Some(true);
        Self {
            proj_lim: p.real("mg projection limit").unwrap_or(0.6),
            lumping: p.logical("mg projection lumping") == Some(true),
            direct,
            direct_limit: p
                .integer("mg direct interpolate limit")
                .map_or(if direct { 0 } else { usize::MAX }, |k| k.max(0) as usize),
            favor_neg: p.real("mg negative connection favor factor").unwrap_or(1.0),
            component1: if components > 1 {
                p.integer("mg determining component").map_or(1, |c| c.clamp(1, components as i32) as usize)
            } else {
                1
            },
        }
    }
}

/// `InterpolateF2C`: the projector (fine nodes x coarse nodes, the
/// determining component's rows) from the C/F split `cf` (1-based coarse
/// indices, 0 for fine) and the rows `fixed` by Dirichlet elimination.
/// A coarse row interpolates itself; a fine one from its coarse
/// neighbours, through its fine neighbours too (unless directly), with the
/// dominant sign's coefficients kept beyond `MG Projection Limit` of the
/// largest and the other sign's while they still count.
#[allow(clippy::needless_range_loop)]
pub(crate) fn interpolate_f2c(a: &CsrMatrix<f64>, cf: &[i32], fixed: &[bool], components: usize, o: &InterpOptions) -> CsrMatrix<f64> {
    let fdofs = a.nrows;
    let comps = components.max(1);
    let same = |x: usize, y: usize| (x + 1) % comps == (y + 1) % comps;
    let val_at = |k: usize| a.values[k];
    let diag_k = |i: usize| (a.row_ptr[i]..a.row_ptr[i + 1]).find(|&k| a.col_idx[k] as usize == i);
    let ncoarse = cf.iter().filter(|&&c| c > 0).count();
    let nrows = fdofs / comps;
    let mut p = CsrMatrix { nrows, ncols: ncoarse, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
    let mut coeffs_inds = vec![0usize; fdofs];
    for ind in (o.component1 - 1..fdofs).step_by(comps) {
        if cf[ind] > 0 {
            p.col_idx.push(cf[ind] as u32 - 1);
            p.values.push(1.0);
            p.row_ptr.push(p.values.len());
            continue;
        }
        let mut inds: Vec<usize> = Vec::new();
        let mut coeffs: Vec<f64> = Vec::new();
        let mut diag = 0usize;
        let mut ncoarse_nb = 0usize;
        if o.direct {
            for k in a.row_ptr[ind]..a.row_ptr[ind + 1] {
                let ci = a.col_idx[k] as usize;
                if !same(ci, ind) || fixed[ci] {
                    continue;
                }
                let v = val_at(k);
                if v.abs() < 1.0e-50 {
                    continue;
                }
                inds.push(ci);
                coeffs.push(v);
                if ci == ind {
                    diag = inds.len();
                } else if cf[ci] > 0 {
                    ncoarse_nb += 1;
                }
            }
        }
        if ncoarse_nb < o.direct_limit {
            inds.clear();
            coeffs.clear();
            diag = 0;
            for k in a.row_ptr[ind]..a.row_ptr[ind + 1] {
                let ci = a.col_idx[k] as usize;
                if !same(ci, ind) || fixed[ci] {
                    continue;
                }
                if cf[ci] > 0 || ci == ind {
                    let v = val_at(k);
                    if v.abs() < 1.0e-50 {
                        continue;
                    }
                    inds.push(ci);
                    coeffs.push(v);
                    coeffs_inds[ci] = inds.len();
                    if ci == ind {
                        diag = inds.len();
                    }
                }
            }
            for k in a.row_ptr[ind]..a.row_ptr[ind + 1] {
                let ci = a.col_idx[k] as usize;
                if !same(ci, ind) || cf[ci] > 0 || ci == ind || fixed[ci] {
                    continue;
                }
                let Some(dk) = diag_k(ci) else { continue };
                for j in a.row_ptr[ci]..a.row_ptr[ci + 1] {
                    let cj = a.col_idx[j] as usize;
                    if cj == ci || !same(cj, ind) || fixed[cj] {
                        continue;
                    }
                    let v = val_at(k) * val_at(j) / val_at(dk);
                    if v.abs() < 1.0e-50 {
                        continue;
                    }
                    let mut slot = coeffs_inds[cj];
                    if slot == 0 {
                        inds.push(cj);
                        coeffs.push(0.0);
                        slot = inds.len();
                        coeffs_inds[cj] = slot;
                    }
                    coeffs[slot - 1] -= v;
                }
            }
            for &i in &inds {
                coeffs_inds[i] = 0;
            }
        }
        let no = inds.len();
        if no <= 1 || diag == 0 {
            p.row_ptr.push(p.values.len());
            continue;
        }
        let d = diag - 1;
        if o.lumping {
            let s: f64 = coeffs.iter().sum();
            coeffs[d] -= s;
        }
        if coeffs[d] < 0.0 {
            coeffs.iter_mut().for_each(|c| *c = -*c);
        }
        let (mut possum, mut negsum, mut posmax, mut negmax) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut posi, mut negi) = (0usize, 0usize);
        let mut diagsum = 0.0;
        for i in 0..no {
            let v = coeffs[i];
            if i == d {
                diagsum = v;
            } else if v > 0.0 {
                possum += v;
                if cf[inds[i]] > 0 {
                    posmax = posmax.max(v);
                    posi += 1;
                }
            } else {
                negsum += v;
                if cf[inds[i]] > 0 {
                    negmax = negmax.min(v);
                    negi += 1;
                }
            }
        }
        if posi == 0 && negi == 0 {
            p.row_ptr.push(p.values.len());
            continue;
        }
        let mut pos: Vec<(usize, f64)> = Vec::new();
        let mut neg: Vec<(usize, f64)> = Vec::new();
        let (mut poscsum, mut negcsum) = (0.0, 0.0);
        if posi == 0 || possum.abs() <= negsum.abs() * o.favor_neg {
            if negi == 0 {
                p.row_ptr.push(p.values.len());
                continue;
            }
            for i in 0..no {
                let (v, ci) = (coeffs[i], inds[i]);
                if cf[ci] == 0 {
                    continue;
                }
                if v < o.proj_lim * negmax {
                    neg.push((ci, v));
                    negcsum += v;
                } else if v > o.proj_lim * posmax {
                    pos.push((ci, v));
                }
            }
            negi = neg.len();
            posi = 0;
            let refbond = -o.proj_lim * negsum * negmax / (o.favor_neg * negcsum);
            if possum > refbond {
                bubble(&mut pos, |a, b| a < b);
                poscsum = 0.0;
                for i in 0..pos.len() {
                    if i == 0 {
                        if possum < refbond {
                            break;
                        }
                    } else {
                        let posbond = possum * pos[i].1 / (poscsum + pos[i].1);
                        if posbond < refbond && pos[i].1 < 0.99 * pos[i - 1].1 {
                            break;
                        }
                    }
                    posi = i + 1;
                    poscsum += pos[i].1;
                }
            }
        } else {
            for i in 0..no {
                let (v, ci) = (coeffs[i], inds[i]);
                if cf[ci] == 0 {
                    continue;
                }
                if v > o.proj_lim * posmax {
                    pos.push((ci, v));
                    poscsum += v;
                } else if v < o.proj_lim * negmax {
                    neg.push((ci, v));
                }
            }
            posi = pos.len();
            negi = 0;
            let refbond = o.proj_lim * possum * posmax / (poscsum * o.favor_neg);
            if -negsum > refbond {
                bubble(&mut neg, |a, b| a > b);
                negcsum = 0.0;
                for i in 0..neg.len() {
                    if i == 0 {
                        if -negsum < refbond {
                            break;
                        }
                    } else {
                        let negbond = negsum * neg[i].1 / (negcsum + neg[i].1);
                        if -negbond < refbond && neg[i].1 > 0.99 * neg[i - 1].1 {
                            break;
                        }
                    }
                    negi = i + 1;
                    negcsum += neg[i].1;
                }
            }
        }
        if posi == 0 {
            diagsum += possum / o.favor_neg;
        }
        if negi == 0 {
            diagsum = (diagsum + negsum) / o.favor_neg;
        }
        for &(ci, v) in &neg[..negi] {
            p.col_idx.push(cf[ci] as u32 - 1);
            p.values.push(-negsum * v / (negcsum * diagsum));
        }
        for &(ci, v) in &pos[..posi] {
            p.col_idx.push(cf[ci] as u32 - 1);
            p.values.push(-possum * v / (poscsum * diagsum * o.favor_neg));
        }
        p.row_ptr.push(p.values.len());
    }
    p
}

/// The Fortran's bubble passes: entries `i` and `i + 1` swap when
/// `swap(v_i, v_i+1)`.
fn bubble(v: &mut [(usize, f64)], swap: impl Fn(f64, f64) -> bool) {
    let n = v.len();
    for _ in 1..n {
        for i in 0..n - 1 {
            if swap(v[i].1, v[i + 1].1) {
                v.swap(i, i + 1);
            }
        }
    }
}

/// `CRS_ProjectMatrixCreate(A, P, R, B, DOFs)` with `P = R^T`: the coarse
/// matrix `B = P A R` from the node projector `r` (fine nodes x coarse
/// nodes, `ProjT`), `dofs` components per node; with several components
/// the entries of `A` below `MG Matrix Create Epsilon` are skipped. The
/// Fortran stores a row's columns in the order first met; they are sorted
/// here (the sums are unchanged).
pub(crate) fn project_matrix_create(a: &CsrMatrix<f64>, r: &CsrMatrix<f64>, dofs: usize, epsilon: f64) -> CsrMatrix<f64> {
    // P = R^T, rows in increasing fine order (CRS_Transpose).
    let nc = r.ncols;
    let mut prow: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nc];
    for i in 0..r.nrows {
        for k in r.row_ptr[i]..r.row_ptr[i + 1] {
            prow[r.col_idx[k] as usize].push((i, r.values[k]));
        }
    }
    let nb = nc * dofs;
    let mut b = CsrMatrix { nrows: nb, ncols: nb, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
    let mut slot = vec![usize::MAX; nb];
    for prow_i in &prow {
        for m in 0..dofs {
            let mut cols: Vec<usize> = Vec::new();
            let mut vals: Vec<f64> = Vec::new();
            for &(pj, pv) in prow_i {
                let cj = dofs * pj + m;
                for k in a.row_ptr[cj]..a.row_ptr[cj + 1] {
                    if dofs > 1 && a.values[k].abs() < epsilon {
                        continue;
                    }
                    let ck = a.col_idx[k] as usize;
                    let (n, node) = (ck % dofs, ck / dofs);
                    for l in r.row_ptr[node]..r.row_ptr[node + 1] {
                        let cl = dofs * r.col_idx[l] as usize + n;
                        let v = pv * a.values[k] * r.values[l];
                        if slot[cl] == usize::MAX {
                            slot[cl] = cols.len();
                            cols.push(cl);
                            vals.push(v);
                        } else {
                            vals[slot[cl]] += v;
                        }
                    }
                }
            }
            let mut row: Vec<(usize, f64)> = cols.iter().copied().zip(vals).collect();
            for &c in &cols {
                slot[c] = usize::MAX;
            }
            row.sort_by_key(|&(c, _)| c);
            for (c, v) in row {
                b.col_idx.push(c as u32);
                b.values.push(v);
            }
            b.row_ptr.push(b.values.len());
        }
    }
    b
}

/// The level hierarchy of `AMGSolve` (coarsest first): from the finest
/// matrix down while the level is above 1 and the matrix has at least `MG
/// Lowest Linear Solver Limit` (20) rows; Dirichlet rows are eliminated on
/// the finest level only. The interpolation operators are expanded to the
/// dofs (`CRS_ProjectVector` applies the node projector per component).
pub(crate) fn build(a: &CsrMatrix<f64>, p: &ParamMap, dofs: usize, levels: usize, ilu: Option<usize>) -> Result<Hierarchy, SolverError> {
    let direct_limit = p.integer("mg lowest linear solver limit").map_or(20, |k| k.max(0) as usize);
    let epsilon = p.real("mg matrix create epsilon").unwrap_or(f64::MIN_POSITIVE.sqrt());
    let mut mats: Vec<(CsrMatrix<f64>, Option<CsrMatrix<f64>>)> = Vec::new();
    let mut cur = a.clone();
    let mut level = levels;
    loop {
        if level <= 1 || cur.nrows < direct_limit {
            mats.push((cur, None));
            break;
        }
        let (cf, _, fixed) = choose_coarse_nodes(&cur, p, dofs, level == levels)?;
        let proj = interpolate_f2c(&cur, &cf, &fixed, dofs, &InterpOptions::from_params(p, dofs));
        let next = project_matrix_create(&cur, &proj, dofs, epsilon);
        mats.push((cur, Some(expand(&proj, dofs))));
        cur = next;
        level -= 1;
    }
    mats.reverse();
    let mut h = Hierarchy::new(mats, ilu)?;
    h.kind = Kind::Algebraic;
    h.dofs = dofs;
    Ok(h)
}

/// The node projector applied per component: row `dofs i + c`, columns
/// `dofs j + c`.
fn expand(p: &CsrMatrix<f64>, dofs: usize) -> CsrMatrix<f64> {
    if dofs == 1 {
        return p.clone();
    }
    let mut e = CsrMatrix { nrows: p.nrows * dofs, ncols: p.ncols * dofs, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
    for i in 0..p.nrows {
        for c in 0..dofs {
            for k in p.row_ptr[i]..p.row_ptr[i + 1] {
                e.col_idx.push(dofs as u32 * p.col_idx[k] + c as u32);
                e.values.push(p.values[k]);
            }
            e.row_ptr.push(e.values.len());
        }
    }
    e
}

/// `AMGSolve` at the finest level of `h`.
pub(crate) fn solve(
    h: &Hierarchy,
    x: &mut [f64],
    b: &mut [f64],
    opts: &MgOptions,
    coarse: &CoarseSolve<'_>,
    log: &dyn elmer_core::Logger,
) -> Result<(), SolverError> {
    h.solve(h.levels.len(), x, b, opts, coarse, log)
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

    #[test]
    fn coarsening_splits_the_isotropic_grid() {
        let a = grid(5, 1.0);
        let (cf, nc, _) = choose_coarse_nodes(&a, &ParamMap::new(), 1, true).unwrap();
        // Every fine node has a strongly connected coarse neighbour and the
        // coarse nodes are numbered 1..=nc.
        assert!(nc > 0 && nc < 25, "{nc}");
        let mut seen: Vec<i32> = cf.iter().copied().filter(|&c| c > 0).collect();
        seen.sort_unstable();
        assert_eq!(seen, (1..=nc as i32).collect::<Vec<_>>());
        for i in 0..25 {
            if cf[i] == 0 {
                let has = (a.row_ptr[i]..a.row_ptr[i + 1]).any(|k| cf[a.col_idx[k] as usize] > 0);
                assert!(has, "fine node {i} without a coarse neighbour");
            }
        }
    }

    #[test]
    fn anisotropic_grid_coarsens_along_the_strong_direction() {
        // Strong y couplings only: each column is coarsened on its own, so
        // no two vertically adjacent nodes are both coarse.
        let a = grid(6, 0.01);
        let (cf, _, _) = choose_coarse_nodes(&a, &ParamMap::new(), 1, true).unwrap();
        for i in 0..30 {
            assert!(!(cf[i] > 0 && cf[i + 6] > 0), "nodes {i} and {} both coarse", i + 6);
        }
    }

    #[test]
    fn interpolation_reproduces_constants() {
        // Rows summing to zero (no Dirichlet rows): a constant is
        // interpolated exactly, so each fine row of the projector sums to
        // one; a coarse row is the identity.
        let m = 6;
        let n = m * m;
        let mut a = grid(m, 1.0);
        for i in 0..n {
            let s: f64 = (a.row_ptr[i]..a.row_ptr[i + 1])
                .filter(|&k| a.col_idx[k] as usize != i)
                .map(|k| a.values[k])
                .sum();
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                if a.col_idx[k] as usize == i {
                    a.values[k] = -s;
                }
            }
        }
        let (cf, nc, fixed) = choose_coarse_nodes(&a, &ParamMap::new(), 1, true).unwrap();
        let p = interpolate_f2c(&a, &cf, &fixed, 1, &InterpOptions::from_params(&ParamMap::new(), 1));
        assert_eq!((p.nrows, p.ncols), (n, nc));
        for i in 0..n {
            let row = &p.values[p.row_ptr[i]..p.row_ptr[i + 1]];
            if cf[i] > 0 {
                assert_eq!(row, &[1.0]);
            } else {
                let s: f64 = row.iter().sum();
                assert!((s - 1.0).abs() < 1e-12, "row {i}: {row:?}");
            }
        }
    }

    fn dense_solve(a: &CsrMatrix<f64>, b: &[f64], x: &mut [f64]) -> Result<(), SolverError> {
        let n = a.nrows;
        let mut m = vec![0.0; n * n];
        for i in 0..n {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                m[i * n + a.col_idx[k] as usize] = a.values[k];
            }
        }
        x.copy_from_slice(b);
        for c in 0..n {
            for r in c + 1..n {
                let f = m[r * n + c] / m[c * n + c];
                for k in c..n {
                    m[r * n + k] -= f * m[c * n + k];
                }
                x[r] -= f * x[c];
            }
        }
        for c in (0..n).rev() {
            let s: f64 = (c + 1..n).map(|k| m[c * n + k] * x[k]).sum();
            x[c] = (x[c] - s) / m[c * n + c];
        }
        Ok(())
    }

    #[test]
    fn amg_cycles_solve_the_grid_laplacian() {
        let m = 12;
        let n = m * m;
        let a = grid(m, 1.0);
        let h = build(&a, &ParamMap::new(), 1, 10, None).unwrap();
        assert!(h.levels.len() >= 3, "{} levels", h.levels.len());
        let mut b = vec![1.0; n];
        let mut x = vec![0.0; n];
        let opts = MgOptions {
            max_iter: 100,
            tolerance: 1.0e-10,
            level_max_iter: 1,
            level_tolerance: f64::MAX,
            sweeps: vec![1],
            pre: vec![2],
            post: vec![2],
            smoother: "sgs".into(),
        };
        solve(&h, &mut x, &mut b, &opts, &dense_solve, &elmer_core::NullLogger).unwrap();
        let mut exact = vec![0.0; n];
        dense_solve(&a, &b, &mut exact).unwrap();
        for i in 0..n {
            assert!((x[i] - exact[i]).abs() < 1e-8 * exact[i].abs().max(1.0), "{i}: {} vs {}", x[i], exact[i]);
        }
    }

    #[test]
    fn coarse_matrix_is_the_galerkin_product() {
        let a = grid(6, 1.0);
        let (cf, _, fixed) = choose_coarse_nodes(&a, &ParamMap::new(), 1, true).unwrap();
        let r = interpolate_f2c(&a, &cf, &fixed, 1, &InterpOptions::from_params(&ParamMap::new(), 1));
        let b = project_matrix_create(&a, &r, 1, 0.0);
        let g = crate::multigrid::galerkin(&a, &r);
        assert_eq!((b.row_ptr.clone(), b.col_idx.clone()), (g.row_ptr.clone(), g.col_idx.clone()));
        for (x, y) in b.values.iter().zip(&g.values) {
            assert!((x - y).abs() < 1e-12 * x.abs().max(1.0));
        }
    }
}
