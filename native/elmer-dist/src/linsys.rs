//! One rank's `Linear System Save` output from a parallel ElmerSolver run
//! (`SaveLinearSystem`, SolveCore.F90:6938, `Parallel` branch).
//!
//! Rank `R` writes `<prefix>_a.dat.R` (`row col value`) and `<prefix>_b.dat.R`
//! (`row value`), with rows and columns given as Elmer's `GlobalDOFs`. Rows
//! follow the local row order, one `_b` line per local row. The matrix is the
//! local assembled one: shared rows carry only this rank's contributions, so
//! the global system is the sum over ranks, which is what [`crate::DistCsr`]
//! takes.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use elmer_core::CsrMatrix;
use elmer_sparse::save::parse_fortran_f64;

use crate::DistError;

/// One rank's local system with the global id of each local row.
#[derive(Clone, Debug)]
pub struct RankSystem {
    pub gids: Vec<u64>,
    pub matrix: CsrMatrix<f64>,
    pub rhs: Vec<f64>,
}

fn file_err(path: &Path, msg: String) -> DistError {
    DistError::File {
        path: path.to_path_buf(),
        msg,
    }
}

fn read(path: &Path) -> Result<String, DistError> {
    std::fs::read_to_string(path).map_err(|e| file_err(path, e.to_string()))
}

/// Read `<dir>/<prefix>_a.dat.<rank>` and `<dir>/<prefix>_b.dat.<rank>`
/// (stiffness only, as saved without `Save Mass`/`Save Damp`).
pub fn read_rank_system(dir: &Path, prefix: &str, rank: usize) -> Result<RankSystem, DistError> {
    let pa: PathBuf = dir.join(format!("{prefix}_a.dat.{rank}"));
    let pb: PathBuf = dir.join(format!("{prefix}_b.dat.{rank}"));

    let mut gids = Vec::new();
    let mut rhs = Vec::new();
    for (ln, line) in read(&pb)?.lines().enumerate() {
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.is_empty() {
            continue;
        }
        let bad = || file_err(&pb, format!("line {}: expected `row value`", ln + 1));
        if t.len() != 2 {
            return Err(bad());
        }
        gids.push(t[0].parse().map_err(|_| bad())?);
        rhs.push(parse_fortran_f64(t[1]).ok_or_else(bad)?);
    }
    let local: HashMap<u64, u32> = gids
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i as u32))
        .collect();
    if local.len() != gids.len() {
        return Err(file_err(&pb, "a row appears twice".into()));
    }

    let n = gids.len();
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::new();
    let mut values = Vec::new();
    let mut cur = 0usize;
    for (ln, line) in read(&pa)?.lines().enumerate() {
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.is_empty() {
            continue;
        }
        let bad = |what: &str| file_err(&pa, format!("line {}: {what}", ln + 1));
        if t.len() != 3 {
            return Err(bad("expected `row col value`"));
        }
        let lookup = |tok: &str| -> Result<usize, DistError> {
            let g: u64 = tok.parse().map_err(|_| bad("bad index"))?;
            local
                .get(&g)
                .map(|&l| l as usize)
                .ok_or_else(|| bad("index not among the rhs rows"))
        };
        let (i, j) = (lookup(t[0])?, lookup(t[1])?);
        if i < cur {
            return Err(bad("rows out of local order"));
        }
        cur = i;
        row_ptr[i + 1] += 1;
        col_idx.push(j as u32);
        values.push(parse_fortran_f64(t[2]).ok_or_else(|| bad("bad value"))?);
    }
    for i in 0..n {
        row_ptr[i + 1] += row_ptr[i];
    }
    let mut matrix = CsrMatrix {
        nrows: n,
        ncols: n,
        row_ptr,
        col_idx,
        values,
        diag: None,
    };
    elmer_sparse::sort_rows(&mut matrix);
    Ok(RankSystem { gids, matrix, rhs })
}
