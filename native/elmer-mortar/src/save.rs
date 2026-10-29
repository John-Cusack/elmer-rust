//! `Save Projector` and `Save Constraint Matrix` output (`SaveProjector`,
//! MortarUtils.F90), serial case.
//!
//! `<prefix>.dat` has one line `row col value` per stored entry,
//! `<prefix>_rsum.dat` one line `row i nnz diag rowsum` per row and, for a
//! matrix with a right-hand side, `<prefix>_rhs.dat` one line `row i rhs`.
//! `row` is the row's `InvPerm` entry and `i` the 1-based matrix row, all
//! written list-directed by gfortran.

use std::fmt::Write as _;
use std::io;
use std::path::{Path, PathBuf};

use elmer_core::NO_DOF;
use elmer_sparse::save::{list_f64, list_int};

use crate::constraint::ConstraintMatrix;
use crate::projector::Projector;

/// A CRS matrix as `SaveProjector` sees it, indices 1-based.
struct Crs<'a> {
    rows: &'a [usize],
    /// `InvPerm(i)`, or `None` for a row `SaveProjector` skips
    /// (`InvPerm(i) == 0`).
    inv_perm: Vec<Option<i64>>,
    cols: Vec<i64>,
    values: &'a [f64],
    rhs: Option<&'a [f64]>,
}

fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

/// `SaveProjector(Projector, SaveRowSum=.TRUE., Prefix)`.
fn save(m: &Crs<'_>, prefix: &Path) -> io::Result<()> {
    let mut dat = String::new();
    let mut rsum = String::new();
    let mut rhs = String::new();
    for i in 0..m.rows.len() - 1 {
        let Some(ii) = m.inv_perm[i] else { continue };
        let r = m.rows[i]..m.rows[i + 1];
        let mut dia = 0.0;
        let mut sum = 0.0;
        for k in r.clone() {
            let (jj, v) = (m.cols[k], m.values[k]);
            writeln!(dat, "{}{}{}", list_int(ii), list_int(jj), list_f64(v)).unwrap();
            if jj == ii {
                dia = v;
            }
            sum += v;
        }
        writeln!(
            rsum,
            "{}{}{}{}{}",
            list_int(ii),
            list_int(i as i64 + 1),
            list_int(r.len() as i64),
            list_f64(dia),
            list_f64(sum)
        )
        .unwrap();
        if let Some(b) = m.rhs {
            writeln!(
                rhs,
                "{}{}{}",
                list_int(ii),
                list_int(i as i64 + 1),
                list_f64(b[i])
            )
            .unwrap();
        }
    }
    std::fs::write(with_suffix(prefix, ".dat"), dat)?;
    std::fs::write(with_suffix(prefix, "_rsum.dat"), rsum)?;
    if m.rhs.is_some() {
        std::fs::write(with_suffix(prefix, "_rhs.dat"), rhs)?;
    }
    Ok(())
}

/// `SaveProjector(Projector, SaveRowSum=.TRUE., Prefix)`: writes
/// `<prefix>.dat` and `<prefix>_rsum.dat` in parent node numbers. Rows
/// without a parent node are skipped, as the Fortran skips `InvPerm == 0`.
pub fn save_projector(p: &Projector, prefix: &Path) -> io::Result<()> {
    save(
        &Crs {
            rows: &p.rows,
            inv_perm: p
                .inv_perm
                .iter()
                .map(|&n| (n != NO_DOF).then_some(n as i64 + 1))
                .collect(),
            cols: p.cols.iter().map(|&c| c as i64 + 1).collect(),
            values: &p.values,
            rhs: None,
        },
        prefix,
    )
}

/// `Save Constraint Matrix` (`SaveProjector(A % ConstraintMatrix,
/// .TRUE., 'cm')`): writes `<prefix>.dat`, `<prefix>_rsum.dat` and
/// `<prefix>_rhs.dat` in matrix row numbers.
pub fn save_constraint_matrix(c: &ConstraintMatrix, prefix: &Path) -> io::Result<()> {
    save(
        &Crs {
            rows: &c.rows,
            inv_perm: c.inv_perm.iter().map(|&n| Some(n as i64 + 1)).collect(),
            cols: c.cols.iter().map(|&c| c as i64 + 1).collect(),
            values: &c.values,
            rhs: Some(&c.rhs),
        },
        prefix,
    )
}
