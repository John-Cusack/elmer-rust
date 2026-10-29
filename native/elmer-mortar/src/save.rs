//! `Save Projector` output (`SaveProjector`, MortarUtils.F90), serial case.
//!
//! `<prefix>.dat` has one line `row col value` per stored entry and
//! `<prefix>_rsum.dat` one line `row i nnz diag rowsum` per row, where `row`
//! and `col` are 1-based parent mesh nodes and `i` the 1-based projector row,
//! all written list-directed by gfortran.

use std::fmt::Write as _;
use std::io;
use std::path::Path;

use elmer_core::NO_DOF;
use elmer_sparse::save::{list_f64, list_int};

use crate::projector::Projector;

/// `SaveProjector(Projector, SaveRowSum=.TRUE., Prefix)`: writes
/// `<prefix>.dat` and `<prefix>_rsum.dat`. Rows without a parent node are
/// skipped, as the Fortran skips `InvPerm == 0`.
pub fn save_projector(p: &Projector, prefix: &Path) -> io::Result<()> {
    let mut dat = String::new();
    let mut rsum = String::new();
    for i in 0..p.n_rows() {
        let node = p.inv_perm[i];
        if node == NO_DOF {
            continue;
        }
        let ii = node as i64 + 1;
        let mut dia = 0.0;
        let mut sum = 0.0;
        for (c, v) in p.row(i) {
            let jj = c as i64 + 1;
            writeln!(dat, "{}{}{}", list_int(ii), list_int(jj), list_f64(v)).unwrap();
            if jj == ii {
                dia = v;
            }
            sum += v;
        }
        let nnz = (p.rows[i + 1] - p.rows[i]) as i64;
        writeln!(
            rsum,
            "{}{}{}{}{}",
            list_int(ii),
            list_int(i as i64 + 1),
            list_int(nnz),
            list_f64(dia),
            list_f64(sum)
        )
        .unwrap();
    }
    let with = |suffix: &str| {
        let mut s = prefix.as_os_str().to_owned();
        s.push(suffix);
        std::path::PathBuf::from(s)
    };
    std::fs::write(with(".dat"), dat)?;
    std::fs::write(with("_rsum.dat"), rsum)
}
