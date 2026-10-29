//! The `Linear System Save` dump format (`SaveLinearSystem`, `PrintMatrix`,
//! `PrintRHS`; SolveCore.F90:890-1046 and 6938), serial case.
//!
//! Files, all with 1-based indices and gfortran list-directed numbers:
//!
//! | file | lines |
//! |---|---|
//! | `<prefix>_a.dat` | `row col v...`, one line per stored entry; `v...` are stiffness, damping, mass in that order, only those saved |
//! | `<prefix>_b.dat` | `row rhs` (or `row rhs rowsum abssum` with `Save Sum`) |
//! | `<prefix>_perm.dat` | `i Perm(i)` (0 = not in the system) |
//! | `<prefix>_sol.dat` | `i value`, `ES15.6` |
//! | `<prefix>_sizes.dat` | rows, nnz, perm size, min col, max col, optional offset |

use elmer_core::{CsrMatrix, NO_DOF};
use std::fmt::Write as _;
use std::io;
use std::path::Path;

use crate::ops::find_diag;

/// Everything one `Linear System Save` writes.
#[derive(Clone, Debug, Default)]
pub struct SavedSystem {
    pub stiff: Option<CsrMatrix<f64>>,
    /// Damping coefficients sharing the structure of `stiff`.
    pub damp: Option<Vec<f64>>,
    /// Mass coefficients sharing the structure of `stiff`.
    pub mass: Option<Vec<f64>>,
    pub rhs: Option<Vec<f64>>,
    /// Node -> dof slot, [`NO_DOF`] where absent.
    pub perm: Option<Vec<u32>>,
    pub solution: Option<Vec<f64>>,
}

/// Which coefficient columns `_a.dat` carries, in file order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Columns {
    pub stiff: bool,
    pub damp: bool,
    pub mass: bool,
}

impl Default for Columns {
    fn default() -> Self {
        Columns {
            stiff: true,
            damp: false,
            mass: false,
        }
    }
}

impl Columns {
    fn count(self) -> usize {
        self.stiff as usize + self.damp as usize + self.mass as usize
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct SaveOptions {
    /// Keep only the columns that exist; stiffness is on unless `no_stiff`.
    pub no_stiff: bool,
    pub save_damp: bool,
    pub save_mass: bool,
    /// `Linear System Save Skip Zeros`: drop entries whose saved values sum to ~0.
    pub skip_zeros: bool,
    /// `Linear System Save Sum`: add row sum and absolute row sum to `_b.dat`.
    pub save_sum: bool,
    /// The sixth line of `_sizes.dat`.
    pub offset_ind: Option<i64>,
}

/// gfortran list-directed `REAL(8)` item, including its leading separator blank
/// and (for the F branch) trailing exponent-width blanks.
pub fn list_f64(x: f64) -> String {
    if x.is_nan() {
        return format!(" {:>25}", "NaN");
    }
    if x.is_infinite() {
        return format!(" {:>25}", if x < 0.0 { "-Infinity" } else { "Infinity" });
    }
    let neg = x.is_sign_negative();
    let mag = x.abs();
    let body = if mag == 0.0 {
        Some("0.0000000000000000".to_string())
    } else {
        let sci = format!("{:.16e}", mag);
        let e: i32 = sci[sci.find('e').unwrap() + 1..].parse().unwrap();
        if (-1..17).contains(&e) {
            let dec = if e >= 0 { (16 - e) as usize } else { 17 };
            let mut s = format!("{:.*}", dec, mag);
            if dec == 0 {
                s.push('.');
            }
            Some(s)
        } else {
            None
        }
    };
    match body {
        Some(b) => {
            let v = if neg { format!("-{b}") } else { b };
            format!(" {v:>20}     ")
        }
        None => {
            let sci = format!("{:.16e}", mag);
            let (m, e) = sci.split_once('e').unwrap();
            let e: i32 = e.parse().unwrap();
            let v = format!(
                "{}{m}E{}{:03}",
                if neg { "-" } else { "" },
                if e < 0 { '-' } else { '+' },
                e.abs()
            );
            format!(" {v:>25}")
        }
    }
}

/// gfortran list-directed `INTEGER(4)` item.
pub fn list_int(i: i64) -> String {
    format!(" {i:>11}")
}

/// `ES15.6`, as `_sol.dat` uses.
fn es15_6(x: f64) -> String {
    let sci = format!("{:.6e}", x);
    let (m, e) = sci.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    let s = format!("{m}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs());
    format!("{s:>15}")
}

fn path(dir: &Path, prefix: &str, suffix: &str) -> std::path::PathBuf {
    dir.join(format!("{prefix}{suffix}"))
}

/// Write the files of a serial `Linear System Save` into `dir`.
pub fn write_linear_system_save(
    dir: &Path,
    prefix: &str,
    sys: &SavedSystem,
    opts: &SaveOptions,
) -> io::Result<()> {
    let a = sys
        .stiff
        .as_ref()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "no matrix to save"))?;
    let do_stiff = !opts.no_stiff;
    let do_damp = opts.save_damp && sys.damp.is_some();
    let do_mass = opts.save_mass && sys.mass.is_some();
    let mut out = String::new();
    for i in 0..a.nrows {
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            let mut vals = [0.0f64; 3];
            let mut cnt = 0;
            if do_stiff {
                vals[cnt] = a.values[j];
                cnt += 1;
            }
            if do_damp {
                vals[cnt] = sys.damp.as_ref().unwrap()[j];
                cnt += 1;
            }
            if do_mass {
                vals[cnt] = sys.mass.as_ref().unwrap()[j];
                cnt += 1;
            }
            if opts.skip_zeros {
                // SUM(ABS(Vals(1:IndMax))) < EPSILON
                let mut s = 0.0;
                for v in &vals[..cnt] {
                    s += v.abs();
                }
                if s < f64::EPSILON {
                    continue;
                }
            }
            let _ = write!(out, "{} {} ", i + 1, a.col_idx[j] + 1);
            if cnt > 0 {
                for v in &vals[..cnt] {
                    out.push_str(&list_f64(*v));
                }
            } else {
                out.push(' ');
            }
            out.push('\n');
        }
    }
    std::fs::write(path(dir, prefix, "_a.dat"), out)?;

    if sys.rhs.is_some() || opts.save_sum {
        let mut out = String::new();
        for i in 0..a.nrows {
            let _ = write!(out, "{} ", i + 1);
            let val = sys.rhs.as_ref().map(|r| r[i]);
            let sums = opts.save_sum.then(|| {
                let (mut rs, mut asum) = (0.0, 0.0);
                for j in a.row_ptr[i]..a.row_ptr[i + 1] {
                    rs += a.values[j];
                    asum += a.values[j].abs();
                }
                (rs, asum)
            });
            match (val, sums) {
                (Some(v), Some((r, s))) => {
                    out.push_str(&list_f64(v));
                    out.push_str(&list_f64(r));
                    out.push_str(&list_f64(s));
                }
                (None, Some((r, s))) => {
                    out.push_str(&list_f64(r));
                    out.push_str(&list_f64(s));
                }
                (Some(v), None) => {
                    if v.abs() <= f64::MIN_POSITIVE {
                        out.push_str("0.0");
                    } else {
                        out.push_str(&list_f64(v));
                    }
                }
                (None, None) => unreachable!(),
            }
            out.push('\n');
        }
        std::fs::write(path(dir, prefix, "_b.dat"), out)?;
    }

    if let Some(perm) = &sys.perm {
        let mut out = String::new();
        for (i, &p) in perm.iter().enumerate() {
            let slot = if p == NO_DOF { 0 } else { p as u64 + 1 };
            let _ = writeln!(out, "{} {}", i + 1, slot);
        }
        std::fs::write(path(dir, prefix, "_perm.dat"), out)?;
    }
    if let Some(sol) = &sys.solution {
        let mut out = String::new();
        for (i, &v) in sol.iter().enumerate() {
            let _ = writeln!(out, "{}{}", i + 1, es15_6(v));
        }
        std::fs::write(path(dir, prefix, "_sol.dat"), out)?;
    }

    let mut out = String::new();
    let min_col = a.col_idx.iter().min().map_or(0, |&c| c as i64 + 1);
    let max_col = a.col_idx.iter().max().map_or(0, |&c| c as i64 + 1);
    for v in [
        a.nrows as i64,
        a.values.len() as i64,
        sys.perm.as_ref().map_or(0, |p| p.len() as i64),
        min_col,
        max_col,
    ] {
        out.push_str(&list_int(v));
        out.push('\n');
    }
    if let Some(o) = opts.offset_ind {
        out.push_str(&list_int(o));
        out.push('\n');
    }
    std::fs::write(path(dir, prefix, "_sizes.dat"), out)
}

/// Parse a Fortran-written real: `D` exponents and a missing `E` before a
/// signed three-digit exponent are accepted.
pub fn parse_fortran_f64(tok: &str) -> Option<f64> {
    if let Ok(v) = tok.parse::<f64>() {
        return Some(v);
    }
    let t = tok.replace(['D', 'd'], "E");
    if let Ok(v) = t.parse::<f64>() {
        return Some(v);
    }
    // 1.234-100 style: exponent sign directly after the mantissa
    let pos = t.rfind(['+', '-'])?;
    if pos == 0 {
        return None;
    }
    format!("{}E{}", &t[..pos], &t[pos..]).parse().ok()
}

fn bad(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// Read a serial `Linear System Save` back. `cols` says which coefficient
/// columns `_a.dat` holds. `_b`, `_perm`, `_sol` are read when present
/// (`_b.dat` only in its plain `row rhs` form).
pub fn read_linear_system_save(dir: &Path, prefix: &str, cols: Columns) -> io::Result<SavedSystem> {
    let text = std::fs::read_to_string(path(dir, prefix, "_a.dat"))?;
    let ncoef = cols.count();
    let mut entries: Vec<(usize, u32, [f64; 3])> = Vec::new();
    let (mut max_row, mut max_col) = (0usize, 0usize);
    for (ln, line) in text.lines().enumerate() {
        let mut it = line.split_whitespace();
        let Some(r) = it.next() else { continue };
        let r: usize = r
            .parse()
            .map_err(|_| bad(format!("_a.dat line {}: row", ln + 1)))?;
        let c: usize = it
            .next()
            .and_then(|t| t.parse().ok())
            .ok_or_else(|| bad(format!("_a.dat line {}: col", ln + 1)))?;
        let mut v = [0.0; 3];
        let mut got = 0;
        for t in it {
            if got >= 3 {
                return Err(bad(format!("_a.dat line {}: too many values", ln + 1)));
            }
            v[got] = parse_fortran_f64(t)
                .ok_or_else(|| bad(format!("_a.dat line {}: value {t}", ln + 1)))?;
            got += 1;
        }
        if got != ncoef {
            return Err(bad(format!(
                "_a.dat line {}: {got} values, expected {ncoef}",
                ln + 1
            )));
        }
        if r == 0 || c == 0 {
            return Err(bad(format!("_a.dat line {}: index below 1", ln + 1)));
        }
        max_row = max_row.max(r);
        max_col = max_col.max(c);
        entries.push((r - 1, (c - 1) as u32, v));
    }
    let mut nrows = max_row;
    if let Ok(s) = std::fs::read_to_string(path(dir, prefix, "_sizes.dat")) {
        if let Some(n) = s
            .split_whitespace()
            .next()
            .and_then(|t| t.parse::<usize>().ok())
        {
            nrows = nrows.max(n);
        }
    }
    let ncols = nrows.max(max_col);
    let mut row_ptr = vec![0usize; nrows + 1];
    for e in &entries {
        row_ptr[e.0 + 1] += 1;
    }
    for i in 0..nrows {
        row_ptr[i + 1] += row_ptr[i];
    }
    let nnz = entries.len();
    let mut col_idx = vec![0u32; nnz];
    let mut stiff = vec![0.0; nnz];
    let mut damp = vec![0.0; nnz];
    let mut mass = vec![0.0; nnz];
    let mut next = row_ptr.clone();
    for (r, c, v) in entries {
        let k = next[r];
        next[r] += 1;
        col_idx[k] = c;
        let mut idx = 0;
        if cols.stiff {
            stiff[k] = v[idx];
            idx += 1;
        }
        if cols.damp {
            damp[k] = v[idx];
            idx += 1;
        }
        if cols.mass {
            mass[k] = v[idx];
        }
    }
    let mut a = CsrMatrix {
        nrows,
        ncols,
        row_ptr,
        col_idx,
        values: stiff,
        diag: None,
    };
    a.diag = Some(find_diag(&a));

    let mut sys = SavedSystem {
        damp: cols.damp.then_some(damp),
        mass: cols.mass.then_some(mass),
        ..Default::default()
    };
    if let Ok(s) = std::fs::read_to_string(path(dir, prefix, "_b.dat")) {
        let mut rhs = vec![0.0; nrows];
        for (ln, line) in s.lines().enumerate() {
            let mut it = line.split_whitespace();
            let (Some(i), Some(v)) = (it.next(), it.next()) else {
                continue;
            };
            let i: usize = i
                .parse()
                .map_err(|_| bad(format!("_b.dat line {}", ln + 1)))?;
            if i == 0 || i > nrows {
                return Err(bad(format!("_b.dat line {}: row {i} out of range", ln + 1)));
            }
            rhs[i - 1] =
                parse_fortran_f64(v).ok_or_else(|| bad(format!("_b.dat line {}", ln + 1)))?;
        }
        sys.rhs = Some(rhs);
    }
    if let Ok(s) = std::fs::read_to_string(path(dir, prefix, "_perm.dat")) {
        let mut perm = Vec::new();
        for (ln, line) in s.lines().enumerate() {
            let mut it = line.split_whitespace();
            let (Some(_), Some(p)) = (it.next(), it.next()) else {
                continue;
            };
            let p: u64 = p
                .parse()
                .map_err(|_| bad(format!("_perm.dat line {}", ln + 1)))?;
            perm.push(if p == 0 { NO_DOF } else { (p - 1) as u32 });
        }
        sys.perm = Some(perm);
    }
    if let Ok(s) = std::fs::read_to_string(path(dir, prefix, "_sol.dat")) {
        let mut sol = Vec::new();
        for (ln, line) in s.lines().enumerate() {
            let mut it = line.split_whitespace();
            let (Some(_), Some(v)) = (it.next(), it.next()) else {
                continue;
            };
            sol.push(parse_fortran_f64(v).ok_or_else(|| bad(format!("_sol.dat line {}", ln + 1)))?);
        }
        sys.solution = Some(sol);
    }
    sys.stiff = Some(a);
    Ok(sys)
}
