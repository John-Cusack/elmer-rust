//! File formats of the `SaveData` solvers: `SaveScalars` (`<file>` + `<file>.names`)
//! and `SaveLine` (`<file>` + `<file>.names`). This module owns the layout
//! only (column names, ordering, number formatting, appending); the values
//! come from the operators that live elsewhere.
//!
//! References: `modules/SaveData/SaveScalars.F90` (names file ~1280-1360, rows
//! ~1380-1405, column names ~840-920) and `SaveLine.F90` (rows ~925 and
//! ~1110, names file ~2260-2320).
//!
//! Not implemented: the `.csv`, `.json` and `.marker` side files of
//! SaveScalars, per-partition files (`Parallel Reduce` off) and the eigen-value
//! column names.

use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::error::{Result, SaveIoError};
use crate::fortran::{fortran_es_e, fortran_i, now_elmer_date};

/// The header lines shared by both `.names` files.
#[derive(Clone, Debug, Default)]
pub struct NamesMeta {
    /// `Simulation` section `Comment`.
    pub simulation_comment: Option<String>,
    /// Solver section `Comment`.
    pub comment: Option<String>,
    /// `GetVersion()`.
    pub version: String,
    pub revision: Option<String>,
    pub compilation_date: Option<String>,
    /// `GetSifName()`.
    pub sif_name: Option<String>,
    /// `File started at:` text; the current UTC time if `None`.
    pub started_at: Option<String>,
}

fn names_preamble(meta: &NamesMeta, what: &str, file: &str) -> String {
    let mut s = String::new();
    for c in [&meta.simulation_comment, &meta.comment]
        .into_iter()
        .flatten()
    {
        s.push_str(c.trim_end());
        s.push('\n');
    }
    s.push_str(&format!("Metadata for {what} file: {file}\n"));
    s.push_str(&format!("Elmer version: {}\n", meta.version));
    if let Some(r) = &meta.revision {
        s.push_str(&format!("Elmer revision: {r}\n"));
    }
    if let Some(d) = &meta.compilation_date {
        s.push_str(&format!("Elmer compilation date: {d}\n"));
    }
    if let Some(n) = &meta.sif_name {
        s.push_str(&format!("Solver input file: {n}\n"));
    }
    let stamp = meta.started_at.clone().unwrap_or_else(now_elmer_date);
    s.push_str(&format!("File started at: {stamp}\n"));
    s
}

/// The scope suffix of an area/boundary operator's column name.
#[derive(Clone, Debug, PartialEq)]
pub enum ScalarScope {
    /// Whole domain: `max: temperature`.
    Domain,
    /// ` over bc <id>`.
    Bc(usize),
    /// ` over bc <mask name>`.
    BcMask(String),
    /// ` over polyline <j>`.
    Polyline(usize),
}

/// Column name of a SaveScalars value: `oper: variable[ over ...]`.
pub fn scalar_column_name(oper: &str, variable: &str, scope: &ScalarScope) -> String {
    let base = format!("{oper}: {variable}");
    match scope {
        ScalarScope::Domain => base,
        ScalarScope::Bc(j) => format!("{base} over bc {j}"),
        ScalarScope::BcMask(m) => format!("{base} over bc {m}"),
        ScalarScope::Polyline(j) => format!("{base} over polyline {j}"),
    }
}

/// Suffix the Fortran appends to a parallel-reduced column.
pub fn parallel_suffix(par_oper: &str) -> String {
    format!(" : mpi_{par_oper}")
}

/// Text of `<file>.names` for SaveScalars. `line_marker` is the `Line Marker`
/// keyword (`Some(n)` with `n != 0` adds a first column).
pub fn scalars_names_text(
    meta: &NamesMeta,
    file: &str,
    names: &[String],
    line_marker: Option<i32>,
) -> String {
    let mut s = names_preamble(meta, "SaveScalars", file);
    s.push_str(" \nVariables in columns of matrix: \n");
    let mut shift = 0;
    if line_marker.is_some_and(|m| m != 0) {
        s.push_str(&format!("{}: Line Marker\n", fortran_i(1, 4)));
        shift = 1;
    }
    for (i, n) in names.iter().enumerate() {
        s.push_str(&format!("{}: {n}\n", fortran_i((i + 1 + shift) as i64, 4)));
    }
    s
}

/// One data row of a SaveScalars file: optional `I10` line marker, then
/// `ES22.12E3` per value (`I10` for integer columns).
pub fn scalars_row_text(values: &[f64], integer: &[bool], line_marker: Option<i32>) -> String {
    let mut s = String::new();
    if let Some(m) = line_marker.filter(|&m| m != 0) {
        s.push_str(&fortran_i(m.into(), 10));
    }
    for (i, v) in values.iter().enumerate() {
        if integer.get(i).copied().unwrap_or(false) {
            s.push_str(&fortran_i(v.round() as i64, 10));
        } else {
            s.push_str(&fortran_es_e(*v, 22, 12, Some(3)));
        }
    }
    s.push('\n');
    s
}

/// Writes a SaveScalars file and its `.names` companion.
///
/// The names file is rewritten when the first row is written and whenever the
/// column set changes; rows after the first always append. `file_append`
/// (`File Append = True`) also appends to a file left by an earlier run.
pub struct ScalarsWriter {
    path: PathBuf,
    names_path: PathBuf,
    meta: NamesMeta,
    /// `./scalars.dat` as shown in the metadata line.
    display_name: String,
    line_marker: Option<i32>,
    file_append: bool,
    last_names: Option<Vec<String>>,
    rows: usize,
}

impl ScalarsWriter {
    /// `display_name` is the file name as the solver was given it
    /// (`Filename`, which the Fortran prefixes with `./` in the metadata).
    pub fn new(path: &Path, display_name: &str, meta: NamesMeta) -> Self {
        let mut np = path.as_os_str().to_owned();
        np.push(".names");
        ScalarsWriter {
            path: path.to_path_buf(),
            names_path: PathBuf::from(np),
            meta,
            display_name: display_name.to_string(),
            line_marker: None,
            file_append: false,
            last_names: None,
            rows: 0,
        }
    }

    pub fn with_line_marker(mut self, marker: i32) -> Self {
        self.line_marker = Some(marker);
        self
    }

    pub fn with_file_append(mut self, append: bool) -> Self {
        self.file_append = append;
        self
    }

    /// Write `<file>.names` for `names` now.
    pub fn write_names(&mut self, names: &[String]) -> Result<()> {
        std::fs::write(
            &self.names_path,
            scalars_names_text(&self.meta, &self.display_name, names, self.line_marker),
        )?;
        self.last_names = Some(names.to_vec());
        Ok(())
    }

    /// Append one row of `(column name, value)`; writes the names file first
    /// when the columns differ from the previous row.
    pub fn row(&mut self, values: &[(String, f64)]) -> Result<()> {
        let flagged: Vec<(String, f64, bool)> =
            values.iter().map(|(n, v)| (n.clone(), *v, false)).collect();
        self.row_flagged(&flagged)
    }

    /// Like [`row`](Self::row), with a flag per column for integer output.
    pub fn row_flagged(&mut self, values: &[(String, f64, bool)]) -> Result<()> {
        if values.is_empty() {
            return Err(SaveIoError::Input("SaveScalars row without values".into()));
        }
        let names: Vec<String> = values.iter().map(|v| v.0.clone()).collect();
        if self.last_names.as_ref() != Some(&names) {
            self.write_names(&names)?;
        }
        let nums: Vec<f64> = values.iter().map(|v| v.1).collect();
        let ints: Vec<bool> = values.iter().map(|v| v.2).collect();
        let text = scalars_row_text(&nums, &ints, self.line_marker);
        let mut f = if self.rows > 0 || self.file_append {
            OpenOptions::new()
                .append(true)
                .create(true)
                .open(&self.path)?
        } else {
            File::create(&self.path)?
        };
        f.write_all(text.as_bytes())?;
        self.rows += 1;
        Ok(())
    }
}

// ---------------------------------------------------------------- SaveLine

/// Node counts reported in the `SaveLine` names file.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct LineNodeCounts {
    pub bc: usize,
    pub polyline: usize,
    pub circle: usize,
    pub isocurve: usize,
}

/// Options that decide the leading columns of a `SaveLine` file.
#[derive(Copy, Clone, Debug, Default)]
pub struct LineColumns {
    /// `Line Marker` is nonzero: first column.
    pub line_marker: bool,
    /// Transient simulation: a `Timestep` column.
    pub timestep: bool,
    /// `Skip Boundary Info`: drops `Call count`, `Boundary condition`, `Node index`.
    pub skip_boundary_info: bool,
}

impl LineColumns {
    /// Number of label columns before the values.
    pub fn n_labels(&self) -> usize {
        if self.skip_boundary_info {
            0
        } else {
            3 + usize::from(self.line_marker) + usize::from(self.timestep)
        }
    }
}

/// Text of `<file>.names` for SaveLine.
pub fn line_names_text(
    meta: &NamesMeta,
    file: &str,
    counts: LineNodeCounts,
    cols: LineColumns,
    value_names: &[String],
) -> String {
    let mut s = names_preamble(meta, "SaveLine", file);
    s.push_str("Number of data nodes for each step\n");
    s.push_str(&format!("  bc nodes: {}\n", counts.bc));
    s.push_str(&format!("  polyline nodes: {}\n", counts.polyline));
    s.push_str(&format!("  circle nodes: {}\n", counts.circle));
    s.push_str(&format!("  isocurve nodes: {}\n", counts.isocurve));
    s.push_str("Data on different columns\n");
    let mut j = 0;
    let mut push = |s: &mut String, name: &str| {
        j += 1;
        s.push_str(&format!("{}: {name}\n", fortran_i(j, 3)));
    };
    if !cols.skip_boundary_info {
        if cols.line_marker {
            push(&mut s, "Line Marker");
        }
        if cols.timestep {
            push(&mut s, "Timestep");
        }
        push(&mut s, "Call count");
        push(&mut s, "Boundary condition");
        push(&mut s, "Node index");
    }
    for n in value_names {
        push(&mut s, n);
    }
    s
}

/// One `SaveLine` row: each label as `<int> `, then `ES20.11E3` per value.
pub fn line_row_text(labels: &[i64], values: &[f64]) -> String {
    let mut s = String::new();
    for l in labels {
        s.push_str(&format!("{l} "));
    }
    for v in values {
        s.push_str(&fortran_es_e(*v, 20, 11, Some(3)));
    }
    s.push('\n');
    s
}

/// Append-or-create a text file the way the SaveData solvers do.
pub fn append_text(path: &Path, text: &str, append: bool) -> Result<()> {
    let mut f = if append {
        OpenOptions::new().append(true).create(true).open(path)?
    } else {
        File::create(path)?
    };
    f.write_all(text.as_bytes())?;
    Ok(())
}
