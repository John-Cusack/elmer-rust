//! Writing `.result` files from live fields (`SaveResult`), and the file-name
//! conventions around it (cyclic and separate files, `_last_nc`).
//!
//! Layout rules reproduced from the Fortran writer:
//!
//! * the variable list holds vector fields first, then scalar fields; only
//!   scalar fields carry data, so a vector field `name` with `D` dofs
//!   contributes `D` scalar components named by [`component_name`], read from
//!   its interleaved `values` (`values[slot * D + c]`);
//! * scalar order is the order of [`FieldSource::fields`], with the components
//!   of a vector field at the position of that field;
//! * a permutation equal to the last *written* table is stored as
//!   `Perm: use previous`;
//! * dof slots in [`FieldRef::perm`] are 0-based ([`NO_DOF`] = absent); the
//!   file stores node and slot numbers 1-based.
//!
//! Deliberate differences from `SaveResult`, all of which only avoid writing
//! a file that the reader cannot parse: `Output Variable i` selection uses the
//! header rule in both passes (the Fortran data pass uses an exact-name test
//! that drops numbered components the header listed), and `Omit unchanged
//! variables` is not modelled (every save is a full save).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use elmer_core::field::{FieldKind, FieldRef, NO_DOF};
use elmer_core::{FieldSource, MeshView};

use crate::error::{Result, SaveIoError};
use crate::fortran::now_elmer_date;
use crate::result::{
    variable_listed, PermRecord, ResultEncoder, ResultFormat, ResultHeader, StepRecord, VarHeader,
    VarRecord,
};

/// The `Simulation` keywords that steer `SaveResult`.
#[derive(Clone, Debug, Default)]
pub struct ResultOptions {
    /// `Output Global Variables`.
    pub save_global: bool,
    /// `Output Coordinates`.
    pub save_coordinates: bool,
    /// `Output IP Variables` (fields of kind [`FieldKind::Ip`]).
    pub save_ip: bool,
    /// `Output Variable 1..n`. Empty selects everything.
    pub output_variables: Vec<String>,
    /// Field name -> solver `Equation`; `no equation` when absent. Looked up
    /// by the full field name, then by the name before `[`.
    pub equations: HashMap<String, String>,
    /// Text of the `!File started at:` stamp; the current UTC time if `None`.
    pub started_at: Option<String>,
}

/// Name of component `c` (1-based) of a vector field: `ComponentName` in
/// `GeneralUtils.F90`. `velocity:3 pressure:1` gives `velocity 1..3`, `pressure`.
pub fn component_name(base: &str, c: usize, dofs: usize) -> Result<String> {
    if base == "flow solution" {
        return Ok(if c == dofs {
            "pressure".into()
        } else {
            format!("velocity {c}")
        });
    }
    let Some(open) = base.find('[') else {
        return Ok(format!("{base} {c}"));
    };
    let mut total = 0usize;
    let mut rest = &base[open + 1..];
    loop {
        let colon = rest.find(':').ok_or_else(|| {
            SaveIoError::Input(format!("syntax error in variable definition '{base}'"))
        })?;
        let part = rest[colon + 1..]
            .chars()
            .next()
            .and_then(|d| d.to_digit(10))
            .ok_or_else(|| {
                SaveIoError::Input(format!("syntax error in variable definition '{base}'"))
            })? as usize;
        total += part;
        if total >= c {
            let label = rest[..colon].trim_start();
            return Ok(if part > 1 {
                format!("{label} {}", c + part - total)
            } else {
                label.to_string()
            });
        }
        rest = rest[colon + 2..].trim_start_matches(' ');
        if rest.is_empty() {
            return Err(SaveIoError::Input(format!(
                "component {c} outside '{base}'"
            )));
        }
    }
}

/// One saved scalar or vector as the file sees it.
struct View<'a> {
    name: String,
    dofs: usize,
    size: usize,
    perm: Option<&'a [u32]>,
    values: &'a [f64],
    offset: usize,
    stride: usize,
    equation: String,
}

impl View<'_> {
    fn at(&self, slot: usize) -> f64 {
        self.values[self.offset + slot * self.stride]
    }
}

struct Selection<'a> {
    vectors: Vec<View<'a>>,
    scalars: Vec<View<'a>>,
}

fn keep(v: &View<'_>, opts: &ResultOptions) -> bool {
    let save = if v.size == v.dofs {
        opts.save_global
    } else {
        !v.name.to_ascii_lowercase().starts_with("coordinate") || opts.save_coordinates
    };
    save && variable_listed(&opts.output_variables, &v.name)
}

fn select<'a>(fields: &[FieldRef<'a>], opts: &ResultOptions) -> Result<Selection<'a>> {
    let mut sel = Selection {
        vectors: Vec::new(),
        scalars: Vec::new(),
    };
    for f in fields {
        if f.dofs == 0 || f.values.len() % f.dofs != 0 {
            return Err(SaveIoError::Input(format!(
                "field '{}': {} values for {} dofs",
                f.name,
                f.values.len(),
                f.dofs
            )));
        }
        if f.kind == FieldKind::Ip {
            if !opts.save_ip {
                continue;
            }
            if f.perm.is_some() {
                return Err(SaveIoError::Unsupported(format!(
                    "IP field '{}' with a permutation",
                    f.name
                )));
            }
        }
        if let Some(p) = f.perm {
            if let Some(&bad) = p
                .iter()
                .find(|&&s| s != NO_DOF && s as usize >= f.values.len() / f.dofs)
            {
                return Err(SaveIoError::Input(format!(
                    "field '{}': dof slot {bad} out of range",
                    f.name
                )));
            }
        }
        let equation = opts
            .equations
            .get(f.name)
            .or_else(|| f.name.split('[').next().and_then(|b| opts.equations.get(b)))
            .cloned()
            .unwrap_or_else(|| "no equation".into());
        let whole = View {
            name: f.name.to_string(),
            dofs: f.dofs,
            size: f.values.len(),
            perm: f.perm,
            values: f.values,
            offset: 0,
            stride: 1,
            equation: equation.clone(),
        };
        if f.dofs == 1 {
            if keep(&whole, opts) {
                sel.scalars.push(whole);
            }
            continue;
        }
        if keep(&whole, opts) {
            sel.vectors.push(View {
                name: whole.name.clone(),
                equation: equation.clone(),
                ..whole
            });
        }
        for c in 0..f.dofs {
            let comp = View {
                name: component_name(f.name, c + 1, f.dofs)?,
                dofs: 1,
                size: f.values.len() / f.dofs,
                perm: f.perm,
                values: f.values,
                offset: c,
                stride: f.dofs,
                equation: equation.clone(),
            };
            if keep(&comp, opts) {
                sel.scalars.push(comp);
            }
        }
    }
    Ok(sel)
}

fn header_of(
    sel: &Selection<'_>,
    n_nodes: usize,
    format: ResultFormat,
    started_at: String,
) -> ResultHeader {
    let line = |v: &View<'_>| VarHeader {
        name: v.name.clone(),
        field_size: v.size,
        perm_size: v.perm.map_or(0, <[u32]>::len),
        dofs: v.dofs,
        equation: v.equation.clone(),
    };
    ResultHeader {
        format,
        started_at,
        vars: sel.vectors.iter().chain(&sel.scalars).map(line).collect(),
        n_nodes,
    }
}

fn step_of(sel: &Selection<'_>, saved_count: i32, timestep: i32, time: f64) -> StepRecord {
    let mut prev: Option<&[u32]> = None;
    let mut vars = Vec::with_capacity(sel.scalars.len());
    for v in &sel.scalars {
        let (perm, values) = match v.perm {
            None => (PermRecord::Null, (0..v.size).map(|s| v.at(s)).collect()),
            Some(p) if prev == Some(p) => (PermRecord::UsePrevious, table_values(v, p)),
            Some(p) => {
                prev = Some(p);
                let entries = p
                    .iter()
                    .enumerate()
                    .filter(|(_, &s)| s != NO_DOF)
                    .map(|(i, &s)| (i as u32 + 1, s + 1))
                    .collect();
                (
                    PermRecord::Table {
                        size: p.len(),
                        entries,
                    },
                    table_values(v, p),
                )
            }
        };
        vars.push(VarRecord {
            name: v.name.clone(),
            perm,
            values,
        });
    }
    StepRecord {
        saved_count,
        timestep,
        time,
        vars,
    }
}

fn table_values(v: &View<'_>, p: &[u32]) -> Vec<f64> {
    p.iter()
        .filter(|&&s| s != NO_DOF)
        .map(|&s| v.at(s as usize))
        .collect()
}

/// Appends saves to one `.result` file (and its `.pos` index for binary),
/// like repeated calls of `SaveResult` on one mesh.
pub struct ResultWriter {
    format: ResultFormat,
    opts: ResultOptions,
    enc: ResultEncoder,
    header: Option<ResultHeader>,
    saves_done: usize,
    main: BufWriter<File>,
    pos: Option<BufWriter<File>>,
}

impl ResultWriter {
    /// Create (truncate) `path`, and `path.pos` for binary.
    pub fn create(path: &Path, format: ResultFormat, opts: ResultOptions) -> Result<Self> {
        let main = BufWriter::new(File::create(path)?);
        let pos = match format {
            ResultFormat::Ascii => None,
            ResultFormat::Binary(_) => Some(BufWriter::new(File::create(pos_path(path))?)),
        };
        Ok(ResultWriter {
            format,
            opts,
            enc: ResultEncoder::new(format),
            header: None,
            saves_done: 0,
            main,
            pos,
        })
    }

    /// Number of saves written so far.
    pub fn saves_done(&self) -> usize {
        self.saves_done
    }

    /// Append the current state of `src` and return the save count.
    pub fn save(&mut self, mesh: &dyn MeshView, src: &dyn FieldSource) -> Result<usize> {
        let fields = src.fields();
        let sel = select(&fields, &self.opts)?;
        let timestep =
            i32::try_from(src.step()).map_err(|_| SaveIoError::Input("step exceeds i32".into()))?;
        match &self.header {
            None => {
                let stamp = self.opts.started_at.clone().unwrap_or_else(now_elmer_date);
                let h = header_of(&sel, mesh.n_nodes(), self.format, stamp);
                let (m, p) = self.enc.header(&h)?;
                self.emit(&m, p.as_deref())?;
                self.header = Some(h);
            }
            Some(h) => {
                let now = header_of(&sel, mesh.n_nodes(), self.format, h.started_at.clone());
                if now != *h {
                    return Err(SaveIoError::Input(
                        "the set of saved variables changed since the first save".into(),
                    ));
                }
            }
        }
        let rec = step_of(&sel, self.saves_done as i32 + 1, timestep, src.time());
        let (m, p) = self.enc.step(&rec)?;
        self.emit(&m, p.as_deref())?;
        self.saves_done += 1;
        Ok(self.saves_done)
    }

    fn emit(&mut self, main: &[u8], pos: Option<&[u8]>) -> Result<()> {
        self.main.write_all(main)?;
        if let (Some(w), Some(p)) = (self.pos.as_mut(), pos) {
            w.write_all(p)?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.main.flush()?;
        if let Some(p) = self.pos.as_mut() {
            p.flush()?;
        }
        Ok(())
    }
}

fn pos_path(path: &Path) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".pos");
    PathBuf::from(s)
}

/// Write one save to a fresh file with default options.
pub fn write_result(
    path: &Path,
    mesh: &dyn MeshView,
    fields: &dyn FieldSource,
    format: ResultFormat,
) -> Result<()> {
    write_result_with(path, mesh, fields, format, ResultOptions::default())
}

pub fn write_result_with(
    path: &Path,
    mesh: &dyn MeshView,
    fields: &dyn FieldSource,
    format: ResultFormat,
    opts: ResultOptions,
) -> Result<()> {
    let mut w = ResultWriter::create(path, format, opts)?;
    w.save(mesh, fields)?;
    w.finish()
}

// ------------------------------------------------------------ file naming

/// `Output File Cycle = n`: file name of save number `saves_done` (0-based).
/// Returns the name and the 1-based file index.
pub fn cycle_file_name(base: &str, file_cycle: usize, saves_done: usize) -> (String, usize) {
    let ind = saves_done % file_cycle + 1;
    (format!("{base}_{ind}nc"), ind)
}

/// `Separate Results Files`: the `calve_ind`-th file (counting from 0).
pub fn separate_file_name(base: &str, calve_ind: usize) -> String {
    format!("{base}_{calve_ind}nc")
}

/// Per-partition name in a parallel run.
pub fn partition_file_name(base: &str, my_pe: usize) -> String {
    format!("{base}.{my_pe}")
}

/// Contents of the `<file>_last_nc` marker written after each cyclic save.
pub fn last_nc_contents(date: &str, file_index: usize) -> String {
    format!("!File created at: {date}\n!Last Saved File Cycle:\n{file_index}\n")
}

/// Cycle index stored in a `_last_nc` file (`Restart File Cycle = 0`).
pub fn parse_last_nc(contents: &str) -> Result<usize> {
    contents
        .lines()
        .nth(2)
        .and_then(|l| l.trim().parse().ok())
        .ok_or_else(|| SaveIoError::Format("bad _last_nc file".into()))
}

/// Restart file name for `Restart File Cycle = j`, with the position rule of
/// `LoadRestartFile` (`index of the first '.' + 6`, which is the end of
/// `case.result`).
pub fn restart_cycle_name(restart_file: &str, j: usize) -> Result<String> {
    let dot = restart_file.find('.').map_or(0, |i| i + 1);
    let at = dot + 6;
    if at > restart_file.len() || !restart_file.is_char_boundary(at) {
        return Err(SaveIoError::Input(format!(
            "'{restart_file}' is too short for Restart File Cycle"
        )));
    }
    Ok(format!(
        "{}_{j}nc{}",
        &restart_file[..at],
        &restart_file[at..]
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn component_names() {
        let n = "flow solution[velocity:3 pressure:1]";
        let got: Vec<_> = (1..=4).map(|c| component_name(n, c, 4).unwrap()).collect();
        assert_eq!(got, ["velocity 1", "velocity 2", "velocity 3", "pressure"]);
        assert_eq!(
            component_name("displacement", 2, 3).unwrap(),
            "displacement 2"
        );
        assert_eq!(component_name("flow solution", 3, 3).unwrap(), "pressure");
        assert_eq!(component_name("a[b:2 c:1 d:2]", 5, 5).unwrap(), "d 2");
        assert!(component_name("a[b:2]", 3, 2).is_err());
    }

    #[test]
    fn naming() {
        assert_eq!(
            cycle_file_name("case.result", 3, 4),
            ("case.result_2nc".into(), 2)
        );
        assert_eq!(
            restart_cycle_name("case.result", 2).unwrap(),
            "case.result_2nc"
        );
        // The Fortran rule cuts at the first '.', so a leading "../" garbles the name.
        assert_eq!(
            restart_cycle_name("../x.result", 2).unwrap(),
            "../x.re_2ncsult"
        );
        assert_eq!(
            parse_last_nc(&last_nc_contents("2026/01/01 00:00:00", 7)).unwrap(),
            7
        );
    }
}
