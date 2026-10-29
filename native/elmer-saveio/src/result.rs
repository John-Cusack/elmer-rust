//! `.result` restart files: the data model, the byte-exact encoder and the
//! parser (Fortran `SaveResult` / `LoadRestartFile`, `ModelDescription.F90`).
//!
//! A file is a header (format line, `!File started at:` stamp, one line per
//! variable, `Total DOFs`, `Number Of Nodes`) followed by time steps. Each step
//! stores every *scalar* variable (vector fields are written component-wise)
//! as name, permutation and values. Binary files also have a `<file>.pos`
//! index which [`ResultEncoder`] regenerates.
//!
//! Versions 1 and 2 of the format (obsolete) are rejected with
//! [`SaveIoError::Unsupported`].

use crate::error::{format_err, Result, SaveIoError};
use crate::fortran::{fortran_es_e, fortran_i, list_i32, list_real, parse_list_real};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Endian {
    Little,
    Big,
}

impl Endian {
    pub const NATIVE: Endian = if cfg!(target_endian = "little") {
        Endian::Little
    } else {
        Endian::Big
    };

    fn letter(self) -> u8 {
        match self {
            Endian::Little => b'L',
            Endian::Big => b'B',
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ResultFormat {
    Ascii,
    Binary(Endian),
}

impl ResultFormat {
    /// Binary in the byte order of this machine, as the Fortran writer does.
    pub const BINARY: ResultFormat = ResultFormat::Binary(Endian::NATIVE);
}

/// One line of the variable list.
#[derive(Clone, Debug, PartialEq)]
pub struct VarHeader {
    /// Full name, including a `[velocity:3 pressure:1]` suffix when present.
    pub name: String,
    /// `SIZE(Var % Values)`.
    pub field_size: usize,
    /// `SIZE(Var % Perm)`, 0 without a permutation.
    pub perm_size: usize,
    pub dofs: usize,
    /// Solver `Equation`, `no equation` when none.
    pub equation: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResultHeader {
    pub format: ResultFormat,
    /// Text after `!File started at: `.
    pub started_at: String,
    /// Vector variables first, then scalar variables, as in the file.
    pub vars: Vec<VarHeader>,
    pub n_nodes: usize,
}

impl ResultHeader {
    /// `Total DOFs`: the number of scalar (single component) variables.
    pub fn total_dofs(&self) -> usize {
        self.vars.iter().filter(|v| v.dofs == 1).count()
    }

    fn scalars(&self) -> impl Iterator<Item = &VarHeader> {
        self.vars.iter().filter(|v| v.dofs == 1)
    }
}

/// A saved permutation.
#[derive(Clone, Debug, PartialEq)]
pub enum PermRecord {
    /// `Perm: NULL`: one value per entry of the field.
    Null,
    /// Same table as the previous variable of this step.
    UsePrevious,
    /// `size` = `SIZE(Perm)`; `(node, slot)` pairs, both 1-based as in the
    /// file, for every nonzero entry.
    Table {
        size: usize,
        entries: Vec<(u32, u32)>,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct VarRecord {
    pub name: String,
    pub perm: PermRecord,
    /// For a table: one value per entry, in entry order (`Values(Perm(i))`).
    /// Otherwise the whole field.
    pub values: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct StepRecord {
    /// 1-based count of saves in this file.
    pub saved_count: i32,
    pub timestep: i32,
    pub time: f64,
    pub vars: Vec<VarRecord>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResultFile {
    pub header: ResultHeader,
    pub steps: Vec<StepRecord>,
}

/// A variable of one step resolved the way `LoadRestartFile` assigns it:
/// `values[i]` belongs to `node[i]`, stored in dof `slot[i]` (all 0-based).
#[derive(Clone, Debug, PartialEq)]
pub struct RestartField {
    pub node: Vec<u32>,
    pub slot: Vec<u32>,
    pub values: Vec<f64>,
}

impl ResultFile {
    /// Variable `name` of step `step`. `None` when the step does not hold it;
    /// an error when it uses `Perm: use previous` with nothing before it.
    pub fn restart_field(&self, step: usize, name: &str) -> Result<Option<RestartField>> {
        let Some(st) = self.steps.get(step) else {
            return Ok(None);
        };
        let mut table: Option<&Vec<(u32, u32)>> = None;
        for v in &st.vars {
            let this = match &v.perm {
                PermRecord::Table { entries, .. } => {
                    table = Some(entries);
                    Some(entries)
                }
                PermRecord::UsePrevious => Some(table.ok_or_else(|| {
                    SaveIoError::Format(format!(
                        "'{}': Perm: use previous with no previous table",
                        v.name
                    ))
                })?),
                PermRecord::Null => None,
            };
            if v.name != name {
                continue;
            }
            return Ok(Some(match this {
                Some(entries) => RestartField {
                    node: entries.iter().map(|e| e.0 - 1).collect(),
                    slot: entries.iter().map(|e| e.1 - 1).collect(),
                    values: v.values.clone(),
                },
                None => RestartField {
                    node: (0..v.values.len() as u32).collect(),
                    slot: (0..v.values.len() as u32).collect(),
                    values: v.values.clone(),
                },
            }));
        }
        Ok(None)
    }

    /// Index of the step whose `saved_count` is `count`, or the last step for
    /// `count == 0` (the `Restart Position = 0` convention).
    pub fn find_step(&self, count: i32) -> Option<usize> {
        if count == 0 {
            self.steps.len().checked_sub(1)
        } else {
            self.steps.iter().position(|s| s.saved_count == count)
        }
    }

    /// Serialize to `(main file, .pos index)`; the index exists for binary only.
    pub fn to_bytes(&self) -> Result<(Vec<u8>, Option<Vec<u8>>)> {
        let mut enc = ResultEncoder::new(self.header.format);
        let (mut main, mut pos) = enc.header(&self.header)?;
        for s in &self.steps {
            let (m, p) = enc.step(s)?;
            main.extend(m);
            if let (Some(pos), Some(p)) = (pos.as_mut(), p) {
                pos.extend(p);
            }
        }
        Ok((main, pos))
    }

    /// Parse a main file. `pos` (the `.pos` bytes) is not needed: every step is
    /// self-delimiting.
    pub fn from_bytes(bytes: &[u8]) -> Result<ResultFile> {
        parse(bytes)
    }
}

// ---------------------------------------------------------------- encoding

fn var_line(h: &VarHeader) -> String {
    let mut s = String::new();
    s.push_str(&h.name);
    if h.name.len() < 25 {
        while s.len() < 24 {
            s.push(' ');
        }
    }
    s.push_str(" : ");
    s.push_str(&fortran_i(h.field_size as i64, 8));
    s.push_str(&fortran_i(h.perm_size as i64, 8));
    s.push_str(&fortran_i(h.dofs as i64, 4));
    s.push_str(" : ");
    s.push_str(&h.equation);
    s
}

/// Stateful encoder: header once, then one call per saved step, mirroring
/// the append behaviour of `SaveResult`.
pub struct ResultEncoder {
    format: ResultFormat,
    /// Bytes of the main file produced so far.
    offset: u64,
    n_scalars: Option<usize>,
}

impl ResultEncoder {
    pub fn new(format: ResultFormat) -> Self {
        ResultEncoder {
            format,
            offset: 0,
            n_scalars: None,
        }
    }

    /// Bytes of the main file written so far.
    pub fn offset(&self) -> u64 {
        self.offset
    }

    pub fn header(&mut self, h: &ResultHeader) -> Result<(Vec<u8>, Option<Vec<u8>>)> {
        if h.format != self.format {
            return Err(SaveIoError::Input(
                "header format differs from encoder format".into(),
            ));
        }
        let mut t = String::new();
        match self.format {
            ResultFormat::Ascii => t.push_str(" ASCII 3\n"),
            ResultFormat::Binary(e) => t.push_str(&format!(" BINARY 3.{}\n", e.letter() as char)),
        }
        t.push_str(&format!("!File started at: {}\n", h.started_at));
        t.push_str(" Degrees of freedom: \n");
        for pass_vector in [true, false] {
            for v in h.vars.iter().filter(|v| (v.dofs > 1) == pass_vector) {
                t.push_str(&var_line(v));
                t.push('\n');
            }
        }
        t.push_str(&format!(
            " Total DOFs: {}\n",
            list_i32(h.total_dofs() as i32)
        ));
        t.push_str(&format!(
            " Number Of Nodes: {}\n",
            list_i32(h.n_nodes as i32)
        ));
        let mut main = t.into_bytes();
        self.n_scalars = Some(h.total_dofs());
        let pos = match self.format {
            ResultFormat::Ascii => None,
            ResultFormat::Binary(e) => {
                main.push(0); // SwitchToBinary: an empty string, i.e. a NUL byte
                let mut p = vec![e.letter()];
                p.extend(put_i32(e, h.total_dofs() as i32));
                for v in h.scalars() {
                    p.extend(v.name.as_bytes());
                    p.push(0);
                }
                Some(p)
            }
        };
        self.offset += main.len() as u64;
        Ok((main, pos))
    }

    pub fn step(&mut self, s: &StepRecord) -> Result<(Vec<u8>, Option<Vec<u8>>)> {
        let n = self
            .n_scalars
            .ok_or_else(|| SaveIoError::Input("step before header".into()))?;
        if s.vars.len() > n {
            return Err(SaveIoError::Input(format!(
                "{} variables, header has {n}",
                s.vars.len()
            )));
        }
        match self.format {
            ResultFormat::Ascii => Ok((self.step_ascii(s)?, None)),
            ResultFormat::Binary(e) => self.step_binary(e, s).map(|(m, p)| (m, Some(p))),
        }
    }

    fn step_ascii(&mut self, s: &StepRecord) -> Result<Vec<u8>> {
        let mut t = format!(
            "Time: {}{}{}\n",
            fortran_i(s.saved_count.into(), 7),
            fortran_i(s.timestep.into(), 7),
            fortran_es_e(s.time, 17, 8, Some(3))
        );
        let mut seen_table = false;
        for v in &s.vars {
            t.push_str(&v.name);
            t.push('\n');
            match &v.perm {
                PermRecord::Null => t.push_str("Perm: NULL\n"),
                PermRecord::UsePrevious => {
                    if !seen_table {
                        return Err(SaveIoError::Input(format!(
                            "'{}': use previous without a table",
                            v.name
                        )));
                    }
                    t.push_str("Perm: use previous\n");
                }
                PermRecord::Table { size, entries } => {
                    seen_table = true;
                    t.push_str(&format!(
                        "Perm: {} {}\n",
                        fortran_i(*size as i64, 12),
                        fortran_i(entries.len() as i64, 12)
                    ));
                    for (i, k) in entries {
                        t.push_str(&format!(
                            "{}{}\n",
                            fortran_i((*i).into(), 11),
                            fortran_i((*k).into(), 11)
                        ));
                    }
                }
            }
            for x in &v.values {
                t.push_str(&list_real(*x));
                t.push('\n');
            }
        }
        self.offset += t.len() as u64;
        Ok(t.into_bytes())
    }

    fn step_binary(&mut self, e: Endian, s: &StepRecord) -> Result<(Vec<u8>, Vec<u8>)> {
        let base = self.offset;
        let mut m: Vec<u8> = Vec::new();
        let mut p: Vec<u8> = Vec::new();
        p.extend(put_i64(e, (base + m.len() as u64) as i64));
        m.extend(b"Time:\0");
        m.extend(put_i32(e, s.saved_count));
        m.extend(put_i32(e, s.timestep));
        m.extend(put_f64(e, s.time));
        let mut prev_table_pos: Option<u64> = None;
        for v in &s.vars {
            p.extend(put_i64(e, (base + m.len() as u64) as i64));
            m.extend(v.name.as_bytes());
            m.push(0);
            match &v.perm {
                PermRecord::Null => m.extend(put_i32(e, 0)),
                PermRecord::UsePrevious => {
                    let at = prev_table_pos.ok_or_else(|| {
                        SaveIoError::Input(format!("'{}': use previous without a table", v.name))
                    })?;
                    m.extend(put_i32(e, -1));
                    m.extend(put_i64(e, at as i64));
                }
                PermRecord::Table { size, entries } => {
                    prev_table_pos = Some(base + m.len() as u64);
                    m.extend(put_i32(e, *size as i32));
                    m.extend(put_i32(e, entries.len() as i32));
                    for (i, k) in entries {
                        m.extend(put_i32(e, *i as i32));
                        m.extend(put_i32(e, *k as i32));
                    }
                }
            }
            for x in &v.values {
                m.extend(put_f64(e, *x));
            }
        }
        self.offset += m.len() as u64;
        Ok((m, p))
    }
}

fn put_i32(e: Endian, x: i32) -> [u8; 4] {
    if e == Endian::Little {
        x.to_le_bytes()
    } else {
        x.to_be_bytes()
    }
}
fn put_i64(e: Endian, x: i64) -> [u8; 8] {
    if e == Endian::Little {
        x.to_le_bytes()
    } else {
        x.to_be_bytes()
    }
}
fn put_f64(e: Endian, x: f64) -> [u8; 8] {
    if e == Endian::Little {
        x.to_le_bytes()
    } else {
        x.to_be_bytes()
    }
}

// ----------------------------------------------------------------- parsing

struct Cursor<'a> {
    b: &'a [u8],
    at: usize,
    e: Endian,
}

impl<'a> Cursor<'a> {
    fn eof(&self) -> bool {
        self.at >= self.b.len()
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.at + n > self.b.len() {
            return format_err("unexpected end of binary data");
        }
        let s = &self.b[self.at..self.at + n];
        self.at += n;
        Ok(s)
    }
    fn i32(&mut self) -> Result<i32> {
        let a: [u8; 4] = self.take(4)?.try_into().expect("4 bytes");
        Ok(if self.e == Endian::Little {
            i32::from_le_bytes(a)
        } else {
            i32::from_be_bytes(a)
        })
    }
    fn i64(&mut self) -> Result<i64> {
        let a: [u8; 8] = self.take(8)?.try_into().expect("8 bytes");
        Ok(if self.e == Endian::Little {
            i64::from_le_bytes(a)
        } else {
            i64::from_be_bytes(a)
        })
    }
    fn f64(&mut self) -> Result<f64> {
        let a: [u8; 8] = self.take(8)?.try_into().expect("8 bytes");
        Ok(if self.e == Endian::Little {
            f64::from_le_bytes(a)
        } else {
            f64::from_be_bytes(a)
        })
    }
    /// NUL-terminated string (the NUL is consumed).
    fn string(&mut self) -> Result<String> {
        let rest = &self.b[self.at..];
        let n = rest
            .iter()
            .position(|&c| c == 0)
            .ok_or_else(|| SaveIoError::Format("unterminated string".into()))?;
        let s = String::from_utf8_lossy(&rest[..n]).into_owned();
        self.at += n + 1;
        Ok(s)
    }
}

fn parse(bytes: &[u8]) -> Result<ResultFile> {
    let first_nl = bytes
        .iter()
        .position(|&c| c == b'\n')
        .ok_or_else(|| SaveIoError::Format("no header line".into()))?;
    let first = String::from_utf8_lossy(&bytes[..first_nl]).into_owned();
    let t = first.trim_start();
    let (format, version) = if let Some(r) = t.strip_prefix("ASCII") {
        (ResultFormat::Ascii, r.trim().to_string())
    } else if let Some(r) = t.strip_prefix("BINARY") {
        let r = r.trim();
        let (v, e) = r
            .split_once('.')
            .ok_or_else(|| SaveIoError::Format(format!("bad binary header '{first}'")))?;
        let e = match e.trim().chars().next() {
            Some('L') => Endian::Little,
            Some('B') => Endian::Big,
            _ => return format_err(format!("bad endianness in '{first}'")),
        };
        (ResultFormat::Binary(e), v.trim().to_string())
    } else {
        return format_err(format!("could not determine file format from '{first}'"));
    };
    if version != "3" {
        return Err(SaveIoError::Unsupported(format!(
            "result file version '{version}' (only 3)"
        )));
    }

    // The text header ends at "Number Of Nodes" (and, for binary, a NUL).
    let (text_end, body_start) = match format {
        ResultFormat::Ascii => (bytes.len(), 0),
        ResultFormat::Binary(_) => {
            let n = bytes
                .iter()
                .position(|&c| c == 0)
                .ok_or_else(|| SaveIoError::Format("no binary section".into()))?;
            (n, n + 1)
        }
    };
    let text = String::from_utf8_lossy(&bytes[..text_end]).into_owned();
    let mut lines = text.split('\n').peekable();
    lines.next(); // format line
    let mut started_at = String::new();
    let mut vars: Vec<VarHeader> = Vec::new();
    let mut n_nodes = 0usize;
    let mut total_dofs = None;
    let mut in_vars = false;
    let mut consumed = first_nl + 1;
    for raw in lines.by_ref() {
        consumed += raw.len() + 1;
        let line = raw.trim().to_lowercase();
        if let Some(rest) = raw.trim_start().strip_prefix("!File started at:") {
            started_at = rest
                .strip_prefix(' ')
                .unwrap_or(rest)
                .trim_end()
                .to_string();
            continue;
        }
        if line.is_empty() || line.starts_with('!') {
            continue;
        }
        if !in_vars {
            if line.contains("degrees of freedom") {
                in_vars = true;
            }
            continue;
        }
        if let Some(k) = line.rfind("total dofs:") {
            total_dofs = Some(
                line[k + 11..]
                    .trim()
                    .parse::<usize>()
                    .map_err(|_| SaveIoError::Format(format!("bad Total DOFs '{raw}'")))?,
            );
            continue;
        }
        if let Some(r) = line.strip_prefix("number of nodes:") {
            n_nodes = r
                .trim()
                .parse()
                .map_err(|_| SaveIoError::Format(format!("bad Number Of Nodes '{raw}'")))?;
            if format == ResultFormat::Ascii {
                break;
            }
            continue;
        }
        if total_dofs.is_some() {
            continue;
        }
        vars.push(parse_var_line(&line)?);
    }
    let header = ResultHeader {
        format,
        started_at,
        vars,
        n_nodes,
    };
    let td = total_dofs.ok_or_else(|| SaveIoError::Format("no 'Total DOFs' line".into()))?;
    if td != header.total_dofs() {
        return format_err(format!(
            "Total DOFs {td} but {} scalar variables listed",
            header.total_dofs()
        ));
    }
    let sizes: Vec<usize> = header.scalars().map(|v| v.field_size).collect();

    let steps = match format {
        ResultFormat::Ascii => parse_ascii_steps(&bytes[consumed.min(bytes.len())..], &sizes)?,
        ResultFormat::Binary(e) => parse_binary_steps(&bytes[body_start..], e, &sizes)?,
    };
    Ok(ResultFile { header, steps })
}

/// `name[...] : size permsize dofs : equation` (already trimmed, lowercase).
fn parse_var_line(line: &str) -> Result<VarHeader> {
    let bad = || SaveIoError::Format(format!("cannot read variable line '{line}'"));
    let (name, rest) = match line.find(']') {
        Some(j) => (&line[..=j], &line[j + 1..]),
        None => {
            let j = line.find(':').filter(|&j| j > 0).ok_or_else(bad)?;
            (line[..j].trim_end(), &line[j..])
        }
    };
    let rest = rest.trim_start().strip_prefix(':').ok_or_else(bad)?;
    let (nums, eq) = rest.split_once(':').ok_or_else(bad)?;
    let n: Vec<usize> = nums
        .split_whitespace()
        .map(|t| t.parse().map_err(|_| bad()))
        .collect::<Result<_>>()?;
    let [field_size, perm_size, dofs] = n[..] else {
        return Err(bad());
    };
    Ok(VarHeader {
        name: name.to_string(),
        field_size,
        perm_size,
        dofs,
        equation: eq.trim().to_string(),
    })
}

fn parse_ascii_steps(body: &[u8], sizes: &[usize]) -> Result<Vec<StepRecord>> {
    let text = String::from_utf8_lossy(body);
    let mut it = text.lines().peekable();
    let mut steps: Vec<StepRecord> = Vec::new();
    while let Some(line) = it.next() {
        if line.trim().is_empty() {
            continue;
        }
        let lower = line.trim_start().to_lowercase();
        let Some(rest) = lower.strip_prefix("time:") else {
            return format_err(format!("expected 'Time:', found '{line}'"));
        };
        let toks: Vec<&str> = rest.split_whitespace().collect();
        if toks.len() < 3 {
            return format_err(format!("bad time line '{line}'"));
        }
        let mut step = StepRecord {
            saved_count: toks[0]
                .parse()
                .map_err(|_| SaveIoError::Format(format!("bad time line '{line}'")))?,
            timestep: toks[1]
                .parse()
                .map_err(|_| SaveIoError::Format(format!("bad time line '{line}'")))?,
            time: parse_list_real(toks[2])
                .ok_or_else(|| SaveIoError::Format(format!("bad time line '{line}'")))?,
            vars: Vec::new(),
        };
        let mut prev_n_positive: Option<usize> = None;
        while let Some(next) = it.peek() {
            if next.trim_start().to_lowercase().starts_with("time:") {
                break;
            }
            let name = it.next().expect("peeked").trim().to_string();
            if name.is_empty() {
                continue;
            }
            let idx = step.vars.len();
            let permline = it
                .next()
                .ok_or_else(|| SaveIoError::Format(format!("'{name}': missing Perm line")))?;
            let pl = permline.trim_start();
            let after = pl
                .strip_prefix("Perm:")
                .or_else(|| pl.strip_prefix("perm:"))
                .ok_or_else(|| {
                    SaveIoError::Format(format!("'{name}': expected Perm line, found '{permline}'"))
                })?;
            let after = after.trim();
            let (perm, count) = if after.eq_ignore_ascii_case("null") {
                (
                    PermRecord::Null,
                    *sizes.get(idx).ok_or_else(|| {
                        SaveIoError::Format(format!("more variables than header lists at '{name}'"))
                    })?,
                )
            } else if after.starts_with("use previous") {
                (
                    PermRecord::UsePrevious,
                    prev_n_positive.ok_or_else(|| {
                        SaveIoError::Format(format!("'{name}': use previous with no table"))
                    })?,
                )
            } else {
                let n: Vec<usize> = after
                    .split_whitespace()
                    .map(|t| {
                        t.parse()
                            .map_err(|_| SaveIoError::Format(format!("bad Perm line '{permline}'")))
                    })
                    .collect::<Result<_>>()?;
                let [size, npos] = n[..] else {
                    return format_err(format!("bad Perm line '{permline}'"));
                };
                let mut entries = Vec::with_capacity(npos);
                for _ in 0..npos {
                    let l = it
                        .next()
                        .ok_or_else(|| SaveIoError::Format("truncated permutation".into()))?;
                    let t: Vec<u32> = l
                        .split_whitespace()
                        .map(|t| {
                            t.parse().map_err(|_| {
                                SaveIoError::Format(format!("bad permutation entry '{l}'"))
                            })
                        })
                        .collect::<Result<_>>()?;
                    let [i, k] = t[..] else {
                        return format_err(format!("bad permutation entry '{l}'"));
                    };
                    entries.push((i, k));
                }
                prev_n_positive = Some(npos);
                (PermRecord::Table { size, entries }, npos)
            };
            let mut values = Vec::with_capacity(count);
            for _ in 0..count {
                let l = it
                    .next()
                    .ok_or_else(|| SaveIoError::Format(format!("'{name}': truncated values")))?;
                values
                    .push(parse_list_real(l).ok_or_else(|| {
                        SaveIoError::Format(format!("'{name}': bad value '{l}'"))
                    })?);
            }
            step.vars.push(VarRecord { name, perm, values });
        }
        steps.push(step);
    }
    Ok(steps)
}

fn parse_binary_steps(body: &[u8], e: Endian, sizes: &[usize]) -> Result<Vec<StepRecord>> {
    let mut c = Cursor { b: body, at: 0, e };
    let mut steps: Vec<StepRecord> = Vec::new();
    while !c.eof() {
        let tag = c.string()?;
        if tag != "Time:" {
            return format_err(format!("expected 'Time:', found '{tag}'"));
        }
        let mut step = StepRecord {
            saved_count: c.i32()?,
            timestep: c.i32()?,
            time: c.f64()?,
            vars: Vec::new(),
        };
        let mut prev_n_positive: Option<usize> = None;
        loop {
            if c.eof() || body[c.at..].starts_with(b"Time:\0") {
                break;
            }
            let name = c.string()?;
            let idx = step.vars.len();
            let np = c.i32()?;
            let (perm, count) = if np < 0 {
                c.i64()?; // position of the table; the previous table is in memory
                (
                    PermRecord::UsePrevious,
                    prev_n_positive.ok_or_else(|| {
                        SaveIoError::Format(format!("'{name}': use previous with no table"))
                    })?,
                )
            } else if np == 0 {
                (
                    PermRecord::Null,
                    *sizes.get(idx).ok_or_else(|| {
                        SaveIoError::Format(format!("more variables than header lists at '{name}'"))
                    })?,
                )
            } else {
                let npos = c.i32()?;
                if npos < 0 {
                    return format_err("negative permutation count");
                }
                let mut entries = Vec::with_capacity(npos as usize);
                for _ in 0..npos {
                    let i = c.i32()?;
                    let k = c.i32()?;
                    entries.push((i as u32, k as u32));
                }
                prev_n_positive = Some(npos as usize);
                (
                    PermRecord::Table {
                        size: np as usize,
                        entries,
                    },
                    npos as usize,
                )
            };
            let mut values = Vec::with_capacity(count);
            for _ in 0..count {
                values.push(c.f64()?);
            }
            step.vars.push(VarRecord { name, perm, values });
        }
        steps.push(step);
    }
    Ok(steps)
}

/// Parsed `.pos` index of a binary file.
#[derive(Clone, Debug, PartialEq)]
pub struct PosFile {
    pub endian: Endian,
    /// Scalar variable names in header order.
    pub names: Vec<String>,
    /// Per step: byte offset of `Time:`, then of each stored variable.
    pub offsets: Vec<Vec<u64>>,
}

impl PosFile {
    pub fn from_bytes(bytes: &[u8]) -> Result<PosFile> {
        let endian = match bytes.first() {
            Some(b'L') => Endian::Little,
            Some(b'B') => Endian::Big,
            _ => return format_err("bad .pos endianness byte"),
        };
        let mut c = Cursor {
            b: bytes,
            at: 1,
            e: endian,
        };
        let n = c.i32()?;
        let names = (0..n).map(|_| c.string()).collect::<Result<Vec<_>>>()?;
        let per = (names.len() + 1) * 8;
        let rest = bytes.len() - c.at;
        if !rest.is_multiple_of(per) {
            return format_err(".pos body is not a whole number of steps");
        }
        let mut offsets = Vec::new();
        while !c.eof() {
            offsets.push(
                (0..=names.len())
                    .map(|_| c.i64().map(|x| x as u64))
                    .collect::<Result<Vec<_>>>()?,
            );
        }
        Ok(PosFile {
            endian,
            names,
            offsets,
        })
    }
}

/// Is `name` selected by `Output Variable i` / `Restart Variable i` entries?
/// A listed vector name also selects its numbered components
/// (`velocity` selects `velocity 1`, `velocity 2`, ...). An empty list selects all.
pub fn variable_listed(list: &[String], name: &str) -> bool {
    if list.is_empty() {
        return true;
    }
    list.iter().any(|entry| {
        name.len() >= entry.len()
            && name.is_char_boundary(entry.len())
            && &name[..entry.len()] == entry
            && name[entry.len()..]
                .chars()
                .all(|c| " 0123456789".contains(c))
    })
}

/// Read a `.result` file (ASCII or binary, format version 3).
pub fn read_result(path: &std::path::Path) -> Result<ResultFile> {
    ResultFile::from_bytes(&std::fs::read(path)?)
}
