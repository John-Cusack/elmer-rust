//! Direct runner: executes each test's `runtest.cmake` steps outside ctest.
//!
//! The runner stages a test directory (source files copied to a run dir),
//! replays the preparation tools (ElmerGrid, ViewFactors, ...), then runs
//! `ElmerSolver` with the environment `RUN_ELMER_TEST` sets, and reads the
//! solver's own `TEST.PASSED` verdict plus the `CompareToReferenceSolution`
//! comparison lines. Only the declarative subset of `runtest.cmake` is
//! covered; anything with custom CMake logic (macros, conditionals, file
//! comparisons) is reported as [`Plan::custom`] with the offending
//! construct named, never silently skipped.

use crate::{io_with, OracleError, Result};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Paths into one Fortran build arm plus the MPI switch.
#[derive(Debug, Clone)]
pub struct Arm {
    pub name: String,
    pub build_dir: PathBuf,
    pub elmergrid: PathBuf,
    pub elmersolver: PathBuf,
    pub viewfactors: PathBuf,
    pub radiators: PathBuf,
    pub mesh2d: Option<PathBuf>,
    pub mpiexec: Option<PathBuf>,
    pub with_mpi: bool,
    /// `elmersolver` is a candidate binary, not the arm's own: cmake
    /// replays pass it as `ELMERSOLVER_BIN` instead of the recorded one.
    pub solver_override: bool,
}

impl Arm {
    pub fn tool_path(&self, bin: &ToolBin) -> Option<PathBuf> {
        match bin {
            ToolBin::ElmerGrid => Some(self.elmergrid.clone()),
            ToolBin::ViewFactors => Some(self.viewfactors.clone()),
            ToolBin::Radiators => Some(self.radiators.clone()),
            ToolBin::Mesh2d => self.mesh2d.clone(),
        }
    }
}

/// Test registration parsed from the test's `CMakeLists.txt`.
#[derive(Debug, Clone, PartialEq)]
pub struct TestSpec {
    pub name: String,
    /// Task counts from `NPROCS`; empty means serial.
    pub nprocs: Vec<u32>,
    pub labels: Vec<String>,
    /// True when `ADD_ELMERTEST_MODULE` compiles a test-local Fortran
    /// plugin: the test needs a Rust fixture port (D12) and stays out of
    /// the golden pool until the fixture exists.
    pub has_test_module: bool,
    /// Files staged by `CONFIGURE_FILE`/`file(COPY ...)`.
    pub staged_files: Vec<String>,
    /// Sif named by `ELMERSOLVER_STARTINFO`, when present.
    pub startinfo_sif: Option<String>,
}

/// Preparation-tool binaries addressable from `execute_process`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ToolBin {
    ElmerGrid,
    ViewFactors,
    Radiators,
    Mesh2d,
}

/// One replayable step of a test plan.
#[derive(Debug, Clone, PartialEq)]
pub enum Step {
    Tool {
        bin: ToolBin,
        args: Vec<String>,
    },
    Solver {
        /// Explicit sif (`EXECUTE_ELMER_SOLVER`); `None` means the
        /// `ELMERSOLVER_STARTINFO` file (`RUN_ELMER_TEST`).
        sif: Option<String>,
        np: u32,
        extra_lib: Vec<String>,
        /// Only `RUN_ELMER_TEST` steps check `TEST.PASSED`; intermediate
        /// `EXECUTE_*` steps only need a zero exit, mirroring ctest.
        verdict: bool,
    },
}

/// The executable plan for one (test, np) pair. `custom` names the first
/// unsupported construct when the file is not declarative.
#[derive(Debug, Clone)]
pub struct Plan {
    pub steps: Vec<Step>,
    pub custom: Option<String>,
    pub stdout_log: String,
    pub stderr_log: String,
    pub passed_file: String,
}

impl Plan {
    fn unsupported(reason: String) -> Self {
        Plan {
            steps: Vec::new(),
            custom: Some(reason),
            stdout_log: String::new(),
            stderr_log: String::new(),
            passed_file: String::new(),
        }
    }
}

/// Join physical lines into logical lines on parenthesis depth, so
/// multi-line `execute_process(...)` calls parse as one command.
fn logical_lines(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut depth: i32 = 0;
    for raw in text.lines() {
        let line = raw.trim();
        if cur.is_empty() && (line.is_empty() || line.starts_with('#')) {
            continue;
        }
        if !cur.is_empty() {
            cur.push(' ');
        }
        cur.push_str(line);
        depth += line.bytes().filter(|&b| b == b'(').count() as i32;
        depth -= line.bytes().filter(|&b| b == b')').count() as i32;
        if depth <= 0 {
            // Strip a trailing `# comment` outside quotes.
            let stripped = strip_comment(&cur);
            if !stripped.trim().is_empty() {
                out.push(stripped);
            }
            cur.clear();
            depth = 0;
        }
    }
    if !cur.trim().is_empty() {
        out.push(strip_comment(&cur));
    }
    out
}

fn strip_comment(line: &str) -> String {
    let mut in_quote = false;
    let bytes = line.as_bytes();
    for i in 0..bytes.len() {
        match bytes[i] {
            b'"' => in_quote = !in_quote,
            b'#' if !in_quote && (i == 0 || bytes[i - 1].is_ascii_whitespace()) => {
                return line[..i].to_string();
            }
            _ => {}
        }
    }
    line.to_string()
}

/// Split a logical line into tokens honoring double quotes. Parentheses
/// are their own tokens, so `RUN_ELMER_TEST()` and `file(COPY ...)` split
/// into a bare command head plus arguments. Quotes are removed; empty
/// quoted strings are dropped.
fn tokenize(line: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut cur = String::new();
    let mut in_quote = false;
    let flush = |cur: &mut String, tokens: &mut Vec<String>| {
        if !cur.is_empty() {
            tokens.push(std::mem::take(cur));
        }
    };
    for c in line.chars() {
        match c {
            '"' => in_quote = !in_quote,
            // An unquoted `;` separates CMake list elements, so
            // `ElmerGrid 1 2 Step;` passes `Step` and nothing else.
            _ if (c.is_whitespace() || c == '(' || c == ')' || c == ';') && !in_quote => {
                flush(&mut cur, &mut tokens);
            }
            _ => cur.push(c),
        }
    }
    flush(&mut cur, &mut tokens);
    tokens
}

fn tool_var(token: &str) -> Option<ToolBin> {
    match token {
        "${ELMERGRID_BIN}" => Some(ToolBin::ElmerGrid),
        "${VIEWFACTORS_BIN}" => Some(ToolBin::ViewFactors),
        "${RADIATORS_BIN}" => Some(ToolBin::Radiators),
        "${MESH2D_BIN}" => Some(ToolBin::Mesh2d),
        _ => None,
    }
}

/// Build the plan for one (test, np) pair from the test's `runtest.cmake`.
pub fn plan_test(spec: &TestSpec, runtest: &str, np: u32, _arm: &Arm) -> Result<Plan> {
    let mut steps = Vec::new();
    let log_suffix = if np > 1 {
        format!("_{np}")
    } else {
        String::new()
    };
    let mut last_solver_log = ("test-stdout.log".to_string(), "test-stderr.log".to_string());
    if np > 1 {
        last_solver_log = (
            format!("test-stdout{log_suffix}.log"),
            format!("test-stderr{log_suffix}.log"),
        );
    }

    for line in logical_lines(runtest) {
        let tokens = tokenize(&line);
        if tokens.is_empty() {
            continue;
        }
        let head = tokens[0].to_lowercase();
        match head.as_str() {
            "include" => continue,
            "execute_process" => {
                let cmd_pos = tokens
                    .iter()
                    .position(|t| t == "COMMAND")
                    .ok_or_else(|| {
                        OracleError::Unsupported(format!("execute_process without COMMAND: {line}"))
                    })?;
                // Options before COMMAND are not supported.
                if cmd_pos != 1 {
                    return Ok(Plan::unsupported(format!(
                        "execute_process options before COMMAND: {line}"
                    )));
                }
                let mut rest = &tokens[cmd_pos + 1..];
                // COMMAND takes exactly the binary plus its args up to the
                // next CMake keyword.
                let bin_token = rest.first().cloned().unwrap_or_default();
                rest = &rest[1.min(rest.len())..];
                let mut args: Vec<String> = Vec::new();
                let mut stdout_file: Option<String> = None;
                let mut stderr_file: Option<String> = None;
                let mut result_var = false;
                let mut i = 0;
                while i < rest.len() {
                    match rest[i].as_str() {
                        "OUTPUT_FILE" => {
                            i += 1;
                            stdout_file = rest.get(i).cloned();
                        }
                        "ERROR_FILE" => {
                            i += 1;
                            stderr_file = rest.get(i).cloned();
                        }
                        "RESULT_VARIABLE" | "OUTPUT_VARIABLE" | "ERROR_VARIABLE" => {
                            result_var = true;
                            i += 1;
                        }
                        "WORKING_DIRECTORY" => {
                            i += 1;
                        }
                        // `cmake -P` runs in the test's binary dir, so
                        // its current binary dir is the run dir itself.
                        _ => args.push(
                            rest[i]
                                .replace("${MPIEXEC_NTASKS}", &np.to_string())
                                .replace("${CMAKE_CURRENT_BINARY_DIR}", "."),
                        ),
                    }
                    i += 1;
                }
                if bin_token == "${CMAKE_COMMAND}" {
                    return Ok(Plan::unsupported(format!(
                        "CMAKE_COMMAND custom step (file compare or helper): {line}"
                    )));
                }
                if bin_token == "${MPIEXEC}" {
                    return Ok(Plan::unsupported(format!(
                        "bare MPIEXEC outside solver macros: {line}"
                    )));
                }
                let Some(bin) = tool_var(&bin_token) else {
                    return Ok(Plan::unsupported(format!(
                        "unknown tool binary {bin_token}: {line}"
                    )));
                };
                if result_var || stdout_file.is_some() || stderr_file.is_some() {
                    // Log redirection is recorded; RESULT_VARIABLE checks
                    // are custom verdict logic.
                    if result_var {
                        return Ok(Plan::unsupported(format!(
                            "RESULT_VARIABLE verdict logic: {line}"
                        )));
                    }
                }
                let _ = (stdout_file, stderr_file);
                steps.push(Step::Tool { bin, args });
            }
            "run_elmer_test" => {
                let mut extra_lib = Vec::new();
                if let Some(pos) = tokens.iter().position(|t| t == "ELMER_LIB") {
                    extra_lib = tokens[pos + 1..].to_vec();
                }
                steps.push(Step::Solver {
                    sif: None,
                    np,
                    extra_lib,
                    verdict: true,
                });
            }
            "execute_elmer_solver" => {
                let sif = tokens.get(1).cloned().filter(|s| !s.starts_with("ELMER_LIB"));
                let mut extra_lib = Vec::new();
                if let Some(pos) = tokens.iter().position(|t| t == "ELMER_LIB") {
                    extra_lib = tokens[pos + 1..].to_vec();
                }
                steps.push(Step::Solver {
                    sif,
                    np: 1,
                    extra_lib,
                    verdict: false,
                });
            }
            "execute_elmer_solver_mpi" => {
                let sif = tokens.get(1).cloned().filter(|s| !s.starts_with("ELMER_LIB"));
                steps.push(Step::Solver {
                    sif,
                    np,
                    extra_lib: Vec::new(),
                    verdict: false,
                });
            }
            _ => {
                return Ok(Plan::unsupported(format!(
                    "non-declarative construct: {line}"
                )));
            }
        }
    }
    if steps.iter().all(|s| !matches!(s, Step::Solver { .. })) {
        return Ok(Plan::unsupported(String::from(
            "no solver step: file never runs ElmerSolver",
        )));
    }
    let _ = spec;
    Ok(Plan {
        steps,
        custom: None,
        stdout_log: last_solver_log.0,
        stderr_log: last_solver_log.1,
        passed_file: if np > 1 {
            format!("TEST.PASSED{log_suffix}")
        } else {
            String::from("TEST.PASSED")
        },
    })
}

/// Parse `CMakeLists.txt` registration for one test directory.
pub fn parse_cmakelists(name: &str, text: &str) -> TestSpec {
    let mut spec = TestSpec {
        name: name.to_string(),
        nprocs: Vec::new(),
        labels: Vec::new(),
        has_test_module: text.contains("ADD_ELMERTEST_MODULE"),
        staged_files: Vec::new(),
        startinfo_sif: None,
    };
    for line in logical_lines(text) {
        let tokens = tokenize(&line);
        if tokens.is_empty() {
            continue;
        }
        match tokens[0].as_str() {
            "ADD_ELMER_TEST" => {
                if let Some(pos) = tokens.iter().position(|t| t == "NPROCS") {
                    for t in &tokens[pos + 1..] {
                        if t == "LABELS" {
                            break;
                        }
                        if let Ok(n) = t.parse::<u32>() {
                            spec.nprocs.push(n);
                        }
                    }
                }
                if let Some(pos) = tokens.iter().position(|t| t == "LABELS") {
                    spec.labels = tokens[pos + 1..].to_vec();
                }
            }
            "CONFIGURE_FILE" => {
                if let Some(src) = tokens.get(1) {
                    spec.staged_files.push(src.clone());
                }
            }
            _ => {}
        }
        // file(COPY a b DESTINATION ...) stages mesh and startinfo files.
        if tokens[0] == "file" && tokens.iter().any(|t| t == "COPY") {
            for t in tokens.iter().skip(1) {
                if t == "COPY" || t == "DESTINATION" {
                    continue;
                }
                if t.contains('.') {
                    spec.staged_files.push(t.clone());
                }
            }
        }
    }
    spec
}

/// A test's `runtest.cmake` (`runTest.cmake` in the ElmerIce suite).
pub fn read_runtest(dir: &Path) -> Result<String> {
    let path = ["runtest.cmake", "runTest.cmake"]
        .iter()
        .map(|n| dir.join(n))
        .find(|p| p.exists())
        .unwrap_or_else(|| dir.join("runtest.cmake"));
    std::fs::read_to_string(&path).map_err(|e| io_with(&path, e))
}

/// Read the sif name from an `ELMERSOLVER_STARTINFO` file (first token).
pub fn read_startinfo(dir: &Path) -> Result<Option<String>> {
    let path = dir.join("ELMERSOLVER_STARTINFO");
    if !path.exists() {
        return Ok(None);
    }
    let text = std::fs::read_to_string(&path).map_err(|e| io_with(&path, e))?;
    Ok(text.split_whitespace().next().map(String::from))
}

/// Outcome of running one plan.
#[derive(Debug)]
pub struct RunOutcome {
    pub passed: bool,
    /// Whether the last executed solver/tool process exited zero.
    /// Intermediate `EXECUTE_*` steps use this instead of the verdict.
    pub exit_ok: bool,
    pub failed_step: Option<String>,
    pub compare_lines: Vec<String>,
    pub cpu_time: Option<f64>,
    pub real_time: Option<f64>,
}

/// Golden-capture augmentation for staged `.sif` files. The keywords use
/// the exact forms the pristine parser accepts (verified against
/// `ModelDescription.F90` and the quick suite): quoted `Output File`
/// needs a `file`-typed keyword, and logicals need the `Logical` word.
/// `Linear System Save Prefix` is deliberately absent: it is not a
/// declared keyword and the parser rejects it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SifInject {
    /// Append `Output File` + `Binary Output` to the Simulation section.
    pub o2: bool,
    /// Append `Linear System Save` + `Save Perm` to every Solver section.
    pub o2m: bool,
    /// Raise (never lower) `Max Output Level` in the Simulation section,
    /// so iteration counts reach the log as diagnostics.
    pub max_output: Option<u32>,
}

impl SifInject {
    pub fn is_active(&self) -> bool {
        self.o2 || self.o2m || self.max_output.is_some()
    }
}

/// Result file the O2 injection names for one sif. The startinfo sif keeps
/// `case.result` (the name every existing golden uses); any other sif of a
/// multi-sif test gets its own `<stem>.result`, so successive solver runs
/// against the same mesh cannot overwrite each other's golden.
pub fn result_name_for(sif: &str, startinfo: Option<&str>) -> String {
    if startinfo.is_none_or(|s| s == sif) {
        return String::from("case.result");
    }
    let stem = sif.strip_suffix(".sif").unwrap_or(sif);
    format!("{stem}.result")
}

/// Insert capture keywords into every root-level `.sif` in `dir`.
/// Returns the modified file names. Sections already carrying the
/// keywords are left untouched.
pub fn inject_sif(dir: &Path, inj: &SifInject) -> Result<Vec<String>> {
    inject_sif_with(dir, inj, &[])
}

/// [`inject_sif`], plus `Reference Norm` overrides `(solver number,
/// value)` in the startinfo sif only (the bitwise grader validation,
/// §4.3 item 2): each named Solver section gains the value with 17
/// significant digits and `Reference Norm Tolerance = 0`.
pub fn inject_sif_with(dir: &Path, inj: &SifInject, refnorms: &[(u32, f64)]) -> Result<Vec<String>> {
    let mut modified = Vec::new();
    let startinfo = read_startinfo(dir)?;
    let entries = std::fs::read_dir(dir).map_err(|e| io_with(dir, e))?;
    for entry in entries {
        let entry = entry.map_err(|e| io_with(dir, e))?;
        let path = entry.path();
        if path.extension().is_none_or(|e| e != "sif") {
            continue;
        }
        // Sifs are bytes, not UTF-8 (several quick tests carry Latin-1
        // comments): decode each byte as one char so the round trip is
        // lossless; the injection itself only adds ASCII.
        let bytes = std::fs::read(&path).map_err(|e| io_with(&path, e))?;
        let text: String = bytes.iter().map(|&b| char::from(b)).collect();
        let name = path.file_name().unwrap_or_default().to_string_lossy().to_string();
        let result_name = result_name_for(&name, startinfo.as_deref());
        let is_startinfo = startinfo.as_deref() == Some(name.as_str());
        let out = inject_text_with(&text, inj, &result_name, if is_startinfo { refnorms } else { &[] });
        if out != text {
            let out: Vec<u8> = out.chars().map(|c| c as u32 as u8).collect();
            std::fs::write(&path, &out).map_err(|e| io_with(&path, e))?;
            modified.push(path.file_name().unwrap_or_default().to_string_lossy().to_string());
        }
    }
    Ok(modified)
}

/// Raise every parseable `Max Output Level` line to `level`; append one
/// when the section has none. Returns the block lines to emit.
fn raise_max_output(block: &[String], level: u32) -> Vec<String> {
    let mut lines = block.to_vec();
    let mut found = false;
    for l in lines.iter_mut() {
        if keyword_of(l).as_deref() != Some("max output level") {
            continue;
        }
        let Some((head, tail)) = strip_sif_comment(l).split_once('=') else {
            continue;
        };
        let head = head.to_string();
        let value = tail.trim();
        let value = value
            .strip_prefix("Integer")
            .or_else(|| value.strip_prefix("integer"))
            .unwrap_or(value)
            .trim();
        if let Ok(v) = value.parse::<u32>() {
            found = true;
            if v < level {
                *l = format!("{head}= {level}");
            }
        }
    }
    if !found {
        lines.push(format!("  Max Output Level = {level}"));
    }
    lines
}

/// Whether a (comment-stripped) sif line opens a multi-line table that
/// its own `End` closes: the value ends in a bare type word, optionally
/// with interpolation modifiers (`Real`, `= Variable Time Real`,
/// `Real cubic`), as `SectionContents` reads it. `Real MATC "..."` and
/// `Integer 20` carry their value and open nothing.
pub(crate) fn opens_table(line: &str) -> bool {
    let value = line.rsplit_once('=').map_or(line, |(_, v)| v).to_lowercase();
    let mut tokens: Vec<&str> = value.split_whitespace().collect();
    while tokens
        .last()
        .is_some_and(|t| matches!(*t, "cubic" | "monotone" | "harmonic" | "linear"))
    {
        tokens.pop();
    }
    tokens
        .last()
        .is_some_and(|t| matches!(*t, "real" | "integer" | "logical" | "string" | "file"))
}

/// A sif line with its `!`/`#` comment removed (quotes respected).
pub(crate) fn strip_sif_comment(line: &str) -> &str {
    let mut in_quote = false;
    for (i, c) in line.char_indices() {
        match c {
            '"' => in_quote = !in_quote,
            '!' | '#' if !in_quote => return &line[..i],
            _ => {}
        }
    }
    line
}

/// Lowercased keyword of a `key = value` line, whitespace collapsed;
/// `None` for comments, table rows and lines without `=`.
pub(crate) fn keyword_of(line: &str) -> Option<String> {
    let (key, _) = strip_sif_comment(line).split_once('=')?;
    let key = key.split_whitespace().collect::<Vec<_>>().join(" ").to_lowercase();
    (!key.is_empty()).then_some(key)
}

/// Insert the capture keywords into the Simulation and Solver sections of
/// one sif, following the top level of `LoadInputFile`: directives
/// (`Check Keywords`, `Echo`, `Include`, `RUN`), `Section :: key = value`
/// one-liners and `/* */` blocks pass through; any other line opens a
/// section closed by its `End` (table `End`s do not count).
#[cfg(test)]
fn inject_text(text: &str, inj: &SifInject, result_name: &str) -> String {
    inject_text_with(text, inj, result_name, &[])
}

fn inject_text_with(text: &str, inj: &SifInject, result_name: &str, refnorms: &[(u32, f64)]) -> String {
    let mut out = String::new();
    let mut section = String::new();
    let mut block: Vec<String> = Vec::new();
    // `Simulation :: Output File = ...` one-liners count as present.
    let oneliner_output = text.lines().any(|l| {
        let code = strip_sif_comment(l).to_lowercase();
        code.split_once("::").is_some_and(|(sec, rest)| {
            sec.trim() == "simulation" && keyword_of(rest).as_deref() == Some("output file")
        })
    });
    let flush = |out: &mut String, section: &str, block: &[String]| {
        let is_simulation = section.starts_with("simulation");
        let is_solver = section.starts_with("solver");
        let keys: Vec<String> = block.iter().filter_map(|l| keyword_of(l)).collect();
        let has_output = oneliner_output || keys.iter().any(|k| k == "output file");
        let has_save = keys.iter().any(|k| k == "linear system save");
        let lines = match inj.max_output {
            Some(level) if is_simulation => raise_max_output(block, level),
            _ => block.to_vec(),
        };
        for line in &lines {
            out.push_str(line);
            out.push('\n');
        }
        if inj.o2 && is_simulation && !has_output {
            out.push_str(&format!("  Output File = \"{result_name}\"\n"));
            out.push_str("  Binary Output = Logical True\n");
        }
        let number: Option<u32> = section
            .strip_prefix("solver")
            .and_then(|n| n.trim().parse().ok());
        for (_, value) in refnorms
            .iter()
            .filter(|(n, _)| Some(*n) == number && !oneliner_refnorm(text, *n))
        {
            out.push_str(&format!(
                "  Reference Norm = Real {}\n",
                crate::norm::format_grader_norm(*value)
            ));
            out.push_str("  Reference Norm Tolerance = Real 0\n");
        }
        if inj.o2m && is_solver && !has_save {
            out.push_str("  Linear System Save = Logical True\n");
            out.push_str("  Linear System Save Perm = Logical True\n");
        }
    };
    let mut in_table = false;
    let mut in_comment_block = false;
    for line in text.lines() {
        let code = strip_sif_comment(line).trim();
        let lower = code.to_lowercase();
        if in_comment_block {
            in_comment_block = lower != "*/";
            out.push_str(line);
            out.push('\n');
        } else if in_table {
            // Elmer ends a table at the first row that starts with `end`.
            in_table = !lower.starts_with("end");
            block.push(line.to_string());
        } else if !section.is_empty() {
            if lower == "end" {
                flush(&mut out, &section, &block);
                block.clear();
                section.clear();
                out.push_str(line);
                out.push('\n');
            } else {
                in_table = opens_table(code);
                block.push(line.to_string());
            }
        } else {
            let directive = lower.is_empty()
                || lower.starts_with('$')
                || lower.contains("::")
                || lower.starts_with("include")
                || lower.starts_with("echo ")
                || lower.starts_with("check ")
                || lower == "run"
                || lower == "stop";
            if lower == "/*" {
                in_comment_block = true;
            } else if !directive {
                section = lower;
            }
            out.push_str(line);
            out.push('\n');
        }
    }
    // Unterminated trailing block: flush without injection.
    for line in &block {
        out.push_str(line);
        out.push('\n');
    }
    // A `Solver N :: Reference Norm...` one-liner overrides the section, so
    // an override for such a solver goes right after its last one, ahead of
    // any `RUN`.
    for (n, value) in refnorms.iter().filter(|(n, _)| oneliner_refnorm(text, *n)) {
        let lines: Vec<&str> = out.lines().collect();
        let Some(at) = lines.iter().rposition(|l| is_refnorm_oneliner(l, *n)) else {
            continue;
        };
        let mut rebuilt: Vec<String> = lines.iter().map(|l| l.to_string()).collect();
        rebuilt.insert(
            at + 1,
            format!(
                "Solver {n} :: Reference Norm = Real {}\nSolver {n} :: Reference Norm Tolerance = Real 0",
                crate::norm::format_grader_norm(*value)
            ),
        );
        out = rebuilt.join("\n");
        out.push('\n');
    }
    out
}

/// `Solver N :: Reference Norm...` (any keyword starting so).
fn is_refnorm_oneliner(line: &str, n: u32) -> bool {
    let code = strip_sif_comment(line).to_lowercase();
    let Some((sec, rest)) = code.split_once("::") else {
        return false;
    };
    let sec_n: Option<u32> = sec.trim().strip_prefix("solver").and_then(|x| x.trim().parse().ok());
    sec_n == Some(n) && keyword_of(rest).is_some_and(|k| k.starts_with("reference norm"))
}

fn oneliner_refnorm(text: &str, n: u32) -> bool {
    text.lines().any(|l| is_refnorm_oneliner(l, n))
}

 /// Stage a test by copying its source tree, then run the plan steps.
 pub fn run_plan(
     arm: &Arm,
     spec: &TestSpec,
     plan: &Plan,
     source_dir: &Path,
     run_dir: &Path,
    staging: &Staging,
 ) -> Result<RunOutcome> {
    stage_test(arm, source_dir, run_dir, staging)?;
    for step in &plan.steps {
        match step {
            Step::Tool { bin, args } => {
                let exe = arm.tool_path(bin).ok_or_else(|| {
                    OracleError::SolverFailed(format!("tool {bin:?} not built in arm {}", arm.name))
                })?;
                let status = Command::new(&exe)
                    .args(args)
                    .current_dir(run_dir)
                    .envs(tool_env(arm))
                    .status()
                    .map_err(|e| OracleError::SolverFailed(format!("spawn {}: {e}", exe.display())))?;
                if !status.success() {
                    return Ok(RunOutcome {
                        passed: false,
                        exit_ok: false,
                        failed_step: Some(format!("tool {bin:?} exited {status}")),
                        compare_lines: Vec::new(),
                        cpu_time: None,
                        real_time: None,
                    });
                }
            }
            Step::Solver {
                sif,
                np,
                extra_lib,
                verdict,
            } => {
                let sif_name = match sif {
                    Some(s) => s.clone(),
                    None => read_startinfo(run_dir)?.ok_or_else(|| {
                        OracleError::SolverFailed(String::from(
                            "RUN_ELMER_TEST without ELMERSOLVER_STARTINFO",
                        ))
                    })?,
                };
                let log_base = if sif.is_some() {
                    sif_name.clone()
                } else {
                    String::from("test")
                };
                let outcome = run_solver(arm, spec, run_dir, &sif_name, *np, extra_lib, &log_base)?;
                if *verdict {
                    if !outcome.passed {
                        return Ok(outcome);
                    }
                    // Multi-solver tests keep going only on success,
                    // mirroring ctest: the last solver step owns the
                    // plan logs.
                    let _ = outcome;
                } else if !outcome.exit_ok {
                    return Ok(RunOutcome {
                        passed: false,
                        failed_step: Some(format!(
                            "intermediate solver on {sif_name} exited nonzero"
                        )),
                        compare_lines: outcome.compare_lines,
                        cpu_time: outcome.cpu_time,
                        real_time: outcome.real_time,
                        exit_ok: false,
                    });
                }
            }
        }
    }
    // Re-read the final verdict files for the graded outcome.
    let stdout = std::fs::read_to_string(run_dir.join(&plan.stdout_log)).unwrap_or_default();
    Ok(RunOutcome {
        passed: std::fs::read_to_string(run_dir.join(&plan.passed_file))
            .map(|s| s.trim() == "1")
            .unwrap_or(false),
        exit_ok: true,
        failed_step: None,
        compare_lines: compare_lines(&stdout),
        cpu_time: parse_total_time(&stdout, 1),
        real_time: parse_total_time(&stdout, 2),
    })
}

/// What the oracle adds to a staged test beyond its sources.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Staging {
    pub inject: Option<SifInject>,
    /// Directories to create in the run dir before any step runs: the
    /// homes of solver-internal meshes (`<mesh>-isomesh`, split meshes)
    /// that exist only in memory, so the injected `Output File`, which
    /// Elmer writes as `<mesh name>/<file>`, has somewhere to land.
    pub mkdirs: Vec<String>,
    /// `Reference Norm` overrides for the startinfo sif (bitwise grader
    /// validation).
    pub refnorms: Vec<(u32, f64)>,
}

/// Copy the test sources into `run_dir`, add the arm's own build of any
/// test-local plugin, then apply the capture injection.
fn stage_test(arm: &Arm, source_dir: &Path, run_dir: &Path, staging: &Staging) -> Result<()> {
    stage_tree(source_dir, run_dir)?;
    if let Some(dir) = source_dir.file_name() {
        // ElmerIce tests live in `elmerice/Tests`, the rest in `fem/tests`.
        let suite = match source_dir.parent().and_then(|p| p.file_name()) {
            Some(p) if p == "Tests" => arm.build_dir.join("elmerice").join("Tests"),
            _ => arm.build_dir.join("fem").join("tests"),
        };
        stage_arm_modules(&suite.join(dir), run_dir)?;
    }
    // Configure-time `FILE(MAKE_DIRECTORY ...)` dirs: git keeps no empty
    // directories, so the source tree cannot carry them.
    let cmake = std::fs::read_to_string(source_dir.join("CMakeLists.txt")).unwrap_or_default();
    for d in make_directories(&cmake) {
        let path = run_dir.join(d);
        std::fs::create_dir_all(&path).map_err(|e| io_with(&path, e))?;
    }
    if staging.inject.is_some() || !staging.refnorms.is_empty() {
        inject_sif_with(run_dir, &staging.inject.unwrap_or_default(), &staging.refnorms)?;
    }
    for d in &staging.mkdirs {
        let path = run_dir.join(d);
        std::fs::create_dir_all(&path).map_err(|e| io_with(&path, e))?;
    }
    Ok(())
}

/// Output directories a failed run could not open: the parent of every
/// `Cannot open file '<dir>/<name>.result...'` the solver logs report.
/// Only result files the injection asked for qualify, so a test whose own
/// inputs are missing is never papered over.
pub fn missing_output_dirs(run_dir: &Path) -> Vec<String> {
    let mut dirs = Vec::new();
    let Ok(entries) = std::fs::read_dir(run_dir) else {
        return dirs;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().is_none_or(|e| e != "log") {
            continue;
        }
        let Ok(bytes) = std::fs::read(&path) else {
            continue;
        };
        let text = String::from_utf8_lossy(&bytes);
        for piece in text.split("Cannot open file '").skip(1) {
            let Some((file, _)) = piece.split_once('\'') else {
                continue;
            };
            let Some((dir, name)) = file.rsplit_once('/') else {
                continue;
            };
            let is_result = name.ends_with(".result")
                || name
                    .rsplit_once(".result.")
                    .is_some_and(|(_, n)| n.bytes().all(|b| b.is_ascii_digit()));
            if is_result
                && !dir.is_empty()
                && !dir.starts_with('/')
                && !dir.contains("..")
                && !dirs.iter().any(|d| d == dir)
            {
                dirs.push(dir.to_string());
            }
        }
    }
    dirs.sort();
    dirs
}

/// Relative directories a test's `CMakeLists.txt` creates with
/// `FILE(MAKE_DIRECTORY ...)` in its binary dir.
pub fn make_directories(cmake: &str) -> Vec<String> {
    let mut out = Vec::new();
    for line in logical_lines(cmake) {
        let tokens = tokenize(&line);
        if !tokens.first().is_some_and(|t| t.eq_ignore_ascii_case("file")) {
            continue;
        }
        if !tokens.get(1).is_some_and(|t| t == "MAKE_DIRECTORY") {
            continue;
        }
        for t in &tokens[2..] {
            let rel = t
                .strip_prefix("${CMAKE_CURRENT_BINARY_DIR}/")
                .unwrap_or(t)
                .trim_end_matches('/');
            if !rel.is_empty() && !rel.contains("${") && !rel.starts_with('/') && !rel.contains("..") {
                out.push(rel.to_string());
            }
        }
    }
    out
}

/// `ADD_ELMERTEST_MODULE` builds test-local plugins into the arm's binary
/// test directory, never the source tree; each arm uses its own build of
/// them (the FMA arm's plugin carries the FMA flags too).
fn stage_arm_modules(build_test_dir: &Path, run_dir: &Path) -> Result<()> {
    let Ok(entries) = std::fs::read_dir(build_test_dir) else {
        return Ok(());
    };
    for entry in entries {
        let entry = entry.map_err(|e| io_with(build_test_dir, e))?;
        let src = entry.path();
        if src.is_file() && src.extension().is_some_and(|e| e == "so") {
            let dst = run_dir.join(entry.file_name());
            std::fs::copy(&src, &dst).map_err(|e| io_with(&src, e))?;
        }
    }
    Ok(())
}

/// Replay a registered test's own `runtest.cmake` through `cmake -P` in a
/// staged run dir, with the arguments ctest would pass (container paths
/// remapped to the host). This is the faithful path for tests whose
/// `runtest.cmake` carries custom logic (`IF`, macros, `file(...)`) that
/// [`plan_test`] refuses to interpret: the verdict is CMake's own, a zero
/// exit plus a passing `TEST.PASSED`.
pub fn run_replay(
    arm: &Arm,
    reg: &crate::registry::RegisteredTest,
    source_dir: &Path,
    run_dir: &Path,
    staging: &Staging,
) -> Result<RunOutcome> {
    stage_test(arm, source_dir, run_dir, staging)?;
    let recorded = reg.define("BINARY_DIR").ok_or_else(|| {
        OracleError::Parse(format!("{}: registration has no BINARY_DIR", reg.name))
    })?;
    let mut args: Vec<String> = Vec::new();
    let mut iter = reg.command.iter().skip(1);
    while let Some(a) = iter.next() {
        if a == "-P" {
            // Run the staged copy, never the source tree's (ElmerIce
            // spells it `runTest.cmake`).
            let script = iter
                .next()
                .and_then(|p| p.rsplit('/').next())
                .unwrap_or("runtest.cmake")
                .to_string();
            args.push(a.clone());
            args.push(run_dir.join(script).display().to_string());
            continue;
        }
        if arm.solver_override && a.starts_with("-DELMERSOLVER_BIN=") {
            args.push(format!("-DELMERSOLVER_BIN={}", arm.elmersolver.display()));
            continue;
        }
        args.push(crate::registry::remap_arg(a, recorded, &arm.build_dir));
    }
    let passed_name = if reg.np > 1 {
        format!("TEST.PASSED_{}", reg.np)
    } else {
        String::from("TEST.PASSED")
    };
    let _ = std::fs::remove_file(run_dir.join(&passed_name));
    let build = arm.build_dir.display().to_string();
    let output = Command::new("cmake")
        .args(&args)
        .current_dir(run_dir)
        .env("LD_LIBRARY_PATH", arm_lib_dirs(&build))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .map_err(|e| OracleError::SolverFailed(format!("spawn cmake replay: {e}")))?;
    std::fs::write(run_dir.join("replay-stdout.log"), &output.stdout)
        .map_err(|e| io_with(&run_dir.join("replay-stdout.log"), e))?;
    std::fs::write(run_dir.join("replay-stderr.log"), &output.stderr)
        .map_err(|e| io_with(&run_dir.join("replay-stderr.log"), e))?;
    // An MPI-enabled build names the logs `_<ntasks>` even for one task.
    let stdout = std::fs::read_to_string(run_dir.join(format!("test-stdout_{}.log", reg.np)))
        .or_else(|_| std::fs::read_to_string(run_dir.join("test-stdout.log")))
        .unwrap_or_default();
    let verdict = std::fs::read_to_string(run_dir.join(&passed_name))
        .map(|s| s.trim() == "1")
        .unwrap_or(false);
    let passed = verdict && output.status.success();
    Ok(RunOutcome {
        passed,
        exit_ok: output.status.success(),
        failed_step: if passed {
            None
        } else {
            Some(format!("cmake replay exited {}", output.status))
        },
        compare_lines: compare_lines(&stdout),
        cpu_time: parse_total_time(&stdout, 1),
        real_time: parse_total_time(&stdout, 2),
    })
}

/// Library directories of one arm build, mapped onto the host. The arm
/// binaries carry container-baked RUNPATH (`/work/...`), so the host
/// loader needs these explicitly; system BLAS/LAPACK/MPI resolve from
/// the host's matching Ubuntu 24.04 stack (recorded per arm).
fn arm_lib_dirs(build: &str) -> String {
    let existing = std::env::var("LD_LIBRARY_PATH").unwrap_or_default();
    let dirs = [
        format!("{build}/fem/src"),
        format!("{build}/matc/src"),
        format!("{build}/fhutiter/src"),
        format!("{build}/mathlibs/src/arpack"),
        format!("{build}/mathlibs/src/parpack"),
        // ElmerIce solver and USF libraries link ElmerIceUtils.
        format!("{build}/elmerice/Utils"),
    ];
    if existing.is_empty() {
        dirs.join(":")
    } else {
        format!("{}:{existing}", dirs.join(":"))
    }
}

fn tool_env(arm: &Arm) -> Vec<(String, String)> {
    let build = arm.build_dir.display().to_string();
    vec![
        (String::from("ELMER_HOME"), format!("{build}/fem/src")),
        (String::from("ELMER_LIB"), format!("{build}/fem/src/modules")),
        (
            String::from("PATH"),
            format!(
                "{build}/meshgen2d/src:{build}/fem/src:{}",
                std::env::var("PATH").unwrap_or_default()
            ),
        ),
        (String::from("LD_LIBRARY_PATH"), arm_lib_dirs(&build)),
    ]
}

fn run_solver(
    arm: &Arm,
    _spec: &TestSpec,
    run_dir: &Path,
    sif: &str,
    np: u32,
    extra_lib: &[String],
    log_base: &str,
) -> Result<RunOutcome> {
    let build = arm.build_dir.display().to_string();
    let mut elmer_lib = format!("{build}/fem/src/modules");
    for extra in extra_lib {
        elmer_lib.push(':');
        elmer_lib.push_str(&format!("{build}/fem/src/modules/{extra}"));
    }
    let suffix = if np > 1 { format!("_{np}") } else { String::new() };
    let stdout_name = format!("{log_base}-stdout{suffix}.log");
    let stderr_name = format!("{log_base}-stderr{suffix}.log");
    // Remove a stale verdict so a no-run cannot pass.
    let passed_name = if np > 1 {
        format!("TEST.PASSED_{np}")
    } else {
        String::from("TEST.PASSED")
    };
    let _ = std::fs::remove_file(run_dir.join(&passed_name));

    let mut cmd = if np > 1 {
        let mpiexec = arm.mpiexec.clone().ok_or_else(|| {
            OracleError::SolverFailed(format!("arm {} has no mpiexec", arm.name))
        })?;
        let mut c = Command::new(mpiexec);
        c.arg("-n").arg(np.to_string());
        c.arg(&arm.elmersolver);
        c
    } else {
        Command::new(&arm.elmersolver)
    };
    // Solver binaries take the sif as argv; ELMER_LIB/ELMER_HOME select
    // the module tree, exactly as RUN_ELMER_TEST sets them.
    let output = cmd
        .arg(sif)
        .current_dir(run_dir)
        .env("ELMER_HOME", format!("{build}/fem/src"))
        .env("ELMER_LIB", &elmer_lib)
        .env(
            "PATH",
            format!(
                "{build}/meshgen2d/src:{build}/fem/src:{}",
                std::env::var("PATH").unwrap_or_default()
            ),
        )
        .env("LD_LIBRARY_PATH", arm_lib_dirs(&build))
        .env("OMP_STACKSIZE", std::env::var("OMP_STACKSIZE").unwrap_or_else(|_| String::from("8M")))
        .env("OMP_WAIT_POLICY", std::env::var("OMP_WAIT_POLICY").unwrap_or_else(|_| String::from("PASSIVE")))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .map_err(|e| OracleError::SolverFailed(format!("spawn solver: {e}")))?;
    std::fs::write(run_dir.join(&stdout_name), &output.stdout)
        .map_err(|e| io_with(&run_dir.join(&stdout_name), e))?;
    std::fs::write(run_dir.join(&stderr_name), &output.stderr)
        .map_err(|e| io_with(&run_dir.join(&stderr_name), e))?;
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let passed = std::fs::read_to_string(run_dir.join(&passed_name))
        .map(|s| s.trim() == "1")
        .unwrap_or(false);
    Ok(RunOutcome {
        passed,
        exit_ok: output.status.success(),
        failed_step: if passed {
            None
        } else {
            Some(format!("solver on {sif} wrote no passing verdict"))
        },
        compare_lines: compare_lines(&stdout),
        cpu_time: parse_total_time(&stdout, 1),
        real_time: parse_total_time(&stdout, 2),
    })
}

/// Lines the solver wrote about its own reference comparison.
pub fn compare_lines(stdout: &str) -> Vec<String> {
    stdout
        .lines()
        .filter(|l| l.contains("CompareToReferenceSolution"))
        .map(String::from)
        .collect()
}

/// Parse `SOLVER TOTAL TIME(CPU,REAL): <cpu> <real>`; `which` selects 1/2.
pub fn parse_total_time(stdout: &str, which: usize) -> Option<f64> {
    stdout.lines().rev().find_map(|l| {
        let (_, tail) = l.split_once("SOLVER TOTAL TIME(CPU,REAL):")?;
        let mut nums = tail.split_whitespace();
        let first = nums.next()?;
        let second = nums.next()?;
        let pick = if which == 1 { first } else { second };
        pick.parse::<f64>().ok()
    })
}

fn stage_tree(source: &Path, dest: &Path) -> Result<()> {
    std::fs::create_dir_all(dest).map_err(|e| io_with(dest, e))?;
    let entries = std::fs::read_dir(source).map_err(|e| io_with(source, e))?;
    for entry in entries {
        let entry = entry.map_err(|e| io_with(source, e))?;
        let name = entry.file_name();
        // Never copy build leftovers into a fresh run dir.
        if name == "TEST.PASSED"
            || name.to_string_lossy().starts_with("TEST.PASSED_")
            || name.to_string_lossy().starts_with("test-stdout")
            || name.to_string_lossy().starts_with("test-stderr")
        {
            continue;
        }
        let src = entry.path();
        let dst = dest.join(&name);
        if src.is_dir() {
            stage_tree(&src, &dst)?;
        } else {
            std::fs::copy(&src, &dst).map_err(|e| io_with(&src, e))?;
        }
    }
    Ok(())
}

/// Discover test directories: every child of `tests_dir` holding a
/// `runtest.cmake` is one test, mirroring ctest registration.
pub fn discover_tests(tests_dir: &Path) -> Result<Vec<String>> {
    let mut names = Vec::new();
    let entries = std::fs::read_dir(tests_dir).map_err(|e| io_with(tests_dir, e))?;
    for entry in entries {
        let entry = entry.map_err(|e| io_with(tests_dir, e))?;
        if entry.path().join("runtest.cmake").exists() {
            names.push(entry.file_name().to_string_lossy().to_string());
        }
    }
    names.sort();
    Ok(names)
}

/// Load known `${VAR} → path` substitutions for plan display.
pub fn known_bins(arm: &Arm) -> HashMap<String, String> {
    let mut m = HashMap::new();
    m.insert(String::from("${ELMERGRID_BIN}"), arm.elmergrid.display().to_string());
    m.insert(String::from("${ELMERSOLVER_BIN}"), arm.elmersolver.display().to_string());
    m.insert(String::from("${VIEWFACTORS_BIN}"), arm.viewfactors.display().to_string());
    m.insert(String::from("${RADIATORS_BIN}"), arm.radiators.display().to_string());
    if let Some(m2) = &arm.mesh2d {
        m.insert(String::from("${MESH2D_BIN}"), m2.display().to_string());
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_arm() -> Arm {
        Arm {
            name: "test".to_string(),
            build_dir: PathBuf::from("/build"),
            elmergrid: PathBuf::from("/build/bin/ElmerGrid"),
            elmersolver: PathBuf::from("/build/fem/src/ElmerSolver"),
            viewfactors: PathBuf::from("/build/fem/src/ViewFactors"),
            radiators: PathBuf::from("/build/fem/src/Radiators"),
            mesh2d: None,
            mpiexec: None,
            with_mpi: false,
            solver_override: false,
        }
    }

    fn spec() -> TestSpec {
        TestSpec {
            name: "poisson".to_string(),
            nprocs: Vec::new(),
            labels: vec![String::from("quick")],
            has_test_module: false,
            staged_files: Vec::new(),
            startinfo_sif: Some(String::from("case.sif")),
        }
    }

    #[test]
    fn poisson_plan_is_declarative() {
        let runtest = "include(test_macros)\nexecute_process(COMMAND ${ELMERGRID_BIN} 1 2 hexwire.grd)\nRUN_ELMER_TEST()\n";
        let plan = plan_test(&spec(), runtest, 1, &test_arm()).unwrap();
        assert!(plan.custom.is_none());
        assert_eq!(plan.steps.len(), 2);
        assert!(matches!(plan.steps[0], Step::Tool { .. }));
        assert_eq!(
            plan.steps[1],
            Step::Solver {
                sif: None,
                np: 1,
                extra_lib: Vec::new(),
                verdict: true
            }
        );
        assert_eq!(plan.passed_file, "TEST.PASSED");
    }

    #[test]
    fn tool_args_split_lists_and_expand_ntasks() {
        let runtest = "include(test_macros)\nexecute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step;)\nexecute_process(COMMAND ${ELMERGRID_BIN} 1 2 a.grd -metis ${MPIEXEC_NTASKS} 3)\nRUN_ELMER_TEST()\n";
        let plan = plan_test(&spec(), runtest, 4, &test_arm()).unwrap();
        let args: Vec<_> = plan
            .steps
            .iter()
            .filter_map(|s| match s {
                Step::Tool { args, .. } => Some(args.clone()),
                _ => None,
            })
            .collect();
        assert_eq!(args[0], vec!["1", "2", "Step"]);
        assert_eq!(args[1], vec!["1", "2", "a.grd", "-metis", "4", "3"]);
    }

    #[test]
    fn max_output_raises_never_lowers() {
        let sif = "Simulation\n  Max Output Level = 3\nEnd\n";
        let inj = SifInject { o2: false, o2m: false, max_output: Some(5) };
        let once = inject_text(sif, &inj, "case.result");
        assert!(once.contains("Max Output Level = 5"));
        assert_eq!(inject_text(&once, &inj, "case.result"), once);
        let high = "Simulation\n  Max Output Level = 30\nEnd\n";
        assert_eq!(inject_text(high, &inj, "case.result"), high);
        let none = inject_text("Simulation\nEnd\n", &inj, "case.result");
        assert!(none.contains("Max Output Level = 5"));
    }

    #[test]
    fn non_startinfo_sifs_get_their_own_result_file() {
        assert_eq!(result_name_for("case.sif", Some("case.sif")), "case.result");
        assert_eq!(result_name_for("case.sif", None), "case.result");
        assert_eq!(result_name_for("quad.sif", Some("tet.sif")), "quad.result");
    }

    #[test]
    fn missing_output_dirs_reads_only_result_targets() {
        let dir = std::env::temp_dir().join(format!("elmer-oracle-mod-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("test-stderr.log"),
            "Fortran runtime error: Cannot open file 'square-tri-isomesh/case.result': No such file\n\
             Fortran runtime error: Cannot open file 'part/case.result.3': No such file\n\
             Fortran runtime error: Cannot open file 'mesh/mesh.nodes': No such file\n",
        )
        .unwrap();
        assert_eq!(missing_output_dirs(&dir), vec!["part", "square-tri-isomesh"]);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn latin1_sif_survives_injection_byte_for_byte() {
        let dir = std::env::temp_dir().join(format!("elmer-oracle-l1-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut sif = b"! Bj\xf6rn's case\nSimulation\n  Simulation Type = Steady\nEnd\n".to_vec();
        std::fs::write(dir.join("case.sif"), &sif).unwrap();
        let inj = SifInject { o2: true, o2m: false, max_output: None };
        inject_sif(&dir, &inj).unwrap();
        let out = std::fs::read(dir.join("case.sif")).unwrap();
        let tail = b"  Output File = \"case.result\"\n  Binary Output = Logical True\nEnd\n";
        sif.truncate(sif.len() - 4);
        sif.extend_from_slice(tail);
        assert_eq!(out, sif);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn table_end_does_not_close_the_section() {
        let sif = "Simulation\n  Density = Variable Coordinate 1\n    Real\n      0 1\n      1 2\n    End\n  Max Output Level = 5\nEnd\n";
        let inj = SifInject { o2: true, o2m: false, max_output: None };
        let out = inject_text(sif, &inj, "case.result");
        let table_end = out.find("    End\n").unwrap();
        let output = out.find("Output File").unwrap();
        assert!(output > table_end, "injected inside the table:\n{out}");
        assert!(out.ends_with("  Binary Output = Logical True\nEnd\n"));
        assert!(opens_table("Real cubic"));
        assert!(opens_table("Name = Variable Time Real"));
        assert!(!opens_table("Real MATC \"1+tx\""));
        assert!(!opens_table("Extruded Mesh Levels = Integer 20"));
    }

    #[test]
    fn directives_and_comments_do_not_hide_sections() {
        let sif = "Check Keywords \"Warn\"\n\nHeader :: Mesh DB \".\" \"square\"\n\nSimulation\n! Output File = \"old.dat\"\n! Max Output Level = 3\n  Max Output Level = 4\nEnd\n/*\nSolver 9\nEnd\n*/\nSolver 1 :: Reference Norm = 1.0\nRUN\n";
        let inj = SifInject { o2: true, o2m: true, max_output: Some(5) };
        let out = inject_text(sif, &inj, "case.result");
        assert!(out.contains("  Output File = \"case.result\"\n"), "{out}");
        assert!(out.contains("! Max Output Level = 3\n"), "comment rewritten: {out}");
        assert!(out.contains("  Max Output Level = 5\n"), "{out}");
        // The commented-out Solver 9 inside /* */ gains nothing.
        assert!(!out.contains("Linear System Save"), "{out}");
        assert!(out.ends_with("Solver 1 :: Reference Norm = 1.0\nRUN\n"));
    }

    #[test]
    fn elmerice_runtest_spelling_is_found() {
        let dir = std::env::temp_dir().join(format!("elmer-oracle-rt-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("runTest.cmake"), "RUN_ELMERICE_TEST()\n").unwrap();
        assert_eq!(read_runtest(&dir).unwrap(), "RUN_ELMERICE_TEST()\n");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn make_directory_targets_are_relative_dirs() {
        let cmake = "IF(HAVE_UMFPACK)\n  FILE(MAKE_DIRECTORY results)\n  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/kulma)\n  FILE(COPY a DESTINATION \"${CMAKE_CURRENT_BINARY_DIR}/\")\nENDIF()\n";
        assert_eq!(make_directories(cmake), vec!["results", "kulma"]);
    }

    #[test]
    fn reference_norm_override_lands_in_the_named_solver() {
        let sif = "Solver 1\n  Equation = A\nEnd\nSolver 2\n  Equation = B\n  Reference Norm = 1.5\nEnd\n";
        let idle = SifInject::default();
        let out = inject_text_with(sif, &idle, "case.result", &[(2, 0.1)]);
        let s2 = out.split("Solver 2").nth(1).unwrap();
        assert!(s2.contains("Reference Norm = Real 1.00000000000000006e-1\n"), "{out}");
        assert!(s2.contains("Reference Norm Tolerance = Real 0\n"));
        assert!(!out.split("Solver 2").next().unwrap().contains("Tolerance"));
    }

    #[test]
    fn reference_norm_override_follows_the_last_oneliner() {
        let sif = "Solver 2\n  Variable = T\nEnd\nSolver 2 :: Reference Norm = 1.5\nSolver 2 :: Reference Norm Tolerance = 1e-5\nRUN\n";
        let out = inject_text_with(sif, &SifInject::default(), "case.result", &[(2, 0.1)]);
        assert!(
            out.ends_with("Solver 2 :: Reference Norm Tolerance = 1e-5\nSolver 2 :: Reference Norm = Real 1.00000000000000006e-1\nSolver 2 :: Reference Norm Tolerance = Real 0\nRUN\n"),
            "{out}"
        );
        // Nothing lands inside the section.
        assert!(!out.split("End").next().unwrap().contains("Real 0"));
    }

    #[test]
    fn mpi_plan_names_suffixed_files() {
        let runtest = "include(test_macros)\nRUN_ELMER_TEST()\n";
        let plan = plan_test(&spec(), runtest, 4, &test_arm()).unwrap();
        assert_eq!(plan.passed_file, "TEST.PASSED_4");
        assert_eq!(plan.stdout_log, "test-stdout_4.log");
    }

    #[test]
    fn explicit_solver_step_uses_named_sif() {
        let runtest = "include(test_macros)\nEXECUTE_ELMER_SOLVER(numrays.sif)\nRUN_ELMER_TEST()\n";
        let plan = plan_test(&spec(), runtest, 1, &test_arm()).unwrap();
        assert_eq!(plan.steps.len(), 2);
        assert_eq!(
            plan.steps[0],
            Step::Solver {
                sif: Some(String::from("numrays.sif")),
                np: 1,
                extra_lib: Vec::new(),
                verdict: false
            }
        );
    }

    #[test]
    fn custom_logic_is_named_not_skipped() {
        let runtest = "include(test_macros)\nMACRO(RUN_CORNER SIF)\nRUN_TOOL(${ELMERGRID_BIN} elmergrid 1 2 box.grd)\nENDMACRO(RUN_CORNER)\n";
        let plan = plan_test(&spec(), runtest, 1, &test_arm()).unwrap();
        assert!(plan.custom.is_some());
        assert!(plan.steps.is_empty());
    }

    #[test]
    fn cmake_command_steps_are_custom() {
        let runtest = "include(test_macros)\nexecute_process(COMMAND ${CMAKE_COMMAND} -E compare_files a b)\nRUN_ELMER_TEST()\n";
        let plan = plan_test(&spec(), runtest, 1, &test_arm()).unwrap();
        assert!(plan.custom.as_ref().unwrap().contains("CMAKE_COMMAND"));
    }

    #[test]
    fn cmakelists_parses_nprocs_labels_and_module() {
        let text = "INCLUDE(test_macros)\nADD_ELMER_TEST(foo NPROCS 2 4 LABELS quick parallel)\nADD_ELMERTEST_MODULE(foo mymod mymod.F90)\n";
        let spec = parse_cmakelists("foo", text);
        assert_eq!(spec.nprocs, vec![2, 4]);
        assert_eq!(spec.labels, vec!["quick", "parallel"]);
        assert!(spec.has_test_module);
    }

    #[test]
    fn total_time_parses_cpu_and_real() {
        let out = "noise\nSOLVER TOTAL TIME(CPU,REAL):   12.340    56.780\ntrailer\n";
        assert_eq!(parse_total_time(out, 1), Some(12.34));
        assert_eq!(parse_total_time(out, 2), Some(56.78));
    }

    #[test]
    fn compare_lines_selects_solver_verdicts() {
        let out = "foo\nCompareToReferenceSolution: Norm 1.0 vs 1.0\nbar\n";
        assert_eq!(compare_lines(out).len(), 1);
    }

    #[test]
    fn injection_targets_sections_and_is_idempotent() {
        let sif = "Header\n  Mesh DB = \".\" \"m\"\nEnd\n\nSimulation\n  Simulation Type = \"Steady State\"\nEnd\n\nSolver 1\n  Equation = \"Heat Equation\"\nEnd\n";
        let inj = SifInject { o2: true, o2m: true, max_output: None };
        let once = inject_text(sif, &inj, "case.result");
        assert!(once.contains("Output File = \"case.result\""));
        assert!(once.contains("Binary Output = Logical True"));
        assert!(once.contains("Linear System Save = Logical True"));
        assert!(once.contains("Linear System Save Perm = Logical True"));
        // Simulation section must not gain solver keywords and vice versa.
        let sim = once.split("Solver 1").next().unwrap();
        assert!(!sim.contains("Linear System Save"));
        let solver = once.split("Solver 1").nth(1).unwrap();
        assert!(!solver.contains("Output File"));
        // Second pass changes nothing.
        assert_eq!(inject_text(&once, &inj, "case.result"), once);
    }

    #[test]
    fn injection_skips_files_without_target_sections() {
        let sif = "Header\n  Mesh DB = \".\" \"m\"\nEnd\n";
        let inj = SifInject { o2: true, o2m: true, max_output: None };
        assert_eq!(inject_text(sif, &inj, "case.result"), sif);
    }

    #[test]
    fn trailing_prefix_and_run_lines_survive() {
        // `Solver 1 ::` prefixes and the `RUN` trailer live outside any
        // section; a bare word like RUN must not eat them (HeatControl2
        // lost its Reference Norm this way and graded vacuous).
        let sif = "Solver 1\n  Equation = \"Heat Equation\"\nEnd\n\nSolver 1 :: Reference Norm = 428.79\nRUN\n";
        let idle = SifInject { o2: false, o2m: false, max_output: None };
        assert_eq!(inject_text(sif, &idle, "case.result"), sif);
        let inj = SifInject { o2: true, o2m: true, max_output: None };
        let once = inject_text(sif, &inj, "case.result");
        assert!(once.contains("Solver 1 :: Reference Norm = 428.79"));
        assert!(once.contains("\nRUN\n"));
    }
}
