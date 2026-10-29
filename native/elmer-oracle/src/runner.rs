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
            _ if (c.is_whitespace() || c == '(' || c == ')') && !in_quote => {
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
                        _ => args.push(rest[i].clone()),
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
}

/// Insert capture keywords into every root-level `.sif` in `dir`.
/// Returns the modified file names. Sections already carrying the
/// keywords are left untouched.
pub fn inject_sif(dir: &Path, inj: &SifInject) -> Result<Vec<String>> {
    let mut modified = Vec::new();
    let entries = std::fs::read_dir(dir).map_err(|e| io_with(dir, e))?;
    for entry in entries {
        let entry = entry.map_err(|e| io_with(dir, e))?;
        let path = entry.path();
        if path.extension().is_none_or(|e| e != "sif") {
            continue;
        }
        let text = std::fs::read_to_string(&path).map_err(|e| io_with(&path, e))?;
        let out = inject_text(&text, inj);
        if out != text {
            std::fs::write(&path, &out).map_err(|e| io_with(&path, e))?;
            modified.push(path.file_name().unwrap_or_default().to_string_lossy().to_string());
        }
    }
    Ok(modified)
}

fn inject_text(text: &str, inj: &SifInject) -> String {
    let mut out = String::new();
    let mut section = String::new();
    let mut block: Vec<String> = Vec::new();
    let flush = |out: &mut String, section: &str, block: &[String]| {
        let is_simulation = section.eq_ignore_ascii_case("simulation");
        let is_solver = section.to_lowercase().starts_with("solver");
        let has_output = block.iter().any(|l| l.to_lowercase().contains("output file"));
        let has_save = block.iter().any(|l| l.to_lowercase().contains("linear system save"));
        for line in block {
            out.push_str(line);
            out.push('\n');
        }
        if inj.o2 && is_simulation && !has_output {
            out.push_str("  Output File = \"case.result\"\n");
            out.push_str("  Binary Output = Logical True\n");
        }
        if inj.o2m && is_solver && !has_save {
            out.push_str("  Linear System Save = Logical True\n");
            out.push_str("  Linear System Save Perm = Logical True\n");
        }
    };
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.eq_ignore_ascii_case("end") && !section.is_empty() {
            flush(&mut out, &section, &block);
            block.clear();
            section.clear();
            out.push_str(line);
            out.push('\n');
        } else if section.is_empty() && !trimmed.is_empty() && !trimmed.starts_with('!') && !trimmed.contains('=') && !trimmed.starts_with('$') {
            // Lines outside any section (leading comments, `Solver 1 ::`
            // prefixes, the `RUN` trailer) pass through verbatim; only
            // then does the new section start.
            for pending in block.drain(..) {
                out.push_str(&pending);
                out.push('\n');
            }
            section = trimmed.to_string();
            out.push_str(line);
            out.push('\n');
        } else {
            block.push(line.to_string());
        }
    }
    // Unterminated trailing block: flush without injection.
    for line in &block {
        out.push_str(line);
        out.push('\n');
    }
    out
}

 /// Stage a test by copying its source tree, then run the plan steps.
 pub fn run_plan(
     arm: &Arm,
     spec: &TestSpec,
     plan: &Plan,
     source_dir: &Path,
     run_dir: &Path,
    inject: Option<&SifInject>,
 ) -> Result<RunOutcome> {
     stage_tree(source_dir, run_dir)?;
    if let Some(inj) = inject {
        inject_sif(run_dir, inj)?;
    }
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
        let inj = SifInject { o2: true, o2m: true };
        let once = inject_text(sif, &inj);
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
        assert_eq!(inject_text(&once, &inj), once);
    }

    #[test]
    fn injection_skips_files_without_target_sections() {
        let sif = "Header\n  Mesh DB = \".\" \"m\"\nEnd\n";
        let inj = SifInject { o2: true, o2m: true };
        assert_eq!(inject_text(sif, &inj), sif);
    }

    #[test]
    fn trailing_prefix_and_run_lines_survive() {
        // `Solver 1 ::` prefixes and the `RUN` trailer live outside any
        // section; a bare word like RUN must not eat them (HeatControl2
        // lost its Reference Norm this way and graded vacuous).
        let sif = "Solver 1\n  Equation = \"Heat Equation\"\nEnd\n\nSolver 1 :: Reference Norm = 428.79\nRUN\n";
        let idle = SifInject { o2: false, o2m: false };
        assert_eq!(inject_text(sif, &idle), sif);
        let inj = SifInject { o2: true, o2m: true };
        let once = inject_text(sif, &inj);
        assert!(once.contains("Solver 1 :: Reference Norm = 428.79"));
        assert!(once.contains("\nRUN\n"));
    }
}
