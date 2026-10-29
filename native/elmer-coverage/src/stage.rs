//! Staging one quick test: a scratch copy of its working directory, the
//! mesh generated with the pristine ElmerGrid, and the sequence of solver
//! runs its `runtest.cmake` performs.

use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use crate::discover::QuickTest;

/// One step of a test's `runtest.cmake`, in order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Step {
    /// `execute_process(COMMAND ${ELMERGRID_BIN} ...)`.
    ElmerGrid(Vec<String>),
    /// `RUN_ELMER_TEST()` (the `ELMERSOLVER_STARTINFO` case) or
    /// `EXECUTE_ELMER_SOLVER(x.sif)`.
    Solve(Option<String>),
    /// Any other tool the test runs first; the native path has no
    /// equivalent (`ViewFactors`, `Radiators`, `Mesh2D`, ...).
    Other(String),
}

/// A CMake command call: lowercased name and its arguments, quotes removed.
fn calls(text: &str) -> Vec<(String, Vec<String>)> {
    let src: String = text
        .lines()
        .map(|l| match l.find('#') {
            Some(p) if l[..p].matches('"').count() % 2 == 0 => &l[..p],
            _ => l,
        })
        .collect::<Vec<_>>()
        .join("\n");
    let b = src.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if !(b[i].is_ascii_alphabetic() || b[i] == b'_') {
            i += 1;
            continue;
        }
        let start = i;
        while i < b.len() && (b[i].is_ascii_alphanumeric() || b[i] == b'_') {
            i += 1;
        }
        let name = src[start..i].to_ascii_lowercase();
        let mut j = i;
        while j < b.len() && (b[j] == b' ' || b[j] == b'\t') {
            j += 1;
        }
        if j >= b.len() || b[j] != b'(' {
            continue;
        }
        let (mut depth, mut in_q, mut k) = (0usize, false, j);
        let mut arg = String::new();
        let mut args = Vec::new();
        let mut quoted = false;
        while k < b.len() {
            let c = b[k] as char;
            k += 1;
            if c == '"' {
                in_q = !in_q;
                quoted = true;
                continue;
            }
            if in_q {
                arg.push(c);
                continue;
            }
            match c {
                '(' => {
                    depth += 1;
                    if depth > 1 {
                        arg.push(c);
                    }
                }
                ')' => {
                    depth -= 1;
                    if depth == 0 {
                        break;
                    }
                    arg.push(c);
                }
                c if c.is_whitespace() => {
                    if !arg.is_empty() || quoted {
                        args.push(std::mem::take(&mut arg));
                        quoted = false;
                    }
                }
                c => arg.push(c),
            }
        }
        if !arg.is_empty() || quoted {
            args.push(arg);
        }
        out.push((name, args));
        i = k;
    }
    out
}

/// Replaces the variables a quick test's mesh commands use.
fn substitute(arg: &str, nprocs: u32, source: &Path) -> String {
    arg.replace("${MPIEXEC_NTASKS}", &nprocs.to_string())
        .replace("${TEST_SOURCE}", &source.to_string_lossy())
        .replace("${CMAKE_CURRENT_BINARY_DIR}", ".")
}

/// The steps of `runtest.cmake` in order.
pub fn steps(runtest: &str, nprocs: u32, source: &Path) -> Vec<Step> {
    let mut steps = Vec::new();
    for (name, args) in calls(runtest) {
        match name.as_str() {
            "execute_process" => {
                let words: Vec<String> = args.iter().map(|a| substitute(a, nprocs, source)).collect();
                let Some(p) = words.iter().position(|w| w.eq_ignore_ascii_case("COMMAND")) else { continue };
                let cmd: Vec<String> = words[p + 1..]
                    .iter()
                    .take_while(|w| !is_process_option(w))
                    .cloned()
                    .collect();
                match cmd.first().map(String::as_str) {
                    // An unquoted `;` separates CMake list items (`Step;`).
                    Some("${ELMERGRID_BIN}") => steps.push(Step::ElmerGrid(
                        cmd[1..].iter().flat_map(|a| a.split(';')).filter(|a| !a.is_empty()).map(str::to_string).collect(),
                    )),
                    Some("${VIEWFACTORS_BIN}") => steps.push(Step::Other("ViewFactors".into())),
                    Some("${RADIATORS_BIN}") => steps.push(Step::Other("Radiators".into())),
                    Some("${MESH2D_BIN}") => steps.push(Step::Other("Mesh2D".into())),
                    Some(other) => steps.push(Step::Other(other.trim_matches(&['$', '{', '}'][..]).to_string())),
                    None => {}
                }
            }
            "run_elmer_test" => steps.push(Step::Solve(None)),
            "execute_elmer_solver" => steps.push(Step::Solve(args.first().cloned())),
            _ => {}
        }
    }
    steps
}

fn is_process_option(w: &str) -> bool {
    matches!(
        w.to_ascii_uppercase().as_str(),
        "WORKING_DIRECTORY"
            | "OUTPUT_QUIET"
            | "ERROR_QUIET"
            | "OUTPUT_VARIABLE"
            | "ERROR_VARIABLE"
            | "RESULT_VARIABLE"
            | "OUTPUT_FILE"
            | "ERROR_FILE"
            | "INPUT_FILE"
            | "TIMEOUT"
            | "COMMAND_ECHO"
            | "COMMAND_ERROR_IS_FATAL"
            | "ENCODING"
    )
}

/// Module names the test compiles from its own Fortran sources
/// (`ADD_ELMERTEST_MODULE(<test> <module> <file>)` in its `CMakeLists.txt`).
pub fn local_modules(cmake_lists: &str) -> Vec<String> {
    calls(cmake_lists)
        .into_iter()
        .filter(|(n, _)| n == "add_elmertest_module")
        .filter_map(|(_, a)| a.get(1).cloned())
        .collect()
}

/// What `runtest.cmake` and `CMakeLists.txt` tell about a test.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Recipe {
    pub steps: Vec<Step>,
    pub local_modules: Vec<String>,
}

/// Reads a test's recipe from its source directory.
pub fn recipe(test: &QuickTest) -> Recipe {
    let read = |f: &str| fs::read_to_string(test.source_dir.join(f)).unwrap_or_default();
    Recipe {
        steps: steps(&read("runtest.cmake"), test.nprocs, &test.source_dir),
        local_modules: local_modules(&read("CMakeLists.txt")),
    }
}

/// Files of the build-tree test directory that are CMake machinery.
fn is_cmake_noise(name: &str) -> bool {
    matches!(name, "CMakeFiles" | "Makefile" | "cmake_install.cmake" | "CTestTestfile.cmake")
}

fn copy_tree(src: &Path, dst: &Path, skip_noise: bool, overwrite: bool) -> io::Result<()> {
    fs::create_dir_all(dst)?;
    for e in fs::read_dir(src)? {
        let e = e?;
        let name = e.file_name();
        if skip_noise && is_cmake_noise(&name.to_string_lossy()) {
            continue;
        }
        let to = dst.join(&name);
        let ty = e.file_type()?;
        if ty.is_dir() {
            copy_tree(&e.path(), &to, false, overwrite)?;
        } else if ty.is_file() && (overwrite || !to.exists()) {
            fs::copy(e.path(), &to)?;
        } else if ty.is_symlink() && !to.exists() {
            if let Ok(real) = fs::canonicalize(e.path()) {
                if real.is_file() {
                    fs::copy(real, &to)?;
                }
            }
        }
    }
    Ok(())
}

/// Runs `cmd` to completion or kills it at `limit`.
pub fn run_limited(mut cmd: Command, limit: Duration) -> io::Result<(Option<i32>, bool)> {
    let mut child = cmd.stdin(Stdio::null()).stdout(Stdio::null()).stderr(Stdio::null()).spawn()?;
    let start = Instant::now();
    loop {
        if let Some(st) = child.try_wait()? {
            return Ok((st.code(), false));
        }
        if start.elapsed() > limit {
            let _ = child.kill();
            let _ = child.wait();
            return Ok((None, true));
        }
        std::thread::sleep(Duration::from_millis(20));
    }
}

/// Copies the test into `scratch/<name>` (configured build files first,
/// then whatever the source directory has that the build copy lacks) and
/// runs its ElmerGrid steps with the pristine binary, serially at nice 10.
pub fn stage(test: &QuickTest, steps: &[Step], scratch: &Path, elmergrid: &Path) -> Result<PathBuf, String> {
    let dir = scratch.join(&test.name);
    let _ = fs::remove_dir_all(&dir);
    copy_tree(&test.work_dir, &dir, true, true).map_err(|e| format!("copy {}: {e}", test.work_dir.display()))?;
    if test.source_dir.is_dir() {
        copy_tree(&test.source_dir, &dir, false, false).map_err(|e| format!("copy {}: {e}", test.source_dir.display()))?;
    }
    for step in steps {
        if let Step::ElmerGrid(args) = step {
            let mut cmd = Command::new("nice");
            cmd.args(["-n", "10"]).arg(elmergrid).args(args).current_dir(&dir);
            match run_limited(cmd, Duration::from_secs(120)) {
                Ok((Some(0), _)) => {}
                Ok((code, timed_out)) => {
                    return Err(format!("ElmerGrid {} failed ({})", args.join(" "), if timed_out { "timeout".into() } else { format!("{code:?}") }))
                }
                Err(e) => return Err(format!("ElmerGrid: {e}")),
            }
        }
    }
    Ok(dir)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn elmergrid_and_run() {
        let text = "include(test_macros)\n# comment\nexecute_process(COMMAND ${ELMERGRID_BIN} 14 2 OnePipe.msh -autoclean )\nRUN_ELMER_TEST()\n";
        assert_eq!(
            steps(text, 1, Path::new("/s")),
            vec![Step::ElmerGrid(vec!["14".into(), "2".into(), "OnePipe.msh".into(), "-autoclean".into()]), Step::Solve(None)]
        );
    }

    #[test]
    fn ntasks_and_options() {
        let text = "execute_process(COMMAND ${ELMERGRID_BIN} 1 2 sq -partdual ${MPIEXEC_NTASKS}\n  OUTPUT_QUIET)\nEXECUTE_ELMER_SOLVER(init.sif)\n";
        assert_eq!(
            steps(text, 2, Path::new("/s")),
            vec![
                Step::ElmerGrid(vec!["1".into(), "2".into(), "sq".into(), "-partdual".into(), "2".into()]),
                Step::Solve(Some("init.sif".into())),
            ]
        );
    }

    #[test]
    fn semicolon_is_a_list_separator() {
        let text = "execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step;)\n";
        assert_eq!(steps(text, 1, Path::new("/s")), vec![Step::ElmerGrid(vec!["1".into(), "2".into(), "Step".into()])]);
    }

    #[test]
    fn other_tools() {
        let text = "execute_process(COMMAND ${VIEWFACTORS_BIN} a)\nexecute_process(COMMAND ${MESH2D_BIN} k.mif .)\n";
        assert_eq!(steps(text, 1, Path::new("/s")), vec![Step::Other("ViewFactors".into()), Step::Other("Mesh2D".into())]);
    }

    #[test]
    fn local_plugin_modules() {
        let text = "ADD_ELMERTEST_MODULE(1dtests Poisson Poisson.f90)\nADD_ELMERTEST_MODULE(1dtests Source Source.f90)\n";
        assert_eq!(local_modules(text), ["Poisson", "Source"]);
    }
}
