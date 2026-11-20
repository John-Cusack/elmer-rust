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
    /// `file(COPY <files> DESTINATION <dir>)`: files of the test copied into
    /// a directory of the run (`reload` puts its restart file in the mesh
    /// directory). Paths are relative to the run directory.
    Copy(Vec<String>, String),
    /// `file(MAKE_DIRECTORY <dir>)`, relative to the run directory.
    MakeDir(String),
    /// `file(REMOVE <files>)`, relative to the run directory.
    Remove(Vec<String>),
    /// `${MESH2D_BIN} <args>` in a serial run: the pristine `Mesh2D`, an
    /// external tool (GPL, not ported) as ElmerGrid is. Adaptive remeshing
    /// runs it again from the solver, through the child's `PATH`.
    Mesh2D(Vec<String>),
    /// `${VIEWFACTORS_BIN}` in a serial run: the heat solver computes the
    /// factors itself, but from the mesh files as they are, not from the
    /// `f20.12` rewrite RadiationFactors makes when it runs `ViewFactors`
    /// (the runner tells the child with `ELMER_VIEWFACTORS_PRESTEP`).
    ViewFactors,
    /// Any other tool the test runs first; the native path has no
    /// equivalent (`Mesh2D`, ...). `ViewFactors` and `Radiators` need no
    /// step in a serial run: the heat solver computes their factors itself.
    Other(String),
}

/// The prestep a partitioned run with `${VIEWFACTORS_BIN}` reports.
const PARTITIONED_VIEWFACTORS: &str = "ViewFactors (partitioned radiation: RadiationParallelMeshDistribute)";

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
                    // `elmer-physics-heat` computes the view and radiator
                    // factors in process (`elmer-viewfactors`), serially
                    // only. A partitioned radiation run also needs
                    // `RadiationParallelMeshDistribute` (MeshBasics.F90:1411):
                    // every radiating boundary element and its nodes imported
                    // into each partition, the Gebhart couplings to those
                    // nodes in the parallel matrix, and the Gebhart factors
                    // solved on one rank and broadcast. Not ported.
                    Some("${VIEWFACTORS_BIN}") if nprocs > 1 => steps.push(Step::Other(PARTITIONED_VIEWFACTORS.into())),
                    Some("${RADIATORS_BIN}") if nprocs > 1 => steps.push(Step::Other("Radiators".into())),
                    Some("${VIEWFACTORS_BIN}") => steps.push(Step::ViewFactors),
                    Some("${RADIATORS_BIN}") => {}
                    Some("${MESH2D_BIN}") if nprocs > 1 => steps.push(Step::Other("Mesh2D".into())),
                    Some("${MESH2D_BIN}") => steps.push(Step::Mesh2D(cmd[1..].to_vec())),
                    Some(other) => steps.push(Step::Other(other.trim_matches(&['$', '{', '}'][..]).to_string())),
                    None => {}
                }
            }
            "file" if args.first().is_some_and(|a| a.eq_ignore_ascii_case("COPY")) => {
                let words: Vec<String> = args[1..].iter().map(|a| substitute(a, nprocs, source)).collect();
                let Some(p) = words.iter().position(|w| w.eq_ignore_ascii_case("DESTINATION")) else { continue };
                if let Some(dest) = words.get(p + 1) {
                    steps.push(Step::Copy(words[..p].to_vec(), dest.clone()));
                }
            }
            "file" if args.first().is_some_and(|a| a.eq_ignore_ascii_case("MAKE_DIRECTORY")) => {
                for d in &args[1..] {
                    steps.push(Step::MakeDir(substitute(d, nprocs, source)));
                }
            }
            "file" if args.first().is_some_and(|a| a.eq_ignore_ascii_case("REMOVE")) => {
                steps.push(Step::Remove(args[1..].iter().map(|a| substitute(a, nprocs, source)).collect()));
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
/// The file that marks a scratch directory as elmer-coverage's own: only
/// in such a directory are per-test directories deleted.
pub const SCRATCH_MARKER: &str = ".elmer-coverage-scratch";

/// Takes `scratch` as the run's scratch directory: one that does not exist
/// is created, an empty one adopted, and either gets [`SCRATCH_MARKER`]. A
/// non-empty directory without the marker was not created by elmer-coverage
/// and is refused, so nothing in it is ever deleted.
pub fn prepare_scratch(scratch: &Path) -> Result<(), String> {
    let marker = scratch.join(SCRATCH_MARKER);
    if marker.is_file() {
        return Ok(());
    }
    if scratch.exists() {
        let empty = fs::read_dir(scratch).map_err(|e| format!("scratch {}: {e}", scratch.display()))?.next().is_none();
        if !empty {
            return Err(format!(
                "scratch {} exists, is not empty and was not created by elmer-coverage (no {SCRATCH_MARKER}); pass a new or empty --scratch",
                scratch.display()
            ));
        }
    } else {
        fs::create_dir_all(scratch).map_err(|e| format!("scratch {}: {e}", scratch.display()))?;
    }
    fs::write(&marker, "elmer-coverage scratch: per-test directories here are deleted and recreated\n")
        .map_err(|e| format!("scratch {}: {e}", marker.display()))
}

/// Whether `scratch` carries [`SCRATCH_MARKER`] (see [`prepare_scratch`]).
pub fn owns_scratch(scratch: &Path) -> bool {
    scratch.join(SCRATCH_MARKER).is_file()
}

/// Stages `test` in `scratch/<name>`, replacing an earlier copy. Refuses
/// a scratch that [`prepare_scratch`] did not take.
pub fn stage(test: &QuickTest, steps: &[Step], scratch: &Path, elmergrid: &Path, mesh2d: &Path) -> Result<PathBuf, String> {
    if !owns_scratch(scratch) {
        return Err(format!("scratch {} is not an elmer-coverage scratch (no {SCRATCH_MARKER})", scratch.display()));
    }
    let dir = scratch.join(&test.name);
    let _ = fs::remove_dir_all(&dir);
    copy_tree(&test.work_dir, &dir, true, true).map_err(|e| format!("copy {}: {e}", test.work_dir.display()))?;
    if test.source_dir.is_dir() {
        copy_tree(&test.source_dir, &dir, false, false).map_err(|e| format!("copy {}: {e}", test.source_dir.display()))?;
    }
    for step in steps {
        if let Step::Copy(files, dest) = step {
            copy_into(&dir, files, dest)?;
        }
        if let Step::MakeDir(d) = step {
            let to = inside(&dir, d).ok_or_else(|| format!("file(MAKE_DIRECTORY): {d} outside the run directory"))?;
            fs::create_dir_all(&to).map_err(|e| format!("file(MAKE_DIRECTORY) {}: {e}", to.display()))?;
        }
        if let Step::Remove(files) = step {
            for f in files {
                let p = inside(&dir, f).ok_or_else(|| format!("file(REMOVE): {f} outside the run directory"))?;
                let _ = fs::remove_file(&p);
            }
        }
        if let Step::Mesh2D(args) = step {
            let mut cmd = Command::new("nice");
            cmd.args(["-n", "10"]).arg(mesh2d).args(args).current_dir(&dir);
            match run_limited(cmd, Duration::from_secs(120)) {
                Ok((Some(0), _)) => {}
                Ok((code, timed_out)) => {
                    return Err(format!("Mesh2D {} failed ({})", args.join(" "), if timed_out { "timeout".into() } else { format!("{code:?}") }))
                }
                Err(e) => return Err(format!("Mesh2D: {e}")),
            }
        }
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

/// A relative path that stays inside the run directory, or `None`.
fn inside(dir: &Path, rel: &str) -> Option<PathBuf> {
    let p = Path::new(rel);
    let ok = p.components().all(|c| matches!(c, std::path::Component::Normal(_) | std::path::Component::CurDir));
    ok.then(|| dir.join(p))
}

/// [`Step::Copy`] in the staged directory `dir`: only paths inside it.
fn copy_into(dir: &Path, files: &[String], dest: &str) -> Result<(), String> {
    let to = inside(dir, dest).ok_or_else(|| format!("file(COPY): destination {dest} outside the run directory"))?;
    fs::create_dir_all(&to).map_err(|e| format!("file(COPY) {}: {e}", to.display()))?;
    for f in files {
        let from = inside(dir, f).ok_or_else(|| format!("file(COPY): {f} outside the run directory"))?;
        let name = from.file_name().ok_or_else(|| format!("file(COPY): {f} names no file"))?;
        fs::copy(&from, to.join(name)).map_err(|e| format!("file(COPY) {}: {e}", from.display()))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("elmer-coverage-test-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        d
    }

    /// A non-empty directory elmer-coverage did not create is refused and
    /// left alone, staging into it included; a new or empty one is taken
    /// and marked.
    #[test]
    fn foreign_scratch_is_never_cleaned() {
        let foreign = tmp("foreign");
        fs::create_dir_all(foreign.join("SomeTest")).unwrap();
        fs::write(foreign.join("SomeTest/keep.txt"), "old lane data").unwrap();
        assert!(prepare_scratch(&foreign).is_err());
        assert!(!owns_scratch(&foreign));
        let test = QuickTest {
            name: "SomeTest".into(),
            work_dir: foreign.join("nowhere"),
            source_dir: foreign.join("nowhere"),
            nprocs: 1,
            labels: vec![],
        };
        assert!(stage(&test, &[], &foreign, Path::new("ElmerGrid"), Path::new("Mesh2D")).is_err());
        assert_eq!(fs::read_to_string(foreign.join("SomeTest/keep.txt")).unwrap(), "old lane data");

        let fresh = tmp("fresh");
        prepare_scratch(&fresh).unwrap();
        assert!(owns_scratch(&fresh));
        prepare_scratch(&fresh).unwrap();
        let empty = tmp("empty");
        fs::create_dir_all(&empty).unwrap();
        prepare_scratch(&empty).unwrap();
        assert!(owns_scratch(&empty));
        for d in [foreign, fresh, empty] {
            fs::remove_dir_all(d).unwrap();
        }
    }

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

    /// `reload`'s recipe: the restart file goes into the mesh directory
    /// ElmerGrid makes, and nothing outside the run directory is touched.
    #[test]
    fn file_copy_into_the_run() {
        let text = "execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step)\nfile(COPY flow.dat flow.dat.pos DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Step/)\nRUN_ELMER_TEST()\n";
        let st = steps(text, 1, Path::new("/s"));
        assert_eq!(st[1], Step::Copy(vec!["flow.dat".into(), "flow.dat.pos".into()], "./Step/".into()));
        let d = tmp("copy");
        fs::create_dir_all(d.join("Step")).unwrap();
        fs::write(d.join("flow.dat"), "x").unwrap();
        fs::write(d.join("flow.dat.pos"), "y").unwrap();
        copy_into(&d, &["flow.dat".into(), "flow.dat.pos".into()], "./Step/").unwrap();
        assert_eq!(fs::read_to_string(d.join("Step/flow.dat.pos")).unwrap(), "y");
        assert!(copy_into(&d, &["flow.dat".into()], "../elsewhere").is_err());
        assert!(copy_into(&d, &["/etc/hostname".into()], "Step").is_err());
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn make_directory_and_remove() {
        let text = "file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/kulma)\nfile(REMOVE mesh.nodes mesh.header)\n";
        assert_eq!(
            steps(text, 1, Path::new("/s")),
            vec![Step::MakeDir("./kulma".into()), Step::Remove(vec!["mesh.nodes".into(), "mesh.header".into()])]
        );
    }

    #[test]
    fn other_tools() {
        let text = "execute_process(COMMAND ${VIEWFACTORS_BIN} a)\nexecute_process(COMMAND ${MESH2D_BIN} k.mif .)\n";
        assert_eq!(steps(text, 1, Path::new("/s")), vec![Step::ViewFactors, Step::Mesh2D(vec!["k.mif".into(), ".".into()])]);
        assert_eq!(steps(text, 2, Path::new("/s")), vec![Step::Other(PARTITIONED_VIEWFACTORS.into()), Step::Other("Mesh2D".into())]);
    }

    #[test]
    fn local_plugin_modules() {
        let text = "ADD_ELMERTEST_MODULE(1dtests Poisson Poisson.f90)\nADD_ELMERTEST_MODULE(1dtests Source Source.f90)\n";
        assert_eq!(local_modules(text), ["Poisson", "Source"]);
    }
}
