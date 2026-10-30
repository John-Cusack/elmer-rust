//! ctest registration of one arm build, read from its generated
//! `CTestTestfile.cmake` files.
//!
//! The registration is the denominator: a test directory whose
//! `CMakeLists.txt` says `quick` but sits behind `IF(HAVE_LUA)` is not a
//! registered test in an arm built without Lua, and an `NPROCS 2 4` test
//! registers as `X_np2` and `X_np4`. Reading the generated files gives
//! exactly what `ctest` would run, without needing `ctest` itself.

use crate::{io_with, OracleError, Result};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// One registered ctest test.
#[derive(Debug, Clone, PartialEq)]
pub struct RegisteredTest {
    /// ctest name, e.g. `AngleMetisLayer_np2`.
    pub name: String,
    /// Test source directory name under `fem/tests` (basename of
    /// `TEST_SOURCE`).
    pub dir: String,
    /// `MPIEXEC_NTASKS`.
    pub np: u32,
    pub labels: Vec<String>,
    /// The `add_test` command verbatim, as the build recorded it
    /// (container paths, before any host remapping).
    pub command: Vec<String>,
    /// `WORKING_DIRECTORY` as recorded.
    pub working_dir: String,
}

impl RegisteredTest {
    /// Value of a `-DNAME=value` argument of the recorded command.
    pub fn define(&self, name: &str) -> Option<&str> {
        let prefix = format!("-D{name}=");
        self.command.iter().find_map(|a| a.strip_prefix(prefix.as_str()))
    }

    pub fn has_label(&self, label: &str) -> bool {
        self.labels.iter().any(|l| l == label)
    }
}

/// Quoted arguments of one CMake call, `name(...)` already stripped.
fn quoted_args(body: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_quote = false;
    let mut chars = body.chars();
    while let Some(c) = chars.next() {
        match c {
            '\\' if in_quote => {
                if let Some(n) = chars.next() {
                    cur.push(n);
                }
            }
            '"' => {
                if in_quote {
                    out.push(std::mem::take(&mut cur));
                }
                in_quote = !in_quote;
            }
            _ if in_quote => cur.push(c),
            _ if c.is_whitespace() => {}
            _ => {
                // Unquoted token (test names, PROPERTIES keywords).
                cur.push(c);
                let mut rest = String::new();
                for n in chars.by_ref() {
                    if n.is_whitespace() {
                        break;
                    }
                    rest.push(n);
                }
                cur.push_str(&rest);
                out.push(std::mem::take(&mut cur));
            }
        }
    }
    out
}

/// Body between the first `(` and the last `)` of a call line.
fn call_body<'a>(line: &'a str, head: &str) -> Option<&'a str> {
    let rest = line.strip_prefix(head)?.strip_prefix('(')?;
    let end = rest.rfind(')')?;
    Some(&rest[..end])
}

/// Parse one `CTestTestfile.cmake`, returning its tests and its `subdirs`.
pub fn parse_ctest_file(text: &str) -> (Vec<RegisteredTest>, Vec<String>) {
    let mut tests: BTreeMap<String, RegisteredTest> = BTreeMap::new();
    let mut subdirs = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if let Some(body) = call_body(line, "subdirs") {
            subdirs.extend(quoted_args(body));
        } else if let Some(body) = call_body(line, "add_test") {
            let args = quoted_args(body);
            let Some((name, command)) = args.split_first() else {
                continue;
            };
            let command = command.to_vec();
            let source = command
                .iter()
                .find_map(|a| a.strip_prefix("-DTEST_SOURCE="))
                .unwrap_or_default();
            let dir = source.rsplit('/').next().unwrap_or_default().to_string();
            let np = command
                .iter()
                .find_map(|a| a.strip_prefix("-DMPIEXEC_NTASKS="))
                .and_then(|v| v.parse().ok())
                .unwrap_or(1);
            tests.insert(
                name.clone(),
                RegisteredTest {
                    name: name.clone(),
                    dir,
                    np,
                    labels: Vec::new(),
                    command,
                    working_dir: String::new(),
                },
            );
        } else if let Some(body) = call_body(line, "set_tests_properties") {
            let args = quoted_args(body);
            let Some(t) = args.first().and_then(|n| tests.get_mut(n)) else {
                continue;
            };
            let mut i = 1;
            while i + 1 < args.len() {
                match args[i].as_str() {
                    "LABELS" => {
                        t.labels = args[i + 1].split(';').map(String::from).collect();
                        i += 2;
                    }
                    "WORKING_DIRECTORY" => {
                        t.working_dir = args[i + 1].clone();
                        i += 2;
                    }
                    "PROPERTIES" => i += 1,
                    _ => i += 2,
                }
            }
        }
    }
    (tests.into_values().collect(), subdirs)
}

/// Every test registered under `<build>/fem/tests`, sorted by name.
pub fn read_registry(build_dir: &Path) -> Result<Vec<RegisteredTest>> {
    read_registry_under(&build_dir.join("fem").join("tests"))
}

/// Every test registered under `root` (a binary-tree directory with a
/// `CTestTestfile.cmake`, e.g. `<build>/elmerice` for the ice suite).
pub fn read_registry_under(root: &Path) -> Result<Vec<RegisteredTest>> {
    let root = root.to_path_buf();
    let mut out = Vec::new();
    let mut queue: Vec<PathBuf> = vec![root.clone()];
    while let Some(dir) = queue.pop() {
        let path = dir.join("CTestTestfile.cmake");
        if !path.exists() {
            continue;
        }
        let text = std::fs::read_to_string(&path).map_err(|e| io_with(&path, e))?;
        let (tests, subdirs) = parse_ctest_file(&text);
        out.extend(tests);
        for s in subdirs {
            queue.push(dir.join(s));
        }
    }
    if out.is_empty() {
        return Err(OracleError::Parse(format!(
            "no ctest registration under {}",
            root.display()
        )));
    }
    out.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(out)
}

/// Rewrite the container prefix the arm was configured under (the parent
/// of the recorded `BINARY_DIR`) to the host prefix of `build_dir`.
pub fn remap_arg(arg: &str, recorded_binary_dir: &str, build_dir: &Path) -> String {
    let container_root = recorded_binary_dir.rsplit_once('/').map(|(p, _)| p).unwrap_or("");
    let host_root = build_dir.parent().unwrap_or(build_dir).display().to_string();
    if container_root.is_empty() || container_root == host_root {
        return arg.to_string();
    }
    arg.replace(&format!("{container_root}/"), &format!("{host_root}/"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"# CMake generated Testfile for
add_test(Angle_np2 "/usr/bin/cmake" "-DCMAKE_MODULE_PATH=/work/b.src/cmake/Modules" "/usr/share/cmake-3.28/Modules" "-DELMERSOLVER_BIN=/work/b/fem/src/ElmerSolver_mpi" "-DTEST_SOURCE=/work/b.src/fem/tests/Angle" "-DBINARY_DIR=/work/b" "-DMPIEXEC_NTASKS=2" "-P" "/work/b.src/fem/tests/Angle/runtest.cmake")
set_tests_properties(Angle_np2 PROPERTIES  DEPENDS "Angle" LABELS "parallel;heateq;quick" PROCESSORS "2" WORKING_DIRECTORY "/work/b/fem/tests/Angle" _BACKTRACE_TRIPLES "x;1;y")
subdirs("Angle")
"#;

    #[test]
    fn parses_name_dir_np_labels_and_subdirs() {
        let (tests, subdirs) = parse_ctest_file(SAMPLE);
        assert_eq!(subdirs, vec!["Angle"]);
        assert_eq!(tests.len(), 1);
        let t = &tests[0];
        assert_eq!(t.name, "Angle_np2");
        assert_eq!(t.dir, "Angle");
        assert_eq!(t.np, 2);
        assert!(t.has_label("quick"));
        assert!(t.has_label("parallel"));
        assert_eq!(t.working_dir, "/work/b/fem/tests/Angle");
        assert_eq!(t.define("BINARY_DIR"), Some("/work/b"));
        // The stray list element CMake split off stays a separate arg.
        assert_eq!(t.command[2], "/usr/share/cmake-3.28/Modules");
    }

    #[test]
    fn remaps_container_prefix_to_host() {
        let host = Path::new("/home/x/data/b");
        assert_eq!(
            remap_arg("-DTEST_SOURCE=/work/b.src/fem/tests/A", "/work/b", host),
            "-DTEST_SOURCE=/home/x/data/b.src/fem/tests/A"
        );
        assert_eq!(remap_arg("/usr/bin/mpiexec", "/work/b", host), "/usr/bin/mpiexec");
    }
}
