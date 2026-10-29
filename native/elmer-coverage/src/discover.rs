//! The quick suite, as ctest registers it: `ctest -N -V -L quick` lists
//! every quick test with its working directory, task count and source
//! directory without running anything.

use std::path::{Path, PathBuf};
use std::process::Command;

/// One quick test.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuickTest {
    pub name: String,
    /// The configured directory ctest runs in (`<build>/fem/tests/<dir>`).
    pub work_dir: PathBuf,
    /// The test's source directory (`-DTEST_SOURCE`), holding `runtest.cmake`.
    pub source_dir: PathBuf,
    /// `-DMPIEXEC_NTASKS`.
    pub nprocs: u32,
    pub labels: Vec<String>,
}

/// `-DKEY=value` out of ctest's quoted test command.
fn define<'a>(command: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"-D{key}=");
    let start = command.find(&needle)? + needle.len();
    let end = command[start..].find('"')?;
    Some(&command[start..start + end])
}

/// Parses the output of `ctest -N -V`. Each test prints `N: Test command:`,
/// `N: Working Directory:` and `Labels:` before its `Test #N: name` line.
pub fn parse_ctest(text: &str) -> Vec<QuickTest> {
    let mut tests = Vec::new();
    let (mut command, mut work, mut labels) = (String::new(), String::new(), Vec::new());
    for line in text.lines() {
        if let Some(pos) = line.find("Test command: ") {
            command = line[pos + "Test command: ".len()..].to_string();
        } else if let Some(pos) = line.find("Working Directory: ") {
            work = line[pos + "Working Directory: ".len()..].trim().to_string();
        } else if let Some(rest) = line.strip_prefix("Labels: ") {
            labels = rest.split_whitespace().map(str::to_string).collect();
        } else if let Some(rest) = line.trim_start().strip_prefix("Test").map(str::trim_start).and_then(|r| r.strip_prefix('#')) {
            let Some((_, name)) = rest.split_once(": ") else { continue };
            let nprocs = define(&command, "MPIEXEC_NTASKS").and_then(|n| n.parse().ok()).unwrap_or(1);
            tests.push(QuickTest {
                name: name.trim().to_string(),
                work_dir: PathBuf::from(&work),
                source_dir: PathBuf::from(define(&command, "TEST_SOURCE").unwrap_or("")),
                nprocs,
                labels: std::mem::take(&mut labels),
            });
        }
    }
    tests
}

/// Lists the tests carrying `label` in the pristine build tree. This reads
/// ctest's registry only; no test is run.
pub fn list(build_dir: &Path, label: &str) -> Result<Vec<QuickTest>, String> {
    let out = Command::new("nice")
        .args(["-n", "10", "ctest", "-N", "-V", "-L", label])
        .current_dir(build_dir)
        .output()
        .map_err(|e| format!("ctest: {e}"))?;
    if !out.status.success() {
        return Err(format!("ctest -N failed: {}", String::from_utf8_lossy(&out.stderr)));
    }
    Ok(parse_ctest(&String::from_utf8_lossy(&out.stdout)))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"
745: Test command: /usr/bin/cmake "-DELMERGRID_BIN=/b/ElmerGrid" "-DTEST_SOURCE=/src/ThermalPipe2D" "-DMPIEXEC_NTASKS=1" "-P" "/src/ThermalPipe2D/runtest.cmake"
745: Working Directory: /b/fem/tests/ThermalPipe2D
Labels: heateq quick serial
  Test #745: ThermalPipe2D
7: Test command: /usr/bin/cmake "-DTEST_SOURCE=/src/Foo" "-DMPIEXEC_NTASKS=2" "-P" "/src/Foo/runtest.cmake"
7: Working Directory: /b/fem/tests/Foo
Labels: quick parallel
  Test    #7: Foo_np2
"#;

    #[test]
    fn parses_name_dirs_and_tasks() {
        let t = parse_ctest(SAMPLE);
        assert_eq!(t.len(), 2);
        assert_eq!(t[0].name, "ThermalPipe2D");
        assert_eq!(t[0].work_dir, PathBuf::from("/b/fem/tests/ThermalPipe2D"));
        assert_eq!(t[0].source_dir, PathBuf::from("/src/ThermalPipe2D"));
        assert_eq!(t[0].nprocs, 1);
        assert_eq!(t[0].labels, ["heateq", "quick", "serial"]);
        assert_eq!((t[1].name.as_str(), t[1].nprocs), ("Foo_np2", 2));
    }
}
