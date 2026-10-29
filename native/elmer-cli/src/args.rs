//! Command line of `ElmerSolver` (ElmerSolver.F90:190-240, 335-352).

use std::fmt;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CliArgs {
    /// `-s` / `--silent`.
    pub silent: bool,
    /// `-v` / `--version`.
    pub version: bool,
    /// First argument when it does not start with `-`.
    pub sif: Option<String>,
    /// Second argument after a model name (`eq`).
    pub equation: Option<String>,
    /// `-rpar n v1 … vn`.
    pub rpar: Vec<f64>,
    /// `-ipar n v1 … vn`.
    pub ipar: Vec<i32>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct ArgError(pub String);

impl fmt::Display for ArgError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ArgError {}

/// Mirrors the Fortran scan: every argument is examined for `-rpar`,
/// `-ipar`, `-s`/`--silent`, `-v`/`--version`; the model name is argument 1
/// only. Unknown options are ignored, as in Fortran.
pub fn parse(args: &[String]) -> Result<CliArgs, ArgError> {
    let mut out = CliArgs::default();
    if let Some(first) = args.first() {
        if !first.starts_with('-') {
            out.sif = Some(first.clone());
            out.equation = args.get(1).cloned();
        }
    }
    let mut i = 0;
    while i < args.len() {
        let a = args[i].as_str();
        match a {
            "-rpar" | "-ipar" => {
                let n: usize = args
                    .get(i + 1)
                    .and_then(|s| s.trim().parse().ok())
                    .ok_or_else(|| ArgError(format!("{a}: missing or invalid parameter count")))?;
                let vals = args
                    .get(i + 2..i + 2 + n)
                    .ok_or_else(|| ArgError(format!("{a}: expected {n} values")))?;
                if a == "-rpar" {
                    for v in vals {
                        out.rpar.push(parse_real(v).ok_or_else(|| ArgError(format!("-rpar: bad value `{v}`")))?);
                    }
                } else {
                    for v in vals {
                        out.ipar.push(v.trim().parse().map_err(|_| ArgError(format!("-ipar: bad value `{v}`")))?);
                    }
                }
                i += 2 + n;
                continue;
            }
            "-s" | "--silent" => out.silent = true,
            "-v" | "--version" => out.version = true,
            _ => {}
        }
        i += 1;
    }
    Ok(out)
}

/// Fortran list-directed real read accepts `1.d0` style exponents.
fn parse_real(s: &str) -> Option<f64> {
    let t = s.trim().replace(['d', 'D'], "e");
    t.parse().ok()
}

/// First line of `ELMERSOLVER_STARTINFO` in `dir` (`READ(1,'(a)') ModelName`),
/// trailing blanks dropped. `None` when the file cannot be opened.
pub fn read_startinfo(dir: &Path) -> Option<String> {
    let text = std::fs::read_to_string(dir.join("ELMERSOLVER_STARTINFO")).ok()?;
    let line = text.lines().next().unwrap_or("");
    Some(line.trim_end().to_string())
}

/// Resolves the model file as `Fortran` does; relative names are relative to
/// `cwd`.
pub fn resolve_sif(args: &CliArgs, cwd: &Path) -> Option<PathBuf> {
    let name = match &args.sif {
        Some(s) => s.clone(),
        None => read_startinfo(cwd)?,
    };
    Some(cwd.join(name))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn plain_case_and_flags() {
        let a = parse(&v(&["case.sif", "heat", "-s"])).unwrap();
        assert_eq!(a.sif.as_deref(), Some("case.sif"));
        assert_eq!(a.equation.as_deref(), Some("heat"));
        assert!(a.silent && !a.version);
    }

    #[test]
    fn probe_has_no_model() {
        for f in ["-v", "--version"] {
            let a = parse(&v(&[f])).unwrap();
            assert!(a.version && a.sif.is_none());
        }
    }

    #[test]
    fn rpar_ipar_consume_values() {
        let a = parse(&v(&["x.sif", "-rpar", "2", "1.5", "2.d0", "-ipar", "1", "7", "-v"])).unwrap();
        assert_eq!(a.rpar, vec![1.5, 2.0]);
        assert_eq!(a.ipar, vec![7]);
        assert!(a.version);
        assert!(parse(&v(&["-rpar", "3", "1"])).is_err());
        assert!(parse(&v(&["-ipar", "x"])).is_err());
    }

    #[test]
    fn value_that_looks_like_flag_is_not_a_flag() {
        // `-rpar 1 -v` consumes `-v` as a (bad) value, as Fortran's READ would fail.
        assert!(parse(&v(&["-rpar", "1", "-v"])).is_err());
    }

    #[test]
    fn startinfo_first_line() {
        let dir = std::env::temp_dir().join(format!("elmer-cli-args-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("ELMERSOLVER_STARTINFO"), "case.sif  \nignored\n").unwrap();
        assert_eq!(read_startinfo(&dir).as_deref(), Some("case.sif"));
        let p = resolve_sif(&CliArgs::default(), &dir).unwrap();
        assert_eq!(p, dir.join("case.sif"));
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(read_startinfo(&dir), None);
    }
}
