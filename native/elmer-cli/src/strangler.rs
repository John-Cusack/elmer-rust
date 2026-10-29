//! Strangler launcher: hands a case to the pristine Fortran binary named by
//! `ELMER_FORTRAN_SOLVER`. The Fortran code is never linked into this crate; it
//! is only executed as a separate program.

use std::ffi::OsString;
use std::fmt;
use std::io;
use std::os::unix::process::{CommandExt, ExitStatusExt};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Path of the pristine serial `ElmerSolver`.
pub const ENV_FORTRAN: &str = "ELMER_FORTRAN_SOLVER";
/// Path of the pristine `ElmerSolver_mpi`; defaults to [`ENV_FORTRAN`].
pub const ENV_FORTRAN_MPI: &str = "ELMER_FORTRAN_SOLVER_MPI";
/// Set on the delegated process so a mis-set path cannot loop back into us.
pub const ENV_GUARD: &str = "ELMER_CLI_DELEGATED";

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Variant {
    Serial,
    Mpi,
}

#[derive(Debug)]
pub enum DelegateError {
    NotConfigured,
    Missing(PathBuf),
    /// The configured path is this very executable, or we were already
    /// delegated to.
    Recursion(PathBuf),
}

impl fmt::Display for DelegateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DelegateError::NotConfigured => write!(f, "{ENV_FORTRAN} is not set"),
            DelegateError::Missing(p) => write!(f, "Fortran solver {} is not an executable file", p.display()),
            DelegateError::Recursion(p) => write!(f, "Fortran solver {} would re-enter elmer-cli", p.display()),
        }
    }
}

fn is_executable(p: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    std::fs::metadata(p).is_ok_and(|m| m.is_file() && m.permissions().mode() & 0o111 != 0)
}

/// Picks and validates the Fortran executable for `variant`. `get` reads the
/// environment; `current_exe` is this process's own image, for the recursion
/// guard (canonical paths are compared).
pub fn resolve_fortran(
    variant: Variant,
    get: &dyn Fn(&str) -> Option<OsString>,
    current_exe: Option<&Path>,
) -> Result<PathBuf, DelegateError> {
    let serial = get(ENV_FORTRAN).filter(|v| !v.is_empty());
    let path = match variant {
        Variant::Serial => serial,
        Variant::Mpi => get(ENV_FORTRAN_MPI).filter(|v| !v.is_empty()).or(serial),
    }
    .map(PathBuf::from)
    .ok_or(DelegateError::NotConfigured)?;
    if get(ENV_GUARD).is_some() {
        return Err(DelegateError::Recursion(path));
    }
    if !is_executable(&path) {
        return Err(DelegateError::Missing(path));
    }
    if let (Ok(a), Some(Ok(b))) = (path.canonicalize(), current_exe.map(Path::canonicalize)) {
        if a == b {
            return Err(DelegateError::Recursion(path));
        }
    }
    Ok(path)
}

/// Replaces this process with `path`, keeping cwd, environment (plus
/// [`ENV_GUARD`]) and signal dispositions. Returns only on failure.
pub fn exec_fortran(path: &Path, args: &[OsString]) -> io::Error {
    Command::new(path).args(args).env(ENV_GUARD, "1").exec()
}

/// Runs `path` as a child and returns its exit status (`128 + signal` when it
/// was killed). Used by shadow mode.
pub fn run_child(path: &Path, args: &[OsString]) -> io::Result<i32> {
    let st = Command::new(path).args(args).env(ENV_GUARD, "1").status()?;
    Ok(st.code().unwrap_or_else(|| 128 + st.signal().unwrap_or(0)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn envmap(pairs: &[(&str, &str)]) -> HashMap<String, OsString> {
        pairs.iter().map(|(k, v)| (k.to_string(), OsString::from(v))).collect()
    }

    #[test]
    fn unset_and_missing() {
        let e = envmap(&[]);
        assert!(matches!(
            resolve_fortran(Variant::Serial, &|k| e.get(k).cloned(), None),
            Err(DelegateError::NotConfigured)
        ));
        let e = envmap(&[(ENV_FORTRAN, "/nonexistent/ElmerSolver")]);
        assert!(matches!(
            resolve_fortran(Variant::Serial, &|k| e.get(k).cloned(), None),
            Err(DelegateError::Missing(_))
        ));
    }

    #[test]
    fn guard_blocks_second_delegation() {
        let e = envmap(&[(ENV_FORTRAN, "/bin/sh"), (ENV_GUARD, "1")]);
        assert!(matches!(
            resolve_fortran(Variant::Serial, &|k| e.get(k).cloned(), None),
            Err(DelegateError::Recursion(_))
        ));
    }

    #[test]
    fn self_path_is_recursion() {
        let e = envmap(&[(ENV_FORTRAN, "/bin/sh")]);
        let me = Path::new("/bin/sh");
        assert!(matches!(
            resolve_fortran(Variant::Serial, &|k| e.get(k).cloned(), Some(me)),
            Err(DelegateError::Recursion(_))
        ));
        assert!(resolve_fortran(Variant::Serial, &|k| e.get(k).cloned(), Some(Path::new("/bin/ls"))).is_ok());
    }

    #[test]
    fn mpi_variant_prefers_its_own_path() {
        let e = envmap(&[(ENV_FORTRAN, "/bin/sh"), (ENV_FORTRAN_MPI, "/bin/ls")]);
        let get = |k: &str| e.get(k).cloned();
        assert_eq!(resolve_fortran(Variant::Mpi, &get, None).unwrap(), PathBuf::from("/bin/ls"));
        assert_eq!(resolve_fortran(Variant::Serial, &get, None).unwrap(), PathBuf::from("/bin/sh"));
    }
}
