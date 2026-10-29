//! Running a solver executable the way ElmerGUI does.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

#[derive(Clone, Debug)]
pub struct Outcome {
    /// `None` when the process was killed by a signal or by our timeout.
    pub exit: Option<i32>,
    pub timed_out: bool,
    /// stdout followed by stderr (the GUI reads both channels).
    pub log: String,
    pub started: bool,
}

/// Variables that make `elmer-cli` hand a case to the Fortran solver.
/// `scrub_fortran` removes them so the executable has to answer by itself.
pub const FORTRAN_VARS: [&str; 2] = ["ELMER_FORTRAN_SOLVER", "ELMER_FORTRAN_SOLVER_MPI"];

/// Run `exe` with `args` in `cwd`, killing it after `timeout`. The
/// environment is inherited except for the variables that keep Elmer serial
/// and quiet about hwloc. `Err` means the process could not be started,
/// which is the one condition the GUI's own probe tests (`waitForStarted`).
pub fn run(
    exe: &Path,
    args: &[&str],
    cwd: &Path,
    timeout: Duration,
    scrub_fortran: bool,
) -> std::io::Result<Outcome> {
    let mut cmd = Command::new(exe);
    if scrub_fortran {
        for var in FORTRAN_VARS {
            cmd.env_remove(var);
        }
    }
    let mut child = cmd
        .args(args)
        .current_dir(cwd)
        .env("OMP_NUM_THREADS", "1")
        .env("HWLOC_COMPONENTS", "-gl")
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    let mut so = child.stdout.take().expect("piped");
    let mut se = child.stderr.take().expect("piped");
    let t_out = std::thread::spawn(move || {
        let mut b = Vec::new();
        let _ = so.read_to_end(&mut b);
        b
    });
    let t_err = std::thread::spawn(move || {
        let mut b = Vec::new();
        let _ = se.read_to_end(&mut b);
        b
    });
    let start = Instant::now();
    let (status, timed_out) = loop {
        if let Some(s) = child.try_wait()? {
            break (Some(s), false);
        }
        if start.elapsed() >= timeout {
            let _ = child.kill();
            let _ = child.wait();
            break (None, true);
        }
        std::thread::sleep(Duration::from_millis(20));
    };
    let mut log = String::from_utf8_lossy(&t_out.join().unwrap_or_default()).into_owned();
    log.push_str(&String::from_utf8_lossy(&t_err.join().unwrap_or_default()));
    Ok(Outcome {
        exit: status.and_then(|s| s.code()),
        timed_out,
        log,
        started: true,
    })
}

/// A solver executable named on the command line must be absolute; relative
/// names would resolve through `PATH` and silently pick another install.
pub fn require_absolute(p: &str) -> Result<PathBuf, String> {
    let path = PathBuf::from(p);
    if !path.is_absolute() {
        return Err(format!("`{p}` is not an absolute path"));
    }
    if !path.is_file() {
        return Err(format!("`{p}` is not a file"));
    }
    Ok(path)
}
