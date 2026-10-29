//! Process-level policy of `ElmerSolver` / `ElmerSolver_mpi`, as a pure function
//! from (arguments, environment) to an [`Action`], so it is testable without
//! `exec`.
//!
//! Policy (swarm plan C08, L11):
//! * `-v` is always answered natively (the ElmerGUI only checks that it starts).
//! * Default is Fortran-first: when `ELMER_FORTRAN_SOLVER` names a valid
//!   executable the case is delegated without printing anything, so its stdout is
//!   the Fortran stdout. Native is opt-in with `ELMER_NATIVE=1`.
//! * Native mode asks the [`Engine`]; a refusal delegates to Fortran when
//!   configured, else exits [`ExitCode::Refuse`] with the reasons on stderr.
//! * `ELMER_SHADOW=1` runs Fortran as a child (shadow comparison is a stub).
//! * Under `mpirun` with more than one rank native mode is disabled.

use std::collections::HashMap;
use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime};

use crate::args::{self, CliArgs};
use crate::logger::ElmerLogger;
use crate::strangler::{self, DelegateError, Variant};
use crate::version::{print_probe, total_output_silenced, VersionInfo};
use crate::ExitCode;

pub const ENV_NATIVE: &str = "ELMER_NATIVE";
pub const ENV_SHADOW: &str = "ELMER_SHADOW";
const ENV_TOTAL: &str = "ELMERSOLVER_OUTPUT_TOTAL";

/// Environment variables that carry the MPI world size.
const RANK_VARS: [&str; 3] = ["OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_SIZE"];

#[derive(Clone, Debug)]
pub struct Environment {
    pub cwd: PathBuf,
    pub current_exe: Option<PathBuf>,
    vars: HashMap<String, OsString>,
}

impl Environment {
    pub fn from_process() -> Self {
        Environment {
            cwd: std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
            current_exe: std::env::current_exe().ok(),
            vars: std::env::vars_os()
                .filter_map(|(k, v)| Some((k.into_string().ok()?, v)))
                .collect(),
        }
    }

    pub fn new(cwd: PathBuf) -> Self {
        Environment {
            cwd,
            current_exe: None,
            vars: HashMap::new(),
        }
    }

    pub fn with_var(mut self, key: &str, value: &str) -> Self {
        self.vars.insert(key.to_string(), OsString::from(value));
        self
    }

    pub fn get(&self, key: &str) -> Option<OsString> {
        self.vars.get(key).cloned()
    }

    fn flag(&self, key: &str) -> bool {
        self.get(key)
            .and_then(|v| v.into_string().ok())
            .is_some_and(|v| matches!(v.to_ascii_lowercase().as_str(), "1" | "true" | "on" | "yes"))
    }

    fn mpi_ranks(&self) -> usize {
        RANK_VARS
            .iter()
            .filter_map(|k| self.get(k)?.into_string().ok()?.trim().parse().ok())
            .next()
            .unwrap_or(1)
    }
}

/// What the native engine got to work with.
pub struct RunRequest<'a> {
    pub sif: &'a Path,
    pub args: &'a CliArgs,
    pub env: &'a Environment,
}

pub enum EngineOutcome {
    Done(ExitCode),
    /// Not supported natively; the strings are human-readable reasons.
    Refused(Vec<String>),
}

/// The native solver as the CLI sees it. L12/S1 supply the real one.
pub trait Engine {
    fn run(&self, req: &RunRequest<'_>, log: &ElmerLogger) -> EngineOutcome;
}

/// Refuses every case, until a real engine lands.
pub struct NoEngine;

impl Engine for NoEngine {
    fn run(&self, _req: &RunRequest<'_>, _log: &ElmerLogger) -> EngineOutcome {
        EngineOutcome::Refused(vec!["native engine not available".to_string()])
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum Action {
    Exit(ExitCode),
    /// `exec` the Fortran solver with these (original) arguments.
    Delegate { path: PathBuf, args: Vec<OsString> },
    /// Run Fortran as a child, native shadow is a stub.
    Shadow { path: PathBuf, args: Vec<OsString> },
}

/// Decides what the process does. All user-visible text goes to `log` (stdout
/// stream) or `diag` (stderr stream).
pub fn run(
    variant: Variant,
    os_args: &[OsString],
    env: &Environment,
    engine: &dyn Engine,
    log: &ElmerLogger,
    diag: &mut dyn Write,
) -> Action {
    let started = Instant::now();
    let strings: Vec<String> = os_args.iter().map(|a| a.to_string_lossy().into_owned()).collect();
    let cli = match args::parse(&strings) {
        Ok(c) => c,
        Err(e) => {
            log.fatal("MAIN", &e.to_string());
            return Action::Exit(ExitCode::Input);
        }
    };

    if cli.version {
        let total = env.get(ENV_TOTAL).and_then(|v| v.into_string().ok());
        print_probe(
            log,
            &VersionInfo::default(),
            cli.silent,
            total_output_silenced(total.as_deref()),
            SystemTime::now(),
            (0.0, started.elapsed().as_secs_f64()),
        );
        return Action::Exit(ExitCode::Pass);
    }

    let get = |k: &str| env.get(k);
    let fortran = strangler::resolve_fortran(variant, &get, env.current_exe.as_deref());
    let multi_rank = env.mpi_ranks() > 1;
    let native_first = env.flag(ENV_NATIVE) && !multi_rank;

    if env.flag(ENV_SHADOW) {
        if let Ok(path) = fortran {
            return Action::Shadow { path, args: os_args.to_vec() };
        }
    }
    if !native_first {
        match &fortran {
            Ok(path) => return Action::Delegate { path: path.clone(), args: os_args.to_vec() },
            // A configured-but-unusable path is an operator error, not a refusal.
            Err(e @ (DelegateError::Missing(_) | DelegateError::Recursion(_))) => {
                let _ = writeln!(diag, "elmer-cli: {e}");
                return Action::Exit(ExitCode::Input);
            }
            Err(DelegateError::NotConfigured) => {}
        }
    }

    // Native attempt (opt-in, or no Fortran to fall back on).
    let reasons = if multi_rank {
        vec!["native solver is serial-only; refusing under multi-rank MPI".to_string()]
    } else {
        let Some(sif) = args::resolve_sif(&cli, &env.cwd) else {
            log.fatal("MAIN", "Unable to find ELMERSOLVER_STARTINFO, can not execute.");
            return Action::Exit(ExitCode::Input);
        };
        if !sif.is_file() {
            let shown = cli.sif.clone().or_else(|| args::read_startinfo(&env.cwd)).unwrap_or_default();
            log.fatal("MAIN", &format!("Unable to find input file [{shown}], can not execute."));
            return Action::Exit(ExitCode::Input);
        }
        match engine.run(&RunRequest { sif: &sif, args: &cli, env }, log) {
            EngineOutcome::Done(code) => return Action::Exit(code),
            EngineOutcome::Refused(r) => r,
        }
    };

    match fortran {
        Ok(path) => Action::Delegate { path, args: os_args.to_vec() },
        Err(e) => {
            for r in &reasons {
                let _ = writeln!(diag, "elmer-cli: refused: {r}");
            }
            if !matches!(e, DelegateError::NotConfigured) {
                let _ = writeln!(diag, "elmer-cli: {e}");
            }
            let _ = writeln!(diag, "elmer-cli: set {} to run this case with the Fortran solver", strangler::ENV_FORTRAN);
            Action::Exit(ExitCode::Refuse)
        }
    }
}

/// Entry point shared by both binaries.
pub fn main_entry(variant: Variant) -> ! {
    let env = Environment::from_process();
    let os_args: Vec<OsString> = std::env::args_os().skip(1).collect();
    let log = ElmerLogger::stdout();
    let action = run(variant, &os_args, &env, &NoEngine, &log, &mut std::io::stderr());
    match action {
        Action::Exit(code) => std::process::exit(code.code()),
        Action::Delegate { path, args } => {
            let err = strangler::exec_fortran(&path, &args);
            eprintln!("elmer-cli: cannot execute {}: {err}", path.display());
            std::process::exit(ExitCode::Input.code());
        }
        Action::Shadow { path, args } => {
            eprintln!("elmer-cli: shadow mode: native engine not available, comparison skipped");
            match strangler::run_child(&path, &args) {
                Ok(code) => std::process::exit(code),
                Err(err) => {
                    eprintln!("elmer-cli: cannot run {}: {err}", path.display());
                    std::process::exit(ExitCode::Input.code());
                }
            }
        }
    }
}
