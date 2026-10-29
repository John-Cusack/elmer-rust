//! `mpirun -np N elmer-mpi-run [case.sif]`: a partitioned run of the
//! Rust-native Elmer, the native counterpart of `ElmerSolver_mpi`. Without
//! an argument the case is the one `ELMERSOLVER_STARTINFO` names. Rank 0
//! prints Elmer's messages; exit codes are `elmer-cli`'s (0 solved, 1
//! failed, 2 bad input, 78 refused). Run with `HWLOC_COMPONENTS=-gl`.

use std::path::PathBuf;

use elmer_cli::driver::EngineOutcome;
use elmer_cli::{ElmerLogger, ExitCode};
use elmer_core::{Communicator, Logger, NullLogger};
use elmer_mpi::MpiRuntime;
use elmer_mpi_run::{MpiEngine, WorldExchange};

fn case_path() -> Option<PathBuf> {
    let cwd = std::env::current_dir().ok()?;
    match std::env::args_os().nth(1) {
        Some(a) => Some(cwd.join(a)),
        None => elmer_cli::args::read_startinfo(&cwd).map(|s| cwd.join(s.trim())),
    }
}

fn main() {
    let code = {
        let runtime = match MpiRuntime::init() {
            Ok(r) => r,
            Err(e) => {
                eprintln!("elmer-mpi-run: {e}");
                std::process::exit(ExitCode::Input.code());
            }
        };
        let world = runtime.world();
        let rank = world.rank();
        let elmer_log = ElmerLogger::stdout();
        let log: &dyn Logger = if rank == 0 { &elmer_log } else { &NullLogger };
        let code = match case_path() {
            Some(sif) if sif.is_file() => {
                let engine = MpiEngine::new(Box::new(WorldExchange(world)));
                let fatal = |caller: &str, msg: &str| {
                    if rank == 0 {
                        elmer_log.fatal(caller, msg);
                    }
                };
                match engine.run_sif(&sif, log, &fatal) {
                    EngineOutcome::Done(c) => c,
                    EngineOutcome::Refused(reasons) => {
                        if rank == 0 {
                            for r in reasons {
                                eprintln!("elmer-mpi-run: refused: {r}");
                            }
                        }
                        ExitCode::Refuse
                    }
                }
            }
            other => {
                if rank == 0 {
                    eprintln!("elmer-mpi-run: no case file ({other:?})");
                }
                ExitCode::Input
            }
        };
        // Every rank reaches here with the same outcome (errors are
        // collective), so MPI_Finalize runs on all of them.
        drop(runtime);
        code
    };
    std::process::exit(code.code());
}
