//! `ElmerSolver`'s process behaviour with [`elmer_run::NativeEngine`] behind
//! `ELMER_NATIVE=1`. Until a `.sif`/mesh front end is on main the engine
//! refuses every case, and the Fortran fallback (`ELMER_FORTRAN_SOLVER`)
//! applies as in `elmer-cli`.

fn main() {
    elmer_cli::driver::main_entry_with(elmer_cli::strangler::Variant::Serial, &elmer_run::NativeEngine::default())
}
