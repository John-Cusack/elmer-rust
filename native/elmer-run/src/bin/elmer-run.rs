//! `ElmerSolver`'s process behaviour with [`elmer_run::NativeEngine`] behind
//! `ELMER_NATIVE=1`: `elmer-sif` reads the case and [`elmer_run::mesh_db`]
//! its mesh DB. Cases the engine refuses fall back to the Fortran solver
//! (`ELMER_FORTRAN_SOLVER`) as in `elmer-cli`.

fn main() {
    let engine = elmer_run::NativeEngine::new(Box::new(elmer_run::SifFrontend { mesh: elmer_run::mesh_db::load }));
    elmer_cli::driver::main_entry_with(elmer_cli::strangler::Variant::Serial, &engine)
}
