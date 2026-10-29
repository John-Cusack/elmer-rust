//! Which provider serves which seam of `elmer-solver`.

use elmer_core::{DirectSolver, ParamMap, Params};
use elmer_direct::{kind_for_method, BandCholesky, BandLu, DirectKind, FaerCholesky, FaerLu};
use elmer_solver::{Simulation, SolverError, SolverRegistry};

/// `Linear System Direct Method` values `elmer-direct` serves
/// (`DirectSolve.F90:4266`).
const DIRECT_METHODS: [&str; 14] = [
    "banded",
    "symmetric banded",
    "umfpack",
    "big umfpack",
    "smumps",
    "cmumps",
    "mumps",
    "dmumps",
    "zmumps",
    "mumpslocal",
    "superlu",
    "pardiso",
    "cpardiso",
    "cholmod",
];

/// The procedures this build can run.
pub fn solver_registry() -> SolverRegistry {
    let mut reg = SolverRegistry::new();
    elmer_physics_heat::register(&mut reg);
    reg
}

/// Adds what `elmer-solver`'s default backends (assembly, time) leave out:
/// `elmer-direct` for the direct methods. Iterative methods stay unregistered,
/// and so refused, until `elmer-solver` wires `elmer-krylov` and
/// `elmer-itermethods` with `Linear System Scaling` (L12).
pub fn install_backends(sim: &mut Simulation) {
    for m in DIRECT_METHODS {
        sim.backends.linear.register_direct(m, direct);
    }
}

/// `elmer-direct`'s choice for the section. The method is read here, not
/// taken from the registry key, because Elmer's default without the keyword
/// is `banded` (DirectSolve.F90:4311).
fn direct(params: &ParamMap) -> Result<Box<dyn DirectSolver<f64> + Send>, SolverError> {
    let kind = match params.string(elmer_direct::select::METHOD_KEY) {
        Some(m) => kind_for_method(m)?,
        None => DirectKind::BandLu,
    };
    Ok(match kind {
        DirectKind::BandLu => Box::new(BandLu::<f64>::new()),
        DirectKind::BandCholesky => Box::new(BandCholesky::<f64>::new()),
        DirectKind::SparseLu => Box::new(FaerLu::<f64>::new()),
        DirectKind::SparseCholesky => Box::new(FaerCholesky::<f64>::new()),
    })
}
