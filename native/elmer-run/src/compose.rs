//! Which provider serves which seam of `elmer-solver`. Assembly, time and the
//! linear solvers are `elmer-solver`'s default backends; this crate adds the
//! physics and output procedures, the element library they and the
//! framework's own integrals share and the MATC engine.

use std::sync::OnceLock;

use elmer_core::{ElementBasis, MeshView, Params, Quadrature};
use elmer_fe::{CoordinateSystem, ElmerQuadrature, LagrangeBasis};
use elmer_physics_heat::FeProvider;
use elmer_solver::{FeBackend, Model, SolverError, SolverRegistry};

/// The procedures this build can run.
pub fn solver_registry() -> SolverRegistry {
    let mut reg = SolverRegistry::new();
    elmer_physics_heat::register(&mut reg);
    elmer_physics_flux::register(&mut reg);
    elmer_physics_elastic::register(&mut reg);
    elmer_physics_elasticsolve::register(&mut reg);
    elmer_physics_flow::register(&mut reg);
    elmer_physics_incompressible::register(&mut reg);
    elmer_physics_scalar::register(&mut reg);
    elmer_physics_advdiff::register(&mut reg);
    elmer_physics_meshtools::register(&mut reg);
    elmer_physics_thermoelectric::register(&mut reg);
    elmer_physics_magneto::register(&mut reg);
    elmer_physics_datatofield::register(&mut reg);
    elmer_savedata::register(&mut reg);
    reg
}

/// [`solver_registry`] plus the test-local procedures `elmer-fixtures`
/// ports for the `fem/tests` case named `case` (the case directory's name),
/// as ElmerSolver loads a test's own compiled modules.
pub fn solver_registry_for(case: &str) -> SolverRegistry {
    let mut reg = solver_registry();
    elmer_fixtures::registry::register(case, &mut reg);
    reg
}

/// The case's coordinate system (`SetCoordinateSystem`); its `dim` is
/// `CoordinateSystemDimension()`.
pub fn coordinate_system(model: &Model, mesh: &dyn MeshView) -> Result<CoordinateSystem, SolverError> {
    let keyword = model.simulation.string("coordinate system");
    CoordinateSystem::from_keyword(keyword, mesh.dim())
        .ok_or_else(|| SolverError::fatal("SetCoordinateSystem", format!("Unknown Coordinate System: {}", keyword.unwrap_or(""))))
}

fn lagrange<const CDIM: usize>() -> (Box<dyn ElementBasis + Send>, Box<dyn Quadrature + Send>) {
    (Box::new(LagrangeBasis::new(CDIM)), Box::new(ElmerQuadrature))
}

/// `elmer-fe`'s Lagrange basis and Gauss rules for `CoordinateSystemDimension() = cdim`.
pub fn fe_provider(cdim: usize) -> FeProvider {
    match cdim {
        1 => lagrange::<1>,
        2 => lagrange::<2>,
        _ => lagrange::<3>,
    }
}

/// `elmer-fe` for the integrals `elmer-solver` computes itself (`Coupled
/// Solver` integral constraints, `Apply Source Control` sources).
pub fn fe_backend(cdim: usize) -> FeBackend {
    FeBackend { basis: Box::new(LagrangeBasis::new(cdim)), quadrature: Box::new(ElmerQuadrature) }
}

static INSTALLED: OnceLock<usize> = OnceLock::new();

/// Installs `elmer-fe` into every physics crate that takes it from outside
/// (`elmer-physics-advdiff` builds its own) for `cdim`. The crates keep
/// one process-wide basis, so a second case needing another dimension is
/// refused rather than evaluated with the wrong one.
pub fn install_fe(cdim: usize) -> Result<(), SolverError> {
    let first = *INSTALLED.get_or_init(|| {
        let provider = fe_provider(cdim);
        let _ = elmer_physics_heat::install_fe(provider);
        let _ = elmer_physics_flux::install_fe(provider);
        let _ = elmer_physics_flow::install_fe(provider);
        let _ = elmer_physics_incompressible::install_fe(provider);
        let _ = elmer_physics_scalar::install_fe(provider);
        let _ = elmer_physics_thermoelectric::install_fe(provider);
        let _ = elmer_savedata::install_fe(provider);
        let services = || elmer_physics_elastic::FeServices {
            basis: Box::new(LagrangeBasis::new(cdim)),
            quad: Box::new(ElmerQuadrature),
        };
        let _ = elmer_physics_elastic::install_fe(services());
        let _ = elmer_physics_elasticsolve::install_fe(services());
        cdim
    });
    if first != cdim {
        return Err(SolverError::NotWired(format!(
            "element basis already installed for coordinate dimension {first}, case needs {cdim}"
        )));
    }
    Ok(())
}
