//! Which provider serves which seam of `elmer-solver`. Assembly, time and the
//! linear solvers are `elmer-solver`'s default backends; this crate adds the
//! physics and output procedures, the element library they and the
//! framework's own integrals share and the MATC engine.

use std::path::Path;
use std::sync::OnceLock;

use elmer_core::{ElementBasis, Logger, MeshView, Params, Quadrature};
use elmer_fe::{CoordinateSystem, ElmerQuadrature, LagrangeBasis};
use elmer_mortar::{generate_periodic_projectors, PeriodicAssembly, ProjectorEnv};
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
    elmer_savedata::register(&mut reg);
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

/// `GeneratePeriodicProjectors`: with a `Periodic BC`, `elmer-mortar`'s
/// projectors and its `PeriodicAssembly` (the default assembly with the
/// projector couplings and `SetPeriodicBoundaries`) as the assembly
/// backend. `Save Projector` files go to `save_dir`.
pub fn periodic_assembly(
    model: &Model,
    mesh: &dyn MeshView,
    cdim: usize,
    save_dir: &Path,
    log: &dyn Logger,
) -> Result<Option<PeriodicAssembly>, SolverError> {
    if !model.boundary_conditions.iter().any(|bc| bc.get("periodic bc").is_some()) {
        return Ok(None);
    }
    // The Fortran stops the run once the projector is saved; refused before
    // any file is written.
    if model.boundary_conditions.iter().any(|bc| bc.logical("save projector and stop") == Some(true)) {
        return Err(SolverError::Unsupported("Save Projector And Stop".into()));
    }
    let bcs: Vec<&dyn Params> = model.boundary_conditions.iter().map(|b| b as &dyn Params).collect();
    let bc_of = |tag: i32| model.bc_of_tag(tag);
    let env = ProjectorEnv {
        mesh,
        bcs: &bcs,
        bc_of: &bc_of,
        simulation: &model.simulation,
        solver: model.solvers.first().map(|s| s as &dyn Params),
        coord_dim: cdim,
        log,
    };
    let (projectors, _) = generate_periodic_projectors(&env, save_dir, 0)?;
    Ok(Some(PeriodicAssembly::new(model.boundary_conditions.clone(), projectors, mesh.n_nodes())))
}
