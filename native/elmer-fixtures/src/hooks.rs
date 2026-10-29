//! The ported kernels behind `elmer-solver`'s hooks: [`PoissonAssembly`] as
//! the [`ElementAssembly`] routines `Bulk Assembly Procedure` and `Boundary
//! Assembly Procedure` name, fed the way the Fortran's `DefUtils` calls feed
//! it, and integrated with `elmer-fe`.

use elmer_core::{ElemId, LocalSystem, MeshView, ParamMap, Params};
use elmer_fe::{CoordinateSystem, ElmerQuadrature, LagrangeBasis};
use elmer_solver::{ElementAssembly, SolverCtx, SolverError};

use crate::poisson_assembly::{BoundaryProps, BulkProps, Element, PoissonAssembly};

/// `ElementInfo` for `CoordinateSystemDimension()`: the case's `Coordinate
/// System` on its mesh, as `elmer-run` installs `elmer-fe` for the physics
/// crates.
pub(crate) fn lagrange_basis(ctx: &SolverCtx<'_>) -> Result<LagrangeBasis, SolverError> {
    Ok(LagrangeBasis::new(coordinate_system(ctx)?.dim))
}

/// `SetCoordinateSystem`: the case's `Coordinate System` on its mesh.
pub(crate) fn coordinate_system(ctx: &SolverCtx<'_>) -> Result<CoordinateSystem, SolverError> {
    let keyword = ctx.simulation_params().string("coordinate system");
    CoordinateSystem::from_keyword(keyword, ctx.mesh.dim()).ok_or_else(|| {
        SolverError::fatal("SetCoordinateSystem", format!("Unknown Coordinate System: {}", keyword.unwrap_or("")))
    })
}

/// `GetElementNodes`: the element's node coordinates in Elmer order.
pub(crate) fn element_nodes(mesh: &dyn MeshView, e: ElemId) -> Vec<[f64; 3]> {
    mesh.element(e).nodes.iter().map(|&n| mesh.coord(n)).collect()
}

/// `GetReal(list, key)` at element `e`'s nodes; `None` for an absent keyword
/// or list, which `GetReal` answers with zeros.
fn real(ctx: &SolverCtx<'_>, list: Option<&ParamMap>, key: &str, e: ElemId) -> Result<Option<Vec<f64>>, SolverError> {
    match list {
        Some(p) => ctx.real_at_nodes(p, key, e),
        None => Ok(None),
    }
}

/// `PoissonAssembly.f90`'s `BulkAssembly`: `Conductivity`, `Source` and
/// `Reaction` of the element's material (`GetMaterial`).
#[derive(Default)]
pub struct PoissonBulk(PoissonAssembly);

impl ElementAssembly for PoissonBulk {
    fn assemble(
        &mut self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        _dt: f64,
        _transient: bool,
        local: &mut LocalSystem,
    ) -> Result<(), SolverError> {
        let el = ctx.mesh.element(e);
        let xyz = element_nodes(ctx.mesh, e);
        let material = ctx.body_params(e).and_then(|b| b.material);
        let cond = real(ctx, material, "conductivity", e)?;
        let source = real(ctx, material, "source", e)?;
        let reaction = real(ctx, material, "reaction", e)?;
        let props = BulkProps { conductivity: cond.as_deref(), source: source.as_deref(), reaction: reaction.as_deref() };
        self.0
            .bulk(
                Element { etype: el.etype, nodes: &xyz },
                props,
                &lagrange_basis(ctx)?,
                &ElmerQuadrature,
                ctx.log,
                local,
            )
            .map_err(|err| SolverError::fatal("PoissonBulkAssembly", err.to_string()))
    }
}

/// `PoissonAssembly.f90`'s `BoundaryAssembly`: `Flux`, `Coefficient` and
/// `External` of the element's boundary condition (`GetBC`).
#[derive(Default)]
pub struct PoissonBoundary(PoissonAssembly);

impl ElementAssembly for PoissonBoundary {
    fn assemble(
        &mut self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        _dt: f64,
        _transient: bool,
        local: &mut LocalSystem,
    ) -> Result<(), SolverError> {
        let el = ctx.mesh.element(e);
        let xyz = element_nodes(ctx.mesh, e);
        let bc = ctx.bc_params(e);
        let flux = real(ctx, bc, "flux", e)?;
        let coefficient = real(ctx, bc, "coefficient", e)?;
        let external = real(ctx, bc, "external", e)?;
        let props = BoundaryProps {
            flux: flux.as_deref(),
            coefficient: coefficient.as_deref(),
            external: external.as_deref(),
        };
        self.0
            .boundary(
                Element { etype: el.etype, nodes: &xyz },
                props,
                &lagrange_basis(ctx)?,
                &ElmerQuadrature,
                ctx.log,
                local,
            )
            .map(|_| ())
            .map_err(|err| SolverError::fatal("PoissonBoundaryAssembly", err.to_string()))
    }
}
