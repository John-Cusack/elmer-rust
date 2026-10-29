//! `SUBROUTINE PoissonSolver` of the `Poisson.f90` plugin that
//! `Procedure = "Poisson" "PoissonSolver"` names in `1dtests`,
//! `1dtests_MeshLevels`, `ExtrudedMesh` and `ExtrudedMeshBaseline1`.
//!
//! Two sources serve those cases. `1dtests/Poisson.f90` (also
//! `1dtests_MeshLevels`) sizes its arrays by nodes; `ExtrudedMesh/Poisson.f90`
//! (also `ExtrudedMeshBaseline1`, and the p-element cases this crate does not
//! register) sizes them by `GetElementNOFDOFs() + GetElementNOFBDOFs()` and
//! condenses the bubble dofs away (`LCondensate`). On Lagrange `n:1`
//! elements there are no bubbles and both reduce to the same element
//! matrix, with the same association order; the bubble path is refused.
//!
//! The solver does one linear solve per visit: `DefaultInitialize`, the
//! active-element loop with the body force `Source` at the element nodes,
//! `DefaultFinishAssembly`, `DefaultDirichletBCs`, `DefaultSolve`. It has no
//! boundary loop and no time term, so a transient case solves the same
//! steady problem with each step's source.

use elmer_core::{BasisEval, ElementBasis, LocalSystem, ParamMap, Params, Quadrature};
use elmer_fe::ElmerQuadrature;
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::{element_nodes, lagrange_basis};

const CALLER: &str = "PoissonSolve";

/// Which `Poisson.f90` a case compiled. They differ only off the `n:1`
/// path (see the module note), so this records provenance and decides
/// nothing numerically.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PoissonSource {
    /// `1dtests/Poisson.f90`.
    Nodal,
    /// `ExtrudedMesh/Poisson.f90`, with static condensation of bubbles.
    Condensing,
}

pub struct PoissonSolver {
    source: PoissonSource,
    local: LocalSystem,
}

impl PoissonSolver {
    pub fn new(source: PoissonSource) -> Self {
        Self { source, local: LocalSystem::zeros(0) }
    }

    pub fn source(&self) -> PoissonSource {
        self.source
    }
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for the `1dtests` source.
pub fn nodal_factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    check_section(params)?;
    Ok(Box::new(PoissonSolver::new(PoissonSource::Nodal)))
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for the `ExtrudedMesh`
/// source.
pub fn condensing_factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    check_section(params)?;
    Ok(Box::new(PoissonSolver::new(PoissonSource::Condensing)))
}

/// Element definitions with more dofs than nodes (p-, bubble, edge or DG
/// elements) are not the `n:1` path ported here.
fn check_section(params: &ParamMap) -> Result<(), SolverError> {
    if let Some(e) = params.string("element") {
        if elmer_core::normalize_key(e) != "n:1" {
            return Err(SolverError::Unsupported(format!(
                "Poisson PoissonSolver: Element = {e} (only the nodal n:1 path is ported)"
            )));
        }
    }
    if params.logical("bubbles") == Some(true) {
        return Err(SolverError::Unsupported("Poisson PoissonSolver: Bubbles (LCondensate) is not ported".into()));
    }
    Ok(())
}

/// `LocalMatrix(STIFF, FORCE, LOAD, Element, n)`: over the element's Gauss
/// points, `STIFF += s*detJ * MATMUL(dBasisdx, TRANSPOSE(dBasisdx))` and
/// `FORCE += s*detJ * LoadAtIP * Basis`, products taken left to right as
/// the Fortran writes them.
pub fn local_matrix(
    etype: elmer_core::ElementType,
    xyz: &[[f64; 3]],
    load: &[f64],
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    local: &mut LocalSystem,
) -> Result<(), String> {
    let n = xyz.len();
    *local = LocalSystem::zeros(n);
    let mut ev = BasisEval::default();
    for qp in quad.gauss_points(etype, None) {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        if ev.basis.len() < n || ev.dbasis_dx.len() < n {
            return Err(format!("ElementInfo gave {} basis functions for {n} nodes", ev.basis.len()));
        }
        let load_at_ip = ev.basis[..n].iter().zip(load).fold(0.0, |s, (b, l)| s + b * l);
        let w = qp.weight * ev.det_j;
        let g = &ev.dbasis_dx;
        for p in 0..n {
            for q in 0..n {
                let m = g[p][0] * g[q][0] + g[p][1] * g[q][1] + g[p][2] * g[q][2];
                local.stiff[p * n + q] += w * m;
            }
            local.force[p] += w * load_at_ip * ev.basis[p];
        }
    }
    Ok(())
}

impl Solver for PoissonSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        crate::hooks::check_element_definitions(ctx)?;
        if ctx.primary().is_none_or(|v| v.dofs != 1) {
            return Err(SolverError::Unsupported("Poisson PoissonSolver: a scalar Variable (Variable DOFs = 1)".into()));
        }
        // `GetActiveElement` also walks boundary elements that a BC's `Body
        // ID` puts into a body (ExtrudedMeshBaseline1's baseline);
        // `active_elements` has bulk elements only, so the solve would be
        // empty.
        if let Some(i) = ctx.sim.model.boundary_conditions.iter().position(|bc| bc.get("body id").is_some()) {
            return Err(SolverError::Unsupported(format!(
                "Poisson PoissonSolver: Boundary Condition {} Body ID (boundary elements as a lower-dimensional body)",
                i + 1
            )));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let basis = lagrange_basis(ctx)?;
        ctx.default_initialize()?;
        let active: Vec<_> = ctx.active_elements().collect();
        for e in active {
            let etype = ctx.mesh.element(e).etype;
            let xyz = element_nodes(ctx.mesh, e);
            let n = xyz.len();
            // LOAD = 0; Load(1:n) = GetReal(BodyForce, 'Source', Found)
            let body_force = ctx.body_params(e).and_then(|b| b.body_force);
            let load = match body_force {
                Some(bf) => ctx.real_at_nodes(bf, "source", e)?.unwrap_or_else(|| vec![0.0; n]),
                None => vec![0.0; n],
            };
            local_matrix(etype, &xyz, &load, &basis, &ElmerQuadrature, &mut self.local)
                .map_err(|m| SolverError::fatal(CALLER, m))?;
            ctx.default_update(&self.local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}
