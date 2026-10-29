//! `elmer-fixtures`' `PoissonAssembly` behind the solver's assembly hook.

use elmer_core::{ElemId, LocalSystem};
use elmer_fixtures::poisson_assembly::{BulkProps, Element, PoissonAssembly};
use elmer_solver::{ElementAssembly, SolverCtx, SolverError};

/// `BulkAssembly` as the routine sees `DefUtils`: the element's material
/// values at its nodes (`GetReal`) and its coordinates (`GetElementNodes`).
#[derive(Default)]
pub struct Bulk(PoissonAssembly);

impl ElementAssembly for Bulk {
    fn assemble(
        &mut self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        _dt: f64,
        _transient: bool,
        local: &mut LocalSystem,
    ) -> Result<(), SolverError> {
        let el = ctx.mesh.element(e);
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
        let mat = ctx.body_params(e).and_then(|b| b.material).expect("material");
        let get = |k: &str| ctx.real_at_nodes(mat, k, e);
        let (cond, source, reaction) = (get("conductivity")?, get("source")?, get("reaction")?);
        let props = BulkProps {
            conductivity: cond.as_deref(),
            source: source.as_deref(),
            reaction: reaction.as_deref(),
        };
        self.0
            .bulk(
                Element { etype: el.etype, nodes: &xyz },
                props,
                &super::basis::TestBasis,
                &super::quad::TestQuadrature,
                ctx.log,
                local,
            )
            .map_err(|err| SolverError::fatal("PoissonBulkAssembly", err.to_string()))
    }
}
