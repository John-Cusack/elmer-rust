//! `fem/tests/BlockLinElast1/BlockLinElast.f90`'s `BulkAssembly`, the
//! case's own element routine for a `Block Solver` of linear elasticity
//! split by displacement component (strategy D12: a test-local plugin, not
//! product code). Block `(i, j)` (`Block Matrix Row` / `Column`) gets
//! `lambda dN_q/dx_j dN_p/dx_i + mu dN_q/dx_i dN_p/dx_j`, plus on the
//! diagonal `mu grad N_q . grad N_p` and the source `Source N_p`, on the
//! element's p basis (`GetElementNOFDOFs` functions).
//!
//! `Source` is read in the block's namespace first (`Block ij: Source`,
//! the `ListPushNamespace` of `BlockSystemAssembly`), then plainly.

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params};
use elmer_solver::{ElementAssembly, SolverCtx, SolverError};

use crate::hooks::{element, nodal};

const CALLER: &str = "BlockLinElast BulkAssembly";

#[derive(Default)]
pub struct BlockLinElastBulk;

/// `GetCReal(Material, key)`: the constant value, 0 when absent.
fn constant(material: Option<&ParamMap>, key: &str) -> Result<f64, SolverError> {
    match material.and_then(|m| m.get(key)) {
        None => Ok(0.0),
        Some(_) => material
            .and_then(|m| m.real(key))
            .ok_or_else(|| SolverError::Unsupported(format!("{CALLER}: a non-constant {key}"))),
    }
}

impl ElementAssembly for BlockLinElastBulk {
    fn assemble(
        &mut self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        _dt: f64,
        _transient: bool,
        local: &mut LocalSystem,
    ) -> Result<(), SolverError> {
        let p = ctx.params();
        let i = p.integer("block matrix row").unwrap_or(1) as usize;
        let j = p.integer("block matrix column").unwrap_or(1) as usize;
        let el = element(ctx, e)?;
        let (n, nd) = (el.nodes.len(), el.nd);
        if local.force.len() != nd {
            return Err(SolverError::fatal(CALLER, "local system and element dofs differ"));
        }
        let material = ctx.body_params(e).and_then(|b| b.material);
        let block_key = format!("block {i}{j}: source");
        let key = if material.is_some_and(|m| m.get(&block_key).is_some()) { block_key.as_str() } else { "source" };
        let source = nodal(ctx, material, key, e, &el.nodes)?;
        let young = constant(material, "youngs modulus")?;
        let poisson = constant(material, "poisson ratio")?;
        let lame1 = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
        let lame2 = young / (2.0 * (1.0 + poisson));
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
        let mut be = BasisEval::default();
        for qp in &el.rule {
            el.basis
                .eval(el.etype, &xyz, qp.uvw, &mut be)
                .map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let w = qp.weight * be.det_j;
            let (b, g) = (&be.basis[..nd], &be.dbasis_dx[..nd]);
            let s_ip = (0..n).fold(0.0, |acc, k| acc + b[k] * source[k]);
            if i == j {
                for pp in 0..nd {
                    for q in 0..nd {
                        // MATMUL(dBasisdx, TRANSPOSE(dBasisdx)) over the three columns.
                        let mut d = 0.0;
                        for c in 0..3 {
                            d += g[pp][c] * g[q][c];
                        }
                        local.stiff[pp * nd + q] += w * lame2 * d;
                    }
                    local.force[pp] += w * s_ip * b[pp];
                }
            }
            for pp in 0..nd {
                for q in 0..nd {
                    let a = &mut local.stiff[pp * nd + q];
                    *a += w * lame1 * g[q][j - 1] * g[pp][i - 1];
                    *a += w * lame2 * g[q][i - 1] * g[pp][j - 1];
                }
            }
        }
        Ok(())
    }
}
