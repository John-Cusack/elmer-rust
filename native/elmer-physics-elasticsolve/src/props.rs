//! Keyword evaluation for the solver: constants, tables, MATC and `Procedure`
//! values through `elmer-expr`, at the nodes of an element (`ListGetReal`,
//! `ListGetRealArray`) or at its integration points (`ListGetElementReal`).

use elmer_core::{ElemId, ParamMap, Params, Value};
use elmer_expr::{compile, MatcEngine};
use elmer_solver::expr::NoMatc;
use elmer_solver::{SolverCtx, SolverError};

const CALLER: &str = "ElasticSolver";

/// A keyword's value at each node of an element: `rows x cols` numbers per
/// node, node-major and row-major within a node.
#[derive(Clone, Debug)]
pub(crate) struct NodeArray {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
}

impl NodeArray {
    pub fn ncomp(&self) -> usize {
        self.rows * self.cols
    }

    /// The `rows x cols` block of node `k`.
    pub fn node(&self, k: usize) -> &[f64] {
        &self.data[k * self.ncomp()..(k + 1) * self.ncomp()]
    }
}

fn with_matc<R>(ctx: &SolverCtx<'_>, f: impl FnOnce(&dyn MatcEngine) -> R) -> R {
    let none = NoMatc;
    match ctx.sim.backends.matc.as_deref() {
        Some(m) => f(m),
        None => f(&none),
    }
}

fn not_real(key: &str) -> SolverError {
    SolverError::fatal(CALLER, format!("'{key}' is not a real"))
}

/// `ListGetRealArray` at the nodes of `e`; `None` when the section lacks `key`.
pub(crate) fn node_array(
    ctx: &SolverCtx<'_>,
    section: &ParamMap,
    key: &str,
    e: ElemId,
) -> Result<Option<NodeArray>, SolverError> {
    let nodes = ctx.mesh.element(e).nodes;
    match section.get(key) {
        None => Ok(None),
        Some(Value::Real { rows, cols, data }) => {
            // `elmer-sif` reports `Name(r,c)` with `rows = r*c`; the data length is what counts.
            let rows = if *cols > 0 && rows * cols != data.len() && data.len() % cols == 0 { data.len() / cols } else { *rows };
            Ok(Some(NodeArray {
                rows,
                cols: *cols,
                data: data.iter().copied().cycle().take(data.len() * nodes.len()).collect(),
            }))
        }
        Some(v @ Value::Dependent(_)) => {
            let c = compile(v);
            let mut data = vec![0.0; nodes.len() * c.ncomp()];
            with_matc(ctx, |m| c.eval_nodes(nodes, &ctx.fields(), m, &mut data))?;
            Ok(Some(NodeArray { rows: c.rows(), cols: c.cols(), data }))
        }
        Some(_) => Err(not_real(key)),
    }
}

/// `ListGetReal` at the nodes of `e`.
pub(crate) fn scalar_nodes(
    ctx: &SolverCtx<'_>,
    section: &ParamMap,
    key: &str,
    e: ElemId,
) -> Result<Option<Vec<f64>>, SolverError> {
    match node_array(ctx, section, key, e)? {
        None => Ok(None),
        Some(a) if a.ncomp() == 1 => Ok(Some(a.data)),
        Some(_) => Err(SolverError::fatal(CALLER, format!("'{key}' must be a scalar"))),
    }
}

/// `ListGetElementReal` at each integration point of `e`. `shape` holds the
/// shape function values at the points (`n_ip` rows of `n` values).
pub(crate) fn scalar_ips(
    ctx: &SolverCtx<'_>,
    section: &ParamMap,
    key: &str,
    e: ElemId,
    shape: &[f64],
) -> Result<Option<Vec<f64>>, SolverError> {
    let Some(v) = section.get(key) else { return Ok(None) };
    if !matches!(v, Value::Real { .. } | Value::Dependent(_)) {
        return Err(not_real(key));
    }
    let nodes = ctx.mesh.element(e).nodes;
    let c = compile(v);
    if c.ncomp() != 1 {
        return Err(SolverError::fatal(CALLER, format!("'{key}' must be a scalar")));
    }
    let n_ip = shape.len() / nodes.len();
    let mut out = vec![0.0; n_ip];
    with_matc(ctx, |m| c.eval_ips(nodes, shape, &ctx.fields(), m, &mut out))?;
    Ok(Some(out))
}

/// `GetCReal`: a keyword that must not vary in space.
pub(crate) fn const_real(ctx: &SolverCtx<'_>, section: &ParamMap, key: &str) -> Result<Option<f64>, SolverError> {
    let Some(v) = section.get(key) else { return Ok(None) };
    with_matc(ctx, |m| elmer_solver::expr::eval_const(v, Some(m))).map(Some)
}
