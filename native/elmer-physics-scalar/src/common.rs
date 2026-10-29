//! Helpers shared by the two solvers: nodal interpolation, the coordinate
//! system switch and the finite-element provider the registry factories use.

use std::sync::OnceLock;

use elmer_core::{normalize_key, ElemId, ElementBasis, ParamMap, Params, Quadrature};
use elmer_solver::{SolverCtx, SolverError};

use crate::error::ScalarError;

/// Builds the element basis and quadrature the solvers integrate with. The
/// registry's factories are plain `fn` pointers, so the front end (`elmer-fe`)
/// hands its implementations over through [`install_fe`] instead.
pub type FeProvider = fn() -> (Box<dyn ElementBasis + Send>, Box<dyn Quadrature + Send>);

static FE: OnceLock<FeProvider> = OnceLock::new();

/// Installs the basis and quadrature provider used by the factories of both
/// solvers. Returns `false` if one is already installed.
pub fn install_fe(provider: FeProvider) -> bool {
    FE.set(provider).is_ok()
}

pub(crate) fn fe_provider() -> Result<FeProvider, SolverError> {
    FE.get()
        .copied()
        .ok_or_else(|| SolverError::NotWired("element basis and quadrature (elmer-fe)".into()))
}

impl From<ScalarError> for SolverError {
    fn from(e: ScalarError) -> Self {
        match e {
            ScalarError::Unsupported { .. } => SolverError::Unsupported(e.to_string()),
            e => SolverError::fatal("elmer-physics-scalar", e.to_string()),
        }
    }
}

pub(crate) fn unsupported(what: impl Into<String>, reason: &'static str) -> SolverError {
    ScalarError::Unsupported {
        what: what.into(),
        reason,
    }
    .into()
}

/// Value at an integration point from nodal values and basis values.
pub(crate) fn interp(phi: &[f64], nodal: &[f64]) -> f64 {
    phi.iter().zip(nodal).map(|(b, v)| b * v).sum()
}

pub(crate) fn check_len(what: &'static str, n: usize, got: usize) -> Result<(), ScalarError> {
    if got == n {
        Ok(())
    } else {
        Err(ScalarError::NodalLength {
            what,
            expected: n,
            got,
        })
    }
}

pub(crate) fn check_opt(
    what: &'static str,
    n: usize,
    v: &Option<Vec<f64>>,
) -> Result<(), ScalarError> {
    v.as_ref().map_or(Ok(()), |v| check_len(what, n, v.len()))
}

/// `CurrentCoordinateSystem() /= Cartesian`: axisymmetric (`r` first, `z`
/// second) or cylindric symmetric. Other systems are refused.
pub(crate) fn axisymmetric(ctx: &SolverCtx<'_>) -> Result<bool, SolverError> {
    match ctx.simulation_params().string("coordinate system") {
        None => Ok(false),
        Some(s) => {
            let s = normalize_key(s);
            if s.starts_with("cartesian") {
                Ok(false)
            } else if s == "axi symmetric" || s == "cylindric symmetric" {
                Ok(true)
            } else {
                Err(SolverError::Unsupported(format!("Coordinate System = {s}")))
            }
        }
    }
}

/// `ListGetReal(section, key, n, Indexes)` for an optional section: `None`
/// when the section or the keyword is absent.
pub(crate) fn real_at_nodes(
    ctx: &SolverCtx<'_>,
    p: Option<&ParamMap>,
    key: &str,
    e: ElemId,
) -> Result<Option<Vec<f64>>, SolverError> {
    match p {
        Some(p) => ctx.real_at_nodes(p, key, e),
        None => Ok(None),
    }
}

/// The running solver's slot, for the keywords `_init` adds to its section.
pub(crate) fn slot_mut<'a>(ctx: &'a mut SolverCtx<'_>) -> &'a mut elmer_solver::SolverSlot {
    let id = ctx.solver_id();
    ctx.sim
        .slots
        .iter_mut()
        .find(|s| s.id == id)
        .expect("the running solver has a slot")
}

/// A material coefficient at the nodes of element `e`, by the shape of the
/// keyword's array (`ListGetElementReal`'s `Rdim`): 1 row is a scalar, 1
/// column a diagonal, else a full tensor. A scalar may depend on fields; a
/// tensor must be constant. `None` when the keyword is absent.
pub(crate) fn coefficient_at_nodes(
    ctx: &SolverCtx<'_>,
    p: Option<&ParamMap>,
    key: &str,
    e: ElemId,
    n: usize,
    dim: usize,
) -> Result<Option<crate::potential::Coefficient>, SolverError> {
    use crate::potential::Coefficient;
    use elmer_core::Value;
    let Some(p) = p else { return Ok(None) };
    let Some(v) = p.get(key) else { return Ok(None) };
    match v {
        Value::Real { rows, cols, data } if rows * cols > 1 => {
            let (r, c) = (*rows, *cols);
            if c == 1 {
                if r < dim {
                    return Err(SolverError::fatal(
                        "ListGetElementReal",
                        format!("\"{key}\" has {r} rows, the problem has {dim} dimensions"),
                    ));
                }
                let mut d = [0.0; 3];
                d[..dim].copy_from_slice(&data[..dim]);
                Ok(Some(Coefficient::Diagonal(vec![d; n])))
            } else {
                if r < dim || c < dim {
                    return Err(SolverError::fatal(
                        "ListGetElementReal",
                        format!("\"{key}\" is {r}x{c}, the problem has {dim} dimensions"),
                    ));
                }
                let mut k = [[0.0; 3]; 3];
                for (i, row) in k.iter_mut().enumerate().take(dim) {
                    for (j, kij) in row.iter_mut().enumerate().take(dim) {
                        *kij = data[i * c + j];
                    }
                }
                Ok(Some(Coefficient::Full(vec![k; n])))
            }
        }
        Value::Dependent(d) if d.rows * d.cols > 1 => Err(SolverError::Unsupported(format!(
            "\"{key}\": a non-constant tensor"
        ))),
        _ => Ok(ctx.real_at_nodes(p, key, e)?.map(Coefficient::Scalar)),
    }
}

/// `Initialize Dirichlet Conditions` (ElmerSolver.F90, `SetInitialConditions`):
/// unless `Simulation` turns it off, the solver variable starts at the
/// Dirichlet values on its boundary nodes, so the first time step of a
/// transient run sees them in the history. `elmer-solver` does not do this
/// for any solver yet.
pub(crate) fn init_dirichlet(ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
    if ctx
        .simulation_params()
        .logical("initialize dirichlet conditions")
        == Some(false)
    {
        return Ok(());
    }
    let fixed = ctx.dirichlet_values()?;
    if let Some(v) = ctx.primary_mut() {
        for (dof, x) in fixed {
            v.values[dof as usize] = x;
        }
    }
    Ok(())
}
