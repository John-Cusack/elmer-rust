//! `OdeSolver` (fem/src/modules/OdeSolver.F90): an ordinary differential
//! equation system `M u'' + D u' + K u = f` (or `D u' + K u = f` with `Time
//! Derivative Order = 1`) on a global variable of `Variable Dofs` unknowns,
//! its coefficients `Stiffness Matrix ij`, `Damping Matrix ij`, `Mass Matrix
//! ij`, `Force i` from the `Active Components` component (or the solver
//! section), solved with Elmer's global-mass-matrix time integration
//! (`Default1stOrderTimeGlobal` / `Default2ndOrderTimeGlobal`).

use elmer_core::{ParamMap, Params, Value};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

pub const PROCEDURE: (&str, &str) = ("OdeSolver", "OdeSolver");
const CALLER: &str = "OdeSolver";

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(OdeSolver))
}

/// Whether `params` is an OdeSolver section.
pub fn is_section(params: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(params)
        .is_some_and(|(lib, entry)| lib.eq_ignore_ascii_case(PROCEDURE.0) && entry.eq_ignore_ascii_case(PROCEDURE.1))
}

/// `OdeSolver_init`, before the solver's variable is created.
pub fn init_section(params: &mut ParamMap) {
    if params.get("time derivative order").is_none() {
        params.insert("Time derivative order", Value::Integer(vec![2]));
    }
    params.insert("Use Global Mass Matrix", Value::Logical(true));
    params.insert("Ode Matrix", Value::Logical(true));
    params.insert("Variable Global", Value::Logical(true));
}

struct OdeSolver;

impl Solver for OdeSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "Solving ordinary differential equation", 4);
        let max_iter = ctx.params().integer("nonlinear system max iterations").unwrap_or(1).max(0);
        let time_order = ctx.params().integer("time derivative order").unwrap_or(0);
        for iter in 1..=max_iter {
            if max_iter > 1 {
                ctx.info(CALLER, &format!("Nonlinear iteration: {iter}"), 5);
            }
            ctx.default_initialize()?;
            ode_matrix(ctx, time_order)?;
            ctx.default_finish_bulk_assembly()?;
            ctx.default_finish_boundary_assembly()?;
            ctx.default_solve()?;
            if ctx.convergence().is_some_and(|c| c.nonlin_converged == 1) {
                break;
            }
        }
        ctx.info(CALLER, "All done", 5);
        Ok(())
    }
}

/// `OdeMatrix`: `CRS_SetMatrixElement` of every coefficient, from the
/// active component or the solver section (`ListGetCReal`, 0 when absent).
fn ode_matrix(ctx: &mut SolverCtx<'_>, time_order: i32) -> Result<(), SolverError> {
    let list: ParamMap = match ctx.params().get("active components") {
        Some(Value::Integer(ids)) => {
            if ids.len() > 1 {
                return Err(SolverError::fatal(CALLER, "Currently implemented only for one component!"));
            }
            let i = ids[0];
            if i < 1 || i as usize > ctx.sim.components.len() {
                return Err(SolverError::fatal(CALLER, "Active Component index out of range"));
            }
            ctx.info(CALLER, &format!("Using active component: {i}"), 10);
            ctx.sim.components[i as usize - 1].clone()
        }
        _ => {
            ctx.info(CALLER, "Using solver section", 10);
            ctx.params().clone()
        }
    };
    let dofs = ctx.primary().map(|v| v.dofs).ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    let get = |key: String| -> Result<f64, SolverError> { Ok(ctx.creal(&list, &key)?.unwrap_or(0.0)) };
    let second = time_order >= 2;
    let mut stiff = vec![0.0; dofs * dofs];
    let mut mass = vec![0.0; dofs * dofs];
    let mut damp = vec![0.0; dofs * dofs];
    let mut force = vec![0.0; dofs];
    for i in 0..dofs {
        for j in 0..dofs {
            stiff[i * dofs + j] = get(format!("stiffness matrix {}{}", i + 1, j + 1))?;
            // The first-order term of a first-order system is `Damping Matrix`.
            let m = if second { "mass matrix" } else { "damping matrix" };
            mass[i * dofs + j] = get(format!("{m} {}{}", i + 1, j + 1))?;
            if second {
                damp[i * dofs + j] = get(format!("damping matrix {}{}", i + 1, j + 1))?;
            }
        }
        force[i] = get(format!("force {}", i + 1))?;
    }
    let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no matrix (call default_initialize)"))?;
    set_dense(&mut sys.matrix, dofs, &stiff);
    if let Some(m) = sys.mass.as_mut() {
        set_dense(m, dofs, &mass);
    }
    if second {
        if let Some(d) = sys.damp.as_mut() {
            set_dense(d, dofs, &damp);
        }
    }
    sys.rhs[..dofs].copy_from_slice(&force);
    Ok(())
}

/// `CRS_SetMatrixElement (A, i, j, v)` for every `(i, j)` of the dense ODE
/// matrix (CreateOdeMatrix's full pattern).
fn set_dense(a: &mut elmer_core::CsrMatrix<f64>, dofs: usize, v: &[f64]) {
    for i in 0..dofs {
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            let j = a.col_idx[k] as usize;
            a.values[k] = v[i * dofs + j];
        }
    }
}
