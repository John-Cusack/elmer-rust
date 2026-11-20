//! `OdeSolver` (fem/src/modules/OdeSolver.F90): a system of ordinary
//! differential equations on a global variable, its matrices and force
//! given as constants in the solver section or in the `Active Components`'
//! `Component` section:
//!
//! ```text
//!   D du/dt + K u = f                  (Time Derivative Order = 1)
//!   M d2u/dt2 + D du/dt + K u = f      (Time Derivative Order = 2)
//! ```
//!
//! `OdeSolver_init` asks for a full `Dofs x Dofs` system (`Ode Matrix`) with
//! the time terms added to the assembled matrices (`Use Global Mass
//! Matrix`). Refused: more than one active component.

use elmer_core::{ParamMap, Params, Value};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

/// `Procedure = "OdeSolver" "OdeSolver"`.
pub const PROCEDURE: (&str, &str) = ("OdeSolver", "OdeSolver");

const CALLER: &str = "OdeSolver";

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| Ok(Box::new(OdeSolver)));
}

/// `OdeSolver_init`: the keywords it adds to the solver section, before
/// the variable and the matrix are created.
pub fn solver_defaults(p: &mut ParamMap) {
    if p.get("time derivative order").is_none() {
        p.insert("Time derivative order", Value::Integer(vec![2]));
    }
    p.insert("Use Global Mass Matrix", Value::Logical(true));
    p.insert("Ode Matrix", Value::Logical(true));
    p.insert("Variable Global", Value::Logical(true));
}

struct OdeSolver;

/// `ListGetCReal(OdeList, key, Found)`: 0 when absent.
fn creal(list: &ParamMap, key: &str) -> Result<f64, SolverError> {
    match list.get(key) {
        None => Ok(0.0),
        Some(Value::Real { data, .. }) if data.len() == 1 => Ok(data[0]),
        Some(Value::Integer(v)) if v.len() == 1 => Ok(f64::from(v[0])),
        Some(v) => Err(SolverError::Unsupported(format!("OdeSolver: \"{key}\" is not a constant ({v:?})"))),
    }
}

impl Solver for OdeSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "Solving ordinary differential equation", 4);
        let max_iter = ctx.params().integer("nonlinear system max iterations").unwrap_or(1).max(1);
        let order = ctx.params().integer("time derivative order").unwrap_or(2);
        // OdeMatrix: the section holding the matrices.
        let list = match ctx.params().get("active components") {
            Some(Value::Integer(v)) => {
                if v.len() > 1 {
                    return Err(SolverError::fatal(CALLER, "Currently implemented only for one component!"));
                }
                let i = v[0];
                let comps = &ctx.sim.model.components;
                if i < 1 || i as usize > comps.len() {
                    return Err(SolverError::fatal(CALLER, "Active Component index out of range"));
                }
                comps[i as usize - 1].clone()
            }
            _ => ctx.params().clone(),
        };
        let dofs = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?.dofs;
        let mut stiff = vec![0.0; dofs * dofs];
        let mut damp = vec![0.0; dofs * dofs];
        let mut mass = vec![0.0; dofs * dofs];
        let mut force = vec![0.0; dofs];
        for i in 0..dofs {
            for j in 0..dofs {
                stiff[i * dofs + j] = creal(&list, &format!("stiffness matrix {}{}", i + 1, j + 1))?;
                damp[i * dofs + j] = creal(&list, &format!("damping matrix {}{}", i + 1, j + 1))?;
                mass[i * dofs + j] = creal(&list, &format!("mass matrix {}{}", i + 1, j + 1))?;
            }
            force[i] = creal(&list, &format!("force {}", i + 1))?;
        }
        for iter in 1..=max_iter {
            if max_iter > 1 {
                ctx.info(CALLER, &format!("Nonlinear iteration: {iter}"), 5);
            }
            ctx.default_initialize()?;
            {
                let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no linear system"))?;
                // CRS_SetMatrixElement on the full pattern.
                sys.matrix.values.copy_from_slice(&stiff);
                // MassValues hold the damping of a first-order system, the
                // mass of a second-order one (with DampValues its damping).
                let m = sys.mass.as_mut().ok_or_else(|| SolverError::fatal(CALLER, "no global mass matrix"))?;
                if order >= 2 {
                    m.values.copy_from_slice(&mass);
                    let d = sys.damp.as_mut().ok_or_else(|| SolverError::fatal(CALLER, "no global damping matrix"))?;
                    d.values.copy_from_slice(&damp);
                } else {
                    m.values.copy_from_slice(&damp);
                }
                sys.rhs.copy_from_slice(&force);
            }
            ctx.default_finish_bulk_assembly()?;
            ctx.default_finish_boundary_assembly()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        ctx.info(CALLER, "All done", 5);
        Ok(())
    }
}
