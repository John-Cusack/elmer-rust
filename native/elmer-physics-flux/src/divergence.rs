//! `DivergenceSolver` (DivergenceSolver.F90): the Galerkin projection of
//! `div v` of a vector field `v` onto the nodal basis, `M d = (N, div v)`.
//!
//! Ported: `Divergence Variable` / `Target Variable` (default `Velocity`),
//! `Relative Divergence` and the axisymmetric `v_r / r` term. Refused:
//! `Constant Bulk Matrix`.

use elmer_core::{normalize_key, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

/// `Procedure = "DivergenceSolver" "DivergenceSolver"`.
pub const DIVERGENCE_PROCEDURE: (&str, &str) = ("DivergenceSolver", "DivergenceSolver");

const CALLER: &str = "DivergenceSolver";

pub(crate) fn register(registry: &mut SolverRegistry) {
    registry.register(
        DIVERGENCE_PROCEDURE.0,
        DIVERGENCE_PROCEDURE.1,
        divergence_factory,
    );
}

pub fn divergence_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(DivergenceSolver))
}

fn add_new(p: &mut ParamMap, key: &str, v: Value) {
    if p.get(key).is_none() {
        p.insert(key, v);
    }
}

/// `DivergenceSolver_Init`: the `Div <var>` (or `Divergence`) variable and
/// the default linear system (CG, ILU0, 500 iterations, 1e-8). Call it on
/// the solver section before the variables are created.
pub fn divergence_defaults(p: &mut ParamMap) {
    let name = match p.string("divergence variable") {
        Some(v) => format!("Div {}", v.trim()),
        None => "Divergence".to_string(),
    };
    if p.get("variable").is_none() {
        p.insert("Variable DOFs", Value::Integer(vec![1]));
        p.insert("Variable", Value::Str(name));
    }
    p.insert("Time derivative order", Value::Integer(vec![0]));
    add_new(p, "Linear System Solver", Value::Str("Iterative".into()));
    add_new(p, "Linear System Iterative Method", Value::Str("cg".into()));
    add_new(
        p,
        "Linear System Preconditioning",
        Value::Str("ILU0".into()),
    );
    add_new(p, "Linear System Max Iterations", Value::Integer(vec![500]));
    add_new(
        p,
        "Linear System Convergence Tolerance",
        Value::Real {
            rows: 1,
            cols: 1,
            data: vec![1.0e-8],
        },
    );
}

pub struct DivergenceSolver;

/// `GetScalarLocalSolution( x, name )`: the field at the element's nodes,
/// zero where the field or its permutation entry is missing.
fn local_component(ctx: &SolverCtx<'_>, name: &str, nodes: &[u32]) -> Vec<f64> {
    let mut x = vec![0.0; nodes.len()];
    let Some(f) = ctx.component(name) else {
        return x;
    };
    for (xi, &n) in x.iter_mut().zip(nodes) {
        let slot = match f.perm {
            Some(p) => match p.get(n as usize) {
                Some(&s) if s != NO_DOF => s as usize,
                _ => continue,
            },
            None => n as usize,
        };
        if let Some(v) = f.values.get(slot * f.stride + f.offset) {
            *xi = *v;
        }
    }
    x
}

impl Solver for DivergenceSolver {
    fn init(&mut self, _ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "-------------------------------------", 4);
        ctx.info(CALLER, "Computing the divergence field", 4);
        ctx.info(CALLER, "-------------------------------------", 4);
        let Some(var) = ctx.primary() else {
            return Err(SolverError::fatal(CALLER, "Variable does not exist!"));
        };
        if var.perm.is_some_and(|p| p.iter().all(|&s| s == NO_DOF)) {
            ctx.warn(CALLER, "Size is zero, nothing to compute");
            return Ok(());
        }
        if var.dofs != 1 {
            return Err(SolverError::fatal(
                CALLER,
                "Divergence should have 1 component",
            ));
        }
        let p = ctx.params().clone();
        if p.logical("constant bulk matrix") == Some(true) {
            return Err(SolverError::Unsupported(
                "DivergenceSolver: Constant Bulk Matrix".into(),
            ));
        }
        let cs = ctx
            .simulation_params()
            .string("coordinate system")
            .map(normalize_key)
            .unwrap_or_else(|| "cartesian".into());
        let csymmetry = cs.starts_with("axi") || cs.starts_with("cylindric symmetric");
        let dim = if csymmetry || cs.contains("2d") {
            2
        } else if cs.contains("3d") {
            3
        } else {
            ctx.mesh.dim()
        };
        let var_name = p
            .string("divergence variable")
            .or_else(|| p.string("target variable"))
            .unwrap_or("Velocity")
            .trim()
            .to_string();
        let relative = p.logical("relative divergence") == Some(true);

        ctx.default_initialize()?;
        if ctx.sim.backends.fe.is_none() {
            return Err(SolverError::NotWired(
                "basis and quadrature (elmer-fe)".into(),
            ));
        }
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut ev = elmer_core::BasisEval::default();
        for &e in &elems {
            let el = ctx.mesh.element(e);
            let nodes = el.nodes.to_vec();
            let etype = el.etype;
            let n = nodes.len();
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let comp = |k: usize| format!("{var_name} {k}");
            let vx = local_component(ctx, &comp(1), &nodes);
            let vy = local_component(ctx, &comp(2), &nodes);
            let vz = if dim == 3 {
                local_component(ctx, &comp(3), &nodes)
            } else {
                vec![0.0; n]
            };
            let mut vabs = 0.0;
            if relative {
                let mut velo = [0.0f64; 3];
                velo[0] = vx.iter().sum::<f64>() / n as f64;
                velo[1] = vy.iter().sum::<f64>() / n as f64;
                if dim == 3 {
                    velo[2] = vz.iter().sum::<f64>() / n as f64;
                }
                let s: f64 = velo[..dim].iter().map(|v| v * v).sum();
                vabs = s.sqrt().max(f64::EPSILON);
            }
            let mut ls = LocalSystem::zeros(n);
            let fe = ctx.sim.backends.fe.as_ref().expect("checked above");
            for qp in fe.quadrature.gauss_points(etype, None) {
                fe.basis
                    .eval(etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let mut weight = qp.weight * ev.det_j;
                let mut x = 0.0;
                if csymmetry {
                    x = sum_prod(&ev.basis, xyz.iter().map(|c| c[0]));
                    weight *= x;
                }
                for p in 0..n {
                    for q in 0..n {
                        ls.stiff[p * n + q] += weight * ev.basis[q] * ev.basis[p];
                    }
                }
                let mut source = sum_prod_d(&ev.dbasis_dx, 0, &vx);
                source += sum_prod_d(&ev.dbasis_dx, 1, &vy);
                if dim == 3 {
                    source += sum_prod_d(&ev.dbasis_dx, 2, &vz);
                }
                if csymmetry {
                    source += sum_prod(&ev.basis, vx.iter().copied()) / x;
                }
                if relative {
                    source /= vabs;
                }
                for p in 0..n {
                    ls.force[p] += ev.basis[p] * weight * source;
                }
            }
            ctx.default_update(&ls, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_dirichlet()?;
        let norm = ctx.default_solve()?;
        ctx.info(CALLER, &format!(" Result Norm: {norm}"), 4);
        Ok(())
    }
}

/// `SUM( a(1:n) * b(1:n) )`, in order.
fn sum_prod(a: &[f64], b: impl Iterator<Item = f64>) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// `SUM( dBasisdx(1:n, i) * v(1:n) )`.
fn sum_prod_d(db: &[[f64; 3]], i: usize, v: &[f64]) -> f64 {
    let mut s = 0.0;
    for (d, x) in db.iter().zip(v) {
        s += d[i] * x;
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_names_the_divergence_variable() {
        let mut p = ParamMap::new();
        p.insert("Divergence Variable", Value::Str("CoilCurrent".into()));
        p.insert("Linear System Iterative Method", Value::Str("BiCGStab".into()));
        divergence_defaults(&mut p);
        assert_eq!(p.string("variable"), Some("Div CoilCurrent"));
        assert_eq!(p.integer("variable dofs"), Some(1));
        // ListAddNew: the section's own choice stays.
        assert_eq!(p.string("linear system iterative method"), Some("BiCGStab"));
        assert_eq!(p.string("linear system preconditioning"), Some("ILU0"));
        let mut q = ParamMap::new();
        divergence_defaults(&mut q);
        assert_eq!(q.string("variable"), Some("Divergence"));
    }
}
