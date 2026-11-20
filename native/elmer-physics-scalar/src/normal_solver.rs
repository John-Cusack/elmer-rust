//! `NormalSolver` (fem/src/modules/NormalSolver.F90): the nodal normals of
//! the solver's (lower-dimensional) elements as the L2 projection of
//! `NormalVector`, one right-hand side per component on the same mass
//! matrix, then scaled to unit length at each node.
//!
//! The result goes to `Normals Result Variable` (`NormalSolver_Init` makes
//! it `Normals` when the case names none), with `DOFs = dim`.
//! Refused: `Constant Bulk Matrix`, `Point on Negative Side`, the vector
//! field projection (`Vector Field to Project`, `Project Displacement`,
//! `Project Velocity`), `Enforce Symmetry` and `Set Dirichlet Conditions`.

use elmer_core::{ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::common::{fe_provider, unsupported};

/// `Procedure = "NormalSolver" "NormalSolver"`.
pub const PROCEDURE: (&str, &str) = ("NormalSolver", "NormalSolver");

const CALLER: &str = "NormalSolver";

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| {
        let (basis, quad) = fe_provider()?();
        Ok(Box::new(NormalSolver { basis, quad }))
    });
}

/// `NormalSolver_Init`: the keywords it adds to the solver section, before
/// the variables are created.
pub fn solver_defaults(p: &mut ParamMap, dim: usize) {
    if p.get("variable").is_none() {
        p.insert("Variable", Value::Str("-nooutput nrm_temp".into()));
    }
    if p.get("normals result variable").is_none() {
        p.insert("Normals Result Variable", Value::Str("Normals".into()));
        let mut l = 1;
        while p.get(&format!("exported variable {l}")).is_some() {
            l += 1;
        }
        p.insert(&format!("Exported Variable {l}"), Value::Str(format!("Normals[Normals:{dim}]")));
    }
    p.insert("Time derivative order", Value::Integer(vec![0]));
}

struct NormalSolver {
    basis: Box<dyn ElementBasis + Send>,
    quad: Box<dyn Quadrature + Send>,
}

fn flag(p: &ParamMap, key: &str) -> bool {
    p.logical(key) == Some(true)
}

impl Solver for NormalSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "-------------------------------------", 4);
        ctx.info(CALLER, "Computing the normals", 4);
        ctx.info(CALLER, "-------------------------------------", 4);
        let p = ctx.params().clone();
        for k in ["constant bulk matrix", "project displacement", "project velocity", "enforce symmetry", "set dirichlet conditions"] {
            if flag(&p, k) {
                return Err(unsupported(format!("NormalSolver keyword \"{k}\""), "not ported"));
            }
        }
        if p.get("vector field to project").is_some() {
            return Err(unsupported("NormalSolver keyword \"Vector Field to Project\"", "not ported"));
        }
        let dim = ctx.mesh.dim();
        let Some(prim) = ctx.primary() else { return Ok(()) };
        let Some(perm) = prim.perm.map(|p| p.to_vec()) else { return Ok(()) };
        if perm.iter().all(|&k| k == NO_DOF) {
            return Ok(());
        }
        let n_rows = prim.values.len();
        let var_name = p
            .string("normals result variable")
            .ok_or_else(|| SolverError::fatal(CALLER, "> Normals Result Variable < not found!"))?
            .to_string();
        match ctx.variable(&var_name) {
            Some(v) if v.dofs != dim => return Err(SolverError::fatal(CALLER, "The normals should have DOFs equal to DIM")),
            Some(_) => {}
            None => return Err(unsupported(format!("Normals Result Variable \"{var_name}\" not exported"), "DefaultVariableAdd is not ported here")),
        }
        let coord = ctx.simulation_params().string("coordinate system").map(elmer_core::normalize_key).unwrap_or_default();
        let csym = coord == "axi symmetric" || coord == "cylindric symmetric";

        ctx.default_initialize()?;
        let mut forces = vec![vec![0.0f64; n_rows]; dim];
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut ev = elmer_core::BasisEval::default();
        let mut dofs = Vec::new();
        for &e in &elems {
            if let Some(bp) = ctx.body_params(e) {
                if bp.body.get("point on negative side").is_some() {
                    return Err(unsupported("Body keyword \"Point on Negative Side\"", "not ported"));
                }
            }
            let el = ctx.mesh.element(e);
            let etype = el.etype;
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let nd = xyz.len();
            let mut local = LocalSystem::zeros(nd);
            let mut force = vec![vec![0.0f64; nd]; dim];
            for qp in self.quad.gauss_points(etype, None) {
                self.basis
                    .eval(etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let b = &ev.basis[..nd];
                let mut weight = qp.weight * ev.det_j;
                if csym {
                    let mut r = 0.0;
                    for (bi, x) in b.iter().zip(&xyz) {
                        r += bi * x[0];
                    }
                    weight *= r;
                }
                for p in 0..nd {
                    for q in 0..nd {
                        local.stiff[p * nd + q] += weight * b[q] * b[p];
                    }
                }
                let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [qp.uvw[0], qp.uvw[1]])?;
                for i in 0..dim {
                    for q in 0..nd {
                        force[i][q] += weight * normal[i] * b[q];
                    }
                }
            }
            local.force.copy_from_slice(&force[0]);
            ctx.default_update(&local, e)?;
            // `DefaultUpdateForce( FORCE(i,1:nd) )` into `ForceVector(:,i)`.
            ctx.element_dofs(e, &mut dofs);
            for i in 0..dim {
                for (q, &d) in dofs.iter().enumerate() {
                    if d != NO_DOF {
                        forces[i][d as usize] += force[i][q];
                    }
                }
            }
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;

        let mut result = vec![0.0f64; n_rows * dim];
        for (i, f) in forces.iter().enumerate() {
            let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no linear system"))?;
            sys.rhs.copy_from_slice(f);
            ctx.default_dirichlet()?;
            {
                let v = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
                v.values.iter_mut().for_each(|x| *x = 0.0);
            }
            ctx.default_solve()?;
            let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
            for (r, x) in v.values.iter().enumerate() {
                result[dim * r + i] = *x;
            }
        }
        for r in 0..n_rows {
            let mut nrm = 0.0;
            for j in 0..dim {
                nrm += result[dim * r + j].powi(2);
            }
            if nrm > 0.0 {
                let nrm = 1.0 / nrm.sqrt();
                for j in 0..dim {
                    result[dim * r + j] *= nrm;
                }
            }
        }
        let nv = ctx.variable_mut(&var_name).ok_or_else(|| SolverError::fatal(CALLER, "no normals variable"))?;
        if nv.perm.map(|q| q.to_vec()) != Some(perm.clone()) || nv.values.len() != result.len() {
            return Err(unsupported("Normals Result Variable on another permutation", "not ported"));
        }
        nv.values.copy_from_slice(&result);
        Ok(())
    }
}
