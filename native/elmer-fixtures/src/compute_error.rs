//! `ComputeError.F90` of the `RaviartThomas2D` and `BDM2D` tests
//! (`Procedure = "ComputeError" "IntegrateError"`): the relative L2 errors
//! of the potential and the flux of the `Target Equation`'s mixed Poisson
//! solution against `u = 16 (x - x^2)(y - y^2)`, printed as
//! `L2 Scalar Error = ...` and `L2 Flux Error = ...`.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{ElemId, Params, Value};
use elmer_fe_hp::face::face_element_info_triangle;
use elmer_pdofs::quadrature::gauss_points_triangle_p;
use elmer_pdofs::PDofsBackend;
use elmer_solver::{DofsBackend, Solver, SolverCtx, SolverError};

use crate::poisson_ml::fortran_e;

const CALLER: &str = "ComputeError";

pub fn factory(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(IntegrateError))
}

struct IntegrateError;

impl Solver for IntegrateError {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let target = ctx
            .params()
            .string("target equation")
            .ok_or_else(|| SolverError::fatal("ListGetString", "Requested property: [Target Equation], not found"))?
            .trim()
            .to_string();
        let model = &ctx.sim.model;
        let tid = model
            .solvers
            .iter()
            .position(|s| s.string("equation").is_some_and(|e| e.trim() == target))
            .ok_or_else(|| SolverError::fatal(CALLER, "Target Equation does not exist"))?;
        let tp = model.solvers[tid].clone();
        let second = tp.logical("second kind basis").unwrap_or(false);
        // The target's active elements: bulk elements whose equation lists it.
        let active: Vec<ElemId> = (0..ctx.mesh.n_bulk() as u32)
            .map(ElemId)
            .filter(|&e| {
                model
                    .body_of_tag(ctx.mesh.element(e).tag)
                    .and_then(|b| model.body_params(b))
                    .and_then(|b| b.equation)
                    .and_then(|eq| match eq.get("active solvers") {
                        Some(Value::Integer(v)) => Some(v.contains(&(tid as i32 + 1))),
                        _ => None,
                    })
                    .unwrap_or(false)
            })
            .collect();
        let layout = PDofsBackend::new()
            .for_solver(ctx.mesh, model, &tp, &active, 1)
            .ok_or_else(|| SolverError::fatal(CALLER, "the target has no face element dofs"))??;
        let var = tp.string("variable").unwrap_or("").to_ascii_lowercase();
        let sol = ctx
            .variable(&var)
            .ok_or_else(|| SolverError::fatal(CALLER, format!("no variable {var}")))?
            .values
            .to_vec();
        let rule = gauss_points_triangle_p(3).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
        let (mut pres_err, mut flux_err, mut sol_norm, mut flux_norm) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut slots = Vec::new();
        for &e in &active {
            let el = ctx.mesh.element(e);
            if el.etype.0 != 303 {
                return Err(SolverError::fatal(CALLER, "A simplicial mesh assumed currently"));
            }
            layout.element_slots(e, &mut slots);
            let nd = slots.len();
            let dofs: Vec<f64> = slots.iter().map(|&s| sol[s as usize]).collect();
            let (mut x, mut y, mut z, mut g) = ([0.0; 3], [0.0; 3], [0.0; 3], [0u32; 3]);
            for i in 0..3 {
                let c = ctx.mesh.coord(el.nodes[i]);
                (x[i], y[i], z[i]) = (c[0], c[1], c[2]);
                g[i] = el.nodes[i];
            }
            let mut ek = [0.0f64; 2];
            for t in 0..rule.len() {
                let fi = face_element_info_triangle(&x, &y, &z, &g, rule.u[t], rule.v[t], second, true);
                let (mut xq, mut yq) = (0.0, 0.0);
                for i in 0..3 {
                    xq += x[i] * fi.basis[i];
                    yq += y[i] * fi.basis[i];
                }
                let w = fi.det_f * rule.s[t];
                let solv = 16.0 * (xq - xq * xq) * (yq - yq * yq);
                let grad = [16.0 * (1.0 - 2.0 * xq) * (yq - yq * yq), 16.0 * (1.0 - 2.0 * yq) * (xq - xq * xq)];
                sol_norm += solv * solv * w;
                ek[0] += (solv - dofs[nd - 1]) * (solv - dofs[nd - 1]) * w;
                flux_norm += (grad[0] * grad[0] + grad[1] * grad[1]) * w;
                let mut fs = [0.0f64; 2];
                for p in 0..nd - 1 {
                    for k in 0..2 {
                        fs[k] += fi.fbasis[p][k] * dofs[p];
                    }
                }
                ek[1] += ((grad[0] - fs[0]) * (grad[0] - fs[0]) + (grad[1] - fs[1]) * (grad[1] - fs[1])) * w;
            }
            pres_err += ek[0];
            flux_err += ek[1];
        }
        println!("L2 Scalar Error = {}", fortran_e(pres_err.sqrt() / sol_norm.sqrt(), 16, 8));
        println!("L2 Flux Error ={}", fortran_e(flux_err.sqrt() / flux_norm.sqrt(), 16, 8));
        Ok(())
    }
}
