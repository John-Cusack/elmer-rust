//! `ComputeError.F90` of the `RaviartThomas2D`, `BDM2D` and
//! `RaviartThomas3D` tests (`Procedure = "ComputeError" "IntegrateError"`):
//! the relative L2 errors of the potential and the flux of the `Target
//! Equation`'s mixed Poisson solution against `u = 16 (x - x^2)(y - y^2)`
//! (2D, triangles) or `u = -64 (x^2 - 1/4)(y^2 - 1/4)(z^2 - 1/4)` (3D,
//! tetrahedra), printed as `L2 Scalar Error = ...` and `L2 Flux Error =
//! ...`.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{ElemId, Params, Value};
use elmer_fe_hp::face::{face_element_info_tetra, face_element_info_triangle};
use elmer_pdofs::quadrature::{gauss_points_tetra_p, gauss_points_triangle_p};
use elmer_pdofs::PDofsBackend;
use elmer_solver::{DofsBackend, Solver, SolverCtx, SolverError};

use crate::poisson_ml::fortran_e;

const CALLER: &str = "ComputeError";

/// The 2D `ComputeError.F90` (`RaviartThomas2D`, `BDM2D`).
pub fn factory(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(IntegrateError { three_d: false }))
}

/// The 3D `ComputeError.F90` (`RaviartThomas3D`).
pub fn factory_3d(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(IntegrateError { three_d: true }))
}

struct IntegrateError {
    three_d: bool,
}

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
        let (rule, n, code) = if self.three_d {
            (gauss_points_tetra_p(4), 4, 504)
        } else {
            (gauss_points_triangle_p(3), 3, 303)
        };
        let rule = rule.map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
        let dim = if self.three_d { 3 } else { 2 };
        let (mut pres_err, mut flux_err, mut sol_norm, mut flux_norm) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut slots = Vec::new();
        for &e in &active {
            let el = ctx.mesh.element(e);
            if el.etype.0 != code {
                return Err(SolverError::fatal(CALLER, "A simplicial mesh assumed currently"));
            }
            layout.element_slots(e, &mut slots);
            let nd = slots.len();
            let dofs: Vec<f64> = slots.iter().map(|&s| sol[s as usize]).collect();
            // `np`: nodal dofs before the face dofs (an `n:1` definition).
            let nf = (if self.three_d { 4 } else { 3 }) * if second { if self.three_d { 3 } else { 2 } } else { 1 };
            let np = nd - nf - 1;
            let (mut x, mut y, mut z, mut g) = ([0.0; 4], [0.0; 4], [0.0; 4], [0u32; 4]);
            for i in 0..n {
                let c = ctx.mesh.coord(el.nodes[i]);
                (x[i], y[i], z[i]) = (c[0], c[1], c[2]);
                g[i] = el.nodes[i];
            }
            let mut ek = [0.0f64; 2];
            for t in 0..rule.len() {
                let fi = if self.three_d {
                    face_element_info_tetra(&x, &y, &z, &g, rule.u[t], rule.v[t], rule.w[t], second, true)
                } else {
                    let (x3, y3, z3, g3) = ([x[0], x[1], x[2]], [y[0], y[1], y[2]], [z[0], z[1], z[2]], [g[0], g[1], g[2]]);
                    face_element_info_triangle(&x3, &y3, &z3, &g3, rule.u[t], rule.v[t], second, true)
                };
                let (mut xq, mut yq, mut zq) = (0.0, 0.0, 0.0);
                for i in 0..n {
                    xq += x[i] * fi.basis[i];
                    yq += y[i] * fi.basis[i];
                    zq += z[i] * fi.basis[i];
                }
                let (detj, ws) = (fi.det_f, rule.s[t]);
                let (solv, grad) = if self.three_d {
                    let (a, b, c) = (xq * xq - 1.0 / 4.0, yq * yq - 1.0 / 4.0, zq * zq - 1.0 / 4.0);
                    (-64.0 * a * b * c, [-64.0 * 2.0 * xq * b * c, -64.0 * 2.0 * yq * a * c, -64.0 * 2.0 * zq * a * b])
                } else {
                    let (a, b) = (xq - xq * xq, yq - yq * yq);
                    (16.0 * a * b, [16.0 * (1.0 - 2.0 * xq) * b, 16.0 * (1.0 - 2.0 * yq) * a, 0.0])
                };
                sol_norm += solv * solv * detj * ws;
                ek[0] += (solv - dofs[nd - 1]) * (solv - dofs[nd - 1]) * detj * ws;
                let mut gg = 0.0;
                for k in 0..dim {
                    gg += grad[k] * grad[k];
                }
                flux_norm += gg * detj * ws;
                let mut fs = [0.0f64; 3];
                for p in 0..nf {
                    for k in 0..dim {
                        fs[k] += fi.fbasis[p][k] * dofs[np + p];
                    }
                }
                let mut d = 0.0;
                for k in 0..dim {
                    d += (grad[k] - fs[k]) * (grad[k] - fs[k]);
                }
                ek[1] += d * detj * ws;
            }
            pres_err += ek[0];
            flux_err += ek[1];
        }
        println!("L2 Scalar Error = {}", fortran_e(pres_err.sqrt() / sol_norm.sqrt(), 16, 8));
        println!("L2 Flux Error ={}", fortran_e(flux_err.sqrt() / flux_norm.sqrt(), 16, 8));
        Ok(())
    }
}
