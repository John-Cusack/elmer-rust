//! `FaceFEInterpolationTest.F90` of the `FaceElement2D_BCs` and
//! `FaceElement3D_BCs` tests (`Procedure = "FaceFEInterpolationTest"
//! "BestApproximationSolver"`): the L2 best approximation of the constant
//! field `(1, 1)` (`Test Mode = 1`) or `(1, 1, 1)` (`Test Mode = 2`) in a
//! face element space (RT(0)/BDM(1) on triangles, the Nedelec tetrahedra of
//! the first and second kind), the boundary dofs fixed by `<Variable> {f}
//! j` Dirichlet conditions; with `Error Computation` (the default) it
//! prints `Error Norm = ...`, the relative L2 error.
//!
//! The Gauss rule is `GaussPoints(Element, EdgeBasis = .TRUE.,
//! PReferenceElement = .TRUE., EdgeBasisDegree = 1)`: the p-element
//! triangle's 3 points and tetrahedron's 4, as `FaceCell::rule` gives
//! them. Refused: bricks, elementwise (bubble) dofs after the face ones,
//! and `Compute Div Error` with a `Test Mode` above 2 (which the Fortran
//! does not reach either).

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{ElemId, LocalSystem, Params, NO_DOF};
use elmer_pdofs::face_element::{face_dirichlet, nodal_dofs, FaceCell};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::poisson_ml::fortran_e;

const CALLER: &str = "BestApproximationSolver";

pub fn factory(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(BestApproximation))
}

struct BestApproximation;

/// The field the test approximates, `sol(1:3)`.
fn target(test_mode: i32) -> Result<[f64; 3], SolverError> {
    match test_mode {
        1 => Ok([1.0, 1.0, 0.0]),
        2 => Ok([1.0, 1.0, 1.0]),
        _ => Err(SolverError::fatal("FaceFEInterpolationTest", "This Test Mode needs an update")),
    }
}

/// The cell of active element `e` (bricks refused) and its nodal dof count.
fn cell_of(ctx: &SolverCtx<'_>, e: ElemId, nd: usize, second: bool) -> Result<(FaceCell, usize), SolverError> {
    let cell = FaceCell::of(ctx, e)?;
    if cell.brick() {
        return Err(SolverError::Unsupported("BestApproximationSolver on bricks".into()));
    }
    let np = nodal_dofs(&cell, nd, second, 0)?;
    Ok((cell, np))
}

impl Solver for BestApproximation {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::fatal(CALLER, "the face element dofs (no DofsBackend)"));
        }
        let p = ctx.params();
        let second = p.logical("second kind basis").unwrap_or(false);
        let error_estimation = p.logical("error computation").unwrap_or(true);
        let test_mode = p.integer("test mode").unwrap_or(1);
        let name = ctx.variable_name().unwrap_or("").to_ascii_lowercase();
        let sol = target(test_mode)?;
        ctx.default_initialize()?;

        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut slots = Vec::new();
        for &e in &active {
            ctx.element_slots(e, &mut slots);
            if slots.contains(&NO_DOF) {
                return Err(SolverError::fatal(CALLER, format!("element {} lacks dofs", e.0)));
            }
            let nd = slots.len();
            let (cell, np) = cell_of(ctx, e, nd, second)?;
            let (n, dim) = (cell.n(), cell.dim());
            let bp = ctx.body_params(e);
            let acoef = match bp.as_ref().and_then(|b| b.material) {
                Some(m) => ctx.real_at_nodes(m, "material param", e)?.unwrap_or_else(|| vec![0.0; n]),
                None => vec![0.0; n],
            };
            let rule = cell.rule()?;
            let mut local = LocalSystem::zeros(nd);
            for t in 0..rule.len() {
                let fi = cell.info(&rule, t, second);
                let (detj, w) = (fi.det_f, rule.s[t]);
                let mut matpar = 0.0;
                for i in 0..n {
                    matpar += acoef[i] * fi.basis[i];
                }
                if np > 0 {
                    for p in 0..n {
                        for q in 0..n {
                            local.stiff[p * nd + q] += fi.basis[p] * fi.basis[q] * detj * w;
                        }
                        for q in 0..nd - np {
                            let j = np + q;
                            local.stiff[p * nd + j] -= fi.fbasis[q][0] * fi.basis[p] * detj * w;
                        }
                    }
                }
                for p in 0..nd - np {
                    let i = np + p;
                    for q in 0..nd - np {
                        let j = np + q;
                        let mut d = 0.0;
                        for k in 0..dim {
                            d += fi.fbasis[q][k] * fi.fbasis[p][k];
                        }
                        local.stiff[i * nd + j] += 1.0 * d * detj * w + matpar * fi.div[q] * fi.div[p] * detj * w;
                    }
                    let fb = fi.fbasis[p];
                    local.force[i] += if test_mode == 1 {
                        1.0 * fb[0] * detj * w + 1.0 * fb[1] * detj * w
                    } else {
                        1.0 * fb[0] * detj * w + 1.0 * fb[1] * detj * w + 1.0 * fb[2] * detj * w
                    };
                }
            }
            ctx.default_update(&local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        let fixed = face_dirichlet(ctx, &name, second, 0)?;
        ctx.default_dirichlet_after(&fixed)?;
        ctx.default_solve()?;

        if error_estimation {
            let values = ctx
                .variable(&name)
                .ok_or_else(|| SolverError::fatal(CALLER, format!("no variable {name}")))?
                .values
                .to_vec();
            let (mut err, mut sol_norm) = (0.0f64, 0.0f64);
            for &e in &active {
                ctx.element_slots(e, &mut slots);
                let nd = slots.len();
                let (cell, np) = cell_of(ctx, e, nd, second)?;
                let dim = cell.dim();
                let load: Vec<f64> = slots.iter().map(|&s| values[s as usize]).collect();
                let rule = cell.rule()?;
                for t in 0..rule.len() {
                    let fi = cell.info(&rule, t, second);
                    let (detj, w) = (fi.det_f, rule.s[t]);
                    let mut u = [0.0f64; 3];
                    for i in 0..dim {
                        let mut s = 0.0;
                        for q in 0..nd - np {
                            s += load[np + q] * fi.fbasis[q][i];
                        }
                        u[i] = s;
                    }
                    let mut ev = [0.0f64; 3];
                    for i in 0..dim {
                        ev[i] = sol[i] - u[i];
                    }
                    sol_norm += (sol[0] * sol[0] + sol[1] * sol[1] + sol[2] * sol[2]) * detj * w;
                    err += (ev[0] * ev[0] + ev[1] * ev[1] + ev[2] * ev[2]) * detj * w;
                }
            }
            // ParallelReduction.
            let mut buf = [err, sol_norm];
            ctx.comm.sum_f64(&mut buf);
            let [err, sol_norm] = buf;
            println!("Error Norm = {}", fortran_e(err.sqrt() / sol_norm.sqrt(), 16, 8));
        }
        Ok(())
    }
}
