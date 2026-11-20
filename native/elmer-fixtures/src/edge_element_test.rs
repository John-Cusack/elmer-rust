//! `EdgeElementTest.F90` of the `SecondOrderEdgeElement2D` and
//! `SecondOrderEdgeElement3D` tests (`Procedure = "EdgeElementTest"
//! "EdgeElementSolver"`): the H(curl) projection `(u, v) + (curl u, curl v)`
//! of a manufactured field in the second-order Piola edge element space
//! (`Element = "n:0 e:2 b:2"` on triangles, `"n:0 e:2 -tri_face b:2"` on
//! tetrahedra, bubbles in the global system), and with `Error Computation`
//! the relative energy-norm error printed as `Error Norm = ...`.
//!
//! The two plugins share the assembly and differ in the field: in 2D
//! `u = (-y^3/3, x^3/3)`, `curl u = x^2 + y^2`; in 3D `u = (0, 0, (y x^3 +
//! x y^3)/6)` with `MatPar = 1` on the curl term. The Gauss rule is
//! `GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=.TRUE.,
//! EdgeBasisDegree=2)`: the p triangle's 6 points, the p tetrahedron's 11.
//! Refused: `Tabulate Basis` (the pretabulated variant of the 2D plugin),
//! other element families.

use elmer_core::{ElemId, LocalSystem, Params, NO_DOF};
use elmer_fe_hp::edge::{edge_element_info, p_reference_nodal_basis};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "EdgeElementSolver";

pub fn factory_2d(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(EdgeElementSolver { three_d: false }))
}

pub fn factory_3d(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(EdgeElementSolver { three_d: true }))
}

struct EdgeElementSolver {
    /// The `SecondOrderEdgeElement3D` plugin.
    three_d: bool,
}

/// One integration point: `detJ * IP % s(t)` split as the Fortran keeps
/// them, the nodal basis, the edge basis and its curl.
struct Point {
    det: f64,
    s: f64,
    basis: Vec<f64>,
    e: Vec<[f64; 3]>,
    c: Vec<[f64; 3]>,
}

fn points(ctx: &SolverCtx<'_>, e: ElemId, nd: usize) -> Result<(Vec<Point>, Vec<[f64; 3]>), SolverError> {
    let el = ctx.mesh.element(e);
    let etype = el.etype;
    let rule = match etype.0 {
        303 => elmer_pdofs::quadrature::gauss_points_triangle_p(6),
        504 => elmer_pdofs::quadrature::gauss_points_tetra_p(11),
        c => return Err(SolverError::Unsupported(format!("{CALLER} on element type {c}"))),
    }
    .map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let gindex: Vec<i64> = el.nodes.iter().map(|&k| k as i64 + 1).collect();
    let cdim = ctx.mesh.dim();
    let mut out = Vec::with_capacity(rule.len());
    for t in 0..rule.len() {
        let uvw = [rule.u[t], rule.v[t], rule.w[t]];
        let ev = edge_element_info(etype, &xyz, &gindex, uvw, 2, cdim).map_err(|m| SolverError::fatal(CALLER, m))?;
        if ev.dofs != nd {
            return Err(SolverError::fatal(CALLER, format!("{} edge basis functions for {nd} element dofs", ev.dofs)));
        }
        out.push(Point { det: ev.det, s: rule.s[t], basis: p_reference_nodal_basis(etype, uvw), e: ev.edge, c: ev.curl });
    }
    Ok((out, xyz))
}

/// `SUM(a(1:n) * b(1:n))`.
fn sum(a: impl Iterator<Item = f64>) -> f64 {
    let mut s = 0.0;
    for x in a {
        s += x;
    }
    s
}

impl Solver for EdgeElementSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::fatal(CALLER, "the edge element dofs (no DofsBackend)"));
        }
        let p = ctx.params();
        let error_estimation = p.logical("error computation").unwrap_or(false);
        let energy_norm = p.logical("use energy norm").unwrap_or(false);
        if !self.three_d && p.logical("tabulate basis") == Some(true) {
            return Err(SolverError::Unsupported(format!("{CALLER}: Tabulate Basis")));
        }
        let dim = ctx.mesh.dim();
        let mat_par = 1.0;
        ctx.default_initialize()?;
        let nb = ctx.mesh.n_bulk();
        let active: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nb).collect();
        let mut slots = Vec::new();
        for &e in &active {
            ctx.element_slots(e, &mut slots);
            if slots.contains(&NO_DOF) {
                return Err(SolverError::fatal(CALLER, format!("element {} lacks dofs", e.0)));
            }
            let nd = slots.len();
            if !self.three_d && ctx.mesh.element(e).etype.family() == 5 && nd != 6 {
                return Err(SolverError::fatal(CALLER, "Indices for a tetrahedron erratic"));
            }
            let (pts, xyz) = points(ctx, e, nd)?;
            let n = xyz.len();
            let mut stiff = vec![0.0; nd * nd];
            let mut force = vec![0.0; nd];
            for pt in &pts {
                let (dj, s) = (pt.det, pt.s);
                let xq = sum((0..n).map(|i| xyz[i][0] * pt.basis[i]));
                let yq = sum((0..n).map(|i| xyz[i][1] * pt.basis[i]));
                for pp in 0..nd {
                    for q in 0..nd {
                        // STIFF(i,j) = STIFF(i,j) + mass term + curl term, left to right.
                        let t1 = 1.0 * sum((0..dim).map(|k| pt.e[q][k] * pt.e[pp][k])) * dj * s;
                        let t2 = if self.three_d {
                            mat_par * sum((0..dim).map(|k| pt.c[q][k] * pt.c[pp][k])) * dj * s
                        } else {
                            pt.c[q][2] * pt.c[pp][2] * dj * s
                        };
                        let a = &mut stiff[pp * nd + q];
                        *a = (*a + t1) + t2;
                    }
                    let (e, c) = (&pt.e[pp], &pt.c[pp]);
                    force[pp] = if self.three_d {
                        force[pp]
                            + (0.0) * e[0] * dj * s
                            + (0.0) * e[1] * dj * s
                            + 1.0 / 6.0 * (yq * (xq * xq * xq) + xq * (yq * yq * yq)) * e[2] * dj * s
                            + mat_par * (1.0 / 6.0 * (xq * xq * xq) + xq * (yq * yq) / 2.0) * c[0] * dj * s
                            + mat_par * (-(1.0 / 6.0 * (yq * yq * yq)) - yq * (xq * xq) / 2.0) * c[1] * dj * s
                            + mat_par * (0.0) * c[2] * dj * s
                    } else {
                        force[pp] - (yq * yq * yq) / 3.0 * e[0] * dj * s
                            + (xq * xq * xq) / 3.0 * e[1] * dj * s
                            + (xq * xq + yq * yq) * c[2] * dj * s
                    };
                }
            }
            let ls = LocalSystem { n: nd, stiff, mass: None, damp: None, force };
            ctx.default_update(&ls, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;

        if error_estimation {
            let vals = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?.values.to_vec();
            let mut err = 0.0;
            let mut sol_norm = 0.0;
            for &e in &active {
                ctx.element_slots(e, &mut slots);
                let nd = slots.len();
                let load: Vec<f64> = slots.iter().map(|&s| vals[s as usize]).collect();
                let (pts, xyz) = points(ctx, e, nd)?;
                let n = xyz.len();
                for pt in &pts {
                    let xq = sum((0..n).map(|i| xyz[i][0] * pt.basis[i]));
                    let yq = sum((0..n).map(|i| xyz[i][1] * pt.basis[i]));
                    let comp = |f: &dyn Fn(usize) -> f64| sum((0..nd).map(|j| load[j] * f(j)));
                    let (sol, rotsol, u, rotu) = if self.three_d {
                        let u = [comp(&|j| pt.e[j][0]), comp(&|j| pt.e[j][1]), comp(&|j| pt.e[j][2])];
                        let rotu = [comp(&|j| pt.c[j][0]), comp(&|j| pt.c[j][1]), comp(&|j| pt.c[j][2])];
                        let sol = [0.0, 0.0, 1.0 / 6.0 * (yq * (xq * xq * xq) + xq * (yq * yq * yq))];
                        let rotsol = [
                            0.0 + 1.0 / 6.0 * (xq * xq * xq) + xq * (yq * yq) / 2.0,
                            0.0 - 1.0 / 6.0 * (yq * yq * yq) - yq * (xq * xq) / 2.0,
                            0.0,
                        ];
                        (sol, rotsol, u, rotu)
                    } else {
                        let u = [comp(&|j| pt.e[j][0]), comp(&|j| pt.e[j][1]), 0.0];
                        let rotu = [0.0, 0.0, comp(&|j| pt.c[j][2])];
                        let sol = [-(yq * yq * yq) / 3.0, (xq * xq * xq) / 3.0, 0.0];
                        let rotsol = [0.0, 0.0, xq * xq + yq * yq];
                        (sol, rotsol, u, rotu)
                    };
                    let e3: [f64; 3] = if self.three_d {
                        std::array::from_fn(|k| sol[k] - u[k])
                    } else {
                        [sol[0] - u[0], sol[1] - u[1], 0.0]
                    };
                    let rote: [f64; 3] = std::array::from_fn(|k| rotsol[k] - rotu[k]);
                    let ss = |v: &[f64; 3]| sum(v.iter().map(|x| x * x));
                    if energy_norm {
                        sol_norm += (ss(&sol) + mat_par * ss(&rotsol)) * pt.det * pt.s;
                        err += (ss(&e3) + mat_par * ss(&rote)) * pt.det * pt.s;
                    } else {
                        sol_norm += ss(&sol) * pt.det * pt.s;
                        err += ss(&e3) * pt.det * pt.s;
                    }
                }
            }
            let msg = format!("Error Norm = {}", crate::poisson_ml::fortran_e(err.sqrt() / sol_norm.sqrt(), 16, 8));
            ctx.info(CALLER, &msg, 1);
        }
        Ok(())
    }
}
