//! `bentonite/BentoniteSolver.F90`: coupled moisture (`Xi`), temperature
//! (`T`) and vapour (`Eta`) in bentonite, one 3-dof variable, with the
//! material model written into the element routine and a Robin heat
//! exchange on boundaries with `Robin Condition = True`.
//!
//! On the first visit the third component is replaced at every node by
//! its equilibrium value from the first two (`CompEta`). Each nonlinear
//! iteration assembles, solves, and stops when the relative change of the
//! solution norm, `2|Np - N|/(Np + N)`, falls below `Nonlinear System
//! Convergence Tolerance`. Integration weights carry `CoordinateSqrtMetric`
//! outside Cartesian coordinates.
//!
//! The transient heat capacity term of the Fortran reads `1.38*T + 355.6`
//! where `T` is, by Fortran's case-insensitivity, the integer Gauss-point
//! index `t`, in single precision. [`c_eff`] keeps that.

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementType, LocalSystem, ParamMap, Params, Quadrature, NO_DOF};
use elmer_fe::{CoordinateKind, CoordinateSystem, ElmerQuadrature, LagrangeBasis};
use elmer_solver::fortran_fmt::list_int;
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::{coordinate_system, element_nodes};
use crate::udf::list_real;

const CALLER: &str = "BentoniteSolve";
const DOFS: usize = 3;

const NU_0: f64 = 0.45;
const P_0: f64 = 1.013e5;
const RHO_1: f64 = 998.0;
const M_V: f64 = 0.018;
const R: f64 = 8.314351;
const K_U: f64 = 2.0e-21;
const MU_1: f64 = 1.0e-3;
const A: f64 = 0.4;
const D_0: f64 = 0.216e-4;
const N_D: f64 = 1.8;
const C_VP: f64 = 1.866e3;
const C_LP: f64 = 4.182e3;
const ETA_0: f64 = 0.023;
const XI_MAX: f64 = 0.999;
const T_0: f64 = 293.15;
const RHO_DRY: f64 = 1690.0;
const SMALL: f64 = 1.0e-12;

fn l_0() -> f64 {
    2454.3e3 - (C_VP - C_LP) * T_0
}

/// `CompEta`: the vapour content in equilibrium with `xi` and `temp`.
pub fn comp_eta(xi: f64, temp: f64) -> f64 {
    (M_V / (R * temp) * (l_0() * (temp - T_0) / T_0 + (C_VP - C_LP) * temp * (temp / T_0).ln())
        - A * (1.0 / xi.powi(2) - 1.0 / XI_MAX.powi(2)))
    .exp()
        * ETA_0
}

/// `c_eff` of `LocalMatrix` at Gauss point `t` (1-based): `1.38*T + 355.6`
/// is the single-precision product with the integer `t`.
pub fn c_eff(xi: f64, t: usize) -> f64 {
    let heat = 1.38_f32 * t as f32 + 355.6_f32;
    NU_0 * xi * RHO_1 * C_LP + RHO_DRY * f64::from(heat)
}

/// `LocalMatrix`: mass, stiffness and force for the 3-dof system.
#[allow(clippy::too_many_arguments)]
pub fn local_matrix(
    etype: ElementType,
    xyz: &[[f64; 3]],
    nodal: [&[f64]; 3],
    coords: CoordinateSystem,
    basis: &dyn ElementBasis,
    local: &mut LocalSystem,
) -> Result<(), String> {
    let n = xyz.len();
    let nn = n * DOFS;
    *local = LocalSystem::zeros(nn);
    let mut mass = vec![0.0; nn * nn];
    let dim = coords.dim;
    let mut ev = BasisEval::default();
    for (t, qp) in ElmerQuadrature.gauss_points(etype, None).into_iter().enumerate() {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        let (b, g) = (&ev.basis[..n], &ev.dbasis_dx[..n]);
        let mut s = qp.weight * ev.det_j;
        if coords.kind != CoordinateKind::Cartesian {
            let at = |k: usize| xyz.iter().zip(b).fold(0.0, |acc, (x, bp)| acc + x[k] * bp);
            s *= coords.info(at(0), at(1), at(2)).sqrt_metric;
        }
        let interp = |v: &[f64]| b.iter().zip(v).fold(0.0, |acc, (bp, x)| acc + bp * x);
        let xi = interp(nodal[0]);
        let temp = interp(nodal[1]);
        let eta = interp(nodal[2]).max(SMALL);
        let l0 = l_0();
        let ce = c_eff(xi, t + 1);
        let ar1 = 0.57;
        let ar2 = 1.28;
        let x_0 = 0.65;
        let dx = 0.1;
        let delta_eff = ar2 + (ar1 - ar2) / (1.0 + ((xi - x_0) / dx).exp());

        let a0 = NU_0 * (RHO_1 - eta * M_V / (R * temp) * P_0);
        let a1 = -NU_0 * (1.0 - xi) * eta * M_V / (R * temp.powi(2)) * P_0;
        let a2 = NU_0 * (1.0 - xi) * M_V / (R * temp) * P_0;
        let a3 = RHO_1 * K_U / MU_1 * 2.0 * A * RHO_1 * R / M_V * temp;
        let a4 = -RHO_1 * K_U / MU_1 * 2.0 * A * RHO_1 * R / M_V * (xi / XI_MAX) * (XI_MAX - xi);
        let a5 = M_V / (R * temp) * P_0 * NU_0 * (1.0 - xi) * D_0 * (temp / 273.0).powf(N_D);
        let b0 = -(l0 + (C_VP - C_LP) * temp) * NU_0 * RHO_1;
        let b1 = ce;
        let b2 = -(l0 + (C_VP - C_LP) * temp) * RHO_1 * K_U / MU_1 * 2.0 * A * RHO_1 * R * temp / M_V;
        let b3 = delta_eff;
        let c0 = -2.0 * A / xi.powi(3);
        let c1 = -(M_V * l0 / (R * temp.powi(2)) + M_V * (C_VP - C_LP) / (R * temp));
        let c2 = 1.0 / eta;

        let mut f = 0.0;
        f = f - A * (1.0 / xi.powi(2) - 1.0 / XI_MAX.powi(2)) - 2.0 * A / xi.powi(2);
        f = f + M_V * l0 / (R * T_0) - 2.0 * M_V * l0 / (R * temp);
        f = f + M_V * (C_VP - C_LP) / R * (temp / T_0).ln() - M_V * (C_VP - C_LP) / R;
        f = f - (eta / ETA_0).ln() + 1.0;
        let f = s * f;

        for p in 0..n {
            for q in 0..n {
                let base = s * b[p] * b[q];
                let grad = s * (0..dim).fold(0.0, |acc, k| acc + g[q][k] * g[p][k]);
                let mut lm = [[0.0; 3]; 3];
                let mut ls = [[0.0; 3]; 3];
                lm[0][0] += a0 * base;
                lm[0][1] += a1 * base;
                lm[0][2] += a2 * base;
                ls[0][0] += a3 * grad;
                ls[0][1] += a4 * grad;
                ls[0][2] += a5 * grad;
                lm[1][0] += b0 * base;
                lm[1][1] += b1 * base;
                ls[1][0] += b2 * grad;
                ls[1][1] += b3 * grad;
                ls[2][0] += c0 * base;
                ls[2][1] += c1 * base;
                ls[2][2] += c2 * base;
                for i in 0..DOFS {
                    for j in 0..DOFS {
                        let at = (DOFS * p + i) * nn + DOFS * q + j;
                        mass[at] += lm[i][j];
                        local.stiff[at] += ls[i][j];
                    }
                }
            }
            local.force[DOFS * p + 2] += f * b[p];
        }
    }
    local.mass = Some(mass);
    Ok(())
}

/// `LocalMatrixBoundary`: `H (T - Text)` heat exchange, `H = 1.9`,
/// `Text = 285.15`.
pub fn local_matrix_boundary(
    etype: ElementType,
    xyz: &[[f64; 3]],
    coords: CoordinateSystem,
    basis: &dyn ElementBasis,
    local: &mut LocalSystem,
) -> Result<(), String> {
    const H: f64 = 1.9;
    const TEXT: f64 = 285.15;
    let n = xyz.len();
    let nn = n * DOFS;
    *local = LocalSystem::zeros(nn);
    let mut ev = BasisEval::default();
    for qp in ElmerQuadrature.gauss_points(etype, None) {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        let b = &ev.basis[..n];
        let mut s = qp.weight * ev.det_j;
        if coords.kind != CoordinateKind::Cartesian {
            let at = |k: usize| xyz.iter().zip(b).fold(0.0, |acc, (x, bp)| acc + x[k] * bp);
            s *= coords.sqrt_metric(at(0), at(1), at(2));
        }
        for p in 0..n {
            for q in 0..n {
                local.stiff[(DOFS * p + 1) * nn + DOFS * q + 1] += s * H * b[q] * b[p];
            }
            local.force[DOFS * p + 1] += s * H * TEXT * b[p];
        }
    }
    Ok(())
}

/// `AllocationsDone` (`SAVE`d): `CompEta` has initialised the field.
#[derive(Default)]
pub struct BentoniteSolver {
    initialised: bool,
}

pub fn factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    if params.string("element").is_some_and(|e| elmer_core::normalize_key(e) != "n:1") {
        return Err(SolverError::Unsupported("BentoniteSolver: non-nodal elements".into()));
    }
    Ok(Box::new(BentoniteSolver::default()))
}

/// `Xi`, `T`, `Eta` at the element's nodes.
fn nodal_fields(ctx: &SolverCtx<'_>, e: ElemId) -> Result<[Vec<f64>; 3], SolverError> {
    let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
    let mut out = [Vec::new(), Vec::new(), Vec::new()];
    for &node in ctx.mesh.element(e).nodes {
        let slot = v.perm.and_then(|p| p.get(node as usize).copied()).filter(|&s| s != NO_DOF);
        let slot = slot.ok_or_else(|| SolverError::fatal(CALLER, format!("node {} has no dof", node + 1)))? as usize;
        for (k, o) in out.iter_mut().enumerate() {
            o.push(v.values[DOFS * slot + k]);
        }
    }
    Ok(out)
}

impl Solver for BentoniteSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        crate::hooks::check_element_definitions(ctx)?;
        if ctx.primary().is_none_or(|v| v.dofs != DOFS) {
            return Err(SolverError::Unsupported("BentoniteSolver: Variable DOFs = 3".into()));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        let coords = coordinate_system(ctx)?;
        let basis = LagrangeBasis::new(coords.dim);
        let mut norm = ctx.convergence().map_or(0.0, |c| c.norm);
        if !self.initialised {
            let v = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
            for cell in v.values.as_chunks_mut::<DOFS>().0 {
                cell[2] = comp_eta(cell[0], cell[1]);
            }
            self.initialised = true;
        }
        let tol = ctx.params().real("nonlinear system convergence tolerance").unwrap_or(0.0);
        let iters = ctx.params().integer("nonlinear system max iterations").unwrap_or(1);
        let bulk: Vec<ElemId> = ctx.active_elements().collect();
        let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut local = LocalSystem::zeros(0);
        let mut dofs = Vec::new();
        let rule = "-------------------------------------";
        for iter in 1..=iters {
            ctx.info(CALLER, " ", 4);
            ctx.info(CALLER, rule, 4);
            ctx.info(CALLER, &format!(" Bentonite iteration{}", list_int(iter)), 4);
            ctx.info(CALLER, rule, 4);
            ctx.info(CALLER, " ", 4);
            ctx.info(CALLER, "Starting Assembly", 4);
            ctx.default_initialize()?;
            for &e in &bulk {
                let etype = ctx.mesh.element(e).etype;
                let xyz = element_nodes(ctx.mesh, e);
                let [xi, temp, eta] = nodal_fields(ctx, e)?;
                local_matrix(etype, &xyz, [&xi, &temp, &eta], coords, &basis, &mut local)
                    .map_err(|m| SolverError::fatal(CALLER, m))?;
                if transient {
                    ctx.default_first_order_time(&mut local, e, dt)?;
                }
                ctx.default_update(&local, e)?;
            }
            for &e in &boundary {
                ctx.element_dofs(e, &mut dofs);
                if dofs.contains(&NO_DOF) {
                    continue;
                }
                if ctx.bc_params(e).and_then(|bc| bc.logical("robin condition")) != Some(true) {
                    continue;
                }
                let etype = ctx.mesh.element(e).etype;
                let xyz = element_nodes(ctx.mesh, e);
                local_matrix_boundary(etype, &xyz, coords, &basis, &mut local)
                    .map_err(|m| SolverError::fatal(CALLER, m))?;
                if transient {
                    local.mass = Some(vec![0.0; local.n * local.n]);
                    ctx.default_first_order_time(&mut local, e, dt)?;
                }
                ctx.default_update(&local, e)?;
            }
            // DefaultFinishAssembly, DefaultDirichletBCs
            ctx.default_dirichlet()?;
            let prev = norm;
            norm = ctx.default_solve()?;
            let change = if prev + norm != 0.0 { 2.0 * (prev - norm).abs() / (prev + norm) } else { 0.0 };
            ctx.info(CALLER, " ", 4);
            ctx.info(CALLER, &format!(" Result Norm    : {}", list_real(norm)), 4);
            ctx.info(CALLER, &format!(" Relative Change: {}", list_real(change)), 4);
            if change < tol {
                break;
            }
        }
        Ok(())
    }
}
