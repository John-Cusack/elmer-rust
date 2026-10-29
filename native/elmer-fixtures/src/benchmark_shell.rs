//! `BenchmarkShell.F90` (`p-FEM_two_solvers`, `p-FEM_with_varying_p`,
//! `SD_Classical2DShell` and their `SD_` variants): `SUBROUTINE
//! ShellSolver`, a Reissner-Mindlin model of a cylindrical shell of unit
//! radius in the 2D parameter plane (`y1` angular, `y2` axial), on p
//! elements.
//!
//! The unknowns per dof are `u1, u2, u3` (in-surface and normal
//! displacement) and the rotations `b1, b2`, plus an out-of-plane strain
//! with `Variable DOFs = 6`. Membrane (`BM`), shear (`BS`) and bending
//! (`BB`) strains go through `StrainEnergyDensity` with the plane stress
//! matrix of `ElasticityMatrix`. The load is the normal traction
//! `h 1e9 cos(2 y1)`; the `Normal Pressure` keyword is read and then
//! overwritten, as in the Fortran. After the solve the routine prints the
//! relative energy error against the benchmark's reference work.
//!
//! The local system covers `nd + nb` dofs; Elmer's `DefaultUpdateEquations`
//! condenses the local bubbles, which here happens before
//! [`SolverCtx::default_update`]. Bases, rules and dof counts come from the
//! solver's dof layout (`elmer-pdofs`).

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementType, LocalSystem, ParamMap, Params, QuadPoint};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::{element, nodal};
use crate::linalg::condensate;
use crate::udf::{list_real, stdout_sink, PrintSink};

const CALLER: &str = "ShellSolver";

pub struct ShellSolver {
    out: PrintSink,
    local: LocalSystem,
}

/// [`SolverFactory`](elmer_solver::SolverFactory): prints to standard output.
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(ShellSolver { out: stdout_sink(), local: LocalSystem::zeros(0) }))
}

/// `ElasticityMatrix(CMat, GMat, 1, 1, E, nu)` (row-major 4x4, 2x2).
pub fn elasticity_matrix(e: f64, nu: f64) -> ([[f64; 4]; 4], [[f64; 2]; 2]) {
    let (a1, a2) = (1.0_f64, 1.0_f64);
    let mut c = [[0.0; 4]; 4];
    c[0][0] = 1.0;
    c[0][1] = nu;
    c[1][0] = nu;
    c[1][1] = 1.0;
    c[2][2] = (1.0 - nu) / 2.0;
    let f = e / (1.0 - nu.powi(2));
    for row in &mut c {
        for v in row.iter_mut() {
            *v *= f;
        }
    }
    c[0][0] /= a1.powi(4);
    c[0][1] /= a1.powi(2) * a2.powi(2);
    c[1][0] /= a1.powi(2) * a2.powi(2);
    c[1][1] /= a2.powi(4);
    c[2][2] /= a1.powi(2) * a2.powi(2);
    c[3][3] = (1.0 - nu) * e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mut g = [[0.0; 2]; 2];
    g[0][0] = e / (2.0 * (1.0 + nu) * a1.powi(2));
    g[1][1] = e / (2.0 * (1.0 + nu) * a2.powi(2));
    (c, g)
}

/// `StrainEnergyDensity(A, B, C, m, n, s)`: `A(i,j) += C(k,i) B(k,l) C(l,j) s`
/// in the Fortran's loop order; `c` is `m x n` row-major, `b` `m x m`.
fn strain_energy_density(a: &mut [f64], b: &[f64], c: &[f64], m: usize, n: usize, s: f64) {
    for i in 0..n {
        for j in 0..n {
            let mut acc = a[i * n + j];
            for k in 0..m {
                for l in 0..m {
                    acc += c[k * n + i] * b[k * m + l] * c[l * n + j] * s;
                }
            }
            a[i * n + j] = acc;
        }
    }
}

/// Nodal material values of one element.
pub struct ShellNodal<'a> {
    pub poisson_ratio: &'a [f64],
    pub youngs_modulus: &'a [f64],
    pub thickness: &'a [f64],
}

/// `ShellLocalMatrix` over `nd` dofs of `m` components (`csize` membrane
/// strains), into `local` (`m*nd` square). The strain matrices are filled
/// entry by entry at the Fortran's indices.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
pub fn shell_local_matrix(
    etype: ElementType,
    xyz: &[[f64; 3]],
    nd: usize,
    m: usize,
    csize: usize,
    mat: &ShellNodal<'_>,
    basis: &dyn ElementBasis,
    rule: &[QuadPoint],
    local: &mut LocalSystem,
) -> Result<(), String> {
    const R: f64 = 1.0;
    let n = xyz.len();
    let nn = m * nd;
    *local = LocalSystem::zeros(nn);
    let mut bm = vec![0.0; csize * nn];
    let mut bs = vec![0.0; 2 * nn];
    let mut bb = vec![0.0; 3 * nn];
    let mut ev = BasisEval::default();
    let interp = |v: &[f64], b: &[f64]| v.iter().zip(b).fold(0.0, |s, (x, bp)| s + x * bp);
    for qp in rule {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        if ev.basis.len() < nd {
            return Err(format!("ElementInfo gave {} basis functions for {nd} dofs", ev.basis.len()));
        }
        let (b, g) = (&ev.basis[..nd], &ev.dbasis_dx[..nd]);
        let bn = &b[..n];
        let nu = interp(mat.poisson_ratio, bn);
        let e = interp(mat.youngs_modulus, bn);
        let h = interp(mat.thickness, bn);
        let y1 = xyz.iter().zip(bn).fold(0.0, |s, (x, bp)| s + x[0] * bp);
        let traction = h * 1.0e9 * (2.0 * y1).cos();
        let (cmat, gmat) = elasticity_matrix(e, nu);
        // IP % s(t)
        let s = qp.weight;

        let weight = h * ev.det_j * s;
        for p in 0..nd {
            let col = |k: usize| p * m + k;
            bm[col(0)] = g[p][0];
            bm[col(2)] = b[p] / R;
            bm[nn + col(1)] = g[p][1];
            bm[2 * nn + col(0)] = g[p][1];
            bm[2 * nn + col(1)] = g[p][0];
            if m == 6 {
                bm[3 * nn + col(5)] = b[p];
                bm[3 * nn + col(0)] = -nu / (1.0 - nu) * g[p][0];
                bm[3 * nn + col(1)] = -nu / (1.0 - nu) * g[p][1];
                bm[3 * nn + col(2)] = -nu / (1.0 - nu) * b[p] / R;
            }
        }
        let c: Vec<f64> = (0..csize).flat_map(|i| (0..csize).map(move |j| cmat[i][j])).collect();
        strain_energy_density(&mut local.stiff, &c, &bm, csize, nn, weight);

        for p in 0..nd {
            let col = |k: usize| p * m + k;
            bs[col(0)] = -b[p] / R;
            bs[col(2)] = g[p][0];
            bs[col(3)] = -b[p];
            bs[nn + col(2)] = g[p][1];
            bs[nn + col(4)] = -b[p];
        }
        let gm = [gmat[0][0], gmat[0][1], gmat[1][0], gmat[1][1]];
        strain_energy_density(&mut local.stiff, &gm, &bs, 2, nn, weight);

        let weight = h.powi(3) / 12.0 * ev.det_j * s;
        for p in 0..nd {
            let col = |k: usize| p * m + k;
            bb[col(3)] = g[p][0];
            bb[nn + col(4)] = g[p][1];
            bb[2 * nn + col(3)] = g[p][1];
            bb[2 * nn + col(4)] = g[p][0];
            bb[2 * nn + col(0)] = -1.0 / R * g[p][1];
        }
        let c3: Vec<f64> = (0..3).flat_map(|i| (0..3).map(move |j| cmat[i][j])).collect();
        strain_energy_density(&mut local.stiff, &c3, &bb, 3, nn, weight);

        let weight = ev.det_j * s;
        for p in 0..nd {
            local.force[m * p + 2] += traction * b[p] * weight;
        }
    }
    Ok(())
}

impl Solver for ShellSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        crate::hooks::check_element_definitions(ctx)?;
        let dofs = ctx.params().integer("variable dofs");
        if !matches!(dofs, Some(5 | 6)) {
            return Err(SolverError::fatal(CALLER, format!("Variable DOFs must be 5 or 6, got {dofs:?}")));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let m = ctx.params().integer("variable dofs").unwrap_or(5) as usize;
        let csize = if m == 6 { 4 } else { 3 };
        let fatal = |msg: String| SolverError::fatal(CALLER, msg);
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        for e in active {
            let el = element(ctx, e)?;
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let material = ctx.body_params(e).and_then(|b| b.material);
            let o = &el.nodes;
            let pr = nodal(ctx, material, "poisson ratio", e, o)?;
            let ym = nodal(ctx, material, "youngs modulus", e, o)?;
            let th = nodal(ctx, material, "shell thickness", e, o)?;
            // `Normal Pressure` is read into LOAD and never used.
            let bf = ctx.body_params(e).and_then(|b| b.body_force);
            nodal(ctx, bf, "normal pressure", e, o)?;
            let mat = ShellNodal { poisson_ratio: &pr, youngs_modulus: &ym, thickness: &th };
            let total = el.nd + el.nb;
            shell_local_matrix(el.etype, &xyz, total, m, csize, &mat, el.basis.as_ref(), &el.rule, &mut self.local)
                .map_err(fatal)?;
            if el.nb > 0 {
                let l = &mut self.local;
                condensate(m * el.nd, m * el.nb, &mut l.stiff, &mut l.force).map_err(fatal)?;
                l.n = m * el.nd;
            }
            ctx.default_update(&self.local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;

        // Work = 8 SUM(Values * RHS)
        let values: Vec<f64> = ctx.primary().map(|v| v.values.to_vec()).unwrap_or_default();
        let rhs = ctx.system().map(|s| s.rhs.clone()).unwrap_or_default();
        let work = 8.0 * values.iter().zip(&rhs).fold(0.0, |s, (v, r)| s + v * r);
        let t = 1.0e-1;
        // `1.828629366566552` has no `_dp` suffix: a default (single
        // precision) real.
        let ref_work =
            12.0 * (1.0 - (1.0_f64 / 3.0).powi(2)) * 1.0e9_f64.powi(2) / 7.0e10 * f64::from(1.828_629_4_f32) * t;
        (self.out)(&format!(" Relative energy error = {}", list_real(((ref_work - work).abs() / ref_work).sqrt())));
        Ok(())
    }
}
