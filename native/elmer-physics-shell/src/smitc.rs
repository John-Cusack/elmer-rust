//! `SmitcSolver` (Smitc.F90): the Reissner-Mindlin plate equations with
//! MITC shear interpolation on linear triangles and quadrilaterals. The
//! unknown is `Deflection` = (w, theta_x, theta_y); the plate lies in the
//! xy-plane (only `Nodes % x`, `Nodes % y` enter the shear interpolation).

use elmer_core::{ElemId, LocalSystem, ParamMap, Params, Value};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

pub const SMITC_PROCEDURE: (&str, &str) = ("Smitc", "SmitcSolver");

const CALLER: &str = "SmitcSolver";
const DOFS: usize = 3;

pub fn register_smitc(registry: &mut SolverRegistry) {
    registry.register(SMITC_PROCEDURE.0, SMITC_PROCEDURE.1, smitc_factory);
}

pub fn smitc_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(SmitcSolver))
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("SmitcSolver: {what}"))
}

/// `SmitcSolver_Init0`: `Variable = -dofs 3 Deflection`, `Time Derivative
/// Order = 2`, `Use Global Mass Matrix`, `Plate Solver`.
pub fn smitc_defaults(p: &mut ParamMap) {
    let mut add = |k: &str, v: Value| {
        if p.get(k).is_none() {
            p.insert(k, v);
        }
    };
    add("Variable", Value::Str("-dofs 3 Deflection".into()));
    add("Time derivative order", Value::Integer(vec![2]));
    add("Use Global Mass Matrix", Value::Logical(true));
    add("Plate Solver", Value::Logical(true));
}

pub struct SmitcSolver;

/// `SUM( a(1:n) * Basis(1:n) )`.
pub(crate) fn sum(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// `StrainEnergyDensity( A, B, C, m, n, s )`: `A += s C^T (B C)`.
pub(crate) fn strain_energy(
    a: &mut [f64],
    b: &[Vec<f64>],
    c: &[Vec<f64>],
    m: usize,
    n: usize,
    s: f64,
) {
    let mut t = vec![vec![0.0f64; n]; m];
    for i in 0..m {
        for j in 0..n {
            let mut acc = 0.0;
            for k in 0..m {
                acc += b[i][k] * c[k][j];
            }
            t[i][j] = acc;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0;
            for k in 0..m {
                acc += c[k][i] * t[k][j];
            }
            a[i * n + j] += s * acc;
        }
    }
}

/// `ShearCorrectionFactor( Kappa, h, x, y, n, alpha )` (its triangle edge
/// `x13 = x(1) - x(1)` included).
pub(crate) fn shear_correction(
    h: f64,
    x: &[f64],
    y: &[f64],
    stab: f64,
) -> Result<f64, SolverError> {
    let l = |dx: f64, dy: f64| (dx.powi(2) + dy.powi(2)).sqrt();
    let (alpha, hh) = match x.len() {
        3 => {
            let l21 = l(x[1] - x[0], y[1] - y[0]);
            let l32 = l(x[2] - x[1], y[2] - y[1]);
            // The Fortran's `x13 = x(1) - x(1)`, `y13 = y(1) - y(1)`.
            let l13 = 0.0f64;
            (0.20 * stab, l21.max(l32).max(l13))
        }
        4 => {
            let l21 = l(x[1] - x[0], y[1] - y[0]);
            let l32 = l(x[2] - x[1], y[2] - y[1]);
            let l43 = l(x[3] - x[2], y[3] - y[2]);
            let l14 = l(x[0] - x[3], y[0] - y[3]);
            (0.10 * stab, l21.max(l32).max(l43).max(l14))
        }
        n => {
            return Err(SolverError::fatal(
                "ShearCorrectionFactor",
                format!("Illegal number of nodes for Smitc elements: {n}"),
            ))
        }
    };
    Ok(h.powi(2) / (h.powi(2) + alpha * hh.powi(2)))
}

pub(crate) type J2 = ([[f64; 2]; 2], [[f64; 2]; 2]);

/// `Jacobi3( Jmat, invJ, detJ, x, y )`.
pub(crate) fn jacobi3(x: &[f64], y: &[f64]) -> J2 {
    let j = [[x[1] - x[0], y[1] - y[0]], [x[2] - x[0], y[2] - y[0]]];
    let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    (
        j,
        [
            [j[1][1] / det, -j[0][1] / det],
            [-j[1][0] / det, j[0][0] / det],
        ],
    )
}

/// `Jacobi4( Jmat, invJ, detJ, xi, eta, x, y )`.
pub(crate) fn jacobi4(xi: f64, eta: f64, x: &[f64], y: &[f64]) -> J2 {
    let dxi = [
        -(1.0 - eta) / 4.0,
        (1.0 - eta) / 4.0,
        (1.0 + eta) / 4.0,
        -(1.0 + eta) / 4.0,
    ];
    let deta = [
        -(1.0 - xi) / 4.0,
        -(1.0 + xi) / 4.0,
        (1.0 + xi) / 4.0,
        (1.0 - xi) / 4.0,
    ];
    let mut j = [[0.0f64; 2]; 2];
    for i in 0..4 {
        j[0][0] += dxi[i] * x[i];
        j[0][1] += dxi[i] * y[i];
        j[1][0] += deta[i] * x[i];
        j[1][1] += deta[i] * y[i];
    }
    let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    (
        j,
        [
            [j[1][1] / det, -j[0][1] / det],
            [-j[1][0] / det, j[0][0] / det],
        ],
    )
}

/// `Sdofs(k) = (Tau(1)*Jmat(1,c) + Tau(2)*Jmat(2,c)) / d` for column `c`.
pub(crate) fn tj(tau: [f64; 2], j: &[[f64; 2]; 2], c: usize, d: f64) -> f64 {
    (tau[0] * j[0][c] + tau[1] * j[1][c]) / d
}

/// `CovariantInterpolation`: the MITC shear strain `(2, 3n)`.
fn covariant_interpolation(x: &[f64], y: &[f64], u: f64, v: f64, kernel: bool) -> Vec<Vec<f64>> {
    let n = x.len();
    let nc = 3 * n;
    let mut sr = vec![vec![0.0f64; nc]; 2];
    let inv = match n {
        3 => {
            let (jm, inv) = jacobi3(x, y);
            // Edge 1-2: Tau = (1, 0).
            let tau = [1.0, 0.0];
            let mut sd = vec![0.0f64; nc];
            sd[1] = tj(tau, &jm, 0, 2.0);
            sd[2] = tj(tau, &jm, 1, 2.0);
            sd[4] = tj(tau, &jm, 0, 2.0);
            sd[5] = tj(tau, &jm, 1, 2.0);
            for j in 0..9 {
                sr[0][j] += (1.0 - v) * sd[j];
                sr[1][j] += u * sd[j];
            }
            let r2 = 2.0f64.sqrt();
            let tau = [-1.0 / r2, 1.0 / r2];
            let mut sd = vec![0.0f64; nc];
            sd[4] = tj(tau, &jm, 0, r2);
            sd[5] = tj(tau, &jm, 1, r2);
            sd[7] = tj(tau, &jm, 0, r2);
            sd[8] = tj(tau, &jm, 1, r2);
            for j in 0..9 {
                sr[0][j] += (-v) * sd[j];
                sr[1][j] += u * sd[j];
            }
            let tau = [0.0, -1.0];
            let mut sd = vec![0.0f64; nc];
            sd[1] = tj(tau, &jm, 0, 2.0);
            sd[2] = tj(tau, &jm, 1, 2.0);
            sd[7] = tj(tau, &jm, 0, 2.0);
            sd[8] = tj(tau, &jm, 1, 2.0);
            for j in 0..9 {
                sr[0][j] += (-v) * sd[j];
                sr[1][j] += (-1.0 + u) * sd[j];
            }
            inv
        }
        4 => {
            if kernel {
                let ab = [
                    [(1.0 - v) / 4.0, 0.0],
                    [0.0, (1.0 + u) / 4.0],
                    [(1.0 + v) / 4.0, 0.0],
                    [0.0, (1.0 - u) / 4.0],
                ];
                let a1 = 0.25 * (-x[0] + x[1] + x[2] - x[3]);
                let a2 = 0.25 * (-x[0] - x[1] + x[2] + x[3]);
                let a12 = 0.25 * (x[0] - x[1] + x[2] - x[3]);
                let b1 = 0.25 * (-y[0] + y[1] + y[2] - y[3]);
                let b2 = 0.25 * (-y[0] - y[1] + y[2] + y[3]);
                let b12 = 0.25 * (y[0] - y[1] + y[2] - y[3]);
                // ShearRef(:, col) +/-= abf0basis(:, k), 1-based col and k.
                let addp = |sr: &mut Vec<Vec<f64>>, col: usize, k: usize, sgn: f64| {
                    for r in 0..2 {
                        if sgn > 0.0 {
                            sr[r][col - 1] += ab[k - 1][r];
                        } else {
                            sr[r][col - 1] -= ab[k - 1][r];
                        }
                    }
                };
                addp(&mut sr, 1, 1, 1.0);
                addp(&mut sr, 4, 1, -1.0);
                addp(&mut sr, 4, 2, 1.0);
                addp(&mut sr, 7, 2, -1.0);
                addp(&mut sr, 7, 3, -1.0);
                addp(&mut sr, 10, 3, 1.0);
                addp(&mut sr, 10, 4, -1.0);
                addp(&mut sr, 1, 4, 1.0);
                // ShearRef(:, col) += c * abf0basis(:, k).
                let mut add = |col: usize, c: f64, k: usize| {
                    for r in 0..2 {
                        sr[r][col - 1] += c * ab[k - 1][r];
                    }
                };
                add(2, a1 - a12, 1);
                add(3, b1 - b12, 1);
                add(5, a1 - a12, 1);
                add(6, b1 - b12, 1);
                add(5, a2 + a12, 2);
                add(6, b2 + b12, 2);
                add(8, a2 + a12, 2);
                add(9, b2 + b12, 2);
                add(8, a1 + a12, 3);
                add(9, b1 + b12, 3);
                add(11, a1 + a12, 3);
                add(12, b1 + b12, 3);
                add(2, a1 - a12, 4);
                add(3, b1 - b12, 4);
                add(5, a1 + a2, 4);
                add(6, b1 + b2, 4);
                add(8, a2 - a1, 4);
                add(9, b2 - b1, 4);
                add(11, -a1 - a12, 4);
                add(12, -b1 - b12, 4);
            } else {
                // (xi, eta, Tau, Sdofs columns, row, factor)
                let edges: [(f64, f64, [f64; 2], [usize; 4], usize); 4] = [
                    (0.0, -1.0, [1.0, 0.0], [2, 3, 5, 6], 0),
                    (1.0, 0.0, [0.0, 1.0], [5, 6, 8, 9], 1),
                    (0.0, 1.0, [-1.0, 0.0], [8, 9, 11, 12], 0),
                    (-1.0, 0.0, [0.0, -1.0], [2, 3, 11, 12], 1),
                ];
                for (k, (xi, eta, tau, cols, row)) in edges.into_iter().enumerate() {
                    let (jm, _) = jacobi4(xi, eta, x, y);
                    let mut sd = vec![0.0f64; nc];
                    sd[cols[0] - 1] = tj(tau, &jm, 0, 1.0);
                    sd[cols[1] - 1] = tj(tau, &jm, 1, 1.0);
                    sd[cols[2] - 1] = tj(tau, &jm, 0, 1.0);
                    sd[cols[3] - 1] = tj(tau, &jm, 1, 1.0);
                    let f = match k {
                        0 => (1.0 - v) / 4.0,
                        1 => (1.0 + u) / 4.0,
                        2 => (-1.0 - v) / 4.0,
                        _ => (-1.0 + u) / 4.0,
                    };
                    for j in 0..12 {
                        sr[row][j] += f * sd[j];
                    }
                }
            }
            jacobi4(u, v, x, y).1
        }
        _ => return vec![vec![0.0; nc]; 2],
    };
    // ShearStrain = -MATMUL( invJ, ShearRef ).
    (0..2)
        .map(|i| {
            (0..nc)
                .map(|j| -(0.0 + inv[i][0] * sr[0][j] + inv[i][1] * sr[1][j]))
                .collect()
        })
        .collect()
}

struct ElementData {
    load: Vec<f64>,
    poisson: Vec<f64>,
    density: Vec<f64>,
    thickness: Vec<f64>,
    young: Vec<f64>,
    tension: Vec<f64>,
    spring: Vec<f64>,
    damping: Vec<f64>,
}

impl Solver for SmitcSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        for k in ["material index", "body force index"] {
            if p.get(k).is_some() {
                return Err(unsupported(k));
            }
        }
        if p.logical("hole correction") == Some(true) {
            return Err(unsupported("Hole Correction (perforated plates)"));
        }
        if p.logical("harmonic analysis") == Some(true) {
            return Err(unsupported("harmonic analysis (damping matrix)"));
        }
        if p.logical("harmonic mode") == Some(true) && p.logical("eigen analysis") != Some(true) {
            return Err(unsupported("Harmonic Mode outside an eigen analysis"));
        }
        if ctx.comm.size() > 1 {
            return Err(unsupported("more than one partition"));
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        transient: bool,
    ) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("transient plates (Default2ndOrderTime)"));
        }
        ctx.info(
            CALLER,
            "--------------------------------------------------",
            4,
        );
        ctx.info(
            CALLER,
            "Solving the Reissner-Mindlin equations for plates",
            4,
        );
        ctx.info(
            CALLER,
            "--------------------------------------------------",
            4,
        );
        let p = ctx.params().clone();
        let alpha = ctx.creal(&p, "shear relaxation alpha")?.unwrap_or(1.0);
        let eigen = p.logical("eigen analysis") == Some(true);
        let kernel = p.logical("kernel interpolation") == Some(true);
        let mut grav = 0.0;
        if p.logical("use gravity") == Some(true) {
            match ctx.constants().get("gravity") {
                Some(Value::Real { data, .. }) if !data.is_empty() => {
                    grav = data[data.len() - 1].abs()
                }
                _ => {
                    return Err(SolverError::fatal(
                        CALLER,
                        "Gravity requested but not given as constant!",
                    ));
                }
            }
        }
        let max_iter = p
            .integer("nonlinear system max iterations")
            .unwrap_or(1)
            .max(1);
        let name = ctx
            .variable_name()
            .ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?
            .to_string();
        let id = ctx
            .sim
            .vars
            .id(&name)
            .ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        for _ in 0..max_iter {
            ctx.default_initialize()?;
            for &e in &elems {
                let n = ctx.mesh.element(e).nodes.len();
                let d = self.element_data(ctx, e, n)?;
                let (ls, mass, damp) = local_matrix(ctx, e, &d, alpha, kernel, grav)?;
                ctx.default_update(&ls, e)?;
                if eigen {
                    ctx.default_update_mass(&mass, e)?;
                    ctx.default_update_damp(&damp, e)?;
                }
            }
            ctx.default_finish_bulk_assembly()?;
            ctx.default_finish_boundary_assembly()?;
            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.sim.vars.by_id(id).conv.nonlin_converged == 1 {
                break;
            }
        }
        Ok(())
    }
}

impl SmitcSolver {
    fn element_data(
        &self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        n: usize,
    ) -> Result<ElementData, SolverError> {
        let bp = ctx.body_params(e);
        let bf = bp.as_ref().and_then(|b| b.body_force);
        let mat = bp.as_ref().and_then(|b| b.material);
        let get = |s: Option<&ParamMap>, k: &str| -> Result<Option<Vec<f64>>, SolverError> {
            match s {
                Some(s) => ctx.real_at_nodes(s, k, e),
                None => Ok(None),
            }
        };
        let zero = || vec![0.0; n];
        let mut load = zero();
        if let Some(l) = get(bf, "pressure")? {
            load = l;
            if let Some(l2) = get(bf, "pressure b")? {
                for (a, b) in load.iter_mut().zip(&l2) {
                    *a += b;
                }
                if let Some(l3) = get(bf, "pressure c")? {
                    for (a, b) in load.iter_mut().zip(&l3) {
                        *a += b;
                    }
                }
            }
        }
        // `Damping` (body force, else material) only enters DAMP.
        let damping = match get(bf, "damping")? {
            Some(v) => v,
            None => get(mat, "damping")?.unwrap_or_else(zero),
        };
        let spring = match get(bf, "spring")? {
            Some(v) => v,
            None => get(mat, "spring")?.unwrap_or_else(zero),
        };
        Ok(ElementData {
            load,
            poisson: get(mat, "poisson ratio")?.unwrap_or_else(zero),
            density: get(mat, "density")?.unwrap_or_else(zero),
            thickness: get(mat, "thickness")?.unwrap_or_else(zero),
            young: get(mat, "youngs modulus")?.unwrap_or_else(zero),
            tension: get(mat, "tension")?.unwrap_or_else(zero),
            spring,
            damping,
        })
    }
}

/// `LocalMatrix`: the element stiffness and force, its mass and damping.
fn local_matrix(
    ctx: &SolverCtx<'_>,
    e: ElemId,
    d: &ElementData,
    alpha: f64,
    kernel: bool,
    grav: f64,
) -> Result<(LocalSystem, Vec<f64>, Vec<f64>), SolverError> {
    let el = ctx.mesh.element(e);
    let n = el.nodes.len();
    let nd = DOFS * n;
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let x: Vec<f64> = xyz.iter().map(|c| c[0]).collect();
    let y: Vec<f64> = xyz.iter().map(|c| c[1]).collect();
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
    let mut stiff = vec![0.0f64; nd * nd];
    let mut mass = vec![0.0f64; nd * nd];
    let mut damp = vec![0.0f64; nd * nd];
    let mut force = vec![0.0f64; nd];
    let mut curv = vec![vec![0.0f64; nd]; 3];
    let mut ev = elmer_core::BasisEval::default();
    for qp in fe.quadrature.gauss_points(el.etype, None) {
        fe.basis
            .eval(el.etype, &xyz, qp.uvw, &mut ev)
            .map_err(|m| SolverError::fatal("ElementInfo", m))?;
        let b = &ev.basis[..n];
        let s = qp.weight * ev.det_j;
        let pressure = sum(&d.load, b);
        let h = sum(&d.thickness, b);
        let winkler = sum(&d.spring, b);
        let damp_coef = sum(&d.damping, b);
        // IsotropicElasticity.
        let eu = sum(&d.young, b);
        let pu = sum(&d.poisson, b);
        let tu = sum(&d.thickness, b);
        let gu = eu / (2.0 * (1.0 + pu));
        let mut em = vec![vec![0.0f64; 3]; 3];
        em[0][0] = 1.0;
        em[0][1] = pu;
        em[1][0] = pu;
        em[1][1] = 1.0;
        em[2][2] = (1.0 - pu) / 2.0;
        for r in em.iter_mut() {
            for v in r.iter_mut() {
                *v = *v * eu * tu.powi(3) / (12.0 * (1.0 - pu.powi(2)));
            }
        }
        let gm = vec![vec![gu * tu, 0.0], vec![0.0, gu * tu]];
        let rho = sum(&d.density, b);
        let t = sum(&d.tension, b) * h;
        let tm = vec![vec![t, 0.0], vec![0.0, sum(&d.tension, b) * h]];
        for p in 0..n {
            let db = ev.dbasis_dx[p];
            curv[0][3 * p + 1] = db[0];
            curv[1][3 * p + 2] = db[1];
            curv[2][3 * p + 1] = db[1];
            curv[2][3 * p + 2] = db[0];
        }
        strain_energy(&mut stiff, &em, &curv, 3, nd, s);
        let mut shear = covariant_interpolation(&x, &y, qp.uvw[0], qp.uvw[1], kernel);
        let kappa = shear_correction(h, &x, &y, alpha)?;
        if !kernel {
            for p in 0..n {
                shear[0][3 * p] = ev.dbasis_dx[p][0];
                shear[1][3 * p] = ev.dbasis_dx[p][1];
            }
        }
        strain_energy(&mut stiff, &gm, &shear, 2, nd, kappa * s);
        let mut shear = vec![vec![0.0f64; nd]; 2];
        for p in 0..n {
            shear[0][3 * p] = ev.dbasis_dx[p][0];
            shear[1][3 * p] = ev.dbasis_dx[p][1];
        }
        strain_energy(&mut stiff, &tm, &shear, 2, nd, s);
        for p in 0..n {
            for q in 0..n {
                stiff[(DOFS * p) * nd + DOFS * q] += winkler * b[p] * b[q] * s;
            }
        }
        for p in 0..n {
            force[DOFS * p] += pressure * b[p] * s;
        }
        for p in 0..n {
            force[DOFS * p] -= grav * rho * h * b[p] * s;
        }
        for p in 0..n {
            for q in 0..n {
                mass[(DOFS * p) * nd + DOFS * q] += rho * h * b[p] * b[q] * s;
            }
        }
        for p in 0..n {
            for q in 0..n {
                damp[(DOFS * p) * nd + DOFS * q] += damp_coef * b[p] * b[q] * s;
            }
        }
    }
    Ok((
        LocalSystem {
            n: nd,
            stiff,
            mass: None,
            damp: None,
            force,
        },
        mass,
        damp,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The square `[0, 2]^2` maps one to one onto the reference square.
    #[test]
    fn jacobi4_of_the_doubled_square_is_the_identity() {
        let (x, y) = ([0.0, 2.0, 2.0, 0.0], [0.0, 0.0, 2.0, 2.0]);
        let (j, inv) = jacobi4(0.3, -0.2, &x, &y);
        assert_eq!(j, [[1.0, 0.0], [0.0, 1.0]]);
        assert_eq!(inv, [[1.0, 0.0], [0.0, 1.0]]);
    }

    /// Kappa = h^2 / (h^2 + alpha l^2) with the longest edge; the triangle
    /// leaves out its third edge, as the Fortran does.
    #[test]
    fn shear_correction_uses_the_longest_counted_edge() {
        let h2 = 0.1f64.powi(2);
        let k = shear_correction(0.1, &[0.0, 1.0, 1.0, 0.0], &[0.0, 0.0, 1.0, 1.0], 1.0).unwrap();
        assert_eq!(k, h2 / (h2 + 0.1));
        // Edges 1-2 and 2-3 have length 1; 3-1 (sqrt 2) is not counted.
        let k = shear_correction(0.1, &[0.0, 1.0, 1.0], &[0.0, 0.0, 1.0], 1.0).unwrap();
        assert_eq!(k, h2 / (h2 + 0.2));
        assert!(shear_correction(0.1, &[0.0; 6], &[0.0; 6], 1.0).is_err());
    }
}
