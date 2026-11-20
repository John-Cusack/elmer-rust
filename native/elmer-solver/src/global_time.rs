//! `Use Global Mass Matrix`: the element mass matrices go to the global
//! `Matrix % MassValues` (`DefaultUpdateMass` in place of the local
//! `Default1stOrderTime`), and `DefaultFinishAssembly` adds the first-order
//! time discretisation to the assembled system (`Default1stOrderTimeGlobal`
//! -> `Add1stOrderTime_CRS`, SolverBasics.F90:417) with `BDF_CRS`,
//! `VBDF_CRS` or `NewmarkBeta_CRS` (TimeIntegrate.F90), in their arithmetic
//! order. MATRIX_CRS only, consistent (not lumped) mass.

use elmer_core::CsrMatrix;
use elmer_timeint::Scheme;

use crate::error::SolverError;

/// The previous solutions, newest first (`PrevValues(:, 1..)`), and the
/// previous step sizes, newest first (`Timestep size` `PrevValues`).
pub(crate) struct History<'a> {
    pub prev: &'a [Vec<f64>],
    pub dts: &'a [f64],
}

/// `Add1stOrderTime_CRS(Matrix, Force, dt, Solver)` for `order =
/// MIN(DoneTime, Order)` previous solutions.
pub(crate) fn add_1st_order_time_crs(
    a: &mut CsrMatrix<f64>,
    mass: &CsrMatrix<f64>,
    force: &mut [f64],
    dt: f64,
    scheme: Scheme,
    order: usize,
    h: &History<'_>,
) -> Result<(), SolverError> {
    if order == 0 {
        return Ok(());
    }
    if h.prev.len() < order {
        return Err(SolverError::fatal(
            "Add1stOrderTime_CRS",
            "the variable has fewer previous solutions than the time order",
        ));
    }
    match scheme {
        Scheme::None => Ok(()),
        Scheme::Bdf(_) => {
            let mut dts = vec![dt];
            let mut constant = true;
            for i in 1..order {
                let d = h.dts.get(i - 1).copied().unwrap_or(dt);
                if (d - dt).abs() > 1.0e-6 * dt {
                    constant = false;
                }
                dts.push(d);
            }
            if constant {
                bdf_crs(dt, a, mass, force, h.prev, order)
            } else {
                vbdf_crs(&dts, a, mass, force, h.prev, order);
                Ok(())
            }
        }
        Scheme::NewmarkBeta(beta) => {
            newmark_beta_crs(dt, a, mass, force, &h.prev[0], beta);
            Ok(())
        }
        other => Err(SolverError::Unsupported(format!(
            "Use Global Mass Matrix with the time scheme {other:?}"
        ))),
    }
}

/// `BDF_CRS`, consistent mass.
fn bdf_crs(
    dt: f64,
    a: &mut CsrMatrix<f64>,
    mass: &CsrMatrix<f64>,
    force: &mut [f64],
    prev: &[Vec<f64>],
    order: usize,
) -> Result<(), SolverError> {
    let mut c = [0.0f64; 6];
    match order {
        1 => {
            c[0] = 1.0;
            c[1] = -1.0;
        }
        2 => {
            c[0] = 1.5;
            c[1] = -2.0;
            c[2] = 0.5;
        }
        3 => {
            c[0] = 11.0 / 6.0;
            c[1] = -3.0;
            c[2] = 3.0 / 2.0;
            c[3] = -1.0 / 3.0;
        }
        4 => {
            c[0] = 25.0 / 12.0;
            c[1] = -4.0;
            c[2] = 3.0;
            c[3] = -4.0 / 3.0;
            c[4] = 1.0 / 4.0;
        }
        5 => {
            c[0] = 137.0 / 60.0;
            c[1] = -5.0;
            c[2] = 5.0;
            c[3] = -10.0 / 3.0;
            c[4] = 5.0 / 4.0;
            c[5] = -1.0 / 5.0;
        }
        _ => {
            return Err(SolverError::fatal(
                "BDF_CRS",
                "Constant timestep BDF implemented only to order 5",
            ))
        }
    }
    for x in &mut c {
        *x /= dt;
    }
    for i in 0..a.nrows {
        let mut mu = 0.0;
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            let col = a.col_idx[j] as usize;
            for k in 0..order {
                let uj = prev[k][col];
                mu -= mass.values[j] * c[k + 1] * uj;
            }
            a.values[j] += c[0] * mass.values[j];
        }
        force[i] += mu;
    }
    Ok(())
}

/// `VBDF_CRS`, consistent mass, to order 3.
fn vbdf_crs(
    d: &[f64],
    a: &mut CsrMatrix<f64>,
    mass: &CsrMatrix<f64>,
    force: &mut [f64],
    prev: &[Vec<f64>],
    order: usize,
) {
    let mut c = [0.0f64; 4];
    c[0] = 1.0 / d[0];
    c[1] = -1.0 / d[0];
    if order >= 2 {
        c[0] += 1.0 / (d[0] + d[1]);
        c[1] -= (1.0 + d[0] / d[1]) / (d[0] + d[1]);
        c[2] = (d[0] / d[1]) / (d[0] + d[1]);
    }
    if order >= 3 {
        let s = d[0] + d[1] + d[2];
        c[0] += 1.0 / s;
        c[1] -= (1.0 + d[0] / d[1] * (1.0 + (d[0] + d[1]) / (d[1] + d[2]))) / s;
        c[2] += (d[0] / d[1] * (1.0 + (d[0] + d[1]) / (d[1] + d[2]))
            + d[0] / d[2] * (d[0] + d[1]) / (d[1] + d[2]))
            / s;
        c[3] = -(d[0] / d[2]) * (d[0] + d[1]) / (d[1] + d[2]) / s;
    }
    let kmax = order.min(3);
    for i in 0..a.nrows {
        let mut mu = 0.0;
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            let col = a.col_idx[j] as usize;
            for k in 0..kmax {
                let uj = prev[k][col];
                mu -= mass.values[j] * c[k + 1] * uj;
            }
        }
        force[i] += mu;
    }
    for (s, m) in a.values.iter_mut().zip(&mass.values) {
        *s += c[0] * m;
    }
}

/// `NewmarkBeta_CRS`, consistent mass, no bulk residual.
fn newmark_beta_crs(
    dt: f64,
    a: &mut CsrMatrix<f64>,
    mass: &CsrMatrix<f64>,
    force: &mut [f64],
    prev: &[f64],
    beta: f64,
) {
    for i in 0..a.nrows {
        let (mut su, mut mu) = (0.0, 0.0);
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            let uj = prev[a.col_idx[j] as usize];
            su += a.values[j] * uj;
            mu += mass.values[j] * uj;
        }
        force[i] = force[i] - (1.0 - beta) * su + (1.0 / dt) * mu;
    }
    for (s, m) in a.values.iter_mut().zip(&mass.values) {
        *s = beta * *s + (1.0 / dt) * m;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn system() -> (CsrMatrix<f64>, CsrMatrix<f64>) {
        // 2x2 dense pattern.
        let a = CsrMatrix {
            nrows: 2,
            ncols: 2,
            row_ptr: vec![0, 2, 4],
            col_idx: vec![0, 1, 0, 1],
            values: vec![2.0, -1.0, -1.0, 2.0],
            diag: None,
        };
        let mut m = a.clone();
        m.values = vec![2.0, 1.0, 1.0, 2.0];
        (a, m)
    }

    /// BDF(1) and implicit Euler are `A + M/dt`, `b + M x_prev / dt`.
    #[test]
    fn bdf1_and_implicit_euler_add_m_over_dt() {
        let prev = vec![vec![1.0, 3.0]];
        for scheme in [Scheme::Bdf(1), Scheme::NewmarkBeta(1.0)] {
            let (mut a, m) = system();
            let mut f = vec![0.5, 0.5];
            let h = History { prev: &prev, dts: &[] };
            add_1st_order_time_crs(&mut a, &m, &mut f, 0.5, scheme, 1, &h).unwrap();
            assert_eq!(a.values, vec![6.0, 1.0, 1.0, 6.0]);
            assert_eq!(f, vec![0.5 + (2.0 + 3.0) / 0.5, 0.5 + (1.0 + 6.0) / 0.5]);
        }
    }

    /// Constant-step BDF(2) and the variable-step form with equal steps
    /// give the same system up to rounding.
    #[test]
    fn bdf2_constant_and_variable_agree() {
        let prev = vec![vec![1.0, 3.0], vec![0.5, 2.0]];
        let (mut a1, m) = system();
        let mut f1 = vec![0.0; 2];
        bdf_crs(0.25, &mut a1, &m, &mut f1, &prev, 2).unwrap();
        let (mut a2, _) = system();
        let mut f2 = vec![0.0; 2];
        vbdf_crs(&[0.25, 0.25], &mut a2, &m, &mut f2, &prev, 2);
        for (x, y) in a1.values.iter().zip(&a2.values).chain(f1.iter().zip(&f2)) {
            assert!((x - y).abs() < 1e-12, "{x} vs {y}");
        }
    }
}

/// `Time2ndOrder_CRS(dt, Matrix, Force, Alpha, Beta, PrevValues)`
/// (TimeIntegrate.F90): the Bossak (`param_b < 0`) or generalized-alpha
/// (`param_b` = rho_inf) terms of the global mass and damping matrices,
/// added to the stiffness and the force, with `X, V, A, A2` the previous
/// values' columns 3, 4, 5 and 7 (5 for Bossak).
#[allow(clippy::too_many_arguments)]
pub(crate) fn time_2nd_order_crs(
    dt: f64,
    a: &mut CsrMatrix<f64>,
    mass: Option<&CsrMatrix<f64>>,
    damp: Option<&CsrMatrix<f64>>,
    force: &mut [f64],
    param_a: f64,
    param_b: f64,
    prev: &[Vec<f64>],
) -> Result<(), SolverError> {
    if prev.len() < 7 {
        return Err(SolverError::fatal("Time2ndOrder_CRS", "the variable lacks its seven previous-value columns"));
    }
    let (x, v, acc) = (&prev[2], &prev[3], &prev[4]);
    let (alpha, beta, gamma, delta, a2);
    if param_b >= 0.0 {
        // rho_inf for the generalized-alpha method.
        let b = param_b;
        alpha = (2.0 * b - 1.0) / (b + 1.0);
        beta = 1.0 / ((b + 1.0) * (b + 1.0));
        gamma = (3.0 - b) / (2.0 * (b + 1.0));
        delta = b / (b + 1.0);
        a2 = &prev[6];
    } else {
        // Bossak.
        alpha = param_a;
        gamma = 0.5 - alpha;
        beta = (1.0 - alpha) * (1.0 - alpha) / 4.0;
        delta = 0.0;
        a2 = &prev[4];
    }
    let c1 = (1.0 - alpha) / (1.0 - delta) / (beta * (dt * dt));
    let c2 = (1.0 - alpha) / (1.0 - delta) / (beta * dt);
    let c3 = delta / (1.0 - delta);
    let c4 = -((1.0 - alpha) * (1.0 - 1.0 / (2.0 * beta)) + alpha) / (1.0 - delta);
    let c5 = gamma / (beta * dt);
    let c6 = gamma / beta - 1.0;
    let c7 = -((1.0 - gamma) + gamma * (1.0 - 1.0 / (2.0 * beta))) * dt;
    let n = a.nrows;
    if let Some(m) = mass {
        for i in 0..n {
            let mut s = 0.0;
            for j in m.row_ptr[i]..m.row_ptr[i + 1] {
                let k = m.col_idx[j] as usize;
                s += m.values[j] * (c1 * x[k] + c2 * v[k] + c3 * acc[k] + c4 * a2[k]);
            }
            force[i] += s;
        }
        for (st, mv) in a.values.iter_mut().zip(&m.values) {
            *st += c1 * mv;
        }
    }
    if let Some(d) = damp {
        for i in 0..n {
            let mut s = 0.0;
            for j in d.row_ptr[i]..d.row_ptr[i + 1] {
                let k = d.col_idx[j] as usize;
                s += d.values[j] * (c5 * x[k] + c6 * v[k] + c7 * a2[k]);
            }
            force[i] += s;
        }
        for (st, dv) in a.values.iter_mut().zip(&d.values) {
            *st += c5 * dv;
        }
    }
    Ok(())
}

#[cfg(test)]
mod second_order_tests {
    use super::*;

    fn one(v: f64) -> CsrMatrix<f64> {
        CsrMatrix { nrows: 1, ncols: 1, row_ptr: vec![0, 1], col_idx: vec![0], values: vec![v], diag: Some(vec![0]) }
    }

    /// Bossak on `m x'' + d x' + k x = f`, one dof: the stiffness gains
    /// `c1 m + c5 d` and the force the history terms, with Elmer's
    /// coefficients.
    #[test]
    fn bossak_terms_of_a_single_dof() {
        let (m, d, k, f, dt, al) = (2.0, 0.5, 3.0, 1.0, 0.1, -0.05);
        let prev = vec![vec![0.0], vec![0.0], vec![0.4], vec![0.3], vec![0.2], vec![0.0], vec![9.0]];
        let mut a = one(k);
        let mut rhs = vec![f];
        time_2nd_order_crs(dt, &mut a, Some(&one(m)), Some(&one(d)), &mut rhs, al, -1.0, &prev).unwrap();
        let (gamma, beta) = (0.5 - al, (1.0 - al) * (1.0 - al) / 4.0);
        let c1 = (1.0 - al) / (beta * dt * dt);
        let c2 = (1.0 - al) / (beta * dt);
        let c4 = -((1.0 - al) * (1.0 - 1.0 / (2.0 * beta)) + al);
        let c5 = gamma / (beta * dt);
        let c6 = gamma / beta - 1.0;
        let c7 = -((1.0 - gamma) + gamma * (1.0 - 1.0 / (2.0 * beta))) * dt;
        let want_k = k + c1 * m + c5 * d;
        // A2 is column 5 for Bossak, not column 7.
        let want_f = f + m * (c1 * 0.4 + c2 * 0.3 + c4 * 0.2) + d * (c5 * 0.4 + c6 * 0.3 + c7 * 0.2);
        assert!((a.values[0] - want_k).abs() < 1e-9 * want_k.abs());
        assert!((rhs[0] - want_f).abs() < 1e-9 * want_f.abs());
    }
}
