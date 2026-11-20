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
