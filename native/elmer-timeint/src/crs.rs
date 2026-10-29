//! Global (`_CRS`) time discretisation: `RungeKutta_CRS`, `NewmarkBeta_CRS`,
//! `BDF_CRS`, `VBDF_CRS`, `FractionalStep_CRS` (TimeIntegrate.F90:~700-1090)
//! and `Add1stOrderTime_CRS` (SolverBasics.F90:~420).
//!
//! The stiffness matrix is the `CsrMatrix` (`Matrix % Values`). The other
//! `Matrix_t` fields the routines read arrive in [`CrsTimeData`]; mass,
//! damping and the lumped diagonal share the stiffness sparsity pattern.
//!
//! Two Fortran branches read a stale loop index after the row loop
//! (`MassL(i)` in the lumped `BulkResidual` cases of `RungeKutta_CRS` and
//! `NewmarkBeta_CRS`). This port adds `MassL(row)` to each diagonal entry,
//! which is what the neighbouring branches compute.

use crate::error::TimeError;
use crate::local::{vbdf_coeffs, Aux1st};
use crate::scheme::{FsCoeffs, Scheme};
use elmer_core::{CsrMatrix, Logger};

/// The `Matrix_t` fields besides `Values` that time integration touches.
#[derive(Clone, Copy, Debug, Default)]
pub struct CrsTimeData<'a> {
    /// `MassValues`, one per stored entry.
    pub mass: Option<&'a [f64]>,
    /// `DampValues`, one per stored entry.
    pub damp: Option<&'a [f64]>,
    /// `MassValuesLumped`, one per row.
    pub mass_lumped: Option<&'a [f64]>,
    /// `Matrix % Lumped`. `BDF_CRS` and `VBDF_CRS` branch on this flag; the
    /// other routines branch on `mass_lumped.is_some()`, as in the Fortran.
    pub lumped: bool,
    /// `BulkResidual` (`A x - f` of the previous iteration).
    pub bulk_residual: Option<&'a [f64]>,
}

pub(crate) fn diag_check(a: &CsrMatrix, nforce: usize, caller: &'static str) -> Result<(), TimeError> {
    if a.row_ptr.len() != a.nrows + 1
        || a.values.len() != *a.row_ptr.last().unwrap_or(&0)
        || a.col_idx.len() != a.values.len()
        || nforce < a.nrows
    {
        return Err(TimeError::Dim(format!("{caller}: inconsistent CSR / force sizes")));
    }
    Ok(())
}

/// Index of `A[i,i]` in `values`: `diag` when present, else a search of the
/// sorted row.
fn diag_index(a: &CsrMatrix, i: usize) -> Result<usize, TimeError> {
    if let Some(d) = &a.diag {
        return Ok(d[i]);
    }
    let (lo, hi) = (a.row_ptr[i], a.row_ptr[i + 1]);
    a.col_idx[lo..hi]
        .binary_search(&(i as u32))
        .map(|p| lo + p)
        .map_err(|_| TimeError::Dim(format!("row {i} has no diagonal entry")))
}

fn mass_vals<'a>(d: &CrsTimeData<'a>, a: &CsrMatrix, caller: &'static str) -> Result<&'a [f64], TimeError> {
    let m = d.mass.ok_or(TimeError::MissingMass(caller))?;
    if m.len() != a.values.len() {
        return Err(TimeError::Dim(format!("{caller}: MassValues length")));
    }
    Ok(m)
}

fn lumped_vals<'a>(d: &CrsTimeData<'a>, a: &CsrMatrix, caller: &'static str) -> Result<&'a [f64], TimeError> {
    let m = d.mass_lumped.ok_or(TimeError::MissingLumpedMass(caller))?;
    if m.len() < a.nrows {
        return Err(TimeError::Dim(format!("{caller}: MassValuesLumped length")));
    }
    Ok(m)
}

fn row_dot(a: &CsrMatrix, vals: &[f64], i: usize, x: &[f64]) -> f64 {
    let mut s = 0.0;
    for j in a.row_ptr[i]..a.row_ptr[i + 1] {
        s += vals[j] * x[a.col_idx[j] as usize];
    }
    s
}

/// `RungeKutta_CRS`.
pub fn runge_kutta_crs(
    dt: f64,
    a: &mut CsrMatrix,
    d: &CrsTimeData<'_>,
    force: &mut [f64],
    prev: &[f64],
    cur: &[f64],
) -> Result<(), TimeError> {
    diag_check(a, force.len(), "RungeKutta_CRS")?;
    let n = a.nrows;
    if let Some(res) = d.bulk_residual {
        if d.mass_lumped.is_some() {
            let ml = lumped_vals(d, a, "RungeKutta_CRS")?;
            for i in 0..n {
                force[i] = res[i] + (1.0 / dt) * ml[i] * prev[i];
            }
            a.values.fill(0.0);
            for i in 0..n {
                let j = diag_index(a, i)?;
                a.values[j] += (1.0 / dt) * ml[i];
            }
        } else {
            let mass = mass_vals(d, a, "RungeKutta_CRS")?;
            for i in 0..n {
                let mu = row_dot(a, mass, i, prev);
                force[i] = res[i] + (1.0 / dt) * mu;
            }
            for (s, m) in a.values.iter_mut().zip(mass) {
                *s = (1.0 / dt) * m;
            }
        }
    } else if d.mass_lumped.is_some() {
        let ml = lumped_vals(d, a, "RungeKutta_CRS")?;
        for i in 0..n {
            let su = row_dot(a, &a.values, i, cur);
            let mu = ml[i] * prev[i];
            force[i] = force[i] - su + (1.0 / dt) * mu;
        }
        a.values.fill(0.0);
        for i in 0..n {
            let j = diag_index(a, i)?;
            a.values[j] += (1.0 / dt) * ml[i];
        }
    } else {
        let mass = mass_vals(d, a, "RungeKutta_CRS")?;
        for i in 0..n {
            let mut su = 0.0;
            let mut mu = 0.0;
            for j in a.row_ptr[i]..a.row_ptr[i + 1] {
                let k = a.col_idx[j] as usize;
                su += a.values[j] * cur[k];
                mu += mass[j] * prev[k];
            }
            force[i] = force[i] - su + (1.0 / dt) * mu;
        }
        for (s, m) in a.values.iter_mut().zip(mass) {
            *s = (1.0 / dt) * m;
        }
    }
    Ok(())
}

/// `NewmarkBeta_CRS`: implicit Euler for `beta = 1`, Crank-Nicolson for `0.5`.
pub fn newmark_beta_crs(
    dt: f64,
    a: &mut CsrMatrix,
    d: &CrsTimeData<'_>,
    force: &mut [f64],
    prev: &[f64],
    beta: f64,
) -> Result<(), TimeError> {
    diag_check(a, force.len(), "NewmarkBeta_CRS")?;
    let n = a.nrows;
    if let Some(res) = d.bulk_residual {
        if d.mass_lumped.is_some() {
            let ml = lumped_vals(d, a, "NewmarkBeta_CRS")?;
            for i in 0..n {
                force[i] = beta * force[i] - (1.0 - beta) * res[i] + (1.0 / dt) * ml[i] * prev[i];
            }
            for s in a.values.iter_mut() {
                *s = beta * *s;
            }
            for i in 0..n {
                let j = diag_index(a, i)?;
                a.values[j] += (1.0 / dt) * ml[i];
            }
        } else {
            let mass = mass_vals(d, a, "NewmarkBeta_CRS")?;
            for i in 0..n {
                let mu = row_dot(a, mass, i, prev);
                force[i] = beta * force[i] - (1.0 - beta) * res[i] + (1.0 / dt) * mu;
            }
            for (s, m) in a.values.iter_mut().zip(mass) {
                *s = beta * *s + (1.0 / dt) * m;
            }
        }
    } else if d.mass_lumped.is_some() {
        let ml = lumped_vals(d, a, "NewmarkBeta_CRS")?;
        for i in 0..n {
            let su = row_dot(a, &a.values, i, prev);
            let mu = ml[i] * prev[i];
            force[i] = force[i] - (1.0 - beta) * su + (1.0 / dt) * mu;
        }
        for s in a.values.iter_mut() {
            *s = beta * *s;
        }
        for i in 0..n {
            let j = diag_index(a, i)?;
            a.values[j] += (1.0 / dt) * ml[i];
        }
    } else {
        let mass = mass_vals(d, a, "NewmarkBeta_CRS")?;
        for i in 0..n {
            let mut su = 0.0;
            let mut mu = 0.0;
            for j in a.row_ptr[i]..a.row_ptr[i + 1] {
                let uj = prev[a.col_idx[j] as usize];
                su += a.values[j] * uj;
                mu += mass[j] * uj;
            }
            force[i] = force[i] - (1.0 - beta) * su + (1.0 / dt) * mu;
        }
        for (s, m) in a.values.iter_mut().zip(mass) {
            *s = beta * *s + (1.0 / dt) * m;
        }
    }
    Ok(())
}

/// Constant-step BDF weights `a(1..6)` of `BDF_CRS`, already divided by `dt`.
fn bdf_weights(dt: f64, order: usize) -> Result<[f64; 6], TimeError> {
    let mut a = [0.0; 6];
    match order {
        1 => {
            a[0] = 1.0;
            a[1] = -1.0;
        }
        2 => {
            a[0] = 1.5;
            a[1] = -2.0;
            a[2] = 0.5;
        }
        3 => {
            a[0] = 11.0 / 6.0;
            a[1] = -3.0;
            a[2] = 3.0 / 2.0;
            a[3] = -1.0 / 3.0;
        }
        4 => {
            a[0] = 25.0 / 12.0;
            a[1] = -4.0;
            a[2] = 3.0;
            a[3] = -4.0 / 3.0;
            a[4] = 1.0 / 4.0;
        }
        5 => {
            a[0] = 137.0 / 60.0;
            a[1] = -5.0;
            a[2] = 5.0;
            a[3] = -10.0 / 3.0;
            a[4] = 5.0 / 4.0;
            a[5] = -1.0 / 5.0;
        }
        o => return Err(TimeError::InvalidBdfOrder(o as i64)),
    }
    for w in a.iter_mut() {
        *w /= dt;
    }
    Ok(a)
}

fn check_prev(prev: &[&[f64]], order: usize, n: usize, caller: &'static str) -> Result<(), TimeError> {
    if prev.len() < order {
        return Err(TimeError::MissingHistory {
            caller,
            needed: order,
            got: prev.len(),
        });
    }
    if prev.iter().take(order).any(|p| p.len() < n) {
        return Err(TimeError::Dim(format!("{caller}: previous solution shorter than {n}")));
    }
    Ok(())
}

/// `BDF_CRS`: constant-step BDF of `order` 1..=5 (0 does nothing). `prev[k]` is
/// the solution `k + 1` steps back.
pub fn bdf_crs(
    dt: f64,
    a: &mut CsrMatrix,
    d: &CrsTimeData<'_>,
    force: &mut [f64],
    prev: &[&[f64]],
    order: usize,
) -> Result<(), TimeError> {
    if order == 0 {
        return Ok(());
    }
    diag_check(a, force.len(), "BDF_CRS")?;
    let w = bdf_weights(dt, order)?;
    check_prev(prev, order, a.nrows, "BDF_CRS")?;
    sum_history_crs(a, d, force, prev, order, &w)
}

/// `VBDF_CRS`: variable-step BDF. `dts[0]` is the current step, `dts[1..]` the
/// earlier ones, newest first. Orders above 3 warn and use three terms.
pub fn vbdf_crs(
    dts: &[f64],
    a: &mut CsrMatrix,
    d: &CrsTimeData<'_>,
    force: &mut [f64],
    prev: &[&[f64]],
    order: usize,
    logger: &dyn Logger,
) -> Result<(), TimeError> {
    diag_check(a, force.len(), "VBDF_CRS")?;
    if order > 3 {
        logger.warn("VBDF_CRS", "Variable timestep BDF implemented only to order 3");
    }
    let used = order.min(3);
    if dts.len() < used.max(1) {
        return Err(TimeError::MissingHistory {
            caller: "VBDF_CRS",
            needed: used.max(1),
            got: dts.len(),
        });
    }
    check_prev(prev, used, a.nrows, "VBDF_CRS")?;
    let c = vbdf_coeffs(dts, order);
    let mut w = [0.0; 6];
    w[..4].copy_from_slice(&c);
    sum_history_crs(a, d, force, prev, used, &w)
}

/// Shared body of `BDF_CRS` and `VBDF_CRS`: `Force -= M * sum(w[k] u_k)`,
/// `K += w[0] M`.
///
/// The two Fortran routines add `a(1) * M` at different points (inside the
/// row loop versus after it), which does not change any value.
fn sum_history_crs(
    a: &mut CsrMatrix,
    d: &CrsTimeData<'_>,
    force: &mut [f64],
    prev: &[&[f64]],
    terms: usize,
    w: &[f64; 6],
) -> Result<(), TimeError> {
    let n = a.nrows;
    if d.lumped {
        let ml = lumped_vals(d, a, "BDF_CRS")?;
        for i in 0..n {
            for k in 1..=terms {
                let ui = prev[k - 1][i];
                force[i] = force[i] - ml[i] * w[k] * ui;
            }
            let j = diag_index(a, i)?;
            a.values[j] += w[0] * ml[i];
        }
    } else {
        let mass = mass_vals(d, a, "BDF_CRS")?;
        for i in 0..n {
            let mut mu = 0.0;
            for j in a.row_ptr[i]..a.row_ptr[i + 1] {
                let col = a.col_idx[j] as usize;
                for k in 1..=terms {
                    let uj = prev[k - 1][col];
                    mu = mu - mass[j] * w[k] * uj;
                }
                a.values[j] += w[0] * mass[j];
            }
            force[i] += mu;
        }
    }
    Ok(())
}

/// `FractionalStep_CRS`. In the lumped branch the Fortran scales only the
/// diagonal entry by `MassCoeff`; off-diagonal stiffness is left as is.
pub fn fractional_step_crs(
    dt: f64,
    a: &mut CsrMatrix,
    d: &CrsTimeData<'_>,
    force: &mut [f64],
    prev: &[f64],
    fs: FsCoeffs,
) -> Result<(), TimeError> {
    diag_check(a, force.len(), "FractionalStep_CRS")?;
    let n = a.nrows;
    if d.mass_lumped.is_some() {
        let ml = lumped_vals(d, a, "FractionalStep_CRS")?;
        for i in 0..n {
            let su = row_dot(a, &a.values, i, prev);
            let j = diag_index(a, i)?;
            let ui = prev[a.col_idx[j] as usize];
            force[i] = force[i] - fs.force_coeff * su + (1.0 / dt) * ml[i] * ui;
            a.values[j] = fs.mass_coeff * a.values[j] + (1.0 / dt) * ml[i];
        }
    } else {
        let mass = mass_vals(d, a, "FractionalStep_CRS")?;
        for i in 0..n {
            let mut su = 0.0;
            let mut mu = 0.0;
            for j in a.row_ptr[i]..a.row_ptr[i + 1] {
                let uj = prev[a.col_idx[j] as usize];
                su += a.values[j] * uj;
                mu += mass[j] * uj;
            }
            force[i] = force[i] - fs.force_coeff * su + (1.0 / dt) * mu;
        }
        for (s, m) in a.values.iter_mut().zip(mass) {
            *s = fs.mass_coeff * *s + (1.0 / dt) * m;
        }
    }
    Ok(())
}

/// `Add1stOrderTime_CRS`. `prev[k]` is the solution `k + 1` steps back and
/// `prev.len()` plays the role of `Order = min(DoneTime, Solver % Order)`;
/// nothing happens when it is empty. `aux.cur` is the current iterate and
/// `aux.dts` the previous step sizes, newest first.
pub fn add_1st_order_time_crs(
    a: &mut CsrMatrix,
    d: &CrsTimeData<'_>,
    force: &mut [f64],
    dt: f64,
    scheme: &Scheme,
    prev: &[&[f64]],
    aux: &Aux1st<'_>,
) -> Result<(), TimeError> {
    if prev.is_empty() || matches!(scheme, Scheme::None) {
        return Ok(());
    }
    let n = a.nrows;
    check_prev(prev, 1, n, "Add1stOrderTime_CRS")?;
    match *scheme {
        Scheme::None => Ok(()),
        Scheme::FractionalStep(fs) => fractional_step_crs(dt, a, d, force, prev[0], fs),
        Scheme::Bdf(o) | Scheme::Vbdf(o) => {
            let order = (o as usize).min(prev.len());
            let constant = !matches!(scheme, Scheme::Vbdf(_)) && {
                order <= 1 || {
                    let tail = aux.dts.get(..order - 1).ok_or(TimeError::MissingHistory {
                        caller: "Add1stOrderTime_CRS",
                        needed: order - 1,
                        got: aux.dts.len(),
                    })?;
                    tail.iter().all(|x| (x - dt).abs() <= 1.0e-6 * dt)
                }
            };
            if constant {
                bdf_crs(dt, a, d, force, prev, order)
            } else {
                let mut dts = vec![dt];
                dts.extend_from_slice(aux.dts.get(..order.saturating_sub(1)).ok_or(
                    TimeError::MissingHistory {
                        caller: "Add1stOrderTime_CRS",
                        needed: order.saturating_sub(1),
                        got: aux.dts.len(),
                    },
                )?);
                vbdf_crs(&dts, a, d, force, prev, order, aux.logger)
            }
        }
        Scheme::RungeKutta => {
            let cur = aux.cur.ok_or(TimeError::MissingCurrent("RungeKutta_CRS"))?;
            runge_kutta_crs(dt, a, d, force, prev[0], cur)
        }
        Scheme::AdamsBashforth { .. } | Scheme::AdamsMoulton { .. } => {
            Err(TimeError::NotImplemented("Add1stOrderTime_CRS: Adams methods"))
        }
        Scheme::NewmarkBeta(beta) => newmark_beta_crs(dt, a, d, force, prev[0], beta),
        Scheme::Bossak { .. } | Scheme::GeneralizedAlpha { .. } => Err(TimeError::NotImplemented(
            "Add1stOrderTime_CRS for a second-order scheme",
        )),
    }
}
