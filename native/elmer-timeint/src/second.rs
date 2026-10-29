//! Second-order-in-time schemes: `Time2ndOrder`, `Time2ndOrder_CRS`,
//! `Update2ndOrder` (TimeIntegrate.F90:431-565) and `Add2ndOrderTime`
//! (SolverBasics.F90:494).
//!
//! `param_b >= 0` selects generalized-alpha with `rho_inf = param_b`;
//! otherwise Bossak with `Alpha = param_a`, exactly as the Fortran branches.

use crate::crs::{CrsTimeData, diag_check};
use crate::error::TimeError;
use crate::local::lump_matrix;
use elmer_core::{CsrMatrix, LocalSystem};

/// Column indices (0-based) of the `Solver % Variable % PrevValues` array that
/// the second-order schemes use.
pub mod col {
    /// Predicted velocity written by [`super::update_2nd_order`] (`PrevValues(:,1)`).
    pub const VEL_PRED: usize = 0;
    /// Predicted acceleration (`PrevValues(:,2)`).
    pub const ACC_PRED: usize = 1;
    /// Displacement at the previous step (`PrevValues(:,3)`).
    pub const X: usize = 2;
    /// Velocity (`PrevValues(:,4)`).
    pub const V: usize = 3;
    /// Acceleration (`PrevValues(:,5)`).
    pub const A: usize = 4;
    /// Generalized-alpha work column (`PrevValues(:,6)`).
    pub const A_ALPHA: usize = 5;
    /// Generalized-alpha acceleration used as `A2` (`PrevValues(:,7)`).
    pub const A2: usize = 6;
}

/// The seven `PrevValues` columns of a second-order solver, each of length `n`.
#[derive(Clone, Debug)]
pub struct SecondOrderState {
    pub cols: [Vec<f64>; 7],
}

impl SecondOrderState {
    pub fn zeros(n: usize) -> Self {
        Self {
            cols: std::array::from_fn(|_| vec![0.0; n]),
        }
    }

    pub fn len(&self) -> usize {
        self.cols[0].len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The column `Time2ndOrder` calls `A2`: column 7 for generalized-alpha,
    /// column 5 (the acceleration) for Bossak.
    pub fn a2(&self, param_b: f64) -> &[f64] {
        if param_b >= 0.0 {
            &self.cols[col::A2]
        } else {
            &self.cols[col::A]
        }
    }
}

/// `Alpha, Beta, Gamma, Delta` of the scheme.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Newmark2 {
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
    pub delta: f64,
}

impl Newmark2 {
    pub fn new(param_a: f64, param_b: f64) -> Self {
        if param_b >= 0.0 {
            Newmark2 {
                alpha: (2.0 * param_b - 1.0) / (param_b + 1.0),
                beta: 1.0 / ((param_b + 1.0) * (param_b + 1.0)),
                gamma: (3.0 - param_b) / (2.0 * (param_b + 1.0)),
                delta: param_b / (param_b + 1.0),
            }
        } else {
            let alpha = param_a;
            Newmark2 {
                alpha,
                beta: (1.0 - alpha) * (1.0 - alpha) / 4.0,
                gamma: 0.5 - alpha,
                delta: 0.0,
            }
        }
    }
}

/// `Time2ndOrder` (TimeIntegrate.F90:431) on row-major `n * n` matrices.
/// `a2` is the already chosen `A2b` column.
#[allow(clippy::too_many_arguments)]
pub fn time_2nd_order(
    n: usize,
    dt: f64,
    mass: &[f64],
    damp: &[f64],
    stiff: &mut [f64],
    force: &mut [f64],
    x: &[f64],
    v: &[f64],
    a: &[f64],
    a2b: &[f64],
    param_a: f64,
    param_b: f64,
) {
    let Newmark2 {
        alpha,
        beta,
        gamma,
        delta,
    } = Newmark2::new(param_a, param_b);
    let dt2 = dt * dt;
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            let m = mass[i * n + j];
            let d = damp[i * n + j];
            s = s + (1.0 - alpha) / (1.0 - delta) / (beta * dt2) * m * x[j];
            s = s + (1.0 - alpha) / (1.0 - delta) / (beta * dt) * m * v[j];
            s = s + delta / (1.0 - delta) * m * a[j];
            s = s
                - ((1.0 - alpha) * (1.0 - 1.0 / (2.0 * beta)) + alpha) / (1.0 - delta)
                    * m
                    * a2b[j];

            s = s + gamma / (beta * dt) * d * x[j];
            s = s + (gamma / beta - 1.0) * d * v[j];
            s = s - ((1.0 - gamma) + gamma * (1.0 - 1.0 / (2.0 * beta))) * dt * d * a2b[j];

            stiff[i * n + j] = stiff[i * n + j]
                + (1.0 - alpha) / (1.0 - delta) / (beta * dt2) * m
                + (gamma / (beta * dt)) * d;
        }
        force[i] += s;
    }
}

/// Adds the second-order time terms to `local` (`Add2ndOrderTime`).
///
/// `x, v, a, a2` are the element-local `PrevValues` columns 3, 4, 5 and 7.
/// `lumped` lumps mass and damping first. A missing damping matrix counts as
/// zero. As in the Fortran, `UpdateGlobalForce` of the untouched force is the
/// caller's job.
#[allow(clippy::too_many_arguments)]
pub fn add_2nd_order_time(
    local: &mut LocalSystem,
    dt: f64,
    x: &[f64],
    v: &[f64],
    a: &[f64],
    a2: &[f64],
    param_a: f64,
    param_b: f64,
    lumped: bool,
) -> Result<(), TimeError> {
    let n = local.n;
    if local.stiff.len() != n * n || local.force.len() != n {
        return Err(TimeError::Dim("Add2ndOrderTime: local system".into()));
    }
    if [x, v, a, a2].iter().any(|c| c.len() < n) {
        return Err(TimeError::Dim("Add2ndOrderTime: history shorter than n".into()));
    }
    let mut mass = local
        .mass
        .clone()
        .ok_or(TimeError::MissingMass("Time2ndOrder"))?;
    let mut damp = local.damp.clone().unwrap_or_else(|| vec![0.0; n * n]);
    if mass.len() != n * n || damp.len() != n * n {
        return Err(TimeError::Dim("Add2ndOrderTime: mass/damp length".into()));
    }
    if lumped {
        lump_matrix(&mut mass, n);
        lump_matrix(&mut damp, n);
        local.mass = Some(mass.clone());
        if local.damp.is_some() {
            local.damp = Some(damp.clone());
        }
    }
    let a2b = if param_b >= 0.0 { a2 } else { a };
    time_2nd_order(
        n,
        dt,
        &mass,
        &damp,
        &mut local.stiff,
        &mut local.force,
        x,
        v,
        a,
        a2b,
        param_a,
        param_b,
    );
    Ok(())
}

/// `Time2ndOrder_CRS` (TimeIntegrate.F90:487) on the global CRS system.
/// `state` holds the seven `PrevValues` columns.
pub fn time_2nd_order_crs(
    dt: f64,
    matrix: &mut CsrMatrix,
    data: &CrsTimeData<'_>,
    force: &mut [f64],
    param_a: f64,
    param_b: f64,
    state: &SecondOrderState,
) -> Result<(), TimeError> {
    let n = matrix.nrows;
    diag_check(matrix, force.len(), "Time2ndOrder_CRS")?;
    if data.lumped {
        // Add2ndOrderTime_CRS
        return Err(TimeError::NotImplemented(
            "Add2ndOrderTime_CRS: matrix lumping",
        ));
    }
    if state.len() < n {
        return Err(TimeError::Dim("Time2ndOrder_CRS: state shorter than n".into()));
    }
    let Newmark2 {
        alpha,
        beta,
        gamma,
        delta,
    } = Newmark2::new(param_a, param_b);
    let x = &state.cols[col::X];
    let v = &state.cols[col::V];
    let a = &state.cols[col::A];
    let a2 = state.a2(param_b);

    let dt2 = dt * dt;
    let c1 = (1.0 - alpha) / (1.0 - delta) / (beta * dt2);
    let c2 = (1.0 - alpha) / (1.0 - delta) / (beta * dt);
    let c3 = delta / (1.0 - delta);
    let c4 = -((1.0 - alpha) * (1.0 - 1.0 / (2.0 * beta)) + alpha) / (1.0 - delta);

    let c5 = gamma / (beta * dt);
    let c6 = gamma / beta - 1.0;
    let c7 = -((1.0 - gamma) + gamma * (1.0 - 1.0 / (2.0 * beta))) * dt;

    let rows = &matrix.row_ptr;
    let cols = &matrix.col_idx;
    if let Some(mass) = data.mass {
        for i in 0..n {
            let mut s = 0.0;
            for j in rows[i]..rows[i + 1] {
                let k = cols[j] as usize;
                s += mass[j] * (c1 * x[k] + c2 * v[k] + c3 * a[k] + c4 * a2[k]);
            }
            force[i] += s;
        }
        for (st, m) in matrix.values.iter_mut().zip(mass) {
            *st += c1 * m;
        }
    }
    if let Some(damp) = data.damp {
        for i in 0..n {
            let mut s = 0.0;
            for j in rows[i]..rows[i + 1] {
                let k = cols[j] as usize;
                s += damp[j] * (c5 * x[k] + c6 * v[k] + c7 * a2[k]);
            }
            force[i] += s;
        }
        for (st, d) in matrix.values.iter_mut().zip(damp) {
            *st += c5 * d;
        }
    }
    Ok(())
}

/// `Update2ndOrder` (TimeIntegrate.F90:~545): after the solve, refreshes the
/// predicted velocity and acceleration columns from the new solution `x`.
pub fn update_2nd_order(
    dt: f64,
    x: &[f64],
    prev: &mut SecondOrderState,
    param_a: f64,
    param_b: f64,
) {
    let Newmark2 {
        alpha,
        beta,
        gamma,
        delta,
    } = Newmark2::new(param_a, param_b);
    let n = x.len().min(prev.len());
    let dt2 = dt * dt;
    for i in 0..n {
        let dx = x[i] - prev.cols[col::X][i];
        let v = prev.cols[col::V][i];
        let a = prev.cols[col::A][i];
        if delta == 0.0 {
            prev.cols[col::VEL_PRED][i] = gamma / (beta * dt) * dx
                - (gamma - beta) / beta * v
                - dt * (gamma - 2.0 * beta) / (2.0 * beta) * a;
            prev.cols[col::ACC_PRED][i] =
                1.0 / (beta * dt2) * dx - 1.0 / (beta * dt) * v - (1.0 - 2.0 * beta) / (2.0 * beta) * a;
        } else {
            let a2 = prev.cols[col::A2][i];
            prev.cols[col::VEL_PRED][i] = gamma / (beta * dt) * dx
                - (gamma - beta) / beta * v
                - dt * (gamma - 2.0 * beta) / (2.0 * beta) * a2;
            prev.cols[col::ACC_PRED][i] = (1.0 - alpha) / (1.0 - delta) / (beta * dt2) * dx
                - (1.0 - alpha) / (1.0 - delta) / (beta * dt) * v
                - delta / (1.0 - delta) * a
                - (1.0 - alpha - 2.0 * beta) / (1.0 - delta) / (2.0 * beta) * a2;
            prev.cols[col::A_ALPHA][i] = (1.0 - delta) / (1.0 - alpha) * prev.cols[col::ACC_PRED][i]
                + delta / (1.0 - alpha) * a
                - alpha / (1.0 - alpha) * a2;
        }
    }
}
