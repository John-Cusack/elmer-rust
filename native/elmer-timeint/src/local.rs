//! Element-local time discretisation: the first-order routines of
//! TimeIntegrate.F90 (`NewmarkBeta`, `BDFLocal`, `VBDFLocal`, `RungeKutta`,
//! Adams-Bashforth/Moulton, `FractionalStep`) and `Add1stOrderTime`
//! (SolverBasics.F90:245).
//!
//! The arithmetic keeps the Fortran evaluation order (left-to-right sums over
//! `j`, the same grouping of every product) so results are bitwise equal to the
//! Fortran compiled without fused multiply-add.

use crate::error::TimeError;
use crate::scheme::{effective_pred_corr_order, FsCoeffs, Scheme};
use elmer_core::{LocalSystem, Logger, NullLogger};

/// State `AdamsBashforth` keeps in `Element % propertydata`: `2n` values, the
/// previous residual and the previous force. Keep one per element.
#[derive(Clone, Debug, Default)]
pub struct AdamsState {
    pub values: Option<Vec<f64>>,
}

/// Everything `Add1stOrderTime` reads from the solver besides the scheme and
/// the previous solutions.
pub struct Aux1st<'a> {
    /// `Solver % DoneTime`; only the Adams schemes read it.
    pub done_time: usize,
    /// Previous step sizes, newest first: `Dts(2..Order)` of the Fortran, i.e.
    /// `Timestep size` `PrevValues(1, 1..Order-1)`. Needed by BDF when
    /// `Order > 1`.
    pub dts: &'a [f64],
    /// Current iterate `Solver % Variable % Values` (Runge-Kutta).
    pub cur: Option<&'a [f64]>,
    /// Per-element state for the Adams schemes.
    pub adams: Option<&'a mut AdamsState>,
    /// `Solver % Matrix % Lumped`: lump the mass matrix first.
    pub lumped: bool,
    pub logger: &'a dyn Logger,
}

impl Default for Aux1st<'_> {
    fn default() -> Self {
        Aux1st {
            done_time: 0,
            dts: &[],
            cur: None,
            adams: None,
            lumped: false,
            logger: &NullLogger,
        }
    }
}

/// The row-sum lumping of `Add1stOrderTime` (SolverBasics.F90:276): off-diagonals
/// are dropped and the diagonal is scaled so the total mass is kept. `m` is
/// row-major `n * n`.
pub fn lump_matrix(m: &mut [f64], n: usize) {
    let mut s = 0.0;
    let mut t = 0.0;
    for i in 0..n {
        for j in 0..n {
            s += m[i * n + j];
            if i != j {
                m[i * n + j] = 0.0;
            }
        }
        t += m[i * n + i];
    }
    if t != 0.0 {
        for k in 0..n {
            m[k * n + k] = m[k * n + k] * s / t;
        }
    }
}

fn check_len(local: &LocalSystem, caller: &'static str) -> Result<(), TimeError> {
    let n = local.n;
    if local.stiff.len() != n * n || local.force.len() != n {
        return Err(TimeError::Dim(format!(
            "{caller}: stiff {} / force {} for n = {n}",
            local.stiff.len(),
            local.force.len()
        )));
    }
    Ok(())
}

fn mass_ref<'a>(
    mass: &'a Option<Vec<f64>>,
    n: usize,
    caller: &'static str,
) -> Result<&'a [f64], TimeError> {
    let m = mass.as_deref().ok_or(TimeError::MissingMass(caller))?;
    if m.len() != n * n {
        return Err(TimeError::Dim(format!("{caller}: mass length {}", m.len())));
    }
    Ok(m)
}

fn need_prev(
    prev: &[&[f64]],
    needed: usize,
    n: usize,
    caller: &'static str,
) -> Result<(), TimeError> {
    if prev.len() < needed {
        return Err(TimeError::MissingHistory {
            caller,
            needed,
            got: prev.len(),
        });
    }
    if prev.iter().take(needed).any(|p| p.len() < n) {
        return Err(TimeError::Dim(format!("{caller}: previous solution shorter than {n}")));
    }
    Ok(())
}

/// Adds the first-order time terms to `local` (`Add1stOrderTime`).
///
/// `prev[k]` is the solution `k + 1` steps back in the element's local dof
/// order, `Order = min(DoneTime, Solver % Order)` vectors in all. The caller is
/// responsible for the Fortran's `UpdateGlobalForce` of the *untouched* force
/// before this call (clone `local.force` first): it feeds the bulk residual
/// that Crank-Nicolson uses later.
pub fn add_1st_order_time(
    local: &mut LocalSystem,
    dt: f64,
    scheme: &Scheme,
    prev: &[&[f64]],
    aux: &mut Aux1st<'_>,
) -> Result<(), TimeError> {
    check_len(local, "Add1stOrderTime")?;
    let n = local.n;
    if matches!(scheme, Scheme::None) {
        return Ok(());
    }
    if aux.lumped {
        if let Some(m) = local.mass.as_mut() {
            lump_matrix(m, n);
        }
    }
    match *scheme {
        Scheme::None => Ok(()),
        Scheme::NewmarkBeta(beta) => {
            need_prev(prev, 1, n, "NewmarkBeta")?;
            let m = mass_ref(&local.mass, n, "NewmarkBeta")?;
            newmark_beta(n, dt, m, &mut local.stiff, &mut local.force, prev[0], beta);
            Ok(())
        }
        Scheme::RungeKutta => {
            need_prev(prev, 1, n, "RungeKutta")?;
            let cur = aux.cur.ok_or(TimeError::MissingCurrent("RungeKutta"))?;
            let m = mass_ref(&local.mass, n, "RungeKutta")?;
            runge_kutta(n, dt, m, &mut local.stiff, &mut local.force, prev[0], cur);
            Ok(())
        }
        Scheme::FractionalStep(fs) => {
            need_prev(prev, 1, n, "FractionalStep")?;
            let m = mass_ref(&local.mass, n, "FractionalStep")?;
            fractional_step(n, dt, m, &mut local.stiff, &mut local.force, prev[0], fs);
            Ok(())
        }
        Scheme::Bdf(order) | Scheme::Vbdf(order) => {
            let order = (order as usize).min(prev.len());
            if order == 0 {
                return Err(TimeError::InvalidBdfOrder(0));
            }
            need_prev(prev, order, n, "BDFLocal")?;
            let constant = !matches!(scheme, Scheme::Vbdf(_))
                && dts_constant(dt, aux.dts, order, "BDFLocal")?;
            let m = mass_ref(&local.mass, n, "BDFLocal")?;
            if constant {
                if order > 5 {
                    return Err(TimeError::InvalidBdfOrder(order as i64));
                }
                bdf_local(n, dt, m, &mut local.stiff, &mut local.force, prev, order);
                Ok(())
            } else {
                let mut dts = Vec::with_capacity(order);
                dts.push(dt);
                dts.extend_from_slice(aux.dts.get(..order.saturating_sub(1)).ok_or(
                    TimeError::MissingHistory {
                        caller: "VBDFLocal",
                        needed: order.saturating_sub(1),
                        got: aux.dts.len(),
                    },
                )?);
                vbdf_local(n, &dts, m, &mut local.stiff, &mut local.force, prev, order, aux.logger)
            }
        }
        Scheme::AdamsBashforth { zeta, order } => {
            need_prev(prev, 1, n, "AdamsBashforth")?;
            let pc = effective_pred_corr_order(order, aux.done_time);
            let state = aux
                .adams
                .as_deref_mut()
                .ok_or(TimeError::AdamsMoultonBeforeBashforth)?;
            let m = mass_ref(&local.mass, n, "AdamsBashforth")?;
            adams_bashforth(n, dt, m, &mut local.stiff, &mut local.force, prev[0], zeta, pc, state);
            Ok(())
        }
        Scheme::AdamsMoulton { order } => {
            need_prev(prev, 2, n, "AdamsMoulton")?;
            let pc = effective_pred_corr_order(order, aux.done_time);
            let state = aux
                .adams
                .as_deref()
                .ok_or(TimeError::AdamsMoultonBeforeBashforth)?;
            let m = mass_ref(&local.mass, n, "AdamsMoulton")?;
            adams_moulton(n, dt, m, &mut local.stiff, &mut local.force, prev, pc, state)
        }
        Scheme::Bossak { .. } | Scheme::GeneralizedAlpha { .. } => {
            Err(TimeError::NotImplemented("Add1stOrderTime for a second-order scheme"))
        }
    }
}

/// `ConstantDt` test of `Add1stOrderTime`: all previous sizes within 1e-6
/// relative of `dt`.
fn dts_constant(dt: f64, dts: &[f64], order: usize, caller: &'static str) -> Result<bool, TimeError> {
    if order <= 1 {
        return Ok(true);
    }
    let prev = dts.get(..order - 1).ok_or(TimeError::MissingHistory {
        caller,
        needed: order - 1,
        got: dts.len(),
    })?;
    Ok(prev.iter().all(|d| (d - dt).abs() <= 1.0e-6 * dt))
}

/// `NewmarkBeta` (TimeIntegrate.F90:87): `Beta = 1` is implicit Euler.
pub fn newmark_beta(
    n: usize,
    dt: f64,
    mass: &[f64],
    stiff: &mut [f64],
    force: &mut [f64],
    prev: &[f64],
    beta: f64,
) {
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            s = s + (1.0 / dt) * mass[i * n + j] * prev[j]
                - (1.0 - beta) * stiff[i * n + j] * prev[j];
        }
        for j in 0..n {
            stiff[i * n + j] = beta * stiff[i * n + j] + (1.0 / dt) * mass[i * n + j];
        }
        force[i] += s;
    }
}

/// `RungeKutta` stage (TimeIntegrate.F90:56).
pub fn runge_kutta(
    n: usize,
    dt: f64,
    mass: &[f64],
    stiff: &mut [f64],
    force: &mut [f64],
    prev: &[f64],
    cur: &[f64],
) {
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            s = s + (1.0 / dt) * mass[i * n + j] * prev[j] - stiff[i * n + j] * cur[j];
        }
        for j in 0..n {
            stiff[i * n + j] = (1.0 / dt) * mass[i * n + j];
        }
        force[i] += s;
    }
}

/// `AdamsBashforth` (TimeIntegrate.F90:~140). `pc_order` is the already
/// capped `PredCorrOrder`. The state is allocated on first use; the Fortran
/// leaves it uninitialised, so the zeros here only matter when `pc_order == 0`.
#[allow(clippy::too_many_arguments)]
pub fn adams_bashforth(
    n: usize,
    dt: f64,
    mass: &[f64],
    stiff: &mut [f64],
    force: &mut [f64],
    prev: &[f64],
    zeta: f64,
    pc_order: usize,
    state: &mut AdamsState,
) {
    let pd = state.values.get_or_insert_with(|| vec![0.0; 2 * n]);
    for i in 0..n {
        let mut s_curr = 0.0;
        let mut m_curr = 0.0;
        for j in 0..n {
            s_curr += stiff[i * n + j] * prev[j];
            m_curr += (1.0 / dt) * mass[i * n + j] * prev[j];
        }
        let curr_res = -s_curr;
        let (residual, pre_force) = if pc_order == 1 {
            (curr_res, force[i])
        } else {
            (pd[i], pd[i + n])
        };
        pd[i] = curr_res;
        pd[i + n] = force[i];

        force[i] = force[i] - s_curr
            + m_curr
            + 0.5 * zeta * (force[i] - pre_force)
            + 0.5 * zeta * (curr_res - residual);

        for j in 0..n {
            stiff[i * n + j] = (1.0 / dt) * mass[i * n + j];
        }
    }
}

/// `AdamsMoulton` (TimeIntegrate.F90:~205). `prev[0]` is the corrector,
/// `prev[1]` the converged previous solution.
#[allow(clippy::too_many_arguments)]
pub fn adams_moulton(
    n: usize,
    dt: f64,
    mass: &[f64],
    stiff: &mut [f64],
    force: &mut [f64],
    prev: &[&[f64]],
    pc_order: usize,
    state: &AdamsState,
) -> Result<(), TimeError> {
    let pd = state
        .values
        .as_deref()
        .ok_or(TimeError::AdamsMoultonBeforeBashforth)?;
    for i in 0..n {
        let mut s_curr = 0.0;
        let mut m_curr = 0.0;
        for j in 0..n {
            s_curr += stiff[i * n + j] * prev[0][j];
            m_curr += (1.0 / dt) * mass[i * n + j] * prev[1][j];
        }
        for j in 0..n {
            stiff[i * n + j] = (1.0 / dt) * mass[i * n + j];
        }
        let residual = pd[i];
        let pre_force = pd[i + n];
        if pc_order == 1 {
            force[i] = force[i] + m_curr - s_curr;
        } else {
            force[i] = 0.5 * (force[i] + pre_force) + m_curr + 0.5 * (-s_curr + residual);
        }
    }
    Ok(())
}

/// `BDFLocal` (TimeIntegrate.F90:~260), constant step, order 1..=5.
pub fn bdf_local(
    n: usize,
    dt: f64,
    mass: &[f64],
    stiff: &mut [f64],
    force: &mut [f64],
    prev: &[&[f64]],
    order: usize,
) {
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            let m = mass[i * n + j];
            s += match order {
                1 => (1.0 / dt) * m * prev[0][j],
                2 => (1.0 / dt) * m * (2.0 * prev[0][j] - 0.5 * prev[1][j]),
                3 => {
                    (1.0 / dt)
                        * m
                        * (3.0 * prev[0][j] - (3.0 / 2.0) * prev[1][j] + (1.0 / 3.0) * prev[2][j])
                }
                4 => {
                    (1.0 / dt)
                        * m
                        * (4.0 * prev[0][j] - 3.0 * prev[1][j] + (4.0 / 3.0) * prev[2][j]
                            - (1.0 / 4.0) * prev[3][j])
                }
                _ => {
                    (1.0 / dt)
                        * m
                        * (5.0 * prev[0][j] - 5.0 * prev[1][j] + (10.0 / 3.0) * prev[2][j]
                            - (5.0 / 4.0) * prev[3][j]
                            + (1.0 / 5.0) * prev[4][j])
                }
            };
        }
        force[i] += s;
        for j in 0..n {
            let m = mass[i * n + j];
            let k = stiff[i * n + j];
            stiff[i * n + j] = match order {
                1 => (1.0 / dt) * m + k,
                2 => (1.5 / dt) * m + k,
                3 => (11.0 / (6.0 * dt)) * m + k,
                4 => (25.0 / (12.0 * dt)) * m + k,
                _ => (137.0 / (60.0 * dt)) * m + k,
            };
        }
    }
}

/// Variable-step BDF coefficients `a(1..4)` shared by `VBDFLocal` and `VBDF_CRS`.
/// `dts[0]` is the current step, `dts[1..]` the earlier ones, newest first.
pub fn vbdf_coeffs(dts: &[f64], order: usize) -> [f64; 4] {
    let mut a = [0.0; 4];
    a[0] = 1.0 / dts[0];
    a[1] = -1.0 / dts[0];
    if order >= 2 {
        a[0] += 1.0 / (dts[0] + dts[1]);
        a[1] -= (1.0 + dts[0] / dts[1]) / (dts[0] + dts[1]);
        a[2] = (dts[0] / dts[1]) / (dts[0] + dts[1]);
    }
    if order >= 3 {
        a[0] += 1.0 / (dts[0] + dts[1] + dts[2]);
        a[1] -= (1.0 + dts[0] / dts[1] * (1.0 + (dts[0] + dts[1]) / (dts[1] + dts[2])))
            / (dts[0] + dts[1] + dts[2]);
        a[2] += (dts[0] / dts[1] * (1.0 + (dts[0] + dts[1]) / (dts[1] + dts[2]))
            + dts[0] / dts[2] * (dts[0] + dts[1]) / (dts[1] + dts[2]))
            / (dts[0] + dts[1] + dts[2]);
        a[3] = -(dts[0] / dts[2]) * (dts[0] + dts[1]) / (dts[1] + dts[2])
            / (dts[0] + dts[1] + dts[2]);
    }
    a
}

/// `VBDFLocal` (TimeIntegrate.F90:~372). Orders above 3 warn and use order 3
/// terms in the sum, as the Fortran does.
#[allow(clippy::too_many_arguments)]
pub fn vbdf_local(
    n: usize,
    dts: &[f64],
    mass: &[f64],
    stiff: &mut [f64],
    force: &mut [f64],
    prev: &[&[f64]],
    order: usize,
    logger: &dyn Logger,
) -> Result<(), TimeError> {
    if order > 3 {
        logger.warn("VBDFLocal", "Variable timestep BDF implemented only to order 3");
    }
    let used = order.min(3);
    if dts.len() < used {
        return Err(TimeError::MissingHistory {
            caller: "VBDFLocal",
            needed: used,
            got: dts.len(),
        });
    }
    let a = vbdf_coeffs(dts, order);
    for i in 0..n {
        let mut s = 0.0;
        for k in 0..used {
            for j in 0..n {
                s = s - a[k + 1] * mass[i * n + j] * prev[k][j];
            }
        }
        force[i] += s;
        for j in 0..n {
            stiff[i * n + j] = a[0] * mass[i * n + j] + stiff[i * n + j];
        }
    }
    Ok(())
}

/// `FractionalStep` (TimeIntegrate.F90:~520).
pub fn fractional_step(
    n: usize,
    dt: f64,
    mass: &[f64],
    stiff: &mut [f64],
    force: &mut [f64],
    prev: &[f64],
    fs: FsCoeffs,
) {
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            s = s + (1.0 / dt) * mass[i * n + j] * prev[j]
                - fs.force_coeff * stiff[i * n + j] * prev[j];
        }
        force[i] += s;
        for j in 0..n {
            stiff[i * n + j] = fs.mass_coeff * stiff[i * n + j] + (1.0 / dt) * mass[i * n + j];
        }
    }
}
