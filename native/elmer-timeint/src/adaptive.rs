//! Adaptive sub-stepping and predictor-corrector time-step control as pure
//! state machines. The caller solves; these decide the step sizes.
//!
//! Sources: ExecSimulation (ElmerSolver.F90:3187-3410) for the halving /
//! divergence control, MainUtils.F90:5958-6200 for the PI controller.

use crate::error::TimeError;
use elmer_core::Params;

/// Why a sub-step controller exists: `Adaptive Timestepping` or `Convergence Control`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SubstepMode {
    /// Step-doubling error control (`Adaptive Timestepping = True`).
    Error,
    /// Shrink on nonlinear divergence (`Convergence Control = True`).
    Divergence,
}

/// Outcome of one trial sub-step.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Decision {
    /// The sub-step of size `used` is kept; `next` is the proposed `ddt`.
    Accept { used: f64, next: f64 },
    /// Restore the saved solution and redo with `next`.
    Reject { next: f64 },
}

/// A sub-step to solve: size and the physical time at its end.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Substep {
    pub dt: f64,
    pub time: f64,
}

/// Settings read from the `Simulation` section.
#[derive(Clone, Debug, PartialEq)]
pub struct SubstepConfig {
    pub mode: SubstepMode,
    /// `Adaptive Time Error` (required for [`SubstepMode::Error`]).
    pub limit: f64,
    pub keep_smallest: i32,
    pub max_timestep: Option<f64>,
    pub min_timestep: Option<f64>,
    pub increase: f64,
    pub decrease: f64,
    pub rough: bool,
    pub smart: bool,
}

impl SubstepConfig {
    /// `None` when neither `Adaptive Timestepping` nor `Convergence Control`
    /// is set. Adaptive stepping wins when both are.
    pub fn from_params(sim: &dyn Params) -> Result<Option<Self>, TimeError> {
        let adaptive = sim.logical("adaptive timestepping").unwrap_or(false);
        let divergence = sim.logical("convergence control").unwrap_or(false);
        if !adaptive && !divergence {
            return Ok(None);
        }
        let mut limit = 0.0;
        let mut keep_smallest = 0;
        if adaptive {
            limit = sim.real("adaptive time error").ok_or_else(|| {
                TimeError::Config(
                    "Adaptive Time Error must be given for adaptive stepping scheme.".into(),
                )
            })?;
            keep_smallest = sim.integer("adaptive keep smallest").unwrap_or(0).max(0);
        }
        Ok(Some(SubstepConfig {
            mode: if adaptive {
                SubstepMode::Error
            } else {
                SubstepMode::Divergence
            },
            limit,
            keep_smallest,
            max_timestep: sim.real("adaptive max timestep"),
            min_timestep: sim.real("adaptive min timestep"),
            increase: sim.real("adaptive increase coefficient").unwrap_or(2.0),
            decrease: sim.real("adaptive decrease coefficient").unwrap_or(0.5),
            rough: sim.logical("adaptive rough timestep").unwrap_or(false),
            smart: sim.logical("adaptive smart timestep").unwrap_or(false),
        }))
    }
}

/// Splits one nominal step of size `dt` into sub-steps.
///
/// Per nominal step: [`begin_step`](Self::begin_step), then repeat
/// [`next_substep`](Self::next_substep) / solve /
/// [`finish_error`](Self::finish_error) (or
/// [`finish_divergence`](Self::finish_divergence)) until `next_substep`
/// returns `None`.
#[derive(Clone, Debug)]
pub struct SubstepController {
    cfg: SubstepConfig,
    ddt: f64,
    cum: f64,
    dt: f64,
    start: f64,
    max_dt: f64,
    min_dt: f64,
    smallest_count: i32,
    /// `StepControl` of ExecSimulation, saved between nominal steps.
    step_control: i32,
    pending: Option<f64>,
}

impl SubstepController {
    pub fn new(cfg: SubstepConfig) -> Self {
        Self {
            cfg,
            ddt: -1.0,
            cum: 0.0,
            dt: 0.0,
            start: 0.0,
            max_dt: 0.0,
            min_dt: 0.0,
            smallest_count: 0,
            step_control: -1,
            pending: None,
        }
    }

    /// Current proposed sub-step size (`ddt`).
    pub fn ddt(&self) -> f64 {
        self.ddt
    }

    /// Time covered so far within the nominal step.
    pub fn cumulative(&self) -> f64 {
        self.cum
    }

    /// Starts a nominal step of size `dt` that begins at `start_time`.
    pub fn begin_step(&mut self, dt: f64, start_time: f64) {
        self.dt = dt;
        self.start = start_time;
        self.cum = 0.0;
        self.pending = None;
        self.max_dt = match self.cfg.max_timestep {
            Some(m) => m.min(dt),
            None => dt,
        };
        self.min_dt = self.cfg.min_timestep.unwrap_or(1.0e-8 * self.max_dt);
        if self.cfg.mode == SubstepMode::Divergence {
            self.ddt *= self.cfg.increase;
        } else {
            self.smallest_count = 0;
        }
        if self.ddt < 0.0 || self.ddt > self.max_dt {
            self.ddt = self.max_dt;
        }
    }

    /// The next trial, or `None` when the nominal step is complete.
    pub fn next_substep(&mut self) -> Option<Substep> {
        if self.cum >= self.dt - 1.0e-12 {
            return None;
        }
        if !self.cfg.rough {
            self.ddt = self.ddt.min(self.dt - self.cum);
            if self.cfg.smart && self.dt - self.cum - self.ddt > 1.0e-12 {
                self.ddt = self.ddt.min((self.dt - self.cum) / 2.0);
            }
        }
        self.pending = Some(self.ddt);
        Some(Substep {
            dt: self.ddt,
            time: self.start + self.cum + self.ddt,
        })
    }

    /// Result of the step-doubling comparison. `max_err` is the externally
    /// given `Adaptive Error Measure` difference or the largest relative norm
    /// change.
    pub fn finish_error(&mut self, max_err: f64) -> Decision {
        let used = self.ddt;
        if max_err < self.cfg.limit || used <= self.min_dt {
            self.cum += used;
            if self.smallest_count >= self.cfg.keep_smallest || self.step_control > 0 {
                self.ddt = (self.cfg.increase * self.ddt).min(self.max_dt);
                self.step_control = 1;
                self.smallest_count = 0;
            } else {
                self.step_control = 0;
                self.smallest_count += 1;
            }
            Decision::Accept {
                used,
                next: self.ddt,
            }
        } else {
            self.ddt *= self.cfg.decrease;
            self.step_control = -1;
            Decision::Reject { next: self.ddt }
        }
    }

    /// Result of a solve under `Convergence Control`.
    pub fn finish_divergence(&mut self, diverged: bool) -> Result<Decision, TimeError> {
        let used = self.ddt;
        if !diverged {
            self.cum += used;
            self.ddt *= self.cfg.increase;
            if self.ddt > self.max_dt {
                self.ddt = self.max_dt;
                self.step_control = 0;
            }
            Ok(Decision::Accept {
                used,
                next: self.ddt,
            })
        } else {
            if self.ddt < self.min_dt * (1.0 + 1.0e-8) {
                return Err(TimeError::NoStableTimestep);
            }
            self.ddt = (self.cfg.decrease * self.ddt).max(self.min_dt);
            self.step_control = 1;
            Ok(Decision::Reject { next: self.ddt })
        }
    }
}

/// `Predictor-Corrector Control` settings (`ReadPredCorrParams`, MainUtils.F90:6133).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PredCorrConfig {
    pub order: i32,
    pub epsilon: f64,
    pub beta1: f64,
    pub beta2: f64,
}

impl PredCorrConfig {
    /// Looks in `solver` first, then `sim`, then the defaults.
    pub fn from_params(solver: &dyn Params, sim: &dyn Params) -> Self {
        let order = solver
            .integer("predictor-corrector scheme order")
            .or_else(|| sim.integer("predictor-corrector scheme order"))
            .unwrap_or(2);
        let r = |k: &str| solver.real(k).or_else(|| sim.real(k));
        PredCorrConfig {
            order,
            epsilon: r("predictor-corrector control tolerance").unwrap_or(1.0e-6),
            beta1: r("predictor-corrector control beta 1").unwrap_or(0.6 / (order as f64 + 1.0)),
            beta2: r("predictor-corrector control beta 2").unwrap_or(-0.2 / (order as f64 + 1.0)),
        }
    }
}

/// `PredCorrErrorEstimate`: local truncation error of the AB/AM pair.
pub fn pred_corr_error_estimate(dt: f64, order: i32, time_error: f64, zeta: f64) -> f64 {
    if dt > 0.0 {
        if order == 2 {
            time_error * zeta / dt / (zeta + 1.0) / 3.0
        } else {
            time_error / dt / 2.0
        }
    } else {
        0.0
    }
}

/// `TimeStepLimiter`: `dt = dtOld * (1 + k atan((g - 1)/k))`, `k = 2` by default.
pub fn time_step_limiter(dt_old: f64, gfactor: f64) -> f64 {
    dt_old * (1.0 + 2.0 * ((gfactor - 1.0) * 0.5).atan())
}

/// `TimeStepController`: PI controller. Keeps `dt_old` when an error is zero.
pub fn time_step_controller(
    dt_old: f64,
    eta: f64,
    eta_old: f64,
    epsilon: f64,
    beta1: f64,
    beta2: f64,
) -> f64 {
    if eta != 0.0 && eta_old != 0.0 {
        let g = (epsilon / eta).powf(beta1) * (epsilon / eta_old).powf(beta2);
        time_step_limiter(dt_old, g)
    } else {
        dt_old
    }
}

/// `PredictorCorrectorControl` (MainUtils.F90:5958) with its `SAVE`d locals.
#[derive(Clone, Debug)]
pub struct PredCorrControl {
    cfg: PredCorrConfig,
    dt_old: f64,
    eta_old: f64,
    zeta: f64,
}

impl PredCorrControl {
    pub fn new(cfg: PredCorrConfig) -> Self {
        Self {
            cfg,
            dt_old: 0.0,
            eta_old: 0.0,
            zeta: 1.0,
        }
    }

    /// `Adams Zeta` for the predictor of the coming step.
    pub fn zeta(&self) -> f64 {
        self.zeta
    }

    /// The step size to use. `timestep` is the 1-based step within the
    /// interval (the Fortran passes it as `RealTimestep`); `dt` the nominal
    /// size; `time_error` the max-norm of `Values - PrevValues(:,1)` over all
    /// partitions, read after the previous step (unused for `timestep <= 2`).
    /// Step 1 returns 0, the Fortran's "do nothing" first predictor step.
    pub fn control(&mut self, timestep: usize, dt: f64, time_error: f64) -> f64 {
        match timestep {
            1 => {
                self.dt_old = dt;
                self.zeta = 1.0;
                0.0
            }
            2 => {
                self.zeta = 1.0;
                self.dt_old
            }
            _ => {
                let order = if timestep == 3 { 1 } else { self.cfg.order };
                let eta = pred_corr_error_estimate(self.dt_old, order, time_error, self.zeta);
                if timestep == 3 {
                    self.eta_old = eta;
                }
                let new_dt = time_step_controller(
                    self.dt_old,
                    eta,
                    self.eta_old,
                    self.cfg.epsilon,
                    self.cfg.beta1,
                    self.cfg.beta2,
                );
                self.zeta = new_dt / self.dt_old;
                self.eta_old = eta;
                self.dt_old = new_dt;
                new_dt
            }
        }
    }
}
