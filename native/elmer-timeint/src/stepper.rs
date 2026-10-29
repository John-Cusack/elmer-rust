//! The `Timestep Intervals` / `Timestep Sizes` / `Output Intervals` loop of
//! `ExecSimulation` (ElmerSolver.F90:2703) and `InitializeIntervals`
//! (ElmerSolver.F90:1016) as an iterator over steps.
//!
//! Not modelled: time-parallel / sliced runs (`Parallel Timestepping`,
//! `Parallel Slices`), `Real Time Max`, and the `Periodic Time*` keyword
//! prefix (only `Time Period` and `Periodic Timesteps` are read).

use crate::adaptive::PredCorrControl;
use crate::error::TimeError;
use elmer_core::{Params, Value};

/// `Simulation Type`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum RunKind {
    Steady,
    Transient,
    Scanning,
}

impl RunKind {
    /// Maps the `Simulation Type` string; anything but `transient` and
    /// `scanning` is steady, as in ElmerSolver.F90:522.
    pub fn from_params(sim: &dyn Params) -> RunKind {
        match sim
            .string("simulation type")
            .map(|s| s.trim().to_ascii_lowercase())
            .as_deref()
        {
            Some("transient") => RunKind::Transient,
            Some("scanning") => RunKind::Scanning,
            _ => RunKind::Steady,
        }
    }

    fn is_time(self) -> bool {
        self != RunKind::Steady
    }
}

/// One step handed to the solvers.
#[derive(Clone, Debug, PartialEq)]
pub struct Step {
    /// 1-based over the whole run (`cum_Timestep`, `Timestep` variable).
    pub cum_step: usize,
    /// 1-based interval.
    pub interval: usize,
    /// 1-based step within the interval (`timestep`).
    pub timestep: usize,
    /// Step size (`sSize`).
    pub dt: f64,
    /// Physical time at the end of the step (`sTime`).
    pub time: f64,
    /// Total number of steps in the run (`stepcount`).
    pub total_steps: usize,
    /// `sPrevSizes`: earlier step sizes, newest first.
    pub prev_sizes: [f64; 5],
    /// `OutputIntervals(interval)` asks for a save after this step
    /// (`MOD(timestep-1, OutputIntervals) == 0`). A steady-state hit also saves.
    pub output_due: bool,
    /// `(cycle, time in cycle)` when the run is periodic.
    pub periodic: Option<(f64, f64)>,
    kind: RunKind,
}

impl Step {
    /// The banner line, e.g. `Time: 3/10:   6.000E+00` (as printed after the `MAIN: ` prefix); steady runs get the
    /// iteration line.
    pub fn time_line(&self) -> String {
        if self.kind.is_time() {
            format!(
                "Time: {}/{}:{}",
                self.cum_step,
                self.total_steps,
                fortran_es(self.time, 12, 3)
            )
        } else {
            format!("Steady state iteration: {:>12}", self.cum_step)
        }
    }
}

/// `ES<width>.<digits>` edit descriptor.
pub fn fortran_es(x: f64, width: usize, digits: usize) -> String {
    let body = if x == 0.0 {
        format!("{:.*}E+00", digits, 0.0)
    } else if !x.is_finite() {
        if x.is_nan() { "NaN".to_string() } else if x > 0.0 { "Infinity".to_string() } else { "-Infinity".to_string() }
    } else {
        let s = format!("{:.*e}", digits, x);
        let (mant, exp) = s.split_once('e').unwrap_or((&s, "0"));
        let e: i32 = exp.parse().unwrap_or(0);
        let sign = if e < 0 { '-' } else { '+' };
        if e.abs() >= 100 {
            format!("{mant}{sign}{:03}", e.abs())
        } else {
            format!("{mant}E{sign}{:02}", e.abs())
        }
    };
    format!("{body:>width$}")
}

/// A `Timestep Size` that cannot be read as a constant.
pub type DtFunction = Box<dyn FnMut(f64) -> f64>;

enum SizeSpec {
    Absent,
    Const(f64),
    Function,
}

fn size_spec(sim: &dyn Params, key: &str) -> SizeSpec {
    match sim.get(key) {
        None => SizeSpec::Absent,
        Some(Value::Dependent(_)) => SizeSpec::Function,
        Some(_) => sim.real(key).map_or(SizeSpec::Absent, SizeSpec::Const),
    }
}

/// First column of a `Real` array value, as `ListGetConstRealArray(...)(i,1)`.
fn real_column(v: &Value) -> Option<Vec<f64>> {
    match v {
        Value::Real { rows, cols, data } if *rows * *cols == data.len() && *cols > 0 => {
            Some((0..*rows).map(|r| data[r * cols]).collect())
        }
        _ => None,
    }
}

/// `GeometricTimesteps` (ElmerSolver.F90:978): redistributes steps so the
/// geometric ratio is met, keeping the interval length. Updates the three
/// arrays in place.
pub fn geometric_timesteps(counts: &mut [i32], sizes: &mut [f64], ratios: &mut [f64]) {
    let eps = f64::EPSILON;
    for i in 0..counts.len() {
        if ratios[i].abs() < eps {
            ratios[i] = 1.0;
        }
        if (ratios[i] - 1.0).abs() < eps {
            continue;
        }
        let n0 = counts[i] as f64;
        let q = 1.0 + (ratios[i] - 1.0) / n0;
        let n = ((1.0 + (q - 1.0) * n0).ln() / q.ln()).round() as i32;
        let dt = n0 * sizes[i] * (1.0 - q) / (1.0 - q.powi(n));
        ratios[i] = q;
        sizes[i] = dt;
        counts[i] = n;
    }
}

/// Iterates the steps of a run. Feed back solver results with
/// [`set_steady_state_reached`](Self::set_steady_state_reached) and
/// [`stop`](Self::stop) between `next()` calls.
pub struct TimeStepper {
    kind: RunKind,
    counts: Vec<i32>,
    sizes: Vec<f64>,
    ratios: Option<Vec<f64>>,
    output: Vec<i32>,
    coupled_min_iter: i32,
    total_steps: usize,
    time_period: Option<f64>,
    n_periodic: i32,
    start_zero: bool,
    first_size: SizeSpec,
    size: SizeSpec,
    size_iterations: i32,
    dt_fn: Option<DtFunction>,
    pred_corr: Option<PredCorrControl>,
    time_error: f64,
    // state
    interval: usize,
    timestep: usize,
    cum: usize,
    time: f64,
    dt: f64,
    prev_sizes: [f64; 5],
    size_prev: f64,
    steady_reached: bool,
    stopped: bool,
    pending_exit: bool,
    last_saved: bool,
    error: Option<TimeError>,
}

impl TimeStepper {
    /// Reads the run description from the `Simulation` section.
    pub fn new(sim: &dyn Params, kind: RunKind) -> Result<TimeStepper, TimeError> {
        let cfg = |m: String| TimeError::Config(m);
        let mut ratios = None;
        let (counts, sizes) = if kind.is_time() {
            let counts = match sim.get("timestep intervals") {
                Some(Value::Integer(v)) => v.clone(),
                _ => {
                    return Err(cfg(
                        "Keyword > Timestep Intervals < MUST be defined for transient and scanning simulations"
                            .into(),
                    ))
                }
            };
            let sizes = match sim.get("timestep sizes").and_then(real_column) {
                Some(s) => s,
                None => {
                    if kind == RunKind::Scanning || sim.get("timestep size").is_some() {
                        vec![1.0; counts.len()]
                    } else {
                        return Err(cfg(
                            "Keyword [Timestep Sizes] MUST be defined for time dependent simulations"
                                .into(),
                        ));
                    }
                }
            };
            if sizes.len() < counts.len() {
                return Err(cfg(format!(
                    "Timestep Sizes has {} entries for {} intervals",
                    sizes.len(),
                    counts.len()
                )));
            }
            ratios = sim.get("timestep ratios").and_then(real_column);
            if let Some(r) = &ratios {
                if r.len() < counts.len() {
                    return Err(cfg("Timestep Ratios shorter than Timestep Intervals".into()));
                }
            }
            (counts, sizes)
        } else {
            let n = sim.integer("steady state max iterations").unwrap_or(1);
            if n < 1 {
                return Err(cfg("Steady State Max Iterations must be at least 1".into()));
            }
            (vec![n], vec![1.0])
        };
        let mut counts = counts;
        let mut sizes = sizes;
        if let Some(r) = ratios.as_mut() {
            let m = counts.len();
            geometric_timesteps(&mut counts, &mut sizes[..m], &mut r[..m]);
        }
        let coupled_min_iter = sim.integer("steady state min iterations").unwrap_or(1);
        let output = match sim.get("output intervals") {
            Some(Value::Integer(v)) => {
                if v.len() != counts.len() {
                    return Err(cfg(
                        "> Output Intervals < should have the same size as > Timestep Intervals < !"
                            .into(),
                    ));
                }
                v.clone()
            }
            _ => vec![1; counts.len()],
        };
        let total_steps = counts.iter().map(|&c| c.max(0) as usize).sum();
        let size = if sim.get("timestep function").is_some() {
            size_spec(sim, "timestep function")
        } else {
            size_spec(sim, "timestep size")
        };
        Ok(TimeStepper {
            kind,
            counts,
            sizes,
            ratios,
            output,
            coupled_min_iter,
            total_steps,
            time_period: sim.real("time period"),
            n_periodic: sim.integer("periodic timesteps").unwrap_or(0),
            start_zero: sim.logical("timestep start zero").unwrap_or(false),
            first_size: size_spec(sim, "first timestep size"),
            size,
            size_iterations: sim.integer("timestep size iterations").unwrap_or(0),
            dt_fn: None,
            pred_corr: None,
            time_error: 0.0,
            interval: 0,
            timestep: 0,
            cum: 0,
            time: 0.0,
            dt: 1.0,
            prev_sizes: [0.0; 5],
            size_prev: 0.0,
            steady_reached: false,
            stopped: false,
            pending_exit: false,
            last_saved: false,
            error: None,
        })
    }

    /// Evaluator for a `Timestep Size` (or `First Timestep Size`) that depends
    /// on `time`. It is called with the value of `time` to evaluate at.
    pub fn set_dt_function(&mut self, f: DtFunction) {
        self.dt_fn = Some(f);
    }

    /// Turns on `Predictor-Corrector Control`.
    pub fn set_predictor_corrector(&mut self, c: PredCorrControl) {
        self.pred_corr = Some(c);
    }

    pub fn predictor_corrector(&self) -> Option<&PredCorrControl> {
        self.pred_corr.as_ref()
    }

    /// Max-norm of `Values - PrevValues(:,1)` after the last step, for the
    /// predictor-corrector controller.
    pub fn set_time_error(&mut self, e: f64) {
        self.time_error = e;
    }

    /// `SteadyStateReached` after the last step.
    pub fn set_steady_state_reached(&mut self, reached: bool) {
        self.steady_reached = reached;
    }

    /// Ends the run after the last step (exit conditions, `Real Time Max`).
    pub fn stop(&mut self) {
        self.stopped = true;
    }

    pub fn total_steps(&self) -> usize {
        self.total_steps
    }

    /// Steps per interval after any `Timestep Ratios` redistribution.
    pub fn intervals(&self) -> &[i32] {
        &self.counts
    }

    /// Step sizes per interval after redistribution.
    pub fn interval_sizes(&self) -> &[f64] {
        &self.sizes
    }

    /// Whether a save is due after `step`, given the last
    /// [`set_steady_state_reached`](Self::set_steady_state_reached). Records it
    /// for [`needs_final_save`](Self::needs_final_save).
    pub fn should_save(&mut self, step: &Step) -> bool {
        let oi = self.output.get(step.interval - 1).copied().unwrap_or(1);
        let due = oi != 0 && (step.output_due || self.steady_reached);
        self.last_saved = due;
        due
    }

    /// The last step is always saved: true when the final step was not.
    pub fn needs_final_save(&self) -> bool {
        !self.last_saved
    }

    /// The error that ended the iteration early, if any.
    pub fn take_error(&mut self) -> Option<TimeError> {
        self.error.take()
    }

    fn eval_fn(&mut self, t: f64) -> Result<f64, TimeError> {
        match self.dt_fn.as_mut() {
            Some(f) => Ok(f(t)),
            None => Err(TimeError::Config(
                "Timestep Size depends on a variable: supply TimeStepper::set_dt_function".into(),
            )),
        }
    }

    fn pick_dt(&mut self) -> Result<(), TimeError> {
        let mut got = false;
        if self.cum == 1 {
            match self.first_size {
                SizeSpec::Const(v) => {
                    self.dt = v;
                    got = true;
                }
                SizeSpec::Function => {
                    self.dt = self.eval_fn(self.time)?;
                    got = true;
                }
                SizeSpec::Absent => {}
            }
        }
        if self.kind.is_time() && !got {
            let t0 = self.time;
            let dtfunc = match self.size {
                SizeSpec::Absent => None,
                SizeSpec::Const(v) => Some(v),
                SizeSpec::Function => {
                    let mut d = self.eval_fn(t0)?;
                    for _ in 0..self.size_iterations.max(0) {
                        d = self.eval_fn(t0 + d)?;
                    }
                    Some(d)
                }
            };
            match dtfunc {
                Some(d) => {
                    if d < f64::EPSILON {
                        return Err(TimeError::TimestepTooSmall(d));
                    }
                    self.dt = d;
                }
                None => {
                    let mut dt = self.sizes[self.interval];
                    if let Some(r) = &self.ratios {
                        let q = r[self.interval];
                        if (1.0 - q).abs() > f64::EPSILON {
                            dt *= q.powi(self.timestep as i32 - 1);
                        }
                    }
                    self.dt = dt;
                }
            }
        }
        Ok(())
    }
}

fn fortran_modulo(a: f64, p: f64) -> f64 {
    a - (a / p).floor() * p
}

impl Iterator for TimeStepper {
    type Item = Step;

    fn next(&mut self) -> Option<Step> {
        if self.stopped || self.error.is_some() {
            return None;
        }
        // `IF (SteadyStateReached .AND. .NOT.(Transient.OR.Scanning)) THEN
        //    IF (Timestep >= CoupledMinIter) EXIT`
        if self.cum > 0
            && self.steady_reached
            && !self.kind.is_time()
            && self.timestep as i32 >= self.coupled_min_iter
        {
            self.pending_exit = true;
        }
        if self.pending_exit {
            self.pending_exit = false;
            self.interval += 1;
            self.timestep = 0;
        }
        while self.interval < self.counts.len() && self.timestep as i32 >= self.counts[self.interval] {
            self.interval += 1;
            self.timestep = 0;
        }
        if self.interval >= self.counts.len() {
            return None;
        }
        self.timestep += 1;
        self.cum += 1;

        if let Err(e) = self.pick_dt() {
            self.error = Some(e);
            return None;
        }
        if let Some(pc) = self.pred_corr.as_mut() {
            self.dt = pc.control(self.timestep, self.dt, self.time_error);
        }

        if !(self.cum == 1 && self.start_zero) {
            self.time += self.dt;
        }

        let mut periodic = None;
        let mut period = self.time_period.unwrap_or(f64::MAX);
        if self.n_periodic > 0 {
            period = self.n_periodic as f64 * self.dt;
        }
        if self.time_period.is_some() || self.n_periodic > 0 {
            periodic = Some((self.time / period, fortran_modulo(self.time, period)));
        }

        if self.timestep > 1 || self.interval > 0 {
            for i in (1..self.prev_sizes.len()).rev() {
                self.prev_sizes[i] = self.prev_sizes[i - 1];
            }
            self.prev_sizes[0] = self.size_prev;
        }
        self.size_prev = self.dt;

        let oi = self.output[self.interval];
        let output_due = oi != 0 && (self.timestep as i64 - 1) % oi as i64 == 0;

        Some(Step {
            cum_step: self.cum,
            interval: self.interval + 1,
            timestep: self.timestep,
            dt: self.dt,
            time: self.time,
            total_steps: self.total_steps,
            prev_sizes: self.prev_sizes,
            output_due,
            periodic,
            kind: self.kind,
        })
    }
}
