//! Choice of time discretisation: `InitializeTimestep` (SolverBasics.F90:2623)
//! and the solver setup in MainUtils.F90:1431 reduced to a value type.

use crate::error::TimeError;
use elmer_core::Params;

/// Coefficients of one fractional-step stage (`FractionalStep`, TimeIntegrate.F90:~520).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FsCoeffs {
    pub mass_coeff: f64,
    pub force_coeff: f64,
}

impl FsCoeffs {
    /// Reads `fsstep`, `fsTheta`, `fsdTheta`, `fsAlpha`, `fsBeta` from the solver
    /// section. A missing key reads as 0 like `ListGetConstReal`. The caller
    /// re-reads this whenever the solver advances `fsstep`.
    pub fn from_params(solver: &dyn Params) -> Result<Self, TimeError> {
        let g = |k: &str| solver.real(k).unwrap_or(0.0);
        let (step, theta, dtheta, alpha, beta) = (
            g("fsstep"),
            g("fstheta"),
            g("fsdtheta"),
            g("fsalpha"),
            g("fsbeta"),
        );
        // `INT(fsstep)` truncates toward zero.
        match step as i64 {
            1 | 3 => Ok(Self {
                mass_coeff: alpha * theta,
                force_coeff: beta * theta,
            }),
            2 => Ok(Self {
                mass_coeff: beta * dtheta,
                force_coeff: alpha * dtheta,
            }),
            _ => Err(TimeError::InvalidFsStep(step)),
        }
    }
}

/// A time discretisation. The scheme is data; the numerical work is in
/// [`crate::add_1st_order_time`], [`crate::add_2nd_order_time`] and the `_crs`
/// functions.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Scheme {
    /// `Timestepping Method = none`: nothing is added.
    None,
    /// Implicit Euler is `NewmarkBeta(1.0)`, Crank-Nicolson `NewmarkBeta(0.5)`,
    /// explicit Euler `NewmarkBeta(0.0)`.
    NewmarkBeta(f64),
    /// Constant-step BDF of order 1..=5. When the previous step sizes differ
    /// from the current one the first-order dispatch switches to variable-step
    /// BDF, as `Add1stOrderTime` does.
    Bdf(u8),
    /// Variable-step BDF, forced. Implemented to order 3 only (Fortran warns above that).
    Vbdf(u8),
    /// One stage of the Runge-Kutta scheme. The stage loop belongs to the scheduler.
    RungeKutta,
    FractionalStep(FsCoeffs),
    /// `order` is `Predictor-Corrector Scheme Order`; it is capped at `done_time / 2`.
    AdamsBashforth { zeta: f64, order: u8 },
    AdamsMoulton { order: u8 },
    /// Second order in time with `Alpha = alpha`, `Beta = (1-alpha)^2/4`, `Gamma = 1/2 - alpha`.
    /// `alpha = 0` is Newmark's average acceleration scheme.
    Bossak { alpha: f64 },
    /// Second order in time, `rho_inf` is `Generalized-alpha rinf`.
    GeneralizedAlpha { rho_inf: f64 },
}

const METHOD_KEY: &str = "timestepping method";

const KNOWN_METHODS: &[&str] = &[
    "implicit euler",
    "explicit euler",
    "runge-kutta",
    "crank-nicolson",
    "fs",
    "adams-bashforth",
    "adams-moulton",
    "newmark",
    "bdf",
    "generalized-alpha",
    "bossak",
    "none",
];

fn lookup_real(solver: &dyn Params, sim: &dyn Params, key: &str) -> Option<f64> {
    solver.real(key).or_else(|| sim.real(key))
}

fn lookup_int(solver: &dyn Params, sim: &dyn Params, key: &str) -> Option<i32> {
    solver.integer(key).or_else(|| sim.integer(key))
}

impl Scheme {
    /// Resolves the scheme for a solver from its section and the `Simulation`
    /// section, with `time_order` the solver's `Time Derivative Order` (1 or 2).
    ///
    /// Defaults follow the Fortran: no method and no `Newmark Beta` gives
    /// implicit Euler (first order) or Bossak with `alpha = -0.05` (second
    /// order); `BDF Order` defaults to 2 (Fortran warns); `Generalized-alpha rinf`
    /// defaults to 0.75; `Adams Zeta` to 1; the predictor-corrector order to 2.
    pub fn from_params(
        solver: &dyn Params,
        simulation: &dyn Params,
        time_order: u8,
    ) -> Result<Scheme, TimeError> {
        let method = solver
            .string(METHOD_KEY)
            .or_else(|| simulation.string(METHOD_KEY))
            .map(|m| m.trim().to_ascii_lowercase());

        if let Some(m) = &method {
            if !KNOWN_METHODS.contains(&m.as_str()) {
                return Err(TimeError::UnknownMethod(m.clone()));
            }
            if m == "none" {
                return Ok(Scheme::None);
            }
        }

        if time_order >= 2 {
            return Ok(match method.as_deref() {
                Some("generalized-alpha") => Scheme::GeneralizedAlpha {
                    rho_inf: lookup_real(solver, simulation, "generalized-alpha rinf")
                        .unwrap_or(0.75),
                },
                // Every other method name selects Bossak at second order.
                _ => Scheme::Bossak {
                    alpha: lookup_real(solver, simulation, "bossak alpha").unwrap_or(-0.05),
                },
            });
        }

        let newmark_beta = || lookup_real(solver, simulation, "newmark beta");
        Ok(match method.as_deref() {
            None => Scheme::NewmarkBeta(newmark_beta().unwrap_or(1.0)),
            Some("implicit euler") | Some("generalized-alpha") | Some("bossak") => {
                Scheme::NewmarkBeta(1.0)
            }
            Some("explicit euler") => Scheme::NewmarkBeta(0.0),
            Some("crank-nicolson") => Scheme::NewmarkBeta(0.5),
            // A missing `Newmark Beta` reads as 0 (ListGetConstReal); Fortran only warns outside [0,1].
            Some("newmark") => Scheme::NewmarkBeta(newmark_beta().unwrap_or(0.0)),
            Some("runge-kutta") => Scheme::RungeKutta,
            Some("fs") => Scheme::FractionalStep(FsCoeffs::from_params(solver)?),
            Some("adams-bashforth") => Scheme::AdamsBashforth {
                zeta: solver.real("adams zeta").unwrap_or(1.0),
                order: pred_corr_order(solver, simulation),
            },
            Some("adams-moulton") => Scheme::AdamsMoulton {
                order: pred_corr_order(solver, simulation),
            },
            Some("bdf") => {
                let o = lookup_int(solver, simulation, "bdf order").unwrap_or(2);
                if !(1..=5).contains(&o) {
                    return Err(TimeError::InvalidBdfOrder(o as i64));
                }
                Scheme::Bdf(o as u8)
            }
            Some(other) => return Err(TimeError::UnknownMethod(other.to_string())),
        })
    }

    /// `InitializeTimestep` forces `Beta = 1` on the first step of a first-order
    /// scheme unless `Beta = 0`. Call with the solver's `DoneTime` for the step.
    pub fn for_step(self, done_time: usize) -> Scheme {
        match self {
            Scheme::NewmarkBeta(b) if done_time == 1 && b != 0.0 => Scheme::NewmarkBeta(1.0),
            s => s,
        }
    }

    /// The `(param_a, param_b)` pair `Time2ndOrder` receives: `(alpha, -1)` for
    /// Bossak, `(unused, rho_inf)` for generalized-alpha.
    pub fn second_order_params(&self) -> Option<(f64, f64)> {
        match *self {
            Scheme::Bossak { alpha } => Some((alpha, -1.0)),
            Scheme::GeneralizedAlpha { rho_inf } => Some((0.0, rho_inf)),
            _ => None,
        }
    }
}

fn pred_corr_order(solver: &dyn Params, sim: &dyn Params) -> u8 {
    lookup_int(solver, sim, "predictor-corrector scheme order").unwrap_or(2) as u8
}

/// `MIN(PredCorrOrder, DoneTime/2)`, with integer division.
pub fn effective_pred_corr_order(order: u8, done_time: usize) -> usize {
    (order as usize).min(done_time / 2)
}
