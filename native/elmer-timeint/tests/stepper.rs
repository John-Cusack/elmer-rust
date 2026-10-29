//! Time-step sequence, adaptive control and scheme selection.
//!
//! `fixtures/*.time.txt` are the `Time:` lines printed by the pristine
//! ElmerSolver for `fixtures/*.sif.txt` (adv_diff1 from fem/tests, and the same
//! case with three intervals, `Timestep Ratios` and `Output Intervals`).

use elmer_core::{ParamMap, Value};
use elmer_timeint::adaptive::{pred_corr_error_estimate, time_step_controller, time_step_limiter};
use elmer_timeint::{
    fortran_es, geometric_timesteps, Decision, PredCorrConfig, PredCorrControl, RunKind, Scheme,
    SubstepConfig, SubstepController, SubstepMode, TimeError, TimeStepper,
};

fn reals(v: &[f64]) -> Value {
    Value::Real {
        rows: v.len(),
        cols: 1,
        data: v.to_vec(),
    }
}

fn int(v: &[i32]) -> Value {
    Value::Integer(v.to_vec())
}

fn transient(counts: &[i32], sizes: &[f64]) -> ParamMap {
    let mut p = ParamMap::new();
    p.insert("Simulation Type", Value::Str("Transient".into()));
    p.insert("Timestep Intervals", int(counts));
    p.insert("Timestep Sizes", reals(sizes));
    p
}

fn lines(file: &str) -> Vec<String> {
    file.lines()
        .map(|l| l.trim_start_matches("MAIN: ").to_string())
        .collect()
}

#[test]
fn uniform_steps_match_the_fortran_time_lines() {
    let p = transient(&[10], &[2.0]);
    let got: Vec<String> = TimeStepper::new(&p, RunKind::Transient)
        .unwrap()
        .map(|s| s.time_line())
        .collect();
    assert_eq!(got, lines(include_str!("fixtures/adv_diff1.time.txt")));
}

#[test]
fn intervals_ratios_and_output_match_the_fortran_time_lines() {
    let mut p = transient(&[3, 4, 2], &[1.0, 0.5, 2.0]);
    p.insert("Timestep Ratios", reals(&[1.0, 1.5, 0.0]));
    p.insert("Output Intervals", int(&[1, 2, 0]));
    let mut st = TimeStepper::new(&p, RunKind::Transient).unwrap();
    // 4 steps at ratio 1.5 become 3, so the run has 3 + 3 + 2 = 8 steps.
    assert_eq!(st.intervals(), &[3, 3, 2]);
    let mut got = Vec::new();
    let mut saves = Vec::new();
    while let Some(s) = st.next() {
        got.push(s.time_line());
        saves.push(st.should_save(&s));
    }
    assert_eq!(got, lines(include_str!("fixtures/multi.time.txt")));
    // Interval 1 saves every step, interval 2 every second (1st, 3rd), interval 3 never.
    assert_eq!(saves, [true, true, true, true, false, true, false, false]);
    assert!(st.needs_final_save());
}

#[test]
fn geometric_redistribution_keeps_the_interval_length() {
    let (mut n, mut dt, mut r) = ([4], [0.5], [1.5]);
    geometric_timesteps(&mut n, &mut dt, &mut r);
    assert_eq!(n, [3]);
    let q: f64 = r[0];
    let total = dt[0] * (1.0 + q + q * q);
    assert!((total - 2.0).abs() < 1e-12, "{total}");
    // Zero ratio means one; ratio one leaves everything alone.
    let (mut n, mut dt, mut r) = ([5, 5], [1.0, 1.0], [0.0, 1.0]);
    geometric_timesteps(&mut n, &mut dt, &mut r);
    assert_eq!((n, dt, r), ([5, 5], [1.0, 1.0], [1.0, 1.0]));
}

#[test]
fn fortran_es_format() {
    assert_eq!(fortran_es(20.0, 12, 3), "   2.000E+01");
    assert_eq!(fortran_es(0.0, 12, 3), "   0.000E+00");
    assert_eq!(fortran_es(-1.5e-5, 12, 3), "  -1.500E-05");
    assert_eq!(fortran_es(9.9996, 12, 3), "   1.000E+01");
    assert_eq!(fortran_es(1.0e100, 12, 3), "   1.000+100");
}

#[test]
fn step_fields_prev_sizes_and_errors() {
    let p = transient(&[2, 2], &[0.5, 1.0]);
    let steps: Vec<_> = TimeStepper::new(&p, RunKind::Transient).unwrap().collect();
    let dts: Vec<f64> = steps.iter().map(|s| s.dt).collect();
    assert_eq!(dts, [0.5, 0.5, 1.0, 1.0]);
    let times: Vec<f64> = steps.iter().map(|s| s.time).collect();
    assert_eq!(times, [0.5, 1.0, 2.0, 3.0]);
    assert_eq!(steps[3].prev_sizes, [1.0, 0.5, 0.5, 0.0, 0.0]);
    assert_eq!(steps[2].interval, 2);
    assert_eq!(steps[2].timestep, 1);

    let mut bad = ParamMap::new();
    bad.insert("Simulation Type", Value::Str("Transient".into()));
    assert!(TimeStepper::new(&bad, RunKind::Transient).is_err());
    let mut p = transient(&[2], &[1.0]);
    p.insert("Output Intervals", int(&[1, 1]));
    assert!(TimeStepper::new(&p, RunKind::Transient).is_err());
}

#[test]
fn first_timestep_size_and_start_zero() {
    let mut p = transient(&[3], &[1.0]);
    p.insert("First Timestep Size", Value::real(0.25));
    let t: Vec<f64> = TimeStepper::new(&p, RunKind::Transient)
        .unwrap()
        .map(|s| s.time)
        .collect();
    assert_eq!(t, [0.25, 1.25, 2.25]);

    let mut p = transient(&[3], &[1.0]);
    p.insert("Timestep Start Zero", Value::Logical(true));
    let t: Vec<f64> = TimeStepper::new(&p, RunKind::Transient)
        .unwrap()
        .map(|s| s.time)
        .collect();
    assert_eq!(t, [0.0, 1.0, 2.0]);
}

#[test]
fn time_dependent_timestep_size_iterates_for_consistency() {
    use elmer_core::{Dependent, DependentKind};
    let mut p = ParamMap::new();
    p.insert("Timestep Intervals", int(&[3]));
    p.insert(
        "Timestep Size",
        Value::Dependent(Dependent {
            variables: vec!["time".into()],
            kind: DependentKind::Matc("dt".into()),
            rows: 1,
            cols: 1,
        }),
    );
    p.insert("Timestep Size Iterations", Value::Integer(vec![2]));
    let mut st = TimeStepper::new(&p, RunKind::Transient).unwrap();
    assert!(st.next().is_none() && st.take_error().is_some(), "needs a function");

    let mut st = TimeStepper::new(&p, RunKind::Transient).unwrap();
    st.set_dt_function(Box::new(|t| 0.1 + 0.1 * t));
    let steps: Vec<_> = st.by_ref().collect();
    // dt = 0.1 + 0.1 (t0 + dt) solved by two fixed-point sweeps from dt(t0).
    let mut t = 0.0f64;
    for s in &steps {
        let mut d = 0.1 + 0.1 * t;
        d = 0.1 + 0.1 * (t + d);
        d = 0.1 + 0.1 * (t + d);
        t += d;
        assert_eq!(s.time.to_bits(), t.to_bits());
    }
}

#[test]
fn steady_state_stops_an_interval_after_min_iterations() {
    let mut p = ParamMap::new();
    p.insert("Steady State Max Iterations", Value::Integer(vec![10]));
    p.insert("Steady State Min Iterations", Value::Integer(vec![3]));
    let mut st = TimeStepper::new(&p, RunKind::Steady).unwrap();
    let mut n = 0;
    while let Some(s) = st.next() {
        n += 1;
        assert_eq!(s.dt, 1.0);
        st.set_steady_state_reached(s.timestep >= 2);
    }
    assert_eq!(n, 3, "reached at step 2 but exits only once timestep >= min");

    let st = TimeStepper::new(&p, RunKind::Steady).unwrap();
    assert_eq!(st.count(), 10);
}

#[test]
fn stop_ends_the_run() {
    let p = transient(&[10], &[1.0]);
    let mut st = TimeStepper::new(&p, RunKind::Transient).unwrap();
    st.next().unwrap();
    st.stop();
    assert!(st.next().is_none());
}

#[test]
fn periodic_time_reports_cycle_and_phase() {
    let mut p = transient(&[5], &[0.4]);
    p.insert("Time Period", Value::real(1.0));
    let steps: Vec<_> = TimeStepper::new(&p, RunKind::Transient).unwrap().collect();
    let (cycle, phase) = steps[2].periodic.unwrap();
    assert!((cycle - 1.2).abs() < 1e-12 && (phase - 0.2).abs() < 1e-12);
}

// ---- adaptive sub-stepping ---------------------------------------------------

fn adaptive_params(extra: impl FnOnce(&mut ParamMap)) -> SubstepConfig {
    let mut p = ParamMap::new();
    p.insert("Adaptive Timestepping", Value::Logical(true));
    p.insert("Adaptive Time Error", Value::real(1e-3));
    extra(&mut p);
    SubstepConfig::from_params(&p).unwrap().unwrap()
}

#[test]
fn adaptive_requires_an_error_limit_and_is_off_by_default() {
    assert!(SubstepConfig::from_params(&ParamMap::new()).unwrap().is_none());
    let mut p = ParamMap::new();
    p.insert("Adaptive Timestepping", Value::Logical(true));
    assert!(SubstepConfig::from_params(&p).is_err());
}

#[test]
fn adaptive_rejects_halves_then_accepts_and_grows() {
    let mut c = SubstepController::new(adaptive_params(|p| {
        p.insert("Adaptive Keep Smallest", Value::Integer(vec![1]));
    }));
    c.begin_step(2.0, 10.0);
    // First trial uses the whole step.
    let s = c.next_substep().unwrap();
    assert_eq!((s.dt, s.time), (2.0, 12.0));
    // Error too large: halve and retry.
    assert_eq!(c.finish_error(1.0), Decision::Reject { next: 1.0 });
    let s = c.next_substep().unwrap();
    assert_eq!((s.dt, s.time), (1.0, 11.0));
    // Accepted, but KeepSmallest = 1 holds the size for one accepted sub-step.
    assert_eq!(c.finish_error(1e-6), Decision::Accept { used: 1.0, next: 1.0 });
    let s = c.next_substep().unwrap();
    assert_eq!((s.dt, s.time), (1.0, 12.0));
    // Second acceptance finishes the step; ddt then grows (capped at 2).
    assert_eq!(c.finish_error(1e-6), Decision::Accept { used: 1.0, next: 2.0 });
    assert!(c.next_substep().is_none());
    assert_eq!(c.cumulative(), 2.0);

    // Next nominal step starts from the grown ddt, not from scratch.
    c.begin_step(2.0, 12.0);
    assert_eq!(c.next_substep().unwrap().dt, 2.0);
}

#[test]
fn adaptive_accepts_at_the_minimum_step_whatever_the_error() {
    let mut c = SubstepController::new(adaptive_params(|p| {
        p.insert("Adaptive Min Timestep", Value::real(0.5));
        p.insert("Adaptive Decrease Coefficient", Value::real(0.25));
    }));
    c.begin_step(1.0, 0.0);
    assert_eq!(c.next_substep().unwrap().dt, 1.0);
    assert_eq!(c.finish_error(5.0), Decision::Reject { next: 0.25 });
    assert_eq!(c.next_substep().unwrap().dt, 0.25);
    assert!(matches!(c.finish_error(5.0), Decision::Accept { used, .. } if used == 0.25));
}

#[test]
fn smart_and_rough_substepping() {
    // Smart: 1.0 left and ddt 0.75 would leave 0.25, so split into 0.5 + 0.5.
    let mut c = SubstepController::new(adaptive_params(|p| {
        p.insert("Adaptive Smart Timestep", Value::Logical(true));
        p.insert("Adaptive Max Timestep", Value::real(0.75));
    }));
    c.begin_step(1.0, 0.0);
    assert_eq!(c.next_substep().unwrap().dt, 0.5);
    // Rough: the last sub-step may overshoot the step.
    let mut c = SubstepController::new(adaptive_params(|p| {
        p.insert("Adaptive Rough Timestep", Value::Logical(true));
        p.insert("Adaptive Max Timestep", Value::real(0.75));
        p.insert("Adaptive Keep Smallest", Value::Integer(vec![0]));
    }));
    c.begin_step(1.0, 0.0);
    assert_eq!(c.next_substep().unwrap().dt, 0.75);
    c.finish_error(0.0);
    assert_eq!(c.next_substep().unwrap().dt, 0.75);
}

#[test]
fn divergence_control_shrinks_then_gives_up_at_the_minimum() {
    let mut p = ParamMap::new();
    p.insert("Convergence Control", Value::Logical(true));
    p.insert("Adaptive Min Timestep", Value::real(0.2));
    let cfg = SubstepConfig::from_params(&p).unwrap().unwrap();
    assert_eq!(cfg.mode, SubstepMode::Divergence);
    let mut c = SubstepController::new(cfg);
    c.begin_step(1.0, 0.0);
    assert_eq!(c.next_substep().unwrap().dt, 1.0);
    assert_eq!(c.finish_divergence(true).unwrap(), Decision::Reject { next: 0.5 });
    assert_eq!(c.next_substep().unwrap().dt, 0.5);
    assert_eq!(c.finish_divergence(false).unwrap(), Decision::Accept { used: 0.5, next: 1.0 });
    assert_eq!(c.next_substep().unwrap().dt, 0.5);
    assert!(matches!(c.finish_divergence(false), Ok(Decision::Accept { .. })));
    assert!(c.next_substep().is_none());

    // Diverging all the way down to the minimum is fatal.
    let mut c = SubstepController::new(SubstepConfig::from_params(&p).unwrap().unwrap());
    c.begin_step(1.0, 0.0);
    let mut last = Ok(Decision::Reject { next: 1.0 });
    for _ in 0..10 {
        if c.next_substep().is_none() {
            break;
        }
        last = c.finish_divergence(true);
        if last.is_err() {
            break;
        }
    }
    assert!(matches!(last, Err(TimeError::NoStableTimestep)));
}

// ---- predictor-corrector ---------------------------------------------------

#[test]
fn pi_controller_pieces() {
    // atan limiter: g = 1 keeps the step, large g saturates at 1 + pi.
    assert_eq!(time_step_limiter(0.2, 1.0), 0.2);
    assert!((time_step_limiter(1.0, 1e12) - (1.0 + std::f64::consts::PI)).abs() < 1e-9);
    assert_eq!(pred_corr_error_estimate(0.0, 2, 1.0, 1.0), 0.0);
    assert_eq!(pred_corr_error_estimate(0.5, 1, 1.0, 1.0), 1.0);
    let e = pred_corr_error_estimate(0.5, 2, 3.0, 2.0);
    assert_eq!(e, 3.0 * 2.0 / 0.5 / 3.0 / 3.0);
    assert_eq!(time_step_controller(0.1, 0.0, 1.0, 1e-6, 0.2, -0.07), 0.1);
}

#[test]
fn predictor_corrector_sequence() {
    let cfg = PredCorrConfig {
        order: 2,
        epsilon: 1e-6,
        beta1: 0.2,
        beta2: -0.4 / 6.0,
    };
    let mut c = PredCorrControl::new(cfg);
    assert_eq!(c.control(1, 0.1, 0.0), 0.0);
    assert_eq!(c.control(2, 0.1, 0.0), 0.1);
    assert_eq!(c.zeta(), 1.0);
    // Step 3 uses a first-order estimate and etaOld = eta, so g = (eps/eta)^(b1+b2).
    let err = 4e-6;
    let eta = err / 0.1 / 2.0;
    let g = (1e-6f64 / eta).powf(cfg.beta1) * (1e-6f64 / eta).powf(cfg.beta2);
    let want = time_step_limiter(0.1, g);
    let got = c.control(3, 0.1, err);
    assert_eq!(got.to_bits(), want.to_bits());
    assert_eq!(c.zeta().to_bits(), (want / 0.1).to_bits());
    // Step 4 uses the configured order and the saved eta.
    let eta2 = pred_corr_error_estimate(want, 2, 2e-6, want / 0.1);
    let g2 = (1e-6f64 / eta2).powf(cfg.beta1) * (1e-6f64 / eta).powf(cfg.beta2);
    assert_eq!(c.control(4, 0.1, 2e-6).to_bits(), time_step_limiter(want, g2).to_bits());
}

#[test]
fn stepper_applies_the_predictor_corrector() {
    let p = transient(&[4], &[0.1]);
    let mut st = TimeStepper::new(&p, RunKind::Transient).unwrap();
    st.set_predictor_corrector(PredCorrControl::new(PredCorrConfig {
        order: 2,
        epsilon: 1e-6,
        beta1: 0.2,
        beta2: -0.07,
    }));
    let s1 = st.next().unwrap();
    assert_eq!((s1.dt, s1.time), (0.0, 0.0));
    let s2 = st.next().unwrap();
    assert_eq!(s2.dt, 0.1);
    st.set_time_error(1e-6);
    let s3 = st.next().unwrap();
    assert!(s3.dt > 0.0 && s3.time > s2.time);
}

// ---- scheme selection --------------------------------------------------------

fn method(m: &str) -> ParamMap {
    let mut p = ParamMap::new();
    p.insert("Timestepping Method", Value::Str(m.into()));
    p
}

#[test]
fn scheme_from_params() {
    let none = ParamMap::new();
    let s = |solver: &ParamMap, sim: &ParamMap, order| Scheme::from_params(solver, sim, order).unwrap();
    assert_eq!(s(&none, &none, 1), Scheme::NewmarkBeta(1.0));
    assert_eq!(s(&none, &none, 2), Scheme::Bossak { alpha: -0.05 });
    assert_eq!(s(&method("implicit euler"), &none, 1), Scheme::NewmarkBeta(1.0));
    assert_eq!(s(&method("Crank-Nicolson"), &none, 1), Scheme::NewmarkBeta(0.5));
    assert_eq!(s(&method("explicit euler"), &none, 1), Scheme::NewmarkBeta(0.0));
    // Method from the Simulation section when the solver has none.
    assert_eq!(s(&none, &method("bdf"), 1), Scheme::Bdf(2));
    let mut bdf = method("bdf");
    bdf.insert("BDF Order", Value::Integer(vec![4]));
    assert_eq!(s(&bdf, &none, 1), Scheme::Bdf(4));
    bdf.insert("BDF Order", Value::Integer(vec![6]));
    assert!(matches!(
        Scheme::from_params(&bdf, &none, 1),
        Err(TimeError::InvalidBdfOrder(6))
    ));
    let mut nm = method("newmark");
    nm.insert("Newmark Beta", Value::real(0.3));
    assert_eq!(s(&nm, &none, 1), Scheme::NewmarkBeta(0.3));
    let mut sim = ParamMap::new();
    sim.insert("Newmark Beta", Value::real(0.6));
    assert_eq!(s(&none, &sim, 1), Scheme::NewmarkBeta(0.6));
    assert_eq!(s(&method("generalized-alpha"), &none, 2), Scheme::GeneralizedAlpha { rho_inf: 0.75 });
    let mut ba = method("bossak");
    ba.insert("Bossak Alpha", Value::real(-0.1));
    assert_eq!(s(&ba, &none, 2), Scheme::Bossak { alpha: -0.1 });
    assert_eq!(s(&method("none"), &none, 1), Scheme::None);
    assert_eq!(s(&method("runge-kutta"), &none, 1), Scheme::RungeKutta);
    assert!(matches!(
        Scheme::from_params(&method("leapfrog"), &none, 1),
        Err(TimeError::UnknownMethod(_))
    ));
    let mut ab = method("adams-bashforth");
    ab.insert("Adams Zeta", Value::real(1.4));
    assert_eq!(s(&ab, &none, 1), Scheme::AdamsBashforth { zeta: 1.4, order: 2 });
}

#[test]
fn fractional_step_coefficients() {
    let mut p = method("fs");
    for (k, v) in [("fsTheta", 0.3), ("fsdTheta", 0.4), ("fsAlpha", 0.5), ("fsBeta", 0.6)] {
        p.insert(k, Value::real(v));
    }
    let fs = |step: f64| {
        let mut q = p.clone();
        q.insert("fsstep", Value::real(step));
        Scheme::from_params(&q, &ParamMap::new(), 1)
    };
    let Scheme::FractionalStep(c1) = fs(1.0).unwrap() else { panic!() };
    assert_eq!((c1.mass_coeff, c1.force_coeff), (0.5 * 0.3, 0.6 * 0.3));
    let Scheme::FractionalStep(c2) = fs(2.0).unwrap() else { panic!() };
    assert_eq!((c2.mass_coeff, c2.force_coeff), (0.6 * 0.4, 0.5 * 0.4));
    assert!(fs(7.0).is_err());
}

#[test]
fn first_step_of_a_newmark_scheme_is_implicit_euler() {
    assert_eq!(Scheme::NewmarkBeta(0.5).for_step(1), Scheme::NewmarkBeta(1.0));
    assert_eq!(Scheme::NewmarkBeta(0.5).for_step(2), Scheme::NewmarkBeta(0.5));
    assert_eq!(Scheme::NewmarkBeta(0.0).for_step(1), Scheme::NewmarkBeta(0.0));
    assert_eq!(Scheme::Bdf(2).for_step(1), Scheme::Bdf(2));
}
