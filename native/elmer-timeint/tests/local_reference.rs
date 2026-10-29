//! Bitwise comparison of the local routines against scalar transliterations
//! of TimeIntegrate.F90 written with Fortran's 1-based `(i,j)` indexing and
//! the Fortran expression text, evaluated left to right.

#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]

mod support;

use elmer_core::LocalSystem;
use elmer_timeint::local::{AdamsState, Aux1st};
use elmer_timeint::scheme::FsCoeffs;
use elmer_timeint::second::{add_2nd_order_time, update_2nd_order, SecondOrderState};
use elmer_timeint::{add_1st_order_time, lump_matrix, Scheme, TimeError};
use support::{assert_bits, Lcg};

const N: usize = 6;

/// Column-major-free 1-based accessor over a row-major slice.
struct Mat<'a>(&'a [f64], usize);
impl Mat<'_> {
    fn at(&self, i: usize, j: usize) -> f64 {
        self.0[(i - 1) * self.1 + (j - 1)]
    }
}
fn set(m: &mut [f64], n: usize, i: usize, j: usize, v: f64) {
    m[(i - 1) * n + (j - 1)] = v;
}

fn system(rng: &mut Lcg, with_damp: bool) -> LocalSystem {
    let mut l = LocalSystem::zeros(N);
    l.stiff = rng.vec(N * N);
    l.mass = Some(rng.vec(N * N).into_iter().map(|x| x + 1.5).collect());
    if with_damp {
        l.damp = Some(rng.vec(N * N));
    }
    l.force = rng.vec(N);
    l
}

fn ref_newmark(dt: f64, mass: &[f64], stiff: &mut [f64], force: &mut [f64], prev: &[f64], beta: f64) {
    let (m, n) = (Mat(mass, N), N);
    for i in 1..=n {
        let mut s = 0.0;
        for j in 1..=n {
            s = s + (1.0 / dt) * m.at(i, j) * prev[j - 1]
                - (1.0 - beta) * Mat(stiff, N).at(i, j) * prev[j - 1];
        }
        for j in 1..=n {
            let k = Mat(stiff, N).at(i, j);
            set(stiff, N, i, j, beta * k + (1.0 / dt) * m.at(i, j));
        }
        force[i - 1] = force[i - 1] + s;
    }
}

fn ref_bdf(dt: f64, mass: &[f64], stiff: &mut [f64], force: &mut [f64], p: &[&[f64]], order: usize) {
    let m = Mat(mass, N);
    for i in 1..=N {
        let mut s = 0.0;
        for j in 1..=N {
            let pj = |k: usize| p[k - 1][j - 1];
            let term = match order {
                1 => (1.0e0 / dt) * m.at(i, j) * pj(1),
                2 => (1.0 / dt) * m.at(i, j) * (2.0 * pj(1) - 0.5 * pj(2)),
                3 => {
                    (1.0e0 / dt) * m.at(i, j) * (3.0 * pj(1) - (3.0 / 2.0) * pj(2) + (1.0 / 3.0) * pj(3))
                }
                4 => {
                    (1.0 / dt)
                        * m.at(i, j)
                        * (4.0 * pj(1) - 3.0 * pj(2) + (4.0 / 3.0) * pj(3) - (1.0 / 4.0) * pj(4))
                }
                5 => {
                    (1.0 / dt)
                        * m.at(i, j)
                        * (5.0 * pj(1) - 5.0 * pj(2) + (10.0 / 3.0) * pj(3) - (5.0 / 4.0) * pj(4)
                            + (1.0 / 5.0) * pj(5))
                }
                _ => unreachable!(),
            };
            s = s + term;
        }
        force[i - 1] = force[i - 1] + s;
        for j in 1..=N {
            let c = match order {
                1 => 1.0e0 / dt,
                2 => 1.5e0 / dt,
                3 => 11.0 / (6.0 * dt),
                4 => 25.0 / (12.0 * dt),
                _ => 137.0 / (60.0 * dt),
            };
            let k = Mat(stiff, N).at(i, j);
            set(stiff, N, i, j, c * m.at(i, j) + k);
        }
    }
}

#[test]
fn newmark_beta_is_bitwise_the_fortran() {
    let mut rng = Lcg(1);
    for &beta in &[1.0, 0.5, 0.0, 0.37] {
        let sys = system(&mut rng, false);
        let prev = rng.vec(N);
        let dt = 0.013;
        let mut want = sys.clone();
        ref_newmark(dt, sys.mass.as_ref().unwrap(), &mut want.stiff, &mut want.force, &prev, beta);
        let mut got = sys.clone();
        add_1st_order_time(&mut got, dt, &Scheme::NewmarkBeta(beta), &[&prev], &mut Aux1st::default())
            .unwrap();
        assert_bits("stiff", &got.stiff, &want.stiff);
        assert_bits("force", &got.force, &want.force);
        assert_eq!(got.mass, sys.mass, "mass is untouched");
    }
}

#[test]
fn constant_step_bdf_orders_one_to_five() {
    let mut rng = Lcg(2);
    for order in 1..=5usize {
        let sys = system(&mut rng, false);
        let hist: Vec<Vec<f64>> = (0..5).map(|_| rng.vec(N)).collect();
        let prev: Vec<&[f64]> = hist.iter().map(|v| v.as_slice()).take(order).collect();
        let dt = 0.021;
        let mut want = sys.clone();
        ref_bdf(dt, sys.mass.as_ref().unwrap(), &mut want.stiff, &mut want.force, &prev, order);
        let mut got = sys.clone();
        let dts = [dt; 4];
        let mut aux = Aux1st {
            dts: &dts,
            ..Aux1st::default()
        };
        add_1st_order_time(&mut got, dt, &Scheme::Bdf(order as u8), &prev, &mut aux).unwrap();
        assert_bits("stiff", &got.stiff, &want.stiff);
        assert_bits("force", &got.force, &want.force);
    }
}

#[test]
fn bdf_startup_uses_the_available_history() {
    // Order 3 requested, one vector available (DoneTime = 1): BDF1.
    let mut rng = Lcg(3);
    let sys = system(&mut rng, false);
    let p1 = rng.vec(N);
    let mut a = sys.clone();
    add_1st_order_time(&mut a, 0.1, &Scheme::Bdf(3), &[&p1], &mut Aux1st::default()).unwrap();
    let mut b = sys.clone();
    add_1st_order_time(&mut b, 0.1, &Scheme::Bdf(1), &[&p1], &mut Aux1st::default()).unwrap();
    assert_bits("stiff", &a.stiff, &b.stiff);
    assert_bits("force", &a.force, &b.force);
}

fn ref_vbdf(dts: &[f64], mass: &[f64], stiff: &mut [f64], force: &mut [f64], p: &[&[f64]], order: usize) {
    let d = |k: usize| dts[k - 1];
    let mut a = [0.0f64; 4];
    a[0] = 1.0 / d(1);
    a[1] = -1.0 / d(1);
    if order >= 2 {
        a[0] = a[0] + 1.0 / (d(1) + d(2));
        a[1] = a[1] - (1.0 + d(1) / d(2)) / (d(1) + d(2));
        a[2] = (d(1) / d(2)) / (d(1) + d(2));
    }
    if order >= 3 {
        a[0] = a[0] + 1.0 / (d(1) + d(2) + d(3));
        a[1] = a[1]
            - (1.0 + d(1) / d(2) * (1.0 + (d(1) + d(2)) / (d(2) + d(3)))) / (d(1) + d(2) + d(3));
        a[2] = a[2]
            + (d(1) / d(2) * (1.0 + (d(1) + d(2)) / (d(2) + d(3)))
                + d(1) / d(3) * (d(1) + d(2)) / (d(2) + d(3)))
                / (d(1) + d(2) + d(3));
        a[3] = -(d(1) / d(3)) * (d(1) + d(2)) / (d(2) + d(3)) / (d(1) + d(2) + d(3));
    }
    let m = Mat(mass, N);
    for i in 1..=N {
        let mut s = 0.0;
        for k in 1..=order.min(3) {
            for j in 1..=N {
                s = s - a[k] * m.at(i, j) * p[k - 1][j - 1];
            }
        }
        force[i - 1] = force[i - 1] + s;
        for j in 1..=N {
            let kk = Mat(stiff, N).at(i, j);
            set(stiff, N, i, j, a[0] * m.at(i, j) + kk);
        }
    }
}

#[test]
fn variable_step_bdf_orders_one_to_three() {
    let mut rng = Lcg(4);
    for order in 1..=3usize {
        let sys = system(&mut rng, false);
        let hist: Vec<Vec<f64>> = (0..3).map(|_| rng.vec(N)).collect();
        let prev: Vec<&[f64]> = hist.iter().map(|v| v.as_slice()).take(order).collect();
        let dt = 0.05;
        let older = [0.03, 0.07];
        let mut dts = vec![dt];
        dts.extend_from_slice(&older);
        let mut want = sys.clone();
        ref_vbdf(&dts, sys.mass.as_ref().unwrap(), &mut want.stiff, &mut want.force, &prev, order);
        // `Bdf` switches to VBDF because the sizes differ.
        let mut got = sys.clone();
        let mut aux = Aux1st {
            dts: &older,
            ..Aux1st::default()
        };
        add_1st_order_time(&mut got, dt, &Scheme::Bdf(order as u8), &prev, &mut aux).unwrap();
        assert_bits("stiff", &got.stiff, &want.stiff);
        assert_bits("force", &got.force, &want.force);
    }
}

#[test]
fn bdf_sizes_within_a_millionth_stay_constant_step() {
    let mut rng = Lcg(5);
    let sys = system(&mut rng, false);
    let h = (rng.vec(N), rng.vec(N));
    let dt = 0.1;
    let near = [dt * (1.0 + 5.0e-7)];
    let mut a = sys.clone();
    let mut aux = Aux1st {
        dts: &near,
        ..Aux1st::default()
    };
    add_1st_order_time(&mut a, dt, &Scheme::Bdf(2), &[&h.0, &h.1], &mut aux).unwrap();
    let mut want = sys.clone();
    ref_bdf(dt, sys.mass.as_ref().unwrap(), &mut want.stiff, &mut want.force, &[&h.0, &h.1], 2);
    assert_bits("stiff", &a.stiff, &want.stiff);
}

#[test]
fn runge_kutta_stage() {
    let mut rng = Lcg(6);
    let sys = system(&mut rng, false);
    let (prev, cur) = (rng.vec(N), rng.vec(N));
    let dt = 0.2;
    let m = Mat(sys.mass.as_ref().unwrap(), N);
    let mut want = sys.clone();
    for i in 1..=N {
        let mut s = 0.0;
        for j in 1..=N {
            s = s + (1.0 / dt) * m.at(i, j) * prev[j - 1] - Mat(&want.stiff, N).at(i, j) * cur[j - 1];
        }
        for j in 1..=N {
            set(&mut want.stiff, N, i, j, (1.0 / dt) * m.at(i, j));
        }
        want.force[i - 1] += s;
    }
    let mut got = sys.clone();
    let mut aux = Aux1st {
        cur: Some(&cur),
        ..Aux1st::default()
    };
    add_1st_order_time(&mut got, dt, &Scheme::RungeKutta, &[&prev], &mut aux).unwrap();
    assert_bits("stiff", &got.stiff, &want.stiff);
    assert_bits("force", &got.force, &want.force);
}

#[test]
fn fractional_step_stage() {
    let mut rng = Lcg(7);
    let sys = system(&mut rng, false);
    let prev = rng.vec(N);
    let dt = 0.3;
    let fs = FsCoeffs {
        mass_coeff: 0.2929,
        force_coeff: 0.7172,
    };
    let m = Mat(sys.mass.as_ref().unwrap(), N);
    let mut want = sys.clone();
    for i in 1..=N {
        let mut s = 0.0;
        for j in 1..=N {
            s = s + (1.0 / dt) * m.at(i, j) * prev[j - 1]
                - fs.force_coeff * Mat(&want.stiff, N).at(i, j) * prev[j - 1];
        }
        want.force[i - 1] += s;
        for j in 1..=N {
            let k = Mat(&want.stiff, N).at(i, j);
            set(&mut want.stiff, N, i, j, fs.mass_coeff * k + (1.0 / dt) * m.at(i, j));
        }
    }
    let mut got = sys.clone();
    add_1st_order_time(&mut got, dt, &Scheme::FractionalStep(fs), &[&prev], &mut Aux1st::default())
        .unwrap();
    assert_bits("stiff", &got.stiff, &want.stiff);
    assert_bits("force", &got.force, &want.force);
}

#[test]
fn adams_bashforth_then_moulton_share_element_state() {
    let mut rng = Lcg(8);
    let dt = 0.04;
    let zeta = 1.3;
    let mut state = AdamsState::default();
    // Reference state (Element % propertydata % values).
    let mut pd = [0.0; 2 * N];
    for step in 0..3 {
        let done_time = 4 + step; // PredCorrOrder = min(2, DoneTime/2) = 2
        let sys = system(&mut rng, false);
        let (p1, p2) = (rng.vec(N), rng.vec(N));
        let m = Mat(sys.mass.as_ref().unwrap(), N);

        let mut want = sys.clone();
        for i in 1..=N {
            let (mut sc, mut mc) = (0.0, 0.0);
            for j in 1..=N {
                sc = sc + Mat(&want.stiff, N).at(i, j) * p1[j - 1];
                mc = mc + (1.0 / dt) * m.at(i, j) * p1[j - 1];
            }
            let curr_res = -sc;
            let residual = pd[i - 1];
            let pre_force = pd[i - 1 + N];
            pd[i - 1] = curr_res;
            pd[i - 1 + N] = want.force[i - 1];
            let f = want.force[i - 1];
            want.force[i - 1] = f - sc + mc + 0.5 * zeta * (f - pre_force) + 0.5 * zeta * (curr_res - residual);
            for j in 1..=N {
                set(&mut want.stiff, N, i, j, (1.0 / dt) * m.at(i, j));
            }
        }
        let mut got = sys.clone();
        let mut aux = Aux1st {
            done_time,
            adams: Some(&mut state),
            ..Aux1st::default()
        };
        add_1st_order_time(&mut got, dt, &Scheme::AdamsBashforth { zeta, order: 2 }, &[&p1], &mut aux)
            .unwrap();
        assert_bits("AB stiff", &got.stiff, &want.stiff);
        assert_bits("AB force", &got.force, &want.force);

        // Corrector on a fresh local system, using the stored state.
        let sys2 = system(&mut rng, false);
        let m2 = Mat(sys2.mass.as_ref().unwrap(), N);
        let mut want2 = sys2.clone();
        for i in 1..=N {
            let (mut sc, mut mc) = (0.0, 0.0);
            for j in 1..=N {
                sc = sc + Mat(&want2.stiff, N).at(i, j) * p1[j - 1];
                mc = mc + (1.0 / dt) * m2.at(i, j) * p2[j - 1];
            }
            for j in 1..=N {
                set(&mut want2.stiff, N, i, j, (1.0 / dt) * m2.at(i, j));
            }
            let f = want2.force[i - 1];
            want2.force[i - 1] =
                0.5 * (f + pd[i - 1 + N]) + mc + 0.5 * (-sc + pd[i - 1]);
        }
        let mut got2 = sys2.clone();
        let mut aux = Aux1st {
            done_time,
            adams: Some(&mut state),
            ..Aux1st::default()
        };
        add_1st_order_time(&mut got2, dt, &Scheme::AdamsMoulton { order: 2 }, &[&p1, &p2], &mut aux)
            .unwrap();
        assert_bits("AM stiff", &got2.stiff, &want2.stiff);
        assert_bits("AM force", &got2.force, &want2.force);
    }
}

#[test]
fn adams_moulton_before_bashforth_is_refused() {
    let mut rng = Lcg(9);
    let mut sys = system(&mut rng, false);
    let (p1, p2) = (rng.vec(N), rng.vec(N));
    let mut state = AdamsState::default();
    let mut aux = Aux1st {
        done_time: 4,
        adams: Some(&mut state),
        ..Aux1st::default()
    };
    let e = add_1st_order_time(&mut sys, 0.1, &Scheme::AdamsMoulton { order: 2 }, &[&p1, &p2], &mut aux);
    assert!(matches!(e, Err(TimeError::AdamsMoultonBeforeBashforth)));
}

#[test]
fn lumping_keeps_total_mass() {
    let mut rng = Lcg(10);
    let mut m: Vec<f64> = rng.vec(N * N).into_iter().map(|x| x + 2.0).collect();
    let total: f64 = m.iter().sum();
    lump_matrix(&mut m, N);
    for i in 0..N {
        for j in 0..N {
            if i != j {
                assert_eq!(m[i * N + j], 0.0);
            }
        }
    }
    let lumped: f64 = (0..N).map(|i| m[i * N + i]).sum();
    assert!((lumped - total).abs() < 1e-12 * total.abs());
}

#[test]
fn lumped_flag_lumps_before_integrating() {
    let mut rng = Lcg(11);
    let sys = system(&mut rng, false);
    let prev = rng.vec(N);
    let mut a = sys.clone();
    let mut aux = Aux1st {
        lumped: true,
        ..Aux1st::default()
    };
    add_1st_order_time(&mut a, 0.1, &Scheme::NewmarkBeta(1.0), &[&prev], &mut aux).unwrap();
    let mut b = sys.clone();
    lump_matrix(b.mass.as_mut().unwrap(), N);
    add_1st_order_time(&mut b, 0.1, &Scheme::NewmarkBeta(1.0), &[&prev], &mut Aux1st::default())
        .unwrap();
    assert_bits("stiff", &a.stiff, &b.stiff);
    assert_bits("force", &a.force, &b.force);
}

#[test]
fn missing_mass_and_history_are_errors() {
    let mut sys = LocalSystem::zeros(2);
    let p = [0.0, 0.0];
    assert!(matches!(
        add_1st_order_time(&mut sys, 0.1, &Scheme::NewmarkBeta(1.0), &[&p], &mut Aux1st::default()),
        Err(TimeError::MissingMass(_))
    ));
    sys.mass = Some(vec![1.0; 4]);
    assert!(matches!(
        add_1st_order_time(&mut sys, 0.1, &Scheme::Bdf(2), &[], &mut Aux1st::default()),
        Err(TimeError::InvalidBdfOrder(0))
    ));
    assert!(matches!(
        add_1st_order_time(&mut sys, 0.1, &Scheme::Bdf(6), &[&p[..]; 6], &mut Aux1st { dts: &[0.1; 5], ..Aux1st::default() }),
        Err(TimeError::InvalidBdfOrder(6))
    ));
    assert!(matches!(
        add_1st_order_time(&mut sys, 0.1, &Scheme::RungeKutta, &[&p], &mut Aux1st::default()),
        Err(TimeError::MissingCurrent(_))
    ));
}

// ---- second order -------------------------------------------------------

fn params(case: u8) -> (f64, f64) {
    match case {
        0 => (-0.05, -1.0), // Bossak
        1 => (0.0, -1.0),   // Newmark average acceleration
        _ => (0.0, 0.75),   // generalized-alpha
    }
}

#[allow(clippy::too_many_arguments)]
fn ref_second(
    dt: f64,
    sys: &LocalSystem,
    x: &[f64],
    v: &[f64],
    a: &[f64],
    a2: &[f64],
    pa: f64,
    pb: f64,
) -> LocalSystem {
    let (alpha, beta, gamma, delta, a2b);
    if pb >= 0.0 {
        alpha = (2.0 * pb - 1.0) / (pb + 1.0);
        beta = 1.0 / (pb + 1.0).powi(2);
        gamma = (3.0 - pb) / (2.0 * (pb + 1.0));
        delta = pb / (pb + 1.0);
        a2b = a2;
    } else {
        alpha = pa;
        gamma = 0.5e0 - alpha;
        beta = (1.0 - alpha).powi(2) / 4.0;
        delta = 0.0;
        a2b = a;
    }
    let (m, d) = (Mat(sys.mass.as_ref().unwrap(), N), Mat(sys.damp.as_ref().unwrap(), N));
    let mut out = sys.clone();
    for i in 1..=N {
        let mut s = 0.0e0;
        for j in 1..=N {
            let (jj, mij, dij) = (j - 1, m.at(i, j), d.at(i, j));
            s = s + (1.0 - alpha) / (1.0 - delta) / (beta * dt.powi(2)) * mij * x[jj];
            s = s + (1.0 - alpha) / (1.0 - delta) / (beta * dt) * mij * v[jj];
            s = s + delta / (1.0 - delta) * mij * a[jj];
            s = s - ((1.0 - alpha) * (1.0 - 1.0 / (2.0 * beta)) + alpha) / (1.0 - delta) * mij * a2b[jj];
            s = s + gamma / (beta * dt) * dij * x[jj];
            s = s + (gamma / beta - 1.0) * dij * v[jj];
            s = s - ((1.0 - gamma) + gamma * (1.0 - 1.0 / (2.0 * beta))) * dt * dij * a2b[jj];
            let k = Mat(&out.stiff, N).at(i, j);
            set(
                &mut out.stiff,
                N,
                i,
                j,
                k + (1.0 - alpha) / (1.0 - delta) / (beta * dt.powi(2)) * mij + (gamma / (beta * dt)) * dij,
            );
        }
        out.force[i - 1] += s;
    }
    out
}

#[test]
fn second_order_bossak_newmark_and_generalized_alpha() {
    let mut rng = Lcg(12);
    for case in 0..3 {
        let (pa, pb) = params(case);
        let sys = system(&mut rng, true);
        let (x, v, a, a2) = (rng.vec(N), rng.vec(N), rng.vec(N), rng.vec(N));
        let dt = 0.015;
        let want = ref_second(dt, &sys, &x, &v, &a, &a2, pa, pb);
        let mut got = sys.clone();
        add_2nd_order_time(&mut got, dt, &x, &v, &a, &a2, pa, pb, false).unwrap();
        assert_bits("stiff", &got.stiff, &want.stiff);
        assert_bits("force", &got.force, &want.force);
    }
}

#[test]
fn update_2nd_order_matches_the_fortran_formulas() {
    let mut rng = Lcg(13);
    let dt = 0.02;
    for case in [0u8, 2] {
        let (pa, pb) = params(case);
        let x = rng.vec(N);
        let mut st = SecondOrderState::zeros(N);
        for c in 0..7 {
            st.cols[c] = rng.vec(N);
        }
        let before = st.clone();
        update_2nd_order(dt, &x, &mut st, pa, pb);
        let (alpha, beta, gamma, delta) = if pb >= 0.0 {
            ((2.0 * pb - 1.0) / (pb + 1.0), 1.0 / (pb + 1.0).powi(2), (3.0 - pb) / (2.0 * (pb + 1.0)), pb / (pb + 1.0))
        } else {
            (pa, (1.0 - pa).powi(2) / 4.0, 0.5e0 - pa, 0.0)
        };
        for i in 0..N {
            let dx = x[i] - before.cols[2][i];
            let (v, a, a2) = (before.cols[3][i], before.cols[4][i], before.cols[6][i]);
            if delta == 0.0 {
                let p1 = gamma / (beta * dt) * dx - (gamma - beta) / beta * v - dt * (gamma - 2.0 * beta) / (2.0 * beta) * a;
                let p2 = 1.0 / (beta * dt.powi(2)) * dx - 1.0 / (beta * dt) * v - (1.0 - 2.0 * beta) / (2.0 * beta) * a;
                assert_eq!(st.cols[0][i].to_bits(), p1.to_bits());
                assert_eq!(st.cols[1][i].to_bits(), p2.to_bits());
                assert_eq!(st.cols[5][i].to_bits(), before.cols[5][i].to_bits());
            } else {
                let p1 = gamma / (beta * dt) * dx - (gamma - beta) / beta * v - dt * (gamma - 2.0 * beta) / (2.0 * beta) * a2;
                let p2 = (1.0 - alpha) / (1.0 - delta) / (beta * dt.powi(2)) * dx
                    - (1.0 - alpha) / (1.0 - delta) / (beta * dt) * v
                    - delta / (1.0 - delta) * a
                    - (1.0 - alpha - 2.0 * beta) / (1.0 - delta) / (2.0 * beta) * a2;
                let p6 = (1.0 - delta) / (1.0 - alpha) * p2 + delta / (1.0 - alpha) * a - alpha / (1.0 - alpha) * a2;
                assert_eq!(st.cols[0][i].to_bits(), p1.to_bits());
                assert_eq!(st.cols[1][i].to_bits(), p2.to_bits());
                assert_eq!(st.cols[5][i].to_bits(), p6.to_bits());
            }
        }
    }
}

#[test]
fn newmark_average_acceleration_integrates_an_oscillator() {
    // m x'' + k x = 0, m = 1, k = 4: x = cos(2t). With alpha = 0 Newmark is
    // the trapezoidal rule; the energy error is tiny over a few periods.
    let (dt, steps) = (0.01, 300);
    let (mut x, mut v, mut a) = (1.0, 0.0, -4.0);
    let mut st = SecondOrderState::zeros(1);
    for _ in 0..steps {
        st.cols[2][0] = x;
        st.cols[3][0] = v;
        st.cols[4][0] = a;
        let mut sys = LocalSystem::zeros(1);
        sys.stiff[0] = 4.0;
        sys.mass = Some(vec![1.0]);
        add_2nd_order_time(&mut sys, dt, &[x], &[v], &[a], &[a], 0.0, -1.0, false).unwrap();
        let xn = sys.force[0] / sys.stiff[0];
        update_2nd_order(dt, &[xn], &mut st, 0.0, -1.0);
        // After the update: col0 = new velocity, col1 = new acceleration.
        x = xn;
        v = st.cols[0][0];
        a = st.cols[1][0];
    }
    let t = dt * steps as f64;
    assert!((x - (2.0 * t).cos()).abs() < 2e-3, "x = {x}");
}
