//! `_CRS` routines against scalar transliterations of TimeIntegrate.F90 on a
//! random sparse pattern (bitwise).

#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]

mod support;

use elmer_core::{CsrMatrix, NullLogger};
use elmer_timeint::local::Aux1st;
use elmer_timeint::scheme::FsCoeffs;
use elmer_timeint::second::{time_2nd_order_crs, SecondOrderState};
use elmer_timeint::{
    add_1st_order_time_crs, bdf_crs, fractional_step_crs, newmark_beta_crs, runge_kutta_crs,
    vbdf_crs, CrsTimeData, Scheme, TimeError,
};
use support::{assert_bits, Lcg};

const N: usize = 8;

struct Fix {
    a: CsrMatrix,
    mass: Vec<f64>,
    damp: Vec<f64>,
    ml: Vec<f64>,
    res: Vec<f64>,
    force: Vec<f64>,
}

fn fixture(rng: &mut Lcg) -> Fix {
    let mut row_ptr = vec![0usize];
    let mut col_idx = Vec::new();
    let mut diag = Vec::new();
    for i in 0..N {
        for j in 0..N {
            if i == j || rng.next_f64() > 0.2 {
                if i == j {
                    diag.push(col_idx.len());
                }
                col_idx.push(j as u32);
            }
        }
        row_ptr.push(col_idx.len());
    }
    let nnz = col_idx.len();
    Fix {
        a: CsrMatrix {
            nrows: N,
            ncols: N,
            row_ptr,
            col_idx,
            values: rng.vec(nnz),
            diag: Some(diag),
        },
        mass: rng.vec(nnz).into_iter().map(|x| x + 1.4).collect(),
        damp: rng.vec(nnz),
        ml: rng.vec(N).into_iter().map(|x| x + 1.4).collect(),
        res: rng.vec(N),
        force: rng.vec(N),
    }
}

impl Fix {
    fn rows(&self, i: usize) -> std::ops::Range<usize> {
        self.a.row_ptr[i]..self.a.row_ptr[i + 1]
    }
    fn col(&self, j: usize) -> usize {
        self.a.col_idx[j] as usize
    }
    fn data(&self) -> CrsTimeData<'_> {
        CrsTimeData {
            mass: Some(&self.mass),
            ..Default::default()
        }
    }
}

#[test]
fn newmark_beta_all_four_branches() {
    let mut rng = Lcg(21);
    let dt = 0.07;
    for &beta in &[1.0, 0.5, 0.25] {
        for (lumped, resid) in [(false, false), (true, false), (false, true), (true, true)] {
            let f = fixture(&mut rng);
            let prev = rng.vec(N);
            // Reference.
            let (mut st, mut fo) = (f.a.values.clone(), f.force.clone());
            if resid {
                if lumped {
                    let old = st.clone();
                    let _ = old;
                    for i in 0..N {
                        let ui = prev[i];
                        fo[i] = beta * fo[i] - (1.0 - beta) * f.res[i] + (1.0 / dt) * f.ml[i] * ui;
                    }
                    for v in st.iter_mut() {
                        *v = beta * *v;
                    }
                    for i in 0..N {
                        let j = f.a.diag.as_ref().unwrap()[i];
                        st[j] = st[j] + (1.0 / dt) * f.ml[i];
                    }
                } else {
                    for i in 0..N {
                        let mut mu = 0.0;
                        for j in f.rows(i) {
                            mu = mu + f.mass[j] * prev[f.col(j)];
                        }
                        fo[i] = beta * fo[i] - (1.0 - beta) * f.res[i] + (1.0 / dt) * mu;
                    }
                    for j in 0..st.len() {
                        st[j] = beta * st[j] + (1.0 / dt) * f.mass[j];
                    }
                }
            } else if lumped {
                for i in 0..N {
                    let mut su = 0.0;
                    for j in f.rows(i) {
                        su = su + st[j] * prev[f.col(j)];
                    }
                    let mu = f.ml[i] * prev[i];
                    fo[i] = fo[i] - (1.0 - beta) * su + (1.0 / dt) * mu;
                }
                for v in st.iter_mut() {
                    *v = beta * *v;
                }
                for i in 0..N {
                    let j = f.a.diag.as_ref().unwrap()[i];
                    st[j] = st[j] + (1.0 / dt) * f.ml[i];
                }
            } else {
                for i in 0..N {
                    let (mut su, mut mu) = (0.0, 0.0);
                    for j in f.rows(i) {
                        let uj = prev[f.col(j)];
                        su = su + st[j] * uj;
                        mu = mu + f.mass[j] * uj;
                    }
                    fo[i] = fo[i] - (1.0 - beta) * su + (1.0 / dt) * mu;
                }
                for j in 0..st.len() {
                    st[j] = beta * st[j] + (1.0 / dt) * f.mass[j];
                }
            }
            // Implementation.
            let mut a = f.a.clone();
            let mut force = f.force.clone();
            let d = CrsTimeData {
                mass: Some(&f.mass),
                mass_lumped: lumped.then_some(&f.ml[..]),
                bulk_residual: resid.then_some(&f.res[..]),
                ..Default::default()
            };
            newmark_beta_crs(dt, &mut a, &d, &mut force, &prev, beta).unwrap();
            assert_bits(&format!("values b={beta} l={lumped} r={resid}"), &a.values, &st);
            assert_bits("force", &force, &fo);
        }
    }
}

#[test]
fn runge_kutta_branches() {
    let mut rng = Lcg(22);
    let dt = 0.11;
    for (lumped, resid) in [(false, false), (true, false), (false, true), (true, true)] {
        let f = fixture(&mut rng);
        let (prev, cur) = (rng.vec(N), rng.vec(N));
        let (mut st, mut fo) = (f.a.values.clone(), f.force.clone());
        match (lumped, resid) {
            (false, false) => {
                for i in 0..N {
                    let (mut su, mut mu) = (0.0, 0.0);
                    for j in f.rows(i) {
                        su = su + st[j] * cur[f.col(j)];
                        mu = mu + f.mass[j] * prev[f.col(j)];
                    }
                    fo[i] = fo[i] - su + (1.0 / dt) * mu;
                }
                for j in 0..st.len() {
                    st[j] = (1.0 / dt) * f.mass[j];
                }
            }
            (true, false) => {
                for i in 0..N {
                    let mut su = 0.0;
                    for j in f.rows(i) {
                        su = su + st[j] * cur[f.col(j)];
                    }
                    let mu = f.ml[i] * prev[i];
                    fo[i] = fo[i] - su + (1.0 / dt) * mu;
                }
                st.iter_mut().for_each(|v| *v = 0.0);
                for i in 0..N {
                    let j = f.a.diag.as_ref().unwrap()[i];
                    st[j] = st[j] + (1.0 / dt) * f.ml[i];
                }
            }
            (false, true) => {
                for i in 0..N {
                    let mut mu = 0.0;
                    for j in f.rows(i) {
                        mu = mu + f.mass[j] * prev[f.col(j)];
                    }
                    fo[i] = f.res[i] + (1.0 / dt) * mu;
                }
                for j in 0..st.len() {
                    st[j] = (1.0 / dt) * f.mass[j];
                }
            }
            (true, true) => {
                for i in 0..N {
                    fo[i] = f.res[i] + (1.0 / dt) * f.ml[i] * prev[i];
                }
                st.iter_mut().for_each(|v| *v = 0.0);
                for i in 0..N {
                    let j = f.a.diag.as_ref().unwrap()[i];
                    st[j] = st[j] + (1.0 / dt) * f.ml[i];
                }
            }
        }
        let mut a = f.a.clone();
        let mut force = f.force.clone();
        let d = CrsTimeData {
            mass: Some(&f.mass),
            mass_lumped: lumped.then_some(&f.ml[..]),
            bulk_residual: resid.then_some(&f.res[..]),
            ..Default::default()
        };
        runge_kutta_crs(dt, &mut a, &d, &mut force, &prev, &cur).unwrap();
        assert_bits(&format!("values l={lumped} r={resid}"), &a.values, &st);
        assert_bits("force", &force, &fo);
    }
}

fn bdf_a(order: usize) -> [f64; 6] {
    let mut a = [0.0; 6];
    match order {
        1 => { a[0] = 1.0; a[1] = -1.0; }
        2 => { a[0] = 1.5; a[1] = -2.0; a[2] = 0.5; }
        3 => { a[0] = 11.0 / 6.0; a[1] = -3.0; a[2] = 3.0 / 2.0; a[3] = -1.0 / 3.0; }
        4 => { a[0] = 25.0 / 12.0; a[1] = -4.0; a[2] = 3.0; a[3] = -4.0 / 3.0; a[4] = 1.0 / 4.0; }
        _ => { a[0] = 137.0 / 60.0; a[1] = -5.0; a[2] = 5.0; a[3] = -10.0 / 3.0; a[4] = 5.0 / 4.0; a[5] = -1.0 / 5.0; }
    }
    a
}

#[test]
fn bdf_crs_orders_one_to_five_lumped_and_consistent() {
    let mut rng = Lcg(23);
    let dt = 0.05;
    for order in 1..=5usize {
        for lumped in [false, true] {
            let f = fixture(&mut rng);
            let hist: Vec<Vec<f64>> = (0..5).map(|_| rng.vec(N)).collect();
            let prev: Vec<&[f64]> = hist.iter().map(|v| &v[..]).take(order).collect();
            let mut a_w = bdf_a(order);
            for w in a_w.iter_mut() {
                *w = *w / dt;
            }
            let (mut st, mut fo) = (f.a.values.clone(), f.force.clone());
            if lumped {
                for i in 0..N {
                    for k in 1..=order {
                        let ui = prev[k - 1][i];
                        fo[i] = fo[i] - f.ml[i] * a_w[k] * ui;
                    }
                    let j = f.a.diag.as_ref().unwrap()[i];
                    st[j] = st[j] + a_w[0] * f.ml[i];
                }
            } else {
                for i in 0..N {
                    let mut mu = 0.0;
                    for j in f.rows(i) {
                        for k in 1..=order {
                            let uj = prev[k - 1][f.col(j)];
                            mu = mu - f.mass[j] * a_w[k] * uj;
                        }
                        st[j] = st[j] + a_w[0] * f.mass[j];
                    }
                    fo[i] = fo[i] + mu;
                }
            }
            let mut a = f.a.clone();
            let mut force = f.force.clone();
            let d = CrsTimeData {
                mass: Some(&f.mass),
                mass_lumped: Some(&f.ml),
                lumped,
                ..Default::default()
            };
            bdf_crs(dt, &mut a, &d, &mut force, &prev, order).unwrap();
            assert_bits(&format!("values o={order} l={lumped}"), &a.values, &st);
            assert_bits("force", &force, &fo);
        }
    }
    let f = fixture(&mut rng);
    let mut a = f.a.clone();
    let mut force = f.force.clone();
    assert!(matches!(
        bdf_crs(dt, &mut a, &f.data(), &mut force, &[&[0.0; N][..]; 6], 6),
        Err(TimeError::InvalidBdfOrder(6))
    ));
    bdf_crs(dt, &mut a, &f.data(), &mut force, &[], 0).unwrap();
    assert_bits("order 0 is a no-op", &force, &f.force);
}

#[test]
fn vbdf_crs_matches_bdf_when_steps_are_equal_and_the_formula_otherwise() {
    let mut rng = Lcg(24);
    let f = fixture(&mut rng);
    let hist: Vec<Vec<f64>> = (0..3).map(|_| rng.vec(N)).collect();
    let prev: Vec<&[f64]> = hist.iter().map(|v| &v[..]).collect();
    // Equal steps: coefficients agree with constant-step BDF to rounding.
    let dt = 0.1;
    let (mut a1, mut f1) = (f.a.clone(), f.force.clone());
    let (mut a2, mut f2) = (f.a.clone(), f.force.clone());
    bdf_crs(dt, &mut a1, &f.data(), &mut f1, &prev, 3).unwrap();
    vbdf_crs(&[dt, dt, dt], &mut a2, &f.data(), &mut f2, &prev, 3, &NullLogger).unwrap();
    for (x, y) in a1.values.iter().zip(&a2.values) {
        assert!((x - y).abs() < 1e-12, "{x} {y}");
    }
    for (x, y) in f1.iter().zip(&f2) {
        assert!((x - y).abs() < 1e-12, "{x} {y}");
    }
    // Unequal steps, order 2, against the literal formulas.
    let dts = [0.05, 0.03];
    let a1c = 1.0 / dts[0] + 1.0 / (dts[0] + dts[1]);
    let a2c = -1.0 / dts[0] - (1.0 + dts[0] / dts[1]) / (dts[0] + dts[1]);
    let a3c = (dts[0] / dts[1]) / (dts[0] + dts[1]);
    let (mut st, mut fo) = (f.a.values.clone(), f.force.clone());
    for i in 0..N {
        let mut mu = 0.0;
        for j in f.rows(i) {
            for (k, ak) in [a2c, a3c].iter().enumerate() {
                mu = mu - f.mass[j] * ak * prev[k][f.col(j)];
            }
        }
        fo[i] = fo[i] + mu;
    }
    for j in 0..st.len() {
        st[j] = st[j] + a1c * f.mass[j];
    }
    let (mut a, mut force) = (f.a.clone(), f.force.clone());
    vbdf_crs(&dts, &mut a, &f.data(), &mut force, &prev, 2, &NullLogger).unwrap();
    assert_bits("values", &a.values, &st);
    assert_bits("force", &force, &fo);
}

#[test]
fn fractional_step_crs_both_branches() {
    let mut rng = Lcg(25);
    let dt = 0.09;
    let fs = FsCoeffs { mass_coeff: 0.3, force_coeff: 0.8 };
    for lumped in [false, true] {
        let f = fixture(&mut rng);
        let prev = rng.vec(N);
        let (mut st, mut fo) = (f.a.values.clone(), f.force.clone());
        if lumped {
            for i in 0..N {
                let mut su = 0.0;
                for j in f.rows(i) {
                    su = su + st[j] * prev[f.col(j)];
                }
                let j = f.a.diag.as_ref().unwrap()[i];
                let ui = prev[f.col(j)];
                fo[i] = fo[i] - fs.force_coeff * su + (1.0 / dt) * f.ml[i] * ui;
                st[j] = fs.mass_coeff * st[j] + (1.0 / dt) * f.ml[i];
            }
        } else {
            for i in 0..N {
                let (mut su, mut mu) = (0.0, 0.0);
                for j in f.rows(i) {
                    let uj = prev[f.col(j)];
                    su = su + st[j] * uj;
                    mu = mu + f.mass[j] * uj;
                }
                fo[i] = fo[i] - fs.force_coeff * su + (1.0 / dt) * mu;
            }
            for j in 0..st.len() {
                st[j] = fs.mass_coeff * st[j] + (1.0 / dt) * f.mass[j];
            }
        }
        let (mut a, mut force) = (f.a.clone(), f.force.clone());
        let d = CrsTimeData {
            mass: Some(&f.mass),
            mass_lumped: lumped.then_some(&f.ml[..]),
            ..Default::default()
        };
        fractional_step_crs(dt, &mut a, &d, &mut force, &prev, fs).unwrap();
        assert_bits("values", &a.values, &st);
        assert_bits("force", &force, &fo);
    }
}

#[test]
fn time_2nd_order_crs_matches_formulas() {
    let mut rng = Lcg(26);
    let dt: f64 = 0.02;
    for (pa, pb) in [(-0.05f64, -1.0f64), (0.0, 0.75)] {
        let f = fixture(&mut rng);
        let mut state = SecondOrderState::zeros(N);
        for c in 0..7 {
            state.cols[c] = rng.vec(N);
        }
        let (alpha, beta, gamma, delta) = if pb >= 0.0 {
            ((2.0 * pb - 1.0) / (pb + 1.0), 1.0 / (pb + 1.0).powi(2), (3.0 - pb) / (2.0 * (pb + 1.0)), pb / (pb + 1.0))
        } else {
            (pa, (1.0 - pa).powi(2) / 4.0, 0.5e0 - pa, 0.0)
        };
        let a2 = if pb >= 0.0 { &state.cols[6] } else { &state.cols[4] };
        let (x, v, a) = (&state.cols[2], &state.cols[3], &state.cols[4]);
        let c1 = (1.0 - alpha) / (1.0 - delta) / (beta * dt.powi(2));
        let c2 = (1.0 - alpha) / (1.0 - delta) / (beta * dt);
        let c3 = delta / (1.0 - delta);
        let c4 = -((1.0 - alpha) * (1.0 - 1.0 / (2.0 * beta)) + alpha) / (1.0 - delta);
        let c5 = gamma / (beta * dt);
        let c6 = gamma / beta - 1.0;
        let c7 = -((1.0 - gamma) + gamma * (1.0 - 1.0 / (2.0 * beta))) * dt;
        let (mut st, mut fo) = (f.a.values.clone(), f.force.clone());
        for i in 0..N {
            let mut s = 0.0;
            for j in f.rows(i) {
                let k = f.col(j);
                s = s + f.mass[j] * (c1 * x[k] + c2 * v[k] + c3 * a[k] + c4 * a2[k]);
            }
            fo[i] = fo[i] + s;
        }
        for j in 0..st.len() {
            st[j] = st[j] + c1 * f.mass[j];
        }
        for i in 0..N {
            let mut s = 0.0;
            for j in f.rows(i) {
                let k = f.col(j);
                s = s + f.damp[j] * (c5 * x[k] + c6 * v[k] + c7 * a2[k]);
            }
            fo[i] = fo[i] + s;
        }
        for j in 0..st.len() {
            st[j] = st[j] + c5 * f.damp[j];
        }
        let (mut am, mut force) = (f.a.clone(), f.force.clone());
        let d = CrsTimeData { mass: Some(&f.mass), damp: Some(&f.damp), ..Default::default() };
        time_2nd_order_crs(dt, &mut am, &d, &mut force, pa, pb, &state).unwrap();
        assert_bits("values", &am.values, &st);
        assert_bits("force", &force, &fo);
    }
}

#[test]
fn dispatcher_picks_bdf_or_vbdf_and_refuses_adams() {
    let mut rng = Lcg(27);
    let f = fixture(&mut rng);
    let (p1, p2) = (rng.vec(N), rng.vec(N));
    let dt = 0.1;
    let run = |scheme: Scheme, dts: &[f64]| {
        let (mut a, mut force) = (f.a.clone(), f.force.clone());
        let aux = Aux1st { dts, ..Aux1st::default() };
        let r = add_1st_order_time_crs(&mut a, &f.data(), &mut force, dt, &scheme, &[&p1, &p2], &aux);
        (r, a, force)
    };
    let (r, a, force) = run(Scheme::Bdf(2), &[dt]);
    r.unwrap();
    let (mut a_ref, mut f_ref) = (f.a.clone(), f.force.clone());
    bdf_crs(dt, &mut a_ref, &f.data(), &mut f_ref, &[&p1, &p2], 2).unwrap();
    assert_bits("bdf values", &a.values, &a_ref.values);
    assert_bits("bdf force", &force, &f_ref);

    let (r, a, force) = run(Scheme::Bdf(2), &[0.2]);
    r.unwrap();
    let (mut a_ref, mut f_ref) = (f.a.clone(), f.force.clone());
    vbdf_crs(&[dt, 0.2], &mut a_ref, &f.data(), &mut f_ref, &[&p1, &p2], 2, &NullLogger).unwrap();
    assert_bits("vbdf values", &a.values, &a_ref.values);
    assert_bits("vbdf force", &force, &f_ref);

    let (r, _, _) = run(Scheme::AdamsBashforth { zeta: 1.0, order: 2 }, &[]);
    assert!(matches!(r, Err(TimeError::NotImplemented(_))));
}

#[test]
fn diag_is_found_by_search_when_absent() {
    let mut rng = Lcg(28);
    let mut f = fixture(&mut rng);
    let prev = rng.vec(N);
    let (mut a1, mut f1) = (f.a.clone(), f.force.clone());
    let d = CrsTimeData { mass: Some(&f.mass), mass_lumped: Some(&f.ml), ..Default::default() };
    newmark_beta_crs(0.1, &mut a1, &d, &mut f1, &prev, 1.0).unwrap();
    f.a.diag = None;
    let (mut a2, mut f2) = (f.a.clone(), f.force.clone());
    newmark_beta_crs(0.1, &mut a2, &d, &mut f2, &prev, 1.0).unwrap();
    assert_bits("values", &a1.values, &a2.values);
}
