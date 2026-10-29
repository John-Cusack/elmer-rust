//! Node distributions on the unit segment (MeshBasics.F90:8973-9146).

use crate::MeshLevelsError;

/// `x**n` for an integer `n` as gfortran evaluates it: libgcc's
/// `__powidf2`, square-and-multiply from the low bit.
pub fn powi(x: f64, m: i32) -> f64 {
    let mut n = m.unsigned_abs();
    let mut x = x;
    let mut y = if n % 2 == 1 { x } else { 1.0 };
    loop {
        n >>= 1;
        if n == 0 {
            break;
        }
        x *= x;
        if n % 2 == 1 {
            y *= x;
        }
    }
    if m < 0 {
        1.0 / y
    } else {
        y
    }
}

/// `GeometricUnitDivision`: `n+1` points on `[0,1]`, the element lengths in
/// ratio `q` between the last and the first, or symmetric about the middle
/// for negative `q`.
pub fn geometric(n: usize, q: f64) -> Result<Vec<f64>, MeshLevelsError> {
    if n < 1 {
        return Err(MeshLevelsError::Invalid(format!(
            "cannot create division for {n} element"
        )));
    }
    let ni = n as i32;
    let nf = n as f64;
    let mut w = vec![0.0; n + 1];
    if (q.abs() - 1.0).abs() < 1.0e-6 || (q < 0.0 && n <= 2) || n == 1 {
        for (i, wi) in w.iter_mut().enumerate() {
            *wi = i as f64 / (1.0 * nf);
        }
    } else if q > 0.0 {
        let r = q.powf(1.0 / (ni - 1) as f64);
        let h1 = (1.0 - r) / (1.0 - powi(r, ni));
        for (i, wi) in w.iter_mut().enumerate().take(n).skip(1) {
            *wi = h1 * (1.0 - powi(r, i as i32)) / (1.0 - r);
        }
        w[n] = 1.0;
    } else {
        let q = -q;
        let (r, h1);
        if ni % 2 == 0 {
            r = q.powf(1.0 / (ni / 2 - 1) as f64);
            h1 = 0.5 * (1.0 - r) / (1.0 - powi(r, ni / 2));
        } else {
            r = q.powf(1.0 / ((ni - 1) / 2) as f64);
            h1 = 0.5 / ((1.0 - powi(r, (ni + 1) / 2)) / (1.0 - r) - 0.5 * powi(r, (ni - 1) / 2));
        }
        for i in 1..=ni {
            w[i as usize] = if i <= ni / 2 {
                h1 * (1.0 - powi(r, i)) / (1.0 - r)
            } else {
                1.0 - h1 * (1.0 - powi(r, ni - i)) / (1.0 - r)
            };
        }
        w[n] = 1.0;
    }
    Ok(w)
}

/// `FunctionUnitDivision`: element lengths proportional to the mesh density
/// `h(x)` at element midpoints, by Gauss-Seidel sweeps from an even start.
pub fn function(
    n: usize,
    h: &mut dyn FnMut(f64) -> Result<f64, String>,
) -> Result<Vec<f64>, MeshLevelsError> {
    if n < 1 {
        return Err(MeshLevelsError::Invalid(format!(
            "cannot create division for {n} element"
        )));
    }
    let mut w: Vec<f64> = (0..=n).map(|i| i as f64 / (1.0 * n as f64)).collect();
    if n == 1 {
        return Ok(w);
    }
    let mut hv = vec![0.0; n + 1];
    for _ in 0..10000 {
        let mut minhn = f64::MAX;
        let wold = w.clone();
        for i in 1..=n {
            let xn = (w[i] + w[i - 1]) / 2.0;
            minhn = minhn.min(w[i] - w[i - 1]);
            hv[i] = h(xn).map_err(MeshLevelsError::Invalid)?;
            if hv[i] < f64::EPSILON {
                return Err(MeshLevelsError::Invalid(
                    "given value for h(i) was negative".into(),
                ));
            }
        }
        for i in 1..n {
            w[i] = (w[i - 1] * hv[i + 1] + w[i + 1] * hv[i]) / (hv[i] + hv[i + 1]);
        }
        for i in (1..n).rev() {
            w[i] = (w[i - 1] * hv[i + 1] + w[i + 1] * hv[i]) / (hv[i] + hv[i + 1]);
        }
        let err = w
            .iter()
            .zip(&wold)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max)
            / minhn;
        if err < 1.0e-6 {
            break;
        }
    }
    Ok(w)
}
