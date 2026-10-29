//! Moving the target boundary onto this boundary (`OverlayIntefaceMeshes`,
//! MortarUtils.F90).
//!
//! Products are written as the Fortran's `MATMUL` sums, term by term from
//! zero, so the transformed coordinates match the pristine binary bit for bit.

use elmer_core::{Logger, Params, Value};

use crate::interface::InterfaceMesh;
use crate::MortarError;

const CALLER: &str = "OverlayInterfaceMeshes";

type M4 = [[f64; 4]; 4];

fn identity() -> M4 {
    let mut m = [[0.0; 4]; 4];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    m
}

/// `MATMUL(a, b)` for 4x4 matrices.
fn matmul(a: &M4, b: &M4) -> M4 {
    let mut c = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut s = 0.0;
            for k in 0..4 {
                s += a[i][k] * b[k][j];
            }
            c[i][j] = s;
        }
    }
    c
}

/// `MATMUL(m, x)`: matrix times column vector.
fn mat_vec(m: &M4, x: [f64; 4]) -> [f64; 4] {
    let mut y = [0.0; 4];
    for (i, yi) in y.iter_mut().enumerate() {
        let mut s = 0.0;
        for k in 0..4 {
            s += m[i][k] * x[k];
        }
        *yi = s;
    }
    y
}

/// `MATMUL(x, m)`: row vector times matrix.
fn vec_mat(x: [f64; 4], m: &M4) -> [f64; 4] {
    let mut y = [0.0; 4];
    for (j, yj) in y.iter_mut().enumerate() {
        let mut s = 0.0;
        for k in 0..4 {
            s += x[k] * m[k][j];
        }
        *yj = s;
    }
    y
}

/// A constant real array as `ListGetConstRealArray` returns it:
/// `(rows, cols, row-major data)`.
fn real_array<'a>(p: &'a dyn Params, key: &str) -> Option<(usize, usize, &'a [f64])> {
    match p.get(key)? {
        Value::Real { rows, cols, data } => Some((*rows, *cols, data.as_slice())),
        _ => None,
    }
}

/// Column 1 of an array keyword (`Parray(1:n,1)`), as many rows as it has.
fn first_column(rows: usize, cols: usize, data: &[f64]) -> Vec<f64> {
    (0..rows).map(|r| data[r * cols]).collect()
}

fn bbox(b: &InterfaceMesh) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for x in &b.coords {
        for d in 0..3 {
            lo[d] = lo[d].min(x[d]);
            hi[d] = hi[d].max(x[d]);
        }
    }
    (lo, hi)
}

fn es3(v: [f64; 3]) -> String {
    v.iter().map(|&x| es15_6(x)).collect()
}

/// `ES15.6`.
fn es15_6(x: f64) -> String {
    let s = format!("{:.6e}", x.abs());
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    let sign = if x.is_sign_negative() { "-" } else { "" };
    format!(
        "{:>15}",
        format!("{sign}{m}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
    )
}

/// `OverlayIntefaceMeshes`: transform the nodes of `b2` (the target side)
/// so that they overlay `b1`, from the boundary condition `bc`.
pub fn overlay_interface_meshes(
    b1: &InterfaceMesh,
    b2: &mut InterfaceMesh,
    bc: &dyn Params,
    log: &dyn Logger,
) -> Result<(), MortarError> {
    let (x1_min, x1_max) = bbox(b1);
    log.info(
        CALLER,
        &format!("Minimum values for this periodic BC:  {}", es3(x1_min)),
        8,
    );
    log.info(
        CALLER,
        &format!("Maximum values for this periodic BC:  {}", es3(x1_max)),
        8,
    );
    let (x2_min, x2_max) = bbox(b2);
    log.info(
        CALLER,
        &format!("Minimum values for target periodic BC:{}", es3(x2_min)),
        8,
    );
    log.info(
        CALLER,
        &format!("Maximum values for target periodic BC:{}", es3(x2_max)),
        8,
    );

    let trf = if let Some((rows, cols, data)) = real_array(bc, "periodic bc matrix") {
        // TrfMatrix(i,j) = Parray(j,i); a smaller array would leave the rest
        // of the Fortran's matrix undefined.
        if rows != 4 || cols != 4 {
            return Err(MortarError::Unsupported(format!(
                "Periodic BC Matrix of size {rows}x{cols} (only 4x4 is defined)"
            )));
        }
        let mut t = [[0.0; 4]; 4];
        for (i, row) in t.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v = data[j * cols + i];
            }
        }
        t
    } else {
        let mut trs = identity();
        let mut rot = identity();
        let mut scl = identity();

        let mut angles = [0.0; 3];
        let got_rotate = match real_array(bc, "periodic bc rotate") {
            Some((rows, cols, data)) => {
                let c = first_column(rows, cols, data);
                if c.len() < 3 {
                    return Err(MortarError::Fatal(format!(
                        "{CALLER}: Periodic BC Rotate needs 3 angles, got {}",
                        c.len()
                    )));
                }
                angles.copy_from_slice(&c[..3]);
                true
            }
            None => {
                if bc.logical("periodic bc rotate automatic").unwrap_or(false) {
                    return Err(MortarError::Unsupported(
                        "Periodic BC Rotate Automatic".into(),
                    ));
                }
                false
            }
        };

        let (x2r_min, x2r_max) = if got_rotate {
            log.info(CALLER, &format!("Rotating target with: {}", es3(angles)), 8);
            for (i, &a) in angles.iter().enumerate() {
                let alpha = a * std::f64::consts::PI / 180.0;
                if alpha.abs() < f64::MIN_POSITIVE {
                    continue;
                }
                let mut t = identity();
                let (p, q) = match i {
                    0 => (1, 2),
                    1 => (0, 2),
                    _ => (0, 1),
                };
                t[p][p] = alpha.cos();
                t[p][q] = -alpha.sin();
                t[q][p] = alpha.sin();
                t[q][q] = alpha.cos();
                rot = matmul(&rot, &t);
            }
            for x in &mut b2.coords {
                let y = mat_vec(&rot, [x[0], x[1], x[2], 1.0]);
                *x = [y[0], y[1], y[2]];
            }
            let (lo, hi) = bbox(b2);
            log.info(
                CALLER,
                &format!("Minimum values for rotated target:{}", es3(lo)),
                8,
            );
            log.info(
                CALLER,
                &format!("Maximum values for rotated target:{}", es3(hi)),
                8,
            );
            (lo, hi)
        } else {
            (x2_min, x2_max)
        };

        if let Some((rows, cols, data)) = real_array(bc, "periodic bc scale") {
            for (i, v) in first_column(rows, cols, data)
                .into_iter()
                .enumerate()
                .take(4)
            {
                scl[i][i] = v;
            }
        } else {
            let mut s1 = 0.0;
            let mut s2 = 0.0;
            for d in 0..3 {
                let (a, b) = (x1_max[d] - x1_min[d], x2r_max[d] - x2r_min[d]);
                s1 += a * a;
                s2 += b * b;
            }
            let s = if s2 > f64::EPSILON {
                (s1 / s2).sqrt()
            } else {
                1.0
            };
            log.info(CALLER, &format!("Scaling with: {}", es3([s; 3])), 4);
            for (d, row) in scl.iter_mut().enumerate().take(3) {
                row[d] = s;
            }
        }

        if let Some((rows, cols, data)) = real_array(bc, "periodic bc translate") {
            for (i, v) in first_column(rows, cols, data)
                .into_iter()
                .enumerate()
                .take(4)
            {
                trs[3][i] = v;
            }
        } else {
            for i in 0..3 {
                trs[3][i] = x1_min[i] - scl[i][i] * x2r_min[i];
            }
        }
        log.info(
            CALLER,
            &format!("Translation: {}", es3([trs[3][0], trs[3][1], trs[3][2]])),
            4,
        );
        matmul(&scl, &trs)
    };

    for x in &mut b2.coords {
        let y = vec_mat([x[0], x[1], x[2], 1.0], &trf);
        *x = [y[0] / y[3], y[1] / y[3], y[2] / y[3]];
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::es15_6;

    #[test]
    fn es15_6_matches_fortran() {
        assert_eq!(es15_6(0.0), "   0.000000E+00");
        assert_eq!(es15_6(-1.5e-3), "  -1.500000E-03");
        assert_eq!(es15_6(2.0), "   2.000000E+00");
    }
}
