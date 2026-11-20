//! `BeamStiffnessMatrix` (SolidMechanicsUtils.F90) as ShellSolver calls it
//! for its two-node beam sections: Timoshenko beam with a condensed bubble
//! (`nd = 3`, `nb = 1`), combined with the shell's rotation DOFs
//! (`CombineWithShell`), the moment-free axis and the rotation to global
//! components.

use elmer_fe::invert_matrix;

use crate::geom::{cross, dot, norm, V3};
use crate::ShellError;

const DOFS: usize = 6;
/// `AEPS = 10 * EPSILON(1.0_dp)`.
const AEPS: f64 = 10.0 * f64::EPSILON;

/// Nodal material and load data of a beam element.
pub(crate) struct BeamInput {
    pub youngs: Vec<f64>,
    pub shear: Vec<f64>,
    pub form_factor: Vec<f64>,
    pub area: Vec<f64>,
    pub torsion: Vec<f64>,
    pub moment2: Vec<f64>,
    pub moment3: Vec<f64>,
    /// `Body Force 1..3` at the nodes.
    pub load: [Vec<f64>; 3],
    /// `Director` or `Principal Direction 2` of the material.
    pub director: Option<V3>,
    pub principal2: Option<V3>,
    /// `Moment-free Axis` (default 3).
    pub moment_free_axis: usize,
}

/// `TangentDirections( Normal, Tangent1, Tangent2 )`.
fn tangent_directions(n: V3) -> (V3, V3) {
    let (n1, n2, n3) = (n[0].abs(), n[1].abs(), n[2].abs());
    let mut t1 = if n1 <= n3 && n2 <= n3 {
        [0.0, -n[2], n[1]]
    } else {
        [-n[1], n[0], 0.0]
    };
    let s = norm(&t1);
    for x in &mut t1 {
        *x /= s;
    }
    let mut t2 = [
        n[1] * t1[2] - n[2] * t1[1],
        n[2] * t1[0] - n[0] * t1[2],
        n[0] * t1[1] - n[1] * t1[0],
    ];
    let s = norm(&t2);
    for x in &mut t2 {
        *x /= s;
    }
    (t1, t2)
}

/// `R^T (K R)` and `R^T F` for `n x n` row-major `K`, `R`.
fn rotate(k: &mut [f64], f: &mut [f64], r: &[f64], n: usize) {
    let mut kr = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for l in 0..n {
                s += k[i * n + l] * r[l * n + j];
            }
            kr[i * n + j] = s;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for l in 0..n {
                s += r[l * n + i] * kr[l * n + j];
            }
            k[i * n + j] = s;
        }
    }
    let fo: Vec<f64> = (0..n)
        .map(|i| {
            let mut s = 0.0;
            for l in 0..n {
                s += r[l * n + i] * f[l];
            }
            s
        })
        .collect();
    f[..n].copy_from_slice(&fo);
}

/// The beam's local system (`12 x 12`, row-major), its force and
/// `RHSForce`. `basis(u)` gives the 202 basis and its derivative at `u` on
/// the element `[0, L]`; `points` are the three-point rule; `prev` the
/// previous solution (large deflection).
pub(crate) fn beam_stiffness(
    xyz: [V3; 2],
    inp: &BeamInput,
    points: &[(f64, f64)],
    basis: &dyn Fn(f64, f64) -> Result<([f64; 2], [f64; 2], f64), ShellError>,
    large_deflection: bool,
    prev: &[f64],
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), ShellError> {
    let n = 2;
    let nd = 3;
    let nt = nd * DOFS;
    let mut stiff = vec![0.0f64; nt * nt];
    let mut force = vec![0.0f64; nt];
    let mut e1 = [
        xyz[1][0] - xyz[0][0],
        xyz[1][1] - xyz[0][1],
        xyz[1][2] - xyz[0][2],
    ];
    let l = norm(&e1);
    e1 = [1.0 / l * e1[0], 1.0 / l * e1[1], 1.0 / l * e1[2]];
    let (e2, e3) = if let Some(d) = inp.director {
        let nr = norm(&d);
        let e3 = [1.0 / nr * d[0], 1.0 / nr * d[1], 1.0 / nr * d[2]];
        if dot(&e1, &e3).abs() > 100.0 * AEPS {
            return Err(ShellError::Fatal(
                "BeamStiffnessMatrix: Director should be orthogonal to the beam axis".into(),
            ));
        }
        (cross(e3, e1), e3)
    } else {
        let (e2, e3) = if let Some(d) = inp.principal2 {
            let nr = norm(&d);
            let e2 = [1.0 / nr * d[0], 1.0 / nr * d[1], 1.0 / nr * d[2]];
            (e2, cross(e1, e2))
        } else {
            if (0..n).any(|i| (inp.moment3[i] - inp.moment2[i]).abs() > 5.0 * AEPS) {
                return Err(ShellError::Fatal(
                    "BeamStiffnessMatrix: The default principal directions need the same moments of area".into(),
                ));
            }
            tangent_directions(e1)
        };
        if dot(&e1, &e2).abs() > 100.0 * AEPS {
            return Err(ShellError::Fatal(
                "BeamStiffnessMatrix: Principal Direction 2 should be orthogonal to the beam axis"
                    .into(),
            ));
        }
        (e2, e3)
    };
    // SUM( Basis(1:n) * a(1:n) ), n = 2.
    let sum2 = |b: &[f64; 3], a: &[f64]| {
        let mut acc = 0.0;
        for i in 0..2 {
            acc += b[i] * a[i];
        }
        acc
    };
    for &(u, s) in points {
        let (b2, db2, det_j) = basis(u, l)?;
        // The bubble: Basis(3) = B1 B2, its derivative by the product rule.
        let b = [b2[0], b2[1], b2[0] * b2[1]];
        let db = [db2[0], db2[1], db2[0] * b2[1] + b2[0] * db2[1]];
        let f = [
            sum2(&b, &inp.load[0]),
            sum2(&b, &inp.load[1]),
            sum2(&b, &inp.load[2]),
        ];
        let e = sum2(&b, &inp.youngs);
        let g = sum2(&b, &inp.shear);
        let ff = sum2(&b, &inp.form_factor);
        let a = sum2(&b, &inp.area);
        let ed = [
            g * sum2(&b, &inp.torsion),
            e * sum2(&b, &inp.moment2),
            e * sum2(&b, &inp.moment3),
        ];
        let ga = ff * g * a;
        let ea = e * a;
        let w = s * det_j;
        let idx = |p: usize, q: usize, i: usize, j: usize| (p * DOFS + i) * nt + q * DOFS + j;
        for p in 0..nd {
            for q in 0..nd {
                stiff[idx(p, q, 0, 0)] += ea * db[q] * db[p] * w;
                stiff[idx(p, q, 1, 1)] += ga * db[q] * db[p] * w;
                stiff[idx(p, q, 2, 2)] += ga * db[q] * db[p] * w;
                if q >= n {
                    continue;
                }
                stiff[idx(p, q, 1, 5)] -= ga * b[q] * db[p] * w;
                stiff[idx(p, q, 2, 4)] += ga * b[q] * db[p] * w;
            }
            force[p * DOFS] += w * dot(&f, &e1) * b[p];
            force[p * DOFS + 1] += w * dot(&f, &e2) * b[p];
            force[p * DOFS + 2] += w * dot(&f, &e3) * b[p];
            if p >= n {
                continue;
            }
            for q in 0..nd {
                stiff[idx(p, q, 4, 2)] += ga * b[p] * db[q] * w;
                stiff[idx(p, q, 5, 1)] -= ga * b[p] * db[q] * w;
                if q >= n {
                    continue;
                }
                stiff[idx(p, q, 3, 3)] += ed[0] * db[q] * db[p] * w;
                stiff[idx(p, q, 4, 4)] =
                    stiff[idx(p, q, 4, 4)] + ed[1] * db[q] * db[p] * w + ga * b[p] * b[q] * w;
                stiff[idx(p, q, 5, 5)] =
                    stiff[idx(p, q, 5, 5)] + ed[2] * db[q] * db[p] * w + ga * b[p] * b[q] * w;
            }
        }
    }

    // BeamCondensate( 2, 1, 6, 3, Stiff, Force ): the bubble's three
    // displacement dofs.
    let m = n * DOFS;
    let bd = [m, m + 1, m + 2];
    let mut kbb = [0.0f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            // column-major for InvertMatrix
            kbb[i + 3 * j] = stiff[bd[i] * nt + bd[j]];
        }
    }
    invert_matrix(&mut kbb, 3)
        .ok_or_else(|| ShellError::Fatal("BeamCondensate: singular bubble block".into()))?;
    let kbbm = |i: usize, j: usize| kbb[i + 3 * j];
    let fb: [f64; 3] = std::array::from_fn(|i| force[bd[i]]);
    let kfb: [f64; 3] = std::array::from_fn(|i| {
        let mut s = 0.0;
        for j in 0..3 {
            s += kbbm(i, j) * fb[j];
        }
        s
    });
    let mut kbl = vec![[0.0f64; 3]; m];
    for c in 0..m {
        for i in 0..3 {
            let mut s = 0.0;
            for j in 0..3 {
                s += kbbm(i, j) * stiff[bd[j] * nt + c];
            }
            kbl[c][i] = s;
        }
    }
    let mut k = vec![0.0f64; m * m];
    let mut fr = vec![0.0f64; m];
    for r in 0..m {
        let mut s = 0.0;
        for j in 0..3 {
            s += stiff[r * nt + bd[j]] * kfb[j];
        }
        fr[r] = force[r] - s;
        for c in 0..m {
            let mut s = 0.0;
            for j in 0..3 {
                s += stiff[r * nt + bd[j]] * kbl[c][j];
            }
            k[r * m + c] = stiff[r * nt + c] - s;
        }
    }

    // CombineWithShell: the rotation DOFs of the shell (no offsets).
    let mut r = vec![0.0f64; m * m];
    for i in 0..n {
        let i0 = i * DOFS;
        r[i0 * m + i0] = 1.0;
        r[(i0 + 1) * m + i0 + 1] = 1.0;
        r[(i0 + 2) * m + i0 + 2] = 1.0;
        r[(i0 + 3) * m + i0 + 4] = 1.0;
        r[(i0 + 4) * m + i0 + 3] = -1.0;
        r[(i0 + 5) * m + i0 + 5] = 1.0;
    }
    rotate(&mut k, &mut fr, &r, m);
    for p in 0..n {
        let i = p * DOFS + 2 + inp.moment_free_axis;
        for c in 0..m {
            k[i * m + c] = 0.0;
            k[c * m + i] = 0.0;
        }
        k[i * m + i] = 1.0e1 * AEPS;
        fr[i] = 0.0;
    }
    // The rotation to global components.
    let mut r = vec![0.0f64; m * m];
    let blk = [e1, e2, e3];
    for i in 0..n {
        let i0 = i * DOFS;
        for a in 0..3 {
            for c in 0..3 {
                r[(i0 + a) * m + i0 + c] = blk[a][c];
                r[(i0 + 3 + a) * m + i0 + 3 + c] = blk[a][c];
            }
        }
    }
    rotate(&mut k, &mut fr, &r, m);
    let rhs = fr.clone();
    if large_deflection {
        for i in 0..m {
            let mut s = 0.0;
            for j in 0..m {
                s += k[i * m + j] * prev[j];
            }
            fr[i] -= s;
        }
    }
    Ok((k, fr, rhs))
}
