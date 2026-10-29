//! The constitutive laws `ElasticSolve.F90` evaluates through `Constitutive.F90`
//! (`IsotropicLinearStress`, `AnisotropicLinearStress`) and the plane condensation
//! of `Stress.F90`. The linear crate keeps the same routines private to its
//! small-strain kernel, so they are restated here for the finite strain one.

use crate::tensor::{T3, ZERO};

pub type Voigt6 = [f64; 6];

/// Elmer's Voigt order as index pairs: 11, 22, 33, 12, 23, 13.
pub(crate) const I1: [usize; 6] = [0, 1, 2, 0, 1, 0];
pub(crate) const I2: [usize; 6] = [0, 1, 2, 1, 2, 2];

/// `Stress.F90:1531` `CondensePlaneElasticityMatrix`.
pub(crate) fn condense_plane(mut c: [[f64; 6]; 6], plane_stress: bool) -> [[f64; 6]; 6] {
    if plane_stress {
        c[0][0] -= c[0][2] * c[2][0] / c[2][2];
        c[0][1] -= c[0][2] * c[1][2] / c[2][2];
        c[1][0] -= c[1][2] * c[0][2] / c[2][2];
        c[1][1] -= c[1][2] * c[2][1] / c[2][2];
    } else {
        c[3][0] = c[2][0];
        c[3][1] = c[2][1];
        c[3][2] = c[2][3];
    }
    c[2][2] = c[3][3];
    c[0][2] = 0.0;
    c[2][0] = 0.0;
    c[1][2] = 0.0;
    c[2][1] = 0.0;
    c
}

/// `Stress.F90:2123` `RotateElasticityMatrix`: `C'(a,b,c,d) = T(a,p) T(b,q)
/// T(c,r) T(d,s) C(p,q,r,s)` on the rank-4 tensor behind `c`. In two dimensions
/// only the in-plane 2x2 block of Voigt pairs is rotated.
pub(crate) fn rotate_c(c: &[[f64; 6]; 6], t: &T3, dim: usize) -> [[f64; 6]; 6] {
    let nv = if dim == 2 { 2 } else { 6 };
    let mut full = [[[[0.0; 3]; 3]; 3]; 3];
    for i in 0..nv {
        let (p, q) = (I1[i], I2[i]);
        for j in 0..nv {
            let (r, s) = (I1[j], I2[j]);
            let v = c[i][j];
            full[p][q][r][s] = v;
            full[p][q][s][r] = v;
            full[q][p][r][s] = v;
            full[q][p][s][r] = v;
        }
    }
    let mut rot = [[[[0.0; 3]; 3]; 3]; 3];
    for a in 0..dim {
        for b in 0..dim {
            for cc in 0..dim {
                for d in 0..dim {
                    let mut s = 0.0;
                    for p in 0..dim {
                        for q in 0..dim {
                            for r in 0..dim {
                                for u in 0..dim {
                                    s += t[a][p] * t[b][q] * t[cc][r] * t[d][u] * full[p][q][r][u];
                                }
                            }
                        }
                    }
                    rot[a][b][cc][d] = s;
                }
            }
        }
    }
    let mut out = *c;
    for i in 0..nv {
        for j in 0..nv {
            out[i][j] = rot[I1[i]][I2[i]][I1[j]][I2[j]];
        }
    }
    out
}

/// `Stress.F90:1908` `Vector62Tensor`.
pub(crate) fn voigt_to_tensor(v: &Voigt6, cdim: usize, axi: bool) -> T3 {
    let n = match (cdim, axi) {
        (2, true) => 4,
        (2, false) => 3,
        (3, _) => 6,
        _ => 1,
    };
    // The 2D packings are 11, 22, 12 and, axisymmetric, 11, 22, 33, 12.
    let (i1, i2): ([usize; 6], [usize; 6]) = match (cdim, axi) {
        (2, true) => ([0, 1, 2, 0, 0, 0], [0, 1, 2, 1, 0, 0]),
        (2, false) => ([0, 1, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0]),
        _ => (I1, I2),
    };
    let mut x = ZERO;
    for i in 0..n {
        x[i1[i]][i2[i]] = v[i];
        x[i2[i]][i1[i]] = v[i];
    }
    x
}

/// The linear constitutive law at one integration point.
pub(crate) enum PointLaw {
    Iso { lame1: f64, lame2: f64 },
    /// Condensed (2D) or full (3D) `C`.
    Aniso { g: Box<[[f64; 6]; 6]> },
}

impl PointLaw {
    /// The stress of a strain tensor for a state of dimension `dim` (3 for
    /// axisymmetry). `transposed` evaluates with `C^T`, which is what the
    /// stiffness batch does (`MatPropsT`).
    pub(crate) fn stress(&self, dim: usize, transposed: bool, e: &T3) -> T3 {
        match self {
            PointLaw::Iso { lame1, lame2 } => {
                let tr: f64 = (0..dim).map(|i| e[i][i]).sum();
                let mut s = ZERO;
                for i in 0..3 {
                    for j in 0..3 {
                        s[i][j] = 2.0 * lame2 * e[i][j];
                    }
                }
                for (i, row) in s.iter_mut().enumerate().take(dim) {
                    row[i] += lame1 * tr;
                }
                s
            }
            PointLaw::Aniso { g } => {
                let mut s = ZERO;
                let (n, v): (usize, Voigt6) = if dim == 2 {
                    (3, [e[0][0], e[1][1], 2.0 * e[0][1], 0.0, 0.0, 0.0])
                } else {
                    (6, [e[0][0], e[1][1], e[2][2], 2.0 * e[0][1], 2.0 * e[1][2], 2.0 * e[0][2]])
                };
                // The reduced plane packing has the shear in slot 3.
                let (p, q): (&[usize], &[usize]) = if n == 3 { (&[0, 1, 0], &[0, 1, 1]) } else { (&I1, &I2) };
                for i in 0..n {
                    let mut csum = 0.0;
                    for j in 0..n {
                        csum += if transposed { g[j][i] } else { g[i][j] } * v[j];
                    }
                    s[p[i]][q[i]] = csum;
                    s[q[i]][p[i]] = csum;
                }
                s
            }
        }
    }
}
