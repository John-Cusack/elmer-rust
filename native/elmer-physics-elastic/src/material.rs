//! Constitutive laws of ElasticSolve's small-strain path (`Constitutive.F90`
//! `IsotropicLinearStress` / `AnisotropicLinearStress`) and the plane
//! condensation of `Stress.F90`.

/// Elmer's Voigt order as index pairs (1-based in Fortran, 0-based here):
/// 11, 22, 33, 12, 23, 13.
pub(crate) const I1: [usize; 6] = [0, 1, 2, 0, 1, 0];
pub(crate) const I2: [usize; 6] = [0, 1, 2, 1, 2, 2];

pub type Tensor3 = [[f64; 3]; 3];
pub type Voigt6 = [f64; 6];

pub(crate) const ZERO3: Tensor3 = [[0.0; 3]; 3];

/// How the elasticity moduli are given, per element node.
#[derive(Clone, Debug)]
pub enum ElasticLaw {
    /// `Youngs Modulus` and `Poisson Ratio`, both scalar.
    Isotropic { young: Vec<f64>, poisson: Vec<f64> },
    /// `Youngs Modulus` as a tensor: the full 6x6 stiffness `C` in Elmer's
    /// Voigt order (xx, yy, zz, xy, yz, xz) with engineering shear, per node.
    /// A scalar, a vector or a matrix in the sif all arrive here already
    /// expanded as `InputTensor` does (`Stress.F90:1990`): scalar and vector
    /// forms fill the diagonal only.
    Anisotropic {
        c: Vec<[[f64; 6]; 6]>,
        /// Rows are the material-frame unit vectors
        /// (`Rotate Elasticity Tensor`, `Material Coordinates Unit Vector i`).
        rotation: Option<Tensor3>,
    },
}

/// Fill a diagonal `C` the way `InputTensor` does for a scalar or a vector
/// `Youngs Modulus`: `values[i]` on `C(i,i)`, at most six of them.
pub fn diagonal_c(values: &[f64]) -> [[f64; 6]; 6] {
    let mut c = [[0.0; 6]; 6];
    for (i, v) in values.iter().take(6).enumerate() {
        c[i][i] = *v;
    }
    c
}

/// `Stress.F90:1531` `CondensePlaneElasticityMatrix`. Returns the condensed
/// matrix (reduced plane packing in `c[0..3][0..3]`) and `EzzC`.
pub(crate) fn condense_plane(mut c: [[f64; 6]; 6], plane_stress: bool) -> ([[f64; 6]; 6], [f64; 3]) {
    let mut ezz = [0.0; 3];
    if plane_stress {
        ezz[0] = -c[2][0] / c[2][2];
        ezz[1] = -c[2][1] / c[2][2];
        ezz[2] = -c[2][3] / c[2][2];

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
    (c, ezz)
}

/// `Stress.F90:2123` `RotateElasticityMatrix`: `C'(a,b,c,d) = T(a,p) T(b,q)
/// T(c,r) T(d,s) C(p,q,r,s)` on the rank-4 tensor behind `c`.
///
/// In two dimensions Elmer rotates only the in-plane 2x2 block (it loops
/// `i,j = 1,2` over the Voigt pairs 11 and 22), leaving every other entry of
/// `c` as given. This reproduces that, including the dropped shear coupling.
pub(crate) fn rotate_c(c: &[[f64; 6]; 6], t: &Tensor3, dim: usize) -> [[f64; 6]; 6] {
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

/// The constitutive law evaluated at one integration point.
pub(crate) enum PointLaw {
    Iso { lame1: f64, lame2: f64 },
    /// Condensed (2D) or full (3D) `C`.
    Aniso { g: Box<[[f64; 6]; 6]> },
}

impl PointLaw {
    /// `sigma(strain)` for a state of stress of dimension `dim` (3 for
    /// axisymmetry). `transposed` evaluates with `C^T`, which is what
    /// ElasticSolve's stiffness batch does (`MatPropsT`).
    pub(crate) fn stress(&self, dim: usize, transposed: bool, e: &Tensor3) -> Tensor3 {
        match self {
            PointLaw::Iso { lame1, lame2 } => {
                let tr: f64 = (0..dim).map(|i| e[i][i]).sum();
                let mut s = ZERO3;
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
                let mut s = ZERO3;
                let (n, v): (usize, Voigt6) = if dim == 2 {
                    (3, [e[0][0], e[1][1], 2.0 * e[0][1], 0.0, 0.0, 0.0])
                } else {
                    (
                        6,
                        [
                            e[0][0],
                            e[1][1],
                            e[2][2],
                            2.0 * e[0][1],
                            2.0 * e[1][2],
                            2.0 * e[0][2],
                        ],
                    )
                };
                // The reduced plane packing has the shear in slot 3.
                let (p, q): (&[usize], &[usize]) = if n == 3 {
                    (&[0, 1, 0], &[0, 1, 1])
                } else {
                    (&I1, &I2)
                };
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
