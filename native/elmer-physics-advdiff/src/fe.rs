//! What the kernels need from the element layer beyond `elmer-core`'s
//! [`ElementBasis`] and [`Quadrature`]: the residual-free bubble basis of
//! `ElementInfo(..., Bubbles=.TRUE.)` (ElemInfo.F90:571, `EvalBubbleBasis`),
//! its quadrature rule, and a hook for second derivatives.

use elmer_core::{BasisEval, ElementBasis, ElementType, Quadrature};

use crate::error::AdvDiffError;

/// Second derivatives of the nodal basis in physical space, `ddBasisddx` of
/// Elmer's `ElementInfo(..., SecondDerivatives=.TRUE.)`.
///
/// `elmer-core`'s [`ElementBasis`] returns first derivatives only, and the SUPG
/// residual `SU(p)` needs `d2 N_p / dx_i dx_j`. This trait is the interim shape
/// of that request: when `elmer-fe` provides second derivatives it should
/// implement this for its basis (or `BasisEval` should grow the field), and
/// the SUPG path of [`crate::advdiff_bulk_local`] is then usable from the
/// solver. Until then the solver refuses `Stabilize = True`.
pub trait BasisHessian {
    /// `out[p][i][j]` is `d2 N_p / dx_i dx_j` for the `xyz.len()` nodal basis
    /// functions at reference point `uvw`.
    fn eval_hessian(
        &self,
        etype: ElementType,
        xyz: &[[f64; 3]],
        uvw: [f64; 3],
        out: &mut Vec<[[f64; 3]; 3]>,
    ) -> Result<(), String>;
}

/// The element layer a kernel call evaluates with.
#[derive(Copy, Clone)]
pub struct Fe<'a> {
    pub basis: &'a dyn ElementBasis,
    pub quad: &'a dyn Quadrature,
    /// Needed only for SUPG stabilisation.
    pub hessian: Option<&'a dyn BasisHessian>,
}

impl<'a> Fe<'a> {
    pub fn new(basis: &'a dyn ElementBasis, quad: &'a dyn Quadrature) -> Self {
        Self {
            basis,
            quad,
            hessian: None,
        }
    }
}

/// `Element % TYPE % GaussPoints2`, the second quadrature rule of each
/// element type in `elements.def` (`Gauss Points <rule> <rule2> [<rule0>]`),
/// which `GaussPoints( element, GaussPoints2 )` selects for the bubble basis.
pub fn bubble_gauss_points(etype: ElementType) -> Option<usize> {
    Some(match etype.0 {
        202 => 4,
        203 => 5,
        204 => 6,
        303 => 11,
        306 => 17,
        310 => 36,
        404 => 16,
        408 | 409 => 25,
        412 | 416 => 36,
        504 => 11,
        510 => 64,
        605 | 706 => 64,
        613 | 614 | 715 | 718 => 125,
        808 => 64,
        820 | 827 => 125,
        _ => return None,
    })
}

/// Evaluates `ElementInfo` with `Bubbles=.TRUE.` into `out`: the `n` nodal
/// functions followed by `n` bubble functions `N_i * b`, where `b` is the
/// product of the linear nodal functions of the element family (`L1 L2` on
/// lines, `L1 L2 L3` on triangles, `L1 L3` on quadrilaterals, `L1 L2 L3 L4`
/// on tetrahedra, `L1 L7` on bricks).
///
/// As in Elmer, a non-linear element takes its Jacobian determinant from the
/// linear element of the family (`detJ` is overwritten by the inner
/// `ElementInfo` call).
pub fn eval_with_bubbles(
    basis: &dyn ElementBasis,
    etype: ElementType,
    xyz: &[[f64; 3]],
    uvw: [f64; 3],
    dim: usize,
    out: &mut BasisEval,
) -> Result<(), AdvDiffError> {
    let n = xyz.len();
    basis.eval(etype, xyz, uvw, out).map_err(AdvDiffError::Basis)?;
    let (lin_code, picks): (u16, &[usize]) = match etype.0 / 100 {
        2 => (202, &[0, 1]),
        3 => (303, &[0, 1, 2]),
        4 => (404, &[0, 2]),
        5 => (504, &[0, 1, 2, 3]),
        8 => (808, &[0, 6]),
        _ => {
            return Err(AdvDiffError::Basis(format!(
                "bubbles for element {} are not implemented",
                etype.0
            )))
        }
    };
    let lin_n = (lin_code % 100) as usize;
    let mut lin = BasisEval::default();
    if etype.0 == lin_code {
        lin.det_j = out.det_j;
        lin.basis = out.basis.clone();
        lin.dbasis_dx = out.dbasis_dx.clone();
    } else {
        basis
            .eval(ElementType(lin_code), &xyz[..lin_n], uvw, &mut lin)
            .map_err(AdvDiffError::Basis)?;
        out.det_j = lin.det_j;
    }

    let value: f64 = picks.iter().map(|&k| lin.basis[k]).product();
    out.basis.truncate(n);
    out.dbasis_dx.truncate(n);
    for i in 0..n {
        let bi = out.basis[i];
        let mut d = [0.0; 3];
        for (j, dj) in d.iter_mut().enumerate().take(dim) {
            *dj = out.dbasis_dx[i][j] * value;
            for &a in picks {
                let others: f64 = picks
                    .iter()
                    .filter(|&&b| b != a)
                    .map(|&b| lin.basis[b])
                    .product();
                *dj += bi * lin.dbasis_dx[a][j] * others;
            }
        }
        out.basis.push(bi * value);
        out.dbasis_dx.push(d);
    }
    Ok(())
}

/// `CondensateP( n, nb, K, F, F1 )` (MatrixAssembly.F90:937): eliminates the
/// trailing `nb` bubble rows and columns of the `(n+nb)`-square row-major `k`
/// by static condensation, updating the leading `n x n` block of `k` and the
/// first `n` entries of `f` and `f1`. The leading block is left in place in the
/// `(n+nb)`-strided layout.
pub fn condensate(n: usize, nb: usize, k: &mut [f64], f: &mut [f64], f1: Option<&mut [f64]>) -> Result<(), AdvDiffError> {
    if nb == 0 {
        return Ok(());
    }
    let m = n + nb;
    let mut kbb = vec![0.0; nb * nb];
    for i in 0..nb {
        for j in 0..nb {
            kbb[i * nb + j] = k[(n + i) * m + n + j];
        }
    }
    let inv = invert(&kbb, nb)?;
    // Klb * Kbb^-1, n x nb.
    let mut a = vec![0.0; n * nb];
    for i in 0..n {
        for j in 0..nb {
            a[i * nb + j] = (0..nb).map(|l| k[i * m + n + l] * inv[l * nb + j]).sum();
        }
    }
    let fb: Vec<f64> = f[n..m].to_vec();
    for i in 0..n {
        f[i] -= (0..nb).map(|j| a[i * nb + j] * fb[j]).sum::<f64>();
    }
    if let Some(f1) = f1 {
        let fb: Vec<f64> = f1[n..m].to_vec();
        for i in 0..n {
            f1[i] -= (0..nb).map(|j| a[i * nb + j] * fb[j]).sum::<f64>();
        }
    }
    for i in 0..n {
        for j in 0..n {
            let s: f64 = (0..nb).map(|l| a[i * nb + l] * k[(n + l) * m + j]).sum();
            k[i * m + j] -= s;
        }
    }
    Ok(())
}

/// Dense inverse by Gauss-Jordan with partial pivoting.
fn invert(a: &[f64], n: usize) -> Result<Vec<f64>, AdvDiffError> {
    let mut m = a.to_vec();
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    for c in 0..n {
        let p = (c..n)
            .max_by(|&x, &y| m[x * n + c].abs().total_cmp(&m[y * n + c].abs()))
            .unwrap_or(c);
        if m[p * n + c] == 0.0 {
            return Err(AdvDiffError::Fatal(
                "Condensate: the bubble block is singular".into(),
            ));
        }
        if p != c {
            for j in 0..n {
                m.swap(c * n + j, p * n + j);
                inv.swap(c * n + j, p * n + j);
            }
        }
        let d = m[c * n + c];
        for j in 0..n {
            m[c * n + j] /= d;
            inv[c * n + j] /= d;
        }
        for i in 0..n {
            if i != c {
                let f = m[i * n + c];
                if f != 0.0 {
                    for j in 0..n {
                        m[i * n + j] -= f * m[c * n + j];
                        inv[i * n + j] -= f * inv[c * n + j];
                    }
                }
            }
        }
    }
    Ok(inv)
}
