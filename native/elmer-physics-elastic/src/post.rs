//! Stress and strain postprocessing of `ElasticSolve.F90` `ComputeStressAndStrain`
//! (`:5759`): the element contribution to the L2 nodal projection of the
//! `Strain` and `Stress` fields, the projection solve, von Mises, and the
//! principal values. Small strain by default; with `Large Deflection` the
//! Green-Lagrange strain and the Cauchy stress pushed forward from the
//! second Piola-Kirchhoff one, for the linear laws and the neo-Hookean
//! material (`Constitutive.F90` `NeoHookeanStress`).
//!
//! Output layout follows `SymTensorOutputComponents` (`Stress.F90:1623`): the
//! components `(11, 22, 33, 12)` in two dimensions and `(11, 22, 33, 12, 23, 13)`
//! in three, interleaved per node.

use elmer_core::{BasisEval, ElementBasis, ElementType, Quadrature};

use crate::bulk::{CoordSystem, IpProps};
use crate::error::StressError;
use crate::material::{condense_plane, rotate_c, ElasticLaw, PointLaw, Tensor3, ZERO3};

/// `SymTensorOutputComponents`.
pub fn output_components(cdim: usize) -> usize {
    if cdim <= 2 {
        4
    } else {
        6
    }
}

/// `SYMTENSOR_IND`: the output slot (0-based) of the tensor entry `(i, j)`.
const SLOT: [[usize; 3]; 3] = [[0, 3, 5], [3, 1, 4], [5, 4, 2]];

/// `OutputVector2Tensor`.
pub fn vector_to_tensor(v: &[f64], ncomp: usize) -> Tensor3 {
    let mut t = ZERO3;
    for i in 0..3 {
        for j in 0..3 {
            let k = SLOT[i][j];
            if k < ncomp {
                t[i][j] = v[k];
            }
        }
    }
    t
}

/// What `ComputeStressAndStrain` reads for one element, at its nodes.
#[derive(Clone, Debug)]
pub struct StressPostInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    pub cdim: usize,
    pub coord: CoordSystem,
    pub plane_stress: bool,
    pub law: &'a ElasticLaw,
    /// Nodal displacement, `cdim` components used.
    pub disp: &'a [[f64; 3]],
    /// Temperature minus `Reference Temperature`, per node.
    pub temperature: Option<&'a [f64]>,
    /// Diagonal of `Heat Expansion Coefficient`, per node.
    pub heat_expansion: Option<&'a [[f64; 3]]>,
    pub ip: Option<&'a IpProps>,
    pub want_strain: bool,
    pub want_stress: bool,
    pub np: Option<usize>,
    /// ElasticSolve's `Large Deflection`: `DefG = I + Grad`.
    pub large_deflection: bool,
    /// ElasticSolve's `Neo-Hookean Material` (isotropic `law` only; no
    /// thermal strain, no integration-point moduli).
    pub neo_hookean: bool,
    /// ElasticSolve's `Incompressible`: the pressure unknown at the nodes.
    /// The moduli become `Lame1 = 0`, `Lame2 = E/3` (no plane stress) and the
    /// pressure is taken off the normal stresses.
    pub pressure: Option<&'a [f64]>,
}

/// Element mass matrix (`n x n`, row-major) and right-hand sides
/// (`n * ncomp`, node-major) of the projection.
#[derive(Clone, Debug, Default)]
pub struct StressPostLocal {
    pub mass: Vec<f64>,
    pub strain: Vec<f64>,
    pub stress: Vec<f64>,
}

fn add_tensor(force: &mut [f64], bas: &[f64], ncomp: usize, w: f64, t: &Tensor3) {
    for (p, b) in bas.iter().enumerate() {
        for i in 0..3 {
            for j in i..3 {
                let k = SLOT[i][j];
                if k < ncomp {
                    force[ncomp * p + k] += w * t[i][j] * b;
                }
            }
        }
    }
}

/// `LUDecomp` (LinearAlgebra.F90): LU of the leading `n x n` block with
/// column pivoting by rows, the diagonal replaced by its reciprocal. `None`
/// on a zero pivot. Index loops as in the Fortran.
#[allow(clippy::needless_range_loop)]
fn lu_decomp(a: &mut Tensor3, n: usize) -> Option<[usize; 3]> {
    let mut pivot = [0usize; 3];
    for i in 0..n {
        let mut j = i;
        for k in i + 1..n {
            if a[i][k].abs() > a[i][j].abs() {
                j = k;
            }
        }
        if a[i][j].abs() == 0.0 {
            return None;
        }
        pivot[i] = j;
        if j != i {
            for row in a.iter_mut().take(i + 1) {
                row.swap(i, j);
            }
        }
        for k in i + 1..n {
            a[i][k] /= a[i][i];
        }
        for k in i + 1..n {
            if j != i {
                a[k].swap(i, j);
            }
            for l in i + 1..n {
                a[k][l] -= a[k][i] * a[i][l];
            }
        }
    }
    pivot[n - 1] = n - 1;
    for i in 0..n {
        if a[i][i].abs() == 0.0 {
            return None;
        }
        a[i][i] = 1.0 / a[i][i];
    }
    Some(pivot)
}

/// `InvertMatrix( A, n )` (LinearAlgebra.F90): in place, through
/// [`lu_decomp`], in its operation order.
fn invert_matrix(a: &mut Tensor3, n: usize) -> Result<(), StressError> {
    let pivot = lu_decomp(a, n)
        .ok_or_else(|| StressError::Input("InvertMatrix: inversion needs successful LU-decomposition".into()))?;
    for i in (0..n.saturating_sub(1)).rev() {
        for j in (i + 1..n).rev() {
            let mut s = -a[i][j];
            for k in i + 1..j {
                s -= a[i][k] * a[k][j];
            }
            a[i][j] = s;
        }
    }
    for i in (0..n.saturating_sub(1)).rev() {
        for j in (i + 1..n).rev() {
            let mut s = 0.0;
            for k in i + 1..=j {
                s -= a[j][k] * a[k][i];
            }
            a[j][i] = a[i][i] * s;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in i.max(j)..n {
                if k != i {
                    s += a[i][k] * a[k][j];
                } else {
                    s += a[k][j];
                }
            }
            a[i][j] = s;
        }
    }
    for i in (0..n).rev() {
        if pivot[i] != i {
            a.swap(i, pivot[i]);
        }
    }
    Ok(())
}

/// Fortran's inlined `MATMUL` of two 3x3 matrices: each entry summed over
/// `k` in order from zero.
fn matmul3(a: &Tensor3, b: &Tensor3) -> Tensor3 {
    let mut c = ZERO3;
    for i in 0..3 {
        for j in 0..3 {
            let mut s = 0.0;
            for k in 0..3 {
                s += a[i][k] * b[k][j];
            }
            c[i][j] = s;
        }
    }
    c
}

fn transpose3(a: &Tensor3) -> Tensor3 {
    let mut t = ZERO3;
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}

/// Element contribution to the L2 projection of strain and stress.
pub fn stress_strain_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &StressPostInput<'_>,
) -> Result<StressPostLocal, StressError> {
    let points = quad.gauss_points(inp.etype, inp.np);
    stress_strain_local_p(basis, &points, inp.etype.n_nodes(), inp)
}

/// [`stress_strain_local`] on a p element (`ComputeStress` with `n` nodes
/// and `nd` dofs): `inp.disp` has one entry per dof and the displacement
/// gradient, the mass matrix and the right-hand sides run over the `nd`
/// basis functions; the material and the temperature are interpolated with
/// the nodal ones. `points` is the element's rule.
pub fn stress_strain_local_p(
    basis: &dyn ElementBasis,
    points: &[elmer_core::QuadPoint],
    nd: usize,
    inp: &StressPostInput<'_>,
) -> Result<StressPostLocal, StressError> {
    let n = inp.etype.n_nodes();
    let cdim = inp.cdim;
    let axi = inp.coord == CoordSystem::AxiSymmetric;
    if cdim != 2 && cdim != 3 {
        return Err(StressError::Unsupported(format!("coordinate system dimension {cdim}")));
    }
    if axi && cdim != 2 {
        return Err(StressError::Unsupported("axisymmetric elasticity needs a 2D mesh".into()));
    }
    if inp.xyz.len() != n || inp.disp.len() != nd || nd < n {
        return Err(StressError::Input(format!("element {} needs {n} nodes and {nd} displacements", inp.etype.0)));
    }
    let dim = if axi { 3 } else { cdim };
    let plane_stress = inp.plane_stress && !axi && cdim == 2 && inp.pressure.is_none();
    let ncomp = output_components(cdim);
    let isotropic = matches!(inp.law, ElasticLaw::Isotropic { .. });
    if inp.neo_hookean && !isotropic {
        return Err(StressError::Unsupported("the neo-Hookean material with an anisotropic Youngs Modulus".into()));
    }
    // `DefG`, `DetDefG` and the push-forward: the general path of
    // ComputeStressAndStrain. Small strain without the neo-Hookean material
    // keeps the masked form below, which is the same arithmetic with
    // `DefG = I`.
    let general = inp.large_deflection || inp.neo_hookean;

    // The identity of the state of stress: `ReducedIdentity`.
    let mut ident = [0.0; 3];
    for v in ident.iter_mut().take(cdim) {
        *v = 1.0;
    }
    if axi || !plane_stress && isotropic || cdim == 3 {
        ident[2] = 1.0;
    }

    let (lame1_n, lame2_n): (Vec<f64>, Vec<f64>) = match inp.law {
        ElasticLaw::Isotropic { young, poisson } => (0..n)
            .map(|k| {
                let (e, nu) = (young[k], poisson[k]);
                let l1 = if plane_stress {
                    e * nu / (1.0 - nu * nu)
                } else {
                    e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
                };
                if inp.pressure.is_some() {
                    (0.0, e / 3.0)
                } else {
                    (l1, e / (2.0 * (1.0 + nu)))
                }
            })
            .unzip(),
        ElasticLaw::Anisotropic { .. } => (Vec::new(), Vec::new()),
    };

    let need_heat = !inp.neo_hookean
        && inp.temperature.is_some_and(|t| t.iter().any(|x| *x != 0.0))
        && (inp.heat_expansion.is_some() || inp.ip.is_some_and(|p| p.beta.is_some()));

    if let Some(p) = inp.ip {
        let np = points.len();
        let lens = [p.young.as_ref().map(Vec::len), p.poisson.as_ref().map(Vec::len), p.beta.as_ref().map(Vec::len)];
        if lens.iter().flatten().any(|&l| l != np) {
            return Err(StressError::Input(format!("integration point data must have {np} entries")));
        }
    }

    let mut out = StressPostLocal {
        mass: vec![0.0; nd * nd],
        strain: if inp.want_strain { vec![0.0; nd * ncomp] } else { Vec::new() },
        stress: if inp.want_stress { vec![0.0; nd * ncomp] } else { Vec::new() },
    };
    let mut be = BasisEval::default();

    for (ipt, qp) in points.iter().enumerate() {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be).map_err(StressError::Basis)?;
        let bas = &be.basis[..n];
        let full = &be.basis[..nd];
        let db = &be.dbasis_dx;
        let mut w = be.det_j * qp.weight;
        let mut r = 0.0;
        if axi {
            r = (0..n).map(|k| bas[k] * inp.xyz[k][0]).sum();
            w *= r;
        }

        // Interpolated moduli, in the order `ComputeStressAndStrain` takes them.
        let mut nu_nodal = 0.0;
        let mut ezz_c = [0.0; 3];
        let law = match inp.law {
            ElasticLaw::Isotropic { poisson, .. } => {
                nu_nodal = (0..n).map(|k| poisson[k] * bas[k]).sum();
                let mut lame1: f64 = (0..n).map(|k| lame1_n[k] * bas[k]).sum();
                let mut lame2: f64 = (0..n).map(|k| lame2_n[k] * bas[k]).sum();
                // The at-IP moduli are not read for the neo-Hookean material.
                if let Some(p) = inp.ip.filter(|p| !inp.neo_hookean && (p.young.is_some() || p.poisson.is_some())) {
                    let ElasticLaw::Isotropic { young, .. } = inp.law else { unreachable!() };
                    let e = p.young.as_ref().map_or_else(|| (0..n).map(|k| young[k] * bas[k]).sum(), |v| v[ipt]);
                    let nu: f64 = p.poisson.as_ref().map_or(nu_nodal, |v| v[ipt]);
                    lame1 = if plane_stress {
                        e * nu / (1.0 - nu * nu)
                    } else {
                        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
                    };
                    lame2 = e / (2.0 * (1.0 + nu));
                }
                PointLaw::Iso { lame1, lame2 }
            }
            ElasticLaw::Anisotropic { c, rotation } => {
                let mut g = [[0.0; 6]; 6];
                for (i, row) in g.iter_mut().enumerate() {
                    for (j, v) in row.iter_mut().enumerate() {
                        *v = (0..n).map(|k| bas[k] * c[k][i][j]).sum();
                    }
                }
                if let Some(t) = rotation {
                    g = rotate_c(&g, t, dim);
                }
                if dim == 2 {
                    let (gc, ezz) = condense_plane(g, plane_stress);
                    g = gc;
                    ezz_c = ezz;
                }
                PointLaw::Aniso { g: Box::new(g) }
            }
        };

        // Displacement gradient and small strain.
        let mut grad = ZERO3;
        for a in 0..cdim {
            for b in 0..cdim {
                grad[a][b] = (0..nd).map(|k| inp.disp[k][a] * db[k][b]).sum();
            }
        }
        if axi {
            grad[2][2] = (0..nd).map(|k| inp.disp[k][0] * full[k]).sum::<f64>() / r;
        }
        let mut strain = ZERO3;
        for i in 0..3 {
            for j in 0..3 {
                strain[i][j] = (grad[j][i] + grad[i][j]) / 2.0;
            }
        }
        // DefG = Identity (+ Grad) and its determinant over `dim`.
        let mut defg = ZERO3;
        for i in 0..3 {
            defg[i][i] = ident[i];
        }
        if inp.large_deflection {
            for i in 0..3 {
                for j in 0..3 {
                    defg[i][j] += grad[i][j];
                }
            }
            let gtg = matmul3(&transpose3(&grad), &grad);
            for i in 0..3 {
                for j in 0..3 {
                    strain[i][j] += gtg[i][j] / 2.0;
                }
            }
        }
        let f = &defg;
        let det_defg = match dim {
            2 => f[0][0] * f[1][1] - f[0][1] * f[1][0],
            _ => {
                f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
                    + f[0][1] * (f[1][2] * f[2][0] - f[1][0] * f[2][2])
                    + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0])
            }
        };

        // The elastic strain: thermal part out, before the out-of-plane recovery.
        if need_heat {
            let t = inp.temperature.unwrap();
            let t_ip: f64 = (0..n).map(|k| t[k] * bas[k]).sum();
            for i in 0..dim {
                let beta = match inp.ip.and_then(|p| p.beta.as_ref()) {
                    Some(b) => b[ipt],
                    None => inp.heat_expansion.map_or(0.0, |b| (0..n).map(|k| b[k][i] * bas[k]).sum()),
                };
                strain[i][i] -= beta * t_ip;
            }
        }

        if plane_stress {
            // `nu` is the nodal interpolation, as in the Fortran, which takes it
            // before the at-IP override.
            strain[2][2] = if isotropic {
                -nu_nodal / (1.0 - nu_nodal) * (strain[0][0] + strain[1][1])
            } else {
                ezz_c[0] * strain[0][0] + ezz_c[1] * strain[1][1] + ezz_c[2] * 2.0 * strain[0][1]
            };
        }

        let s2 = match &law {
            // NeoHookeanStress: S = p C^-1 + mu (I - C^-1), C = F^T F over
            // `dim`, p = lambda/2 (J - 1)(J + 1).
            PointLaw::Iso { lame1, lame2 } if inp.neo_hookean => {
                let pres = lame1 / 2.0 * (det_defg - 1.0) * (det_defg + 1.0);
                let mut red = ZERO3;
                for (i, row) in red.iter_mut().enumerate().take(cdim) {
                    row[i] = 1.0;
                }
                if axi || !plane_stress {
                    red[2][2] = 1.0;
                }
                let mut inv_c = matmul3(&transpose3(&defg), &defg);
                invert_matrix(&mut inv_c, dim)?;
                let mut s = ZERO3;
                for i in 0..3 {
                    for j in 0..3 {
                        s[i][j] = pres * inv_c[i][j] + lame2 * (red[i][j] - inv_c[i][j]);
                    }
                }
                s
            }
            PointLaw::Iso { lame1, lame2 } => {
                let tr: f64 = (0..dim).map(|i| strain[i][i]).sum();
                let mut s = ZERO3;
                for i in 0..3 {
                    for j in 0..3 {
                        s[i][j] = 2.0 * lame2 * strain[i][j];
                    }
                    s[i][i] += lame1 * tr * ident[i];
                }
                s
            }
            PointLaw::Aniso { .. } => law.stress(dim, false, &strain),
        };
        // `Stress2(i,i) = Stress2(i,i) - Pressure`, the nodal interpolation.
        let mut s2 = s2;
        if let Some(pn) = inp.pressure {
            let pres: f64 = (0..n).map(|k| pn[k] * bas[k]).sum();
            for (i, row) in s2.iter_mut().enumerate().take(cdim) {
                row[i] -= pres;
            }
        }
        // `DefG` is the identity of the state of stress here; the push-forward
        // only masks the rows and columns it does not carry.
        let mut stress = ZERO3;
        if general {
            // Stress = 1/DetDefG * MATMUL( MATMUL(DefG, Stress2), TRANSPOSE(DefG) )
            let pushed = matmul3(&matmul3(&defg, &s2), &transpose3(&defg));
            let rdet = 1.0 / det_defg;
            for i in 0..3 {
                for j in 0..3 {
                    stress[i][j] = rdet * pushed[i][j];
                }
            }
        } else {
            for i in 0..3 {
                for j in 0..3 {
                    stress[i][j] = ident[i] * s2[i][j] * ident[j];
                }
            }
        }
        if let (PointLaw::Aniso { g }, 2, false) = (&law, dim, plane_stress) {
            stress[2][2] += g[3][0] * strain[0][0] + g[3][1] * strain[1][1] + g[3][2] * 2.0 * strain[0][1];
        }

        for p in 0..nd {
            for q in 0..nd {
                out.mass[p * nd + q] += w * full[q] * full[p];
            }
        }
        if inp.want_strain {
            add_tensor(&mut out.strain, full, ncomp, w, &strain);
        }
        if inp.want_stress {
            add_tensor(&mut out.stress, full, ncomp, w, &stress);
        }
    }
    Ok(out)
}

/// Symmetric positive definite sparse matrix by rows, for the projection's
/// consistent mass matrix.
#[derive(Clone, Debug, Default)]
pub struct MassSystem {
    rows: Vec<std::collections::BTreeMap<u32, f64>>,
}

impl MassSystem {
    pub fn new(n: usize) -> Self {
        Self { rows: vec![Default::default(); n] }
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Add an element's `n x n` matrix at the given slots; `None` slots (nodes
    /// outside the system) are dropped.
    pub fn add(&mut self, slots: &[Option<usize>], local: &[f64]) {
        let n = slots.len();
        for (p, sp) in slots.iter().enumerate() {
            let Some(sp) = sp else { continue };
            for (q, sq) in slots.iter().enumerate() {
                let Some(sq) = sq else { continue };
                *self.rows[*sp].entry(*sq as u32).or_insert(0.0) += local[p * n + q];
            }
        }
    }

    /// The matrix in CSR form (columns sorted, with the diagonal positions).
    pub fn to_csr(&self) -> elmer_core::CsrMatrix<f64> {
        let n = self.rows.len();
        let (mut row_ptr, mut col_idx, mut values, mut diag) = (vec![0], Vec::new(), Vec::new(), Vec::with_capacity(n));
        for (i, row) in self.rows.iter().enumerate() {
            let mut d = None;
            for (&j, &v) in row {
                if j as usize == i {
                    d = Some(col_idx.len());
                }
                col_idx.push(j);
                values.push(v);
            }
            diag.push(d.unwrap_or(usize::MAX));
            row_ptr.push(col_idx.len());
        }
        let diag = diag.iter().all(|&d| d != usize::MAX).then_some(diag);
        elmer_core::CsrMatrix { nrows: n, ncols: n, row_ptr, col_idx, values, diag }
    }

    fn mul(&self, x: &[f64], y: &mut [f64]) {
        for (yi, row) in y.iter_mut().zip(&self.rows) {
            *yi = row.iter().map(|(&j, &v)| v * x[j as usize]).sum();
        }
    }

    /// Jacobi-preconditioned conjugate gradients to a relative residual of
    /// `1e-14`. Returns the iteration count.
    pub fn solve(&self, b: &[f64], x: &mut [f64]) -> Result<usize, StressError> {
        let n = self.rows.len();
        let diag: Vec<f64> = (0..n)
            .map(|i| self.rows[i].get(&(i as u32)).copied().filter(|d| *d != 0.0).unwrap_or(1.0))
            .collect();
        let bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        x.iter_mut().for_each(|v| *v = 0.0);
        if bnorm == 0.0 {
            return Ok(0);
        }
        let mut r = b.to_vec();
        let mut z: Vec<f64> = r.iter().zip(&diag).map(|(r, d)| r / d).collect();
        let mut p = z.clone();
        let mut ap = vec![0.0; n];
        let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        for it in 1..=(10 * n + 100) {
            self.mul(&p, &mut ap);
            let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
            if pap <= 0.0 {
                return Err(StressError::Input("projection mass matrix is not positive definite".into()));
            }
            let alpha = rz / pap;
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            if r.iter().map(|v| v * v).sum::<f64>().sqrt() <= 1e-14 * bnorm {
                return Ok(it);
            }
            for i in 0..n {
                z[i] = r[i] / diag[i];
            }
            let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
            let beta = rz_new / rz;
            rz = rz_new;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
        }
        Err(StressError::Input("projection solve did not converge".into()))
    }
}

/// `sqrt(3/2 s:s)` of the deviatoric part of a nodal stress.
pub fn von_mises(s: &Tensor3) -> f64 {
    let tr = s[0][0] + s[1][1] + s[2][2];
    let mut acc = 0.0;
    for (i, row) in s.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            let d = v - if i == j { tr / 3.0 } else { 0.0 };
            acc += d * d;
        }
    }
    (1.5 * acc).sqrt()
}

/// Eigenvalues of a symmetric 3x3 tensor, descending (`Principal Stress 1..3`),
/// by cyclic Jacobi rotations.
pub fn principal_values(t: &Tensor3) -> [f64; 3] {
    let mut a = *t;
    for _ in 0..64 {
        let off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        let scale = a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2] + 2.0 * off;
        if off <= 1e-30 * scale || off == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q] == 0.0 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let tt = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let tt = if theta == 0.0 { 1.0 } else { tt };
            let c = 1.0 / (tt * tt + 1.0).sqrt();
            let s = tt * c;
            let mut b = a;
            for k in 0..3 {
                b[k][p] = c * a[k][p] - s * a[k][q];
                b[k][q] = s * a[k][p] + c * a[k][q];
            }
            let c2 = b;
            for k in 0..3 {
                b[p][k] = c * c2[p][k] - s * c2[q][k];
                b[q][k] = s * c2[p][k] + c * c2[q][k];
            }
            a = b;
        }
    }
    let mut w = [a[0][0], a[1][1], a[2][2]];
    w.sort_by(|x, y| y.total_cmp(x));
    w
}

/// `Tresca`: half the largest difference of principal stresses.
pub fn tresca(w: &[f64; 3]) -> f64 {
    let m = (w[0] - w[1]) / 2.0;
    let m = m.max((w[1] - w[2]) / 2.0);
    m.max((w[0] - w[2]) / 2.0)
}

/// The solve of the stress projection. On p elements it is
/// `NodalProjectorSolve`'s, with the section's own linear system (its
/// `StressTemp` solver); the nodal projection keeps [`MassSystem::solve`].
pub enum Projector<'m> {
    Cg(&'m MassSystem),
    Direct(Box<dyn elmer_core::DirectSolver<f64> + Send>),
    Iterative(
        Box<dyn elmer_core::IterativeSolver<f64> + Send>,
        Box<dyn elmer_core::Preconditioner<f64> + Send>,
        elmer_core::CsrMatrix<f64>,
    ),
}

impl<'m> Projector<'m> {
    /// The projector of `mass` for the running solver; `p_elements` selects
    /// the section's linear system.
    pub fn new(
        ctx: &elmer_solver::SolverCtx<'_>,
        mass: &'m MassSystem,
        p_elements: bool,
    ) -> Result<Self, elmer_solver::SolverError> {
        if !p_elements {
            return Ok(Projector::Cg(mass));
        }
        let a = mass.to_csr();
        let params = ctx.params().clone();
        Ok(match ctx.sim.backends.linear.select_for(&params, &a, ctx.log)? {
            elmer_solver::backend::LinearSolverChoice::Direct(mut d) => {
                d.factor(&a).map_err(|e| elmer_solver::SolverError::fatal("NodalProjectorSolve", e.to_string()))?;
                Projector::Direct(d)
            }
            elmer_solver::backend::LinearSolverChoice::Iterative { solver, precond } => {
                let m: Box<dyn elmer_core::Preconditioner<f64> + Send> = match precond {
                    Some(f) => f(&params, &a)?,
                    None => Box::new(elmer_core::IdentityPrecond),
                };
                Projector::Iterative(solver, m, a)
            }
        })
    }

    /// `x = M^-1 b`.
    pub fn solve(&mut self, comm: &dyn elmer_core::Communicator, b: &[f64], x: &mut [f64]) -> Result<(), elmer_solver::SolverError> {
        let fatal = |e: String| elmer_solver::SolverError::fatal("NodalProjectorSolve", e);
        match self {
            Projector::Cg(m) => m.solve(b, x).map(|_| ()).map_err(|e| fatal(e.to_string())),
            Projector::Direct(d) => {
                x.copy_from_slice(b);
                d.solve_in_place(x).map_err(|e| fatal(e.to_string()))
            }
            Projector::Iterative(s, m, a) => {
                x.iter_mut().for_each(|v| *v = 0.0);
                s.solve(&*a, &**m, b, x, comm).map(|_| ()).map_err(elmer_solver::SolverError::from)
            }
        }
    }
}
