//! Stress and strain postprocessing of `ElasticSolve.F90` `ComputeStressAndStrain`
//! (`:5759`): the element contribution to the L2 nodal projection of the
//! small-strain `Strain` and `Stress` fields, the projection solve, von Mises,
//! and the principal values.
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

/// Element contribution to the L2 projection of strain and stress.
pub fn stress_strain_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
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
    if inp.xyz.len() != n || inp.disp.len() != n {
        return Err(StressError::Input(format!("element {} needs {n} nodal values", inp.etype.0)));
    }
    let dim = if axi { 3 } else { cdim };
    let plane_stress = inp.plane_stress && !axi && cdim == 2;
    let ncomp = output_components(cdim);
    let isotropic = matches!(inp.law, ElasticLaw::Isotropic { .. });

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
                (l1, e / (2.0 * (1.0 + nu)))
            })
            .unzip(),
        ElasticLaw::Anisotropic { .. } => (Vec::new(), Vec::new()),
    };

    let need_heat = inp.temperature.is_some_and(|t| t.iter().any(|x| *x != 0.0))
        && (inp.heat_expansion.is_some() || inp.ip.is_some_and(|p| p.beta.is_some()));

    let points = quad.gauss_points(inp.etype, inp.np);
    if let Some(p) = inp.ip {
        let np = points.len();
        let lens = [p.young.as_ref().map(Vec::len), p.poisson.as_ref().map(Vec::len), p.beta.as_ref().map(Vec::len)];
        if lens.iter().flatten().any(|&l| l != np) {
            return Err(StressError::Input(format!("integration point data must have {np} entries")));
        }
    }

    let mut out = StressPostLocal {
        mass: vec![0.0; n * n],
        strain: if inp.want_strain { vec![0.0; n * ncomp] } else { Vec::new() },
        stress: if inp.want_stress { vec![0.0; n * ncomp] } else { Vec::new() },
    };
    let mut be = BasisEval::default();

    for (ipt, qp) in points.iter().enumerate() {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be).map_err(StressError::Basis)?;
        let bas = &be.basis[..n];
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
                if let Some(p) = inp.ip.filter(|p| p.young.is_some() || p.poisson.is_some()) {
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
                grad[a][b] = (0..n).map(|k| inp.disp[k][a] * db[k][b]).sum();
            }
        }
        if axi {
            grad[2][2] = (0..n).map(|k| inp.disp[k][0] * bas[k]).sum::<f64>() / r;
        }
        let mut strain = ZERO3;
        for i in 0..3 {
            for j in 0..3 {
                strain[i][j] = (grad[j][i] + grad[i][j]) / 2.0;
            }
        }

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
        // `DefG` is the identity of the state of stress here; the push-forward
        // only masks the rows and columns it does not carry.
        let mut stress = ZERO3;
        for i in 0..3 {
            for j in 0..3 {
                stress[i][j] = ident[i] * s2[i][j] * ident[j];
            }
        }
        if let (PointLaw::Aniso { g }, 2, false) = (&law, dim, plane_stress) {
            stress[2][2] += g[3][0] * strain[0][0] + g[3][1] * strain[1][1] + g[3][2] * 2.0 * strain[0][1];
        }

        for p in 0..n {
            for q in 0..n {
                out.mass[p * n + q] += w * bas[q] * bas[p];
            }
        }
        if inp.want_strain {
            add_tensor(&mut out.strain, bas, ncomp, w, &strain);
        }
        if inp.want_stress {
            add_tensor(&mut out.stress, bas, ncomp, w, &stress);
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
