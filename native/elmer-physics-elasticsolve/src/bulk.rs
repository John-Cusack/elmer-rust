//! Bulk element system of `ElasticSolve.F90` `LocalMatrix` (`:3716`) for a linear
//! elastic material, with or without `Large Deflection`.
//!
//! The assembly is a Newton step written for the new displacement: the stiffness
//! is the tangent at the current iterate and the right-hand side carries the
//! residual plus the tangent applied to that iterate (`-Stress1 : grad(v) +
//! dStress1U : grad(v)`). Under `Large Deflection` the strain is the
//! Green-Lagrange one, so the material is St. Venant-Kirchhoff; without it the
//! deformation gradient stays the identity and one step is the linear solution.
//!
//! DOF layout is node-major, component-minor: row `cdim*p + i` is component `i`
//! of node `p` (Elmer's `DOFs*(p-1)+i`).

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};
use elmer_physics_elastic::{CoordSystem, ElasticLaw, IpProps};

use crate::error::ElasticError;
use crate::law::{condense_plane, rotate_c, voigt_to_tensor, PointLaw, Voigt6};
use crate::tensor::{self as t3, identity, T3, ZERO};

/// Everything `LocalMatrix` reads for one bulk element, evaluated at the element
/// nodes (`n` values each; Lagrange elements, so `ntot == n`).
#[derive(Clone, Debug)]
pub struct BulkInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    /// `CoordinateSystemDimension()`: 2 or 3 (2 for axisymmetric as well).
    pub cdim: usize,
    pub coord: CoordSystem,
    /// `Plane Stress` of the Equation section; Cartesian 2D only.
    pub plane_stress: bool,
    /// `Large Deflection`.
    pub large_deflection: bool,
    pub law: &'a ElasticLaw,
    /// Current displacement at the nodes, `cdim` components used.
    pub disp: &'a [[f64; 3]],
    pub density: &'a [f64],
    pub damping: Option<&'a [f64]>,
    /// `Stress Bodyforce 1..3`.
    pub load: &'a [[f64; 3]],
    /// `Inertial Bodyforce 1..3`, scaled by the density.
    pub inertial: Option<&'a [[f64; 3]]>,
    /// `Stress Pressure`: enters as `p div(v)`.
    pub pressure_load: Option<&'a [f64]>,
    /// `Stress Load` in Voigt form (xx yy zz xy yz xz).
    pub stress_load: Option<&'a [Voigt6]>,
    /// `Strain Load` in Voigt form with engineering shear.
    pub strain_load: Option<&'a [Voigt6]>,
    /// Temperature minus `Reference Temperature`.
    pub temperature: Option<&'a [f64]>,
    /// Diagonal of `Heat Expansion Coefficient`.
    pub heat_expansion: Option<&'a [[f64; 3]]>,
    pub want_mass: bool,
    /// Gauss points per element; `None` is the Elmer default rule.
    pub np: Option<usize>,
    pub ip: Option<&'a IpProps>,
}

/// The displacement gradient `Grad(a,b) = d u_a / d x_b` (with the hoop term
/// `u_r / r` at `(3,3)` for axisymmetry).
pub(crate) fn displacement_gradient(
    disp: &[[f64; 3]],
    bas: &[f64],
    db: &[[f64; 3]],
    dim: usize,
    axi: bool,
    r: f64,
) -> T3 {
    let n = disp.len();
    let mut g = ZERO;
    if axi {
        g[0][0] = (0..n).map(|k| disp[k][0] * db[k][0]).sum();
        g[0][1] = (0..n).map(|k| disp[k][0] * db[k][1]).sum();
        g[2][2] = (0..n).map(|k| disp[k][0] * bas[k]).sum::<f64>() / r;
        g[1][0] = (0..n).map(|k| disp[k][1] * db[k][0]).sum();
        g[1][1] = (0..n).map(|k| disp[k][1] * db[k][1]).sum();
    } else {
        for a in 0..dim {
            for b in 0..dim {
                g[a][b] = (0..n).map(|k| disp[k][a] * db[k][b]).sum();
            }
        }
    }
    g
}

/// The gradient of the trial function "node `p`, component `i`".
pub(crate) fn trial_gradient(i: usize, p: usize, bas: &[f64], db: &[[f64; 3]], axi: bool, r: f64) -> T3 {
    let mut d = ZERO;
    if axi {
        match i {
            0 => {
                d[0][0] = db[p][0];
                d[0][1] = db[p][1];
                d[2][2] = bas[p] / r;
            }
            _ => {
                d[1][0] = db[p][0];
                d[1][1] = db[p][1];
            }
        }
    } else {
        d[i] = db[p];
    }
    d
}

/// Bulk element matrix and force. Overwrites `out`.
pub fn elastic_bulk_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &BulkInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), ElasticError> {
    let n = inp.etype.n_nodes();
    let cdim = inp.cdim;
    let axi = inp.coord == CoordSystem::AxiSymmetric;
    if cdim != 2 && cdim != 3 {
        return Err(ElasticError::Unsupported(format!("coordinate system dimension {cdim}")));
    }
    if axi && cdim != 2 {
        return Err(ElasticError::Unsupported("axisymmetric elasticity needs a 2D mesh".into()));
    }
    if inp.xyz.len() != n || inp.load.len() != n || inp.density.len() != n || inp.disp.len() != n {
        return Err(ElasticError::Input(format!("element {} needs {n} nodal values", inp.etype.0)));
    }
    // Dimension of the state of stress: three for axisymmetry (r, z, hoop).
    let dim = if axi { 3 } else { cdim };
    let plane_stress = inp.plane_stress && !axi && cdim == 2;
    let ndof = n * cdim;
    let large = inp.large_deflection;

    *out = LocalSystem::zeros(ndof);
    if inp.want_mass {
        out.mass = Some(vec![0.0; ndof * ndof]);
        if inp.damping.is_some() {
            out.damp = Some(vec![0.0; ndof * ndof]);
        }
    }

    // Lame parameters are interpolated from nodal values, as Elmer does.
    let (lame1_n, lame2_n): (Vec<f64>, Vec<f64>) = match inp.law {
        ElasticLaw::Isotropic { young, poisson } => {
            if young.len() != n || poisson.len() != n {
                return Err(ElasticError::Input("nodal moduli length".into()));
            }
            (0..n)
                .map(|k| {
                    let (e, nu) = (young[k], poisson[k]);
                    let l1 = if plane_stress {
                        e * nu / (1.0 - nu * nu)
                    } else {
                        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
                    };
                    (l1, e / (2.0 * (1.0 + nu)))
                })
                .unzip()
        }
        ElasticLaw::Anisotropic { c, .. } => {
            if c.len() != n {
                return Err(ElasticError::Input("nodal moduli length".into()));
            }
            (Vec::new(), Vec::new())
        }
    };

    let ip = inp.ip;
    let beta_ip = ip.is_some_and(|p| p.beta.is_some());
    let has_heat = match (inp.temperature, inp.heat_expansion) {
        (Some(t), b) => {
            t.iter().any(|x| *x != 0.0) && (beta_ip || b.is_some_and(|b| b.iter().any(|v| v.iter().any(|x| *x != 0.0))))
        }
        _ => false,
    };
    let has_voigt = inp.stress_load.is_some_and(|s| s.iter().any(|v| v.iter().any(|x| *x != 0.0)))
        || inp.strain_load.is_some_and(|s| s.iter().any(|v| v.iter().any(|x| *x != 0.0)));
    let got_offset = has_voigt || has_heat;

    let points = quad.gauss_points(inp.etype, inp.np);
    if let Some(p) = ip {
        let np = points.len();
        let lens = [
            p.young.as_ref().map(Vec::len),
            p.poisson.as_ref().map(Vec::len),
            p.beta.as_ref().map(Vec::len),
            p.load.as_ref().map(Vec::len),
            p.pressure.as_ref().map(Vec::len),
        ];
        if lens.iter().flatten().any(|&l| l != np) {
            return Err(ElasticError::Input(format!("integration point data must have {np} entries")));
        }
    }

    let eye = identity(dim);
    let mut be = BasisEval::default();
    let mut ddefg: Vec<T3> = vec![ZERO; ndof];
    let mut dstresses: Vec<T3> = vec![ZERO; ndof];

    for (ipt, qp) in points.iter().enumerate() {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be).map_err(ElasticError::Basis)?;
        let bas = &be.basis;
        let db = &be.dbasis_dx;
        let mut s = be.det_j * qp.weight;
        let mut r = 0.0;
        if axi {
            r = (0..n).map(|k| bas[k] * inp.xyz[k][0]).sum();
            s *= r;
        }

        let mut force = [0.0; 3];
        let mut inertial = [0.0; 3];
        for i in 0..cdim {
            force[i] = (0..n).map(|k| inp.load[k][i] * bas[k]).sum();
            if let Some(w) = inp.inertial {
                inertial[i] = (0..n).map(|k| w[k][i] * bas[k]).sum();
            }
        }
        if let Some(l) = ip.and_then(|p| p.load.as_ref()) {
            force[..cdim].copy_from_slice(&l[ipt][..cdim]);
        }
        let density: f64 = (0..n).map(|k| inp.density[k] * bas[k]).sum();
        let damping: f64 = inp.damping.map_or(0.0, |d| (0..n).map(|k| d[k] * bas[k]).sum());
        let mut pressure: f64 = inp.pressure_load.map_or(0.0, |p| (0..n).map(|k| p[k] * bas[k]).sum());
        if let Some(p) = ip.and_then(|p| p.pressure.as_ref()) {
            pressure = p[ipt];
        }

        let mut heat_exp = ZERO;
        let mut temperature = 0.0;
        if has_heat {
            let tt = inp.temperature.unwrap();
            temperature = (0..n).map(|k| tt[k] * bas[k]).sum();
            if let Some(b) = ip.and_then(|p| p.beta.as_ref()) {
                for (i, row) in heat_exp.iter_mut().enumerate() {
                    row[i] = b[ipt];
                }
            } else if let Some(b) = inp.heat_expansion {
                for i in 0..3 {
                    heat_exp[i][i] = (0..n).map(|k| b[k][i] * bas[k]).sum();
                }
            }
        }
        let law = match inp.law {
            ElasticLaw::Isotropic { young, poisson } => {
                let mut lame1 = (0..n).map(|k| lame1_n[k] * bas[k]).sum();
                let mut lame2 = (0..n).map(|k| lame2_n[k] * bas[k]).sum();
                if let Some(p) = ip.filter(|p| p.young.is_some() || p.poisson.is_some()) {
                    let e = p.young.as_ref().map_or_else(|| (0..n).map(|k| young[k] * bas[k]).sum(), |v| v[ipt]);
                    let nu: f64 =
                        p.poisson.as_ref().map_or_else(|| (0..n).map(|k| poisson[k] * bas[k]).sum(), |v| v[ipt]);
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
                if let Some(rot) = rotation {
                    g = rotate_c(&g, rot, dim);
                }
                if has_heat && dim == 2 && !plane_stress {
                    let det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
                    let b33 = heat_exp[2][2];
                    heat_exp[0][0] += b33 * (g[1][1] * g[0][2] - g[0][1] * g[1][2]) / det;
                    heat_exp[1][1] += b33 * (g[0][0] * g[1][2] - g[0][1] * g[0][2]) / det;
                }
                if dim == 2 {
                    g = condense_plane(g, plane_stress);
                }
                PointLaw::Aniso { g: Box::new(g) }
            }
        };

        // The affine part of the response: sigma = C(eps - eps0) + sigma0.
        let mut stress_offset = ZERO;
        if got_offset {
            let interp = |v: Option<&[Voigt6]>| -> Voigt6 {
                let mut o = [0.0; 6];
                if let Some(v) = v {
                    for (i, oi) in o.iter_mut().enumerate() {
                        *oi = (0..n).map(|k| v[k][i] * bas[k]).sum();
                    }
                }
                o
            };
            stress_offset = voigt_to_tensor(&interp(inp.stress_load), cdim, axi);
            let mut eigen = voigt_to_tensor(&interp(inp.strain_load), cdim, axi);
            for (i, row) in eigen.iter_mut().enumerate() {
                for (j, v) in row.iter_mut().enumerate() {
                    if i != j {
                        *v /= 2.0;
                    }
                }
            }
            if has_heat {
                for i in 0..dim {
                    eigen[i][i] += heat_exp[i][i] * temperature;
                }
            }
            let resp = law.stress(dim, false, &eigen);
            stress_offset = t3::sub(&stress_offset, &resp);
        }

        let grad = displacement_gradient(inp.disp, bas, db, dim, axi, r);
        let defg = if large { t3::add(&eye, &grad) } else { eye };
        let det_defg = t3::det(&defg, dim);
        let mut strain = t3::symmetrize(&grad);
        if large {
            strain = t3::add(&strain, &t3::scale(&t3::mul(&t3::transpose(&grad), &grad), 0.5));
        }
        let mut stress2 = law.stress(dim, false, &strain);
        if got_offset {
            stress2 = t3::add(&stress2, &stress_offset);
        }
        let stress1 = t3::mul(&defg, &stress2);

        let dstrain_of = |d: &T3| -> T3 {
            t3::scale(&t3::add(&t3::mul(&t3::transpose(&defg), d), &t3::mul(&t3::transpose(d), &defg)), 0.5)
        };
        let dstress2_u = law.stress(dim, false, &dstrain_of(&grad));
        let mut dstress1_u = t3::mul(&defg, &dstress2_u);
        if large {
            dstress1_u = t3::add(&dstress1_u, &t3::mul(&grad, &stress2));
        }

        for p in 0..n {
            for i in 0..cdim {
                let k = cdim * p + i;
                ddefg[k] = trial_gradient(i, p, bas, db, axi, r);
                dstresses[k] = law.stress(dim, true, &dstrain_of(&ddefg[k]));
            }
        }

        for p in 0..n {
            for i in 0..cdim {
                let k = cdim * p + i;
                let d_defg = &ddefg[k];
                let mut dstress1 = t3::mul(&defg, &dstresses[k]);
                if large {
                    dstress1 = t3::add(&dstress1, &t3::mul(d_defg, &stress2));
                }

                out.force[k] += (bas[p] * force[i] * det_defg
                    + bas[p] * inertial[i] * density
                    + pressure * db[p][i]
                    - t3::ddot(d_defg, &stress1, 3)
                    + t3::ddot(d_defg, &dstress1_u, 3))
                    * s;

                for q in 0..n {
                    for j in 0..cdim {
                        let v = if axi {
                            match j {
                                0 => db[q][0] * dstress1[0][0] + db[q][1] * dstress1[0][1] + bas[q] / r * dstress1[2][2],
                                _ => db[q][0] * dstress1[1][0] + db[q][1] * dstress1[1][1],
                            }
                        } else {
                            db[q][0] * dstress1[j][0] + db[q][1] * dstress1[j][1] + db[q][2] * dstress1[j][2]
                        };
                        out.stiff[k * ndof + cdim * q + j] += v * s;
                    }
                }
            }
        }

        if inp.want_mass {
            for p in 0..n {
                for q in 0..n {
                    for i in 0..cdim {
                        let idx = (cdim * p + i) * ndof + cdim * q + i;
                        let mm = out.mass.as_mut().unwrap();
                        mm[idx] += bas[p] * bas[q] * density * s;
                        if let Some(dm) = out.damp.as_mut() {
                            dm[idx] += bas[p] * bas[q] * damping * s;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}
