//! Bulk element matrix and force of the linear small-strain path of
//! `ElasticSolve.F90` `LocalMatrix` (`:3716`), which is what `StressSolve`
//! delegates to unless `Legacy Assembly` is set (`StressSolve.F90:213`).
//!
//! DOF layout is node-major, component-minor: row `cdim*p + i` is component
//! `i` of node `p` (Elmer's `DOFs*(p-1)+i`).

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, QuadPoint, Quadrature};

use crate::error::StressError;
use crate::material::{
    condense_plane, rotate_c, ElasticLaw, PointLaw, Tensor3, Voigt6, I1, I2, ZERO3,
};

/// `CurrentCoordinateSystem()`: only the two systems the linear solver supports.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CoordSystem {
    Cartesian,
    /// `Coordinate System = Axi Symmetric`: x is the radius, y the axis.
    AxiSymmetric,
}

/// Material and load data asked for at the integration points rather than
/// interpolated from the nodes (`Youngs Modulus at IP`, `Poisson Ratio at IP`,
/// `Heat Expansion Coefficient IP`, `Stress Bodyforce at IP`). Each vector has one
/// entry per Gauss point, in the order of `Quadrature::gauss_points`; a quantity
/// left `None` is interpolated from the nodes as usual. Isotropic materials only,
/// as in `ElasticSolve.F90`.
#[derive(Clone, Debug, Default)]
pub struct IpProps {
    pub young: Option<Vec<f64>>,
    pub poisson: Option<Vec<f64>>,
    /// One value, applied to the three diagonal entries of the expansion tensor.
    pub beta: Option<Vec<f64>>,
    /// `Stress Bodyforce 1..3` per point; replaces the nodal load.
    pub load: Option<Vec<[f64; 3]>>,
    /// `Stress Pressure` per point; read with the load, as `LoadIP_h(4)` is.
    pub pressure: Option<Vec<f64>>,
}

/// Everything `LocalMatrix` reads for one bulk element, already evaluated at
/// the element nodes (`n` values each; Lagrange elements, so `ntot == n`).
#[derive(Clone, Debug)]
pub struct StressBulkInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    /// `CoordinateSystemDimension()`: 2 or 3 (2 for axisymmetric as well).
    pub cdim: usize,
    pub coord: CoordSystem,
    /// `Plane Stress` of the Equation section. Meaningful for `cdim == 2`
    /// Cartesian only, as in Elmer.
    pub plane_stress: bool,
    pub law: &'a ElasticLaw,
    /// `Density`, per node. Used for the mass matrix only.
    pub density: &'a [f64],
    /// `Damping`, per node.
    pub damping: Option<&'a [f64]>,
    /// `Stress Bodyforce 1..3`, per node.
    pub load: &'a [[f64; 3]],
    /// `Stress Pressure`: enters as `p div(v)`, per node.
    pub pressure_load: Option<&'a [f64]>,
    /// `Stress Load` in Voigt form (xx yy zz xy yz xz), per node.
    pub stress_load: Option<&'a [Voigt6]>,
    /// `Strain Load` in Voigt form with engineering shear, per node.
    pub strain_load: Option<&'a [Voigt6]>,
    /// Temperature minus `Reference Temperature`, per node.
    pub temperature: Option<&'a [f64]>,
    /// Diagonal of `Heat Expansion Coefficient`, per node.
    pub heat_expansion: Option<&'a [[f64; 3]]>,
    /// Assemble `out.mass` (and `out.damp` when `damping` is given).
    pub want_mass: bool,
    /// Gauss points per element; `None` is the Elmer default rule.
    pub np: Option<usize>,
    /// Quantities given at the integration points.
    pub ip: Option<&'a IpProps>,
}

/// Shape function values at each Gauss point of the element's rule, `n_ip` rows
/// of `n` values, row-major: what `elmer_expr::Compiled::eval_ips` interpolates
/// its dependencies with.
pub fn shape_at_points(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    np: Option<usize>,
) -> Result<Vec<f64>, StressError> {
    let mut be = BasisEval::default();
    let mut out = Vec::new();
    for qp in quad.gauss_points(etype, np) {
        basis.eval(etype, xyz, qp.uvw, &mut be).map_err(StressError::Basis)?;
        out.extend_from_slice(&be.basis[..xyz.len()]);
    }
    Ok(out)
}

fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `Stress.F90:1908` `Vector62Tensor`.
fn voigt_to_tensor(v: &Voigt6, cdim: usize, axi: bool) -> Tensor3 {
    let n = match (cdim, axi) {
        (2, true) => 4,
        (2, false) => 3,
        (3, _) => 6,
        _ => 1,
    };
    // The axisymmetric 2D packing is 11, 22, 33, 12.
    let (i1, i2): ([usize; 6], [usize; 6]) = if cdim == 2 && axi {
        ([0, 1, 2, 0, 0, 0], [0, 1, 2, 1, 0, 0])
    } else {
        (I1, I2)
    };
    let mut x = ZERO3;
    for i in 0..n {
        x[i1[i]][i2[i]] = v[i];
        x[i2[i]][i1[i]] = v[i];
    }
    x
}

fn symmetrize(g: &Tensor3) -> Tensor3 {
    let mut e = ZERO3;
    for i in 0..3 {
        for j in 0..3 {
            e[i][j] = (g[i][j] + g[j][i]) / 2.0;
        }
    }
    e
}

/// Bulk element matrix and force. Overwrites `out`.
pub fn stress_bulk_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &StressBulkInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), StressError> {
    let points = quad.gauss_points(inp.etype, inp.np);
    bulk_impl(basis, &points, inp.etype.n_nodes(), inp, out)
}

/// [`stress_bulk_local`] for `nd >= n` basis functions per component (the
/// nodal ones first, then the edge, face and bubble functions of a p
/// element) on the given rule, as `StressCompose(..., n, ntot, ...)`:
/// material values and loads are interpolated with the nodal functions only.
/// `out` is `nd * cdim` square, interleaved by component. Integration-point
/// properties (`inp.ip`) are refused: they are tabulated on the default rule.
pub fn stress_bulk_local_p(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    nd: usize,
    inp: &StressBulkInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), StressError> {
    if inp.ip.is_some() {
        return Err(StressError::Unsupported(
            "integration point properties on p elements".into(),
        ));
    }
    if nd < inp.etype.n_nodes() {
        return Err(StressError::Input(format!("{nd} basis functions for {} nodes", inp.etype.n_nodes())));
    }
    bulk_impl(basis, points, nd, inp, out)
}

fn bulk_impl(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    nd: usize,
    inp: &StressBulkInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), StressError> {
    let n = inp.etype.n_nodes();
    let cdim = inp.cdim;
    let axi = inp.coord == CoordSystem::AxiSymmetric;
    if cdim != 2 && cdim != 3 {
        return Err(StressError::Unsupported(format!(
            "coordinate system dimension {cdim}"
        )));
    }
    if axi && cdim != 2 {
        return Err(StressError::Unsupported(
            "axisymmetric elasticity needs a 2D mesh".into(),
        ));
    }
    if inp.xyz.len() != n || inp.load.len() != n || inp.density.len() != n {
        return Err(StressError::Input(format!(
            "element {} needs {n} nodal values, got xyz {} load {} density {}",
            inp.etype.0,
            inp.xyz.len(),
            inp.load.len(),
            inp.density.len()
        )));
    }
    // Dimension of the state of stress: three for axisymmetry (r, z, hoop).
    let dim = if axi { 3 } else { cdim };
    let plane_stress = inp.plane_stress && !axi && cdim == 2;
    let ndof = nd * cdim;

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
                return Err(StressError::Input("nodal moduli length".into()));
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
                return Err(StressError::Input("nodal moduli length".into()));
            }
            (Vec::new(), Vec::new())
        }
    };

    let ip = inp.ip;
    let beta_ip = ip.is_some_and(|p| p.beta.is_some());
    let has_heat = match (inp.temperature, inp.heat_expansion) {
        (Some(t), b) => {
            t.iter().any(|x| *x != 0.0)
                && (beta_ip || b.is_some_and(|b| b.iter().any(|v| v.iter().any(|x| *x != 0.0))))
        }
        _ => false,
    };
    let has_voigt = inp
        .stress_load
        .is_some_and(|s| s.iter().any(|v| v.iter().any(|x| *x != 0.0)))
        || inp
            .strain_load
            .is_some_and(|s| s.iter().any(|v| v.iter().any(|x| *x != 0.0)));
    let got_offset = has_voigt || has_heat;

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
            return Err(StressError::Input(format!("integration point data must have {np} entries")));
        }
    }
    let mut be = BasisEval::default();
    // Reused across points and trial functions.
    let mut dstress: Vec<Tensor3> = vec![ZERO3; ndof];
    let mut ddefg: Vec<Tensor3> = vec![ZERO3; ndof];

    for (ipt, qp) in points.iter().enumerate() {
        basis
            .eval(inp.etype, inp.xyz, qp.uvw, &mut be)
            .map_err(StressError::Basis)?;
        if be.basis.len() < nd {
            return Err(StressError::Input(format!("{} basis functions for {nd} dofs", be.basis.len())));
        }
        let bas = &be.basis;
        let db = &be.dbasis_dx;
        let mut s = be.det_j * qp.weight;
        let mut r = 0.0;
        if axi {
            r = (0..n).map(|k| bas[k] * inp.xyz[k][0]).sum();
            s *= r;
        }

        let mut force = [0.0; 3];
        for (i, f) in force.iter_mut().enumerate().take(cdim) {
            *f = (0..n).map(|k| inp.load[k][i] * bas[k]).sum();
        }
        if let Some(l) = ip.and_then(|p| p.load.as_ref()) {
            force[..cdim].copy_from_slice(&l[ipt][..cdim]);
        }
        let density: f64 = (0..n).map(|k| inp.density[k] * bas[k]).sum();
        let damping: f64 = inp
            .damping
            .map_or(0.0, |d| (0..n).map(|k| d[k] * bas[k]).sum());
        let mut pressure: f64 = inp
            .pressure_load
            .map_or(0.0, |p| (0..n).map(|k| p[k] * bas[k]).sum());
        if let Some(p) = ip.and_then(|p| p.pressure.as_ref()) {
            pressure = p[ipt];
        }

        // Constitutive law at the point.
        let mut heat_exp = ZERO3;
        let mut temperature = 0.0;
        if has_heat {
            let t = inp.temperature.unwrap();
            temperature = (0..n).map(|k| t[k] * bas[k]).sum();
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
                if let Some(t) = rotation {
                    g = rotate_c(&g, t, dim);
                }
                if has_heat && dim == 2 && !plane_stress {
                    let det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
                    let b33 = heat_exp[2][2];
                    heat_exp[0][0] += b33 * (g[1][1] * g[0][2] - g[0][1] * g[1][2]) / det;
                    heat_exp[1][1] += b33 * (g[0][0] * g[1][2] - g[0][1] * g[0][2]) / det;
                }
                if dim == 2 {
                    g = condense_plane(g, plane_stress).0;
                }
                PointLaw::Aniso { g: Box::new(g) }
            }
        };

        // Stress offset: sigma = C(eps - eps0) + sigma0 at zero displacement.
        let mut stress_offset = ZERO3;
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
            let sl = interp(inp.stress_load);
            let el = interp(inp.strain_load);
            stress_offset = voigt_to_tensor(&sl, cdim, axi);
            let mut eigen = voigt_to_tensor(&el, cdim, axi);
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
            for i in 0..3 {
                for j in 0..3 {
                    stress_offset[i][j] -= resp[i][j];
                }
            }
        }

        // Trial strain fields, one per (node, component).
        for p in 0..nd {
            for i in 0..cdim {
                let k = cdim * p + i;
                let mut d = ZERO3;
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
                ddefg[k] = d;
                dstress[k] = law.stress(dim, true, &symmetrize(&d));
            }
        }

        for p in 0..nd {
            for i in 0..cdim {
                let k = cdim * p + i;
                let offset_term = if axi {
                    (0..3)
                        .map(|a| (0..3).map(|b| ddefg[k][a][b] * stress_offset[a][b]).sum::<f64>())
                        .sum::<f64>()
                } else {
                    dot3(&db[p], &stress_offset[i])
                };
                out.force[k] +=
                    (bas[p] * force[i] + pressure * db[p][i] - offset_term) * s;
                let ds = &dstress[k];
                for q in 0..nd {
                    for j in 0..cdim {
                        let col = cdim * q + j;
                        let v = if axi {
                            match j {
                                0 => db[q][0] * ds[0][0] + db[q][1] * ds[0][1] + bas[q] / r * ds[2][2],
                                _ => db[q][0] * ds[1][0] + db[q][1] * ds[1][1],
                            }
                        } else {
                            dot3(&db[q], &ds[j])
                        };
                        out.stiff[k * ndof + col] += v * s;
                    }
                }
            }
        }

        if inp.want_mass {
            for p in 0..nd {
                for q in 0..nd {
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
