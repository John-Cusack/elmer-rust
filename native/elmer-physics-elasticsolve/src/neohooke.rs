//! Bulk element system of `ElasticSolve.F90` `NeoHookeanLocalMatrix` (`:4863`),
//! the compressible neo-Hookean material with the second Piola-Kirchhoff stress
//! `S = Lame1/2 (J^2 - 1) C^-1 + Lame2 (I - C^-1)`, linearised exactly.
//!
//! The assembly has the same Newton form as [`crate::bulk`]. Temperature,
//! `Stress Pressure`, `Stress Load` and `Strain Load` play no part in it; the
//! solver refuses the keywords the Fortran refuses. The mixed (pressure) form is
//! not ported.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};
use elmer_physics_elastic::CoordSystem;

use crate::bulk::{displacement_gradient, trial_gradient};
use crate::error::ElasticError;
use crate::tensor::{self as t3, identity, T3};

#[derive(Clone, Debug)]
pub struct NeoHookeInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    pub cdim: usize,
    pub coord: CoordSystem,
    pub plane_stress: bool,
    pub young: &'a [f64],
    pub poisson: &'a [f64],
    pub disp: &'a [[f64; 3]],
    pub density: &'a [f64],
    pub damping: Option<&'a [f64]>,
    pub load: &'a [[f64; 3]],
    pub inertial: Option<&'a [[f64; 3]]>,
    pub want_mass: bool,
    pub np: Option<usize>,
}

/// Bulk element matrix and force. Overwrites `out`.
pub fn neo_hookean_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &NeoHookeInput<'_>,
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
    if [inp.xyz.len(), inp.young.len(), inp.poisson.len(), inp.disp.len(), inp.density.len(), inp.load.len()]
        .iter()
        .any(|&l| l != n)
    {
        return Err(ElasticError::Input(format!("element {} needs {n} nodal values", inp.etype.0)));
    }
    let dim = if axi { 3 } else { cdim };
    let ndof = n * cdim;
    *out = LocalSystem::zeros(ndof);
    if inp.want_mass {
        out.mass = Some(vec![0.0; ndof * ndof]);
        if inp.damping.is_some() {
            out.damp = Some(vec![0.0; ndof * ndof]);
        }
    }

    let (lame1_n, lame2_n): (Vec<f64>, Vec<f64>) = (0..n)
        .map(|k| {
            let (e, nu) = (inp.young[k], inp.poisson[k]);
            let l1 = if inp.plane_stress && !axi {
                e * nu / (1.0 - nu * nu)
            } else {
                e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
            };
            (l1, e / (2.0 * (1.0 + nu)))
        })
        .unzip();

    let eye = identity(dim);
    let singular = || ElasticError::Numerical("singular deformation gradient".into());
    let mut be = BasisEval::default();

    for qp in quad.gauss_points(inp.etype, inp.np) {
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
        let lame1: f64 = (0..n).map(|k| lame1_n[k] * bas[k]).sum();
        let lame2: f64 = (0..n).map(|k| lame2_n[k] * bas[k]).sum();
        let density: f64 = (0..n).map(|k| inp.density[k] * bas[k]).sum();
        let damping: f64 = inp.damping.map_or(0.0, |d| (0..n).map(|k| d[k] * bas[k]).sum());

        let grad = displacement_gradient(inp.disp, bas, db, dim, axi, r);
        let defg = t3::add(&eye, &grad);
        let det_defg = t3::det(&defg, dim);
        let c = t3::mul(&t3::transpose(&defg), &defg);
        let invc = t3::invert(&c, dim).ok_or_else(singular)?;
        let invdefg = t3::invert(&defg, dim).ok_or_else(singular)?;
        let jj = 0.5 * (det_defg - 1.0) * (det_defg + 1.0);
        let stress2 = t3::add(&t3::scale(&invc, lame1 * jj), &t3::scale(&t3::sub(&eye, &invc), lame2));
        let stress1 = t3::mul(&defg, &stress2);

        // dS for a displacement gradient `g`.
        let dstress2 = |g: &T3| -> T3 {
            let dc = t3::add(&t3::mul(&t3::transpose(&defg), g), &t3::mul(&t3::transpose(g), &defg));
            let sandwich = t3::mul(&invc, &t3::mul(&dc, &invc));
            let tr = t3::trace(&t3::mul(g, &invdefg), dim);
            let a = t3::scale(&invc, lame1 * det_defg * det_defg * tr);
            let b = t3::scale(&sandwich, -lame1 * jj);
            let cc = t3::scale(&sandwich, lame2);
            t3::add(&t3::add(&a, &b), &cc)
        };
        let dstress1 = |g: &T3| -> T3 { t3::add(&t3::mul(g, &stress2), &t3::mul(&defg, &dstress2(g))) };
        let dstress1_u = dstress1(&grad);

        for p in 0..n {
            for i in 0..cdim {
                let k = cdim * p + i;
                let g = trial_gradient(i, p, bas, db, axi, r);
                let ds1 = dstress1(&g);
                out.force[k] += (bas[p] * force[i] * det_defg
                    + bas[p] * inertial[i] * density
                    - t3::ddot(&g, &stress1, 3)
                    + t3::ddot(&g, &dstress1_u, 3))
                    * s;
                for q in 0..n {
                    for j in 0..cdim {
                        let v = if axi {
                            match j {
                                0 => db[q][0] * ds1[0][0] + db[q][1] * ds1[0][1] + bas[q] / r * ds1[2][2],
                                _ => db[q][0] * ds1[1][0] + db[q][1] * ds1[1][1],
                            }
                        } else {
                            db[q][0] * ds1[j][0] + db[q][1] * ds1[j][1] + db[q][2] * ds1[j][2]
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
