//! Bulk element system of `ElasticSolve.F90` `NeoHookeanLocalMatrix` (`:4863`),
//! the compressible neo-Hookean material with the second Piola-Kirchhoff stress
//! `S = Lame1/2 (J^2 - 1) C^-1 + Lame2 (I - C^-1)`, linearised exactly.
//!
//! The assembly has the same Newton form as [`crate::bulk`]. Temperature,
//! `Stress Pressure`, `Stress Load` and `Strain Load` play no part in it; the
//! solver refuses the keywords the Fortran refuses. `Mixed Formulation` adds the
//! pressure as an unknown, last at each node.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, QuadPoint};
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
    /// The displacement at the `ntot` basis functions (zero at the bubbles).
    pub disp: &'a [[f64; 3]],
    pub density: &'a [f64],
    pub damping: Option<&'a [f64]>,
    pub load: &'a [[f64; 3]],
    pub inertial: Option<&'a [[f64; 3]]>,
    pub want_mass: bool,
    /// `ntot`: the basis functions, the nodes first (`n + nb` of a p element
    /// whose bubbles stay out of the global system).
    pub ntot: usize,
    /// `Mixed Formulation`: the pressure unknown at the nodes. The local
    /// system then has `cdim + 1` dofs per basis function, the pressure last.
    pub mixed: Option<&'a [f64]>,
}

/// Bulk element matrix and force over `points`. Overwrites `out`.
pub fn neo_hookean_local(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    inp: &NeoHookeInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), ElasticError> {
    let n = inp.etype.n_nodes();
    let ntot = inp.ntot;
    let cdim = inp.cdim;
    let axi = inp.coord == CoordSystem::AxiSymmetric;
    if cdim != 2 && cdim != 3 {
        return Err(ElasticError::Unsupported(format!("coordinate system dimension {cdim}")));
    }
    if axi && cdim != 2 {
        return Err(ElasticError::Unsupported("axisymmetric elasticity needs a 2D mesh".into()));
    }
    if [inp.xyz.len(), inp.young.len(), inp.poisson.len(), inp.density.len(), inp.load.len()]
        .iter()
        .any(|&l| l != n)
        || inp.disp.len() != ntot
        || ntot < n
        || inp.mixed.is_some_and(|p| p.len() != n)
    {
        return Err(ElasticError::Input(format!("element {} needs {n} nodal values", inp.etype.0)));
    }
    let mixed = inp.mixed.is_some();
    let dim = if axi { 3 } else { cdim };
    // DOFs: the displacement, and the pressure last in the mixed form.
    let dofs = if mixed { cdim + 1 } else { cdim };
    let ndof = ntot * dofs;
    *out = LocalSystem::zeros(ndof);
    if inp.want_mass {
        out.mass = Some(vec![0.0; ndof * ndof]);
        if inp.damping.is_some() {
            out.damp = Some(vec![0.0; ndof * ndof]);
        }
    }

    // ElasticSolve.F90:4922-4945: the mixed form has no Lame1; its pressure
    // constraint is weighted by (1 + nu)(1 - 2 nu) / (E nu).
    let mut pressure_par_n = Vec::new();
    if mixed {
        if inp.plane_stress {
            return Err(ElasticError::Input("Mixed formulation does not support plane stress".into()));
        }
        if inp.poisson.iter().all(|nu| nu.abs() < f64::EPSILON) {
            return Err(ElasticError::Input("Mixed formulation with the zero Poisson ratio is not allowed".into()));
        }
        pressure_par_n = (0..n)
            .map(|k| {
                let (e, nu) = (inp.young[k], inp.poisson[k]);
                (1.0 + nu) * (1.0 - 2.0 * nu) / (e * nu)
            })
            .collect();
    }
    let (lame1_n, lame2_n): (Vec<f64>, Vec<f64>) = (0..n)
        .map(|k| {
            let (e, nu) = (inp.young[k], inp.poisson[k]);
            let l1 = if mixed {
                0.0
            } else if inp.plane_stress && !axi {
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

    for qp in points {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be).map_err(ElasticError::Basis)?;
        if be.basis.len() < ntot {
            return Err(ElasticError::Input(format!("the basis has {} functions, not {ntot}", be.basis.len())));
        }
        let bas = &be.basis[..ntot];
        let db = &be.dbasis_dx[..ntot];
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
        let mut lame2: f64 = (0..n).map(|k| lame2_n[k] * bas[k]).sum();
        // `Pressure = SUM(LocalDisplacement(DOFs,1:n) * Basis(1:n))`.
        let pressure: f64 = inp.mixed.map_or(0.0, |pn| (0..n).map(|k| pn[k] * bas[k]).sum());
        if mixed {
            lame2 += pressure;
        }
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

        for p in 0..ntot {
            for i in 0..cdim {
                let k = dofs * p + i;
                let g = trial_gradient(i, p, bas, db, axi, r);
                let ds1 = dstress1(&g);
                out.force[k] += (bas[p] * force[i] * det_defg
                    + bas[p] * inertial[i] * density
                    - t3::ddot(&g, &stress1, 3)
                    + t3::ddot(&g, &dstress1_u, 3))
                    * s;
                for q in 0..ntot {
                    for j in 0..cdim {
                        let v = if axi {
                            match j {
                                0 => db[q][0] * ds1[0][0] + db[q][1] * ds1[0][1] + bas[q] / r * ds1[2][2],
                                _ => db[q][0] * ds1[1][0] + db[q][1] * ds1[1][1],
                            }
                        } else {
                            db[q][0] * ds1[j][0] + db[q][1] * ds1[j][1] + db[q][2] * ds1[j][2]
                        };
                        out.stiff[k * ndof + dofs * q + j] += v * s;
                    }
                }
            }
        }

        if inp.want_mass {
            for p in 0..ntot {
                for q in 0..ntot {
                    for i in 0..cdim {
                        let idx = (dofs * p + i) * ndof + dofs * q + i;
                        let mm = out.mass.as_mut().unwrap();
                        mm[idx] += bas[p] * bas[q] * density * s;
                        if let Some(dm) = out.damp.as_mut() {
                            dm[idx] += bas[p] * bas[q] * damping * s;
                        }
                    }
                }
            }
        }

        if mixed {
            // ElasticSolve.F90:5076-5155: the pressure constraint and the
            // pressure's part of the displacement equations.
            let pressure_par: f64 = (0..n).map(|k| pressure_par_n[k] * bas[k]).sum();
            let gm = t3::sub(&defg, &eye);
            let j2 = det_defg * det_defg;
            let tr = t3::trace(&t3::mul(&gm, &invdefg), dim);
            for p in 0..n {
                let row = dofs * p + cdim;
                out.force[row] = out.force[row] - 0.5 * (j2 - 1.0) * bas[p] * s + j2 * tr * bas[p] * s;
                for q in 0..ntot {
                    for i in 0..cdim {
                        let v = if axi {
                            match i {
                                0 => j2 * (db[q][0] * invdefg[0][0] + db[q][1] * invdefg[1][0] + bas[q] / r * invdefg[2][2]),
                                _ => j2 * (db[q][0] * invdefg[0][1] + db[q][1] * invdefg[1][1]),
                            }
                        } else {
                            let mut sum = 0.0;
                            for a in 0..cdim {
                                sum += db[q][a] * invdefg[a][i];
                            }
                            j2 * sum
                        };
                        out.stiff[row * ndof + dofs * q + i] += v * bas[p] * s;
                    }
                    if q < n {
                        out.stiff[row * ndof + dofs * q + cdim] += pressure_par * bas[p] * bas[q] * s;
                    }
                }
            }
            let invdefg_t = t3::transpose(&invdefg);
            for p in 0..ntot {
                for i in 0..cdim {
                    let g = trial_gradient(i, p, bas, db, axi, r);
                    let k = dofs * p + i;
                    let dd = t3::ddot(&invdefg_t, &g, dim);
                    out.force[k] = out.force[k] + pressure * db[p][i] * s - pressure * dd * s;
                    if axi && i == 0 {
                        out.force[k] += pressure * bas[p] / r * s;
                    }
                    for q in 0..ntot {
                        for j in 0..cdim {
                            let v = if axi {
                                match j {
                                    0 => db[q][0] * g[0][0] + db[q][1] * g[0][1] + bas[q] / r * g[2][2],
                                    _ => db[q][0] * g[1][0] + db[q][1] * g[1][1],
                                }
                            } else {
                                db[q][0] * g[j][0] + db[q][1] * g[j][1] + db[q][2] * g[j][2]
                            };
                            out.stiff[k * ndof + dofs * q + j] -= pressure * v * s;
                        }
                        if q < n {
                            out.stiff[k * ndof + dofs * q + cdim] -= bas[q] * dd * s;
                        }
                    }
                }
                // `+ Basis(p) * Force(DOFs) * s`: no `Stress Volume Source`.
                out.force[dofs * p + cdim] += bas[p] * 0.0 * s;
            }
        }
    }
    if mixed {
        // ElasticSolve.F90:5301-5309: the pressure only at the nodes.
        for p in n..ntot {
            let i = dofs * p + cdim;
            out.force[i] = 0.0;
            for j in 0..ndof {
                out.stiff[i * ndof + j] = 0.0;
                out.stiff[j * ndof + i] = 0.0;
            }
            out.stiff[i * ndof + i] = 1.0;
        }
    }
    Ok(())
}
