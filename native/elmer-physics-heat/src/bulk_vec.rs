// Index loops spell out the Fortran loop order, which fixes the rounding.
#![allow(clippy::needless_range_loop)]

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

use crate::bulk::{check_nodal_inputs, interp_diag, interp_full, reset};
use crate::error::HeatError;
use crate::linear_forms::{grad_u_dot_grad_u, grad_u_dot_u, u_dot_f, u_dot_v};
use crate::material::{
    CompressModel, Conductivity, HeatBulkConfig, HeatMaterialAtNodes, Stabilization,
};

/// [`crate::heat_bulk_local`] in the summation order of `LocalMatrixVec`
/// (HeatSolve.F90:1006), the kernel of the HeatSolveVec path.
///
/// The basis is evaluated at every integration point first, as
/// `ElementInfoVec` does. Each term is then added over all points at once
/// with the forms of `LinearForms.F90` (see `linear_forms`), in the order
/// the Fortran adds them:
/// * stiffness: diffusion, convection, SUPG;
/// * force: pressure work, then source;
/// * mass: SUPG, then `rho*cp`.
///
/// A keyword is taken at the points as `ListGetElementRealVec` does. When its
/// values at the element nodes all agree, that value is used exactly, as
/// the Fortran returns a constant keyword. Otherwise it is interpolated with
/// the basis. A dependent keyword whose nodal values happen to agree is
/// interpolated by the Fortran and so can differ from this in the last bit.
/// Solution fields (the temperature, the flow velocity and pressure, the
/// mesh velocity) are always interpolated.
///
/// Inputs, outputs and unsupported features are as for
/// [`crate::heat_bulk_local`].
pub fn heat_bulk_local_vec(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    etype: ElementType,
    xyz: &[[f64; 3]],
    mat: &HeatMaterialAtNodes,
    cfg: &HeatBulkConfig,
    out: &mut LocalSystem,
) -> Result<(), HeatError> {
    let n = xyz.len();
    check_nodal_inputs(n, mat)?;
    if !(1..=3).contains(&cfg.dim) {
        return Err(HeatError::Invalid(format!("dim {} not in 1..=3", cfg.dim)));
    }
    let dim = cfg.dim;
    let convected = mat.convection_velocity.is_some();
    let transient = cfg.transient;
    let supg = match cfg.stabilization {
        Stabilization::Supg { h_k, m_k } if convected && mat.conductivity.is_some() => {
            Some((h_k, m_k))
        }
        _ => None,
    };

    reset(out, n, transient);

    let ips = quad.gauss_points(etype, cfg.quad_points);
    let ngp = ips.len();
    if ngp == 0 {
        return Err(HeatError::NoQuadrature(etype.0));
    }

    // ElementInfoVec: Basis(t,p), dBasisdx(t,p,k), and the weights
    // DetJVec = s*detJ (times r when axisymmetric).
    let mut bas = vec![0.0; ngp * n];
    let mut dbas = vec![0.0; 3 * ngp * n];
    let mut w = vec![0.0; ngp];
    let mut be = BasisEval::default();
    for (t, qp) in ips.iter().enumerate() {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(HeatError::Basis)?;
        bas[t * n..(t + 1) * n].copy_from_slice(&be.basis[..n]);
        for p in 0..n {
            for k in 0..dim {
                dbas[(k * ngp + t) * n + p] = be.dbasis_dx[p][k];
            }
        }
        w[t] = qp.weight * be.det_j;
    }
    if cfg.axisymmetric {
        for t in 0..ngp {
            w[t] *= field_at(&bas, n, t, xyz.iter().map(|x| x[0]));
        }
    }
    let grad = |t: usize, p: usize, k: usize| dbas[(k * ngp + t) * n + p];

    // ListGetElementRealVec of a keyword, ListGetElementScalarSolutionVec of
    // a field.
    let keyword = |v: &[f64]| -> Vec<f64> {
        if v.iter().all(|&x| x == v[0]) {
            vec![v[0]; ngp]
        } else {
            (0..ngp)
                .map(|t| field_at(&bas, n, t, v.iter().copied()))
                .collect()
        }
    };
    let solution = |v: &[f64]| -> Vec<f64> {
        (0..ngp)
            .map(|t| field_at(&bas, n, t, v.iter().copied()))
            .collect()
    };
    // ListGetElementReal at the first point.
    let first = |v: &[f64]| -> f64 { keyword(v)[0] };
    let zero = || vec![0.0; ngp];

    let mut rho = mat.density.as_deref().map_or_else(zero, keyword);

    // Diffusion: STIFF += (kappa grad u, grad v).
    let mut cond_tau: Option<Vec<f64>> = None;
    match &mat.conductivity {
        None => {}
        Some(Conductivity::Scalar(k)) => {
            let k = keyword(k);
            grad_u_dot_grad_u(ngp, n, dim, &dbas, &w, &mut out.stiff, &k);
            cond_tau = Some(k);
        }
        Some(Conductivity::Diagonal(k)) => {
            let kt = if k.iter().all(|x| *x == k[0]) {
                vec![k[0]; ngp]
            } else {
                (0..ngp)
                    .map(|t| interp_diag(&bas[t * n..(t + 1) * n], k))
                    .collect()
            };
            for t in 0..ngp {
                for q in 0..n {
                    for p in 0..n {
                        let mut a = 0.0;
                        for i in 0..dim {
                            a += kt[t][i] * grad(t, p, i) * grad(t, q, i);
                        }
                        out.stiff[p * n + q] += w[t] * a;
                    }
                }
            }
            cond_tau = Some(kt.iter().map(|k| k[0]).collect());
        }
        Some(Conductivity::Full(k)) => {
            let kt = if k.iter().all(|x| *x == k[0]) {
                vec![k[0]; ngp]
            } else {
                (0..ngp)
                    .map(|t| interp_full(&bas[t * n..(t + 1) * n], k))
                    .collect()
            };
            for t in 0..ngp {
                for q in 0..n {
                    for p in 0..n {
                        let mut a = 0.0;
                        for i in 0..dim {
                            for j in 0..dim {
                                a += kt[t][i][j] * grad(t, p, i) * grad(t, q, j);
                            }
                        }
                        out.stiff[p * n + q] += w[t] * a;
                    }
                }
            }
            cond_tau = Some(kt.iter().map(|k| k[0][0]).collect());
        }
    }

    // rho*cp (TmpVec): the convection and mass coefficient.
    let cp = if convected || transient {
        mat.heat_capacity.as_deref().map_or_else(zero, keyword)
    } else {
        zero()
    };
    let mut rho_cp: Vec<f64> = (0..ngp).map(|t| cp[t] * rho[t]).collect();

    // Compressibility Model: the density, and the pressure-work coefficient.
    let mut pressure_work = None;
    if let Some(c) = &mat.compressibility {
        let temp = solution(&c.temperature);
        match &c.model {
            CompressModel::Thermal {
                heat_expansion,
                reference_temperature,
            } => {
                let t_ref = reference_temperature.as_deref().map_or_else(zero, keyword);
                let beta = heat_expansion.as_deref().map_or_else(zero, keyword);
                for t in 0..ngp {
                    rho[t] *= 1.0 - beta[t] * (temp[t] - t_ref[t]);
                }
                if convected || transient {
                    rho_cp = (0..ngp).map(|t| cp[t] * rho[t]).collect();
                }
                if let Some(p) = c.pressure.as_deref() {
                    let pcoeff: Vec<f64> = (0..ngp)
                        .map(|t| temp[t] * beta[t] / (1.0 - beta[t] * (temp[t] - t_ref[t])))
                        .collect();
                    pressure_work = Some((pcoeff, p));
                }
            }
            CompressModel::PerfectGas {
                specific_heat_ratio,
                reference_pressure,
            } => {
                let Some(p) = c.pressure.as_deref() else {
                    return Err(HeatError::Invalid(
                        "Compressibility Model \"Perfect Gas\" requires computed convection".into(),
                    ));
                };
                let gamma = specific_heat_ratio.as_deref().map_or(5.0 / 3.0, first);
                let p0 = reference_pressure.as_deref().map_or(0.0, first);
                let press = solution(p);
                for t in 0..ngp {
                    rho[t] = (press[t] + p0) / ((gamma - 1.0) / gamma * cp[t] * temp[t]);
                }
                rho_cp = (0..ngp).map(|t| cp[t] * rho[t]).collect();
                pressure_work = Some((vec![1.0; ngp], p));
            }
        }
    }

    // Convection: STIFF += (rho cp v.grad u, v).
    let mut vel = vec![[0.0; 3]; ngp];
    let mut stream = Vec::new();
    let mut tau = Vec::new();
    if let Some(v) = &mat.convection_velocity {
        for i in 0..dim {
            let vi: Vec<f64> = v.iter().map(|x| x[i]).collect();
            let at = if mat.convection_computed {
                solution(&vi)
            } else {
                keyword(&vi)
            };
            for t in 0..ngp {
                vel[t][i] = at[t];
            }
        }
        if let Some(m) = &mat.mesh_velocity {
            for i in 0..dim {
                let mi: Vec<f64> = m.iter().map(|x| x[i]).collect();
                let at = solution(&mi);
                for t in 0..ngp {
                    vel[t][i] -= at[t];
                }
            }
        }
        grad_u_dot_u(ngp, n, dim, &dbas, &bas, &w, &mut out.stiff, &rho_cp, &vel);

        // Pressure work: FORCE += Pcoeff*(dp/dt + v.grad p).
        if let Some((pcoeff, pn)) = &pressure_work {
            if let Some(c) = &mat.compressibility {
                if let Some((prev, dt)) = &c.prev_pressure {
                    let (press, prev) = (solution(pn), solution(prev));
                    let f: Vec<f64> = (0..ngp)
                        .map(|t| pcoeff[t] * (press[t] - prev[t]) / dt)
                        .collect();
                    u_dot_f(ngp, n, &bas, &w, &f, &mut out.force);
                }
            }
            let f: Vec<f64> = (0..ngp)
                .map(|t| {
                    let mut vg = 0.0;
                    for i in 0..dim {
                        let mut gp = 0.0;
                        for (p, pv) in pn.iter().enumerate() {
                            gp += pv * grad(t, p, i);
                        }
                        vg += vel[t][i] * gp;
                    }
                    pcoeff[t] * vg
                })
                .collect();
            u_dot_f(ngp, n, &bas, &w, &f, &mut out.force);
        }

        // SUPG: STIFF += tau (rho cp v.grad u, rho cp v.grad v), and the
        // streamline-weighted mass.
        if let (Some((h_k, m_k)), Some(k)) = (supg, &cond_tau) {
            stream = vec![0.0; ngp * n];
            for i in 0..dim {
                for p in 0..n {
                    for t in 0..ngp {
                        stream[t * n + p] += rho_cp[t] * vel[t][i] * grad(t, p, i);
                    }
                }
            }
            tau = (0..ngp)
                .map(|t| {
                    let vnorm = vel[t][..dim].iter().map(|v| v * v).sum::<f64>().sqrt();
                    if vnorm > 0.0 && k[t] != 0.0 && rho_cp[t] != 0.0 {
                        let tmin = (m_k * h_k * rho_cp[t] * vnorm / (2.0 * k[t].abs())).min(1.0);
                        h_k * tmin / (2.0 * rho_cp[t] * vnorm)
                    } else {
                        0.0
                    }
                })
                .collect();
            let tw: Vec<f64> = (0..ngp).map(|t| tau[t] * w[t]).collect();
            u_dot_v(ngp, n, &stream, &stream, &tw, &mut out.stiff, None);
            if let Some(mass) = out.mass.as_mut() {
                let tw: Vec<f64> = (0..ngp).map(|t| tau[t] * rho_cp[t] * w[t]).collect();
                u_dot_v(ngp, n, &stream, &bas, &tw, mass, None);
            }
        }
    }

    // Time derivative: MASS += (rho cp dT/dt, v).
    if let Some(mass) = out.mass.as_mut() {
        u_dot_v(ngp, n, &bas, &bas, &w, mass, Some(&rho_cp));
    }

    // Source: Volumetric Heat Source wins; otherwise Heat Source * rho.
    let src = match (&mat.volumetric_heat_source, &mat.heat_source) {
        (Some(s), _) => Some(keyword(s)),
        (None, Some(s)) => {
            let s = keyword(s);
            Some((0..ngp).map(|t| s[t] * rho[t]).collect::<Vec<f64>>())
        }
        (None, None) => None,
    };
    if let Some(src) = src {
        u_dot_f(ngp, n, &bas, &w, &src, &mut out.force);
        if !stream.is_empty() {
            let tw: Vec<f64> = (0..ngp).map(|t| tau[t] * w[t]).collect();
            u_dot_f(ngp, n, &stream, &tw, &src, &mut out.force);
        }
    }

    // Heat Transfer Multiplier, at the first point, on the whole system.
    if let Some(m) = mat.heat_transfer_multiplier.as_deref().map(first) {
        out.stiff.iter_mut().for_each(|x| *x *= m);
        out.force.iter_mut().for_each(|x| *x *= m);
        if let Some(mass) = out.mass.as_mut() {
            mass.iter_mut().for_each(|x| *x *= m);
        }
    }
    Ok(())
}

/// `SUM(Basis(t,1:n) * values(1:n))`, the `MATMUL` of the vector getters.
fn field_at(bas: &[f64], n: usize, t: usize, values: impl Iterator<Item = f64>) -> f64 {
    let mut s = 0.0;
    for (b, v) in bas[t * n..(t + 1) * n].iter().zip(values) {
        s += b * v;
    }
    s
}
