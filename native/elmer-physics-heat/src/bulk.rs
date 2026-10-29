use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

use crate::error::HeatError;
use crate::material::{
    check_len, check_opt, interp, Conductivity, HeatBulkConfig, HeatMaterialAtNodes, Stabilization,
};

/// Bulk element matrices of `LocalMatrixVec` (HeatSolve.F90:1006) for an
/// element with `nd = n` nodes (no bubble degrees of freedom).
///
/// `out` is reset to `n x n`: `stiff` holds diffusion plus convection (plus
/// SUPG), `force` the heat source, and `mass` (set only when
/// `cfg.transient`) the `rho*cp` mass matrix. Matrices are row-major with the
/// test function on the row, so `stiff[p*n + q]` couples test `p` to trial
/// `q`, as Fortran's `STIFF(p,q)`.
///
/// Differences from the Fortran that callers should know about:
/// * `Joule Heat`, compressibility, phase change and bubbles are not
///   implemented (see [`crate::check_supported`]).
/// * SUPG with `rho*cp = 0` at a point gives `tau = 0` where Fortran divides
///   by zero.
pub fn heat_bulk_local(
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
    if ips.is_empty() {
        return Err(HeatError::NoQuadrature(etype.0));
    }

    let mut be = BasisEval::default();
    let mut stream = vec![0.0; n];
    // The Heat Transfer Multiplier is read once, at the first integration
    // point (`Basis(1,1:nd)` in the Fortran).
    let mut multiplier = None;

    for (t, qp) in ips.iter().enumerate() {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(HeatError::Basis)?;
        let phi = &be.basis[..n];
        let dphi = &be.dbasis_dx[..n];

        let mut w = qp.weight * be.det_j;
        if cfg.axisymmetric {
            w *= phi.iter().zip(xyz).map(|(b, x)| b * x[0]).sum::<f64>();
        }

        if t == 0 {
            if let Some(m) = &mat.heat_transfer_multiplier {
                multiplier = Some(interp(phi, m));
            }
        }

        let rho = mat.density.as_deref().map_or(0.0, |r| interp(phi, r));

        // Diffusion: (k grad u, grad v).
        let mut k_tau = 0.0;
        match &mat.conductivity {
            None => {}
            Some(Conductivity::Scalar(k)) => {
                let kt = interp(phi, k);
                k_tau = kt;
                for p in 0..n {
                    for q in 0..n {
                        let mut s = 0.0;
                        for i in 0..dim {
                            s += dphi[p][i] * dphi[q][i];
                        }
                        out.stiff[p * n + q] += w * kt * s;
                    }
                }
            }
            Some(Conductivity::Diagonal(k)) => {
                let kt = interp_diag(phi, k);
                k_tau = kt[0];
                for p in 0..n {
                    for q in 0..n {
                        let mut a = 0.0;
                        for i in 0..dim {
                            a += kt[i] * dphi[p][i] * dphi[q][i];
                        }
                        out.stiff[p * n + q] += w * a;
                    }
                }
            }
            Some(Conductivity::Full(k)) => {
                let kt = interp_full(phi, k);
                k_tau = kt[0][0];
                for p in 0..n {
                    for q in 0..n {
                        let mut a = 0.0;
                        for i in 0..dim {
                            for j in 0..dim {
                                a += kt[i][j] * dphi[p][i] * dphi[q][j];
                            }
                        }
                        out.stiff[p * n + q] += w * a;
                    }
                }
            }
        }

        // rho*cp: the convection coefficient and the mass coefficient.
        let rho_cp = if convected || transient {
            let cp = mat.heat_capacity.as_deref().map_or(0.0, |c| interp(phi, c));
            cp * rho
        } else {
            0.0
        };

        // Convection: (rho cp v.grad u, v).
        let mut vel = [0.0; 3];
        if let Some(v) = &mat.convection_velocity {
            for (i, vi) in vel.iter_mut().enumerate().take(dim) {
                *vi = phi.iter().zip(v).map(|(b, vn)| b * vn[i]).sum();
            }
            for p in 0..n {
                for q in 0..n {
                    let mut vg = 0.0;
                    for i in 0..dim {
                        vg += vel[i] * dphi[q][i];
                    }
                    out.stiff[p * n + q] += vg * phi[p] * w * rho_cp;
                }
            }
        }

        // SUPG: tau (rho cp v.grad u, rho cp v.grad v) and, transient, the
        // streamline-weighted mass.
        let mut tau = 0.0;
        if let Some((h_k, m_k)) = supg {
            for p in 0..n {
                let mut s = 0.0;
                for i in 0..dim {
                    s += rho_cp * vel[i] * dphi[p][i];
                }
                stream[p] = s;
            }
            let vnorm = vel[..dim].iter().map(|v| v * v).sum::<f64>().sqrt();
            if vnorm > 0.0 && k_tau != 0.0 && rho_cp != 0.0 {
                let tmin = (m_k * h_k * rho_cp * vnorm / (2.0 * k_tau.abs())).min(1.0);
                tau = h_k * tmin / (2.0 * rho_cp * vnorm);
            }
            for p in 0..n {
                for q in 0..n {
                    out.stiff[p * n + q] += stream[p] * stream[q] * tau * w;
                }
            }
            if transient {
                if let Some(mass) = out.mass.as_mut() {
                    for p in 0..n {
                        for q in 0..n {
                            mass[p * n + q] += stream[p] * phi[q] * tau * rho_cp * w;
                        }
                    }
                }
            }
        }

        // Time derivative: (rho cp dT/dt, v).
        if transient {
            if let Some(mass) = out.mass.as_mut() {
                for p in 0..n {
                    for q in 0..n {
                        mass[p * n + q] += phi[p] * phi[q] * w * rho_cp;
                    }
                }
            }
        }

        // Source: Volumetric Heat Source wins; otherwise Heat Source * rho.
        let src = if let Some(s) = &mat.volumetric_heat_source {
            Some(interp(phi, s))
        } else {
            mat.heat_source.as_deref().map(|s| interp(phi, s) * rho)
        };
        if let Some(src) = src {
            for p in 0..n {
                out.force[p] += phi[p] * w * src;
            }
            if supg.is_some() {
                for p in 0..n {
                    out.force[p] += stream[p] * tau * w * src;
                }
            }
        }
    }

    if let Some(m) = multiplier {
        out.stiff.iter_mut().for_each(|x| *x *= m);
        out.force.iter_mut().for_each(|x| *x *= m);
        if let Some(mass) = out.mass.as_mut() {
            mass.iter_mut().for_each(|x| *x *= m);
        }
    }
    Ok(())
}

fn reset(out: &mut LocalSystem, n: usize, transient: bool) {
    out.n = n;
    out.stiff.clear();
    out.stiff.resize(n * n, 0.0);
    out.force.clear();
    out.force.resize(n, 0.0);
    out.damp = None;
    out.mass = if transient {
        let mut m = out.mass.take().unwrap_or_default();
        m.clear();
        m.resize(n * n, 0.0);
        Some(m)
    } else {
        None
    };
}

fn interp_diag(phi: &[f64], k: &[[f64; 3]]) -> [f64; 3] {
    let mut r = [0.0; 3];
    for (b, kn) in phi.iter().zip(k) {
        for i in 0..3 {
            r[i] += b * kn[i];
        }
    }
    r
}

fn interp_full(phi: &[f64], k: &[[[f64; 3]; 3]]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for (b, kn) in phi.iter().zip(k) {
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] += b * kn[i][j];
            }
        }
    }
    r
}

fn check_nodal_inputs(n: usize, mat: &HeatMaterialAtNodes) -> Result<(), HeatError> {
    if let Some(k) = &mat.conductivity {
        check_len("heat conductivity", n, k.len())?;
    }
    check_opt("density", n, &mat.density)?;
    check_opt("heat capacity", n, &mat.heat_capacity)?;
    check_opt("heat source", n, &mat.heat_source)?;
    check_opt("volumetric heat source", n, &mat.volumetric_heat_source)?;
    check_opt("heat transfer multiplier", n, &mat.heat_transfer_multiplier)?;
    if let Some(v) = &mat.convection_velocity {
        check_len("convection velocity", n, v.len())?;
    }
    Ok(())
}
