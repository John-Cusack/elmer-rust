//! `BoundaryLocalMatrix` (ThermoElectricSolver.F90): heat flux and heat
//! transfer on the temperature, electric flux, anisotropic current density
//! and the anisotropic Robin condition on the potential.

#![allow(clippy::needless_range_loop)]

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};
use elmer_physics_flow::boundary::{element_center, outward_normal};

use crate::bulk::{reset, Tensor};
use crate::error::ThermoElectricError;

/// What `BoundaryLocalMatrix` reads for one boundary element, at its nodes.
pub struct ThermoElectricBoundaryInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    /// `Heat Flux`.
    pub heat_flux: &'a [f64],
    /// `External Temperature`.
    pub ext_temperature: &'a [f64],
    /// `Heat Transfer Coefficient`.
    pub heat_transfer: &'a [f64],
    /// `Electric Flux`.
    pub electric_flux: &'a [f64],
    /// `Electric Flux Vector`, when present.
    pub electric_flux_vector: Option<&'a [[f64; 3]]>,
    /// `External Potential`.
    pub ext_potential: &'a [f64],
    /// `Electric Transfer Coefficient`, when present.
    pub transfer_tensor: Option<&'a [Tensor]>,
    /// A point inside the parent element, orienting the normal outward.
    pub parent_center: Option<[f64; 3]>,
    pub np: Option<usize>,
}

fn check_len(what: &'static str, expected: usize, got: usize) -> Result<(), ThermoElectricError> {
    if got == expected {
        Ok(())
    } else {
        Err(ThermoElectricError::NodalLength { what, expected, got })
    }
}

fn interp(v: &[f64], bas: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, b) in v.iter().zip(bas) {
        s += x * b;
    }
    s
}

/// The local system of one boundary element (no mass matrix).
pub fn thermoelectric_boundary_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &ThermoElectricBoundaryInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), ThermoElectricError> {
    let n = inp.xyz.len();
    for (what, len) in [
        ("heat_flux", inp.heat_flux.len()),
        ("ext_temperature", inp.ext_temperature.len()),
        ("heat_transfer", inp.heat_transfer.len()),
        ("electric_flux", inp.electric_flux.len()),
        ("ext_potential", inp.ext_potential.len()),
    ] {
        check_len(what, n, len)?;
    }
    if let Some(v) = inp.electric_flux_vector {
        check_len("electric_flux_vector", n, v.len())?;
    }
    if let Some(t) = inp.transfer_tensor {
        check_len("transfer_tensor", n, t.len())?;
    }
    reset(out, 2 * n, false);
    let dim2 = 2 * n;
    let mut ev = BasisEval::default();

    for gp in quad.gauss_points(inp.etype, inp.np) {
        basis
            .eval(inp.etype, inp.xyz, gp.uvw, &mut ev)
            .map_err(ThermoElectricError::Basis)?;
        check_len("basis", n, ev.basis.len())?;
        let bas = &ev.basis;
        let s = gp.weight * ev.det_j;

        let hf = interp(inp.heat_flux, bas);
        let et = interp(inp.ext_temperature, bas);
        let ht = interp(inp.heat_transfer, bas);
        let mut ef = interp(inp.electric_flux, bas);
        let ep = interp(inp.ext_potential, bas);

        let mut hte = 0.0;
        if inp.electric_flux_vector.is_some() || inp.transfer_tensor.is_some() {
            let center = inp.parent_center.unwrap_or_else(|| element_center(inp.xyz));
            let normal = outward_normal(inp.etype, inp.xyz, gp.uvw, center)
                .map_err(|e| ThermoElectricError::Unsupported(e.to_string()))?;
            if let Some(efv_n) = inp.electric_flux_vector {
                let efv: [f64; 3] = std::array::from_fn(|i| {
                    let mut s = 0.0;
                    for k in 0..n {
                        s += efv_n[k][i] * bas[k];
                    }
                    s
                });
                // current density vector dotted with the normal
                ef += efv[0] * normal[0] + efv[1] * normal[1] + efv[2] * normal[2];
            }
            if let Some(etc_n) = inp.transfer_tensor {
                let etc: Tensor = std::array::from_fn(|i| {
                    std::array::from_fn(|j| {
                        let mut s = 0.0;
                        for k in 0..n {
                            s += etc_n[k][i][j] * bas[k];
                        }
                        s
                    })
                });
                // normal-normal component of the conductance tensor
                let mut sum = 0.0;
                for i in 0..3 {
                    let row = etc[i][0] * normal[0] + etc[i][1] * normal[1] + etc[i][2] * normal[2];
                    sum += normal[i] * row;
                }
                hte = sum;
            }
        }

        for p in 0..n {
            let i = 2 * p;
            for q in 0..n {
                let j = 2 * q;
                out.stiff[i * dim2 + j] += s * ht * bas[q] * bas[p];
                out.stiff[(i + 1) * dim2 + j + 1] += s * hte * bas[q] * bas[p];
            }
            out.force[i] += s * (hf + ht * et) * bas[p];
            out.force[i + 1] += s * (ef + hte * ep) * bas[p];
        }
    }
    Ok(())
}
