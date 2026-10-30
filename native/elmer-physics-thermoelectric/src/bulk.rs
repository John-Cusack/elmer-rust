//! `LocalMatrix` (ThermoElectricSolver.F90): the element mass, stiffness and
//! force of the coupled temperature (field 1) and electric potential
//! (field 2) system.
//!
//! Dofs are node-major: local row `2 * p + k` is field `k` at basis function
//! `p`, as in Elmer. The number of basis functions `nd` is that of the
//! `ElementBasis` handed in (`solution.len()`): the `n` element nodes first,
//! then any interior bubbles, whose condensation is [`crate::condense`]'s.
//!
//! Ported: tensor `Heat Conductivity`, `Seebeck Coefficient`, `Electric
//! Conductivity` and `Electric Permittivity` (given per node), the lagged
//! Picard form and its Newton linearization of the Peltier and Joule terms,
//! and the SUPG stabilization of the Peltier drift.

// A structural transliteration of Fortran loops over `1:3` that index several
// arrays by the same component number.
#![allow(clippy::needless_range_loop)]

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

use crate::error::ThermoElectricError;

/// A material tensor, `T[i][j]`.
pub type Tensor = [[f64; 3]; 3];

/// `Element % StabilizationMK` and `Element % hK` (`MeshStabParams`).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Stabilization {
    pub mk: f64,
    pub hk: f64,
}

/// What `LocalMatrix` reads for one element, at its `n` nodes.
pub struct ThermoElectricBulkInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    /// `Heat Conductivity`.
    pub hcond: &'a [Tensor],
    /// `Seebeck Coefficient`.
    pub alpha: &'a [Tensor],
    /// `Electric Conductivity`.
    pub sigma: &'a [Tensor],
    /// `Electric Permittivity` (before the `Permittivity of Vacuum` factor).
    pub epsil: &'a [Tensor],
    /// `Permittivity of Vacuum` when `Relative Permittivity` is present, else 1.
    pub epsil0: f64,
    pub rho: &'a [f64],
    pub cp: &'a [f64],
    /// `Heat Source` times `Density`, or `Volumetric Heat Source`.
    pub heat_source: &'a [f64],
    /// `Current Source`.
    pub current_source: &'a [f64],
    /// Temperature and potential of the previous iterate at the `nd` basis
    /// functions (`GetVectorLocalSolution`); bubble values last.
    pub solution: &'a [[f64; 2]],
    /// Newton linearization of the Peltier and Joule terms.
    pub newton: bool,
    /// `Stabilize`: the SUPG form of the Peltier drift.
    pub stabilize: Option<Stabilization>,
    /// Also return the mass matrix (a transient run).
    pub want_mass: bool,
    /// Gauss point count, `None` for the element default.
    pub np: Option<usize>,
}

fn check_len(what: &'static str, expected: usize, got: usize) -> Result<(), ThermoElectricError> {
    if got == expected {
        Ok(())
    } else {
        Err(ThermoElectricError::NodalLength { what, expected, got })
    }
}

/// `MATMUL(t, v)`.
fn mat_vec(t: &Tensor, v: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| t[i][0] * v[0] + t[i][1] * v[1] + t[i][2] * v[2])
}

/// `MATMUL(a, b)`.
fn mat_mat(a: &Tensor, b: &Tensor) -> Tensor {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j]))
}

/// `SUM(a * b)` of two vectors.
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `SUM(t_n(i,j,:) * Basis(1:n))`.
fn interpolate(t: &[Tensor], basis: &[f64]) -> Tensor {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let mut s = 0.0;
            for (tk, b) in t.iter().zip(basis) {
                s += tk[i][j] * b;
            }
            s
        })
    })
}

fn interpolate_scalar(v: &[f64], basis: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, b) in v.iter().zip(basis) {
        s += x * b;
    }
    s
}

/// Sizes `out` to `n` dofs, zeroed, with a mass matrix when asked.
pub(crate) fn reset(out: &mut LocalSystem, n: usize, want_mass: bool) {
    out.n = n;
    out.stiff.clear();
    out.stiff.resize(n * n, 0.0);
    out.force.clear();
    out.force.resize(n, 0.0);
    out.damp = None;
    out.mass = want_mass.then(|| vec![0.0; n * n]);
}

/// The local system of one bulk element.
pub fn thermoelectric_bulk_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &ThermoElectricBulkInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), ThermoElectricError> {
    let n = inp.xyz.len();
    let nd = inp.solution.len();
    for (what, len) in [
        ("hcond", inp.hcond.len()),
        ("alpha", inp.alpha.len()),
        ("sigma", inp.sigma.len()),
        ("epsil", inp.epsil.len()),
        ("rho", inp.rho.len()),
        ("cp", inp.cp.len()),
        ("heat_source", inp.heat_source.len()),
        ("current_source", inp.current_source.len()),
    ] {
        check_len(what, n, len)?;
    }
    if nd < n {
        return Err(ThermoElectricError::NodalLength {
            what: "solution",
            expected: n,
            got: nd,
        });
    }
    let dim2 = 2 * nd;
    reset(out, dim2, inp.want_mass);
    let mut ev = BasisEval::default();

    for gp in quad.gauss_points(inp.etype, inp.np) {
        basis
            .eval(inp.etype, inp.xyz, gp.uvw, &mut ev)
            .map_err(ThermoElectricError::Basis)?;
        // A p basis may evaluate more functions than the element has dofs
        // (a bubble count below the degree's): the first `nd` are its own.
        if ev.basis.len() < nd {
            return Err(ThermoElectricError::NodalLength { what: "basis", expected: nd, got: ev.basis.len() });
        }
        let bas = &ev.basis[..nd];
        let dbas = &ev.dbasis_dx[..nd];
        let b_n = &bas[..n];

        let hcond = interpolate(inp.hcond, b_n);
        let alpha = interpolate(inp.alpha, b_n);
        let sigma = interpolate(inp.sigma, b_n);
        let mut epsil = interpolate(inp.epsil, b_n);
        for row in epsil.iter_mut() {
            for x in row.iter_mut() {
                *x *= inp.epsil0;
            }
        }
        let rho = interpolate_scalar(inp.rho, b_n);
        let c_p = interpolate_scalar(inp.cp, b_n);

        let ic = interpolate_scalar(inp.current_source, b_n); // input current
        let hs = interpolate_scalar(inp.heat_source, b_n); // input heat source

        // Electric field and temperature gradient from the previous solution.
        let mut ef = [0.0; 3];
        let mut tg = [0.0; 3];
        for c in 0..3 {
            let (mut pot, mut temp) = (0.0, 0.0);
            for k in 0..nd {
                pot += inp.solution[k][1] * dbas[k][c];
                temp += inp.solution[k][0] * dbas[k][c];
            }
            ef[c] = -pot;
            tg[c] = temp;
        }
        let a_tg = mat_vec(&alpha, tg);
        let cd = mat_vec(&sigma, std::array::from_fn(|i| ef[i] - a_tg[i])); // current density
        let jh = dot(cd, ef); // joule heating

        let sb = mat_mat(&sigma, &alpha);
        let mut temp = 0.0;
        for k in 0..nd {
            temp += bas[k] * inp.solution[k][0];
        }
        let pelt: Tensor = std::array::from_fn(|i| std::array::from_fn(|j| temp * alpha[i][j])); // Peltier coefficient
        let pelt_sb = mat_mat(&pelt, &sb);
        let pelt_si = mat_mat(&pelt, &sigma);

        // SUPG stabilization of the temperature equation: the Peltier term
        // linearizes into a drift alpha*J acting on grad(T).
        let mut tau = 0.0;
        let mut stream = vec![0.0; nd];
        if let Some(st) = inp.stabilize {
            let v_eff = mat_vec(&alpha, cd);
            let vnorm = dot(v_eff, v_eff).sqrt();
            let kappa_eff = if vnorm > 0.0 {
                dot(v_eff, mat_vec(&hcond, v_eff)) / (vnorm * vnorm)
            } else {
                0.0
            };
            if vnorm > 0.0 && kappa_eff != 0.0 {
                tau = f64::min(1.0, st.mk * st.hk * vnorm / (2.0 * kappa_eff.abs()));
                tau = st.hk * tau / (2.0 * vnorm);
            }
            for p in 0..nd {
                stream[p] = dot(v_eff, dbas[p]);
            }
        }

        let s = gp.weight * ev.det_j;
        for p in 0..nd {
            let i = 2 * p;
            for q in 0..nd {
                let j = 2 * q;
                let at = |r: usize, c: usize| (i + r) * dim2 + j + c;
                if let Some(m) = out.mass.as_mut() {
                    // thermal damping
                    m[at(0, 0)] += s * rho * c_p * bas[q] * bas[p];
                    m[at(0, 0)] += s * tau * rho * c_p * bas[q] * stream[p]; // SUPG
                    // displacement current
                    m[at(1, 1)] += s * dot(mat_vec(&epsil, dbas[q]), dbas[p]);
                }
                let a = &mut out.stiff;
                // thermal diffusion
                a[at(0, 0)] += s * dot(mat_vec(&hcond, dbas[q]), dbas[p]);
                // SUPG streamline-diffusion consistency term for the Peltier drift
                a[at(0, 0)] += s * tau * stream[q] * stream[p];

                if inp.newton {
                    // Newton linearization of div(pelt*J) = div(T*alpha*J):
                    // the Peltier coefficient implicitly,
                    a[at(0, 0)] -= s * bas[q] * dot(mat_vec(&alpha, cd), dbas[p]);
                    // the temperature gradient part of div(pelt_0*J),
                    a[at(0, 0)] -= s * dot(mat_vec(&pelt_sb, dbas[q]), dbas[p]);
                    // the electric field part of div(pelt_0*J),
                    a[at(0, 1)] -= s * dot(mat_vec(&pelt_si, dbas[q]), dbas[p]);
                    // and of the Joule heating -(J,E) ~ -(J_0,E) + (J,E_0) - (J_0,E_0):
                    // -(J_0,E)
                    a[at(0, 1)] += s * dot(cd, dbas[q]) * bas[p];
                    // the temperature gradient part of (J,E_0)
                    a[at(0, 0)] -= s * dot(mat_vec(&sb, dbas[q]), ef) * bas[p];
                    // the electric field part of (J,E_0)
                    a[at(0, 1)] -= s * dot(mat_vec(&sigma, dbas[q]), ef) * bas[p];
                }

                // Seebeck
                a[at(1, 0)] += s * dot(mat_vec(&sb, dbas[q]), dbas[p]);
                // electric diffusion
                a[at(1, 1)] += s * dot(mat_vec(&sigma, dbas[q]), dbas[p]);
            }

            // heat load
            out.force[i] += s * (hs + jh) * bas[p];
            out.force[i] += s * tau * (hs + jh) * stream[p]; // SUPG
            // div(pelt_0*J_0)
            out.force[i] += s * dot(mat_vec(&pelt, cd), dbas[p]);
            // electric load
            out.force[i + 1] += s * ic * bas[p];
        }
    }
    Ok(())
}
