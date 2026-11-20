//! The element matrices of `FilmFlowSolver` (FilmFlowSolver.F90):
//! `LocalBulkMatrix`, `LocalBoundaryMatrix`, the friction law
//! `FrictionLawPraks` and `NSCondensate` (MatrixAssembly.F90), in the
//! Fortran's loop and evaluation order.
//!
//! Local systems are node-major with `mdim + 1` dofs per basis function:
//! the `mdim` velocity components, then the pressure.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, QuadPoint};

use crate::FilmFlowError;

/// `AEPS = 10 * EPSILON(1.0_dp)` (Types.F90).
const AEPS: f64 = 10.0 * f64::EPSILON;

/// `Friction Model`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Friction {
    /// Poiseuille flow resistance (`laminar`, the default).
    Laminar,
    /// `darcy`: the Darcy friction factor of [`friction_law_praks`].
    Darcy,
    /// `darcy2`: the same with the slope of the hydraulic head.
    Darcy2,
    /// `manning`.
    Manning,
    /// `manning2`.
    Manning2,
}

/// Solver-level switches of `LocalBulkMatrix`.
#[derive(Clone, Debug)]
pub struct BulkOptions {
    pub mdim: usize,
    /// `Axi Symmetric`: `geomc = 2` and the radius in the laminar law.
    pub csymmetry: bool,
    pub friction: Friction,
    /// `Darcy Roughness` or `Manning Coefficient` (`nm`).
    pub roughness: f64,
    /// `|Gravity|` when the constants give it (`GotGrav`).
    pub gravity: Option<f64>,
    pub convect: bool,
    pub newton: bool,
    pub lateral_strain: bool,
    pub gradp: bool,
    pub transient: bool,
    /// `Min Speed` (default 1e-6).
    pub min_speed: f64,
    pub calc_friction_heating: bool,
    /// `Use Heating Source`: the material's `Latent Heat`.
    pub melt_heat: Option<f64>,
    /// `FirstRound`: the heating terms are zero on the very first iteration.
    pub first_round: bool,
}

/// Nodal inputs of one bulk element (`n` nodes; `velo[c][k]` is component
/// `c` of the previous iterate at function `k < n`).
pub struct BulkInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    /// Basis functions: the `n` nodal ones, then the bubbles.
    pub ntot: usize,
    pub rho: &'a [f64],
    pub mu: &'a [f64],
    pub gap: &'a [f64],
    /// `Bedrock Height` (`GotHeight`).
    pub height: Option<&'a [f64]>,
    pub pres: &'a [f64],
    pub velo: &'a [Vec<f64>],
    /// `LOAD(1:mdim+4, 1:n)`: body forces, normal and FSI velocity, flow
    /// admittance and external pressure.
    pub load: &'a [Vec<f64>],
}

/// `FilmFriction HeatFlux` contributions (`s * Basis(1:n) * q_f` summed over
/// the points) and the normalisation weight `s * Basis(1:n)`.
#[derive(Clone, Debug, Default)]
pub struct Heating {
    pub flux: Vec<f64>,
    pub weight: Vec<f64>,
}

/// `FrictionLawPraks(v, rho, nu, D, eps)`: the Darcy friction factor of
/// Praks and Brkic (2020), eq. (29), laminar below the turbulent regime.
pub fn friction_law_praks(v: f64, rho: f64, nu: f64, d: f64, eps: f64) -> f64 {
    let re = v * (d / 2.0) * rho / nu;
    let a = re * eps / 8.0897;
    let b = re.ln() - 0.779626;
    let x = a + b;
    if x <= AEPS {
        return 64.0 / re;
    }
    let c = x.ln();
    let lambda = re * (6.94871 * (b - c + c / (x - 0.5588 * c + 1.2079))).powf(-2.0);
    let lambda = lambda.max(1.0);
    64.0 * lambda / re
}

/// A single-precision Fortran literal quotient (`-2.0/3`) as the double the
/// expression is promoted to.
fn single(x: f32) -> f64 {
    f64::from(x)
}

/// `SUM(a(1:n) * b(1:n))` from zero.
fn dot(a: &[f64], b: impl Iterator<Item = f64>) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// `LocalBulkMatrix`. `out` gets `ntot * (mdim + 1)` rows and the mass
/// matrix when transient; the pressure rows of the bubble functions are
/// fixed to zero. Returns the heating contributions when asked.
pub fn local_bulk_matrix(
    basis: &dyn ElementBasis,
    rule: &[QuadPoint],
    inp: &BulkInput<'_>,
    o: &BulkOptions,
    out: &mut LocalSystem,
) -> Result<Heating, FilmFlowError> {
    let n = inp.xyz.len();
    let mdim = o.mdim;
    let dofs = mdim + 1;
    let ntot = inp.ntot;
    let m = ntot * dofs;
    out.n = m;
    out.stiff.clear();
    out.stiff.resize(m * m, 0.0);
    out.force.clear();
    out.force.resize(m, 0.0);
    out.damp = None;
    let mut mass = vec![0.0; m * m];
    let geomc = if o.csymmetry { 2.0 } else { 1.0 };
    let mut heating = Heating {
        flux: vec![0.0; n],
        weight: vec![0.0; n],
    };
    let mut be = BasisEval::default();
    for qp in rule {
        basis
            .eval(inp.etype, inp.xyz, qp.uvw, &mut be)
            .map_err(FilmFlowError::Basis)?;
        if be.basis.len() < ntot {
            return Err(FilmFlowError::Basis(format!(
                "{} basis functions evaluated, {ntot} needed",
                be.basis.len()
            )));
        }
        let bs = &be.basis;
        let db = &be.dbasis_dx;
        let s = qp.weight * be.det_j;

        let mu = dot(&bs[..n], inp.mu.iter().copied());
        let rho = dot(&bs[..n], inp.rho.iter().copied());
        let gap = dot(&bs[..n], inp.gap.iter().copied());
        let mut gap_grad = [0.0f64; 3];
        let mut pres_grad = [0.0f64; 3];
        for i in 0..mdim {
            gap_grad[i] = dot(inp.gap, db[..n].iter().map(|d| d[i]));
            pres_grad[i] = dot(inp.pres, db[..n].iter().map(|d| d[i]));
        }
        let mut h_grad = [0.0f64; 3];
        match inp.height {
            Some(h) => {
                for i in 0..mdim {
                    h_grad[i] = dot(h, db[..n].iter().map(|d| d[i]));
                }
            }
            None => {
                // The bed is the coordinate beyond the model dimension.
                let c = if mdim == 1 { 1 } else { 2 };
                let x: Vec<f64> = inp.xyz.iter().map(|p| p[c]).collect();
                for i in 0..mdim {
                    h_grad[i] = dot(&x, db[..n].iter().map(|d| d[i]));
                }
            }
        }
        // Velo = MATMUL(NodalVelo(1:mdim,1:nd), Basis(1:nd)), likewise the
        // gradient.
        let mut velo = [0.0f64; 3];
        let mut velo_grad = [[0.0f64; 3]; 3];
        for i in 0..mdim {
            velo[i] = dot(&inp.velo[i][..n], bs[..n].iter().copied());
            for j in 0..mdim {
                velo_grad[i][j] = dot(&inp.velo[i][..n], db[..n].iter().map(|d| d[j]));
            }
        }
        let mut speed = {
            let mut s2 = 0.0;
            for i in 0..mdim {
                s2 += velo[i] * velo[i];
            }
            s2.sqrt()
        };
        let mut load = [0.0f64; 7];
        for i in 0..mdim + 4 {
            load[i] = dot(&bs[..n], inp.load[i].iter().copied());
        }
        if o.convect && o.newton {
            for i in 0..mdim {
                let mut g = 0.0;
                for j in 0..mdim {
                    g += velo_grad[i][j] * velo[j];
                }
                load[i] += rho * g;
            }
        }
        load[mdim] *= geomc;
        load[mdim + 1] *= geomc;

        let gap2 = gap * gap;
        let gapi = gap;
        let grav = o.gravity.unwrap_or(0.0);
        let grad_zphi2 = || {
            let mut s2 = 0.0;
            for i in 0..mdim {
                let x = h_grad[i] + pres_grad[i] / (rho * grav);
                s2 += x * x;
            }
            s2.max(single(1.0e-9))
        };
        let mu_coeff = match o.friction {
            Friction::Darcy | Friction::Darcy2 => {
                speed = speed.max(o.min_speed);
                let d = 2.0 * gapi;
                let fd = friction_law_praks(speed, rho, mu, d, o.roughness);
                if o.friction == Friction::Darcy {
                    fd * rho * speed / (2.0 * d)
                } else {
                    rho * (fd * grav).sqrt()
                        * (2.0 * gap).powf(single(-1.0 / 2.0))
                        * grad_zphi2().powf(single(1.0 / 4.0))
                }
            }
            Friction::Manning => {
                o.roughness
                    * rho
                    * grav
                    * (gapi / 2.0).powf(single(-2.0 / 3.0))
                    * grad_zphi2().powf(single(1.0 / 4.0))
            }
            Friction::Manning2 => {
                rho * grav
                    * (o.roughness * o.roughness)
                    * (gapi / 2.0).powf(single(-4.0 / 3.0))
                    * speed
            }
            Friction::Laminar => {
                if o.csymmetry {
                    8.0 * mu / gap2
                } else {
                    12.0 * mu / gap2
                }
            }
        };
        let mut q_f = if o.calc_friction_heating {
            mu_coeff * gapi * (speed * speed)
        } else {
            0.0
        };
        if o.first_round {
            q_f = 0.0;
        }

        let pi = mdim; // the pressure component
        for p in 0..ntot {
            for q in 0..ntot {
                let at = |a: usize, b: usize| (p * dofs + a) * m + q * dofs + b;
                for i in 0..mdim {
                    if o.transient {
                        mass[at(i, i)] += s * rho * bs[q] * bs[p];
                    }
                    out.stiff[at(i, i)] += s * mu_coeff * bs[q] * bs[p];
                    for j in 0..mdim {
                        if o.lateral_strain {
                            out.stiff[at(i, i)] += s * mu * db[q][j] * db[p][j];
                            out.stiff[at(i, j)] += s * mu * db[q][i] * db[p][j];
                        }
                        if o.convect {
                            out.stiff[at(i, i)] += s * rho * velo[j] * db[q][j] * bs[p];
                            if o.newton {
                                out.stiff[at(i, j)] += s * rho * velo_grad[i][j] * bs[q] * bs[p];
                            }
                        }
                    }
                    if o.gradp {
                        out.stiff[at(i, pi)] += s * db[q][i] * bs[p];
                        out.stiff[at(pi, i)] -= s * gap * rho * bs[q] * db[p][i];
                    } else {
                        out.stiff[at(i, pi)] -= s * bs[q] * db[p][i];
                        let x = &mut out.stiff[at(pi, i)];
                        *x = *x
                            + s * gap * rho * db[q][i] * bs[p]
                            + geomc * s * rho * bs[q] * gap_grad[i] * bs[p];
                    }
                }
                out.stiff[at(pi, pi)] += s * rho * load[mdim + 2] * bs[q] * bs[p];
            }
            let f = &mut out.force[p * dofs..(p + 1) * dofs];
            for k in 0..=mdim {
                f[k] -= s * rho * bs[p] * load[k];
            }
            f[mdim] -= s * rho * bs[p] * load[mdim + 1];
            if o.gravity.is_some() {
                for k in 0..mdim {
                    f[k] -= s * rho * grav * bs[p] * h_grad[k];
                }
            }
            f[mdim] += s * rho * bs[p] * load[mdim + 2] * load[mdim + 3];
            if let Some(melt) = o.melt_heat {
                for k in 0..=mdim {
                    f[k] += s * bs[p] * q_f.max(0.0) / (melt * rho);
                }
            }
        }
        if o.calc_friction_heating {
            for k in 0..n {
                heating.flux[k] += s * bs[k] * q_f;
                heating.weight[k] += s * bs[k];
            }
        }
    }
    // p2/p1 elements: the pressure dofs of the bubbles are fixed.
    for p in n..ntot {
        let i = p * dofs + mdim;
        out.force[i] = 0.0;
        for j in 0..m {
            mass[j * m + i] = 0.0;
            mass[i * m + j] = 0.0;
            out.stiff[i * m + j] = 0.0;
            out.stiff[j * m + i] = 0.0;
        }
        out.stiff[i * m + i] = 1.0;
    }
    out.mass = o.transient.then_some(mass);
    Ok(heating)
}

/// `LocalBoundaryMatrix`: the force `s * Basis(p) * Load` of the boundary
/// loads `FilmPressure i` (`i <= mdim`) and `Mass Flux`.
pub fn local_boundary_matrix(
    basis: &dyn ElementBasis,
    rule: &[QuadPoint],
    etype: ElementType,
    xyz: &[[f64; 3]],
    load: &[Vec<f64>],
    mdim: usize,
    out: &mut LocalSystem,
) -> Result<(), FilmFlowError> {
    let n = xyz.len();
    let dofs = mdim + 1;
    let m = n * dofs;
    out.n = m;
    out.stiff.clear();
    out.stiff.resize(m * m, 0.0);
    out.force.clear();
    out.force.resize(m, 0.0);
    out.mass = None;
    out.damp = None;
    let mut be = BasisEval::default();
    for qp in rule {
        basis
            .eval(etype, xyz, qp.uvw, &mut be)
            .map_err(FilmFlowError::Basis)?;
        let s = qp.weight * be.det_j;
        let l: Vec<f64> = (0..dofs)
            .map(|i| dot(&be.basis[..n], load[i].iter().copied()))
            .collect();
        for p in 0..n {
            for k in 0..dofs {
                out.force[p * dofs + k] += s * be.basis[p] * l[k];
            }
        }
    }
    Ok(())
}

/// `LUDecomp` (LinearAlgebra.F90:56), column-major; a copy of elmer-fe's
/// private port.
fn lu_decomp(a: &mut [f64], n: usize, pivot: &mut [usize]) -> Option<()> {
    let ix = |i: usize, j: usize| i + j * n;
    for i in 0..n {
        let mut j = i;
        for k in i + 1..n {
            if a[ix(i, k)].abs() > a[ix(i, j)].abs() {
                j = k;
            }
        }
        if a[ix(i, j)].abs() == 0.0 {
            return None;
        }
        pivot[i] = j;
        if j != i {
            for k in 0..=i {
                a.swap(ix(k, j), ix(k, i));
            }
        }
        for k in i + 1..n {
            a[ix(i, k)] /= a[ix(i, i)];
        }
        for k in i + 1..n {
            if j != i {
                a.swap(ix(k, i), ix(k, j));
            }
            for l in i + 1..n {
                a[ix(k, l)] -= a[ix(k, i)] * a[ix(i, l)];
            }
        }
    }
    pivot[n - 1] = n - 1;
    for i in 0..n {
        if a[ix(i, i)].abs() == 0.0 {
            return None;
        }
        a[ix(i, i)] = 1.0 / a[ix(i, i)];
    }
    Some(())
}

/// `InvertMatrix` (LinearAlgebra.F90:191), column-major.
fn invert_matrix(a: &mut [f64], n: usize) -> Option<()> {
    let ix = |i: usize, j: usize| i + j * n;
    let mut pivot = vec![0usize; n];
    lu_decomp(a, n, &mut pivot)?;
    for i in (0..n.saturating_sub(1)).rev() {
        for j in (i + 1..n).rev() {
            let mut s = -a[ix(i, j)];
            for k in i + 1..j {
                s -= a[ix(i, k)] * a[ix(k, j)];
            }
            a[ix(i, j)] = s;
        }
    }
    for i in (0..n.saturating_sub(1)).rev() {
        for j in (i + 1..n).rev() {
            let mut s = 0.0;
            for k in i + 1..=j {
                s -= a[ix(j, k)] * a[ix(k, i)];
            }
            a[ix(j, i)] = a[ix(i, i)] * s;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in i.max(j)..n {
                if k != i {
                    s += a[ix(i, k)] * a[ix(k, j)];
                } else {
                    s += a[ix(k, j)];
                }
            }
            a[ix(i, j)] = s;
        }
    }
    for i in (0..n).rev() {
        if pivot[i] != i {
            for j in 0..n {
                a.swap(ix(i, j), ix(pivot[i], j));
            }
        }
    }
    Some(())
}

/// `NSCondensate(N = n, Nb = nb, dim = mdim, K, F)`: the velocity
/// components of the `nb` bubble functions are eliminated against the
/// `n * (mdim + 1)` retained rows; `local` keeps those (its mass, when
/// present, is cut to the same leading block, as `Default1stOrderTime` reads
/// it).
pub fn ns_condensate(
    local: &mut LocalSystem,
    n: usize,
    nb: usize,
    mdim: usize,
) -> Result<(), FilmFlowError> {
    if nb == 0 {
        return Ok(());
    }
    let dofs = mdim + 1;
    let m = local.n;
    let nc = n * dofs;
    let b: Vec<usize> = (0..nb)
        .flat_map(|p| (0..mdim).map(move |i| nc + p * dofs + i))
        .collect();
    let nbb = b.len();
    let k = |i: usize, j: usize| local.stiff[i * m + j];
    let mut inv = vec![0.0; nbb * nbb];
    for p in 0..nbb {
        for q in 0..nbb {
            inv[p + q * nbb] = k(b[p], b[q]);
        }
    }
    invert_matrix(&mut inv, nbb)
        .ok_or_else(|| FilmFlowError::Invalid("NSCondensate: a singular bubble block".into()))?;
    let fb: Vec<f64> = b.iter().map(|&i| local.force[i]).collect();
    // MATMUL(Kbb, Fb) and MATMUL(Kbb, Kbl).
    let t: Vec<f64> = (0..nbb)
        .map(|p| {
            let mut s = 0.0;
            for q in 0..nbb {
                s += inv[p + q * nbb] * fb[q];
            }
            s
        })
        .collect();
    let mut w = vec![0.0; nbb * nc];
    for p in 0..nbb {
        for c in 0..nc {
            let mut s = 0.0;
            for q in 0..nbb {
                s += inv[p + q * nbb] * k(b[q], c);
            }
            w[p * nc + c] = s;
        }
    }
    let mut stiff = vec![0.0; nc * nc];
    let mut force = vec![0.0; nc];
    for r in 0..nc {
        let mut s = 0.0;
        for p in 0..nbb {
            s += k(r, b[p]) * t[p];
        }
        force[r] = local.force[r] - s;
        for c in 0..nc {
            let mut s = 0.0;
            for p in 0..nbb {
                s += k(r, b[p]) * w[p * nc + c];
            }
            stiff[r * nc + c] = k(r, c) - s;
        }
    }
    let mass = local.mass.take().map(|mm| {
        let mut out = vec![0.0; nc * nc];
        for r in 0..nc {
            out[r * nc..(r + 1) * nc].copy_from_slice(&mm[r * m..r * m + nc]);
        }
        out
    });
    local.n = nc;
    local.stiff = stiff;
    local.force = force;
    local.mass = mass;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn praks_is_laminar_at_low_reynolds_number() {
        // Re = v (D/2) rho / nu = 1 * 0.5 * 1 / 1 = 0.5: lambda thresholds
        // to 1, f = 64 / Re.
        let f = friction_law_praks(1.0, 1.0, 1.0, 1.0, 0.01);
        assert!((f - 128.0).abs() < 1e-12, "{f}");
    }

    #[test]
    fn single_precision_exponents_are_promoted() {
        assert_eq!(single(-2.0 / 3.0), f64::from(-2.0f32 / 3.0f32));
        assert_ne!(single(-2.0 / 3.0), -2.0 / 3.0);
        assert_eq!(single(-1.0 / 2.0), -0.5);
    }
}
