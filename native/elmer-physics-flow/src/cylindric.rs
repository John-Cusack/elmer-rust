//! `NavierStokesCylindricalCompose` (NavierStokesCylindrical.F90:58) for
//! `Coordinate System = Cylindric Symmetric`: the incompressible
//! Navier-Stokes equations in the (r, z) plane with the swirl velocity as a
//! third component, four dofs per basis function (`u_r`, `u_z`, `u_phi`, p).
//!
//! Ported: the p-bubble formulation (`Stabilization Method = pbubbles` on
//! p elements whose bubbles are condensed locally), Stokes and Navier-Stokes
//! with Picard or Newton linearization, the `Gradp Discretization` variant,
//! body forces; and [`ns_condensate`] (`NSCondensate`, MatrixAssembly.F90:777).
//! The stabilized, legacy-bubble and P2/P1 formulations, compressible and
//! pseudo-compressible models, porous drag, potential and Lorentz forces are
//! refused by the solver.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, QuadPoint};
use elmer_fe::{CoordinateKind, CoordinateSystem};

use crate::error::FlowError;

/// What `NavierStokesCylindricalCompose` reads for one element.
pub struct CylindricBulkInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    pub density: &'a [f64],
    pub viscosity: &'a [f64],
    /// `Ux, Uy, Uz` at the element nodes (the previous nonlinear iterate).
    pub velocity: &'a [[f64; 3]],
    /// `LoadVector(1:4)` at the nodes: the three momentum loads per unit
    /// mass and the continuity load.
    pub load: &'a [[f64; 4]],
    pub convect: bool,
    pub newton: bool,
    /// `Gradp Discretization`.
    pub gradp: bool,
    /// `Element % BDOFs`: the bubble functions after the `n` vertex ones
    /// (`PBubbles`, `NBasis = n + BDOFs`).
    pub bubbles: usize,
}

/// The local stiffness matrix and force vector of one bulk element,
/// `4 * (n + bubbles)` square, the bubbles' pressure rows set to the
/// identity. `basis` and `points` are the element's own (`ElementInfo`,
/// `GaussPointsAdapt(Element)`).
pub fn cylindric_bulk_local(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    inp: &CylindricBulkInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), FlowError> {
    let n = inp.xyz.len();
    let nbasis = n + inp.bubbles;
    const C: usize = 4;
    let nd = C * nbasis;
    for (what, got) in [
        ("density", inp.density.len()),
        ("viscosity", inp.viscosity.len()),
        ("velocity", inp.velocity.len()),
        ("load", inp.load.len()),
    ] {
        if got != n {
            return Err(FlowError::NodalLength { what, expected: n, got });
        }
    }
    if points.is_empty() {
        return Err(FlowError::NoQuadrature(inp.etype.0));
    }
    out.n = nd;
    out.stiff.clear();
    out.stiff.resize(nd * nd, 0.0);
    out.force.clear();
    out.force.resize(nd, 0.0);
    out.mass = None;
    out.damp = None;
    let csys = CoordinateSystem { kind: CoordinateKind::CylindricSymmetric, dim: 2 };
    let mut be = BasisEval::default();
    let sum_n = |f: &dyn Fn(usize) -> f64| {
        let mut s = 0.0;
        for k in 0..n {
            s += f(k);
        }
        s
    };

    for qp in points {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be).map_err(FlowError::Basis)?;
        if be.basis.len() < nbasis {
            return Err(FlowError::Basis(format!("{} basis functions, expected {nbasis}", be.basis.len())));
        }
        let phi = &be.basis;
        let dphi = &be.dbasis_dx;
        let x = sum_n(&|k| inp.xyz[k][0] * phi[k]);
        let y = sum_n(&|k| inp.xyz[k][1] * phi[k]);
        let z = sum_n(&|k| inp.xyz[k][2] * phi[k]);
        let ci = csys.info(x, y, z);
        let metric = [ci.metric[0][0], ci.metric[1][1], ci.metric[2][2]];
        // Symb(i,j,k) of the Fortran is symb(i-1, j-1, k-1).
        let symb = |i: usize, j: usize, k: usize| ci.symbols[i - 1][j - 1][k - 1];
        let s = ci.sqrt_metric * be.det_j * qp.weight;

        let density = sum_n(&|k| inp.density[k] * phi[k]);
        let mut velo = [0.0; 3];
        for (i, v) in velo.iter_mut().enumerate() {
            *v = sum_n(&|k| inp.velocity[k][i] * phi[k]);
        }
        let mut dvelodx = [[0.0; 3]; 3];
        if inp.newton {
            for i in 0..3 {
                for (a, row) in dvelodx.iter_mut().enumerate() {
                    row[i] = sum_n(&|k| inp.velocity[k][a] * dphi[k][i]);
                }
            }
        }
        let mut force = [0.0; 4];
        for (i, f) in force.iter_mut().enumerate() {
            *f = sum_n(&|k| inp.load[k][i] * phi[k]);
        }
        let viscosity = sum_n(&|k| inp.viscosity[k] * phi[k]);

        for p in 0..nbasis {
            for q in 0..nbasis {
                let (bp, bq) = (phi[p], phi[q]);
                let (dp, dq) = (dphi[p], dphi[q]);
                let base_p = bp;
                let mut a = [[0.0f64; 4]; 4];
                // Diffusive terms.
                for i in 0..3 {
                    for j in 0..3 {
                        a[i][i] += viscosity * metric[j] * dq[j] * dp[j];
                        a[i][j] += viscosity * metric[i] * dq[i] * dp[j];
                    }
                }
                a[0][0] += 2.0 * viscosity * metric[2] * bq * bp;
                a[2][2] -= 2.0 * viscosity * symb(3, 1, 3) * dq[0] * bp;
                if inp.convect {
                    for (i, row) in a.iter_mut().enumerate().take(3) {
                        for j in 0..3 {
                            row[i] += density * dq[j] * velo[j] * bp;
                        }
                    }
                    a[0][2] += density * symb(3, 3, 1) * bq * velo[2] * bp;
                    a[2][0] += density * symb(1, 3, 3) * bq * velo[2] * bp;
                    a[2][2] += density * symb(3, 1, 3) * bq * velo[0] * bp;
                    if inp.newton {
                        for (i, row) in a.iter_mut().enumerate().take(3) {
                            for j in 0..3 {
                                row[j] += density * dvelodx[i][j] * bq * bp;
                            }
                        }
                        a[0][2] += density * symb(3, 3, 1) * bq * velo[2] * bp;
                        a[2][0] += density * symb(3, 1, 3) * bq * velo[2] * bp;
                        a[2][2] += density * symb(1, 3, 3) * bq * velo[0] * bp;
                    }
                }
                // Pressure terms.
                if inp.gradp {
                    for (i, row) in a.iter_mut().enumerate().take(3) {
                        row[3] += metric[i] * dq[i] * bp;
                    }
                } else {
                    for (i, row) in a.iter_mut().enumerate().take(3) {
                        row[3] -= metric[i] * bq * dp[i];
                    }
                    a[0][3] += metric[2] * bq * symb(3, 3, 1) * bp;
                }
                // Continuity equation.
                for i in 0..3 {
                    a[3][i] += match (inp.convect, inp.gradp) {
                        (true, true) => -(density * bq * dp[i]),
                        (true, false) => density * dq[i] * base_p,
                        (false, true) => -(bq * dp[i]),
                        (false, false) => dq[i] * base_p,
                    };
                }
                if !inp.gradp {
                    if inp.convect {
                        a[3][0] += density * symb(1, 3, 3) * bq * base_p;
                    } else {
                        a[3][0] += symb(1, 3, 3) * bq * base_p;
                    }
                }
                for (i, row) in a.iter().enumerate() {
                    for (j, &v) in row.iter().enumerate() {
                        out.stiff[(C * p + i) * nd + C * q + j] += s * v;
                    }
                }
            }
        }

        if inp.convect && inp.newton {
            let mut uvelo = [0.0; 3];
            for (i, u) in uvelo.iter_mut().enumerate() {
                *u = sum_n(&|k| phi[k] * inp.velocity[k][i]);
            }
            for (i, f) in force.iter_mut().enumerate().take(3) {
                for (j, &u) in uvelo.iter().enumerate() {
                    *f += dvelodx[i][j] * u;
                }
            }
            force[0] += symb(3, 3, 1) * uvelo[2] * uvelo[2];
            force[2] += symb(3, 1, 3) * uvelo[0] * uvelo[2];
            force[2] += symb(1, 3, 3) * uvelo[2] * uvelo[0];
        }
        for p in 0..nbasis {
            for (i, &f) in force.iter().enumerate() {
                out.force[C * p + i] += s * (density * f * phi[p]);
            }
        }
    }

    // PBubbles: the bubbles carry no pressure.
    for i in n..nbasis {
        let k = C * i + 3;
        for j in 0..nd {
            out.stiff[k * nd + j] = 0.0;
            out.stiff[j * nd + k] = 0.0;
        }
        out.force[k] = 0.0;
        out.stiff[k * nd + k] = 1.0;
    }
    Ok(())
}

/// `NSCondensate( N, Nb, dim, K, F )`: eliminates the `dim` velocity dofs of
/// the `nb` bubbles (local dofs `(dim+1)*(n+p) + i`) from the
/// `(dim+1)*(n+nb)` system, leaving the leading `(dim+1)*n` block. The bubble
/// block is inverted by `InvertMatrix`, the products are `MATMUL`s.
pub fn ns_condensate(n: usize, nb: usize, dim: usize, local: &mut LocalSystem) -> Result<(), FlowError> {
    if nb == 0 {
        return Ok(());
    }
    let c = dim + 1;
    let m = local.n;
    let nl = c * n;
    let nbd = nb * dim;
    let bdofs: Vec<usize> = (0..nb).flat_map(|p| (0..dim).map(move |i| c * p + i + nl)).collect();
    let k = &local.stiff;
    // Kbb column-major, as InvertMatrix takes it.
    let mut kbb = vec![0.0; nbd * nbd];
    for (j, &bj) in bdofs.iter().enumerate() {
        for (i, &bi) in bdofs.iter().enumerate() {
            kbb[i + j * nbd] = k[bi * m + bj];
        }
    }
    elmer_fe::invert_matrix(&mut kbb, nbd)
        .ok_or_else(|| FlowError::Unsupported("NSCondensate: singular bubble block".into()))?;
    let kbb_at = |i: usize, j: usize| kbb[i + j * nbd];
    // MATMUL(Kbb, Fb) and MATMUL(Kbb, Kbl).
    let t: Vec<f64> = (0..nbd)
        .map(|i| {
            let mut s = 0.0;
            for (j, &bj) in bdofs.iter().enumerate() {
                s += kbb_at(i, j) * local.force[bj];
            }
            s
        })
        .collect();
    let mut tk = vec![0.0; nbd * nl];
    for col in 0..nl {
        for i in 0..nbd {
            let mut s = 0.0;
            for (j, &bj) in bdofs.iter().enumerate() {
                s += kbb_at(i, j) * k[bj * m + col];
            }
            tk[i * nl + col] = s;
        }
    }
    let mut stiff = vec![0.0; nl * nl];
    let mut force = vec![0.0; nl];
    for r in 0..nl {
        let mut s = 0.0;
        for (i, &bi) in bdofs.iter().enumerate() {
            s += k[r * m + bi] * t[i];
        }
        force[r] = local.force[r] - s;
        for col in 0..nl {
            let mut s = 0.0;
            for (i, &bi) in bdofs.iter().enumerate() {
                s += k[r * m + bi] * tk[i * nl + col];
            }
            stiff[r * nl + col] = k[r * m + col] - s;
        }
    }
    local.n = nl;
    local.stiff = stiff;
    local.force = force;
    local.mass = None;
    local.damp = None;
    Ok(())
}
