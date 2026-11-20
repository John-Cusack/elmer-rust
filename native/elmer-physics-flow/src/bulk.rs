//! `NavierStokesCompose` (NavierStokes.F90:59): the local matrices of the
//! incompressible Navier-Stokes equations in cartesian coordinates with
//! equal-order (P1/P1, Q1/Q1) velocity and pressure and the residual-based
//! stabilization of Franca et al. (`Stabilization Method = Stabilized`).
//!
//! Dofs are node-major: local row `p * (dim + 1) + i` is velocity component
//! `i` (or the pressure, `i = dim`) at element node `p`, as in Elmer.
//!
//! Ported: Stokes and Navier-Stokes with Picard or Newton linearization of
//! the convection term, the `Div`, `Gradp` and `Laplace Discretization`
//! variants, mesh velocity, and the mass matrix for a transient solve.
//! A rotating frame (`Angular Velocity`) adds its centrifugal force and
//! Coriolis term. Bubbles, P2/P1, VMS, compressible models, porous drag,
//! Lorentz and potential forces, anisotropic and non-newtonian viscosity are
//! not ported; the solver refuses cases that ask for them.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

use crate::error::FlowError;

/// `Element % hK` (already multiplied by `H scale`) and
/// `Element % StabilizationMK` of the element, as `elmer_fe::stab_param`
/// computes them; `mk = 0` when the mesh is not stabilized
/// (`Mesh % Stabilize`), which turns the stabilization off.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct StabParams {
    pub hk: f64,
    pub mk: f64,
}

/// The discretization switches of the solver section.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct Discretization {
    /// `Div Discretization`: the viscous term in divergence form.
    pub div: bool,
    /// `Gradp Discretization`: the pressure gradient is not integrated by parts.
    pub gradp: bool,
    /// `Laplace Discretization`: the viscous term as a vector Laplacian.
    pub laplace: bool,
}

/// What `NavierStokesCompose` reads for one element, at its nodes.
pub struct FlowBulkInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    /// Coordinate system dimension.
    pub dim: usize,
    pub density: &'a [f64],
    pub viscosity: &'a [f64],
    /// Velocity of the previous nonlinear iterate.
    pub velocity: &'a [[f64; 3]],
    /// `Mesh Velocity`, subtracted from the convecting velocity.
    pub mesh_velocity: Option<&'a [[f64; 3]]>,
    /// Body force per unit mass: components `0..dim` load the momentum
    /// equations, component `dim` the continuity equation (Elmer's
    /// `LoadVector`, multiplied by density inside the kernel).
    pub load: &'a [[f64; 4]],
    pub convect: bool,
    pub newton: bool,
    pub stab: StabParams,
    pub discretization: Discretization,
    /// Also return the mass matrix (transient solves).
    pub want_mass: bool,
    /// Gauss point count, `None` for the element default.
    pub np: Option<usize>,
    /// `Angular Velocity` of a rotating frame: the centrifugal force
    /// `-omega x (omega x x)` in the load and the Coriolis term `2 rho omega
    /// x u` in the matrix.
    pub omega: Option<[f64; 3]>,
    /// `Viscosity Model = power law` (EffectiveViscosity,
    /// MaterialModels.F90), `None` for a Newtonian viscosity.
    pub power_law: Option<PowerLaw<'a>>,
    /// Legacy bubbles (`Stabilization Method = bubbles`): `2n` basis
    /// functions, the nodal ones times the element bubble (ElementInfo's
    /// `EvalBubbleBasis`), on the `GaussPoints2` rule, without
    /// stabilization. `out` is `2n * (dim + 1)`; FlowSolve condenses the
    /// bubbles (`NSCondensate`).
    pub bubbles: bool,
    /// `Porous Media`: `Porous Resistivity` at the nodes, the drag `mu *
    /// Drag(i) * u_i` in momentum equation `i`.
    pub porous: Option<&'a [[f64; 3]]>,
}

/// The nodal parameters of the power-law viscosity `mu = Viscosity *
/// s^((c2-1)/2) / c4^(c2-1)`, `s` half the second invariant of the strain
/// rate, clamped below at `c3^2`.
#[derive(Clone, Copy)]
pub struct PowerLaw<'a> {
    /// `Viscosity Exponent` (c2).
    pub exponent: &'a [f64],
    /// `Critical Shear Rate` (c3).
    pub critical: Option<&'a [f64]>,
    /// `Nominal Shear Rate` (c4).
    pub nominal: Option<&'a [f64]>,
}

fn check_len(what: &'static str, expected: usize, got: usize) -> Result<(), FlowError> {
    if got == expected {
        Ok(())
    } else {
        Err(FlowError::NodalLength { what, expected, got })
    }
}

/// The local stiffness matrix, force vector and (if asked) mass matrix of
/// one bulk element. `out` is resized to `n * (dim + 1)`.
pub fn flow_bulk_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &FlowBulkInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), FlowError> {
    let n = inp.xyz.len();
    let dim = inp.dim;
    let c = dim + 1;
    check_len("element nodes", inp.etype.n_nodes(), n)?;
    check_len("density", n, inp.density.len())?;
    check_len("viscosity", n, inp.viscosity.len())?;
    check_len("velocity", n, inp.velocity.len())?;
    check_len("load", n, inp.load.len())?;
    if let Some(m) = inp.mesh_velocity {
        check_len("mesh velocity", n, m.len())?;
    }
    if !(2..=3).contains(&dim) {
        return Err(FlowError::Unsupported(format!("coordinate dimension {dim}")));
    }

    if inp.bubbles && inp.power_law.is_some() {
        return Err(FlowError::Unsupported("a power-law viscosity with legacy bubbles".into()));
    }
    if let Some(d) = inp.porous {
        check_len("porous resistivity", n, d.len())?;
    }
    // NBasis.
    let nbs = if inp.bubbles { 2 * n } else { n };
    let nd = nbs * c;
    out.n = nd;
    out.stiff.clear();
    out.stiff.resize(nd * nd, 0.0);
    out.force.clear();
    out.force.resize(nd, 0.0);
    out.damp = None;
    let mut mass = if inp.want_mass { Some(vec![0.0; nd * nd]) } else { None };

    let at = |p: usize, i: usize| p * c + i;
    let ij = |r: usize, s: usize| r * nd + s;

    // dnodal[k][a]: gradient of basis function a at element node k.
    let mut be = BasisEval::default();
    let def = elmer_fe::element_type(inp.etype.0).map_err(|e| FlowError::Unsupported(e.to_string()))?;
    let mut dnodal = vec![vec![[0.0; 3]; n]; n];
    for k in 0..n {
        let uvw = [def.node_u[k], def.node_v[k], def.node_w[k]];
        basis
            .eval(inp.etype, inp.xyz, uvw, &mut be)
            .map_err(FlowError::Basis)?;
        dnodal[k].copy_from_slice(&be.dbasis_dx[..n]);
    }

    let np = if inp.bubbles { Some(def.gauss_points2) } else { inp.np };
    let points = quad.gauss_points(inp.etype, np);
    if points.is_empty() {
        return Err(FlowError::NoQuadrature(inp.etype.0));
    }
    let (hk, mk) = (inp.stab.hk, inp.stab.mk);
    let mut su = vec![[[0.0; 4]; 4]; n];
    let mut sw = vec![[[0.0; 4]; 4]; n];
    // second[p][i][j] = sum_k (d phi_p / d x_i at node k) * d phi_k / d x_j.
    let mut second = vec![[[0.0; 3]; 3]; n];
    // The viscosity Jacobian (`JacM`) of a Newton iteration on a power-law
    // viscosity, and whether the last integration point asked for it
    // (`ViscNewtonLin` is the last point's value when it is applied).
    let mut jac = vec![0.0; if inp.power_law.is_some() && inp.newton { nd * nd } else { 0 }];
    let mut visc_newton = false;

    for qp in &points {
        basis
            .eval(inp.etype, inp.xyz, qp.uvw, &mut be)
            .map_err(FlowError::Basis)?;
        let s = be.det_j * qp.weight;
        let phi = &be.basis[..n];
        let dphi = &be.dbasis_dx[..n];
        // The basis of the assembly loops: the nodal one, then the bubbles.
        let (bphi, bdphi) = if inp.bubbles {
            bubble_basis(inp.etype, phi, dphi)?
        } else {
            (phi.to_vec(), dphi.to_vec())
        };

        let interp = |v: &dyn Fn(usize) -> f64| (0..n).map(|k| v(k) * phi[k]).sum::<f64>();
        let rho = interp(&|k| inp.density[k]);
        let mut mu = interp(&|k| inp.viscosity[k]);

        let mut velo = [0.0; 3];
        let mut uvelo = [0.0; 3];
        let mut grad = [[0.0; 3]; 3];
        for i in 0..dim {
            for k in 0..n {
                let mv = inp.mesh_velocity.map_or(0.0, |m| m[k][i]);
                velo[i] += (inp.velocity[k][i] - mv) * phi[k];
                uvelo[i] += inp.velocity[k][i] * phi[k];
                for j in 0..3 {
                    grad[i][j] += inp.velocity[k][i] * dphi[k][j];
                }
            }
        }
        if let Some(pl) = inp.power_law {
            let muder = power_law(&pl, phi, &grad, &mut mu);
            visc_newton = inp.newton && muder != 0.0;
            if visc_newton {
                // StrainS(p, j) = SUM(Strain(j,:) * dBasisdx(p,:)).
                let strain: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| (grad[i][j] + grad[j][i]) / 2.0));
                let ss: Vec<[f64; 3]> = (0..n)
                    .map(|p| std::array::from_fn(|j| (0..3).fold(0.0, |acc, k| acc + strain[j][k] * dphi[p][k])))
                    .collect();
                for i in 0..dim {
                    for j in 0..dim {
                        for q in 0..n {
                            for p in 0..n {
                                jac[ij(at(p, j), at(q, i))] += s * 2.0 * muder * 4.0 * ss[q][i] * ss[p][j];
                            }
                        }
                    }
                }
            }
        }

        let mut force = [0.0; 4];
        for (i, f) in force.iter_mut().enumerate().take(c) {
            *f = interp(&|k| inp.load[k][i]);
        }
        if let Some(om) = inp.omega {
            // omega x (omega x coord) = omega (omega . coord) - coord (omega . omega),
            // multiplied by density later on.
            let coord: [f64; 3] = std::array::from_fn(|a| (0..n).map(|k| phi[k] * inp.xyz[k][a]).sum());
            let oc: f64 = (0..3).map(|a| om[a] * coord[a]).sum();
            let oo: f64 = (0..3).map(|a| om[a] * om[a]).sum();
            for i in 0..dim {
                force[i] -= om[i] * oc;
            }
            for i in 0..dim {
                force[i] += coord[i] * oo;
            }
        }
        if inp.convect && inp.newton {
            for i in 0..dim {
                for j in 0..dim {
                    force[i] += grad[i][j] * uvelo[j];
                }
            }
        }

        // Stabilization: Tau, Delta and the residual (SU) and weight (SW)
        // operators of each basis function.
        // Drag(i) = SUM(NodalDrag(i,1:n) * Basis(1:n)).
        let drag: [f64; 3] = match inp.porous {
            Some(d) => std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + d[k][i] * phi[k])),
            None => [0.0; 3],
        };
        let mut dmudx = [0.0; 3];
        for (i, d) in dmudx.iter_mut().enumerate() {
            *d = (0..n).map(|k| inp.viscosity[k] * dphi[k][i]).sum();
        }
        let (tau, delta) = if inp.bubbles {
            (0.0, 0.0)
        } else if inp.convect {
            let vnorm = velo[..dim].iter().map(|v| v * v).sum::<f64>().sqrt().max(1.0e-12);
            let re = (rho * mk * hk * vnorm / (4.0 * mu)).min(1.0);
            (hk * re / (2.0 * rho * vnorm), rho * re * hk * vnorm)
        } else {
            (mk * hk * hk / (8.0 * mu), 0.0)
        };
        for p in 0..n {
            for i in 0..3 {
                for j in 0..3 {
                    second[p][i][j] = (0..n).map(|k| dnodal[k][p][i] * dphi[k][j]).sum();
                }
            }
        }
        for p in 0..n {
            let (u, w) = (&mut su[p], &mut sw[p]);
            *u = [[0.0; 4]; 4];
            *w = [[0.0; 4]; 4];
            let conv: f64 = (0..dim).map(|j| dphi[p][j] * velo[j]).sum();
            for i in 0..dim {
                u[i][dim] += dphi[p][i];
                if inp.porous.is_some() {
                    u[i][i] += mu * drag[i] * phi[p];
                }
                if inp.convect {
                    u[i][i] += rho * conv;
                }
                for j in 0..dim {
                    u[i][i] -= dmudx[j] * dphi[p][j];
                    u[i][j] -= dmudx[j] * dphi[p][i];
                    u[i][i] -= mu * second[p][j][j];
                    u[i][j] -= mu * second[p][i][j];
                }
                if inp.convect && inp.newton {
                    for j in 0..dim {
                        u[i][j] += rho * grad[i][j] * phi[p];
                    }
                }

                if inp.convect {
                    w[dim][i] += rho * dphi[p][i];
                    w[i][i] += rho * conv;
                } else {
                    w[dim][i] += dphi[p][i];
                }
                for j in 0..dim {
                    w[i][i] -= dmudx[j] * dphi[p][j];
                    w[j][i] -= dmudx[j] * dphi[p][i];
                    w[i][i] -= mu * second[p][j][j];
                    w[j][i] -= mu * second[p][i][j];
                }
            }
        }

        let a = &mut out.stiff;
        let disc = inp.discretization;
        for p in 0..nbs {
            for q in 0..nbs {
                // Viscous term.
                for i in 0..dim {
                    for j in 0..dim {
                        a[ij(at(p, i), at(q, i))] += s * mu * bdphi[q][j] * bdphi[p][j];
                        if disc.div {
                            a[ij(at(p, i), at(q, j))] += s * mu * bdphi[q][j] * bdphi[p][i];
                        } else if !disc.laplace {
                            a[ij(at(p, i), at(q, j))] += s * mu * bdphi[q][i] * bdphi[p][j];
                        }
                    }
                }
                // Convection, Picard and Newton parts.
                // NavierStokes.F90 adds the convection term one direction j at
                // a time, the Newton term of the same j after it.
                if inp.convect {
                    for i in 0..dim {
                        for j in 0..dim {
                            a[ij(at(p, i), at(q, i))] += s * rho * bdphi[q][j] * velo[j] * bphi[p];
                            if inp.newton {
                                a[ij(at(p, i), at(q, j))] += s * rho * grad[i][j] * bphi[q] * bphi[p];
                            }
                        }
                    }
                }
                // Pressure gradient and continuity.
                let cont = if inp.convect { rho } else { 1.0 };
                for i in 0..dim {
                    if disc.gradp {
                        a[ij(at(p, i), at(q, dim))] += s * bdphi[q][i] * bphi[p];
                        a[ij(at(p, dim), at(q, i))] -= s * rho * bphi[q] * bdphi[p][i];
                    } else {
                        a[ij(at(p, i), at(q, dim))] -= s * bphi[q] * bdphi[p][i];
                        a[ij(at(p, dim), at(q, i))] += s * cont * bdphi[q][i] * bphi[p];
                    }
                }
                if inp.porous.is_some() {
                    for i in 0..dim {
                        a[ij(at(p, i), at(q, i))] += s * mu * drag[i] * bphi[q] * bphi[p];
                    }
                }
                if let Some(m) = mass.as_mut() {
                    for i in 0..dim {
                        m[ij(at(p, i), at(q, i))] += s * rho * bphi[q] * bphi[p];
                    }
                }
                // Coriolis (NavierStokes.F90:971-1002). The Fortran's last 3D
                // update writes A(1,3) once more where A(3,1) is meant, so
                // A(1,3) nets out and A(3,1) stays unset; kept as it is.
                if let Some(om) = inp.omega {
                    let mc = 2.0 * s * rho * bphi[q] * bphi[p];
                    a[ij(at(p, 0), at(q, 1))] -= mc * om[2];
                    a[ij(at(p, 1), at(q, 0))] += mc * om[2];
                    if dim == 3 {
                        a[ij(at(p, 0), at(q, 2))] += mc * om[1];
                        a[ij(at(p, 1), at(q, 2))] -= mc * om[0];
                        a[ij(at(p, 2), at(q, 1))] += mc * om[0];
                        a[ij(at(p, 0), at(q, 2))] -= mc * om[1];
                    }
                }

                // Stabilization.
                if inp.bubbles {
                    continue;
                }
                for i in 0..c {
                    for j in 0..c {
                        let swsu: f64 = (0..dim).map(|k| sw[p][i][k] * su[q][k][j]).sum();
                        a[ij(at(p, i), at(q, j))] += s * tau * swsu;
                    }
                }
                if let Some(m) = mass.as_mut() {
                    for j in 0..c {
                        for i in 0..dim {
                            m[ij(at(p, j), at(q, i))] += s * tau * rho * bphi[q] * sw[p][j][i];
                        }
                    }
                }
                for i in 0..dim {
                    for j in 0..dim {
                        a[ij(at(p, j), at(q, i))] += s * delta * bdphi[q][i] * bdphi[p][j];
                    }
                }
            }
        }

        for p in 0..nbs {
            for (i, f) in force.iter().enumerate().take(c) {
                out.force[at(p, i)] += s * rho * f * bphi[p];
            }
            if inp.bubbles {
                continue;
            }
            for i in 0..dim {
                for j in 0..c {
                    out.force[at(p, j)] += s * tau * rho * force[i] * sw[p][j][i];
                }
            }
        }
    }
    if visc_newton {
        // StiffMatrix += JacM; FORCE += MATMUL(JacM, SOL), SOL the nodal
        // velocity with zero pressure.
        let mut sol = vec![0.0; nd];
        for p in 0..n {
            sol[at(p, 0)..at(p, 0) + dim].copy_from_slice(&inp.velocity[p][..dim]);
        }
        for (a, j) in out.stiff.iter_mut().zip(&jac) {
            *a += j;
        }
        for r in 0..nd {
            out.force[r] += (0..nd).fold(0.0, |acc, k| acc + jac[ij(r, k)] * sol[k]);
        }
    }
    out.mass = mass;
    Ok(())
}

/// `EvalBubbleBasis` (ElemInfo.F90): the `n` bubble functions `Basis(i) *
/// Bubble` after the nodal ones, `Bubble` the product of the linear
/// element's functions at opposite corners (all of them on a simplex).
#[allow(clippy::type_complexity)]
fn bubble_basis(etype: ElementType, phi: &[f64], dphi: &[[f64; 3]]) -> Result<(Vec<f64>, Vec<[f64; 3]>), FlowError> {
    let n = phi.len();
    // The factors of the bubble: LinBasis indices (0-based).
    let f: &[usize] = match etype.0 {
        202 => &[0, 1],
        303 => &[0, 1, 2],
        404 => &[0, 2],
        504 => &[0, 1, 2, 3],
        808 => &[0, 6],
        t => {
            return Err(FlowError::Unsupported(format!(
                "legacy bubbles on element type {t} (only linear elements; ElementInfo would evaluate the linear sub-element)"
            )))
        }
    };
    let bubble = f.iter().skip(1).fold(phi[f[0]], |acc, &k| acc * phi[k]);
    let mut b = phi.to_vec();
    let mut db = dphi.to_vec();
    for i in 0..n {
        b.push(phi[i] * bubble);
        let mut g = [0.0; 3];
        for (j, gj) in g.iter_mut().enumerate() {
            *gj = dphi[i][j] * bubble;
            // + Basis(i) * dLinBasisdx(f_m, j) * (product of the other factors).
            for m in 0..f.len() {
                let others = f.iter().enumerate().filter(|&(l, _)| l != m).map(|(_, &k)| phi[k]);
                let mut t = phi[i] * dphi[f[m]][j];
                for o in others {
                    t *= o;
                }
                *gj += t;
            }
        }
        db.push(g);
    }
    Ok((b, db))
}

/// `EffectiveViscosity` for `power law` in cartesian coordinates: replaces
/// `mu` (the interpolated `Viscosity`) by the effective viscosity and
/// returns its derivative `muder` with respect to `s`.
fn power_law(pl: &PowerLaw<'_>, phi: &[f64], grad: &[[f64; 3]; 3], mu: &mut f64) -> f64 {
    let at_ip = |v: &[f64]| (0..phi.len()).fold(0.0, |acc, k| acc + phi[k] * v[k]);
    // SecondInvariant (Cartesian): sum of (dV(i,j) + dV(j,i))^2.
    let mut sec = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let s = grad[i][j] + grad[j][i];
            sec += s * s;
        }
    }
    let visc = *mu;
    let c2 = at_ip(pl.exponent);
    let mut s = 0.5 * sec;
    let mut muder = if s != 0.0 { visc * (c2 - 1.0) / 2.0 * s.powf((c2 - 1.0) / 2.0 - 1.0) } else { 0.0 };
    if let Some(c3) = pl.critical.map(at_ip) {
        if s < c3 * c3 {
            s = c3 * c3;
            muder = 0.0;
        }
    }
    *mu = visc * s.powf((c2 - 1.0) / 2.0);
    if let Some(c4) = pl.nominal.map(at_ip) {
        *mu /= c4.powf(c2 - 1.0);
        muder /= c4.powf(c2 - 1.0);
    }
    muder
}
