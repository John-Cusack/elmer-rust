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
//! Bubbles, P2/P1, VMS, compressible models, porous drag, rotating frames,
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

    let nd = n * c;
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

    let points = quad.gauss_points(inp.etype, inp.np);
    if points.is_empty() {
        return Err(FlowError::NoQuadrature(inp.etype.0));
    }
    let (hk, mk) = (inp.stab.hk, inp.stab.mk);
    let mut su = vec![[[0.0; 4]; 4]; n];
    let mut sw = vec![[[0.0; 4]; 4]; n];
    // second[p][i][j] = sum_k (d phi_p / d x_i at node k) * d phi_k / d x_j.
    let mut second = vec![[[0.0; 3]; 3]; n];

    for qp in &points {
        basis
            .eval(inp.etype, inp.xyz, qp.uvw, &mut be)
            .map_err(FlowError::Basis)?;
        let s = be.det_j * qp.weight;
        let phi = &be.basis[..n];
        let dphi = &be.dbasis_dx[..n];

        let interp = |v: &dyn Fn(usize) -> f64| (0..n).map(|k| v(k) * phi[k]).sum::<f64>();
        let rho = interp(&|k| inp.density[k]);
        let mu = interp(&|k| inp.viscosity[k]);

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

        let mut force = [0.0; 4];
        for (i, f) in force.iter_mut().enumerate().take(c) {
            *f = interp(&|k| inp.load[k][i]);
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
        let mut dmudx = [0.0; 3];
        for (i, d) in dmudx.iter_mut().enumerate() {
            *d = (0..n).map(|k| inp.viscosity[k] * dphi[k][i]).sum();
        }
        let (tau, delta) = if inp.convect {
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
        for p in 0..n {
            for q in 0..n {
                // Viscous term.
                for i in 0..dim {
                    for j in 0..dim {
                        a[ij(at(p, i), at(q, i))] += s * mu * dphi[q][j] * dphi[p][j];
                        if disc.div {
                            a[ij(at(p, i), at(q, j))] += s * mu * dphi[q][j] * dphi[p][i];
                        } else if !disc.laplace {
                            a[ij(at(p, i), at(q, j))] += s * mu * dphi[q][i] * dphi[p][j];
                        }
                    }
                }
                // Convection, Picard and Newton parts.
                if inp.convect {
                    let conv: f64 = (0..dim).map(|j| dphi[q][j] * velo[j]).sum();
                    for i in 0..dim {
                        a[ij(at(p, i), at(q, i))] += s * rho * conv * phi[p];
                        if inp.newton {
                            for j in 0..dim {
                                a[ij(at(p, i), at(q, j))] += s * rho * grad[i][j] * phi[q] * phi[p];
                            }
                        }
                    }
                }
                // Pressure gradient and continuity.
                let cont = if inp.convect { rho } else { 1.0 };
                for i in 0..dim {
                    if disc.gradp {
                        a[ij(at(p, i), at(q, dim))] += s * dphi[q][i] * phi[p];
                        a[ij(at(p, dim), at(q, i))] -= s * rho * phi[q] * dphi[p][i];
                    } else {
                        a[ij(at(p, i), at(q, dim))] -= s * phi[q] * dphi[p][i];
                        a[ij(at(p, dim), at(q, i))] += s * cont * dphi[q][i] * phi[p];
                    }
                }
                if let Some(m) = mass.as_mut() {
                    for i in 0..dim {
                        m[ij(at(p, i), at(q, i))] += s * rho * phi[q] * phi[p];
                    }
                }

                // Stabilization.
                for i in 0..c {
                    for j in 0..c {
                        let swsu: f64 = (0..dim).map(|k| sw[p][i][k] * su[q][k][j]).sum();
                        a[ij(at(p, i), at(q, j))] += s * tau * swsu;
                    }
                }
                if let Some(m) = mass.as_mut() {
                    for j in 0..c {
                        for i in 0..dim {
                            m[ij(at(p, j), at(q, i))] += s * tau * rho * phi[q] * sw[p][j][i];
                        }
                    }
                }
                for i in 0..dim {
                    for j in 0..dim {
                        a[ij(at(p, j), at(q, i))] += s * delta * dphi[q][i] * dphi[p][j];
                    }
                }
            }
        }

        for p in 0..n {
            for (i, f) in force.iter().enumerate().take(c) {
                out.force[at(p, i)] += s * rho * f * phi[p];
            }
            for i in 0..dim {
                for j in 0..c {
                    out.force[at(p, j)] += s * tau * rho * force[i] * sw[p][j][i];
                }
            }
        }
    }
    out.mass = mass;
    Ok(())
}
