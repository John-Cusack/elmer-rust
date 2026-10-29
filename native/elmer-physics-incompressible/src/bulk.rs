//! `LocalBulkMatrix` (IncompressibleNS.F90:120) for the equal-order nodal
//! element (`Element = "n:1"`, no bubble) held stable by `Pressure
//! Stabilization`, the Brezzi-Pitkaranta pressure block with the SUPG/PSPG
//! additions for convection and the time derivative.
//!
//! Dofs are node-major: local row `p * (dim + 1) + i` is velocity component
//! `i` (or the pressure, `i = dim`) at element node `p`, as in Elmer.
//!
//! Ported: Stokes and Navier-Stokes (steady and transient mass matrix),
//! Picard and Newton linearization of the convection, the `GradP` and
//! `Div-Curl` forms, element-constant density, nodal viscosity and body
//! forces, the three `Pressure Stabilization Element Size` measures. Not
//! ported (the solver refuses them): bubbles and p-elements, compressibility
//! models, non-newtonian and turbulent viscosity, axisymmetric coordinates.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Quadrature};

use crate::error::IncompressibleError;

/// `Pressure Stabilization Element Size`.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum HMode {
    /// `shortest edge` (default): `Element % hK`.
    #[default]
    Default,
    /// `longest edge`: `ElementDiameter(.., UseLongEdge = .TRUE.)`.
    LongestEdge,
    /// `volume`: `(sum of the weighted detJ) ** (1 / element dimension)`.
    Volume,
}

/// The pressure stabilization of the solver section.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PressureStab {
    /// `Pressure Stabilization Coefficient`.
    pub coeff: f64,
    pub hmode: HMode,
    /// `Element % hK` as `MeshStabParams` fills it.
    pub hk: f64,
}

/// What `LocalBulkMatrix` reads for one element, at its nodes.
pub struct IncompressibleBulkInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    /// Coordinate system dimension.
    pub dim: usize,
    /// `Density`, constant over the element.
    pub density: f64,
    /// `Viscosity` at the nodes.
    pub viscosity: &'a [f64],
    /// `Flow Bodyforce 1..3` at the nodes.
    pub load: &'a [[f64; 3]],
    /// Velocity and pressure of the previous iterate at the nodes
    /// (`GetLocalSolution`); `None` for `Linear Equation = True`.
    pub solution: Option<&'a [[f64; 4]]>,
    pub stokes: bool,
    pub newton: bool,
    /// Also return the mass matrix (a transient, non-Stokes run).
    pub want_mass: bool,
    /// `GradP Discretization`.
    pub gradp: bool,
    /// `Div-Curl Discretization`.
    pub div_curl: bool,
    /// `Specific Load`: the body force is per unit volume, not per mass.
    pub specific_load: bool,
    pub pstab: PressureStab,
    /// Gauss point count, `None` for the element default.
    pub np: Option<usize>,
}

fn check_len(what: &'static str, expected: usize, got: usize) -> Result<(), IncompressibleError> {
    if got == expected {
        Ok(())
    } else {
        Err(IncompressibleError::NodalLength { what, expected, got })
    }
}

/// Element dimension (`Element % Type % Dimension`).
fn element_dim(etype: ElementType) -> usize {
    match etype.family() {
        1 => 0,
        2 => 1,
        3 | 4 => 2,
        _ => 3,
    }
}

/// The per-family stabilization constant of LocalBulkMatrix.
fn family_constant(etype: ElementType) -> f64 {
    match etype.0 / 100 {
        3 => 0.0096,
        4 => 0.0263,
        5 => 0.0029,
        6 => 0.005,
        7 => 0.0090,
        _ => 0.0232,
    }
}

/// `ListGetElementRealVec` of a nodal field at the integration points: the
/// constant itself when the nodal values agree, else the basis interpolation.
fn at_points(vals: &[f64], bas: &[f64], n: usize, ngp: usize) -> Vec<f64> {
    if vals.iter().all(|&v| v == vals[0]) {
        return vec![vals[0]; ngp];
    }
    (0..ngp).map(|t| (0..n).map(|k| bas[t * n + k] * vals[k]).sum()).collect()
}

/// Local blocks `[i][j]` of `n x n` entries.
struct Blocks {
    n: usize,
    dofs: usize,
    v: Vec<f64>,
}

impl Blocks {
    fn new(n: usize, dofs: usize) -> Self {
        Self {
            n,
            dofs,
            v: vec![0.0; n * n * dofs * dofs],
        }
    }

    fn at(&mut self, i: usize, j: usize) -> &mut [f64] {
        let s = self.n * self.n;
        let o = (i * self.dofs + j) * s;
        &mut self.v[o..o + s]
    }

    /// Node-major scatter: `M(i::dofs, j::dofs) = block(i, j)`.
    fn scatter(&self) -> Vec<f64> {
        let (n, dofs) = (self.n, self.dofs);
        let nd = n * dofs;
        let mut out = vec![0.0; nd * nd];
        for i in 0..dofs {
            for j in 0..dofs {
                let o = (i * dofs + j) * n * n;
                for p in 0..n {
                    for q in 0..n {
                        out[(p * dofs + i) * nd + q * dofs + j] = self.v[o + p * n + q];
                    }
                }
            }
        }
        out
    }
}

/// `LinearForms_UdotV`: `G(p, q) += sum_t U(t, p) V(t, q) w(t) [alpha(t)]`.
fn udotv(g: &mut [f64], n: usize, ngp: usize, u: &[f64], v: &[f64], w: &[f64], alpha: Option<&[f64]>) {
    for t in 0..ngp {
        let a = alpha.map_or(1.0, |a| a[t]);
        for p in 0..n {
            for q in 0..n {
                g[p * n + q] += match alpha {
                    None => u[t * n + p] * v[t * n + q] * w[t],
                    Some(_) => u[t * n + p] * v[t * n + q] * w[t] * a,
                };
            }
        }
    }
}

/// `LinearForms_UdotF`: `F(p) += sum_t U(t, p) f(t) w(t)`.
fn udotf(out: &mut [f64], n: usize, ngp: usize, u: &[f64], w: &[f64], f: &[f64]) {
    for t in 0..ngp {
        for p in 0..n {
            out[p] += u[t * n + p] * f[t] * w[t];
        }
    }
}

/// The local stiffness matrix, force vector and (if asked) mass matrix of
/// one bulk element. `out` is resized to `n * (dim + 1)`.
pub fn incompressible_bulk_local(
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    inp: &IncompressibleBulkInput<'_>,
    out: &mut LocalSystem,
) -> Result<(), IncompressibleError> {
    let n = inp.xyz.len();
    let dim = inp.dim;
    let dofs = dim + 1;
    check_len("element nodes", inp.etype.n_nodes(), n)?;
    check_len("viscosity", n, inp.viscosity.len())?;
    check_len("load", n, inp.load.len())?;
    if let Some(s) = inp.solution {
        check_len("solution", n, s.len())?;
    }
    if !(2..=3).contains(&dim) {
        return Err(IncompressibleError::Unsupported(format!("coordinate dimension {dim}")));
    }
    let nd = n * dofs;
    let stokes = inp.stokes;
    let rho = inp.density;

    let points = quad.gauss_points(inp.etype, inp.np);
    let ngp = points.len();
    if ngp == 0 {
        return Err(IncompressibleError::NoQuadrature(inp.etype.0));
    }

    // BasisVec(t, p), dBasisdxVec(t, p, a) as bas / dbs[a], and the weights.
    let mut be = BasisEval::default();
    let mut bas = vec![0.0; ngp * n];
    let mut dbs = [vec![0.0; ngp * n], vec![0.0; ngp * n], vec![0.0; ngp * n]];
    let mut detj = vec![0.0; ngp];
    for (t, qp) in points.iter().enumerate() {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be).map_err(IncompressibleError::Basis)?;
        detj[t] = be.det_j * qp.weight;
        for p in 0..n {
            bas[t * n + p] = be.basis[p];
            for (a, d) in dbs.iter_mut().enumerate() {
                d[t * n + p] = be.dbasis_dx[p][a];
            }
        }
    }

    // Velocity and pressure of the previous iterate at the points.
    let mut velopres = vec![[0.0; 4]; ngp];
    if let Some(sol) = inp.solution {
        for (t, vp) in velopres.iter_mut().enumerate() {
            for (k, x) in vp.iter_mut().enumerate().take(dofs) {
                *x = (0..n).map(|p| bas[t * n + p] * sol[p][k]).sum();
            }
        }
    }

    let mu = at_points(inp.viscosity, &bas, n, ngp);
    let rhovec = vec![rho; ngp];

    // Body force, with the continuity row empty.
    let mut load_at = vec![[0.0; 4]; ngp];
    for i in 0..dim {
        let col: Vec<f64> = inp.load.iter().map(|l| l[i]).collect();
        let lv = at_points(&col, &bas, n, ngp);
        for t in 0..ngp {
            load_at[t][i] = if inp.specific_load { lv[t] } else { rhovec[t] * lv[t] };
        }
    }

    // GradVec(t, i, j) = d u_i / d x_j of the previous iterate.
    let mut grad = vec![[[0.0; 3]; 3]; ngp];
    if inp.newton {
        let sol = inp.solution;
        for (t, g) in grad.iter_mut().enumerate() {
            for i in 0..dim {
                for j in 0..dim {
                    g[i][j] = sol.map_or(0.0, |s| (0..n).map(|q| dbs[j][t * n + q] * s[q][i]).sum());
                }
            }
        }
        if !stokes {
            for t in 0..ngp {
                for i in 0..dim {
                    let s: f64 = (0..dim).map(|j| grad[t][i][j] * velopres[t][j]).sum();
                    load_at[t][i] += rhovec[t] * s;
                }
            }
        }
    }

    let mut so = Blocks::new(n, dofs);
    let mut mass = if inp.want_mass { Some(Blocks::new(n, dofs)) } else { None };

    let weight_a: Vec<f64> = (0..ngp).map(|t| mu[t] * detj[t]).collect();
    // Viscous terms.
    if inp.div_curl {
        let weight_b: Vec<f64> = weight_a.iter().map(|w| -w).collect();
        let weight_c: Vec<f64> = weight_a.iter().map(|w| 4.0 / 3.0 * w).collect();
        for i in 0..dim {
            for j in 0..dim {
                udotv(so.at(i, j), n, ngp, &dbs[i], &dbs[j], &weight_c, None);
                if i == j {
                    for k in 0..dim {
                        if k != i {
                            udotv(so.at(i, j), n, ngp, &dbs[k], &dbs[k], &weight_a, None);
                        }
                    }
                } else {
                    udotv(so.at(i, j), n, ngp, &dbs[j], &dbs[i], &weight_b, None);
                }
            }
        }
    } else {
        for i in 0..dim {
            for j in 0..dim {
                udotv(so.at(i, i), n, ngp, &dbs[j], &dbs[j], &weight_a, None);
                udotv(so.at(i, j), n, ngp, &dbs[j], &dbs[i], &weight_a, None);
            }
        }
    }

    // Pressure gradient and continuity.
    let neg_detj: Vec<f64> = detj.iter().map(|d| -d).collect();
    for i in 0..dim {
        let blk = if inp.gradp {
            udotv(so.at(i, dim), n, ngp, &bas, &dbs[i], &detj, None);
            so.at(i, dim).to_vec()
        } else {
            udotv(so.at(i, dim), n, ngp, &dbs[i], &bas, &neg_detj, None);
            so.at(i, dim).to_vec()
        };
        let tr = so.at(dim, i);
        for p in 0..n {
            for q in 0..n {
                tr[p * n + q] = blk[q * n + p];
            }
        }
    }

    // Pressure stabilization.
    let mut tau = vec![0.0; ngp];
    let mut conv = vec![0.0; ngp * n];
    let kfam = family_constant(inp.etype);
    let h_stab = match inp.pstab.hmode {
        HMode::Default => inp.pstab.hk,
        HMode::LongestEdge => elmer_fe::element_diameter(inp.etype, inp.xyz, true),
        HMode::Volume => detj.iter().sum::<f64>().powf(1.0 / element_dim(inp.etype) as f64),
    };
    if h_stab <= 0.0 {
        return Err(IncompressibleError::ZeroElementSize);
    }
    if stokes {
        for t in 0..ngp {
            tau[t] = inp.pstab.coeff * kfam * h_stab.powi(2) / mu[t];
        }
    } else {
        for t in 0..ngp {
            let vnorm = (0..dim).map(|k| velopres[t][k].powi(2)).sum::<f64>().sqrt();
            tau[t] = if rhovec[t] * vnorm * h_stab > 0.0 {
                let pe = (2.0 * kfam * rhovec[t] * vnorm * h_stab / mu[t]).min(1.0);
                h_stab * pe / (2.0 * rhovec[t] * vnorm)
            } else {
                kfam * h_stab.powi(2) / mu[t]
            };
            tau[t] *= inp.pstab.coeff;
        }
        for k in 0..dim {
            for t in 0..ngp {
                for q in 0..n {
                    conv[t * n + q] += rhovec[t] * velopres[t][k] * dbs[k][t * n + q];
                }
            }
        }
    }
    let wa: Vec<f64> = (0..ngp).map(|t| -tau[t] * detj[t]).collect();
    let wb: Vec<f64> = (0..ngp).map(|t| tau[t] * detj[t]).collect();
    for i in 0..dim {
        udotv(so.at(dim, dim), n, ngp, &dbs[i], &dbs[i], &wa, None);
    }
    if !stokes {
        for j in 0..dim {
            udotv(so.at(dim, j), n, ngp, &dbs[j], &conv, &wa, None);
        }
        for i in 0..dim {
            udotv(so.at(i, i), n, ngp, &conv, &conv, &wb, None);
        }
        for i in 0..dim {
            udotv(so.at(i, dim), n, ngp, &conv, &dbs[i], &wb, None);
        }
        if let Some(m) = mass.as_mut() {
            let mut part = vec![0.0; n * n];
            udotv(&mut part, n, ngp, &conv, &bas, &wb, Some(&rhovec));
            for i in 0..dim {
                for (x, p) in m.at(i, i).iter_mut().zip(&part) {
                    *x += p;
                }
            }
            for i in 0..dim {
                let mut part = vec![0.0; n * n];
                udotv(&mut part, n, ngp, &dbs[i], &bas, &wa, Some(&rhovec));
                for (x, p) in m.at(dim, i).iter_mut().zip(&part) {
                    *x += p;
                }
            }
        }
        if inp.newton {
            for i in 0..dim {
                for j in 0..dim {
                    let al: Vec<f64> = (0..ngp).map(|t| rhovec[t] * grad[t][i][j]).collect();
                    udotv(so.at(i, j), n, ngp, &conv, &bas, &wb, Some(&al));
                    udotv(so.at(dim, j), n, ngp, &dbs[i], &bas, &wa, Some(&al));
                }
            }
        }
    }

    // Galerkin mass and convection.
    if !stokes {
        if let Some(m) = mass.as_mut() {
            let mut vm = vec![0.0; n * n];
            udotv(&mut vm, n, ngp, &bas, &bas, &detj, Some(&rhovec));
            for i in 0..dim {
                for (x, p) in m.at(i, i).iter_mut().zip(&vm) {
                    *x += p;
                }
            }
        }
        let mut vm = vec![0.0; n * n];
        for k in 0..dim {
            let wk: Vec<f64> = (0..ngp).map(|t| rhovec[t] * velopres[t][k]).collect();
            // weight_a is the first weight of the Fortran call: detJ, with
            // rho*u_k as the extra factor.
            udotv(&mut vm, n, ngp, &bas, &dbs[k], &detj, Some(&wk));
        }
        for i in 0..dim {
            for (x, p) in so.at(i, i).iter_mut().zip(&vm) {
                *x += p;
            }
        }
        if inp.newton {
            for i in 0..dim {
                for j in 0..dim {
                    let al: Vec<f64> = (0..ngp).map(|t| rhovec[t] * grad[t][i][j]).collect();
                    udotv(so.at(i, j), n, ngp, &bas, &bas, &detj, Some(&al));
                }
            }
        }
    }

    // Loads.
    let mut force = vec![0.0; nd];
    for i in 0..dofs {
        let f: Vec<f64> = (0..ngp).map(|t| load_at[t][i]).collect();
        let mut part = vec![0.0; n];
        udotf(&mut part, n, ngp, &bas, &detj, &f);
        for p in 0..n {
            force[p * dofs + i] = part[p];
        }
    }
    let mut part = vec![0.0; n];
    for i in 0..dim {
        let f: Vec<f64> = (0..ngp).map(|t| load_at[t][i]).collect();
        udotf(&mut part, n, ngp, &dbs[i], &wa, &f);
    }
    for p in 0..n {
        force[p * dofs + dim] += part[p];
    }
    if !stokes {
        for i in 0..dim {
            let f: Vec<f64> = (0..ngp).map(|t| load_at[t][i]).collect();
            let mut part = vec![0.0; n];
            udotf(&mut part, n, ngp, &conv, &wb, &f);
            for p in 0..n {
                force[p * dofs + i] += part[p];
            }
        }
    }

    out.n = nd;
    out.stiff = so.scatter();
    out.force = force;
    out.damp = None;
    out.mass = mass.map(|m| m.scatter());
    Ok(())
}
