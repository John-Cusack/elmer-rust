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
//! forces, the three `Pressure Stabilization Element Size` measures, and p
//! elements: Taylor-Hood (the edge and face pressure dofs pinned,
//! [`pin_pressure`]) and the MINI element (velocity bubbles condensed,
//! [`lcondensate`], [`lcondensate_bdf1`]). Not ported (the solver refuses
//! them): compressibility models, non-newtonian and turbulent viscosity,
//! axisymmetric coordinates.

// A structural transliteration of Fortran loops over `1:dim` that index several
// arrays by the same component number.
#![allow(clippy::needless_range_loop)]

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, QuadPoint, Quadrature};

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
    /// Velocity and pressure of the previous iterate at the nodes, then at
    /// the other basis functions of a p element (`GetLocalSolution`); `None`
    /// for `Linear Equation = True`.
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
/// `bas` holds `nf` functions per point, the `n = vals.len()` nodal ones first.
fn at_points(vals: &[f64], bas: &[f64], nf: usize, ngp: usize) -> Vec<f64> {
    if vals.iter().all(|&v| v == vals[0]) {
        return vec![vals[0]; ngp];
    }
    (0..ngp).map(|t| vals.iter().enumerate().map(|(k, v)| bas[t * nf + k] * v).sum()).collect()
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
    let points = quad.gauss_points(inp.etype, inp.np);
    incompressible_bulk_local_with(basis, &points, inp.xyz.len(), inp, out)
}

/// [`incompressible_bulk_local`] for `nf >= n` basis functions (the nodal
/// ones first, then the edge, face and bubble functions of a p element:
/// `LocalBulkMatrix(Element, n, nd, ntot)` with `ntot = nf`) on the given
/// rule. `inp.solution` holds `nf` entries; material data and loads are
/// interpolated with the nodal functions. `out` is resized to
/// `nf * (dim + 1)`; see [`pin_pressure`] for the pressure of a
/// Taylor-Hood element.
pub fn incompressible_bulk_local_with(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    nf: usize,
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
        check_len("solution", nf, s.len())?;
    }
    if !(2..=3).contains(&dim) {
        return Err(IncompressibleError::Unsupported(format!("coordinate dimension {dim}")));
    }
    if nf < n {
        return Err(IncompressibleError::Unsupported(format!("{nf} basis functions for {n} nodes")));
    }
    let nd = nf * dofs;
    let stokes = inp.stokes;
    let rho = inp.density;

    let ngp = points.len();
    if ngp == 0 {
        return Err(IncompressibleError::NoQuadrature(inp.etype.0));
    }

    // BasisVec(t, p), dBasisdxVec(t, p, a) as bas / dbs[a], and the weights.
    let mut be = BasisEval::default();
    let mut bas = vec![0.0; ngp * nf];
    let mut dbs = [vec![0.0; ngp * nf], vec![0.0; ngp * nf], vec![0.0; ngp * nf]];
    let mut detj = vec![0.0; ngp];
    for (t, qp) in points.iter().enumerate() {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be).map_err(IncompressibleError::Basis)?;
        if be.basis.len() < nf {
            return Err(IncompressibleError::Unsupported(format!("{} basis functions for {nf} dofs", be.basis.len())));
        }
        detj[t] = be.det_j * qp.weight;
        for p in 0..nf {
            bas[t * nf + p] = be.basis[p];
            for (a, d) in dbs.iter_mut().enumerate() {
                d[t * nf + p] = be.dbasis_dx[p][a];
            }
        }
    }

    // Velocity and pressure of the previous iterate at the points: the DGEMM
    // of LocalBulkMatrix runs over the `n` nodal functions only (the
    // velocity gradient below takes all `nf`).
    let mut velopres = vec![[0.0; 4]; ngp];
    if let Some(sol) = inp.solution {
        for (t, vp) in velopres.iter_mut().enumerate() {
            for (k, x) in vp.iter_mut().enumerate().take(dofs) {
                *x = (0..n).map(|p| bas[t * nf + p] * sol[p][k]).sum();
            }
        }
    }

    let mu = at_points(inp.viscosity, &bas, nf, ngp);
    let rhovec = vec![rho; ngp];

    // Body force, with the continuity row empty.
    let mut load_at = vec![[0.0; 4]; ngp];
    for i in 0..dim {
        let col: Vec<f64> = inp.load.iter().map(|l| l[i]).collect();
        let lv = at_points(&col, &bas, nf, ngp);
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
                    g[i][j] = sol.map_or(0.0, |s| (0..nf).map(|q| dbs[j][t * nf + q] * s[q][i]).sum());
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

    let mut so = Blocks::new(nf, dofs);
    let mut mass = if inp.want_mass { Some(Blocks::new(nf, dofs)) } else { None };

    let weight_a: Vec<f64> = (0..ngp).map(|t| mu[t] * detj[t]).collect();
    // Viscous terms.
    if inp.div_curl {
        let weight_b: Vec<f64> = weight_a.iter().map(|w| -w).collect();
        let weight_c: Vec<f64> = weight_a.iter().map(|w| 4.0 / 3.0 * w).collect();
        for i in 0..dim {
            for j in 0..dim {
                udotv(so.at(i, j), nf, ngp, &dbs[i], &dbs[j], &weight_c, None);
                if i == j {
                    for k in 0..dim {
                        if k != i {
                            udotv(so.at(i, j), nf, ngp, &dbs[k], &dbs[k], &weight_a, None);
                        }
                    }
                } else {
                    udotv(so.at(i, j), nf, ngp, &dbs[j], &dbs[i], &weight_b, None);
                }
            }
        }
    } else {
        for i in 0..dim {
            for j in 0..dim {
                udotv(so.at(i, i), nf, ngp, &dbs[j], &dbs[j], &weight_a, None);
                udotv(so.at(i, j), nf, ngp, &dbs[j], &dbs[i], &weight_a, None);
            }
        }
    }

    // Pressure gradient and continuity.
    let neg_detj: Vec<f64> = detj.iter().map(|d| -d).collect();
    for i in 0..dim {
        let blk = if inp.gradp {
            udotv(so.at(i, dim), nf, ngp, &bas, &dbs[i], &detj, None);
            so.at(i, dim).to_vec()
        } else {
            udotv(so.at(i, dim), nf, ngp, &dbs[i], &bas, &neg_detj, None);
            so.at(i, dim).to_vec()
        };
        let tr = so.at(dim, i);
        for p in 0..nf {
            for q in 0..nf {
                tr[p * nf + q] = blk[q * nf + p];
            }
        }
    }

    // Pressure stabilization.
    let mut tau = vec![0.0; ngp];
    let mut conv = vec![0.0; ngp * nf];
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
                for q in 0..nf {
                    conv[t * nf + q] += rhovec[t] * velopres[t][k] * dbs[k][t * nf + q];
                }
            }
        }
    }
    let wa: Vec<f64> = (0..ngp).map(|t| -tau[t] * detj[t]).collect();
    let wb: Vec<f64> = (0..ngp).map(|t| tau[t] * detj[t]).collect();
    for i in 0..dim {
        udotv(so.at(dim, dim), nf, ngp, &dbs[i], &dbs[i], &wa, None);
    }
    if !stokes {
        for j in 0..dim {
            udotv(so.at(dim, j), nf, ngp, &dbs[j], &conv, &wa, None);
        }
        for i in 0..dim {
            udotv(so.at(i, i), nf, ngp, &conv, &conv, &wb, None);
        }
        for i in 0..dim {
            udotv(so.at(i, dim), nf, ngp, &conv, &dbs[i], &wb, None);
        }
        if let Some(m) = mass.as_mut() {
            let mut part = vec![0.0; nf * nf];
            udotv(&mut part, nf, ngp, &conv, &bas, &wb, Some(&rhovec));
            for i in 0..dim {
                for (x, p) in m.at(i, i).iter_mut().zip(&part) {
                    *x += p;
                }
            }
            for i in 0..dim {
                let mut part = vec![0.0; nf * nf];
                udotv(&mut part, nf, ngp, &dbs[i], &bas, &wa, Some(&rhovec));
                for (x, p) in m.at(dim, i).iter_mut().zip(&part) {
                    *x += p;
                }
            }
        }
        if inp.newton {
            for i in 0..dim {
                for j in 0..dim {
                    let al: Vec<f64> = (0..ngp).map(|t| rhovec[t] * grad[t][i][j]).collect();
                    udotv(so.at(i, j), nf, ngp, &conv, &bas, &wb, Some(&al));
                    udotv(so.at(dim, j), nf, ngp, &dbs[i], &bas, &wa, Some(&al));
                }
            }
        }
    }

    // Galerkin mass and convection.
    if !stokes {
        if let Some(m) = mass.as_mut() {
            let mut vm = vec![0.0; nf * nf];
            udotv(&mut vm, nf, ngp, &bas, &bas, &detj, Some(&rhovec));
            for i in 0..dim {
                for (x, p) in m.at(i, i).iter_mut().zip(&vm) {
                    *x += p;
                }
            }
        }
        let mut vm = vec![0.0; nf * nf];
        for k in 0..dim {
            let wk: Vec<f64> = (0..ngp).map(|t| rhovec[t] * velopres[t][k]).collect();
            // weight_a is the first weight of the Fortran call: detJ, with
            // rho*u_k as the extra factor.
            udotv(&mut vm, nf, ngp, &bas, &dbs[k], &detj, Some(&wk));
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
                    udotv(so.at(i, j), nf, ngp, &bas, &bas, &detj, Some(&al));
                }
            }
        }
    }

    // Loads.
    let mut force = vec![0.0; nd];
    for i in 0..dofs {
        let f: Vec<f64> = (0..ngp).map(|t| load_at[t][i]).collect();
        let mut part = vec![0.0; nf];
        udotf(&mut part, nf, ngp, &bas, &detj, &f);
        for p in 0..nf {
            force[p * dofs + i] = part[p];
        }
    }
    let mut part = vec![0.0; nf];
    for i in 0..dim {
        let f: Vec<f64> = (0..ngp).map(|t| load_at[t][i]).collect();
        udotf(&mut part, nf, ngp, &dbs[i], &wa, &f);
    }
    for p in 0..nf {
        force[p * dofs + dim] += part[p];
    }
    if !stokes {
        for i in 0..dim {
            let f: Vec<f64> = (0..ngp).map(|t| load_at[t][i]).collect();
            let mut part = vec![0.0; nf];
            udotf(&mut part, nf, ngp, &conv, &wb, &f);
            for p in 0..nf {
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

/// The pressure of a p element without `Pressure Stabilization`
/// (`LocalBulkMatrix` for `nb = 0`): only the nodal pressure dofs are kept,
/// each pressure dof of the functions `n..` gets the row and column of a
/// zero Dirichlet condition. `local` has `nf * dofs` rows, node-major.
pub fn pin_pressure(local: &mut LocalSystem, n: usize, dofs: usize) {
    let m = local.n;
    for p in n..m / dofs {
        let i = p * dofs + dofs - 1;
        local.force[i] = 0.0;
        for j in 0..m {
            local.stiff[i * m + j] = 0.0;
            local.stiff[j * m + i] = 0.0;
        }
        local.stiff[i * m + i] = 1.0;
        if let Some(mass) = local.mass.as_mut() {
            for j in 0..m {
                mass[i * m + j] = 0.0;
                mass[j * m + i] = 0.0;
            }
        }
    }
}

/// `LCondensate(N = nd, nb, dim, M, K, F)` of a steady run: the velocity
/// components of the `nb` bubble functions (local rows `(nd + p) * dofs + i`,
/// `i < dim`) are eliminated against the `nd * dofs` retained rows; the
/// bubble pressure rows are dropped. `local` has `(nd + nb) * dofs` rows,
/// node-major, and is left with `nd * dofs`.
pub fn lcondensate(local: &mut LocalSystem, nd: usize, nb: usize, dim: usize) -> Result<(), IncompressibleError> {
    let dofs = dim + 1;
    let m = local.n;
    if m != (nd + nb) * dofs {
        return Err(IncompressibleError::Unsupported(format!("{m} local rows for {nd} + {nb} functions")));
    }
    let (c, b) = condensed_rows(nd, nb, dim);
    let (nc, nbb) = (c.len(), b.len());
    let k = |i: usize, j: usize| local.stiff[i * m + j];
    let inv = invert((0..nbb * nbb).map(|x| k(b[x / nbb], b[x % nbb])).collect(), nbb)?;
    // Klb Kbb^-1, then the Schur complement on the retained rows.
    let mut klb_inv = vec![0.0; nc * nbb];
    for (r, &ci) in c.iter().enumerate() {
        for q in 0..nbb {
            klb_inv[r * nbb + q] = (0..nbb).map(|p| k(ci, b[p]) * inv[p * nbb + q]).sum();
        }
    }
    let fb: Vec<f64> = b.iter().map(|&i| local.force[i]).collect();
    let mut stiff = vec![0.0; nc * nc];
    let mut force = vec![0.0; nc];
    for (r, &ci) in c.iter().enumerate() {
        let row = &klb_inv[r * nbb..(r + 1) * nbb];
        force[r] = local.force[ci] - row.iter().zip(&fb).map(|(x, f)| x * f).sum::<f64>();
        for (s, &cj) in c.iter().enumerate() {
            let kbl: f64 = row.iter().zip(&b).map(|(x, &bi)| x * k(bi, cj)).sum();
            stiff[r * nc + s] = k(ci, cj) - kbl;
        }
    }
    local.n = nc;
    local.stiff = stiff;
    local.force = force;
    local.mass = None;
    local.damp = None;
    Ok(())
}

/// LCondensate's retained (`Cdofs`) and bubble velocity (`Bdofs`) rows.
fn condensed_rows(nd: usize, nb: usize, dim: usize) -> (Vec<usize>, Vec<usize>) {
    let dofs = dim + 1;
    let c = (0..nd * dofs).collect();
    let b = (0..nb).flat_map(|p| (0..dim).map(move |i| (nd + p) * dofs + i)).collect();
    (c, b)
}

/// `a^-1` (`n x n`, row-major) by Gauss-Jordan with partial pivoting
/// (`InvertMatrix`).
fn invert(mut a: Vec<f64>, nbb: usize) -> Result<Vec<f64>, IncompressibleError> {
    let mut inv = vec![0.0; nbb * nbb];
    for i in 0..nbb {
        inv[i * nbb + i] = 1.0;
    }
    for col in 0..nbb {
        let piv = (col..nbb)
            .max_by(|&x, &y| a[x * nbb + col].abs().total_cmp(&a[y * nbb + col].abs()))
            .expect("non-empty");
        if a[piv * nbb + col] == 0.0 {
            return Err(IncompressibleError::Unsupported("a singular bubble block (LCondensate)".into()));
        }
        for j in 0..nbb {
            a.swap(col * nbb + j, piv * nbb + j);
            inv.swap(col * nbb + j, piv * nbb + j);
        }
        let d = a[col * nbb + col];
        for j in 0..nbb {
            a[col * nbb + j] /= d;
            inv[col * nbb + j] /= d;
        }
        for r in 0..nbb {
            if r != col {
                let f = a[r * nbb + col];
                if f != 0.0 {
                    for j in 0..nbb {
                        a[r * nbb + j] -= f * a[col * nbb + j];
                        inv[r * nbb + j] -= f * inv[col * nbb + j];
                    }
                }
            }
        }
    }
    Ok(inv)
}

/// LCondensate of a transient run (`xprev`, `x` and the element's bubble
/// history given; "only works for the BDF(1) method"): `K += M / dt` and
/// `F += M xlprev / dt` over the whole local system, `xlprev` being `xprev`
/// (the previous timestep at the `nd` retained functions, node-major) with
/// `bubble_prev` (`dim` velocity components per bubble) at the bubble
/// velocity rows; then [`lcondensate`]. Returns the current bubble values
/// `Kbb^-1 (Fb - Kbl xl)` for the history, `xl` being the current iterate
/// at the retained rows.
#[allow(clippy::too_many_arguments)]
pub fn lcondensate_bdf1(
    local: &mut LocalSystem,
    nd: usize,
    nb: usize,
    dim: usize,
    dt: f64,
    xprev: &[f64],
    bubble_prev: &[f64],
    xl: &[f64],
) -> Result<Vec<f64>, IncompressibleError> {
    let m = local.n;
    let (c, b) = condensed_rows(nd, nb, dim);
    if m != (nd + nb) * (dim + 1) || xprev.len() != c.len() || xl.len() != c.len() || bubble_prev.len() != b.len() {
        return Err(IncompressibleError::Unsupported("LCondensate: inconsistent local sizes".into()));
    }
    let Some(mass) = local.mass.as_ref() else {
        return Err(IncompressibleError::Unsupported("LCondensate: a transient run without a mass matrix".into()));
    };
    let mut xlprev = vec![0.0; m];
    for (&i, &x) in c.iter().zip(xprev) {
        xlprev[i] = x;
    }
    for (&i, &x) in b.iter().zip(bubble_prev) {
        xlprev[i] = x;
    }
    for i in 0..m {
        let mx: f64 = (0..m).map(|j| mass[i * m + j] * xlprev[j]).sum();
        local.force[i] += mx / dt;
        for j in 0..m {
            local.stiff[i * m + j] += mass[i * m + j] / dt;
        }
    }
    let nbb = b.len();
    let inv = invert((0..nbb * nbb).map(|x| local.stiff[b[x / nbb] * m + b[x % nbb]]).collect(), nbb)?;
    let r: Vec<f64> = b
        .iter()
        .map(|&bi| local.force[bi] - c.iter().zip(xl).map(|(&cj, x)| local.stiff[bi * m + cj] * x).sum::<f64>())
        .collect();
    let cur = (0..nbb).map(|p| (0..nbb).map(|q| inv[p * nbb + q] * r[q]).sum()).collect();
    lcondensate(local, nd, nb, dim)?;
    Ok(cur)
}
