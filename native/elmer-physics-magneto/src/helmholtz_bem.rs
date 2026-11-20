//! `HelmholtzBEM.F90` (fem/src/modules): the 2D Helmholtz equation solved
//! by the boundary element method, `Procedure = "HelmholtzBEM"
//! "HelmholtzBEMSolver"`.
//!
//! The boundary elements the solver is active on (those a boundary
//! condition puts into a body with `Body Id`) carry the unknowns: per node
//! either the potential (`<Variable> 1|2` or `Potential 1|2` given) or the
//! flux (`Flux 1|2`). `IntegrateMatrix` collocates at the boundary nodes
//! with the 2D Green's function `(J0 - i Y0) / (4 i)` (the plugin's own
//! `Bessel`), the outward `NormalVector(..., Check)`, the 4-point line rule;
//! the dense system is solved with `ZGETRF`/`ZGETRS` (`SolveLapack_cmplx`,
//! [`crate::zlapack`]). `ComputePotential` then evaluates the potential at
//! the other nodes of the solver's mesh, and the exported `Flux` gets the
//! boundary fluxes. The norm is `SQRT(SUM(ABS(Potential)**2)) / nNodes`.
//!
//! Refused: 3D, the iterative `FullIterSolver`, boundary elements other
//! than lines, boundaries between two bodies.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use num_complex::Complex64 as C;

use elmer_core::{BasisEval, ElemId, ElementType, ParamMap, Params};
use elmer_fe::metric::MetricStatus;
use elmer_itermethods::IterScalar;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

const CALLER: &str = "HelmholtzBEMSolver";

pub const PROC: (&str, &str) = ("HelmholtzBEM", "HelmholtzBEMSolver");

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROC.0, PROC.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(HelmholtzBem))
}

struct HelmholtzBem;

fn unsupported(m: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {}", m.into()))
}

#[inline]
fn mul(a: C, b: C) -> C {
    C::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

/// `a * z` for a real `a` (componentwise, as gfortran lowers a real
/// operand).
#[inline]
fn rmul(a: f64, z: C) -> C {
    C::new(a * z.re, a * z.im)
}

/// The plugin's `Bessel(x, j0, y0, dj0, dy0)`, with its literal constants.
#[allow(clippy::approx_constant)]
fn bessel(x: f64) -> Result<(f64, f64, f64, f64), SolverError> {
    const MAXROUNDS: usize = 1000;
    const POLYLIMIT: f64 = 10.0;
    const ACCURACY: f64 = 1.0e-8;
    const GAMMA: f64 = 0.577_215_664_901_532_860_606_512_090_082;
    const A: [f64; 7] = [0.79788456, -0.00000077, -0.00552740, 0.00009512, 0.00137237, -0.00072805, 0.00014476];
    const B: [f64; 7] = [-0.78539816, -0.04166397, -0.00003954, 0.00262573, -0.00054125, -0.00029333, 0.00013558];
    const CC: [f64; 7] = [0.79788456, 0.00000156, 0.01659667, 0.00017105, -0.00249511, 0.00113653, -0.0020033];
    const D: [f64; 7] = [-2.35619449, 0.12499612, 0.00005650, -0.00637879, 0.00074348, 0.00079824, -0.00029166];
    if x.abs() > POLYLIMIT {
        let t = 3.0 / x;
        let (mut p, mut f) = (x, 0.0);
        for k in 0..7 {
            f += A[k] * t.powf(k as f64);
            p += B[k] * t.powf(k as f64);
        }
        let j0 = f * p.cos() / x.sqrt();
        let y0 = f * p.sin() / x.sqrt();
        let (mut p, mut f) = (x, 0.0);
        for k in 0..7 {
            f += CC[k] * t.powf(k as f64);
            p += D[k] * t.powf(k as f64);
        }
        let dj0 = -f * p.cos() / x.sqrt();
        let dy0 = -f * p.sin() / x.sqrt();
        return Ok((j0, y0, dj0, dy0));
    }
    let (mut j0, mut y0, mut dj0, mut dy0) = (1.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut hk = 0.0f64;
    let mut done = false;
    for k in 1..=MAXROUNDS {
        let kf = k as f64;
        hk += 1.0 / kf;
        let mut res = 1.0f64;
        for i in 1..=k {
            let q = x / (2.0 * i as f64);
            res *= q * q;
        }
        let sk: i64 = if k % 2 == 0 { 1 } else { -1 };
        let sk1 = -sk;
        j0 += sk as f64 * res;
        y0 += sk1 as f64 * hk * res;
        dj0 += (sk * k as i64) as f64 / (0.5 * x) * res;
        dy0 += sk1 as f64 * hk * kf / (0.5 * x) * res;
        if (kf / (0.5 * x) * res).abs() < ACCURACY {
            done = true;
            break;
        }
    }
    if !done {
        return Err(SolverError::fatal(CALLER, "Error in evaluating Bessel functions"));
    }
    let lg = (0.5 * x).ln() + GAMMA;
    y0 += lg * j0;
    y0 = y0 * 2.0 / std::f64::consts::PI;
    dy0 = dy0 + (1.0 / x) * j0 + lg * dj0;
    dy0 = dy0 * 2.0 / std::f64::consts::PI;
    Ok((j0, y0, dj0, dy0))
}

/// `Green(2, k, x, y, z, W, GradW)`: `W = (J0 - i Y0) / (i 4)`,
/// `dW/dR = k (dJ0 - i dY0) / (i 4)`, `GradW = (x, y, z) dW/dR / R`.
fn green(k: f64, x: f64, y: f64, z: f64) -> Result<(C, [C; 3]), SolverError> {
    let r = ((x * x + y * y) + z * z).sqrt();
    let (j0, y0, dj0, dy0) = bessel(k * r)?;
    let i4 = C::new(0.0, 4.0);
    // J0 - i*Y0: a real minus (0*Y0, Y0).
    let num = C::new(j0 - 0.0 * y0, -y0);
    let w = num.fortran_div(i4);
    let dnum = C::new(k * dj0, -(k * dy0));
    let dwdr = dnum.fortran_div(i4);
    let g = |c: f64| C::new(c * dwdr.re / r, c * dwdr.im / r);
    Ok((w, [g(x), g(y), g(z)]))
}

/// `NormalVector(Boundary, Nodes, u, v, .TRUE.)` of a 2D line with
/// `CheckNormalDirection` (single parent: away from it, unless it is the
/// `Normal Target Body`).
fn line_normal(ctx: &SolverCtx<'_>, e: ElemId, u: f64) -> Result<[f64; 3], SolverError> {
    let el = ctx.mesh.element(e);
    let xs: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let n = xs.len();
    let fe = |m: elmer_fe::FeError| SolverError::fatal("NormalVector", m.to_string());
    let def = elmer_fe::element_type(el.etype.0).map_err(fe)?;
    let mut dl = vec![[0.0; 3]; n];
    elmer_fe::basis::nodal_first_derivatives(def, u, 0.0, 0.0, &mut dl);
    let dot = |a: usize| (0..n).fold(0.0, |s, i| s + dl[i][0] * xs[i][a]);
    let (dxdu, dydu) = (dot(0), dot(1));
    let det = dxdu * dxdu + dydu * dydu;
    if det <= 0.0 {
        return Ok([0.0; 3]);
    }
    let det = 1.0 / det.sqrt();
    let mut normal = [-dydu * det, dxdu * det, 0.0];
    // CheckNormalDirection: the boundary point at u = 0.
    let mut b = vec![0.0; n];
    elmer_fe::basis::nodal_basis(def, 0.0, 0.0, 0.0, &mut b);
    let pt: [f64; 3] = std::array::from_fn(|a| (0..n).fold(0.0, |s, i| s + b[i] * xs[i][a]));
    let parent = match el.parents {
        [Some(_), Some(_)] => return Err(unsupported("a boundary between two bodies")),
        [Some(p), None] | [None, Some(p)] => p,
        [None, None] => return Ok(normal),
    };
    let out_body = ctx.bc_params(e).and_then(|bc| bc.integer("normal target body")).map_or(-1, i64::from);
    let pel = ctx.mesh.element(parent);
    let pdef = elmer_fe::element_type(pel.etype.0).map_err(fe)?;
    let (pu, pv, pw) = match pel.etype.0 / 100 {
        2 | 4 | 8 => (0.0, 0.0, 0.0),
        3 | 7 => (1.0 / 3.0, 1.0 / 3.0, 0.0),
        5 => (0.25, 0.25, 0.25),
        6 => (0.0, 0.0, 1.0 / 3.0),
        f => return Err(SolverError::fatal("CheckNormalDirection", format!("Invalid elementcode for parent element {f}"))),
    };
    let pn = pel.nodes.len();
    let mut pb = vec![0.0; pn];
    elmer_fe::basis::nodal_basis(pdef, pu, pv, pw, &mut pb);
    let pxs: Vec<[f64; 3]> = pel.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let dc: [f64; 3] = std::array::from_fn(|a| (0..pn).fold(0.0, |s, i| s + pb[i] * pxs[i][a]) - pt[a]);
    let body = ctx.sim.model.body_of_tag(pel.tag).map_or(0, |b| b as i64);
    let s = (0..3).fold(0.0, |s, a| s + dc[a] * normal[a]);
    if (s > 0.0 && body != out_body) || (s <= 0.0 && body == out_body) {
        normal = normal.map(|c| -c);
    }
    Ok(normal)
}

struct Line {
    e: ElemId,
    nodes: Vec<u32>,
    xyz: Vec<[f64; 3]>,
}

/// One quadrature point of a boundary line: `s * detJ`, the basis, the
/// point and the normal.
struct Qp {
    s: f64,
    basis: Vec<f64>,
    p: [f64; 3],
    normal: [f64; 3],
}

fn points(ctx: &SolverCtx<'_>, l: &Line, status: &mut MetricStatus) -> Result<Vec<Qp>, SolverError> {
    let n = l.nodes.len();
    let rule = elmer_fe::quadrature::gauss_points(ElementType(202), Some(4)).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let def = elmer_fe::element_type(202).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let mut out = Vec::with_capacity(rule.len());
    for t in 0..rule.len() {
        let mut be = BasisEval::default();
        elmer_fe::info::element_info(def, &l.xyz, [rule.u[t], rule.v[t], rule.w[t]], 2, status, &mut be)
            .map_err(|m| SolverError::fatal("ElementInfo", m.to_string()))?;
        let s = rule.s[t] * be.det_j;
        let normal = line_normal(ctx, l.e, rule.u[t])?;
        let basis = be.basis[..n].to_vec();
        let p: [f64; 3] = std::array::from_fn(|a| {
            let mut acc = 0.0;
            for i in 0..n {
                acc += l.xyz[i][a] * basis[i];
            }
            acc
        });
        out.push(Qp { s, basis, p, normal });
    }
    Ok(out)
}

impl Solver for HelmholtzBem {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.comm.size() > 1 {
            return Err(SolverError::fatal(CALLER, "Unfortunately this module does not support MPI!"));
        }
        let sp = ctx.params().clone();
        let cs = ctx.simulation_params().string("coordinate system").unwrap_or("cartesian").to_ascii_lowercase();
        if !cs.trim().starts_with("cartesian") {
            return Err(SolverError::fatal(CALLER, "This solver is implemented only for cartesian coordinates!"));
        }
        if ctx.mesh.dim() != 2 {
            return Err(unsupported("the 3D Green's function"));
        }
        match sp.string("linear system solver").map(|s| s.trim().to_ascii_lowercase()).as_deref() {
            Some("direct") => {}
            Some("iterative") => return Err(unsupported("the iterative FullIterSolver")),
            _ => return Err(SolverError::fatal("SolveFull", "Unknown solver type.")),
        }
        let var_name = ctx.variable_name().unwrap_or("").to_ascii_lowercase();
        let (vdofs, vperm) = {
            let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
            (v.dofs, v.perm.map(|p| p.to_vec()))
        };
        if vdofs != 2 {
            return Err(SolverError::fatal(CALLER, "the variable needs two components"));
        }
        let nnodes = ctx.mesh.n_nodes();
        let nbulk = ctx.mesh.n_bulk();
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let first = *active.first().ok_or_else(|| SolverError::fatal(CALLER, "no active elements"))?;
        let k = crate::calcfields2d::angular_frequency(ctx, &sp, first)?.unwrap_or(0.0);

        // The boundary lines and their nodes, numbered in element order.
        let mut lines = Vec::new();
        let mut pperm = vec![usize::MAX; nnodes];
        let mut bperm = Vec::new();
        for &e in &active {
            if (e.0 as usize) < nbulk {
                continue;
            }
            let el = ctx.mesh.element(e);
            if el.etype.0 == 101 {
                continue;
            }
            if el.etype.0 != 202 {
                return Err(unsupported(format!("boundary elements of type {}", el.etype.0)));
            }
            for &k in el.nodes {
                if pperm[k as usize] == usize::MAX {
                    pperm[k as usize] = bperm.len();
                    bperm.push(k as usize);
                }
            }
            let nodes = el.nodes.to_vec();
            let xyz = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            lines.push(Line { e, nodes, xyz });
        }
        let nb = bperm.len();
        ctx.info(CALLER, &format!("Number of nodes on boundaries: {nb}"), 4);

        // Known potentials or fluxes.
        let mut load = vec![C::new(0.0, 0.0); nb];
        let mut known = vec![false; nb];
        for l in &lines {
            let Some(bc) = ctx.bc_params(l.e) else { continue };
            let n = l.nodes.len();
            let get = |a: &str, b: &str| -> Result<(Vec<f64>, bool), SolverError> {
                if let Some(v) = ctx.real_at_nodes(bc, a, l.e)? {
                    return Ok((v, true));
                }
                Ok(match ctx.real_at_nodes(bc, b, l.e)? {
                    Some(v) => (v, true),
                    None => (vec![0.0; n], false),
                })
            };
            let (p1, _) = get(&format!("{var_name} 1"), "potential 1")?;
            let (p2, got) = get(&format!("{var_name} 2"), "potential 2")?;
            let (v1, v2) = if got {
                (p1, p2)
            } else {
                let f1 = ctx.real_at_nodes(bc, "flux 1", l.e)?.unwrap_or_else(|| vec![0.0; n]);
                let f2 = ctx.real_at_nodes(bc, "flux 2", l.e)?.unwrap_or_else(|| vec![0.0; n]);
                (f1, f2)
            };
            for (j, &node) in l.nodes.iter().enumerate() {
                let q = pperm[node as usize];
                load[q] = C::new(v1[j], v2[j]);
                known[q] = got;
            }
        }

        // IntegrateMatrix.
        let mut status = MetricStatus::default();
        let mut mat = vec![C::new(0.0, 0.0); nb * nb];
        let mut diag = vec![C::new(0.0, 0.0); nb];
        let mut force = vec![C::new(0.0, 0.0); nb];
        let coord = |k: usize| ctx.mesh.coord(k as u32);
        for l in &lines {
            for qp in points(ctx, l, &mut status)? {
                let s = qp.s;
                for p in 0..nb {
                    let c = coord(bperm[p]);
                    let (x, y, z) = (qp.p[0] - c[0], qp.p[1] - c[1], qp.p[2] - c[2]);
                    let (g, grad) = green(k, x, y, z)?;
                    let mut dgdn = C::new(0.0, 0.0);
                    for a in 0..3 {
                        dgdn += rmul(qp.normal[a], grad[a]);
                    }
                    for (j, &node) in l.nodes.iter().enumerate() {
                        let q = pperm[node as usize];
                        let bj = qp.basis[j];
                        if known[q] {
                            if p != q {
                                force[p] -= mul(rmul(bj, rmul(s, load[q])), dgdn);
                            }
                            mat[p + q * nb] -= rmul(s * bj, g);
                        } else {
                            force[p] += mul(rmul(bj, rmul(s, load[q])), g);
                            if p != q {
                                mat[p + q * nb] += rmul(s * bj, dgdn);
                            }
                        }
                        // R = -k*SIN(k*x)*G*Normal(1) (- COS(k*x)*dGdN).
                        let mut r = rmul(qp.normal[0], rmul(k * (k * x).sin(), g));
                        r = C::new(-r.re, -r.im);
                        if p != q {
                            r -= rmul((k * x).cos(), dgdn);
                        }
                        diag[p] -= rmul(s * bj, r);
                    }
                }
            }
        }
        for i in 0..nb {
            if known[i] {
                force[i] -= mul(load[i], diag[i]);
            } else {
                mat[i + i * nb] = diag[i];
            }
        }

        // SolveFull: SolveLapack_cmplx.
        let ipiv = crate::zlapack::zgetrf(nb, &mut mat).map_err(|m| SolverError::fatal(CALLER, m))?;
        crate::zlapack::zgetrs(nb, &mat, &ipiv, &mut force);
        let (mut pot, mut flx) = (vec![C::new(0.0, 0.0); nb], vec![C::new(0.0, 0.0); nb]);
        for i in 0..nb {
            if known[i] {
                flx[i] = force[i];
                pot[i] = load[i];
            } else {
                flx[i] = load[i];
                pot[i] = force[i];
            }
        }

        // ComputePotential at the other nodes in the solver's perm.
        let mut potential = vec![C::new(0.0, 0.0); nnodes];
        for i in 0..nb {
            potential[bperm[i]] = pot[i];
        }
        let in_perm = |i: usize| vperm.as_ref().is_none_or(|p| p[i] != elmer_core::NO_DOF);
        for l in &lines {
            for qp in points(ctx, l, &mut status)? {
                let s = qp.s;
                for i in 0..nnodes {
                    if pperm[i] != usize::MAX || !in_perm(i) {
                        continue;
                    }
                    let c = coord(i);
                    let (g, grad) = green(k, qp.p[0] - c[0], qp.p[1] - c[1], qp.p[2] - c[2])?;
                    let mut dgdn = C::new(0.0, 0.0);
                    for a in 0..3 {
                        dgdn += rmul(qp.normal[a], grad[a]);
                    }
                    for (j, &node) in l.nodes.iter().enumerate() {
                        let q = pperm[node as usize];
                        let t = mul(pot[q], dgdn) - mul(flx[q], g);
                        potential[i] -= rmul(s * qp.basis[j], t);
                    }
                }
            }
        }

        // The fields and the norm.
        {
            let v = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
            v.values.fill(0.0);
            for (i, z) in potential.iter().enumerate() {
                let j = match &vperm {
                    Some(p) if p[i] == elmer_core::NO_DOF => continue,
                    Some(p) => p[i] as usize,
                    None => i,
                };
                v.values[2 * j] = z.re;
                v.values[2 * j + 1] = z.im;
            }
        }
        if let Some(f) = ctx.variable_mut("flux") {
            f.values.fill(0.0);
            for i in 0..nb {
                let k = bperm[i];
                if 2 * k + 1 < f.values.len() {
                    f.values[2 * k] = flx[i].re;
                    f.values[2 * k + 1] = flx[i].im;
                }
            }
        }
        let mut s = 0.0;
        for z in &potential {
            let a = z.norm();
            s += a * a;
        }
        let norm = s.sqrt() / nnodes as f64;
        let id = ctx.sim.slots.iter().find(|sl| sl.id == ctx.solver_id()).and_then(|sl| sl.variable);
        if let Some(id) = id {
            let v = ctx.sim.vars.by_id_mut(id);
            v.conv.norm = norm;
            v.values_changed = true;
        }
        Ok(())
    }
}
