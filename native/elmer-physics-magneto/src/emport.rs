//! `EMPort.F90` (fem/src/modules) `EMPortSolver`: the eigenmodes of a 2D
//! waveguide port cross-section with the scalar potential (`Use Potential
//! = True`), their power normalization, the non-homogeneous (unit voltage)
//! solution, the port impedance and the nodal field `EF2D`
//! (`Calculate Nodal Field`).
//!
//! [`init_section`] is `EMPortSolver_Init0`: the complex `Eport` variable
//! of `Element = "n:2 e:1 ..."` (`Ez` and `V` per node, the transverse
//! field on the edges), the eigen analysis keywords and, on `Port Ground` /
//! `Unit Voltage` boundaries, the `Eport re|im {n} 1|2` and `{e}` values.
//!
//! The solve:
//! - `LocalMatrix_With_Potentials` on `GaussPoints(EdgeBasis, PReference
//!   =Piola)` with the Piola basis; complex `Nu`, `Eps`.
//! - `DefaultDirichletBCs` with the `{n}` construct and the `{e}` zeros.
//! - `Eigen System Shift Automatic`: `-(Omega SQRT(maxeps maxmu))**2`, then
//!   `DefaultSolve` (the complex shifted eigen system of elmer-solver).
//! - `Beta = SQRT(-lambda)`, normalization by the port power, the
//!   non-homogeneous mode (`K - lambda M` solved and made orthogonal to the
//!   eigenvector), the impedance, `EMPortPost`.
//!
//! Ported: 2D, `Use Potential`, first-order Piola triangles, no port
//! indices. Refused otherwise.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use num_complex::Complex64 as C;

use elmer_core::{ElemId, ElementType, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_itermethods::IterScalar;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::case::has_prefix;

const CALLER: &str = "EMPortSolver";

pub const PROC: (&str, &str) = ("EMPort", "EMPortSolver");

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROC.0, PROC.1, |_| Ok(Box::new(EmPort)));
}

pub fn is_section(p: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(p)
        .is_some_and(|(l, e)| l.trim_start_matches("./").eq_ignore_ascii_case(PROC.0) && e.eq_ignore_ascii_case(PROC.1))
}

fn unsupported(m: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {}", m.into()))
}

fn add_new(p: &mut ParamMap, key: &str, v: Value) {
    if p.get(key).is_none() {
        p.insert(key, v);
    }
}

/// `EMPortSolver_Init0` (no primary `VectorHelmholtzSolver` keywords are
/// copied: a section with one is refused at the solve).
pub fn init_section(p: &mut ParamMap, bcs: &mut [ParamMap]) {
    let use_v = p.logical("use potential") == Some(true);
    add_new(p, "Linear System Complex", Value::Logical(true));
    add_new(p, "Variable DOFs", Value::Integer(vec![2]));
    add_new(p, "Eigen Analysis", Value::Logical(true));
    add_new(p, "Nonlinear System Max Iterations", Value::Integer(vec![1]));
    if p.get("element").is_none() {
        let quadratic = p.logical("quadratic approximation") == Some(true);
        let second = p.logical("second kind basis") == Some(true);
        let piola = p.logical("use piola transform") == Some(true) || quadratic || second;
        let e = if use_v {
            if piola {
                "n:2 e:1 -quad_face b:2 -quad b:2 -brick b:3"
            } else {
                "n:2 e:1"
            }
        } else if quadratic && second {
            "n:1 e:3 -tri b:3 -tri_face b:3"
        } else if quadratic {
            "n:1 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2"
        } else if second {
            "n:1 e:2"
        } else if piola {
            "n:1 e:1 -quad_face b:2 -quad b:2 -brick b:3"
        } else {
            "n:1 e:1"
        };
        p.insert("Element", Value::Str(e.into()));
    }
    for bc in bcs.iter_mut() {
        let ground = bc.logical("port ground") == Some(true);
        let unit = bc.logical("unit voltage") == Some(true);
        if !(ground || unit) {
            continue;
        }
        let c = if ground { 0.0 } else { 1.0 };
        if use_v {
            bc.insert("Eport re {n} 1", Value::real(0.0));
            bc.insert("Eport re {n} 2", Value::real(c));
            bc.insert("Eport im {n} 1", Value::real(0.0));
            bc.insert("Eport im {n} 2", Value::real(0.0));
        } else {
            bc.insert("Eport re", Value::real(0.0));
            bc.insert("Eport im", Value::real(0.0));
        }
        bc.insert("Eport re {e}", Value::real(0.0));
        bc.insert("Eport im {e}", Value::real(0.0));
    }
    add_new(p, "Variable", Value::Str("Eport[Eport re:1 Eport im:1]".into()));
    p.insert("Linear System Refactorize", Value::Logical(true));
    add_new(p, "Skip Compute Nonlinear Change", Value::Logical(true));
    add_new(p, "Skip Compute Steady State Change", Value::Logical(true));
    add_new(p, "Variable Output", Value::Logical(false));
    add_new(p, "post: Skip Compute Nonlinear Change", Value::Logical(true));
    add_new(p, "post: Linear System Complex", Value::Logical(false));
    add_new(p, "post: Variable Output", Value::Logical(false));
    add_new(p, "Eigen System Sorting", Value::Str("smallest real part".into()));
    add_new(p, "Eigenvector Normalization by Power", Value::Logical(true));
    if p.get("eigen system shift").is_none() {
        add_new(p, "Eigen System Shift Automatic", Value::Logical(true));
    }
    add_new(p, "Linear System Nullify Guess", Value::Logical(true));
}

#[inline]
fn mul(a: C, b: C) -> C {
    C::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

/// A real times a complex, componentwise (gfortran's lowering).
#[inline]
fn rmul(a: f64, z: C) -> C {
    C::new(a * z.re, a * z.im)
}

#[inline]
fn div(a: C, b: C) -> C {
    a.fortran_div(b)
}

/// glibc `csqrt` for finite, moderately sized arguments.
fn csqrt(z: C) -> C {
    let (re, im) = (z.re, z.im);
    if im == 0.0 {
        return if re < 0.0 { C::new(0.0, (-re).sqrt().copysign(im)) } else { C::new(re.sqrt().abs(), im) };
    }
    if re == 0.0 {
        let r = if im.abs() >= 2.0 * f64::MIN_POSITIVE { (0.5 * im.abs()).sqrt() } else { 0.5 * (2.0 * im.abs()).sqrt() };
        return C::new(r, r.copysign(im));
    }
    let d = re.hypot(im);
    let (r, s);
    if re > 0.0 {
        r = (0.5 * (d + re)).sqrt();
        s = 0.5 * (im / r);
    } else {
        s = (0.5 * (d - re)).sqrt();
        r = (0.5 * (im / s)).abs();
    }
    C::new(r, s.copysign(im))
}

/// `SUM(a(1:3) * b(1:3))`.
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2]
}

/// `GaussPointsTriangle(3, PReferenceElement=.TRUE.)`.
fn rule() -> Result<elmer_fe::quadrature::GaussRule, SolverError> {
    let mut p = elmer_fe::quadrature::gauss_points(ElementType(303), Some(3)).map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
    let sqrt3 = 3.0f64.sqrt();
    for i in 0..p.len() {
        let (u, v, s) = (p.u[i], p.v[i], p.s[i]);
        p.u[i] = -1.0 + 2.0 * u + v;
        p.v[i] = sqrt3 * v;
        p.s[i] = sqrt3 * 2.0 * s;
    }
    Ok(p)
}

/// `ListGetElementComplex(h, Basis, Element)` of a material keyword (and
/// its `im`): `None` when neither is given.
fn complex_coeff(ctx: &SolverCtx<'_>, m: Option<&ParamMap>, key: &str, e: ElemId, basis: &[f64]) -> Result<Option<C>, SolverError> {
    let Some(m) = m else { return Ok(None) };
    let one = |k: &str| -> Result<Option<f64>, SolverError> {
        Ok(match m.get(k) {
            None => None,
            Some(Value::Real { data, .. }) if data.len() == 1 => Some(data[0]),
            Some(_) => {
                let v = ctx.real_at_nodes(m, k, e)?.unwrap_or_default();
                let mut s = 0.0;
                for (b, x) in basis.iter().zip(&v) {
                    s += b * x;
                }
                Some(s)
            }
        })
    };
    let re = one(key)?;
    let im = one(&format!("{key} im"))?;
    if re.is_none() && im.is_none() {
        return Ok(None);
    }
    Ok(Some(C::new(re.unwrap_or(0.0), im.unwrap_or(0.0))))
}

struct Pt {
    weight: f64,
    basis: Vec<f64>,
    db: Vec<[f64; 3]>,
    w: Vec<[f64; 3]>,
    curl: Vec<[f64; 3]>,
}

fn points(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Vec<Pt>, SolverError> {
    let el = ctx.mesh.element(e);
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let gindex: Vec<i64> = el.nodes.iter().map(|&k| k as i64 + 1).collect();
    let r = rule()?;
    let mut out = Vec::with_capacity(r.len());
    for t in 0..r.len() {
        let (ev, basis, db) = crate::piola::edge_info_cdim(el.etype, &xyz, &gindex, [r.u[t], r.v[t], r.w[t]], 1, ctx.mesh.dim())?;
        out.push(Pt { weight: r.s[t] * ev.det, basis, db, w: ev.edge, curl: ev.curl });
    }
    Ok(out)
}

/// The (complex) eigenvector `mode` at element `e`'s DOFs
/// (`GetScalarLocalEigenmode`, real and imaginary parts).
fn local_mode(ctx: &SolverCtx<'_>, e: ElemId, mode: usize) -> Result<Vec<C>, SolverError> {
    let mut slots = Vec::new();
    ctx.element_slots(e, &mut slots);
    let id = ctx.sim.slots.iter().find(|s| s.id == ctx.solver_id()).and_then(|s| s.variable).ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
    let ev = &ctx.sim.vars.by_id(id).eigen_vectors;
    let vec = ev.get(mode).ok_or_else(|| SolverError::fatal(CALLER, "no such eigenmode"))?;
    Ok(slots.iter().map(|&s| if s == NO_DOF { C::new(0.0, 0.0) } else { vec[s as usize] }).collect())
}

struct EmPort;

impl EmPort {
    /// `LocalMatrix_With_Potentials`.
    #[allow(clippy::too_many_arguments)]
    fn local(
        ctx: &SolverCtx<'_>,
        e: ElemId,
        nd: usize,
        mu0inv: f64,
        eps0: f64,
        omega: f64,
        maxmu: &mut f64,
        maxeps: &mut f64,
    ) -> Result<(Vec<C>, Vec<C>), SolverError> {
        let el = ctx.mesh.element(e);
        let n = el.nodes.len();
        let ndofs = 2;
        let np = n * ndofs;
        let vdofs = nd - np;
        let material = ctx.body_params(e).and_then(|b| b.material);
        let mut stiff = vec![C::new(0.0, 0.0); nd * nd];
        let mut mass = vec![C::new(0.0, 0.0); nd * nd];
        let om2 = omega * omega;
        for pt in points(ctx, e)? {
            if pt.w.len() != vdofs {
                return Err(SolverError::fatal(CALLER, "edge basis and element DOFs differ"));
            }
            let weight = pt.weight;
            let nu = match complex_coeff(ctx, material, "relative reluctivity", e, &pt.basis)? {
                Some(v) => rmul(mu0inv, v),
                None => C::new(mu0inv, 0.0),
            };
            let eps = match complex_coeff(ctx, material, "relative permittivity", e, &pt.basis)? {
                Some(v) => rmul(eps0, v),
                None => C::new(eps0, 0.0),
            };
            *maxmu = maxmu.max(div(C::new(1.0, 0.0), nu).re);
            *maxeps = maxeps.max(eps.re);
            let (b, db, w, cw) = (&pt.basis, &pt.db, &pt.w, &pt.curl);
            let ix = |i: usize, j: usize| i * nd + j;
            for p in 0..n {
                let i = p * ndofs;
                for q in 0..n {
                    let j = q * ndofs;
                    let t = rmul(dot3(&db[q], &db[p]), nu) - rmul(b[p], rmul(b[q], rmul(om2, eps)));
                    stiff[ix(i, j)] += rmul(weight, t);
                }
                for q in 0..vdofs {
                    let j = np + q;
                    stiff[ix(i, j)] -= rmul(weight, rmul(dot3(&w[q], &db[p]), nu));
                }
                for q in 0..n {
                    let j = q * ndofs + 1;
                    stiff[ix(i, j)] += rmul(weight, rmul(dot3(&db[q], &db[p]), nu));
                }
                let i = p * ndofs + 1;
                for q in 0..n {
                    let j = q * ndofs + 1;
                    let t = div(rmul(weight * om2, eps), nu);
                    stiff[ix(i, j)] -= rmul(b[q], rmul(b[p], t));
                }
                for q in 0..vdofs {
                    let j = np + q;
                    let x = dot3(&w[q], &db[p]) * weight;
                    stiff[ix(i, j)] += C::new(x, 0.0);
                }
            }
            for p in 0..vdofs {
                let i = np + p;
                for q in 0..vdofs {
                    let j = np + q;
                    let t = rmul(dot3(&cw[q], &cw[p]), nu) - rmul(dot3(&w[q], &w[p]), rmul(om2, eps));
                    stiff[ix(i, j)] += rmul(weight, t);
                    mass[ix(i, j)] += rmul(dot3(&w[q], &w[p]), rmul(weight, nu));
                }
                for q in 0..n {
                    let j = q * ndofs;
                    let s = dot3(&db[q], &w[p]);
                    mass[ix(i, j)] -= rmul(s, rmul(weight, nu));
                    let j = q * ndofs + 1;
                    stiff[ix(i, j)] += rmul(s, rmul(weight * om2, eps));
                    mass[ix(i, j)] -= rmul(s, rmul(weight, nu));
                }
            }
        }
        Ok((stiff, mass))
    }

    /// `CalculatePortPower` with `UseV`: `P += Nu/Omega (Beta (EF, EF*) -
    /// i (EF, gradEz*)) weight` and `E2 += (EF, EF*) weight`.
    #[allow(clippy::too_many_arguments)]
    fn port_power(ctx: &SolverCtx<'_>, e: ElemId, mode: usize, beta: C, mu0inv: f64, omega: f64, p: &mut C, e2: &mut C) -> Result<(), SolverError> {
        let el = ctx.mesh.element(e);
        let n = el.nodes.len();
        let (ndofs, np) = (2, 2 * n);
        let f = local_mode(ctx, e, mode)?;
        let vdofs = f.len() - np;
        let material = ctx.body_params(e).and_then(|b| b.material);
        let im = C::new(0.0, 1.0);
        for pt in points(ctx, e)? {
            let weight = pt.weight;
            let nu = match complex_coeff(ctx, material, "relative reluctivity", e, &pt.basis)? {
                Some(v) => rmul(mu0inv, v),
                None => C::new(mu0inv, 0.0),
            };
            let mut ef = [C::new(0.0, 0.0); 3];
            for i in 0..vdofs {
                for c in 0..3 {
                    ef[c] += rmul(pt.w[i][c], f[np + i]);
                }
            }
            for c in 0..3 {
                let (mut sr, mut si) = (0.0, 0.0);
                for k in 0..n {
                    sr += f[k * ndofs + 1].re * pt.db[k][c];
                    si += f[k * ndofs + 1].im * pt.db[k][c];
                }
                ef[c] -= C::new(sr, si);
            }
            let mut grad = [C::new(0.0, 0.0); 3];
            for c in 0..3 {
                let (mut sr, mut si) = (0.0, 0.0);
                for k in 0..n {
                    sr += f[k * ndofs].re * pt.db[k][c];
                    si += f[k * ndofs].im * pt.db[k][c];
                }
                grad[c] += C::new(sr, si);
            }
            for c in 0..3 {
                grad[c] = mul(mul(grad[c], im), beta);
            }
            let mut ee = C::new(0.0, 0.0);
            let mut eg = C::new(0.0, 0.0);
            for c in 0..3 {
                ee += mul(ef[c], ef[c].conj());
                eg += mul(ef[c], grad[c].conj());
            }
            let t = mul(beta, ee) - mul(im, eg);
            *p += rmul(weight, mul(div(nu, C::new(omega, 0.0)), t));
            *e2 += rmul(weight, ee);
        }
        Ok(())
    }
}

impl Solver for EmPort {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let sp = ctx.params().clone();
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::fatal(CALLER, "the edge element DOFs (no DofsBackend)"));
        }
        let cs = ctx.simulation_params().string("coordinate system").unwrap_or("cartesian").to_ascii_lowercase();
        if !cs.trim().starts_with("cartesian") {
            return Err(SolverError::fatal(CALLER, "Implemented only for Cartesian problems!"));
        }
        if ctx.primary().map(|v| v.dofs) != Some(2) {
            return Err(SolverError::fatal(CALLER, "Complex field, specify two DOFs"));
        }
        if sp.logical("use potential") != Some(true) {
            return Err(unsupported("the formulation without the potential (Use Potential = False)"));
        }
        if ctx.mesh.dim() != 2 {
            return Err(unsupported("3D"));
        }
        if sp.logical("quadratic approximation") == Some(true) || sp.logical("second kind basis") == Some(true) || sp.logical("use piola transform") != Some(true) {
            return Err(unsupported("other than the first-order Piola basis"));
        }
        for bc in &ctx.sim.model.boundary_conditions {
            for k in ["port index", "port type", "port beta parent", "calculate impedance"] {
                if has_prefix(bc, k) {
                    return Err(unsupported(format!("boundary condition '{k}'")));
                }
            }
        }
        // InitStuff.
        let c = ctx.constants();
        let mu0inv = c.real("permeability of vacuum").map_or(1.0 / (std::f64::consts::PI * 4.0e-7), |m| 1.0 / m);
        let eps0 = c.real("permittivity of vacuum").unwrap_or(8.854187817e-12);
        let nbulk = ctx.mesh.n_bulk();
        let active: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nbulk).collect();
        let first = *active.first().ok_or_else(|| SolverError::fatal(CALLER, "no active elements"))?;
        let omega = crate::calcfields2d::angular_frequency(ctx, &sp, first)?.ok_or_else(|| SolverError::fatal(CALLER, "Angular frequency required"))?;
        let mode = sp.integer("eigenfunction index").unwrap_or(1).max(1) as usize - 1;

        // `Variable Output = False` (EMPortSolver_Init0).
        if sp.logical("variable output") == Some(false) {
            let ix = ctx.sim.slots.iter().position(|s| s.id == ctx.solver_id());
            if let Some(id) = ix.and_then(|i| ctx.sim.slots[i].variable) {
                ctx.sim.vars.by_id_mut(id).output = false;
            }
        }
        ctx.default_initialize()?;
        let (mut maxmu, mut maxeps) = (0.0f64, 0.0f64);
        let mut slots = Vec::new();
        for &e in &active {
            let el = ctx.mesh.element(e);
            if el.etype.0 != 303 {
                return Err(unsupported(format!("element type {}", el.etype.0)));
            }
            ctx.element_slots(e, &mut slots);
            let nd = slots.len();
            let (stiff, mass) = Self::local(ctx, e, nd, mu0inv, eps0, omega, &mut maxmu, &mut maxeps)?;
            let to_real = |a: &[C]| {
                let m = 2 * nd;
                let mut r = vec![0.0; m * m];
                for p in 0..nd {
                    for q in 0..nd {
                        let z = a[p * nd + q];
                        r[2 * p * m + 2 * q] = z.re;
                        r[2 * p * m + 2 * q + 1] = -z.im;
                        r[(2 * p + 1) * m + 2 * q] = z.im;
                        r[(2 * p + 1) * m + 2 * q + 1] = z.re;
                    }
                }
                r
            };
            let ls = LocalSystem { n: 2 * nd, stiff: to_real(&stiff), mass: None, damp: None, force: vec![0.0; 2 * nd] };
            ctx.default_update(&ls, e)?;
            ctx.default_update_mass(&to_real(&mass), e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        let fixed = self.dirichlet(ctx)?;
        ctx.default_dirichlet_after(&fixed)?;
        let ix = ctx.sim.slots.iter().position(|s| s.id == ctx.solver_id()).ok_or_else(|| SolverError::fatal(CALLER, "solver slot"))?;
        if sp.logical("eigen system shift automatic") == Some(true) {
            let betalim = omega * (maxeps * maxmu).sqrt();
            ctx.sim.slots[ix].params.insert("Eigen System Shift", Value::real(-(betalim * betalim)));
        }
        ctx.default_solve()?;
        let id = ctx.sim.slots[ix].variable.ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
        let lambda = *ctx.sim.vars.by_id(id).eigen_values.get(mode).ok_or_else(|| SolverError::fatal(CALLER, "no eigenvalue"))?;
        let beta = csqrt(-lambda);
        elmer_solver::results::add_result(&mut ctx.sim.model.simulation, "res: Port Beta 1", beta.re);
        let beta_sum = beta.re;

        if sp.logical("eigenvector normalization by power") == Some(true) {
            let (mut power, mut e2) = (C::new(0.0, 0.0), C::new(0.0, 0.0));
            for &e in &active {
                Self::port_power(ctx, e, mode, beta, mu0inv, omega, &mut power, &mut e2)?;
            }
            let d = C::new((power.re / 2.0).sqrt(), 0.0);
            let v = &mut ctx.sim.vars.by_id_mut(id).eigen_vectors[mode];
            for x in v.iter_mut() {
                *x = div(*x, d);
            }
        }

        let nonhom = sp.integer("number of nonhomogeneous modes").unwrap_or(0).max(0) as usize;
        if nonhom > 0 {
            ctx.sim.slots[ix].params.insert("Eigen Analysis", Value::Logical(false));
            for j in 0..nonhom {
                let lj = ctx.sim.vars.by_id(id).eigen_values[j];
                {
                    let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no system"))?;
                    let m = sys.mass.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no mass matrix"))?.values.clone();
                    for (v, mv) in sys.matrix.values.iter_mut().zip(&m) {
                        *v -= lj.re * mv;
                    }
                }
                if let Some(v) = ctx.primary_mut() {
                    v.values.fill(0.0);
                }
                ctx.default_solve()?;
                let vals = ctx.primary().map(|v| v.values.to_vec()).unwrap_or_default();
                let m = vals.len() / 2;
                let cv: Vec<C> = (0..m).map(|k| C::new(vals[2 * k], vals[2 * k + 1])).collect();
                let var = ctx.sim.vars.by_id_mut(id);
                let ev = &mut var.eigen_vectors[j];
                let (mut udotv, mut udotu) = (C::new(0.0, 0.0), C::new(0.0, 0.0));
                for k in 0..m {
                    udotv += mul(ev[k].conj(), cv[k]);
                    udotu += mul(ev[k].conj(), ev[k]);
                }
                let r = div(udotv, udotu);
                for k in 0..m {
                    ev[k] = cv[k] - mul(r, ev[k]);
                }
                let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no system"))?;
                let mm = sys.mass.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no mass matrix"))?.values.clone();
                for (v, mv) in sys.matrix.values.iter_mut().zip(&mm) {
                    *v += lj.re * mv;
                }
            }
            ctx.sim.slots[ix].params.insert("Eigen Analysis", Value::Logical(true));
        }

        if sp.logical("calculate impedance") == Some(true) {
            let (mut power, mut e2) = (C::new(0.0, 0.0), C::new(0.0, 0.0));
            for &e in &active {
                Self::port_power(ctx, e, mode, beta, mu0inv, omega, &mut power, &mut e2)?;
            }
            let z = div(C::new(1.0, 0.0), power.conj());
            let s = &mut ctx.sim.model.simulation;
            elmer_solver::results::add_result(s, "res: Port Power 1", 0.5 * power.re);
            elmer_solver::results::add_result(s, "res: Port Power Im 1", 0.5 * power.im);
            elmer_solver::results::add_result(s, "res: Port Impedance 1", z.re);
            elmer_solver::results::add_result(s, "res: Port Impedance Im 1", z.im);
        }
        if sp.logical("calculate nodal field") == Some(true) {
            self.post(ctx, &sp, &active, mode, beta)?;
        }
        ctx.sim.vars.by_id_mut(id).conv.norm = beta_sum;
        Ok(())
    }
}

impl EmPort {
    /// `DefaultDirichletBCs`: the `{n}` construct (`Eport re|im {n} m` at
    /// nodal DOF `m` of the boundary nodes) and the `{e}` values (zero) of
    /// the boundary edges, for both components.
    fn dirichlet(&self, ctx: &SolverCtx<'_>) -> Result<Vec<(u32, f64)>, SolverError> {
        let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut fixed = Vec::new();
        let mut bslots = Vec::new();
        let mut pslots = Vec::new();
        let nb = ctx.mesh.n_bulk();
        // {e}: the parent's edge DOF (VectorElementEdgeDOFs of a zero load).
        for (dof, comp) in ["eport re", "eport im"].iter().enumerate() {
            let key = format!("{comp} {{e}}");
            for &b in &boundary {
                let Some(bc) = ctx.bc_params(b) else { continue };
                if !has_prefix(bc, &key) {
                    continue;
                }
                for (k, _) in bc.0.iter().filter(|(k, _)| k.starts_with(&key)) {
                    if crate::material::const_real(bc, k).map_err(SolverError::from)? != Some(0.0) {
                        return Err(unsupported(format!("a nonzero '{k}'")));
                    }
                }
                let el = ctx.mesh.element(b);
                if el.etype.0 != 202 {
                    continue;
                }
                let Some(parent) = el.parents.iter().flatten().copied().find(|p| (p.0 as usize) < nb) else { continue };
                let pel = ctx.mesh.element(parent);
                let emap = crate::topology::edge_map(pel.etype.0 / 100);
                let Some(l) = emap.iter().position(|&[i, j]| el.nodes.contains(&pel.nodes[i]) && el.nodes.contains(&pel.nodes[j])) else { continue };
                ctx.element_slots(parent, &mut pslots);
                let np = 2 * pel.nodes.len();
                let s = pslots[np + l];
                if s != NO_DOF {
                    fixed.push((2 * s + dof as u32, 0.0));
                }
            }
        }
        // {n}.
        for (dof, comp) in ["eport re", "eport im"].iter().enumerate() {
            for m in 0..2 {
                let key = format!("{comp} {{n}} {}", m + 1);
                for &b in &boundary {
                    let Some(bc) = ctx.bc_params(b) else { continue };
                    if bc.get(&key).is_none() {
                        continue;
                    }
                    let el = ctx.mesh.element(b);
                    let n = el.nodes.len();
                    let work = ctx.real_at_nodes(bc, &key, b)?.unwrap_or_else(|| vec![0.0; n]);
                    ctx.element_slots(b, &mut bslots);
                    for j in 0..n {
                        let k = j * 2 + m;
                        let Some(&s) = bslots.get(k) else { continue };
                        if s == NO_DOF {
                            continue;
                        }
                        fixed.push((2 * s + dof as u32, work[j]));
                    }
                }
            }
        }
        Ok(fixed)
    }

    /// `EMPortPost`: the nodal `EF2D` (`Re 1..3, Im 1..3, Re V, Im V`) by
    /// the post solver's mass projection, its perm `NodalPerm` numbering the
    /// nodes at `2*node-1` (`MaxNDOFs = 2`) in the order of the active
    /// elements, its matrix without bandwidth optimization, each component
    /// solved with the `post:` keywords.
    fn post(&self, ctx: &mut SolverCtx<'_>, sp: &ParamMap, active: &[ElemId], mode: usize, beta: C) -> Result<(), SolverError> {
        const DOFS: usize = 8;
        let nnodes = ctx.mesh.n_nodes();
        let perm_len = ctx.primary().and_then(|v| v.perm.map(|p| p.len())).unwrap_or(2 * nnodes);
        // InitialPermutation over the active elements.
        let mut nodal_perm = vec![0u32; perm_len];
        let mut k = 0u32;
        for &e in active {
            for &nd in ctx.mesh.element(e).nodes {
                let idx = 2 * nd as usize;
                if nodal_perm[idx] == 0 {
                    k += 1;
                    nodal_perm[idx] = k;
                }
            }
        }
        let n = k as usize;
        // MakeListMatrix pattern.
        let mut rows: Vec<std::collections::BTreeSet<u32>> = vec![Default::default(); n];
        let pidx = |nd: u32| nodal_perm[2 * nd as usize] - 1;
        for &e in active {
            let nodes = ctx.mesh.element(e).nodes.to_vec();
            for &a in &nodes {
                for &b in &nodes {
                    rows[pidx(a) as usize].insert(pidx(b));
                }
            }
        }
        let mut mat = elmer_core::CsrMatrix { nrows: n, ncols: n, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
        let mut diag = Vec::with_capacity(n);
        for (i, r) in rows.iter().enumerate() {
            for &c in r {
                if c as usize == i {
                    diag.push(Some(mat.col_idx.len()));
                }
                mat.col_idx.push(c);
                mat.values.push(0.0);
            }
            mat.row_ptr.push(mat.col_idx.len());
        }
        mat.diag = diag.into_iter().collect();
        let mut gforce = vec![vec![0.0; n]; DOFS];
        let im = C::new(0.0, 1.0);
        for &e in active {
            let el = ctx.mesh.element(e);
            let nodes = el.nodes.to_vec();
            let nn = nodes.len();
            let (ndofs, np) = (2, 2 * nn);
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let mut normal = tri_normal(&xyz);
            let mut mi = 0;
            for c in 1..3 {
                if normal[c].abs() > normal[mi].abs() {
                    mi = c;
                }
            }
            if normal[mi] < 0.0 {
                normal = normal.map(|x| -x);
            }
            let f = local_mode(ctx, e, mode)?;
            let vdofs = f.len() - np;
            let mut mass = vec![0.0; nn * nn];
            let mut lf = vec![[0.0f64; DOFS]; nn];
            for pt in points(ctx, e)? {
                let s = pt.weight;
                let (mut rez, mut imz) = (0.0, 0.0);
                for q in 0..nn {
                    rez += f[q * ndofs].re * pt.basis[q];
                    imz += f[q * ndofs].im * pt.basis[q];
                }
                let ez = mul(mul(C::new(rez, imz), im), beta);
                let (rez, imz) = (ez.re, ez.im);
                let (mut re_e, mut im_e) = ([0.0f64; 3], [0.0f64; 3]);
                for p in 0..vdofs {
                    for c in 0..3 {
                        re_e[c] += f[np + p].re * pt.w[p][c];
                        im_e[c] += f[np + p].im * pt.w[p][c];
                    }
                }
                for c in 0..3 {
                    let (mut sr, mut si) = (0.0, 0.0);
                    for q in 0..nn {
                        sr += f[q * ndofs + 1].re * pt.db[q][c];
                        si += f[q * ndofs + 1].im * pt.db[q][c];
                    }
                    re_e[c] -= sr;
                    im_e[c] -= si;
                }
                let (mut rev, mut imv) = (0.0, 0.0);
                for q in 0..nn {
                    rev += f[q * ndofs + 1].re * pt.basis[q];
                    imv += f[q * ndofs + 1].im * pt.basis[q];
                }
                let rel: [f64; 3] = std::array::from_fn(|c| re_e[c] + normal[c] * rez);
                let iml: [f64; 3] = std::array::from_fn(|c| im_e[c] + normal[c] * imz);
                for p in 0..nn {
                    for q in 0..nn {
                        mass[p * nn + q] += s * pt.basis[p] * pt.basis[q];
                    }
                    for c in 0..3 {
                        lf[p][c] += s * rel[c] * pt.basis[p];
                        lf[p][3 + c] += s * iml[c] * pt.basis[p];
                    }
                    lf[p][6] += s * rev * pt.basis[p];
                    lf[p][7] += s * imv * pt.basis[p];
                }
            }
            let idx: Vec<usize> = nodes.iter().map(|&nd| pidx(nd) as usize).collect();
            for p in 0..nn {
                let r = idx[p];
                for q in 0..nn {
                    let c = idx[q] as u32;
                    let j = (mat.row_ptr[r]..mat.row_ptr[r + 1]).find(|&j| mat.col_idx[j] == c).expect("pattern");
                    mat.values[j] += mass[p * nn + q];
                }
                for d in 0..DOFS {
                    gforce[d][r] += lf[p][d];
                }
            }
        }
        // The post solver's keywords.
        let mut post = ParamMap::new();
        for (k, v) in &sp.0 {
            if let Some(rest) = k.strip_prefix("post:") {
                post.insert(rest.trim(), v.clone());
            }
        }
        let ef_perm: Vec<u32> = nodal_perm.iter().map(|&p| if p == 0 { NO_DOF } else { p - 1 }).collect();
        let mut ef_vals = vec![0.0; DOFS * n];
        for (d, mut b) in gforce.into_iter().enumerate() {
            let mut x = vec![0.0; n];
            let mut a = mat.clone();
            elmer_solver::solver::solve_other_matrix(&ctx.sim.backends.linear, &ctx.sim.model.simulation, &post, &mut a, &mut b, &mut x, None, ctx.log, ctx.comm)?;
            for i in 0..n {
                ef_vals[DOFS * i + d] = x[i];
            }
        }
        match ctx.sim.vars.id("ef2d") {
            Some(id) => {
                let v = ctx.sim.vars.by_id_mut(id);
                v.values = ef_vals;
                v.perm = Some(std::sync::Arc::from(ef_perm));
                v.values_changed = true;
            }
            None => {
                // VariableAddVector('EF2D[EF2D Re:3 EF2D Im:3 Re V:1 Im V:1]', 8,
                // Perm = NodalPerm, Secondary, Output).
                let spec = elmer_solver::sim::parse_variable_spec("EF2D[EF2D Re:3 EF2D Im:3 Re V:1 Im V:1]", Some(DOFS as i32))?;
                let mut v = elmer_solver::Variable::new(&spec.name, elmer_core::FieldKind::Nodal, DOFS, Some(std::sync::Arc::from(ef_perm)), ef_vals);
                v.output = true;
                v.values_changed = true;
                ctx.sim.vars.add_with_components(v, &spec.components)?;
            }
        }
        Ok(())
    }
}

/// `NormalVector` of a linear triangle (no direction check).
fn tri_normal(x: &[[f64; 3]]) -> [f64; 3] {
    let du = [-1.0, 1.0, 0.0];
    let dv = [-1.0, 0.0, 1.0];
    let dot = |d: &[f64; 3], a: usize| (0..3).fold(0.0, |s, i| s + d[i] * x[i][a]);
    let (xu, yu, zu) = (dot(&du, 0), dot(&du, 1), dot(&du, 2));
    let (xv, yv, zv) = (dot(&dv, 0), dot(&dv, 1), dot(&dv, 2));
    let auu = xu * xu + yu * yu + zu * zu;
    let auv = xu * xv + yu * yv + zu * zv;
    let avv = xv * xv + yv * yv + zv * zv;
    let det = 1.0 / (auu * avv - auv * auv).sqrt();
    [(yu * zv - yv * zu) * det, (xv * zu - xu * zv) * det, (xu * yv - xv * yu) * det]
}
