//! `PoissonBEM.F90` (fem/src/modules): the 2D Poisson (Laplace) equation
//! solved by the boundary element method, `Procedure = "PoissonBEM"
//! "PoissonBEMSolver"`; the real-valued sibling of
//! [`crate::helmholtz_bem`].
//!
//! The active boundary elements (those a boundary condition puts into a
//! body with `Body Id`) carry the unknowns: per node the potential
//! (`Potential` or `<Variable>` given) or the `Flux`. `IntegrateMatrix`
//! collocates at the boundary nodes with `W = -LOG(R)/(2 PI)` along the
//! checked `NormalVector` and the 4-point line rule, `ComputePotential`
//! then evaluates the potential at the other nodes, and the exported
//! `Flux` gets the boundary fluxes. The norm is
//! `SQRT(SUM(Potential**2) / COUNT(Perm > 0))`.
//!
//! `SolveFull`: `direct` is Gaussian elimination with partial pivoting
//! (`SolveLapack`), `iterative` is `FullIterSolver`'s HUTI CGS, from the
//! `RANDOM_NUMBER` start when the variable is still zero.
//!
//! Refused: 3D, body forces with `Source` (`IntegrateSource`), boundary
//! elements other than lines.

#![allow(clippy::needless_range_loop)]

use elmer_core::{ElemId, ParamMap, Params, NO_DOF};
use elmer_fe::metric::MetricStatus;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::helmholtz_bem::{points, Line};

const CALLER: &str = "PoissonBEMSolver";

pub const PROC: (&str, &str) = ("PoissonBEM", "PoissonBEMSolver");

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROC.0, PROC.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(PoissonBem))
}

struct PoissonBem;

fn unsupported(m: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {}", m.into()))
}

const PI: f64 = std::f64::consts::PI;

/// The 2D kernel at offset `(x, y, z)`: `W = -LOG(R)/(2*PI)`,
/// `GradW = -(x, y, z) / (2*PI*R**2)`.
fn green(x: f64, y: f64, z: f64) -> (f64, [f64; 3]) {
    let r = ((x * x + y * y) + z * z).sqrt();
    let w = -r.ln() / (2.0 * PI);
    let d = 2.0 * PI * (r * r);
    (w, [-x / d, -y / d, -z / d])
}

fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let mut s = 0.0;
    for k in 0..3 {
        s += a[k] * b[k];
    }
    s
}

/// Gaussian elimination with partial pivoting on the column-major `a`.
fn dense_solve(n: usize, a: &mut [f64], b: &mut [f64]) -> Result<(), SolverError> {
    for k in 0..n {
        let mut piv = k;
        for i in k + 1..n {
            if a[i + k * n].abs() > a[piv + k * n].abs() {
                piv = i;
            }
        }
        if a[piv + k * n] == 0.0 {
            return Err(SolverError::fatal("SolveFull", "singular matrix"));
        }
        if piv != k {
            for j in 0..n {
                a.swap(k + j * n, piv + j * n);
            }
            b.swap(k, piv);
        }
        for i in k + 1..n {
            let m = a[i + k * n] / a[k + k * n];
            if m == 0.0 {
                continue;
            }
            for j in k + 1..n {
                a[i + j * n] -= m * a[k + j * n];
            }
            b[i] -= m * b[k];
        }
    }
    for k in (0..n).rev() {
        let mut s = b[k];
        for j in k + 1..n {
            s -= a[k + j * n] * b[j];
        }
        b[k] = s / a[k + k * n];
    }
    Ok(())
}

/// Reference BLAS `DDOT` (unit stride): the `MOD(n,5)` leading terms,
/// then five products per step.
fn ddot(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    let m = n % 5;
    let mut t = 0.0;
    for i in 0..m {
        t += x[i] * y[i];
    }
    let mut i = m;
    while i < n {
        t = t + x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2] + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
        i += 5;
    }
    t
}

fn nrm2(x: &[f64]) -> f64 {
    let mut s = 0.0;
    for &v in x {
        s += v * v;
    }
    s.sqrt()
}

/// The module's `Matvec` on the column-major dense matrix.
fn matvec(n: usize, a: &[f64], u: &[f64], v: &mut [f64]) {
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            s += a[i + j * n] * u[j];
        }
        v[i] = s;
    }
}

/// `FullIterSolver`: `huti_dcgssolv` with the module's Jacobi `Precond`
/// on the left, the identity on the right, the true residual
/// `DNRM2(A x - b)` against `Linear System Convergence Tolerance`, and a
/// random start (`huti_drandvec`, `REAL(4)` harvests) when `x` is zero.
fn full_iter_solver(ctx: &SolverCtx<'_>, sp: &ParamMap, n: usize, a: &[f64], b: &[f64], mut x: Vec<f64>) -> Result<Vec<f64>, SolverError> {
    let tol = sp.real("linear system convergence tolerance").unwrap_or(0.0);
    let maxtol = sp.real("linear system divergence limit").unwrap_or(1.0e20);
    let maxit = sp.integer("linear system max iterations").unwrap_or(0).max(0) as usize;
    let abort = sp.logical("linear system abort not converged").unwrap_or(true);
    if x.iter().all(|&v| v == 0.0) {
        for v in x.iter_mut() {
            *v = f64::from(elmer_itermethods::even_random_r4());
        }
    }
    let diag: Vec<f64> = (0..n).map(|i| a[i + i * n]).collect();
    let precond = |u: &mut [f64], v: &[f64]| {
        for i in 0..n {
            u[i] = v[i] / diag[i];
        }
    };
    let mut r = vec![0.0; n];
    matvec(n, a, &x, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
    }
    let rtld = r.clone();
    let (mut u, mut p, mut q) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let (mut t1v, mut t2v) = (vec![0.0; n], vec![0.0; n]);
    let mut oldrho = 0.0;
    let mut iter = 1usize;
    let converged = loop {
        let rho = ddot(&rtld, &r);
        if rho == 0.0 {
            break false;
        }
        if iter == 1 {
            u.copy_from_slice(&r);
            p.copy_from_slice(&u);
        } else {
            let beta = rho / oldrho;
            for i in 0..n {
                u[i] = r[i] + beta * q[i];
                p[i] = u[i] + beta * q[i] + beta * beta * p[i];
            }
        }
        precond(&mut t2v, &p);
        t1v.copy_from_slice(&t2v);
        matvec(n, a, &t1v, &mut t2v);
        let alpha = rho / ddot(&rtld, &t2v);
        for i in 0..n {
            q[i] = u[i] - alpha * t2v[i];
            t2v[i] = u[i] + q[i];
        }
        precond(&mut u, &t2v);
        t1v.copy_from_slice(&u);
        for i in 0..n {
            x[i] += alpha * t1v[i];
        }
        matvec(n, a, &t1v, &mut t2v);
        for i in 0..n {
            r[i] -= alpha * t2v[i];
        }
        matvec(n, a, &x, &mut t1v);
        for i in 0..n {
            t1v[i] -= b[i];
        }
        let residual = nrm2(&t1v);
        if residual < tol {
            break true;
        }
        if residual.is_nan() || residual > maxtol {
            break false;
        }
        oldrho = rho;
        iter += 1;
        if iter > maxit {
            break false;
        }
    };
    if !converged {
        if abort {
            return Err(SolverError::fatal("IterSolve", "Failed convergence tolerances."));
        }
        ctx.info("IterSolve", "Failed convergence tolerances.", 0);
    }
    Ok(x)
}

impl Solver for PoissonBem {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.comm.size() > 1 {
            return Err(unsupported("parallel run"));
        }
        let sp = ctx.params().clone();
        let cs = ctx.simulation_params().string("coordinate system").unwrap_or("cartesian").to_ascii_lowercase();
        if !cs.trim().starts_with("cartesian") {
            return Err(SolverError::fatal(CALLER, "This solver is implemented only for cartesian coordinates!"));
        }
        if ctx.mesh.dim() != 2 {
            return Err(unsupported("the 3D kernel (surface boundary elements)"));
        }
        match sp.string("linear system solver").map(|s| s.trim().to_ascii_lowercase()).as_deref() {
            Some("direct") | Some("iterative") => {}
            _ => return Err(SolverError::fatal("SolveFull", "Unknown solver type.")),
        }
        if ctx.sim.model.body_forces.iter().any(|bf| bf.get("source").is_some()) {
            return Err(unsupported("a body force Source (IntegrateSource)"));
        }
        let var_name = ctx.variable_name().unwrap_or("").to_ascii_lowercase();
        let (vdofs, vperm, vlen) = {
            let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
            (v.dofs, v.perm.map(|p| p.to_vec()), v.values.len())
        };
        if vdofs != 1 {
            return Err(SolverError::fatal(CALLER, "the variable needs one component"));
        }
        let nnodes = ctx.mesh.n_nodes();
        let nbulk = ctx.mesh.n_bulk();
        let active: Vec<ElemId> = ctx.active_elements().collect();

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

        // Known potentials (`Potential`, then the variable) or fluxes.
        let mut load = vec![0.0f64; nb];
        let mut known = vec![false; nb];
        for l in &lines {
            let Some(bc) = ctx.bc_params(l.e) else { continue };
            let n = l.nodes.len();
            let (vals, got) = match ctx.real_at_nodes(bc, "potential", l.e)? {
                Some(v) => (v, true),
                None => match ctx.real_at_nodes(bc, &var_name, l.e)? {
                    Some(v) => (v, true),
                    None => (ctx.real_at_nodes(bc, "flux", l.e)?.unwrap_or_else(|| vec![0.0; n]), false),
                },
            };
            for (j, &node) in l.nodes.iter().enumerate() {
                let q = pperm[node as usize];
                load[q] = vals[j];
                known[q] = got;
            }
        }

        // IntegrateMatrix.
        let mut status = MetricStatus::default();
        let mut mat = vec![0.0f64; nb * nb];
        let mut diag = vec![0.0f64; nb];
        let mut force = vec![0.0f64; nb];
        let coord = |k: usize| ctx.mesh.coord(k as u32);
        for l in &lines {
            for qp in points(ctx, l, &mut status)? {
                let s = qp.s;
                for p in 0..nb {
                    let c = coord(bperm[p]);
                    let (w, grad) = green(qp.p[0] - c[0], qp.p[1] - c[1], qp.p[2] - c[2]);
                    let dwdn = dot3(&grad, &qp.normal);
                    for (i, &node) in l.nodes.iter().enumerate() {
                        let q = pperm[node as usize];
                        let bi = qp.basis[i];
                        if known[q] {
                            if p != q {
                                force[p] -= s * load[q] * bi * dwdn;
                            }
                            mat[p + q * nb] -= s * bi * w;
                        } else {
                            force[p] += s * load[q] * bi * w;
                            if p != q {
                                mat[p + q * nb] += s * bi * dwdn;
                            }
                        }
                        if p != q {
                            diag[p] -= s * bi * dwdn;
                        }
                    }
                }
            }
        }
        for i in 0..nb {
            if known[i] {
                force[i] -= load[i] * diag[i];
            } else {
                mat[i + i * nb] = diag[i];
            }
        }

        // SolveFull: `SolveLapack` or `FullIterSolver` (HUTI CGS) into
        // `Potential(1:BoundaryNodes)`.
        let iterative = sp.string("linear system solver").is_some_and(|s| s.trim().eq_ignore_ascii_case("iterative"));
        if iterative {
            let x0: Vec<f64> = {
                let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
                (0..nb).map(|i| v.values.get(i).copied().unwrap_or(0.0)).collect()
            };
            let x = full_iter_solver(ctx, &sp, nb, &mat, &force, x0)?;
            force = x;
        } else {
            dense_solve(nb, &mut mat, &mut force)?;
        }
        let (mut pot, mut flx) = (vec![0.0f64; nb], vec![0.0f64; nb]);
        for i in 0..nb {
            if known[i] {
                flx[i] = force[i];
                pot[i] = load[i];
            } else {
                flx[i] = load[i];
                pot[i] = force[i];
            }
        }

        // `Potential = VolumeForce` (zero), the boundary values, then
        // ComputePotential at the other nodes in the solver's perm.
        let slot = |i: usize| -> Option<usize> {
            match &vperm {
                Some(p) if p[i] == NO_DOF => None,
                Some(p) => Some(p[i] as usize),
                None => Some(i),
            }
        };
        let mut potential = vec![0.0f64; vlen];
        for i in 0..nb {
            if let Some(k) = slot(bperm[i]) {
                potential[k] = pot[i];
            }
        }
        if sp.logical("skip volume computing") != Some(true) {
            for l in &lines {
                for qp in points(ctx, l, &mut status)? {
                    let s = qp.s;
                    for i in 0..nnodes {
                        let Some(k) = slot(i) else { continue };
                        if pperm[i] != usize::MAX {
                            continue;
                        }
                        let c = coord(i);
                        let (w, grad) = green(qp.p[0] - c[0], qp.p[1] - c[1], qp.p[2] - c[2]);
                        let dwdn = dot3(&grad, &qp.normal);
                        for (j, &node) in l.nodes.iter().enumerate() {
                            let q = pperm[node as usize];
                            potential[k] -= s * qp.basis[j] * (pot[q] * dwdn - flx[q] * w);
                        }
                    }
                }
            }
        }

        // The fields and the norm.
        let mut sum = 0.0;
        for &x in &potential {
            sum += x * x;
        }
        let count = vperm.as_ref().map_or(nnodes, |p| p.iter().filter(|&&s| s != NO_DOF).count());
        let norm = (sum / count as f64).sqrt();
        {
            let v = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
            v.values.copy_from_slice(&potential);
        }
        if let Some(f) = ctx.variable_mut("flux") {
            let fperm = f.perm.clone();
            f.values.fill(0.0);
            for i in 0..nb {
                let k = match fperm.as_deref() {
                    Some(p) if p[bperm[i]] == NO_DOF => continue,
                    Some(p) => p[bperm[i]] as usize,
                    None => bperm[i],
                };
                if k < f.values.len() {
                    f.values[k] = flx[i];
                }
            }
        }
        let id = ctx.sim.slots.iter().find(|sl| sl.id == ctx.solver_id()).and_then(|sl| sl.variable);
        if let Some(id) = id {
            let v = ctx.sim.vars.by_id_mut(id);
            v.conv.norm = norm;
            v.values_changed = true;
        }
        Ok(())
    }
}
