//! `StokesSolver` of the library module `ParStokes` (fem/src/modules/
//! ParStokes.F90) with its block preconditioners `VelocityPrecond` and
//! `PressurePrecond` (VelocityPrecond.F90, PressurePrecond.F90): Stokes
//! flow in one `dim + 1` dof variable, optionally P2-P1 (`P2-P1
//! Approximation`: pressure on the corner nodes, the linear element's
//! basis), solved by a GCR outer iteration whose preconditioner is a
//! pressure solve with the `-1/mu` pressure mass matrix and a velocity
//! solve with the velocity block (`Block Preconditioning`), each on the
//! matrix of its own solver section (`Equation = Velocity Preconditioning`
//! / `Pressure Preconditioning`).
//!
//! The block-preconditioned path follows the Fortran step by step: an
//! extra nonlinear round for the true residual, `Linear System Row
//! Equilibration` of the main system (kept, not back-scaled) with the
//! preconditioner right-hand sides scaled back by the row sums, the rotated
//! solution carried between rounds and `BackRotateNTSystem` after each,
//! and the P2-P1 mid-edge pressures interpolated at the end.
//!
//! The Fortran turns on `No Precondition Recompute` for the preconditioner
//! solves from the second GCR round; their matrices do not change between
//! rounds, so each solve here recomputes the same ILU instead (equal up to
//! the rounding of their scaling round trip).
//!
//! Refused: `Viscosity Model` (power law, anisotropy), `Convective`,
//! `Block Diagonal A`, `Calculate Loads`, transient runs, p elements.

use elmer_core::{BasisEval, ElemId, ElementType, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "ParStokes";

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("ParStokes StokesSolver: {what}"))
}

/// `ListAddNew*`: set `key` only where the section lacks it.
fn add_new(ctx: &mut SolverCtx<'_>, key: &str, v: Value) {
    let id = ctx.solver_id();
    if let Some(s) = ctx.sim.slots.iter_mut().find(|s| s.id == id) {
        if s.params.get(&elmer_core::normalize_key(key)).is_none() {
            s.params.insert(key, v);
        }
    }
}

/// `VelocityPrecond` / `PressurePrecond`: solver sections whose matrices
/// ParStokes assembles and solves; their own solve does nothing.
pub struct Precond {
    pressure: bool,
}

pub fn velocity_precond_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Precond { pressure: false }))
}

pub fn pressure_precond_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Precond { pressure: true }))
}

impl Solver for Precond {
    /// `VelocityPrecond_Init0` / `PressurePrecond_Init0`: BiCGStabL + ILU0
    /// to 1e-4 for the velocity, BiCGStab2 + Diagonal to 1e-6 for the
    /// pressure, as defaults; the variable is not saved (`Variable Output`,
    /// a default for the velocity, forced for the pressure).
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let id = ctx.solver_id();
        if let Some(s) = ctx.sim.slots.iter_mut().find(|s| s.id == id) {
            s.params.insert("Exec Solver", Value::Str("before simulation".into()));
            s.params.insert("Skip Compute Nonlinear Change", Value::Logical(true));
            s.params.insert("Back Rotate N-T Solution", Value::Logical(false));
        }
        let output = !self.pressure && ctx.params().logical("variable output").unwrap_or(false);
        let (method, precond, tol) = if self.pressure { ("BiCGStab2", "Diagonal", 1.0e-6) } else { ("BiCGStabL", "ILU0", 1.0e-4) };
        add_new(ctx, "Bubbles in Global System", Value::Logical(false));
        add_new(ctx, "Linear System Solver", Value::Str("Iterative".into()));
        add_new(ctx, "Linear System Iterative Method", Value::Str(method.into()));
        add_new(ctx, "Linear System Max Iterations", Value::Integer(vec![1000]));
        add_new(ctx, "Linear System Preconditioning", Value::Str(precond.into()));
        add_new(ctx, "Linear System Convergence Tolerance", Value::real(tol));
        add_new(ctx, "Linear System Abort Not Converged", Value::Logical(false));
        if let Some(name) = ctx.variable_name().map(str::to_string) {
            if let Some(v) = ctx.sim.vars.get_mut(&name) {
                v.output = output;
            }
        }
        Ok(())
    }

    fn solve(&mut self, _ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        Ok(())
    }
}

pub struct ParStokes;

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(ParStokes))
}

/// The local systems of one element: `STIFF` (`(dim+1) nd` square,
/// interleaved), `FORCE`, the pressure block `Mass` (`nd` square).
struct Local {
    stiff: Vec<f64>,
    force: Vec<f64>,
    mass: Vec<f64>,
}

/// The P2-P1 pressure element of a family: its linear type.
fn linear_code(family: u16) -> Result<u16, SolverError> {
    Ok(match family {
        3 => 303,
        4 => 404,
        5 => 504,
        6 => 605,
        7 => 706,
        8 => 808,
        f => return Err(unsupported(&format!("P2-P1 on element family {f}"))),
    })
}

/// `LocalMatrix(STIFF, VeloBlock, Mass, FORCE, LOAD, rho, mu, Vx, Vy, Vz,
/// Element, n, nd, dim, ...)` without convection, a viscosity model or a
/// block-diagonal velocity block.
#[allow(clippy::too_many_arguments)]
fn local_matrix(
    ctx: &SolverCtx<'_>,
    e: ElemId,
    load: &[[f64; 3]],
    rho: &[f64],
    mu: &[f64],
    dim: usize,
    p2p1: bool,
) -> Result<Local, SolverError> {
    let el = ctx.mesh.element(e);
    let (etype, n) = (el.etype, el.nodes.len());
    let nd = n;
    let c = dim + 1;
    let m = c * nd;
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no FE backend"))?;
    let (ltype, nlin) = if p2p1 {
        let lt = ElementType(linear_code(etype.family())?);
        (Some(lt), lt.n_nodes())
    } else {
        (None, n)
    };
    let mut out = Local { stiff: vec![0.0; m * m], force: vec![0.0; m], mass: vec![0.0; nd * nd] };
    let mut be = BasisEval::default();
    let mut lb = BasisEval::default();
    for qp in fe.quadrature.gauss_points(etype, None) {
        fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|e| SolverError::fatal("ElementInfo", e))?;
        if let Some(lt) = ltype {
            // NodalBasisFunctions of the linear element at the same point.
            fe.basis.eval(lt, &xyz[..nlin], qp.uvw, &mut lb).map_err(|e| SolverError::fatal("NodalBasisFunctions", e))?;
        }
        let (b, g) = (&be.basis[..nd], &be.dbasis_dx[..nd]);
        let s = qp.weight * be.det_j;
        let visc = (0..n).fold(0.0, |acc, k| acc + b[k] * mu[k]);
        let rho_ip = (0..n).fold(0.0, |acc, k| acc + b[k] * rho[k]);
        let mu_ip = visc;
        // LoadAtIP = MATMUL(Basis(1:n), LOAD(1:n,1:dim+1)).
        let load_ip: [f64; 4] = std::array::from_fn(|j| if j < 3 { (0..n).fold(0.0, |acc, k| acc + b[k] * load[k][j]) } else { 0.0 });
        let at = |p: usize, i: usize, q: usize, j: usize| (c * p + i) * m + c * q + j;
        for p in 0..nd {
            for q in 0..nd {
                for i in 0..dim {
                    for j in 0..dim {
                        out.stiff[at(p, i, q, i)] += s * mu_ip * g[q][j] * g[p][j];
                        out.stiff[at(p, i, q, j)] += s * mu_ip * g[q][i] * g[p][j];
                    }
                }
            }
        }
        for p in 0..nd {
            for q in 0..nd {
                for i in 0..dim {
                    if p2p1 {
                        let lbas = &lb.basis[..nlin];
                        if q < nlin {
                            out.stiff[at(p, i, q, dim)] -= s * lbas[q] * g[p][i];
                        }
                        if p < nlin {
                            out.stiff[at(p, dim, q, i)] -= s * g[q][i] * lbas[p];
                        }
                    } else {
                        if q < n {
                            out.stiff[at(p, i, q, dim)] -= s * b[q] * g[p][i];
                        }
                        if p < n {
                            out.stiff[at(p, dim, q, i)] -= s * g[q][i] * b[p];
                        }
                    }
                }
            }
            for k in 0..c {
                out.force[c * p + k] += s * rho_ip * load_ip[k] * b[p];
            }
        }
        // The pressure block: -1/mu Basis Basis (the full basis, as in the
        // Fortran, on the first nlin functions).
        for p in 0..nlin {
            for q in 0..nlin {
                out.mass[p * nd + q] -= s * 1.0 / mu_ip * b[p] * b[q];
            }
        }
    }
    for p in nlin..nd {
        let i = c * (p + 1) - 1;
        out.force[i] = 0.0;
        for k in 0..m {
            out.stiff[i * m + k] = 0.0;
            out.stiff[k * m + i] = 0.0;
        }
        out.stiff[i * m + i] = 1.0;
        out.mass[p * nd + p] = 1.0;
    }
    Ok(out)
}

/// `NavierStokesBoundary` of ParStokes: the given traction (`Surface
/// Traction i` / `Pressure i`), a normal one (`Normal Surface Traction` /
/// `External Pressure`) and slip, in cartesian or n-t components.
fn boundary_local(ctx: &SolverCtx<'_>, e: ElemId, bc: &ParamMap, dim: usize, nt: bool) -> Result<Vec<f64>, SolverError> {
    let el = ctx.mesh.element(e);
    let (etype, n) = (el.etype, el.nodes.len());
    let c = dim + 1;
    let m = c * n;
    let nodal = |key: &str| -> Result<Option<Vec<f64>>, SolverError> { ctx.real_at_nodes(bc, key, e) };
    let ext = match nodal("normal surface traction")? {
        Some(v) => v,
        None => nodal("external pressure")?.unwrap_or_else(|| vec![0.0; n]),
    };
    let mut load = vec![[0.0; 3]; n];
    let mut slip = vec![[0.0; 3]; n];
    for i in 0..3 {
        let t = match nodal(&format!("surface traction {}", i + 1))? {
            Some(v) => Some(v),
            None => nodal(&format!("pressure {}", i + 1))?,
        };
        if let Some(t) = t {
            for (l, x) in load.iter_mut().zip(t) {
                l[i] = x;
            }
        }
        if let Some(sc) = nodal(&format!("slip coefficient {}", i + 1))? {
            for (l, x) in slip.iter_mut().zip(sc) {
                l[i] = x;
            }
        }
    }
    let any_slip = slip.iter().any(|s| s.iter().any(|&x| x != 0.0));
    let mut out = vec![0.0; m * m + m];
    let (stiff, force) = out.split_at_mut(m * m);
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no FE backend"))?;
    let mut be = BasisEval::default();
    for qp in fe.quadrature.gauss_points(etype, None) {
        fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|e| SolverError::fatal("ElementInfo", e))?;
        let b = &be.basis[..n];
        let s = be.det_j * qp.weight;
        let mut f = [0.0; 3];
        for (i, fi) in f.iter_mut().enumerate().take(dim) {
            *fi += (0..n).fold(0.0, |acc, k| acc + load[k][i] * b[k]);
        }
        let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [qp.uvw[0], qp.uvw[1]])?;
        let alpha = (0..n).fold(0.0, |acc, k| acc + ext[k] * b[k]);
        if nt {
            f[0] += alpha;
        } else {
            for i in 0..dim {
                f[i] += alpha * normal[i];
            }
        }
        if dim != 2 && (nt || any_slip) {
            return Err(unsupported("n-t or slip boundaries in 3D"));
        }
        let vects = [normal, [normal[1], -normal[0], 0.0], [0.0; 3]];
        if any_slip {
            for p in 0..n {
                for q in 0..n {
                    for i in 0..dim {
                        let sc = (0..n).fold(0.0, |acc, k| acc + slip[k][i] * b[k]);
                        if nt {
                            for j in 0..dim {
                                for k in 0..dim {
                                    stiff[(p * c + j) * m + q * c + k] += s * sc * b[q] * b[p] * vects[i][j] * vects[i][k];
                                }
                            }
                        } else {
                            stiff[(p * c + i) * m + q * c + i] += s * sc * b[q] * b[p];
                        }
                    }
                }
            }
        }
        for q in 0..n {
            for i in 0..dim {
                if nt {
                    for j in 0..dim {
                        force[q * c + j] += s * b[q] * f[i] * vects[i][j];
                    }
                } else {
                    force[q * c + i] += s * b[q] * f[i];
                }
            }
        }
    }
    Ok(out)
}

/// The velocity block of a `(dim+1)`-interleaved local matrix.
fn velocity_block(stiff: &[f64], nd: usize, dim: usize) -> LocalSystem {
    let (c, m, mv) = (dim + 1, (dim + 1) * nd, dim * nd);
    let mut a = LocalSystem::zeros(mv);
    for p in 0..nd {
        for i in 0..dim {
            for q in 0..nd {
                for j in 0..dim {
                    a.stiff[(p * dim + i) * mv + q * dim + j] = stiff[(p * c + i) * m + q * c + j];
                }
            }
        }
    }
    a
}

fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).fold(0.0, |s, (a, b)| s + a * b)
}

fn norm(x: &[f64]) -> f64 {
    dot(x, x).sqrt()
}

fn matvec(a: &elmer_core::CsrMatrix<f64>, x: &[f64], y: &mut [f64]) {
    elmer_solver::crs_matrix_vector_multiply(a, 0, x, y);
}

/// The GCR state carried over the rounds.
struct Gcr {
    s: Vec<Vec<f64>>,
    v: Vec<Vec<f64>>,
    rr: Vec<f64>,
}

impl Gcr {
    /// `GCRUpdate(n, A, M, x, b, r, Snew, S, V, R, Round, res, MParam)`.
    #[allow(clippy::too_many_arguments)]
    fn update(
        &mut self,
        a: &elmer_core::CsrMatrix<f64>,
        x: &mut [f64],
        b: &[f64],
        r: &mut [f64],
        snew: &[f64],
        round: usize,
        mparam: usize,
        ctx: &SolverCtx<'_>,
    ) -> f64 {
        let n = x.len();
        let j = if round.is_multiple_of(mparam) { mparam } else { round % mparam };
        if j == 1 {
            matvec(a, x, r);
            for i in 0..n {
                r[i] = b[i] - r[i];
            }
            let res = norm(r) / norm(b);
            if round == 1 {
                ctx.info("Outer Iteration", &format!("GCR residual for iterate{:>4}{}", 0, es12_3(res)), 4);
            }
        } else {
            r.copy_from_slice(&self.rr);
        }
        let mut t1 = snew.to_vec();
        let mut t2 = vec![0.0; n];
        matvec(a, &t1, &mut t2);
        for i in 0..j - 1 {
            let beta = dot(&self.v[i], &t2);
            for k in 0..n {
                t1[k] -= beta * self.s[i][k];
                t2[k] -= beta * self.v[i][k];
            }
        }
        let alpha = norm(&t2);
        for k in 0..n {
            t1[k] *= 1.0 / alpha;
            t2[k] *= 1.0 / alpha;
        }
        let beta = dot(&t2, r);
        for k in 0..n {
            x[k] += beta * t1[k];
            r[k] -= beta * t2[k];
        }
        if j != mparam {
            if self.s.len() < j {
                self.s.resize(j, vec![0.0; n]);
                self.v.resize(j, vec![0.0; n]);
            }
            self.s[j - 1] = t1;
            self.v[j - 1] = t2;
        }
        self.rr = r.to_vec();
        let res = norm(r) / norm(b);
        ctx.info("Outer Iteration", &format!("GCR residual for iterate{:>4}{}", round, es12_3(res)), 4);
        res
    }
}

/// Fortran `ES12.3`.
fn es12_3(x: f64) -> String {
    let s = format!("{:.3E}", x);
    let (m, e) = s.split_once('E').expect("exponent");
    let e: i32 = e.parse().expect("integer");
    format!("{:>12}", format!("{m}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs()))
}

/// The slot index of the solver whose `Equation` is `eq`.
fn slot_by_equation(ctx: &SolverCtx<'_>, eq: &str) -> Option<usize> {
    ctx.sim.slots.iter().position(|s| s.params.string("equation").map(elmer_core::normalize_key).as_deref() == Some(eq))
}

impl Solver for ParStokes {
    /// `StokesSolver_Init0`.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        add_new(ctx, "Bubbles in Global System", Value::Logical(false));
        add_new(ctx, "Linear System Solver", Value::Str("Iterative".into()));
        add_new(ctx, "Linear System Iterative Method", Value::Str("GCR".into()));
        add_new(ctx, "Linear System GCR Restart", Value::Integer(vec![50]));
        add_new(ctx, "Linear System Max Iterations", Value::Integer(vec![200]));
        add_new(ctx, "Linear System Row Equilibration", Value::Logical(true));
        add_new(ctx, "Linear System Convergence Tolerance", Value::real(1.0e-6));
        add_new(ctx, "Linear System Base Tolerance", Value::real(1.0e-3));
        add_new(ctx, "Linear System Relative Tolerance", Value::real(1.0e-2));
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("transient runs"));
        }
        if ctx.solver_dofs().is_some() {
            return Err(unsupported("p elements (bubbles)"));
        }
        let p = ctx.params().clone();
        let dim = ctx.mesh.dim();
        let c = dim + 1;
        let p2p1 = p.logical("p2-p1 approximation") == Some(true);
        let block = p.logical("block preconditioning") == Some(true);
        for (key, what) in [
            ("convective", "Convective"),
            ("block diagonal a", "Block Diagonal A"),
            ("calculate loads", "Calculate Loads"),
            ("use velocity laplacian", "Use Velocity Laplacian"),
        ] {
            if p.logical(key) == Some(true) {
                return Err(unsupported(what));
            }
        }
        if ctx.sim.model.materials.iter().any(|m| m.get("viscosity model").is_some()) {
            return Err(unsupported("Viscosity Model"));
        }
        let scaling = p.logical("linear system scaling").unwrap_or(false);
        let equil = p.logical("linear system row equilibration").unwrap_or(false);
        if block && scaling && !equil {
            return Err(unsupported("Linear System Scaling without Row Equilibration in the block path"));
        }
        let do_scaling = scaling || equil;
        let adaptive = p.logical("linear system adaptive tolerance") == Some(true);
        let name = ctx.variable_name().unwrap_or_default().to_string();
        let dofs = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?.dofs;
        if dofs != c {
            return Err(unsupported("a variable without dim + 1 dofs"));
        }

        let (pix, vix) = if block {
            let pix = slot_by_equation(ctx, "pressure preconditioning")
                .ok_or_else(|| SolverError::fatal(CALLER, "Pressure preconditioning operation missing..."))?;
            let vix = slot_by_equation(ctx, "velocity preconditioning")
                .ok_or_else(|| SolverError::fatal(CALLER, "Velocity preconditioning operation missing..."))?;
            let my_perm = ctx.primary().and_then(|v| v.perm.map(|p| p.to_vec()));
            for (ix, what) in [(pix, "pressure"), (vix, "velocity")] {
                let other = ctx.sim.slots[ix].variable.and_then(|id| ctx.sim.vars.by_id(id).perm.clone()).map(|p| p.to_vec());
                if other != my_perm {
                    return Err(SolverError::fatal(
                        CALLER,
                        format!("Nonmatching variable permutations, use Optimize Bandwidth in {what} preconditioning"),
                    ));
                }
            }
            (Some(pix), Some(vix))
        } else {
            (None, None)
        };

        let mut nonlinear_iter = p.integer("nonlinear system max iterations").unwrap_or(0).max(0) as usize;
        let nonlinear_tol = p.real("nonlinear system convergence tolerance").unwrap_or(0.0);
        if block {
            nonlinear_iter += 1;
        }
        let max_iterations = p.integer("linear system max iterations").unwrap_or(200).max(1) as usize;
        let restart = p.integer("linear system gcr restart").map_or(max_iterations, |r| r.max(1) as usize);
        let target_tol = p.real("linear system convergence tolerance").unwrap_or(1.0e-6);
        let rel_tol = p.real("linear system relative tolerance").unwrap_or(1.0e-2);
        let base_tol = p.real("linear system base tolerance").unwrap_or(1.0e-4);
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let (mesh, log, comm) = (ctx.mesh, ctx.log, ctx.comm);
        let mut rotated: Vec<f64> = Vec::new();
        let mut non_lin_error = 2.0;

        for iter in 1..=nonlinear_iter {
            ctx.default_initialize()?;
            for ix in [vix, pix].into_iter().flatten() {
                SolverCtx::new(&mut *ctx.sim, ix, mesh, log, comm).default_initialize()?;
            }
            if block && iter == 1 {
                rotated = vec![0.0; ctx.system().map_or(0, |s| s.rhs.len())];
            }
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let nodes = el.nodes.to_vec();
                let n = nodes.len();
                let bp = ctx.body_params(e).ok_or_else(|| SolverError::fatal(CALLER, "an element without a body"))?;
                let (mat, bf) = (bp.material.cloned(), bp.body_force.cloned());
                let mut load = vec![[0.0; 3]; n];
                if let Some(b) = &bf {
                    for i in 0..dim {
                        if let Some(v) = ctx.real_at_nodes(b, &format!("flow bodyforce {}", i + 1), e)? {
                            for (l, x) in load.iter_mut().zip(v) {
                                l[i] = x;
                            }
                        }
                    }
                }
                let get = |key: &str| -> Result<Vec<f64>, SolverError> {
                    match &mat {
                        Some(m) => Ok(ctx.real_at_nodes(m, key, e)?.unwrap_or_else(|| vec![0.0; n])),
                        None => Ok(vec![0.0; n]),
                    }
                };
                let (rho, mu) = (get("density")?, get("viscosity")?);
                let loc = local_matrix(ctx, e, &load, &rho, &mu, dim, p2p1)?;
                let mut main = LocalSystem::zeros(c * n);
                main.stiff.copy_from_slice(&loc.stiff);
                main.force.copy_from_slice(&loc.force);
                ctx.default_update(&main, e)?;
                if let (Some(vix), Some(pix)) = (vix, pix) {
                    let a = velocity_block(&loc.stiff, n, dim);
                    SolverCtx::new(&mut *ctx.sim, vix, mesh, log, comm).default_update(&a, e)?;
                    let mut pl = LocalSystem::zeros(n);
                    pl.stiff.copy_from_slice(&loc.mass);
                    SolverCtx::new(&mut *ctx.sim, pix, mesh, log, comm).default_update(&pl, e)?;
                }
            }
            ctx.default_finish_bulk_assembly()?;

            let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
            let mut slots = Vec::new();
            for e in bnd {
                ctx.element_dofs(e, &mut slots);
                if slots.is_empty() || slots.contains(&NO_DOF) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                if el.etype.family() == 1 {
                    continue;
                }
                let n = el.nodes.len();
                let Some(bc) = ctx.bc_params(e).cloned() else { continue };
                if bc.logical("flow force bc") == Some(false) {
                    continue;
                }
                let nt = bc.logical(&format!("normal-tangential {name}")) == Some(true);
                let both = boundary_local(ctx, e, &bc, dim, nt)?;
                let m = c * n;
                let mut main = LocalSystem::zeros(m);
                main.stiff.copy_from_slice(&both[..m * m]);
                main.force.copy_from_slice(&both[m * m..]);
                ctx.default_update(&main, e)?;
                if let Some(vix) = vix {
                    let a = velocity_block(&both[..m * m], n, dim);
                    SolverCtx::new(&mut *ctx.sim, vix, mesh, log, comm).default_update(&a, e)?;
                }
            }
            ctx.default_dirichlet()?;
            for ix in [pix, vix].into_iter().flatten() {
                SolverCtx::new(&mut *ctx.sim, ix, mesh, log, comm).default_dirichlet()?;
            }

            let (Some(pix), Some(vix)) = (pix, vix) else {
                ctx.default_solve()?;
                if ctx.convergence().is_some_and(|cv| cv.nonlin_converged == 1) {
                    break;
                }
                continue;
            };

            // Row equilibration of the main system; the preconditioners'
            // right-hand sides are scaled back by the row sums.
            let nrows = ctx.system().map_or(0, |s| s.rhs.len());
            let np = nrows / c;
            let (mut row_p, mut row_a) = (vec![1.0; np], vec![1.0; np * dim]);
            if do_scaling {
                let d = ctx.row_equilibrate()?;
                for i in 0..np {
                    row_p[i] = 1.0 / d[c * (i + 1) - 1];
                    for q in 0..dim {
                        row_a[dim * i + q] = 1.0 / d[c * i + q];
                    }
                }
            }
            if iter > 1 {
                let r = rotated.clone();
                if let Some(v) = ctx.primary_mut() {
                    v.values[..nrows].copy_from_slice(&r[..nrows]);
                }
            }
            let (a, b) = {
                let s = ctx.system().ok_or_else(|| SolverError::fatal(CALLER, "no system"))?;
                (s.matrix.clone(), s.rhs.clone())
            };
            let mut x: Vec<f64> = ctx.primary().expect("checked above").values[..nrows].to_vec();
            let mut residual = vec![0.0; nrows];
            matvec(&a, &x, &mut residual);
            for i in 0..nrows {
                residual[i] = b[i] - residual[i];
            }
            non_lin_error = norm(&residual) / norm(&b);
            ctx.info("StokesSolver", &format!("Residual for nonlinear iterate{:>4}{}", iter - 1, es12_3(non_lin_error)), 3);
            if non_lin_error < nonlinear_tol || iter == nonlinear_iter {
                ctx.back_rotate_nt();
                break;
            }
            let min_tol = if adaptive {
                (rel_tol * non_lin_error).min(base_tol).max(target_tol)
            } else {
                target_tol
            };
            let mut gcr = Gcr { s: Vec::new(), v: Vec::new(), rr: vec![0.0; nrows] };
            let mut snew = vec![0.0; nrows];
            let ppermv: Vec<u32> = ctx.primary().and_then(|v| v.perm.map(|p| p.to_vec())).unwrap_or_default();
            for round in 1..=max_iterations {
                matvec(&a, &x, &mut residual);
                for i in 0..nrows {
                    residual[i] = b[i] - residual[i];
                }
                {
                    let mut pc = SolverCtx::new(&mut *ctx.sim, pix, mesh, log, comm);
                    let ps = pc.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no pressure system"))?;
                    for i in 0..np {
                        ps.rhs[i] = residual[c * (i + 1) - 1] * row_p[i];
                    }
                }
                {
                    let mut vc = SolverCtx::new(&mut *ctx.sim, vix, mesh, log, comm);
                    let vs = vc.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no velocity system"))?;
                    for i in 0..np {
                        for q in 0..dim {
                            vs.rhs[dim * i + q] = residual[c * i + q] * row_a[dim * i + q];
                        }
                    }
                }
                ctx.info("Pressure Preconditioning", "Performing linear solve", 4);
                SolverCtx::new(&mut *ctx.sim, pix, mesh, log, comm).default_solve()?;
                let pvals: Vec<f64> = {
                    let id = ctx.sim.slots[pix].variable.ok_or_else(|| SolverError::fatal(CALLER, "no pressure variable"))?;
                    ctx.sim.vars.by_id(id).values.to_vec()
                };
                let mut temps = vec![0.0; nrows];
                for i in 0..np {
                    temps[c * (i + 1) - 1] = pvals[i];
                }
                let mut tempres = vec![0.0; nrows];
                matvec(&a, &temps, &mut tempres);
                {
                    let mut vc = SolverCtx::new(&mut *ctx.sim, vix, mesh, log, comm);
                    let vs = vc.system_mut().expect("checked above");
                    for i in 0..np {
                        for q in 0..dim {
                            vs.rhs[dim * i + q] -= tempres[c * i + q] * row_a[dim * i + q];
                        }
                    }
                }
                ctx.info("Velocity Preconditioning", "Performing linear solve", 4);
                SolverCtx::new(&mut *ctx.sim, vix, mesh, log, comm).default_solve()?;
                let vvals: Vec<f64> = {
                    let id = ctx.sim.slots[vix].variable.ok_or_else(|| SolverError::fatal(CALLER, "no velocity variable"))?;
                    ctx.sim.vars.by_id(id).values.to_vec()
                };
                for &e in &elems {
                    for &node in ctx.mesh.element(e).nodes {
                        let j = ppermv[node as usize] as usize;
                        snew[c * j + dim] = pvals[j];
                        for q in 0..dim {
                            snew[c * j + q] = vvals[dim * j + q];
                        }
                    }
                }
                let res = gcr.update(&a, &mut x, &b, &mut residual, &snew, round, restart, ctx);
                if res < min_tol {
                    break;
                }
            }
            if let Some(v) = ctx.primary_mut() {
                v.values[..nrows].copy_from_slice(&x);
            }
            rotated = x;
            ctx.back_rotate_nt();
        }
        let _ = non_lin_error;

        if p2p1 {
            p2p1_pressures(ctx, &elems, dim)?;
        }
        Ok(())
    }
}

/// The P2-P1 mid-edge (and quad centre) pressures interpolated from the
/// corners.
fn p2p1_pressures(ctx: &mut SolverCtx<'_>, elems: &[ElemId], dim: usize) -> Result<(), SolverError> {
    let c = dim + 1;
    for &e in elems {
        let el = ctx.mesh.element(e);
        let nodes = el.nodes.to_vec();
        let fam = el.etype.family();
        let (edges, mids): (&[[usize; 2]], &[usize]) = match fam {
            3 => (&[[0, 1], [1, 2], [2, 0]], &[3, 4, 5]),
            4 => (&[[0, 1], [1, 2], [2, 3], [3, 0]], &[4, 5, 6, 7]),
            f => return Err(unsupported(&format!("P2-P1 pressure interpolation on family {f}"))),
        };
        if nodes.len() < mids.last().map_or(0, |m| m + 1) {
            continue;
        }
        let v = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
        let perm = v.perm.map(|p| p.to_vec()).ok_or_else(|| SolverError::fatal(CALLER, "no perm"))?;
        let idx = |k: usize| c * perm[nodes[k] as usize] as usize + dim;
        for (q, ed) in edges.iter().enumerate() {
            let m = idx(mids[q]);
            let (i, j) = (idx(ed[0]), idx(ed[1]));
            v.values[m] = 0.5 * (v.values[i] + v.values[j]);
        }
        if fam == 4 && nodes.len() == 9 {
            let mut res = 0.0;
            for q in 0..4 {
                res += v.values[idx(q)];
            }
            let m = idx(8);
            v.values[m] = 0.25 * res;
        }
    }
    Ok(())
}
