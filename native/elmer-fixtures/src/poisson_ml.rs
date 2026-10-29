//! `PoissonML/PoissonML.f90`: `SUBROUTINE PoissonSolver`, a matrix-free
//! ("matrix-less") Poisson solve on p elements (`No Matrix = True`).
//!
//! The routine keeps every element's condensed stiffness matrix and force
//! (bulk elements, then all boundary elements, those without a BC with zero
//! matrices), sums the matrix diagonals per dof, and applies the Dirichlet
//! conditions inside the element matrices: for each dof of a boundary
//! element whose BC gives the variable a constant value, in every element
//! containing it, the column times the value moves to the force (nodal dofs
//! only), the row and column are cleared and the diagonal becomes the dof's
//! averaged diagonal. It then solves with its own Jacobi-preconditioned CG
//! (`PCG`), whose matrix-vector product loops over the element matrices,
//! into the solver variable; nothing goes through `DefaultSolve`.
//!
//! `PCG` stops on `SUM(r*r) < eps**2 * bnrm**2` with `bnrm = SUM(b*b)`
//! (already squared, as the Fortran writes it) and prints every 25th
//! iteration and the final residual. The dofs, basis, Gauss rules and node
//! orders come from the solver's dof layout (`elmer-pdofs`).

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementType, ParamMap, Params, QuadPoint, NO_DOF};
use elmer_solver::fortran_fmt::list_int;
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::linalg::condensate;
use crate::udf::{list_real, stdout_sink, PrintSink};

const CALLER: &str = "PoissonSolver";

/// One stored element (`ed_t`): its dof slots and row-major matrix.
struct Stored {
    dofs: Vec<usize>,
    stiff: Vec<f64>,
    force: Vec<f64>,
}

pub struct PoissonMl {
    out: PrintSink,
}

/// [`SolverFactory`](elmer_solver::SolverFactory): prints to standard output.
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(PoissonMl { out: stdout_sink() }))
}

/// `LocalMatrix` (`with_stiffness`) or `LocalMatrixBC` over `m` dofs:
/// `STIFF += s*detJ * MATMUL(dBasisdx, TRANSPOSE(dBasisdx))`,
/// `FORCE += s*detJ * LoadAtIP * Basis`.
pub fn local_matrix(
    etype: ElementType,
    xyz: &[[f64; 3]],
    load: &[f64],
    m: usize,
    with_stiffness: bool,
    basis: &dyn ElementBasis,
    rule: &[QuadPoint],
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let n = xyz.len();
    let mut stiff = vec![0.0; m * m];
    let mut force = vec![0.0; m];
    let mut ev = BasisEval::default();
    for qp in rule {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        if ev.basis.len() < m {
            return Err(format!("ElementInfo gave {} basis functions for {m} dofs", ev.basis.len()));
        }
        let (b, g) = (&ev.basis[..m], &ev.dbasis_dx[..m]);
        let load_at = b[..n].iter().zip(load).fold(0.0, |s, (bp, l)| s + bp * l);
        let w = qp.weight * ev.det_j;
        if with_stiffness {
            for p in 0..m {
                for q in 0..m {
                    stiff[p * m + q] += w * (g[p][0] * g[q][0] + g[p][1] * g[q][1] + g[p][2] * g[q][2]);
                }
            }
        }
        for (f, bp) in force.iter_mut().zip(b) {
            *f += w * load_at * bp;
        }
    }
    Ok((stiff, force))
}

/// `GetReal(list, key)` at the element's nodes in `order`; zeros when absent.
fn nodal(ctx: &SolverCtx<'_>, list: Option<&ParamMap>, key: &str, e: ElemId, order: &[u32]) -> Result<Vec<f64>, SolverError> {
    let mesh_nodes = ctx.mesh.element(e).nodes;
    let Some(values) = list.map(|p| ctx.real_at_nodes(p, key, e)).transpose()?.flatten() else {
        return Ok(vec![0.0; order.len()]);
    };
    order
        .iter()
        .map(|n| {
            mesh_nodes
                .iter()
                .position(|m| m == n)
                .map(|i| values[i])
                .ok_or_else(|| SolverError::fatal(CALLER, "basis node order is not the element's node set"))
        })
        .collect()
}

/// gfortran's `Ew.d` (`0.ddddE+xx`).
fn fortran_e(x: f64, w: usize, d: usize) -> String {
    if x == 0.0 {
        return format!("{:>w$}", format!("0.{}E+00", "0".repeat(d)));
    }
    let s = format!("{:.*e}", d - 1, x.abs());
    let (m, e) = s.split_once('e').unwrap_or((&s, "0"));
    let e: i32 = e.parse::<i32>().unwrap_or(0) + 1;
    let digits: String = m.chars().filter(char::is_ascii_digit).collect();
    let sign = if x < 0.0 { "-" } else { "" };
    format!("{:>w$}", format!("{sign}0.{digits}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs()))
}

/// `rmv`: `v = sum over elements of stiff * u(dofs)`.
fn matvec(ed: &[Stored], u: &[f64], v: &mut [f64]) {
    v.iter_mut().for_each(|x| *x = 0.0);
    for el in ed {
        let m = el.dofs.len();
        for i in 0..m {
            let x = (0..m).fold(0.0, |s, j| s + el.stiff[i * m + j] * u[el.dofs[j]]);
            v[el.dofs[i]] += x;
        }
    }
}

fn dot(p: &[f64], q: &[f64]) -> f64 {
    p.iter().zip(q).fold(0.0, |s, (a, b)| s + a * b)
}

impl PoissonMl {
    /// `PCG(n, x, b, eps, maxiter)`.
    fn pcg(&self, ed: &[Stored], diag: &[f64], x: &mut [f64], b: &[f64], eps: f64, maxiter: i32) {
        let n = x.len();
        let bnrm = dot(b, b);
        let eps2 = eps.powi(2) * bnrm.powi(2);
        let mut r = vec![0.0; n];
        matvec(ed, x, &mut r);
        r.iter_mut().zip(b).for_each(|(ri, bi)| *ri = bi - *ri);
        let mut residual = dot(&r, &r);
        if residual < eps2 {
            return;
        }
        let (mut p, mut q, mut z) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        let mut oldrho = 0.0;
        let mut iter = 0;
        for it in 1..=maxiter {
            iter = it;
            z.iter_mut().zip(&r).zip(diag).for_each(|((zi, ri), di)| *zi = ri / di);
            let rho = dot(&r, &z);
            if it == 1 {
                p.copy_from_slice(&z);
            } else {
                let beta = rho / oldrho;
                p.iter_mut().zip(&z).for_each(|(pi, zi)| *pi = zi + beta * *pi);
            }
            oldrho = rho;
            matvec(ed, &p, &mut q);
            let alpha = rho / dot(&p, &q);
            x.iter_mut().zip(&p).for_each(|(xi, pi)| *xi += alpha * pi);
            r.iter_mut().zip(&q).for_each(|(ri, qi)| *ri += -alpha * qi);
            residual = dot(&r, &r);
            if it % 25 == 0 {
                (self.out)(&format!("{it:>8}{}", fortran_e((residual / bnrm.powi(2)).sqrt(), 11, 4)));
            }
            if residual < eps2 {
                break;
            }
        }
        // After a loop that ran to the end, Fortran's DO variable is maxiter+1.
        if residual >= eps2 {
            iter = maxiter + 1;
        }
        matvec(ed, x, &mut r);
        r.iter_mut().zip(b).for_each(|(ri, bi)| *ri = bi - *ri);
        let residual = dot(&r, &r).sqrt();
        // `WRITE(*,*) iter, residual/bnrm`: the record's blank and the integer
        // fill twelve columns.
        (self.out)(&format!("{}{}", list_int(iter), list_real(residual / bnrm)));
    }
}

impl Solver for PoissonMl {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if ctx.primary().is_none_or(|v| v.dofs != 1) {
            return Err(SolverError::Unsupported("PoissonML PoissonSolver: a scalar Variable".into()));
        }
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::NotWired(
                "PoissonML PoissonSolver: its p elements need a dof layout (elmer-pdofs in Backends::dofs)".into(),
            ));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let fatal = |m: String| SolverError::fatal(CALLER, m);
        let n_slots = ctx.primary().map_or(0, |v| v.values.len());
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let nb_first = ctx.mesh.n_bulk() as u32;
        let boundary: Vec<ElemId> = (nb_first..nb_first + ctx.mesh.n_boundary() as u32).map(ElemId).collect();
        let mut ed: Vec<Stored> = Vec::with_capacity(active.len() + boundary.len());
        let mut slots = Vec::new();
        let d = ctx.solver_dofs().ok_or_else(|| fatal("no dof layout".into()))?;

        let mut element = |e: ElemId, bulk: bool| -> Result<Stored, SolverError> {
            let etype = ctx.mesh.element(e).etype;
            d.element_slots(e, &mut slots);
            if slots.contains(&NO_DOF) {
                return Err(fatal(format!("element {} has a dof outside the system", e.0 + 1)));
            }
            let nd = slots.len();
            let order = d.element_nodes(e);
            let xyz: Vec<[f64; 3]> = order.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let basis = d.element_basis(e).map_err(fatal)?;
            let rule = d.element_rule(e).map_err(fatal)?;
            let (stiff, force) = if bulk {
                let nb = d.local_bubbles(e);
                let bf = ctx.body_params(e).and_then(|b| b.body_force);
                let load = nodal(ctx, bf, "source", e, &order)?;
                let (mut k, mut f) = local_matrix(etype, &xyz, &load, nd + nb, true, basis.as_ref(), &rule).map_err(fatal)?;
                condensate(nd, nb, &mut k, &mut f).map_err(fatal)?;
                (k, f)
            } else {
                match ctx.bc_params(e) {
                    Some(bc) => {
                        let load = nodal(ctx, Some(bc), "flux", e, &order)?;
                        local_matrix(etype, &xyz, &load, nd, false, basis.as_ref(), &rule).map_err(fatal)?
                    }
                    None => (vec![0.0; nd * nd], vec![0.0; nd]),
                }
            };
            Ok(Stored { dofs: slots.iter().map(|&s| s as usize).collect(), stiff, force })
        };
        for &e in &active {
            ed.push(element(e, true)?);
        }
        for &e in &boundary {
            ed.push(element(e, false)?);
        }

        // np(j): elements containing dof j and the summed diagonal.
        let mut elems: Vec<Vec<usize>> = vec![Vec::new(); n_slots];
        let mut diag = vec![0.0; n_slots];
        for (t, el) in ed.iter().enumerate() {
            let m = el.dofs.len();
            for (i, &j) in el.dofs.iter().enumerate() {
                elems[j].push(t);
                diag[j] += el.stiff[i * m + i];
            }
        }

        // Dirichlet conditions inside the element matrices.
        let name = ctx.variable_name().unwrap_or_default().to_string();
        for (k, &e) in boundary.iter().enumerate() {
            let Some(bc) = ctx.bc_params(e) else { continue };
            let Some(d_val) = ctx.creal(bc, &name)? else { continue };
            let n_nodes = ctx.mesh.element(e).nodes.len();
            let inds = ed[active.len() + k].dofs.clone();
            for (j, &dof) in inds.iter().enumerate() {
                let dg = diag[dof] / elems[dof].len() as f64;
                for &m in &elems[dof] {
                    let el = &mut ed[m];
                    let sz = el.dofs.len();
                    let l = el.dofs.iter().position(|&x| x == dof).ok_or_else(|| fatal("STOP l".into()))?;
                    if j < n_nodes {
                        for r in 0..sz {
                            el.force[r] -= el.stiff[r * sz + l] * d_val;
                        }
                        el.force[l] = d_val * dg;
                    } else {
                        el.force[l] = 0.0;
                    }
                    for r in 0..sz {
                        el.stiff[r * sz + l] = 0.0;
                        el.stiff[l * sz + r] = 0.0;
                    }
                    el.stiff[l * sz + l] = dg;
                }
            }
        }

        let mut b = vec![0.0; n_slots];
        for el in &ed {
            for (j, &k) in el.dofs.iter().enumerate() {
                b[k] += el.force[j];
            }
        }
        let eps = ctx.creal(ctx.params(), "linear system convergence tolerance")?.unwrap_or(0.0);
        let maxi = ctx.params().integer("linear system max iterations").unwrap_or(0);
        let mut x: Vec<f64> = ctx.primary().map(|v| v.values.to_vec()).unwrap_or_default();
        let start = std::time::Instant::now();
        self.pcg(&ed, &diag, &mut x, &b, eps, maxi);
        let t = start.elapsed().as_secs_f64();
        (self.out)(&format!(" Linear System Timing: {}{}", list_real(t), list_real(t)));
        let v = ctx.primary_mut().ok_or_else(|| fatal("no solver variable".into()))?;
        v.values.copy_from_slice(&x);
        Ok(())
    }
}
