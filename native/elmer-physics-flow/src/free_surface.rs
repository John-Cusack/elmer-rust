//! `MoveBoundary` and `PoissonSolve` (fem/src/FreeSurface.F90), FlowSolve's
//! free surfaces: after each nonlinear iteration the nodes of the `Free
//! Surface` boundaries move along one coordinate (`Free Coordinate`, else
//! the direction the boundary's average normal points along most) until the
//! velocity is tangent to the surface, and a Laplace equation in that
//! coordinate, on the mesh as it was, spreads the motion into the bulk with
//! the nodes of every boundary that is not `Free Moving` fixed.
//!
//! Ported: 2D meshes, linear boundary lines (202). Refused: 3D free
//! surfaces and higher-order boundary elements.
//!
//! The solver sees the mesh as `&dyn MeshView`; the runner moves the nodes
//! through [`install_node_mover`].

use std::cell::RefCell;
use std::collections::BTreeSet;
use std::sync::Arc;

use elmer_core::{BasisEval, CsrMatrix, ElemId, ElementBasis, ParamMap, Params, Quadrature, Value, NO_DOF};
use elmer_solver::backend::LinearSolverChoice;
use elmer_solver::{SolverCtx, SolverError};

const CALLER: &str = "MoveBoundary";
/// `AEPS` (Types.F90).
const AEPS: f64 = 10.0 * f64::EPSILON;

/// Replaces the coordinates of every node (`Mesh % Nodes`).
pub type NodeMover = Arc<dyn Fn(&[[f64; 3]])>;

thread_local! {
    static MOVER: RefCell<Option<NodeMover>> = const { RefCell::new(None) };
}

/// Installs the runner's node mover for the solvers run on this thread,
/// returning the previous one.
pub fn install_node_mover(m: Option<NodeMover>) -> Option<NodeMover> {
    MOVER.with(|c| std::mem::replace(&mut *c.borrow_mut(), m))
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("FlowSolver: free surface {what}"))
}

/// Any boundary condition with `Free Surface = True` (FlowSolve.F90:606).
pub fn any_free_surface(ctx: &SolverCtx<'_>) -> bool {
    ctx.sim.model.boundary_conditions.iter().any(|bc| bc.logical("free surface") == Some(true))
}

/// The free-surface boundary elements in mesh order: a BC with `Free
/// Surface`, a parent the solver is active on (`ActiveBoundaryElement`).
fn free_elements(ctx: &SolverCtx<'_>) -> Result<Vec<ElemId>, SolverError> {
    let mut active = vec![false; ctx.mesh.n_bulk()];
    for e in ctx.active_elements() {
        active[e.0 as usize] = true;
    }
    let mut out = Vec::new();
    for e in ctx.boundary_elements() {
        if ctx.bc_params(e).and_then(|b| b.logical("free surface")) != Some(true) {
            continue;
        }
        let el = ctx.mesh.element(e);
        if !el.parents.iter().flatten().any(|p| active.get(p.0 as usize) == Some(&true)) {
            continue;
        }
        if el.etype.0 != 202 {
            return Err(unsupported(&format!("on boundary elements of type {} (only 2D linear lines)", el.etype.0)));
        }
        out.push(e);
    }
    Ok(out)
}

/// `CheckNormalDirection(Boundary, Nrm, x, y, z, Turned)` without `OutBody`:
/// the normal turns when it points towards the centre of the parent (the
/// one with the larger body id of two).
fn turned(ctx: &SolverCtx<'_>, coords: &[[f64; 3]], e: ElemId, nrm: [f64; 3], x: [f64; 3]) -> Result<bool, SolverError> {
    let el = ctx.mesh.element(e);
    let parent = match el.parents {
        [Some(l), Some(r)] => {
            let body = |p: ElemId| ctx.sim.model.body_of_tag(ctx.mesh.element(p).tag);
            match body(l).cmp(&body(r)) {
                std::cmp::Ordering::Greater => l,
                std::cmp::Ordering::Less => r,
                std::cmp::Ordering::Equal => return Err(unsupported("between two elements of one body")),
            }
        }
        [Some(l), None] => l,
        [None, Some(r)] => r,
        [None, None] => return Ok(false),
    };
    let pe = ctx.mesh.element(parent);
    // NodalBasisFunctions at the reference centre.
    let basis: Vec<f64> = match pe.etype.0 {
        404 => vec![0.25; 4],
        303 => {
            let t = 1.0 / 3.0;
            vec![1.0 - t - t, t, t]
        }
        t => return Err(unsupported(&format!("next to a parent element of type {t}"))),
    };
    let mut d = [0.0; 3];
    for (a, da) in d.iter_mut().enumerate() {
        let c = pe.nodes.iter().zip(&basis).fold(0.0, |s, (&n, b)| s + b * coords[n as usize][a]);
        *da = c - x[a];
    }
    Ok(d[0] * nrm[0] + d[1] * nrm[1] + d[2] * nrm[2] > 0.0)
}

/// `MoveBoundary(Model, Relax)`.
pub fn move_boundary(
    ctx: &mut SolverCtx<'_>,
    relax: f64,
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
) -> Result<(), SolverError> {
    if ctx.mesh.dim() != 2 {
        return Err(unsupported("in 3D"));
    }
    let p = ctx.params().clone();
    let feps = p.real("free surface convergence tolerance").unwrap_or(1.0e-12);
    let fiter = p.integer("free surface max iterations").unwrap_or(100);
    let (dofs, perm, velo) = {
        let f = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no flow variable"))?;
        (f.dofs, f.perm.map(|p| p.to_vec()), f.values.to_vec())
    };
    let nn = ctx.mesh.n_nodes();
    // Velocity 1 and 2 at a node.
    let uv = |k: u32| -> Result<(f64, f64), SolverError> {
        let s = perm.as_ref().map_or(k, |p| p[k as usize]);
        if s == NO_DOF {
            return Err(SolverError::fatal(CALLER, format!("free surface node {} outside the flow system", k + 1)));
        }
        let s = s as usize * dofs;
        Ok((velo[s], velo[s + 1]))
    };
    let orig: Vec<[f64; 3]> = (0..nn as u32).map(|n| ctx.mesh.coord(n)).collect();
    let mut coords = orig.clone();
    let elems = free_elements(ctx)?;
    let n_bcs = ctx.sim.model.boundary_conditions.len();
    let mut visited = vec![false; nn];
    let mut turn = vec![false; nn];
    let mut avg = vec![[0.0f64; 3]; n_bcs];
    // FirstDerivative1D of the linear line's basis.
    let dl = [-0.5, 0.5];
    let normal = |coords: &[[f64; 3]], nodes: &[u32]| -> [f64; 3] {
        let sum = |a: usize| nodes.iter().zip(&dl).fold(0.0, |s, (&n, d)| s + coords[n as usize][a] * d);
        [-sum(1), sum(0), 0.0]
    };
    // Whether node i of the line is downstream (skipped) of the other one.
    let skip = |nodes: &[u32], i: usize, ux: f64, uy: f64| -> bool {
        let (k, j) = (nodes[i] as usize, nodes[1 - i] as usize);
        let (x1, y1) = (orig[j][0] - orig[k][0], orig[j][1] - orig[k][1]);
        ux * x1 + uy * y1 > 0.0
    };

    for &e in &elems {
        let nodes = ctx.mesh.element(e).nodes.to_vec();
        let bc = ctx.bc_id(e).expect("free_elements has a BC") - 1;
        for i in 0..nodes.len() {
            let k = nodes[i];
            let (ux, uy) = uv(k)?;
            if skip(&nodes, i, ux, uy) || visited[k as usize] {
                continue;
            }
            visited[k as usize] = true;
            let mut nrm = normal(&coords, &nodes);
            if turned(ctx, &coords, e, nrm, coords[k as usize])? {
                nrm = nrm.map(|c| -c);
                turn[k as usize] = true;
            }
            let r = (nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]).sqrt();
            for a in 0..3 {
                avg[bc][a] += nrm[a].abs() / r;
            }
        }
    }

    let (mut xmoved, mut ymoved) = (false, false);
    for ii in 1..=fiter {
        let mut s = 0.0;
        visited.iter_mut().for_each(|v| *v = false);
        for &e in &elems {
            let bcp = ctx.bc_params(e).expect("free_elements has a BC");
            let which = match bcp.integer("free coordinate") {
                Some(w) if (1..=2).contains(&w) => w as usize,
                Some(w) => return Err(SolverError::fatal("ListGetInteger", format!("Free Coordinate {w} out of range"))),
                None => {
                    let nrm = avg[ctx.bc_id(e).expect("free_elements has a BC") - 1];
                    let mut w = 1;
                    for i in 2..=3 {
                        if nrm[i - 1] > nrm[w - 1] {
                            w = i;
                        }
                    }
                    w
                }
            };
            let nodes = ctx.mesh.element(e).nodes.to_vec();
            // GetElementNodes: the coordinates as the element starts.
            let bn: Vec<[f64; 3]> = nodes.iter().map(|&n| coords[n as usize]).collect();
            for i in 0..nodes.len() {
                let k = nodes[i] as usize;
                let (ux, uy) = uv(nodes[i])?;
                if skip(&nodes, i, ux, uy) || visited[k] {
                    continue;
                }
                visited[k] = true;
                let sum = |a: usize| bn.iter().zip(&dl).fold(0.0, |s, (x, d)| s + x[a] * d);
                let mut nrm = [-sum(1), sum(0), 0.0];
                if turn[k] {
                    nrm = nrm.map(|c| -c);
                }
                let r = ux * nrm[0] + uy * nrm[1];
                if which == 2 {
                    if ux.abs() > AEPS {
                        if turn[k] {
                            coords[k][1] -= r / (ux * dl[i]);
                        } else {
                            coords[k][1] += r / (ux * dl[i]);
                        }
                        ymoved = true;
                    }
                } else if uy.abs() > AEPS {
                    if turn[k] {
                        coords[k][0] += r / (uy * dl[i]);
                    } else {
                        coords[k][0] -= r / (uy * dl[i]);
                    }
                    xmoved = true;
                }
                let q = (ux * nrm[0] + uy * nrm[1]) / (nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]).sqrt();
                s += q * q;
            }
        }
        let s = s.sqrt();
        ctx.info("FreeSurface", &format!("Iter: {ii:>4} Free surface Residual: {}", fortran_e21_12(s)), 4);
        if s < feps {
            break;
        }
    }

    for (moved, a) in [(xmoved, 0), (ymoved, 1)] {
        if !moved {
            continue;
        }
        let mut sol: Vec<f64> = coords.iter().map(|c| c[a]).collect();
        poisson_solve(ctx, &orig, &mut sol, a, basis, quad)?;
        for (c, (x, o)) in coords.iter_mut().zip(sol.iter().zip(&orig)) {
            c[a] = if relax != 1.0 { relax * x + (1.0 - relax) * o[a] } else { *x };
        }
    }
    if coords != orig {
        let m = MOVER.with(|c| c.borrow().clone()).ok_or_else(|| {
            SolverError::NotWired("FlowSolver: moving the free surface (no node mover installed on this thread)".into())
        })?;
        m(&coords);
    }
    Ok(())
}

/// Fortran `E21.12`: `   0.123456789012E-05`.
fn fortran_e21_12(x: f64) -> String {
    if x == 0.0 {
        return format!("{:>21}", "0.000000000000E+00");
    }
    // 12 significant digits d.ddddddddddd e E -> 0.dddddddddddd E(E+1).
    let s = format!("{:.11e}", x.abs());
    let (m, e) = s.split_once('e').expect("exponent");
    let e: i32 = e.parse().expect("integer exponent");
    let digits: String = m.chars().filter(|c| c.is_ascii_digit()).collect();
    let e = e + 1;
    let sign = if x < 0.0 { "-" } else { "" };
    format!("{:>21}", format!("{sign}0.{digits}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs()))
}

/// `PoissonSolve(Model, X, Y, Z, Solution, Moved)`: `d/dx_moved` Laplace
/// matrix of every bulk element on the coordinates `orig`, the identity
/// numbering of all nodes, the rows of the nodes of every boundary that is
/// not `Free Moving` replaced by the identity with the current coordinate,
/// solved by CGS + ILU0 to 1e-9 from `sol`.
fn poisson_solve(
    ctx: &mut SolverCtx<'_>,
    orig: &[[f64; 3]],
    sol: &mut [f64],
    moved: usize,
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
) -> Result<(), SolverError> {
    let nn = orig.len();
    let (nb, nbnd) = (ctx.mesh.n_bulk(), ctx.mesh.n_boundary());
    // CreateMatrix without an equation: every element's nodes.
    let mut rows: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); nn];
    for t in 0..nb + nbnd {
        let nodes = ctx.mesh.element(ElemId(t as u32)).nodes;
        for &r in nodes {
            rows[r as usize].extend(nodes.iter().copied());
        }
    }
    let mut row_ptr = vec![0usize];
    let mut col_idx = Vec::new();
    for r in &rows {
        col_idx.extend(r.iter().copied());
        row_ptr.push(col_idx.len());
    }
    let mut a = CsrMatrix {
        nrows: nn,
        ncols: nn,
        values: vec![0.0; col_idx.len()],
        row_ptr,
        col_idx,
        diag: None,
    };
    let at = |a: &CsrMatrix<f64>, r: usize, c: u32| -> usize {
        let s = &a.col_idx[a.row_ptr[r]..a.row_ptr[r + 1]];
        a.row_ptr[r] + s.binary_search(&c).expect("in the pattern")
    };
    let mut be = BasisEval::default();
    for t in 0..nb {
        let el = ctx.mesh.element(ElemId(t as u32));
        let n = el.nodes.len();
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| orig[k as usize]).collect();
        let mut local = vec![0.0; n * n];
        for qp in quad.gauss_points(el.etype, None) {
            basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let s = be.det_j * qp.weight;
            for p in 0..n {
                for q in 0..n {
                    local[p * n + q] += s * (be.dbasis_dx[p][moved] * be.dbasis_dx[q][moved]);
                }
            }
        }
        for p in 0..n {
            for q in 0..n {
                let ix = at(&a, el.nodes[p] as usize, el.nodes[q]);
                a.values[ix] += local[p * n + q];
            }
        }
    }
    let diag: Vec<usize> = (0..nn).map(|r| at(&a, r, r as u32)).collect();
    a.diag = Some(diag);
    let mut b = vec![0.0; nn];
    for t in nb..nb + nbnd {
        let e = ElemId(t as u32);
        let Some(bc) = ctx.bc_params(e) else { continue };
        if bc.logical("free moving") == Some(true) {
            continue;
        }
        for &k in ctx.mesh.element(e).nodes {
            let k = k as usize;
            b[k] = sol[k];
            for v in &mut a.values[a.row_ptr[k]..a.row_ptr[k + 1]] {
                *v = 0.0;
            }
            let ix = at(&a, k, k as u32);
            a.values[ix] = 1.0;
        }
    }
    // IterSolver is called directly: the iterative path, no scaling.
    let mut p = ParamMap::new();
    p.insert("Linear System Solver", Value::Str("Iterative".into()));
    p.insert("Linear System Iterative Method", Value::Str("CGS".into()));
    p.insert("Linear System Max Iterations", Value::Integer(vec![500]));
    p.insert("Linear System Convergence Tolerance", Value::real(1.0e-9));
    p.insert("Linear System Preconditioning", Value::Str("ILU0".into()));
    p.insert("Linear System Residual Output", Value::Integer(vec![1]));
    match ctx.sim.backends.linear.select_for(&p, &a, ctx.log)? {
        LinearSolverChoice::Iterative { mut solver, precond } => {
            let m: Box<dyn elmer_core::Preconditioner<f64> + Send> = match precond {
                Some(f) => f(&p, &a)?,
                None => Box::new(elmer_core::IdentityPrecond),
            };
            solver.solve(&a, &*m, &b, sol, ctx.comm).map_err(SolverError::from)?;
        }
        LinearSolverChoice::Direct(_) => return Err(SolverError::fatal("PoissonSolve", "CGS selected a direct solver")),
    }
    Ok(())
}
