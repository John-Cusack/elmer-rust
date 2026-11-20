//! `FollowCurvedBoundary(Model, Mesh, SetP = .TRUE.)` (MeshBasics.F90:2893,
//! from a solver's `Follow P Curvature`, e.g. StatElecSolve.F90:363): on
//! each boundary element with p dofs of a BC following a circle, the
//! displacement of its straight geometry onto the circle is fitted in the
//! element's p basis (the nodes held, `LUDecomp` / `LUSolve`), and the edge
//! dofs' coefficients become their coordinates (`Mesh % Nodes` past the
//! nodes). Elements built on those dofs then map isoparametrically.
//!
//! The circle is the one the mesh-load follow kept in the BC (`Cylinder
//! Center`, `Cylinder Radius`). Refused: the other `Follow` modes.

use std::collections::BTreeMap;

use elmer_core::{BasisEval, ElemId, Params, Value};

use crate::error::SolverError;
use crate::solver::SolverCtx;

const CALLER: &str = "SetCurvedBoundary";

/// `LUDecomp` (LinearAlgebra.F90:191), column-major `a[i + j*n]`, the
/// diagonal replaced by its reciprocal.
pub(crate) fn lu_decomp(a: &mut [f64], n: usize, pivot: &mut [usize]) -> bool {
    let ix = |i: usize, j: usize| i + j * n;
    for i in 0..n {
        let mut j = i;
        for k in i + 1..n {
            if a[ix(i, k)].abs() > a[ix(i, j)].abs() {
                j = k;
            }
        }
        if a[ix(i, j)].abs() == 0.0 {
            return false;
        }
        pivot[i] = j;
        if j != i {
            for k in 0..=i {
                a.swap(ix(k, j), ix(k, i));
            }
        }
        for k in i + 1..n {
            a[ix(i, k)] /= a[ix(i, i)];
        }
        for k in i + 1..n {
            if j != i {
                a.swap(ix(k, i), ix(k, j));
            }
            for l in i + 1..n {
                a[ix(k, l)] -= a[ix(k, i)] * a[ix(i, l)];
            }
        }
    }
    pivot[n - 1] = n - 1;
    for i in 0..n {
        if a[ix(i, i)].abs() == 0.0 {
            return false;
        }
        a[ix(i, i)] = 1.0 / a[ix(i, i)];
    }
    true
}

/// `LUSolve(n, A, x, pivot)` (LinearAlgebra.F90:132).
pub(crate) fn lu_solve(a: &[f64], n: usize, x: &mut [f64], pivot: &[usize]) {
    let ix = |i: usize, j: usize| i + j * n;
    for i in 0..n {
        let mut s = x[i];
        for j in 0..i {
            s -= a[ix(i, j)] * x[j];
        }
        x[i] = a[ix(i, i)] * s;
    }
    for i in (0..n).rev() {
        let mut s = x[i];
        for j in i + 1..n {
            s -= a[ix(i, j)] * x[j];
        }
        x[i] = s;
    }
    for i in (0..n).rev() {
        if pivot[i] != i {
            x.swap(i, pivot[i]);
        }
    }
}

/// `FollowCurvedBoundary(.., SetP = .TRUE.)` for the running solver.
pub fn follow_p_curvature(ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
    let model = &ctx.sim.model;
    if !model.boundary_conditions.iter().any(|bc| bc.0.keys().any(|k| k.starts_with("follow"))) {
        return Ok(());
    }
    let Some(d) = ctx.solver_dofs() else { return Ok(()) };
    let dim = ctx.mesh.dim().max(2);
    // `Mesh % Nodes % x(Indexes(n+1:nd)) = x(n+1:nd)`: set by slot here,
    // stored by global index below.
    let mut geometry: BTreeMap<u32, [f64; 3]> = BTreeMap::new();
    let nb = ctx.mesh.n_bulk();
    for (b, bc) in model.boundary_conditions.iter().enumerate() {
        if bc.logical("follow circle boundary") != Some(true) {
            if bc.0.keys().any(|k| k.starts_with("follow") && k.ends_with("boundary")) {
                return Err(SolverError::Unsupported(format!("Boundary Condition {}: Follow P Curvature with a boundary other than a circle", b + 1)));
            }
            continue;
        }
        let center = match bc.get("cylinder center") {
            Some(Value::Real { data, .. }) if data.len() >= 2 => [data[0], data[1]],
            _ => return Err(SolverError::fatal(CALLER, "the circle of the mesh-load follow is missing")),
        };
        let r = bc.real("cylinder radius").ok_or_else(|| SolverError::fatal(CALLER, "no Cylinder Radius"))?;
        let mut extra = Vec::new();
        let mut ev = BasisEval::default();
        for t in nb..nb + ctx.mesh.n_boundary() {
            let e = ElemId(t as u32);
            if ctx.bc_id(e) != Some(b + 1) {
                continue;
            }
            d.boundary_extra_slots(e, &mut extra);
            if extra.is_empty() {
                continue;
            }
            let etype = ctx.mesh.element(e).etype;
            let nodes = d.element_nodes(e);
            let n = nodes.len();
            let nd = n + extra.len();
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let basis = d.element_basis_scalar(e).map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
            let mut mass = vec![0.0f64; nd * nd]; // MASS(i, j) at i + j*nd
            let mut force = vec![vec![0.0f64; nd]; dim];
            for qp in &rule {
                basis.eval(etype, &xyz, qp.uvw, &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let w = qp.weight * ev.det_j;
                let mut c = [0.0f64; 3];
                for (a, ca) in c.iter_mut().enumerate() {
                    *ca = (0..n).fold(0.0, |s, i| s + xyz[i][a] * ev.basis[i]);
                }
                let c0 = c;
                // Mode 1: radially onto the circle.
                let (cx, cy) = (c[0] - center[0], c[1] - center[1]);
                let rat = r / (cx * cx + cy * cy).sqrt();
                c[0] = rat * cx + center[0];
                c[1] = rat * cy + center[1];
                let disp = [c[0] - c0[0], c[1] - c0[1], c[2] - c0[2]];
                for q in 0..nd {
                    for i in 0..nd {
                        mass[i + q * nd] += w * ev.basis[i] * ev.basis[q];
                    }
                }
                for (a, f) in force.iter_mut().enumerate() {
                    for i in 0..nd {
                        f[i] += w * ev.basis[i] * disp[a];
                    }
                }
            }
            for i in 0..n {
                for j in 0..nd {
                    mass[i + j * nd] = 0.0;
                }
                mass[i + i * nd] = 1.0;
                for f in force.iter_mut() {
                    f[i] = 0.0;
                }
            }
            let mut pivot = vec![0usize; nd];
            if !lu_decomp(&mut mass, nd, &mut pivot) {
                return Err(SolverError::fatal(CALLER, "LU-decomposition fails"));
            }
            for (a, f) in force.iter_mut().enumerate() {
                lu_solve(&mass, nd, f, &pivot);
                for (k, &s) in extra.iter().enumerate() {
                    geometry.entry(s).or_insert([0.0; 3])[a] = f[n + k];
                }
            }
        }
    }
    ctx.info(CALLER, &format!("p-dof coordinates set: {}", geometry.len()), 10);
    let perm = d.perm().perm.clone();
    let coords = ctx.sim.p_coords.get_or_insert_with(Vec::new);
    if coords.len() < perm.len() {
        coords.resize(perm.len(), [0.0; 3]);
    }
    for (g, &s) in perm.iter().enumerate() {
        if let Some(c) = geometry.get(&s) {
            coords[g] = *c;
        }
    }
    Ok(())
}
