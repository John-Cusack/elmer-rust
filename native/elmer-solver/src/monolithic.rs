//! `Monolithic Slave` (`MergeSlaveSolvers`, ContactUtils.F90:2619): the
//! solver's system and its first slave solver's (assembled, with Dirichlet
//! conditions, by the slave's own visit with its linear solver disabled)
//! become one system in which a node both variables have is one dof, solved
//! with the master's settings and split back.

use std::collections::BTreeMap;
use std::sync::Arc;

use elmer_core::{CsrMatrix, LinearOperator, Params, Value, NO_DOF};

use crate::error::SolverError;
use crate::sim::LinearSystem;
use crate::solver::SolverCtx;

const CALLER: &str = "MergeSlaveSolvers";

/// `CRS_MergeMatrix(A1, A2, C, PermA, PermB, PermC)` (CRSMatrix.F90:2193):
/// the merged permutation keeps `perm1` and numbers the nodes only `perm2`
/// has after it, in node order; a merged row holds `A1`'s entries plus, for
/// a node both have, `A2`'s (summed where the columns meet), columns sorted
/// (`CRS_SortMatrix`). Rows of nodes only `A2` has follow in node order.
fn merge_matrix(a1: &CsrMatrix<f64>, a2: &CsrMatrix<f64>, perm1: &[u32], perm2: &[u32]) -> (CsrMatrix<f64>, Vec<u32>) {
    let n1 = a1.nrows;
    let mut inv1 = vec![0usize; n1];
    let mut inv2 = vec![0usize; a2.nrows];
    for (i, (&p1, &p2)) in perm1.iter().zip(perm2).enumerate() {
        if p1 != NO_DOF {
            inv1[p1 as usize] = i;
        }
        if p2 != NO_DOF {
            inv2[p2 as usize] = i;
        }
    }
    let mut perm = perm1.to_vec();
    let mut next = n1 as u32;
    for (p, (&p1, &p2)) in perm.iter_mut().zip(perm1.iter().zip(perm2)) {
        if p1 == NO_DOF && p2 != NO_DOF {
            *p = next;
            next += 1;
        }
    }
    let row_of = |a: &CsrMatrix<f64>, inv: &[usize], r: usize, into: &mut BTreeMap<u32, f64>| {
        for k in a.row_ptr[r]..a.row_ptr[r + 1] {
            let c = perm[inv[a.col_idx[k] as usize]];
            *into.entry(c).or_insert(0.0) += a.values[k];
        }
    };
    let mut rows: Vec<BTreeMap<u32, f64>> = Vec::with_capacity(next as usize);
    for (r, &node) in inv1.iter().enumerate() {
        let mut m = BTreeMap::new();
        row_of(a1, &inv1, r, &mut m);
        if perm2[node] != NO_DOF {
            row_of(a2, &inv2, perm2[node] as usize, &mut m);
        }
        rows.push(m);
    }
    for (i, (&p1, &p2)) in perm1.iter().zip(perm2).enumerate() {
        if p1 == NO_DOF && p2 != NO_DOF {
            debug_assert_eq!(perm[i] as usize, rows.len());
            let mut m = BTreeMap::new();
            row_of(a2, &inv2, p2 as usize, &mut m);
            rows.push(m);
        }
    }
    let mut c = CsrMatrix {
        nrows: rows.len(),
        ncols: rows.len(),
        row_ptr: vec![0],
        col_idx: Vec::new(),
        values: Vec::new(),
        diag: None,
    };
    let mut diag = Vec::with_capacity(rows.len());
    for (i, m) in rows.into_iter().enumerate() {
        let mut d = None;
        for (col, v) in m {
            if col as usize == i {
                d = Some(c.values.len());
            }
            c.col_idx.push(col);
            c.values.push(v);
        }
        diag.push(d);
        c.row_ptr.push(c.values.len());
    }
    c.diag = diag.into_iter().collect();
    (c, perm)
}

/// `DefaultSolve` with `Monolithic Slave`: merge (`PreSolve`), solve, split
/// (`UnmergeSolutions`). The nonlinear change is that of the merged vector,
/// as in Elmer.
pub(crate) fn solve_merged(ctx: &mut SolverCtx<'_>) -> Result<f64, SolverError> {
    ctx.info(CALLER, "Monolithic treatment of solvers", 4);
    let p = ctx.params().clone();
    let slaves = match p.get("slave solvers") {
        Some(Value::Integer(v)) if !v.is_empty() => v.clone(),
        _ => return ctx.solve_system(),
    };
    if slaves.len() > 1 {
        ctx.warn(CALLER, "Cannot current only deal with one slave solver!");
    }
    let master = ctx.slot;
    let six = usize::try_from(slaves[0])
        .ok()
        .and_then(|id| ctx.sim.slots.iter().position(|s| s.id == id))
        .ok_or_else(|| SolverError::fatal(CALLER, format!("no slave solver {}", slaves[0])))?;
    let id1 = ctx.sim.slots[master]
        .variable
        .ok_or_else(|| SolverError::fatal(CALLER, "Not applicable without a variable"))?;
    let id2 = ctx.sim.slots[six]
        .variable
        .ok_or_else(|| SolverError::fatal(CALLER, "slave solver has no variable"))?;
    let (v1, v2) = (ctx.sim.vars.by_id(id1), ctx.sim.vars.by_id(id2));
    if v1.dofs != 1 || v2.dofs != 1 {
        return Err(SolverError::fatal(
            "MergeSlaveSolver",
            format!("Enabled for 1 dof only: {},{}", v1.dofs, v2.dofs),
        ));
    }
    let (perm1, perm2) = match (v1.perm.clone(), v2.perm.clone()) {
        (Some(a), Some(b)) if a.len() == b.len() => (a, b),
        _ => return Err(SolverError::fatal("CRS_MergeMatrix", "Mismatch in perm size")),
    };
    let alpha = ctx.creal(&p, "slave field scaling")?.unwrap_or(1.0);
    let inv_alpha = 1.0 / alpha;
    let beta: Option<Vec<f64>> = match p.string("slave field offset variable") {
        None => None,
        Some(name) => {
            let b = ctx.sim.vars.get(name).ok_or_else(|| {
                SolverError::fatal("MergeSlaveSolver", format!("Could not find offset variable: {name}"))
            })?;
            if b.perm.as_deref() != Some(&perm2[..]) {
                return Err(SolverError::fatal("MergeSlaveSolver", "The offset field should have same permutation!"));
            }
            Some(b.values.clone())
        }
    };
    let a1 = ctx.sim.slots[master]
        .system
        .take()
        .ok_or_else(|| SolverError::fatal(CALLER, "Not applicable without a matrix"))?;
    let a2 = match ctx.sim.slots[six].system.as_ref() {
        Some(s) => s.clone(),
        None => {
            ctx.sim.slots[master].system = Some(a1);
            return Err(SolverError::fatal(CALLER, "slave solver has no assembled system"));
        }
    };
    let (merged, perm) = merge_matrix(&a1.matrix, &a2.matrix, &perm1, &perm2);
    let n = merged.nrows;

    // MergeRhsAndSolutions.
    let rhs_add: Option<Vec<f64>> = beta.as_ref().map(|b| {
        let mut r = vec![0.0; a2.matrix.nrows];
        a2.matrix.apply(&b[..a2.matrix.ncols.min(b.len())], &mut r);
        r.iter().map(|x| -inv_alpha * x).collect()
    });
    let (vals1, vals2) = (ctx.sim.vars.by_id(id1).values.clone(), ctx.sim.vars.by_id(id2).values.clone());
    let mut rhs = vec![0.0; n];
    let mut vals = vec![0.0; n];
    for (i, &j) in perm.iter().enumerate() {
        if j == NO_DOF {
            continue;
        }
        let j = j as usize;
        let (j1, j2) = (perm1[i], perm2[i]);
        if j1 != NO_DOF {
            rhs[j] = a1.rhs[j1 as usize];
        }
        if j2 != NO_DOF {
            let j2 = j2 as usize;
            rhs[j] += inv_alpha * a2.rhs[j2];
            if let Some(r) = &rhs_add {
                rhs[j] += r[j2];
            }
        }
        let mut c1 = 0.0;
        if j1 != NO_DOF {
            c1 = 1.0;
            vals[j] = vals1[j1 as usize];
        }
        if j2 != NO_DOF {
            let j2 = j2 as usize;
            let mut x2 = vals2[j2];
            if let Some(b) = &beta {
                x2 -= b[j2];
            }
            vals[j] = (vals[j] + inv_alpha * x2) / (c1 + 1.0);
        }
    }

    let perm_arc: Arc<[u32]> = perm.clone().into();
    let saved_perm = {
        let v = ctx.sim.vars.by_id_mut(id1);
        v.values = vals;
        v.perm.replace(perm_arc)
    };
    ctx.sim.slots[master].system = Some(LinearSystem { matrix: merged, rhs });
    let solved = ctx.solve_system();

    // UnmergeSolutions, then the master's own system and field back.
    let vals = {
        let v = ctx.sim.vars.by_id_mut(id1);
        v.perm = saved_perm;
        std::mem::replace(&mut v.values, vals1)
    };
    ctx.sim.slots[master].system = Some(a1);
    let mut new2 = vals2;
    {
        let v1 = ctx.sim.vars.by_id_mut(id1);
        for (i, &j) in perm.iter().enumerate() {
            if j != NO_DOF && perm1[i] != NO_DOF {
                v1.values[perm1[i] as usize] = vals[j as usize];
            }
        }
        v1.values_changed = true;
    }
    for (i, &j) in perm.iter().enumerate() {
        if j != NO_DOF && perm2[i] != NO_DOF {
            new2[perm2[i] as usize] = alpha * vals[j as usize];
        }
    }
    if let Some(b) = &beta {
        for (x, bi) in new2.iter_mut().zip(b) {
            *x += bi;
        }
    }
    let v2 = ctx.sim.vars.by_id_mut(id2);
    v2.values = new2;
    v2.values_changed = true;
    ctx.sim.publish_live();
    solved
}
