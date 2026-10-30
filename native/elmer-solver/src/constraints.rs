//! Linear constraints `C x = d` solved with Lagrange multipliers.
//!
//! `Apply Integral BCs` (`GenerateProjectors` -> `IntegralProjector`,
//! MeshBasics.F90:3353, then `GenerateConstraintMatrix`,
//! ProjectorUtils.F90:440) turns every boundary condition or body force with
//! `Integral BC` and an `Integral BC <var>` value (`Integral BC
//! <component>` per component of a vector field) into one row per such
//! value: the integral of each basis function over its elements, equal to
//! the value. One with
//! `Flux Integral BC` gets `RobinProjector`'s row, the flux of its Robin
//! terms, and the multiplier is exported (`Lagrange Multiplier <equation>`). A solver
//! may instead hand over its own rows ([`SolverCtx::set_constraint_matrix`],
//! Elmer's initial `Matrix % ConstraintMatrix`). Either way the system is
//! solved as `[A C^T; C 0]` (`SolveWithLinearRestriction`,
//! SolveCore.F90:5561), the multipliers included in the nonlinear norm and
//! change as in Elmer.

use std::collections::BTreeMap;

use elmer_core::{CsrMatrix, ElemId, Logger, ParamMap, Params, NO_DOF};

use crate::error::SolverError;
use crate::fortran_fmt::{es, i2s};
use crate::solver::SolverCtx;

/// Rows of linear constraints over a solver's dofs (Elmer's
/// `ConstraintMatrix`): `matrix` has one row per constraint, its columns are
/// equation indices (`perm[node] * dofs + k`), `rhs` the constrained values.
/// `inv_perm[i]` is the dof row `i` stands for (possibly plus a multiple of
/// the matrix size, the projector's row offset, as Elmer's `InvPerm`), or
/// [`NO_DOF`]: when that dof has a Dirichlet condition the constraint is
/// dropped and its multiplier set to zero.
#[derive(Clone, Debug)]
pub struct ConstraintMatrix {
    pub matrix: CsrMatrix<f64>,
    pub rhs: Vec<f64>,
    pub inv_perm: Vec<u32>,
}

impl SolverCtx<'_> {
    /// Sets (or with `None` clears) this solver's own constraint rows, solved
    /// with the system by every later `default_solve`.
    pub fn set_constraint_matrix(&mut self, c: Option<ConstraintMatrix>) {
        self.sim.slots[self.slot].own_constraints = c;
    }
}

/// One `IntegralProjector`: the summed `weight * Basis(p)` per node, in node
/// order, and the representative node (`InvPerm(1)`).
struct Projector {
    entries: BTreeMap<u32, f64>,
    node: u32,
}

/// A boundary's (or body force's) projector: an integral condition, or
/// `RobinProjector`'s one row per component over equation indices.
enum Proj {
    Integral(Projector),
    Robin(Vec<BTreeMap<u32, f64>>),
}

/// `RobinProjector` (ProjectorUtils.F90:296): for the dofs of the
/// boundary's (body force's) elements, minus the boundary part of their
/// rows (`A % Values - A % BulkValues`) over the same dofs, summed per
/// component; the flux the Robin terms carry.
fn robin_projector(ctx: &SolverCtx<'_>, bc: usize, nbc: usize, created: bool) -> Result<Vec<BTreeMap<u32, f64>>, SolverError> {
    const CALLER: &str = "RobinProjector";
    if created {
        ctx.info(CALLER, &format!("Updating Robin flux constraint matrix for boundary: {}", i2s(bc as i64)), 6);
    } else {
        ctx.info(CALLER, &format!("Creating Robin flux constraint matrix for boundary: {}", i2s(bc as i64)), 6);
    }
    let slot = &ctx.sim.slots[ctx.slot];
    let sys = slot
        .system
        .as_ref()
        .ok_or_else(|| SolverError::fatal(CALLER, "no assembled system"))?;
    let bulk = slot
        .asm
        .bulk
        .as_ref()
        .map(|(v, _)| v)
        .ok_or_else(|| SolverError::fatal(CALLER, "BulkValues are needed to aveluate Robin terms!"))?;
    let id = slot.variable.expect("checked by the caller");
    let v = ctx.sim.vars.by_id(id);
    let dofs = v.dofs;
    let a = &sys.matrix;
    let nb = ctx.mesh.n_bulk() as u32;
    let range = if bc > nbc { 0..nb } else { nb..nb + ctx.mesh.n_boundary() as u32 };
    let mut active = vec![false; a.nrows];
    for e in range.map(ElemId) {
        if !belongs(ctx, e, bc, nbc) {
            continue;
        }
        let slots: Option<Vec<usize>> = ctx.mesh.element(e).nodes.iter().map(|&k| v.slot(k)).collect();
        let Some(slots) = slots else { continue };
        for s in slots {
            for i in 0..dofs {
                active[dofs * s + i] = true;
            }
        }
    }
    let mut rows = vec![BTreeMap::new(); dofs];
    for j in 0..a.nrows {
        if !active[j] {
            continue;
        }
        for k in a.row_ptr[j]..a.row_ptr[j + 1] {
            let col = a.col_idx[k];
            if !active[col as usize] {
                continue;
            }
            *rows[j % dofs].entry(col).or_insert(0.0) += -(a.values[k] - bulk[k]);
        }
    }
    Ok(rows)
}

/// Whether element `e` belongs to boundary condition `bc` or, past the
/// boundary conditions, to body force `bc - nbc`.
fn belongs(ctx: &SolverCtx<'_>, e: ElemId, bc: usize, nbc: usize) -> bool {
    let tag = ctx.mesh.element(e).tag;
    let model = &ctx.sim.model;
    if bc > nbc {
        model
            .body_of_tag(tag)
            .and_then(|b| model.bodies[b - 1].integer("body force"))
            .is_some_and(|k| k >= 1 && k as usize == bc - nbc)
    } else {
        model.bc_of_tag(tag) == Some(bc)
    }
}

/// `IntegralProjector` for boundary condition `bc` (or body force
/// `bc - nbc`).
fn integral_projector(ctx: &SolverCtx<'_>, bc: usize, nbc: usize) -> Result<Projector, SolverError> {
    const CALLER: &str = "IntegralProjector";
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature for integral BCs (elmer-fe)".into()))?;
    if bc > nbc {
        ctx.info(CALLER, &format!("Creating integral constraint matrix for body force: {}", i2s((bc - nbc) as i64)), 6);
    } else {
        ctx.info(CALLER, &format!("Creating integral constraint matrix for boundary: {}", i2s(bc as i64)), 6);
    }
    let axi = ctx
        .sim
        .model
        .simulation
        .string("coordinate system")
        .is_some_and(|c| matches!(elmer_core::normalize_key(c).as_str(), "axi symmetric" | "cylindric symmetric"));
    let nb = ctx.mesh.n_bulk() as u32;
    let range = if bc > nbc { 0..nb } else { nb..nb + ctx.mesh.n_boundary() as u32 };
    let mut entries = BTreeMap::new();
    let mut ev = elmer_core::BasisEval::default();
    for e in range.clone().map(ElemId) {
        if !belongs(ctx, e, bc, nbc) {
            continue;
        }
        let el = ctx.mesh.element(e);
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
        for qp in fe.quadrature.gauss_points(el.etype, None) {
            fe.basis
                .eval(el.etype, &xyz, qp.uvw, &mut ev)
                .map_err(|m| SolverError::fatal(CALLER, m))?;
            let mut w = ev.det_j * qp.weight;
            if axi {
                w *= ev.basis.iter().zip(&xyz).map(|(b, c)| b * c[0]).sum::<f64>();
            }
            for (p, &node) in el.nodes.iter().enumerate() {
                *entries.entry(node).or_insert(0.0) += w * ev.basis[p];
            }
        }
    }
    let area = entries.values().fold(0.0, |s, v| s + v);
    ctx.info(CALLER, &format!("Total area of boundary integral:{}", es(area, 12, 4)), 6);

    // SetInvPermIndex: the largest node of these elements that no other
    // boundary (or body force) touches.
    let mut other = vec![false; ctx.mesh.n_nodes()];
    for e in range.clone().map(ElemId) {
        if !belongs(ctx, e, bc, nbc) {
            for &k in ctx.mesh.element(e).nodes {
                other[k as usize] = true;
            }
        }
    }
    let mut node = None;
    for e in range.map(ElemId) {
        if belongs(ctx, e, bc, nbc) {
            for &k in ctx.mesh.element(e).nodes {
                if !other[k as usize] {
                    node = node.max(Some(k));
                }
            }
        }
    }
    let node = node.ok_or_else(|| SolverError::fatal(CALLER, "Could not determine maximum unset index!"))?;
    ctx.info(CALLER, &format!("Setting the representative node to: {}", i2s(node as i64 + 1)), 8);
    Ok(Projector { entries, node })
}

/// `GenerateProjectors` + `GenerateConstraintMatrix` for `Apply Integral
/// BCs`; refreshes the slot's generated constraint rows.
pub(crate) fn generate(ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
    const CALLER: &str = "GenerateConstraintMatrix";
    let ix = ctx.slot;
    ctx.sim.slots[ix].constraints = None;
    let p = ctx.params();
    if p.logical("apply integral bcs") != Some(true) {
        return Ok(());
    }
    for k in ["apply mortar bcs", "apply contact bcs"] {
        if p.logical(k) == Some(true) {
            return Err(SolverError::Unsupported(format!("{k} (mortar projectors)")));
        }
    }
    for k in ["mortar bc master solver", "mortar diag", "lumped diag"] {
        if p.get(k).is_some() {
            return Err(SolverError::Unsupported(format!("{k} with Apply Integral BCs")));
        }
    }
    if ctx.sim.slots[ix].own_constraints.is_some() {
        return Err(SolverError::Unsupported(
            "Apply Integral BCs together with a solver's own constraint matrix".into(),
        ));
    }
    let id = ctx.sim.slots[ix]
        .variable
        .ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    let (name, dofs, perm) = {
        let v = ctx.sim.vars.by_id(id);
        (v.name.clone(), v.dofs, v.perm.clone())
    };
    let perm = perm.ok_or_else(|| SolverError::fatal(CALLER, "solver variable has no permutation"))?;
    let skip_constrained = p.logical("skip already constrained dofs") == Some(true);
    let eps = p.real("minimum projector value").unwrap_or(1.0e-8);

    let model = &ctx.sim.model;
    let (nbc, nbf) = (model.boundary_conditions.len(), model.body_forces.len());
    let mut projectors: Vec<Option<Proj>> = Vec::with_capacity(nbc + nbf);
    for i in 1..=nbc + nbf {
        let s = if i > nbc {
            &ctx.sim.model.body_forces[i - nbc - 1]
        } else {
            &ctx.sim.model.boundary_conditions[i - 1]
        };
        let robin = s.logical("flux integral bc") == Some(true);
        let mortar = s.get("mortar bc").is_some() || s.get("contact bc").is_some();
        let integral = s.logical("integral bc") == Some(true);
        if robin && integral {
            return Err(SolverError::Unsupported("Integral BC and Flux Integral BC on one boundary".into()));
        }
        if mortar && integral {
            return Err(SolverError::fatal(
                "GenerateProjectors",
                format!("Boundary {} can only be one of mortar, contact and integral!", i2s(i as i64)),
            ));
        }
        if mortar {
            return Err(SolverError::Unsupported("Mortar BC / Contact BC projectors".into()));
        }
        projectors.push(if integral { Some(Proj::Integral(integral_projector(ctx, i, nbc)?)) } else { None });
    }
    // GenerateRobinProjectors (ProjectorUtils.F90:253), at DefaultSolve.
    let mut robins = false;
    for i in 1..=nbc + nbf {
        let s = if i > nbc {
            &ctx.sim.model.body_forces[i - nbc - 1]
        } else {
            &ctx.sim.model.boundary_conditions[i - 1]
        };
        if s.logical("flux integral bc") != Some(true) {
            continue;
        }
        ctx.info("GenerateRobinProjector", &format!("Generating flux integral conditions for BC: {}", i2s(i as i64)), 4);
        let created = ctx.sim.slots[ix].robin_projectors.contains(&i);
        projectors[i - 1] = Some(Proj::Robin(robin_projector(ctx, i, nbc, created)?));
        ctx.sim.slots[ix].robin_projectors.insert(i);
        robins = true;
    }
    if robins {
        // A multiplier of a physical meaning (ListAddNewLogical).
        let p = &mut ctx.sim.slots[ix].params;
        if p.get("export lagrange multiplier").is_none() {
            p.insert("Export Lagrange Multiplier", elmer_core::Value::Logical(true));
        }
    }
    let nrows: usize = projectors
        .iter()
        .flatten()
        .map(|p| match p {
            Proj::Integral(_) => 1,
            Proj::Robin(rows) => rows.len(),
        })
        .sum();
    if nrows == 0 {
        ctx.info(CALLER, "Nothing to do since there are no constrained dofs!", 12);
        return Ok(());
    }
    // ProjectorUtils.F90:638-650 (Apply Contact BCs is refused above).
    let dim = ctx.mesh.dim();
    if !(dofs == 1 || dofs == dim || dofs == dim + 1) {
        return Err(SolverError::fatal(CALLER, format!("Invalid number of dofs for field: {}", i2s(dofs as i64))));
    }
    if dofs != 1 && projectors.iter().flatten().any(|p| matches!(p, Proj::Robin(_))) {
        return Err(SolverError::Unsupported(format!("Flux Integral BC for a field of {dofs} components")));
    }
    let components: Vec<String> = if dofs == 1 { vec![name.clone()] } else { ctx.sim.vars.component_names(id) };
    if components.len() != dofs {
        return Err(SolverError::fatal(CALLER, format!("{name} has no names for its {dofs} components")));
    }
    let n = ctx.sim.vars.by_id(id).values.len();
    ctx.info(CALLER, &format!("There are {} initial rows in constraint matrices", i2s(nrows as i64)), 10);

    let constrained: std::collections::HashSet<u32> = ctx.sim.slots[ix].asm.constrained.iter().map(|&(d, _)| d).collect();
    let mut c = ConstraintMatrix {
        matrix: CsrMatrix {
            nrows: 0,
            ncols: ctx.sim.vars.by_id(id).values.len(),
            row_ptr: vec![0],
            col_idx: Vec::new(),
            values: Vec::new(),
            diag: None,
        },
        rhs: Vec::new(),
        inv_perm: Vec::new(),
    };
    for bc in (1..=nbc + nbf).rev() {
        let Some(proj) = &projectors[bc - 1] else { continue };
        let proj = match proj {
            Proj::Integral(p) => p,
            Proj::Robin(rows) => {
                // AddScalarConstraint for a Robin projector: the row as it
                // is, the diagonal entry -rsum at column `arows + j`
                // (MortarBC % Diag = 1), no right-hand side, no InvPerm.
                ctx.info(CALLER, &format!("Adding flux constraint for BC: {}", i2s(bc as i64)), 8);
                ctx.info(CALLER, &format!("Adding projector rows: {}", i2s(rows.len() as i64)), 12);
                for row in rows {
                    let valsum = row.values().fold(0.0, |s, v: &f64| s + v.abs());
                    let mut rsum = 0.0;
                    for (&col, &v) in row {
                        if v.abs() < eps * valsum {
                            continue;
                        }
                        rsum += v;
                        c.matrix.col_idx.push(col);
                        c.matrix.values.push(v);
                    }
                    c.matrix.col_idx.push(n as u32);
                    c.matrix.values.push(-rsum);
                    c.matrix.row_ptr.push(c.matrix.values.len());
                    c.matrix.nrows += 1;
                    c.rhs.push(0.0);
                    c.inv_perm.push(NO_DOF);
                }
                continue;
            }
        };
        let s = if bc > nbc {
            &ctx.sim.model.body_forces[bc - nbc - 1]
        } else {
            &ctx.sim.model.boundary_conditions[bc - 1]
        };
        // ListGetCReal(BC, 'Integral BC '//component): a component
        // without a value is not constrained (ActiveComponents; SetDefined
        // takes the stale Found of 'Integral BC', so none is "unspecified").
        let mut values = Vec::with_capacity(dofs);
        for comp in &components {
            values.push(ctx.creal(s, &format!("integral bc {comp}"))?);
        }
        if dofs == 1 && values[0].is_none() {
            ctx.info(CALLER, &format!("Skipping component: {}", i2s(1)), 12);
            continue;
        }
        ctx.info(CALLER, &format!("Adding integral constraint for BC: {}", i2s(bc as i64)), 8);
        ctx.info(CALLER, &format!("Adding projector rows: {}", i2s(1)), 12);
        let kk = perm.get(proj.node as usize).copied().unwrap_or(NO_DOF);
        if kk == NO_DOF {
            continue;
        }
        let valsum = proj.entries.values().fold(0.0, |s, v: &f64| s + v.abs());
        // AddScalarConstraint, or AddVectorConstraint's row per component.
        for (j, value) in values.iter().enumerate() {
            let Some(value) = *value else {
                ctx.info(CALLER, &format!("Skipping component: {}", i2s(j as i64 + 1)), 12);
                continue;
            };
            let row_dof = kk * dofs as u32 + j as u32;
            if skip_constrained && constrained.contains(&row_dof) {
                continue;
            }
            for (&node, &v) in &proj.entries {
                if v.abs() < eps * valsum {
                    continue;
                }
                let col = perm.get(node as usize).copied().unwrap_or(NO_DOF);
                if col == NO_DOF {
                    continue;
                }
                c.matrix.col_idx.push(col * dofs as u32 + j as u32);
                c.matrix.values.push(v);
            }
            c.matrix.row_ptr.push(c.matrix.values.len());
            c.matrix.nrows += 1;
            c.rhs.push(value);
            c.inv_perm.push(row_dof);
        }
    }
    let nnz = c.matrix.values.len();
    if nnz == 0 {
        ctx.info(CALLER, "No entries in constraint matrix!", 8);
        return Ok(());
    }
    let (rows, nz) = (i2s(c.matrix.nrows as i64), i2s(nnz as i64));
    ctx.info(CALLER, &format!("Allocating {rows} rows and {nz} nonzeros"), 8);
    ctx.info(CALLER, &format!("Used {rows} rows and {nz} nonzeros"), 7);
    ctx.sim.slots[ix].constraints = Some(c);
    Ok(())
}

/// Options of `SolveWithLinearRestriction` and of the solve around it that
/// are not ported.
pub(crate) fn refuse_options(p: &ParamMap) -> Result<(), SolverError> {
    for k in [
        "no explicit constrained matrix",
        "eliminate linear constraints",
        "constraint scaling",
        "use tree gauge",
        "restricted system norm",
        "restriction system residual mode",
        "linear system residual mode",
        "linear system complex",
        "linear system block mode",
        "nonlinear system convergence without constraints",
        "nonlinear system consistent norm",
        "nonlinear system aitken relaxation",
    ] {
        if p.logical(k) == Some(true) {
            return Err(SolverError::Unsupported(format!("{k} with linear constraints")));
        }
    }
    if let Some(m) = p.string("nonlinear system convergence measure") {
        if m != "norm" {
            return Err(SolverError::Unsupported(format!(
                "Nonlinear System Convergence Measure = {m} with linear constraints"
            )));
        }
    }
    Ok(())
}

/// `SolveWithLinearRestriction`'s collection system `[A C^T; C 0]` and its
/// right-hand side `[b; d]`, built in the Fortran's order. With `Enforce
/// Exact Dirichlet BCs` (default) a Dirichlet dof's column is moved to the
/// right-hand side, and a row whose representative dof has a Dirichlet
/// condition becomes `lambda = 0`.
pub(crate) fn collection(
    a: &CsrMatrix<f64>,
    b: &[f64],
    c: &ConstraintMatrix,
    fixed: &[(u32, f64)],
    p: &ParamMap,
    log: &dyn Logger,
) -> (CsrMatrix<f64>, Vec<f64>) {
    const CALLER: &str = "SolveWithLinearRestriction";
    let (n, m) = (a.nrows, c.matrix.nrows);
    let enforce = p.logical("enforce exact dirichlet bcs").unwrap_or(true);
    let mut dval: Vec<Option<f64>> = vec![None; n];
    if enforce {
        for &(d, v) in fixed {
            if let Some(x) = dval.get_mut(d as usize) {
                *x = Some(v);
            }
        }
    }
    let mut rows: Vec<BTreeMap<u32, f64>> = vec![BTreeMap::new(); n + m];
    let mut rhs = vec![0.0; n + m];
    log.info(CALLER, "Adding ConstraintMatrix into CollectionMatrix", 10);
    let mut empty = 0;
    for i in (0..m).rev() {
        let k = n + i;
        rows[k].insert(k as u32, 0.0);
        let span = c.matrix.row_ptr[i]..c.matrix.row_ptr[i + 1];
        let rowsum = span.clone().rev().fold(0.0, |s, j| s + c.matrix.values[j].abs());
        let mut nonempty = false;
        if rowsum > f64::EPSILON {
            for j in span.rev() {
                let (col, v) = (c.matrix.col_idx[j] as usize, c.matrix.values[j]);
                if v.abs() < f64::EPSILON * rowsum {
                    continue;
                }
                match dval.get(col).copied().flatten() {
                    Some(d) => rhs[k] -= v * d,
                    None => {
                        *rows[col].entry(k as u32).or_insert(0.0) += v;
                        if col < n {
                            *rows[k].entry(col as u32).or_insert(0.0) += v;
                            nonempty |= v != 0.0;
                        }
                    }
                }
            }
        }
        // SolveCore.F90:5864-5867: InvPerm carries the projector's row
        // offset (a multiple of the matrix size); the dof is it modulo n.
        let rep = c.inv_perm.get(i).copied().unwrap_or(NO_DOF);
        let rep = if rep == NO_DOF || n == 0 { NO_DOF } else { (rep as usize % n) as u32 };
        if enforce && rep != NO_DOF && dval.get(rep as usize).copied().flatten().is_some() {
            rows[k].values_mut().for_each(|v| *v = 0.0);
            rows[k].insert(k as u32, 1.0);
            rhs[k] = 0.0;
        } else if !nonempty {
            empty += 1;
            rhs[k] = 0.0;
        } else {
            rhs[k] += c.rhs[i];
        }
    }
    if empty > 0 {
        log.info(
            CALLER,
            &format!("Constraint Matrix has {} empty rows out of {}", i2s(empty), i2s(m as i64)),
            6,
        );
    }
    for i in (0..n).rev() {
        for j in (a.row_ptr[i]..a.row_ptr[i + 1]).rev() {
            *rows[i].entry(a.col_idx[j]).or_insert(0.0) += a.values[j];
        }
        rhs[i] += b[i];
    }
    log.info(CALLER, "Reverting CollectionMatrix back to CRS matrix", 10);
    let mut out = CsrMatrix {
        nrows: n + m,
        ncols: n + m,
        row_ptr: vec![0],
        col_idx: Vec::new(),
        values: Vec::new(),
        diag: None,
    };
    let mut diag = Vec::with_capacity(n + m);
    for (i, r) in rows.into_iter().enumerate() {
        let mut d = None;
        for (col, v) in r {
            if col as usize == i {
                d = Some(out.values.len());
            }
            out.col_idx.push(col);
            out.values.push(v);
        }
        diag.push(d);
        out.row_ptr.push(out.values.len());
    }
    out.diag = diag.into_iter().collect();
    log.info(
        CALLER,
        &format!(
            "Collection matrix increased with {} rows and {} non-zeros",
            i2s(m as i64),
            i2s((out.values.len() - a.values.len()) as i64)
        ),
        8,
    );
    (out, rhs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn csr(rows: &[&[(u32, f64)]], ncols: usize) -> CsrMatrix<f64> {
        let mut m = CsrMatrix { nrows: rows.len(), ncols, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
        for r in rows {
            for &(c, v) in *r {
                m.col_idx.push(c);
                m.values.push(v);
            }
            m.row_ptr.push(m.values.len());
        }
        m
    }

    /// A row whose `InvPerm` is a Dirichlet dof plus the matrix size (the
    /// second projector's offset) becomes `lambda = 0`, as with the plain
    /// dof (SolveCore.F90:5864-5876).
    #[test]
    fn inv_perm_wraps_by_the_matrix_size() {
        let a = csr(&[&[(0, 1.0)], &[(1, 2.0), (2, -1.0)], &[(1, -1.0), (2, 2.0)]], 3);
        let b = [0.5, 1.0, 1.0];
        let fixed = [(0u32, 0.5)];
        let row = |inv: u32| ConstraintMatrix {
            matrix: csr(&[&[(0, 1.0), (1, -1.0)]], 4),
            rhs: vec![0.25],
            inv_perm: vec![inv],
        };
        let log = elmer_core::NullLogger;
        let p = ParamMap::new();
        for inv in [0, 3, 6] {
            let (m, rhs) = collection(&a, &b, &row(inv), &fixed, &p, &log);
            let k = 3;
            let entries: Vec<(u32, f64)> = (m.row_ptr[k]..m.row_ptr[k + 1]).map(|j| (m.col_idx[j], m.values[j])).collect();
            assert!(entries.iter().all(|&(c, v)| if c == 3 { v == 1.0 } else { v == 0.0 }), "inv {inv}: {entries:?}");
            assert_eq!(rhs[k], 0.0, "inv {inv}");
        }
        // A free representative keeps the row.
        let (m, rhs) = collection(&a, &b, &row(4), &fixed, &p, &log);
        let k = 3;
        assert!((m.row_ptr[k]..m.row_ptr[k + 1]).any(|j| m.col_idx[j] == 1 && m.values[j] == -1.0));
        assert_eq!(rhs[k], 0.25 - 0.5);
    }
}
