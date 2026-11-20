//! The constraint matrix of a solver's mortar projectors
//! (`GenerateConstraintMatrix`, ProjectorUtils.F90).

use elmer_core::{Logger, Params, NO_DOF};

use crate::projector::Projector;
use crate::MortarError;

const CALLER: &str = "GenerateConstraintMatrix";

/// `A % ConstraintMatrix`: one row per constrained node of each projector,
/// columns in the solver's matrix numbering (0-based), in the order the
/// Fortran stores them (not sorted).
#[derive(Clone, Debug, PartialEq)]
pub struct ConstraintMatrix {
    /// Length `n_rows + 1`.
    pub rows: Vec<usize>,
    pub cols: Vec<u32>,
    pub values: Vec<f64>,
    pub rhs: Vec<f64>,
    /// `InvPerm`, 0-based: `rowoffset + slot`, where `slot` is the
    /// constrained node's matrix row and `rowoffset` grows by the matrix size
    /// for each projector added before this one.
    pub inv_perm: Vec<usize>,
}

impl ConstraintMatrix {
    pub fn n_rows(&self) -> usize {
        self.rows.len() - 1
    }
}

/// The solver-side inputs of [`generate_constraint_matrix`].
pub struct ConstraintTarget<'a> {
    /// The solver section.
    pub solver: &'a dyn Params,
    /// The variable name (`Solver % Variable % Name`), for the per-component
    /// keywords such as `Mortar BC Temperature`.
    pub name: &'a str,
    pub dofs: usize,
    /// `CoordinateSystemDimension()`: a field of `dim` (or, under `Apply
    /// Contact BCs`, `dim + 1`) components is a vector field.
    pub dim: usize,
    /// Parent node -> matrix row (`Solver % Variable % Perm`), [`NO_DOF`]
    /// when absent.
    pub perm: &'a [u32],
    /// `Solver % Matrix % NumberOfRows`.
    pub n_rows: usize,
}

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// `GenerateConstraintMatrix` for the solver's mortar projectors
/// (`projectors[i]` from boundary condition `i + 1`, as
/// [`crate::generate_mortar_projectors`] returns them). Ported: a scalar
/// field, a real vector field (`AddVectorConstraint`: one row per projector
/// row and active component, the contact's `MortarBC % Active` mask, `Mortar
/// BC <var> <i>` per component, not additive and without `MortarBC % Diag`), `Minimum Projector Value`, nodal (`[D - P]`) and Galerkin
/// projectors, `Mortar BC <var> = False` to skip a boundary, and `Mortar
/// BCs Additive`: one row per constrained node across all projectors, the
/// contributions of every projector summed (`CRS_PackMatrix`), with
/// `Projector Priority` deciding which boundary owns a shared node. Returns
/// `None` when there is nothing to constrain.
pub fn generate_constraint_matrix(
    bcs: &[&dyn Params],
    projectors: &[Option<Projector>],
    t: &ConstraintTarget<'_>,
    log: &dyn Logger,
) -> Result<Option<ConstraintMatrix>, MortarError> {
    let s = t.solver;
    let bcount = projectors.iter().flatten().count();
    log.info(CALLER, &format!("Number of mortar matrices: {bcount}"), 12);
    if projectors.iter().flatten().all(|p| p.n_rows() == 0) {
        log.info(
            CALLER,
            "Nothing to do since there are no constrained dofs!",
            12,
        );
        return Ok(None);
    }
    if s.get("mortar diag").is_some() {
        return Err(MortarError::Unsupported("Mortar Diag".into()));
    }
    if flag(s, "skip already constrained dofs") {
        return Err(MortarError::Unsupported(
            "Skip Already Constrained Dofs".into(),
        ));
    }
    let summed = match s.logical("mortar bcs additive") {
        Some(a) => a,
        None => {
            (bcount > 1 && flag(s, "eliminate linear constraints"))
                || flag(s, "apply conforming bcs")
        }
    };
    // ComplexMatrix and cdofs (ProjectorUtils.F90:616-647).
    let complex = s
        .logical("linear system complex")
        .unwrap_or(t.dofs == 2 * t.dim);
    let contact = flag(s, "apply contact bcs");
    let cdofs = if t.dofs == 1 || t.dofs == t.dim {
        t.dofs
    } else if t.dofs == t.dim + 1 {
        if contact {
            t.dim
        } else {
            t.dofs
        }
    } else {
        return Err(MortarError::Fatal(format!(
            "{CALLER}: Invalid number of dofs for field: {}",
            t.dofs
        )));
    };
    if t.dofs > 1 && complex {
        return Err(MortarError::Unsupported(format!(
            "mortar constraints on a complex {}-component field",
            t.dofs
        )));
    }
    if t.dofs > 1 && summed {
        return Err(MortarError::Unsupported(
            "Mortar BCs Additive on a vector field".into(),
        ));
    }
    let mut any_priority = bcs.iter().any(|bc| bc.get("projector priority").is_some());
    if any_priority && !summed {
        log.warn(CALLER, "Priority has effect only in additive mode!");
        any_priority = false;
    }
    let eps = s.real("minimum projector value").unwrap_or(1.0e-8);

    // The order the Fortran adds the projectors in: from the last boundary
    // condition down, or by descending `Projector Priority * NumberOfBCs -
    // bc` (`SortI` of BCPriority, read from the end).
    let mut order: Vec<usize> = (0..projectors.len()).collect();
    if any_priority {
        let n = bcs.len() as i64;
        let key = |i: usize| {
            -(i as i64 + 1) + bcs[i].integer("projector priority").unwrap_or(0) as i64 * n
        };
        let (mut by_bc, rest): (Vec<usize>, Vec<usize>) =
            order.into_iter().partition(|&i| i < bcs.len());
        by_bc.sort_by_key(|&i| key(i));
        order = by_bc;
        order.extend(rest);
    }

    let mut rows: Vec<Vec<(u32, f64)>> = Vec::new();
    let mut rhs: Vec<f64> = Vec::new();
    let mut inv_perm: Vec<usize> = Vec::new();
    let mut sum_perm: Vec<i64> = vec![0; t.n_rows];
    let (mut eliminated, mut neglected) = (0usize, 0usize);
    let mut rowoffset = 0usize;
    let mut prev_priority = -1i64;
    for &i in order.iter().rev() {
        let Some(p) = &projectors[i] else { continue };
        let bc = *bcs.get(i).ok_or_else(|| {
            MortarError::Unsupported(format!(
                "mortar projector of body force {}",
                i + 1 - bcs.len()
            ))
        })?;
        let priority = if any_priority {
            bc.integer("projector priority").unwrap_or(0) as i64
        } else {
            -1
        };
        if flag(bc, "integral bc") {
            return Err(MortarError::Unsupported("Integral BC".into()));
        }
        log.info(
            CALLER,
            &format!("Adding mortar projector for BC: {}", i + 1),
            8,
        );
        log.info(
            CALLER,
            &format!("Adding projector rows: {}", p.n_rows()),
            12,
        );
        // ActiveComponents: `Mortar BC <component>` (the variable name for
        // a scalar).
        let mut comp_active = vec![true; cdofs];
        let (mut some_set, mut some_skip, mut all_defined) = (false, false, true);
        for (c, a) in comp_active.iter_mut().enumerate() {
            let key = if cdofs > 1 {
                format!("mortar bc {} {}", t.name, c + 1)
            } else {
                format!("mortar bc {}", t.name)
            };
            match bc.logical(&key) {
                Some(v) => {
                    *a = v;
                    if v {
                        some_set = true;
                    } else {
                        some_skip = true;
                    }
                }
                None => all_defined = false,
            }
        }
        if some_set && !all_defined {
            if some_skip {
                return Err(MortarError::Fatal(format!(
                    "{CALLER}: Do not know what to do with all {cdofs} components"
                )));
            }
            for (c, a) in comp_active.iter_mut().enumerate() {
                let key = format!("mortar bc {} {}", t.name, c + 1);
                if bc.logical(&key).is_none() {
                    *a = false;
                }
            }
        }
        // A skipped boundary does not advance the row offset or the
        // previous priority either.
        if t.dofs == 1 && !comp_active[0] {
            continue;
        }
        if t.dofs > 1 {
            if p.diag.is_some() {
                return Err(MortarError::Unsupported(
                    "MortarBC Diag on a vector field".into(),
                ));
            }
            add_vector(p, &comp_active, cdofs, eps, rowoffset, t, &mut rows, &mut rhs, &mut inv_perm)?;
            rowoffset += t.n_rows;
            prev_priority = priority;
            continue;
        }
        // MasterScale and SlaveScale: 1 unless the projector is a
        // Lagrange periodic constraint.
        let create_self = !p.galerkin;
        if summed && create_self {
            return Err(MortarError::Fatal(format!(
                "{CALLER}: It is impossible to sum up nodal projectors!"
            )));
        }
        let slot = |i: usize| -> Option<usize> {
            if p.rows[i] == p.rows[i + 1] {
                return None;
            }
            let k = p.inv_perm[i];
            if k == NO_DOF {
                return None;
            }
            let kk = t.perm[k as usize];
            (kk != NO_DOF).then_some(kk as usize)
        };

        // MortarBC % Perm: the nodes of the projector's rows.
        let slave_nodes: std::collections::HashSet<u32> = if p.mortar_perm {
            p.inv_perm
                .iter()
                .copied()
                .filter(|&k| k != NO_DOF)
                .collect()
        } else {
            Default::default()
        };

        // Number the summed rows.
        if summed {
            for kk in (0..p.n_rows()).filter_map(slot) {
                if sum_perm[kk] == 0 {
                    rows.push(Vec::new());
                    rhs.push(0.0);
                    inv_perm.push(0);
                    sum_perm[kk] = rows.len() as i64;
                } else if priority != prev_priority && sum_perm[kk] < 0 {
                    neglected += 1;
                } else {
                    eliminated += 1;
                }
            }
        }

        for i in 0..p.n_rows() {
            let Some(kk) = slot(i) else { continue };
            let row = if summed {
                let row = sum_perm[kk];
                if any_priority {
                    if row < 0 {
                        continue;
                    }
                    if priority != prev_priority {
                        sum_perm[kk] = -row;
                    }
                }
                if row <= 0 {
                    return Err(MortarError::Fatal(format!(
                        "{CALLER}: Invalid row index: {row}"
                    )));
                }
                row as usize - 1
            } else {
                rows.push(Vec::new());
                rhs.push(0.0);
                inv_perm.push(0);
                rows.len() - 1
            };
            inv_perm[row] = rowoffset + kk;
            let mut valsum = 0.0;
            for v in &p.values[p.rows[i]..p.rows[i + 1]] {
                valsum += v.abs();
            }
            let mut wsum = 0.0;
            for (col, val) in p.row(i) {
                if val.abs() < eps * valsum {
                    continue;
                }
                let col2 = *t.perm.get(col as usize).ok_or_else(|| {
                    MortarError::Fatal(format!(
                        "{CALLER}: col index too large: {} vs {}",
                        col + 1,
                        t.perm.len()
                    ))
                })?;
                if col2 == NO_DOF {
                    continue;
                }
                // A nodal projector's target side takes MasterScale; with
                // MortarBC % Perm a Galerkin row's own nodes take
                // SlaveScale and only they count in the weight.
                let scale = if create_self {
                    wsum += val;
                    p.master_scale
                } else if p.mortar_perm {
                    if slave_nodes.contains(&col) {
                        wsum += val;
                        p.slave_scale
                    } else {
                        p.master_scale
                    }
                } else {
                    wsum += val;
                    1.0
                };
                rows[row].push((col2, scale * val));
            }
            if create_self {
                rows[row].push((kk as u32, p.slave_scale * wsum));
            }
            if let Some(d) = &p.diag {
                if summed {
                    return Err(MortarError::Unsupported(
                        "MortarBC Diag in additive mode".into(),
                    ));
                }
                rows[row].push(((t.n_rows + row) as u32, -d[i] * wsum));
            }
            if let Some(r) = &p.rhs {
                rhs[row] += wsum * r[i];
            }
        }
        if !summed {
            rowoffset += t.n_rows;
        }
        prev_priority = priority;
    }

    if rows.iter().all(Vec::is_empty) {
        log.info(CALLER, "No entries in constraint matrix!", 8);
        return Ok(None);
    }
    if summed {
        log.info(
            CALLER,
            &format!("Number of eliminated rows: {eliminated}"),
            6,
        );
        if eliminated > 0 {
            for r in &mut rows {
                pack_row(r);
            }
        }
    }
    if neglected > 0 {
        log.info(CALLER, &format!("Number of neglected rows: {neglected}"), 6);
    }

    let mut c = ConstraintMatrix {
        rows: vec![0],
        cols: Vec::new(),
        values: Vec::new(),
        rhs,
        inv_perm,
    };
    for r in rows {
        for (col, v) in r {
            c.cols.push(col);
            c.values.push(v);
        }
        c.rows.push(c.cols.len());
    }
    log.info(
        CALLER,
        &format!("Used {} rows and {} nonzeros", c.n_rows(), c.cols.len()),
        7,
    );
    Ok(Some(c))
}

/// `AddVectorConstraint` (ProjectorUtils.F90:1469-1753) of one projector,
/// not summed and without `MortarBC % Diag`: for each projector row `i` and
/// active component `j` (and active `MortarBC % Active(cdofs*i + j)`), a row
/// on the matrix dofs `dofs * Perm(col) + j`. `wsum` counts only the nodal
/// projector's own entries or, with `MortarBC % Perm`, the slave ones.
#[allow(clippy::too_many_arguments)]
fn add_vector(
    p: &Projector,
    comp_active: &[bool],
    cdofs: usize,
    eps: f64,
    rowoffset: usize,
    t: &ConstraintTarget<'_>,
    rows: &mut Vec<Vec<(u32, f64)>>,
    rhs: &mut Vec<f64>,
    inv_perm: &mut Vec<usize>,
) -> Result<(), MortarError> {
    let create_self = !p.galerkin;
    let slave_nodes: std::collections::HashSet<u32> = if p.mortar_perm {
        p.inv_perm
            .iter()
            .copied()
            .filter(|&k| k != NO_DOF)
            .collect()
    } else {
        Default::default()
    };
    let dofs = t.dofs;
    for i in 0..p.n_rows() {
        if p.rows[i] >= p.rows[i + 1] {
            continue;
        }
        for (j, &on) in comp_active.iter().enumerate() {
            if !on {
                continue;
            }
            if let Some(a) = &p.active {
                if !a[cdofs * i + j] {
                    continue;
                }
            }
            let k = p.inv_perm[i];
            if k == NO_DOF {
                continue;
            }
            let kk = t.perm[k as usize];
            if kk == NO_DOF {
                continue;
            }
            let row = rows.len();
            rows.push(Vec::new());
            rhs.push(0.0);
            inv_perm.push(rowoffset + dofs * kk as usize + j);
            let mut valsum = 0.0;
            for v in &p.values[p.rows[i]..p.rows[i + 1]] {
                valsum += v.abs();
            }
            let mut wsum = 0.0;
            for (col, val) in p.row(i) {
                if val.abs() < eps * valsum {
                    continue;
                }
                // A column past the permutation is skipped (with a PRINT).
                let Some(&col2) = t.perm.get(col as usize) else { continue };
                if col2 == NO_DOF {
                    continue;
                }
                let scale = if create_self {
                    wsum += val;
                    p.master_scale
                } else if p.mortar_perm {
                    if slave_nodes.contains(&col) {
                        wsum += val;
                        p.slave_scale
                    } else {
                        p.master_scale
                    }
                } else {
                    1.0
                };
                rows[row].push(((dofs * col2 as usize + j) as u32, scale * val));
            }
            if create_self {
                rows[row].push(((dofs * kk as usize + j) as u32, p.slave_scale * wsum));
            }
            if let Some(r) = &p.rhs {
                rhs[row] = wsum * r[cdofs * i + j];
            }
        }
    }
    Ok(())
}

/// `CRS_PackMatrix` on one row: repeated columns merged into their first
/// occurrence, summed in order.
fn pack_row(r: &mut Vec<(u32, f64)>) {
    let mut out: Vec<(u32, f64)> = Vec::with_capacity(r.len());
    for &(c, v) in r.iter() {
        match out.iter_mut().find(|e| e.0 == c) {
            Some(e) => e.1 += v,
            None => out.push((c, v)),
        }
    }
    *r = out;
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::{NullLogger, ParamMap, Value};

    /// The contact case of `AddVectorConstraint`: `MortarBC % Active` picks
    /// the rows, `MortarBC % Perm` gives the slave columns `SlaveScale` and
    /// the weight, the master ones `MasterScale`, and the rhs is by
    /// component.
    #[test]
    fn vector_rows_follow_the_active_mask() {
        // Two slave nodes (0, 1) against two master nodes (2, 3).
        let p = Projector {
            rows: vec![0, 3, 6],
            cols: vec![0, 1, 2, 0, 1, 3],
            values: vec![0.5, 0.25, -0.75, 0.25, 0.5, -0.75],
            inv_perm: vec![0, 1],
            projector_bc: 1,
            galerkin: true,
            slave_scale: 1.0,
            master_scale: -1.0,
            mortar_perm: true,
            rhs: Some(vec![1.0, 2.0, 3.0, 4.0]),
            diag: None,
            active: Some(vec![true, false, false, true]),
        };
        let mut solver = ParamMap::new();
        solver.insert("Apply Contact BCs", Value::Logical(true));
        let bc = ParamMap::new();
        let bcs: Vec<&dyn Params> = vec![&bc];
        let perm: Vec<u32> = vec![0, 1, 2, 3];
        let t = ConstraintTarget { solver: &solver, name: "displacement", dofs: 2, dim: 2, perm: &perm, n_rows: 8 };
        let c = generate_constraint_matrix(&bcs, &[Some(p)], &t, &NullLogger).unwrap().unwrap();
        assert_eq!(c.n_rows(), 2);
        // Row 0: node 0, x; row 1: node 1, y.
        assert_eq!(c.inv_perm, vec![0, 3]);
        assert_eq!(c.cols, vec![0, 2, 4, 1, 3, 7]);
        assert_eq!(c.values, vec![0.5, 0.25, 0.75, 0.25, 0.5, 0.75]);
        // wsum is the slave weight: 0.75 per row.
        assert_eq!(c.rhs, vec![0.75 * 1.0, 0.75 * 4.0]);
    }
}
