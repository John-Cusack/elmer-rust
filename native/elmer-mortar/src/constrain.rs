//! Using a periodic projector in a linear system: the matrix graph
//! (`ProjectorTopology`, ElementUtils.F90) and the constraint rows
//! (`SetPeriodicBoundariesPass1/2`, BoundaryConditionUtils.F90).
//!
//! The solver calls [`add_projector_topology`] while it builds the node
//! graph, before converting it to CSR, and [`set_periodic_boundaries`] once
//! per component at the start of its Dirichlet step, as
//! `SetDirichletBoundaries` does before any Dirichlet value is set.

use elmer_assembly::list::ListMatrix;
use elmer_core::{CsrMatrix, Params, NO_DOF};

use crate::projector::Projector;
use crate::MortarError;

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// Projectors whose constraint is not assembled into the matrix (explicit or
/// Lagrange-multiplier periodic conditions) add no couplings.
fn implicit(bc: &dyn Params) -> bool {
    !flag(bc, "periodic bc explicit") && !flag(bc, "periodic bc use lagrange coefficient")
}

/// `ProjectorTopology` for every boundary projector (the `DoProjectors`
/// loop of `MakeListMatrix`). `reorder` maps parent nodes to graph rows
/// ([`NO_DOF`] when absent); `bcs[i]` and `projectors[i]` belong to boundary
/// condition `i + 1`.
pub fn add_projector_topology(
    list: &mut ListMatrix,
    reorder: &[u32],
    bcs: &[&dyn Params],
    projectors: &[Option<Projector>],
) {
    for (bc, p) in bcs.iter().zip(projectors) {
        let Some(p) = p else { continue };
        if implicit(*bc) {
            projector_topology(list, reorder, p);
        }
    }
}

fn projector_topology(list: &mut ListMatrix, reorder: &[u32], p: &Projector) {
    let mut ms: Vec<u32> = Vec::new();
    let mut ss: Vec<u32> = Vec::new();
    for i in 0..p.n_rows() {
        let node = p.inv_perm[i];
        if node == NO_DOF || reorder[node as usize] == NO_DOF {
            continue;
        }
        let k = reorder[node as usize];
        // ProjectorRowSlaves
        ms.clear();
        ms.extend(
            p.row(i)
                .map(|(c, _)| reorder[c as usize])
                .filter(|&m| m != NO_DOF),
        );
        if ms.is_empty() {
            continue;
        }
        ms.sort_unstable();
        list.add_sorted(k as usize, &ms);
        // GatherRowWithSelf: row k with k itself merged in, unless no entry
        // of the row reaches k.
        ss.clear();
        let mut got_k = false;
        for &c in list.row(k as usize) {
            if !got_k && c >= k {
                if c > k {
                    ss.push(k);
                }
                got_k = true;
            }
            ss.push(c);
        }
        for &m in &ms {
            list.add_sorted(m as usize, &ss);
        }
        for &s in &ss {
            list.add_sorted(s as usize, &ms);
        }
    }
}

/// The component a periodic condition is applied to.
#[derive(Clone, Copy, Debug)]
pub struct PeriodicTarget<'a> {
    /// The variable or component name of the keywords, e.g. `Temperature`
    /// (`Periodic BC Temperature`) or `Velocity 1`.
    pub name: &'a str,
    /// 0-based component within the node block.
    pub dof: usize,
    pub ndofs: usize,
    /// Parent node -> 0-based node slot, [`NO_DOF`] when absent.
    pub perm: &'a [u32],
}

/// The `ActivePart` test of `SetDirichletBoundaries`: does boundary
/// condition `bc` make `name` periodic?
pub fn periodic_active(bc: &dyn Params, name: &str) -> bool {
    flag(bc, &format!("periodic bc {name}"))
        || flag(bc, &format!("anti periodic bc {name}"))
        || bc.get(&format!("periodic bc scale {name}")).is_some()
}

/// The sign tying the slave to the master side, or `None` when the boundary
/// condition does not constrain `name`.
fn periodic_scale(bc: &dyn Params, name: &str) -> Option<f64> {
    if flag(bc, &format!("periodic bc {name}")) {
        Some(if flag(bc, "antisymmetric bc") {
            1.0
        } else {
            -1.0
        })
    } else if flag(bc, &format!("anti periodic bc {name}")) {
        Some(1.0)
    } else {
        bc.real(&format!("periodic bc scale {name}"))
    }
}

/// `AddToMatrixElement` on the CSR structure `(row_ptr, cols)` with
/// coefficients `values`: a missing entry is an error unless the value is
/// zero.
fn add(
    row_ptr: &[usize],
    cols: &[u32],
    values: &mut [f64],
    i: usize,
    j: u32,
    v: f64,
) -> Result<(), MortarError> {
    let r = row_ptr[i]..row_ptr[i + 1];
    match cols[r.clone()].binary_search(&j) {
        Ok(k) => {
            values[r.start + k] += v;
            Ok(())
        }
        Err(_) if v == 0.0 => Ok(()),
        Err(_) => Err(MortarError::MissingEntry { row: i, col: j }),
    }
}

/// `SetPeriodicBoundariesPass1` and `Pass2` over all boundary conditions
/// for one component: first the slave rows are added to their master rows,
/// then each slave row becomes `P u_master - u_slave = -offset` (with the
/// sign from `Antisymmetric BC` / `Anti Periodic BC` / `Periodic BC Scale`).
///
/// `bulk` is the optional `(BulkValues, BulkRHS)` copy that Pass 1 updates
/// alongside `a` and `b` (it shares `a`'s structure). `n_nodes` sizes the
/// per-node `Done` mask. Returns the number of slave nodes set. Eigen and
/// harmonic mass/damping rows are not touched.
pub fn set_periodic_boundaries(
    a: &mut CsrMatrix<f64>,
    b: &mut [f64],
    mut bulk: (Option<&mut [f64]>, Option<&mut [f64]>),
    bcs: &[&dyn Params],
    projectors: &[Option<Projector>],
    t: &PeriodicTarget<'_>,
    n_nodes: usize,
) -> Result<usize, MortarError> {
    let active: Vec<bool> = bcs.iter().map(|bc| periodic_active(*bc, t.name)).collect();
    if !active.iter().any(|&x| x) {
        return Ok(0);
    }
    for (bc, &on) in bcs.iter().zip(&active) {
        if !on {
            continue;
        }
        for k in [
            "periodic bc explicit",
            "periodic bc use lagrange coefficient",
        ] {
            if flag(*bc, k) {
                return Err(MortarError::Unsupported(k.into()));
            }
        }
        if bc
            .get(&format!("periodic bc coefficient {}", t.name))
            .is_some()
        {
            return Err(MortarError::Unsupported(
                "Periodic BC Coefficient (jump)".into(),
            ));
        }
    }
    let row = |node: u32| -> Option<usize> {
        let s = t.perm[node as usize];
        (s != NO_DOF).then(|| t.ndofs * s as usize + t.dof)
    };

    // Pass 1: move the slave equations onto the master rows.
    let mut done = vec![false; n_nodes];
    for ((bc, p), _) in bcs
        .iter()
        .zip(projectors)
        .zip(&active)
        .filter(|(_, &on)| on)
    {
        let (Some(scale), Some(p)) = (periodic_scale(*bc, t.name), p) else {
            continue;
        };
        for i in 0..p.n_rows() {
            let ii = p.inv_perm[i];
            if ii == NO_DOF {
                continue;
            }
            if let (false, Some(k)) = (done[ii as usize], row(ii)) {
                for (c, w) in p.row(i) {
                    if w == 0.0 {
                        continue;
                    }
                    let Some(m) = row(c) else { continue };
                    for nn in a.row_ptr[k]..a.row_ptr[k + 1] {
                        let (col, v) = (a.col_idx[nn], a.values[nn]);
                        add(
                            &a.row_ptr,
                            &a.col_idx,
                            &mut a.values,
                            m,
                            col,
                            -scale * w * v,
                        )?;
                        if let Some(f) = bulk.0.as_deref_mut() {
                            let fv = f[nn];
                            add(&a.row_ptr, &a.col_idx, f, m, col, -scale * w * fv)?;
                        }
                    }
                    b[m] -= scale * w * b[k];
                    if let Some(f) = bulk.1.as_deref_mut() {
                        f[m] -= scale * w * f[k];
                    }
                }
            }
            done[ii as usize] = true;
        }
    }

    // Pass 2: the constraint rows.
    done.fill(false);
    for ((bc, p), _) in bcs
        .iter()
        .zip(projectors)
        .zip(&active)
        .filter(|(_, &on)| on)
    {
        let (Some(scale), Some(p)) = (periodic_scale(*bc, t.name), p) else {
            continue;
        };
        let offset = bc
            .real(&format!("periodic bc offset {}", t.name))
            .unwrap_or(0.0);
        for i in 0..p.n_rows() {
            let ii = p.inv_perm[i];
            if ii == NO_DOF {
                continue;
            }
            if let (false, Some(k)) = (done[ii as usize], row(ii)) {
                let r = a.row_ptr[k]..a.row_ptr[k + 1];
                a.values[r].fill(0.0);
                b[k] = 0.0;
                for (c, w) in p.row(i) {
                    let Some(m) = row(c) else { continue };
                    add(&a.row_ptr, &a.col_idx, &mut a.values, k, m as u32, w)?;
                }
                b[k] -= offset;
                add(&a.row_ptr, &a.col_idx, &mut a.values, k, k as u32, scale)?;
            }
            done[ii as usize] = true;
        }
    }
    Ok(done.iter().filter(|&&d| d).count())
}
