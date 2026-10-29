//! Field selection: which `Scalar Field n` / `Vector Field n` names become
//! point-data and cell-data arrays, and how their values are read.
//!
//! This follows the `VariableGet` lookups of `WriteVtuFile`
//! (VtuOutputSolver.F90:827-1010 for nodal data, :1280-1450 for elemental
//! data) and the shared checks of `WritePvtuFile`.

use std::io;

use elmer_core::{ElemId, ElementRef, FieldKind, FieldRef, NO_DOF};

use crate::VtuOptions;

pub(crate) fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, msg.into())
}

pub(crate) fn unsupported(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::Unsupported, msg.into())
}

/// A selected field: either one variable with `dofs` components, or 1-3
/// single-component variables `name 1`, `name 2`, `name 3` ("component
/// vector"). As in Fortran, the first variable's `Perm` serves them all.
pub(crate) struct Selected<'a> {
    pub name: String,
    pub comps: Vec<FieldRef<'a>>,
    pub packed: bool,
    pub dofs: usize,
    /// Components written: `max(dofs, 3)` for vectors, 1 for scalars.
    pub sdofs: usize,
}

pub(crate) struct Plan<'a> {
    pub point: Vec<Selected<'a>>,
    pub cell: Vec<Selected<'a>>,
    pub warnings: Vec<String>,
}

fn find<'a, 'b>(fields: &'b [FieldRef<'a>], name: &str) -> Option<&'b FieldRef<'a>> {
    let name = name.trim();
    fields.iter().find(|f| f.name.trim().eq_ignore_ascii_case(name))
}

/// `Name` as one variable, else `Name 1..3` as components.
fn resolve<'a>(
    fields: &[FieldRef<'a>],
    name: &str,
    vector_rank: bool,
    is_cell: bool,
) -> Option<Selected<'a>> {
    let mut sol = find(fields, name).copied();
    // Elemental pass: asking for a vector, a one-component variable will not do.
    if is_cell && vector_rank {
        if let Some(f) = &sol {
            if f.dofs <= 1 {
                sol = None;
            }
        }
    }
    let (comps, packed, dofs) = match sol {
        Some(f) => (vec![f], true, f.dofs),
        None => {
            let f1 = *find(fields, &format!("{name} 1"))?;
            let mut comps = vec![f1];
            if let Some(f2) = find(fields, &format!("{name} 2")) {
                comps.push(*f2);
                if let Some(f3) = find(fields, &format!("{name} 3")) {
                    comps.push(*f3);
                }
            }
            let n = comps.len();
            (comps, false, n)
        }
    };
    let sdofs = if dofs > 1 { dofs.max(3) } else { 1 };
    Some(Selected { name: name.to_string(), comps, packed, dofs, sdofs })
}

/// Select point-data and cell-data arrays for `opts`.
pub(crate) fn plan<'a>(fields: &[FieldRef<'a>], opts: &VtuOptions) -> io::Result<Plan<'a>> {
    let mut p = Plan { point: Vec::new(), cell: Vec::new(), warnings: Vec::new() };
    let order = opts
        .scalar_fields
        .iter()
        .map(|n| (n, false))
        .chain(opts.vector_fields.iter().map(|n| (n, true)));

    for (name, vector_rank) in order.clone() {
        if opts.save_nodal {
            match resolve(fields, name, vector_rank, false) {
                None => p.warnings.push(format!("Nonexistent variable: {name}")),
                Some(s) => match s.comps[0].kind {
                    FieldKind::Nodal => p.point.push(s),
                    // Elemental, DG and integration-point variables are
                    // cell data when DG saving is off; nothing to do here.
                    FieldKind::Elemental | FieldKind::Dg | FieldKind::Ip => {}
                    FieldKind::Global | FieldKind::Edge => p
                        .warnings
                        .push(format!("Variable {name} has no point data layout; skipped")),
                },
            }
        }
    }
    if opts.save_elemental {
        for (name, vector_rank) in order {
            let Some(s) = resolve(fields, name, vector_rank, true) else {
                // Already reported by the point pass when that ran.
                if !opts.save_nodal {
                    p.warnings.push(format!("Nonexistent variable: {name}"));
                }
                continue;
            };
            match s.comps[0].kind {
                FieldKind::Elemental | FieldKind::Ip => p.cell.push(s),
                FieldKind::Dg => {
                    return Err(unsupported(format!(
                        "DG field {name}: elmer-core::MeshView has no per-element DG indexes"
                    )))
                }
                FieldKind::Nodal | FieldKind::Global | FieldKind::Edge => {}
            }
        }
    }
    Ok(p)
}

/// Slot of `node` in `f`'s value array, `None` where the variable is absent.
pub(crate) fn node_slot(f: &FieldRef<'_>, node: usize) -> io::Result<Option<usize>> {
    match f.perm {
        None => Ok(Some(node)),
        Some(p) => match p.get(node) {
            None => Err(invalid(format!(
                "Too big index for Perm of {} at node {}",
                f.name, node
            ))),
            Some(&NO_DOF) => Ok(None),
            Some(&s) => Ok(Some(s as usize)),
        },
    }
}

fn get(vals: &[f64], idx: usize, what: &str) -> io::Result<f64> {
    vals.get(idx)
        .copied()
        .ok_or_else(|| invalid(format!("value index {idx} out of range for {what}")))
}

impl Selected<'_> {
    /// Component `k` (0-based, `< dofs`) of value slot `slot`.
    pub(crate) fn at(&self, slot: usize, k: usize) -> io::Result<f64> {
        if self.packed {
            get(self.comps[0].values, self.dofs * slot + k, &self.name)
        } else {
            get(self.comps[k].values, slot, &self.name)
        }
    }

    /// The `sdofs` values of this field at `node`, zero where the variable is
    /// absent or the component does not exist.
    pub(crate) fn nodal(&self, node: usize, out: &mut [f64]) -> io::Result<()> {
        out.fill(0.0);
        if let Some(slot) = node_slot(&self.comps[0], node)? {
            for (k, o) in out.iter_mut().enumerate().take(self.dofs) {
                *o = self.at(slot, k)?;
            }
        }
        Ok(())
    }

    /// The `sdofs` values of this field on element `e`: the elemental value,
    /// or the mean over integration points (VtuOutputSolver.F90:1383-1475).
    ///
    /// Unlike Fortran, an element with zero integration points yields zeros
    /// (Fortran divides 0 by 0).
    pub(crate) fn elemental(
        &self,
        e: ElemId,
        el: &ElementRef<'_>,
        out: &mut [f64],
    ) -> io::Result<()> {
        out.fill(0.0);
        let f = &self.comps[0];
        let e = e.0 as usize;
        match f.kind {
            FieldKind::Elemental => {
                let slot = match f.perm {
                    None => Some(e),
                    Some(p) => {
                        let mut m = e;
                        if m >= p.len() {
                            // Boundary element: use the parent cell.
                            m = el.parents[0].or(el.parents[1]).map(|id| id.0 as usize).ok_or_else(
                                || invalid(format!("Cannot define parent cell index for element: {}", e + 1)),
                            )?;
                        }
                        match p.get(m) {
                            None => {
                                return Err(invalid(format!(
                                    "Perm of {} too short for element {}",
                                    f.name,
                                    m + 1
                                )))
                            }
                            Some(&NO_DOF) => None,
                            Some(&s) => Some(s as usize),
                        }
                    }
                };
                if let Some(slot) = slot {
                    for (k, o) in out.iter_mut().enumerate().take(self.dofs) {
                        *o = self.at(slot, k)?;
                    }
                }
            }
            FieldKind::Ip => {
                let p = f.perm.ok_or_else(|| {
                    invalid(format!("integration-point field {} has no Perm", f.name))
                })?;
                let (start, n) = match (p.get(e), p.get(e + 1)) {
                    (Some(&a), Some(&b)) if b >= a => (a as usize, (b - a) as usize),
                    _ => (0, 0),
                };
                if n > 0 {
                    for (k, o) in out.iter_mut().enumerate().take(self.dofs) {
                        let mut sum = 0.0;
                        for j in 0..n {
                            sum += self.at(start + j, k)?;
                        }
                        *o = sum / n as f64;
                    }
                }
            }
            _ => unreachable!("planned cell fields are elemental or integration-point"),
        }
        Ok(())
    }
}
