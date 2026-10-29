//! The lumped Dirichlet conditions of `SetDirichletBoundaries`
//! (BoundaryConditionUtils.F90:838-1300): `<var> Single Node` fixes one node
//! of a boundary or body force, `<var> Constant` and `<var> Profile` tie
//! every node of it to one supernode (`SetLumpedRows`), with the optional
//! `Resistance`, `Flux` and `Coefficient`.

use std::collections::BTreeMap;

use elmer_core::{CsrMatrix, ElemId, ParamMap, Params};

use crate::error::SolverError;
use crate::fortran_fmt::{es, i2s};
use crate::solver::SolverCtx;

const CALLER: &str = "SetDirichletBoundaries";

impl SolverCtx<'_> {
    /// The sections of the boundary conditions and then the body forces,
    /// with whether each is a body force.
    fn bc_and_bf_sections(&self) -> Vec<(&ParamMap, bool)> {
        let m = &self.sim.model;
        m.boundary_conditions
            .iter()
            .map(|s| (s, false))
            .chain(m.body_forces.iter().map(|s| (s, true)))
            .collect()
    }

    /// The elements of boundary condition `bc` (1-based) or, past the
    /// boundary conditions, of body force `bc - nbc`, in mesh order.
    fn elements_of(&self, bc: usize) -> Vec<ElemId> {
        let m = &self.sim.model;
        let nbc = m.boundary_conditions.len();
        let nb = self.mesh.n_bulk() as u32;
        if bc > nbc {
            (0..nb)
                .map(ElemId)
                .filter(|&e| {
                    m.body_of_tag(self.mesh.element(e).tag)
                        .and_then(|b| m.bodies[b - 1].integer("body force"))
                        .is_some_and(|k| k >= 1 && k as usize == bc - nbc)
                })
                .collect()
        } else {
            (nb..nb + self.mesh.n_boundary() as u32)
                .map(ElemId)
                .filter(|&e| m.bc_of_tag(self.mesh.element(e).tag) == Some(bc))
                .collect()
        }
    }

    /// `<name> Single Node`: for each boundary condition and body force
    /// giving it, the first node (in element order) with a dof, fixed to the
    /// value, for variable `id`. Returns `(node, value)` in section order.
    pub(crate) fn single_nodes(&self, id: crate::VarId, name: &str) -> Result<Vec<(u32, f64)>, SolverError> {
        let key = format!("{name} single node");
        let mut out = Vec::new();
        let sections = self.bc_and_bf_sections();
        if !sections.iter().any(|(s, _)| s.get(&key).is_some()) {
            return Ok(out);
        }
        let v = self.sim.vars.by_id(id);
        for (i, (s, _)) in sections.iter().enumerate() {
            let Some(val) = self.creal(s, &key)? else { continue };
            let mut node = None;
            'elems: for e in self.elements_of(i + 1) {
                for &k in self.mesh.element(e).nodes {
                    if v.slot(k).is_some() {
                        node = Some(k);
                        break 'elems;
                    }
                }
            }
            match node {
                None => self.warn(CALLER, &format!("Could not find node to set: {name} Single Node")),
                Some(k) => {
                    self.info(CALLER, &format!("Fixing single node {}", i2s(k as i64 + 1)), 8);
                    out.push((k, val));
                }
            }
        }
        Ok(out)
    }

    /// `<name> Constant` and `<name> Profile` on the assembled system,
    /// before the Dirichlet conditions are enforced: the row of each node
    /// of the boundary (body force) is added, times `Resistance`, to the
    /// supernode's row and replaced by `x_k - c x_0 = 0` (`c` = 1, or the
    /// profile over its value at the supernode). The matrix gains the new
    /// entries; the kept bulk system follows its topology.
    pub(crate) fn lumped_rows(&mut self) -> Result<(), SolverError> {
        let ix = self.slot;
        let Some(id) = self.sim.slots[ix].variable else { return Ok(()) };
        let (name, dofs) = {
            let v = self.sim.vars.by_id(id);
            (v.name.clone(), v.dofs)
        };
        let sections = self.bc_and_bf_sections();
        let lumped = |s: &ParamMap, profile: bool| {
            if profile {
                s.get(&format!("{name} profile")).is_some()
            } else {
                s.logical(&format!("{name} constant")) == Some(true)
            }
        };
        if !sections.iter().any(|(s, _)| lumped(s, false) || lumped(s, true)) {
            return Ok(());
        }
        if dofs != 1 {
            return Err(SolverError::Unsupported(format!(
                "{name} Constant / Profile for a field of {dofs} components"
            )));
        }
        let slot_of = |k: u32| self.sim.vars.by_id(id).slot(k);
        let Some(sys) = self.sim.slots[ix].system.as_ref() else {
            return Err(SolverError::fatal("DefaultDirichletBCs", "call default_initialize first"));
        };
        let n = sys.matrix.nrows;
        let mut rows: Vec<BTreeMap<u32, f64>> = (0..n)
            .map(|i| {
                (sys.matrix.row_ptr[i]..sys.matrix.row_ptr[i + 1])
                    .map(|j| (sys.matrix.col_idx[j], sys.matrix.values[j]))
                    .collect()
            })
            .collect();
        let mut b = sys.rhs.clone();
        let nbc = self.sim.model.boundary_conditions.len();
        let mut supernodes: Vec<(usize, usize)> = Vec::new();
        let mut last: Option<&ParamMap> = None;
        for profile in [false, true] {
            let dir = format!("{name} {}", if profile { "Profile" } else { "Constant" });
            let key = dir.to_lowercase();
            let (any_bc, any_bf) = (
                sections.iter().any(|(s, bf)| !bf && lumped(s, profile)),
                sections.iter().any(|(s, bf)| *bf && lumped(s, profile)),
            );
            if any_bc {
                self.info(CALLER, &format!("Found BC constraint: {dir}"), 6);
            }
            if any_bf {
                self.info(CALLER, &format!("Found BodyForce constraint: {dir}"), 6);
            }
            for (i, (s, _)) in sections.iter().enumerate() {
                if !lumped(s, profile) {
                    continue;
                }
                last = Some(s);
                let bc = i + 1;
                let move_coeff = self.creal(s, &format!("{key} resistance"))?.unwrap_or(1.0);
                let elems = self.elements_of(bc);
                let mult = |e: ElemId| -> Result<Vec<f64>, SolverError> {
                    self.real_at_nodes(s, &key, e)?
                        .ok_or_else(|| SolverError::fatal(CALLER, format!("Could not find {key}")))
                };
                // The supernode: the first node with a dof, or the one where
                // the profile is largest in magnitude.
                let (mut ind, mut max_mult) = (None, 0.0f64);
                'elems: for &e in &elems {
                    let m = if profile { mult(e)? } else { Vec::new() };
                    for (p, &k) in self.mesh.element(e).nodes.iter().enumerate() {
                        if slot_of(k).is_none() {
                            continue;
                        }
                        if profile {
                            if m[p].abs() > max_mult.abs() {
                                ind = Some(k);
                                max_mult = m[p];
                            }
                        } else {
                            ind = Some(k);
                            break 'elems;
                        }
                    }
                }
                self.info(CALLER, &format!("{dir} Max Value{}", es(max_mult, 12, 3)), 15);
                let Some(ind) = ind else {
                    self.warn(CALLER, &format!("Could not find node to set: {dir}"));
                    continue;
                };
                let k0 = slot_of(ind).expect("chosen with a dof");
                if (move_coeff - 1.0).abs() > f64::EPSILON {
                    rows[k0].values_mut().for_each(|v| *v *= move_coeff);
                    b[k0] *= move_coeff;
                }
                let mut set = vec![false; self.mesh.n_nodes()];
                set[ind as usize] = true;
                supernodes.push((bc, k0));
                for &e in &elems {
                    let m = if profile { mult(e)? } else { Vec::new() };
                    for (p, &node) in self.mesh.element(e).nodes.iter().enumerate() {
                        if set[node as usize] {
                            continue;
                        }
                        set[node as usize] = true;
                        let Some(k) = slot_of(node) else { continue };
                        let coeff = if profile { m[p] / max_mult } else { 1.0 };
                        // MoveRow(A, k, k0, MoveCoeff): row k0 += c row k, row
                        // k kept as structure with zero values.
                        let moved: Vec<(u32, f64)> = rows[k].iter().map(|(&c, &v)| (c, v)).collect();
                        for (c, v) in moved {
                            *rows[k0].entry(c).or_insert(0.0) += move_coeff * v;
                        }
                        rows[k].values_mut().for_each(|v| *v = 0.0);
                        b[k0] += move_coeff * b[k];
                        *rows[k].entry(k as u32).or_insert(0.0) += 1.0;
                        *rows[k].entry(k0 as u32).or_insert(0.0) += -coeff;
                        b[k] = 0.0;
                    }
                }
                if let Some(f) = self.creal(s, &format!("{key} flux"))? {
                    b[k0] += f;
                }
                if let Some(c) = self.creal(s, &format!("{key} coefficient"))? {
                    *rows[k0].entry(k0 as u32).or_insert(0.0) += c;
                }
                let count = set.iter().filter(|&&x| x).count();
                self.info(CALLER, &format!("Number of lumped nodes set: {}", i2s(count as i64)), 10);
            }
        }

        // The previous flux through each supernode (Elmer reads the
        // coefficient from the section it processed last).
        let x = &self.sim.vars.by_id(id).values;
        let coef = match last {
            Some(s) => self.creal(s, &format!("{} constant coefficient", name.to_lowercase()))?.unwrap_or(0.0),
            None => 0.0,
        };
        for &(bc, k0) in &supernodes {
            let mut flux = rows[k0].iter().fold(0.0, |s, (&c, &v)| s + v * x[c as usize]);
            flux -= coef * x[k0];
            let what = if bc <= nbc {
                format!("Previous bc {} lumped flux: ", i2s(bc as i64))
            } else {
                format!("Previous bf {} lumped flux: ", i2s((bc - nbc) as i64))
            };
            self.info(CALLER, &format!("{what}{}", es(flux, 12, 3)), 4);
        }

        let slot = &mut self.sim.slots[ix];
        let sys = slot.system.as_mut().expect("checked above");
        let old = std::mem::replace(&mut sys.matrix, to_csr(rows));
        sys.rhs = b;
        // CRS_ChangeTopology: the kept bulk values follow the new structure.
        if let Some((bulk, _)) = slot.asm.bulk.as_mut() {
            let mut moved = vec![0.0; sys.matrix.values.len()];
            for i in 0..old.nrows {
                for j in old.row_ptr[i]..old.row_ptr[i + 1] {
                    let c = old.col_idx[j];
                    let span = sys.matrix.row_ptr[i]..sys.matrix.row_ptr[i + 1];
                    let at = span.start + sys.matrix.col_idx[span].binary_search(&c).expect("structure only grows");
                    moved[at] = bulk[j];
                }
            }
            *bulk = moved;
        }
        self.info(
            "DefaultDirichletBCs",
            &format!("Modified matrix non-zeros: {}", i2s(self.sim.slots[ix].system.as_ref().map_or(0, |s| s.matrix.values.len()) as i64)),
            8,
        );
        Ok(())
    }
}

/// Rows with sorted columns as a CSR matrix with its diagonal index.
fn to_csr(rows: Vec<BTreeMap<u32, f64>>) -> CsrMatrix<f64> {
    let n = rows.len();
    let mut a = CsrMatrix {
        nrows: n,
        ncols: n,
        row_ptr: vec![0],
        col_idx: Vec::new(),
        values: Vec::new(),
        diag: None,
    };
    let mut diag = Vec::with_capacity(n);
    for (i, r) in rows.into_iter().enumerate() {
        let mut d = None;
        for (c, v) in r {
            if c as usize == i {
                d = Some(a.values.len());
            }
            a.col_idx.push(c);
            a.values.push(v);
        }
        diag.push(d);
        a.row_ptr.push(a.values.len());
    }
    a.diag = diag.into_iter().collect();
    a
}
