//! Elemental (`-elem`) and integration-point (`-ip`) exported variables and
//! `UpdateExportedVariables` (SolverBasics.F90:6602): an exported variable
//! is filled from the body forces (and boundary conditions) that give a
//! keyword of its name, at the nodes of each element (nodal), as their mean
//! (elemental), or evaluated at each Gauss point with its dependencies
//! interpolated there (integration points, `ListGetElementReal`).

use elmer_core::{ElemId, FieldKind, Logger, MeshView, ParamMap, Params, NO_DOF};
use elmer_expr::FieldLookup;

use crate::error::SolverError;
use crate::expr::SimFields;
use crate::fortran_fmt::i2s;
use crate::sim::{Simulation, SolverSlot};

const CALLER: &str = "UpdateExportedVariables";

/// `InvActiveElements`: element -> 0-based position among the solver's
/// active elements, over every bulk and boundary element.
pub(crate) fn element_perm(mesh: &dyn MeshView, slot: &SolverSlot) -> (Vec<u32>, usize) {
    let mut perm = vec![NO_DOF; mesh.n_bulk() + mesh.n_boundary()];
    for (i, e) in slot.active.iter().enumerate() {
        perm[e.0 as usize] = i as u32;
    }
    (perm, slot.active.len())
}

/// `CreateIpPerm`: the offset of each element's Gauss points (the default
/// rule of an active element, or `Gauss Points on Ip Variables`), one entry
/// more than the bulk and boundary elements.
pub(crate) fn ip_perm(sim: &Simulation, mesh: &dyn MeshView, slot: &SolverSlot) -> Result<Vec<u32>, SolverError> {
    let fe = sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature for IP variables (elmer-fe)".into()))?;
    let n_ip = slot.params.integer("gauss points on ip variables").filter(|&n| n > 0);
    let n = mesh.n_bulk() + mesh.n_boundary();
    let mut active = vec![false; n];
    for e in &slot.active {
        active[e.0 as usize] = true;
    }
    let mut offsets = vec![0u32; n + 1];
    let mut count = 0u32;
    for t in 0..n {
        if active[t] {
            count += match n_ip {
                Some(k) => k as u32,
                None => fe.quadrature.gauss_points(mesh.element(ElemId(t as u32)).etype, None).len() as u32,
            };
        }
        offsets[t + 1] = count;
    }
    Ok(offsets)
}

/// The fields as seen at one point of an element: every variable and the
/// coordinates interpolated with the basis values there.
struct AtPoint<'a> {
    inner: &'a SimFields<'a>,
    nodes: &'a [u32],
    basis: &'a [f64],
}

impl FieldLookup for AtPoint<'_> {
    fn value(&self, var: &str, _node: u32) -> Option<f64> {
        let mut s = 0.0;
        for (&n, &b) in self.nodes.iter().zip(self.basis) {
            s += b * self.inner.value(var, n)?;
        }
        Some(s)
    }

    fn coord(&self, _node: u32) -> [f64; 3] {
        let mut x = [0.0; 3];
        for (&n, &b) in self.nodes.iter().zip(self.basis) {
            let c = self.inner.coord(n);
            for k in 0..3 {
                x[k] += b * c[k];
            }
        }
        x
    }

    fn time(&self) -> f64 {
        self.inner.time()
    }
}

impl Simulation {
    /// `UpdateExportedVariables` for slot `ix`.
    pub(crate) fn update_exported_variables(&mut self, ix: usize, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        log.info(CALLER, "Updating variables, if any!", 20);
        for l in 1.. {
            let Some(text) = self.slots[ix].params.string(&format!("exported variable {l}")).map(str::to_string) else {
                break;
            };
            let spec = crate::sim::parse_variable_spec(&text, None)?;
            let Some(id) = self.vars.id(&spec.name) else { continue };
            let (kind, dofs, len) = {
                let v = self.vars.by_id(id);
                (v.kind, v.dofs, v.values.len())
            };
            let names = self.vars.component_names(id);
            let comp = |j: usize| if dofs > 1 { names.get(j).cloned().unwrap_or_default() } else { spec.name.clone() };
            let m = &self.model;

            // A state variable: one value from the first body force giving it.
            if len == dofs || kind == FieldKind::Global {
                for j in 0..dofs {
                    let key = comp(j);
                    if let Some(bf) = m.body_forces.iter().find(|b| b.get(&key).is_some()) {
                        log.info(CALLER, "Found a proper definition for state variable", 6);
                        if let Some(v) = self.creal(bf, &key, mesh)? {
                            self.vars.by_id_mut(id).values[j] = v;
                        }
                    }
                }
                continue;
            }

            for j in 0..dofs {
                let key = comp(j);
                let cond = format!("{key} condition");
                let parts: Vec<bool> = m.body_forces.iter().map(|b| b.get(&key).is_some()).collect();
                let parts_bc: Vec<bool> = m.boundary_conditions.iter().map(|b| b.get(&key).is_some()).collect();
                let count = parts.iter().chain(&parts_bc).filter(|&&a| a).count();
                if count == 0 {
                    continue;
                }
                log.info(
                    CALLER,
                    &format!("Exported Variable {} defined in {} sections", i2s(l as i64), i2s(count as i64)),
                    8,
                );
                let mut writes: Vec<(usize, f64)> = Vec::new();
                let fields = self.fields(mesh);
                let b = &self.backends;
                let nb = mesh.n_bulk();
                for t in 0..nb + mesh.n_boundary() {
                    let e = ElemId(t as u32);
                    let el = mesh.element(e);
                    let mut section: Option<&ParamMap> = None;
                    if t < nb {
                        if let Some(bf) = m
                            .body_of_tag(el.tag)
                            .and_then(|bd| m.bodies[bd - 1].integer("body force"))
                            .filter(|&k| k >= 1 && parts.get(k as usize - 1) == Some(&true))
                        {
                            section = Some(&m.body_forces[bf as usize - 1]);
                        }
                    } else if let Some(bc) = m.bc_of_tag(el.tag).filter(|&i| parts_bc[i - 1]) {
                        section = Some(&m.boundary_conditions[bc - 1]);
                    }
                    let Some(section) = section else { continue };
                    let conditional = section.get(&cond).is_some();
                    let value = section.get(&key).expect("an active part");
                    let perm = self.vars.by_id(id).perm.clone();
                    let perm = perm.as_deref().unwrap_or(&[]);
                    match kind {
                        FieldKind::Ip => {
                            let (i1, i2) = (perm[t] as usize, perm[t + 1] as usize);
                            if i2 == i1 {
                                continue;
                            }
                            let fe = b.fe.as_ref().ok_or_else(|| {
                                SolverError::NotWired("basis and quadrature for IP variables (elmer-fe)".into())
                            })?;
                            let qps = fe.quadrature.gauss_points(el.etype, None);
                            if qps.len() != i2 - i1 {
                                return Err(SolverError::Unsupported(
                                    "IP variable whose Gauss point count changed (UpdateIpPerm)".into(),
                                ));
                            }
                            if conditional {
                                log.warn(CALLER, "Elemental variable cannot be conditional!");
                            }
                            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
                            let mut ev = elmer_core::BasisEval::default();
                            for (k, qp) in qps.iter().enumerate() {
                                fe.basis
                                    .eval(el.etype, &xyz, qp.uvw, &mut ev)
                                    .map_err(|msg| SolverError::fatal(CALLER, msg))?;
                                let at = AtPoint { inner: &fields, nodes: el.nodes, basis: &ev.basis };
                                let v = crate::expr::eval_at_nodes_with(value, &[0], &at, b.matc.as_deref(), &b.udfs)?;
                                writes.push(((i1 + k) * dofs + j, v[0]));
                            }
                        }
                        FieldKind::Elemental => {
                            if conditional {
                                log.warn(CALLER, "Elemental variables not conditional!");
                            }
                            let local = crate::expr::eval_at_nodes_with(value, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                            let s = perm.get(t).copied().unwrap_or(NO_DOF);
                            if s != NO_DOF {
                                let mean = local.iter().sum::<f64>() / local.len() as f64;
                                writes.push((s as usize * dofs + j, mean));
                            }
                        }
                        FieldKind::Nodal => {
                            let local = crate::expr::eval_at_nodes_with(value, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                            let slot = |n: u32| perm.get(n as usize).copied().filter(|&s| s != NO_DOF);
                            if conditional {
                                let c = crate::expr::eval_at_nodes_with(
                                    section.get(&cond).expect("checked"),
                                    el.nodes,
                                    &fields,
                                    b.matc.as_deref(),
                                    &b.udfs,
                                )?;
                                for (i, &n) in el.nodes.iter().enumerate() {
                                    if let (true, Some(s)) = (c[i] > 0.0, slot(n)) {
                                        writes.push((s as usize * dofs + j, local[i]));
                                    }
                                }
                            } else if el.nodes.iter().all(|&n| slot(n).is_some()) {
                                for (i, &n) in el.nodes.iter().enumerate() {
                                    writes.push((slot(n).expect("checked") as usize * dofs + j, local[i]));
                                }
                            }
                        }
                        k => return Err(SolverError::Unsupported(format!("updating a {k:?} exported variable"))),
                    }
                }
                let v = self.vars.by_id_mut(id);
                for (i, x) in writes {
                    v.values[i] = x;
                }
                v.values_changed = true;
            }
        }
        self.publish_live();
        Ok(())
    }
}
