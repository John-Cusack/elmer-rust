//! Elemental (`-elem`) and integration-point (`-ip`) exported variables and
//! `UpdateExportedVariables` (SolverBasics.F90:6602): an exported variable
//! is filled from the body forces (and boundary conditions) that give a
//! keyword of its name, at the nodes of each element (nodal), as their mean
//! (elemental), or evaluated at each Gauss point with its dependencies
//! interpolated there (integration points, `ListGetElementReal`).

use elmer_core::{ElemId, FieldKind, Logger, MeshView, ParamMap, Params, Value, NO_DOF};
use elmer_expr::FieldLookup;

use crate::error::SolverError;
use crate::expr::SimFields;
use crate::fortran_fmt::i2s;
use crate::sim::{Simulation, SolverSlot};

const CALLER: &str = "UpdateExportedVariables";

/// `Exported Variable n Mask = "[section:] keyword"`: the section (default
/// `body force`) and the logical keyword that select the elements.
pub(crate) fn mask_of(params: &ParamMap, key: &str) -> Option<(String, String)> {
    let text = params.string(&format!("{key} mask"))?;
    let text = elmer_core::normalize_key(text);
    Some(match text.split_once(':') {
        Some((sec, mask)) if !sec.is_empty() => (sec.trim().to_string(), mask.trim().to_string()),
        _ => ("body force".to_string(), text),
    })
}

/// `ListGetSection(Element, name)`: the section of that kind applying to
/// element `e`, through its body (or its boundary condition).
fn section_of<'m>(sim: &'m Simulation, mesh: &dyn MeshView, e: ElemId, name: &str) -> Result<Option<&'m ParamMap>, SolverError> {
    let m = &sim.model;
    let tag = mesh.element(e).tag;
    let is_bulk = (e.0 as usize) < mesh.n_bulk();
    let body = if is_bulk { m.body_of_tag(tag).and_then(|b| m.bodies.get(b - 1)) } else { None };
    let linked = |sections: &'m [ParamMap], key: &str| {
        body.and_then(|b| b.integer(key)).filter(|&k| k >= 1).and_then(|k| sections.get(k as usize - 1))
    };
    Ok(match name {
        "body" => body,
        "material" => linked(&m.materials, "material"),
        "body force" => linked(&m.body_forces, "body force"),
        "initial condition" => linked(&m.initial_conditions, "initial condition"),
        "equation" => linked(&m.equations, "equation"),
        "boundary condition" => (!is_bulk).then(|| m.bc_of_tag(tag)).flatten().map(|i| &m.boundary_conditions[i - 1]),
        other => return Err(SolverError::fatal("ListGetSection", format!("Unknown section name: {other}"))),
    })
}

/// `ListGetLogicalGen(section, mask)`: false without the section or key.
fn mask_on(section: Option<&ParamMap>, mask: &str) -> bool {
    match section.and_then(|s| s.get(mask)) {
        Some(Value::Logical(b)) => *b,
        Some(Value::Real { data, .. }) => data.first().is_some_and(|&x| x > 0.0),
        Some(Value::Integer(v)) => v.first().is_some_and(|&x| x > 0),
        _ => false,
    }
}

/// `CreateElementsPerm` with a mask: active elements whose section sets the
/// mask, numbered in element order. As in the Fortran, an element whose
/// equation is not active (a boundary element, say) keeps the previous
/// element's verdict.
pub(crate) fn element_perm_masked(
    sim: &Simulation,
    mesh: &dyn MeshView,
    slot: &SolverSlot,
    (sec, mask): (&str, &str),
) -> Result<(Vec<u32>, usize), SolverError> {
    let n = mesh.n_bulk() + mesh.n_boundary();
    let mut active = vec![false; n];
    for e in &slot.active {
        active[e.0 as usize] = true;
    }
    let mut perm = vec![NO_DOF; n];
    let (mut on, mut m) = (false, 0usize);
    for t in 0..n {
        if active[t] {
            on = mask_on(section_of(sim, mesh, ElemId(t as u32), sec)?, mask);
        }
        if on {
            perm[t] = m as u32;
            m += 1;
        }
    }
    Ok((perm, m))
}

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
pub(crate) fn ip_perm(
    sim: &Simulation,
    mesh: &dyn MeshView,
    slot: &SolverSlot,
    mask: Option<(&str, &str)>,
) -> Result<Vec<u32>, SolverError> {
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
        let on = active[t]
            && match mask {
                Some((sec, key)) => mask_on(section_of(sim, mesh, ElemId(t as u32), sec)?, key),
                None => true,
            };
        if on {
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
    /// `InitCond` (ElmerSolver.F90:2062) for an elemental or IP exported
    /// variable: from the initial condition of each element's body, the
    /// element mean of the nodal values (scalar elemental fields), or the
    /// value at each Gauss point (IP fields, component by component).
    pub(crate) fn exported_initial_values(
        &self,
        mesh: &dyn MeshView,
        spec: &crate::sim::VarSpec,
        perm: &[u32],
        values: &mut [f64],
    ) -> Result<(), SolverError> {
        let (kind, name, components, dofs) = (spec.kind, spec.name.as_str(), &spec.components, spec.dofs);
        let m = &self.model;
        let names: Vec<(usize, &str)> = if dofs == 1 {
            vec![(0, name)]
        } else {
            components.iter().enumerate().map(|(j, c)| (j, c.as_str())).collect()
        };
        if dofs > 1 && m.initial_conditions.iter().any(|ic| ic.get(name).is_some()) {
            let what = if kind == FieldKind::Ip { "elemental" } else { "elements" };
            return Err(SolverError::fatal("InitCond", format!("Initialization only for scalar {what} fields!")));
        }
        if kind == FieldKind::Elemental && dofs > 1 {
            return Ok(());
        }
        let fields = self.fields(mesh);
        let b = &self.backends;
        for (j, key) in names {
            if !m.initial_conditions.iter().any(|ic| ic.get(key).is_some()) {
                continue;
            }
            for t in 0..mesh.n_bulk() {
                let e = ElemId(t as u32);
                let Some(ic) = section_of(self, mesh, e, "initial condition")? else { continue };
                let Some(value) = ic.get(key) else { continue };
                let el = mesh.element(e);
                match kind {
                    FieldKind::Elemental => {
                        let local = crate::expr::eval_at_nodes_with(value, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                        if let Some(&s) = perm.get(t).filter(|&&s| s != NO_DOF) {
                            values[s as usize] = local.iter().sum::<f64>() / local.len() as f64;
                        }
                    }
                    _ => {
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
                        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
                        let mut ev = elmer_core::BasisEval::default();
                        for (k, qp) in qps.iter().enumerate() {
                            fe.basis
                                .eval(el.etype, &xyz, qp.uvw, &mut ev)
                                .map_err(|msg| SolverError::fatal("InitCond", msg))?;
                            let at = AtPoint { inner: &fields, nodes: el.nodes, basis: &ev.basis };
                            let v = crate::expr::eval_at_nodes_with(value, &[0], &at, b.matc.as_deref(), &b.udfs)?;
                            values[(i1 + k) * dofs + j] = v[0];
                        }
                    }
                }
            }
        }
        Ok(())
    }

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
