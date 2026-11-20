//! `UpdateDependentComponents` (LumpingUtils.F90:759): after a solver with
//! `Update Components`, each listed `Component` section's `Operator k` on
//! its `Variable k` is reduced over the component's `Master Bodies` (or
//! `Master Boundaries`), with `ComponentNodalReduction` and
//! `ComponentIntegralReduction`, and saved as `res: ...` keywords in the
//! component section and in the Simulation. Refused: the `force`, `moment`
//! and `torque` operators and elemental (DG) variables.

use elmer_core::{ElemId, FieldKind, Logger, MeshView, Params, Value};

use crate::error::SolverError;
use crate::expr::SimFields;
use crate::sim::Simulation;


const CALLER: &str = "UpdateDependentComponents";

/// The master entities of a component: `(ids, boundary mode)`.
fn masters(comp: &elmer_core::ParamMap) -> Option<(Vec<i32>, bool)> {
    if let Some(Value::Integer(v)) = comp.get("master bodies") {
        return Some((v.clone(), false));
    }
    if let Some(Value::Integer(v)) = comp.get("master boundaries") {
        return Some((v.clone(), true));
    }
    None
}

impl Simulation {
    /// Whether element `e` belongs to the component's master entities.
    fn in_component(&self, mesh: &dyn MeshView, e: ElemId, ids: &[i32], bc_mode: bool) -> bool {
        let tag = mesh.element(e).tag;
        if bc_mode {
            // Element % BoundaryInfo % Constraint: the boundary condition.
            self.model.bc_of_tag(tag).is_some_and(|b| ids.contains(&(b as i32)))
        } else {
            self.model.body_of_tag(tag).is_some_and(|b| ids.contains(&(b as i32)))
        }
    }

    fn elements(&self, mesh: &dyn MeshView, bc_mode: bool) -> std::ops::Range<u32> {
        let nb = mesh.n_bulk() as u32;
        if bc_mode {
            nb..nb + mesh.n_boundary() as u32
        } else {
            0..nb
        }
    }

    /// `ComponentNodalReduction`.
    fn component_nodal_reduction(&self, mesh: &dyn MeshView, comp: &elmer_core::ParamMap, var: &str, oper: &str) -> Result<f64, SolverError> {
        let Some((ids, bc_mode)) = masters(comp) else { return Ok(0.0) };
        let v = self.vars.get(var).ok_or_else(|| SolverError::fatal(CALLER, format!("no variable {var}")))?;
        if v.kind != FieldKind::Nodal {
            return Err(SolverError::Unsupported(format!("component reduction of the {:?} variable {var}", v.kind)));
        }
        let dofs = v.dofs;
        let mut visited = vec![false; mesh.n_nodes()];
        let (mut sumi, mut sumx, mut sumxx, mut sumabsx) = (0usize, 0.0f64, 0.0f64, 0.0f64);
        let (mut maximum, mut minimum, mut absmax, mut absmin) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut initialized = false;
        for t in self.elements(mesh, bc_mode) {
            let e = ElemId(t);
            if !self.in_component(mesh, e, &ids, bc_mode) {
                continue;
            }
            for &j in mesh.element(e).nodes {
                let Some(k) = v.slot(j) else { continue };
                if visited[j as usize] {
                    continue;
                }
                visited[j as usize] = true;
                let x = if dofs <= 1 {
                    v.values[k]
                } else {
                    let mut x = 0.0;
                    for l in 0..dofs {
                        x += v.values[dofs * k + l] * v.values[dofs * k + l];
                    }
                    x.sqrt()
                };
                if !initialized {
                    initialized = true;
                    (maximum, minimum, absmax, absmin) = (x, x, x, x);
                }
                sumi += 1;
                sumx += x;
                sumxx += x * x;
                sumabsx += x.abs();
                maximum = x.max(maximum);
                minimum = x.min(minimum);
                if x.abs() > absmax.abs() {
                    absmax = x;
                }
                if x.abs() < absmin.abs() {
                    absmin = x;
                }
            }
        }
        if sumi == 0 {
            return Ok(0.0);
        }
        let n = sumi as f64;
        Ok(match oper {
            "sum" => sumx,
            "sum abs" => sumabsx,
            "min" => minimum,
            "max" => maximum,
            "min abs" => absmin,
            "max abs" => absmax,
            "range" => maximum - minimum,
            "mean" => sumx / n,
            "mean abs" => sumabsx / n,
            "variance" => (sumxx / n - (sumx * sumx) / (n * n)).sqrt(),
            _ => 0.0,
        })
    }

    /// `ComponentIntegralReduction`.
    fn component_integral_reduction(
        &self,
        mesh: &dyn MeshView,
        comp: &elmer_core::ParamMap,
        var: &str,
        oper: &str,
        coeff: Option<&str>,
    ) -> Result<f64, SolverError> {
        let Some((ids, bc_mode)) = masters(comp) else { return Ok(0.0) };
        let v = self.vars.get(var).ok_or_else(|| SolverError::fatal(CALLER, format!("no variable {var}")))?;
        if v.kind != FieldKind::Nodal {
            return Err(SolverError::Unsupported(format!("component integral of the {:?} variable {var}", v.kind)));
        }
        let fe = self
            .backends
            .fe
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("basis and quadrature for component integrals (elmer-fe)".into()))?;
        let fields = SimFields { vars: &self.vars, mesh, dg: self.dg_indexes.as_ref() };
        let dofs = v.dofs;
        let (mut integral, mut vol) = (0.0f64, 0.0f64);
        let mut ev = elmer_core::BasisEval::default();
        for t in self.elements(mesh, bc_mode) {
            let e = ElemId(t);
            if !self.in_component(mesh, e, &ids, bc_mode) {
                continue;
            }
            let el = mesh.element(e);
            let n = el.nodes.len();
            let Some(perm) = el.nodes.iter().map(|&j| v.slot(j)).collect::<Option<Vec<usize>>>() else { continue };
            let c: Option<Vec<f64>> = match coeff {
                Some(name) => {
                    let mat = self
                        .model
                        .body_of_tag(el.tag)
                        .and_then(|b| self.model.bodies.get(b - 1))
                        .and_then(|b| b.integer("material"))
                        .and_then(|m| self.model.materials.get((m as usize).wrapping_sub(1)));
                    Some(match mat.and_then(|m| m.get(name)) {
                        Some(val) => crate::expr::eval_at_nodes_with(val, el.nodes, &fields, self.backends.matc.as_deref(), &self.backends.udfs)?,
                        None => vec![0.0; n],
                    })
                }
                None => None,
            };
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
            for qp in fe.quadrature.gauss_points(el.etype, None) {
                fe.basis.eval(el.etype, &xyz, qp.uvw, &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let s = ev.det_j * qp.weight;
                let coef = match &c {
                    Some(c) => {
                        let mut x = 0.0;
                        for i in 0..n {
                            x += c[i] * ev.basis[i];
                        }
                        x
                    }
                    None => 1.0,
                };
                vol += s;
                let interp = |k: usize| -> f64 {
                    let mut f = 0.0;
                    for i in 0..n {
                        f += v.values[dofs * perm[i] + k] * ev.basis[i];
                    }
                    f
                };
                match oper {
                    "volume" => integral += coef * s,
                    "int" | "int mean" => integral += s * coef * interp(0),
                    "int abs" | "int abs mean" => integral += s * coef * interp(0).abs(),
                    "diffusive energy" => {
                        let mut g = [0.0f64; 3];
                        for (j, gj) in g.iter_mut().enumerate() {
                            let mut x = 0.0;
                            for i in 0..n {
                                x += ev.dbasis_dx[i][j] * v.values[perm[i]];
                            }
                            *gj = x;
                        }
                        integral += s * coef * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
                    }
                    "convective energy" => {
                        for k in 0..dofs {
                            let f = interp(k);
                            integral += s * coef * (f * f);
                        }
                    }
                    "potential energy" => integral += s * coef * interp(0),
                    _ => {}
                }
            }
        }
        Ok(match oper {
            "int mean" | "int abs mean" => integral / vol,
            "diffusive energy" | "convective energy" => 0.5 * integral,
            _ => integral,
        })
    }

    /// `UpdateDependentComponents(ComponentList)`.
    pub(crate) fn update_dependent_components(&mut self, list: &[i32], mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        log.info("UpdateDepedentComponents", "Updating Components to reflect new solution", 6);
        for j in 1..=self.model.components.len() {
            if !list.contains(&(j as i32)) {
                continue;
            }
            log.info("UpdateDepedentComponents", &format!("Updating component: {j}"), 4);
            for k in 1.. {
                let comp = self.model.components[j - 1].clone();
                let oper = comp.string(&format!("operator {k}")).map(str::to_string);
                let mut var = comp.string(&format!("variable {k}")).map(str::to_string);
                let mut coeff = comp.string(&format!("coefficient {k}")).map(str::to_string);
                if var.is_none() && oper.as_deref() == Some("electric resistance") {
                    var = Some("Potential".into());
                }
                let (Some(oper), Some(var)) = (oper, var) else { break };
                let oper = oper.trim().to_ascii_lowercase();
                if self.vars.get(&var).is_none() {
                    log.info(CALLER, &format!("Variable not available: {var}"), 4);
                    continue;
                }
                let value = match oper.as_str() {
                    "electric resistance" => {
                        let coeff = coeff.get_or_insert_with(|| "electric conductivity".into()).clone();
                        let power = self.component_integral_reduction(mesh, &comp, &var, "diffusive energy", Some(&coeff))?;
                        let voltage = self.component_nodal_reduction(mesh, &comp, &var, "range")?;
                        let r = voltage * voltage / power;
                        self.model.components[j - 1].insert(&format!("res: {oper}"), Value::real(r));
                        r
                    }
                    "sum" | "sum abs" | "min" | "max" | "min abs" | "max abs" | "range" | "mean" | "mean abs" | "variance" => {
                        let x = self.component_nodal_reduction(mesh, &comp, &var, &oper)?;
                        self.model.components[j - 1].insert(&format!("res: {oper} {var}"), Value::real(x));
                        x
                    }
                    "volume" | "int" | "int abs" | "int mean" | "int abs mean" | "diffusive energy" | "convective energy" | "potential energy" => {
                        let x = self.component_integral_reduction(mesh, &comp, &var, &oper, coeff.as_deref())?;
                        self.model.components[j - 1].insert(&format!("res: {oper} {var}"), Value::real(x));
                        x
                    }
                    "force" | "moment" | "torque" => {
                        return Err(SolverError::Unsupported(format!("the \"{oper}\" component operator (ComponentNodalForceReduction)")));
                    }
                    _ => return Err(SolverError::fatal(CALLER, format!("Uknown operator: {oper}"))),
                };
                log.info(CALLER, &format!("comp {j}: {oper}: {var}: {}", crate::fortran_fmt::es(value, 15, 6)), 5);
                self.model.simulation.insert(&format!("res: comp {j}: {oper} {var}"), Value::real(value));
            }
        }
        Ok(())
    }
}
