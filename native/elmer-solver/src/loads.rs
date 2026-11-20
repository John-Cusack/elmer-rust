//! Nodal loads: `<var> Load` boundary and body-force values added to the
//! right-hand side (`SetNodalLoads`, SolverBasics.F90:1822), the bulk
//! system snapshot (`DefaultFinishBulkAssembly`, DefUtils.F90:7756) and the
//! `<var> Loads` / `<var> Residual` fields computed from it after the solve
//! (`CalculateLoads`, SolveCore.F90:364, called at :1739 and :1769).

use elmer_core::{CsrMatrix, ElemId, LinearOperator, Logger, MeshView, ParamMap, Params, Value};

use crate::error::SolverError;
use crate::sim::{Simulation, VarSpec};
use crate::variable::{VarId, Variable};

/// Per-assembly bookkeeping of one slot, reset by `DefaultInitialize`.
#[derive(Clone, Debug, Default)]
pub struct AssemblyState {
    /// `A % BulkValues`, `A % BulkRHS`.
    pub bulk: Option<(Vec<f64>, Vec<f64>)>,
    /// `DefaultFinishBulkAssembly` has run for this assembly.
    pub bulk_done: bool,
    /// `DefaultFinishBoundaryAssembly` has run for this assembly.
    pub boundary_done: bool,
    /// A boundary element added a nonzero contribution.
    pub boundary_touched: bool,
    /// The Dirichlet dofs of this assembly (`A % ConstrainedDOF`).
    pub constrained: Vec<(u32, f64)>,
    /// `A % PrecValues`: a separate matrix on the system's pattern that the
    /// preconditioner is built from (`DefaultUpdatePrec`). The Dirichlet
    /// conditions do not touch it.
    pub prec: Option<Vec<f64>>,
}

/// Whether `DefaultFinishBulkAssembly` keeps the bulk system: `Calculate
/// Loads` (unless `Calculate Loads Slot` names another stage), `Constant
/// Bulk System`, `Save Bulk System`, `Constant Bulk Matrix`, `Constraint
/// Modes Analysis` or `Control Use Loads` (DefUtils.F90:7797-7808).
pub(crate) fn keeps_bulk(p: &ParamMap) -> Result<bool, SolverError> {
    if [
        "control use loads",
        "constraint modes analysis",
        "constant bulk system",
        "save bulk system",
        "constant bulk matrix",
    ]
    .iter()
    .any(|k| p.logical(k) == Some(true))
    {
        return Ok(true);
    }
    if p.logical("calculate loads") != Some(true) {
        return Ok(false);
    }
    match p.string("calculate loads slot") {
        None => Ok(true),
        Some(s) if elmer_core::normalize_key(s) == "bulk assembly" => Ok(true),
        Some(s) => Err(SolverError::Unsupported(format!("Calculate Loads Slot = {s}"))),
    }
}

/// `A x - b` over the first `n` rows.
pub(crate) fn residual(a: &CsrMatrix<f64>, values: &[f64], x: &[f64], b: &[f64]) -> Vec<f64> {
    let m = CsrMatrix {
        nrows: a.nrows,
        ncols: a.ncols,
        row_ptr: a.row_ptr.clone(),
        col_idx: a.col_idx.clone(),
        values: values.to_vec(),
        diag: None,
    };
    let mut r = vec![0.0; a.nrows];
    m.apply(x, &mut r);
    for (ri, bi) in r.iter_mut().zip(b) {
        *ri -= bi;
    }
    r
}

impl Simulation {
    /// The `<var> Loads` and `<var> Residual` fields of a solver that asks
    /// for them (`Calculate Loads`, `Calculate Residual`; `Calculate
    /// Boundary Fluxes` implies loads), with the solver variable's perm
    /// (AddEquationBasics, MainUtils.F90:2342-2390).
    pub(crate) fn add_load_variables(&mut self, params: &mut ParamMap, spec: &VarSpec, primary: VarId) -> Result<(), SolverError> {
        // MainUtils.F90:2345: fluxes are summed from the nodal loads.
        if params.logical("calculate boundary fluxes") == Some(true) {
            if params.logical("calculate boundary weights") == Some(true) {
                return Err(SolverError::Unsupported("Calculate Boundary Weights".into()));
            }
            params.insert("Calculate Loads", elmer_core::Value::Logical(true));
        }
        for what in ["loads", "residual"] {
            if params.logical(&format!("calculate {what}")) != Some(true) {
                continue;
            }
            let name = format!("{} {what}", spec.name);
            if self.vars.id(&name).is_some() {
                continue;
            }
            let p = self.vars.by_id(primary);
            let mut v = Variable::new(&name, p.kind, p.dofs, p.perm.clone(), vec![0.0; p.values.len()]);
            v.output = params.logical(&format!("save {what}")).unwrap_or(p.output);
            self.vars.add(v)?;
        }
        Ok(())
    }

    /// `SetNodalLoads` for slot `ix`'s variable into its right-hand side:
    /// each node gets the first boundary condition's `<name> Load` that
    /// reaches it, then (separately) the first body force's.
    pub(crate) fn set_nodal_loads(
        &mut self,
        ix: usize,
        mesh: &dyn MeshView,
        log: &dyn Logger,
    ) -> Result<(), SolverError> {
        let Some(id) = self.slots[ix].variable else { return Ok(()) };
        let (name, dofs) = {
            let v = self.vars.by_id(id);
            (v.name.clone(), v.dofs)
        };
        let mut rhs_add: Vec<(usize, f64)> = Vec::new();
        // `ComponentName(name, DOF)`: the declared component names (`U 3` of
        // `Deflection[U:3 DNU:3]`) or `"<name> k"`.
        let names = self.vars.component_names(id);
        for dof in 0..dofs {
            let comp = if dofs == 1 {
                name.clone()
            } else {
                names.get(dof).cloned().unwrap_or_else(|| format!("{name} {}", dof + 1))
            };
            let key = format!("{comp} load");
            let point = |p: &ParamMap| p.get("target nodes").is_some() || p.get("target coordinates").is_some();
            for p in self.model.boundary_conditions.iter().chain(&self.model.body_forces) {
                if p.get(&format!("{key} dofs")).is_some() {
                    return Err(SolverError::Unsupported(format!("{key} DOFs")));
                }
                if p.get(&key).is_some() && point(p) && p.get(&format!("{key} im")).is_some() {
                    return Err(SolverError::Unsupported(format!("{key} im on target nodes")));
                }
            }
            let nslots = self.vars.by_id(id).values.len() / dofs;
            // Boundary conditions (one DoneLoad over all of them), then body forces.
            for bulk in [false, true] {
                let mut done = vec![false; nslots];
                let elements = if bulk {
                    0..mesh.n_bulk() as u32
                } else {
                    mesh.n_bulk() as u32..(mesh.n_bulk() + mesh.n_boundary()) as u32
                };
                let n_sections = if bulk { self.model.body_forces.len() } else { self.model.boundary_conditions.len() };
                for sec in 0..n_sections {
                    let section = if bulk { &self.model.body_forces[sec] } else { &self.model.boundary_conditions[sec] };
                    let Some(val) = section.get(&key) else { continue };
                    // Point loads come after the element loads.
                    if !bulk && point(section) {
                        continue;
                    }
                    if !bulk && section.get("target boundaries").is_none() && section.get("name").is_none() {
                        continue;
                    }
                    for e in elements.clone().map(ElemId) {
                        let el = mesh.element(e);
                        let owner = if bulk {
                            self.model
                                .body_of_tag(el.tag)
                                .and_then(|b| self.model.bodies[b - 1].integer("body force"))
                                .map(|b| b as usize)
                        } else {
                            self.model.bc_of_tag(el.tag)
                        };
                        if owner != Some(sec + 1) {
                            continue;
                        }
                        // `CurrentModel % CurrentElement => Element`.
                        elmer_expr::set_current_element(e.0, el.nodes);
                        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| mesh.coord(n)).collect();
                        elmer_expr::set_current_coords(&xyz);
                        let w = crate::expr::eval_at_nodes_with(
                            val,
                            el.nodes,
                            &self.fields(mesh),
                            self.backends.matc.as_deref(),
                            &self.backends.udfs,
                        )?;
                        let v = self.vars.by_id(id);
                        for (&n, &x) in el.nodes.iter().zip(&w) {
                            let Some(k) = v.slot(n) else { continue };
                            if done[k] {
                                continue;
                            }
                            done[k] = true;
                            rhs_add.push((dofs * k + dof, x));
                        }
                    }
                }
            }
            // SetPointLoads: the BCs' loads at their Target Nodes (or the
            // nodes closest to their Target Coordinates).
            for section in &self.model.boundary_conditions {
                let Some(val) = section.get(&key) else { continue };
                let Some(nodes) = crate::solver::target_nodes_of(mesh, section) else { continue };
                if nodes.is_empty() {
                    continue;
                }
                let w = crate::expr::eval_at_nodes_with(
                    val,
                    &nodes,
                    &self.fields(mesh),
                    self.backends.matc.as_deref(),
                    &self.backends.udfs,
                )?;
                let v = self.vars.by_id(id);
                for (&n, &x) in nodes.iter().zip(&w) {
                    if let Some(k) = v.slot(n) {
                        rhs_add.push((dofs * k + dof, x));
                    }
                }
            }
        }
        if let Some(sys) = self.slots[ix].system.as_mut() {
            for (k, x) in rhs_add {
                if let Some(b) = sys.rhs.get_mut(k) {
                    *b += x;
                }
            }
        }
        self.set_couple_loads(ix, mesh, log)
    }

    /// `SetCoupleLoads` (SolverBasics.F90:2127), which DefaultDirichletBCs
    /// calls for every solver: a `Torque 1..3` at a BC's `Target Nodes`
    /// with a `Director` becomes a force couple on the two matrix neighbours
    /// of the node that lie furthest along and against the director (the
    /// solid nodes across a shell-solid interface). Without a director the
    /// Fortran reads a `Torque` it never set; that is refused.
    fn set_couple_loads(&mut self, ix: usize, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        const CALLER: &str = "SetCoupleLoads";
        const AEPS: f64 = 10.0 * f64::EPSILON;
        let bcs = &self.model.boundary_conditions;
        if !bcs.iter().any(|b| b.0.keys().any(|k| k.starts_with("torque"))) {
            return Ok(());
        }
        let Some(id) = self.slots[ix].variable else { return Ok(()) };
        let dofs = self.vars.by_id(id).dofs;
        let Some(a) = self.slots[ix].system.as_ref().map(|s| &s.matrix) else { return Ok(()) };
        let v = self.vars.by_id(id);
        let mut inv = vec![u32::MAX; a.nrows / dofs.max(1)];
        for n in 0..mesh.n_nodes() as u32 {
            if let Some(k) = v.slot(n) {
                if k < inv.len() {
                    inv[k] = n;
                }
            }
        }
        let fields = self.fields(mesh);
        let eval = |bc: &ParamMap, key: &str, nodes: &[u32]| -> Result<Vec<f64>, SolverError> {
            match bc.get(key) {
                Some(val) => crate::expr::eval_at_nodes_with(
                    val,
                    nodes,
                    &fields,
                    self.backends.matc.as_deref(),
                    &self.backends.udfs,
                ),
                None => Ok(vec![0.0; nodes.len()]),
            }
        };
        let dot = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let e = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut add: Vec<(usize, f64)> = Vec::new();
        for bc in bcs {
            if !["torque 1", "torque 2", "torque 3"].iter().any(|k| bc.get(k).is_some()) {
                continue;
            }
            let Some(Value::Integer(ids)) = bc.get("target nodes") else {
                return Err(SolverError::fatal(CALLER, "Target Nodes not found"));
            };
            let nodes: Vec<u32> = ids.iter().map(|&i| (i - 1) as u32).collect();
            let work: Vec<Vec<f64>> =
                (1..=3).map(|k| eval(bc, &format!("torque {k}"), &nodes)).collect::<Result<_, _>>()?;
            if !["director 1", "director 2", "director 3"].iter().any(|k| bc.get(k).is_some()) {
                return Err(SolverError::Unsupported(
                    "Torque without a Director (SetCoupleLoads reads an unset torque there)".into(),
                ));
            }
            let dirs: Vec<Vec<f64>> =
                (1..=3).map(|k| eval(bc, &format!("director {k}"), &nodes)).collect::<Result<_, _>>()?;
            for (p, &target) in nodes.iter().enumerate() {
                let Some(ti) = v.slot(target) else { continue };
                // The last component's row; every dofs-th entry from its last.
                let row = ti * dofs + dofs - 1;
                let (lo, hi) = (a.row_ptr[row], a.row_ptr[row + 1]);
                let near: Vec<u32> =
                    (lo + dofs - 1..hi).step_by(dofs).map(|i| inv[a.col_idx[i] as usize / dofs]).collect();
                let x0 = mesh.coord(target);
                let rel = |n: u32| -> [f64; 3] {
                    let x = mesh.coord(n);
                    [x[0] - x0[0], x[1] - x0[1], x[2] - x0[2]]
                };
                let d = [dirs[0][p], dirs[1][p], dirs[2][p]];
                let nd = dot(&d, &d).sqrt();
                let e3 = [d[0] / nd, d[1] / nd, d[2] / nd];
                let w = [dot(&e[0], &e3).abs(), dot(&e[1], &e3).abs(), dot(&e[2], &e3).abs()];
                // MINLOC: the first smallest.
                let mut k = 0;
                for i in 1..3 {
                    if w[i] < w[k] {
                        k = i;
                    }
                }
                let vk = e[k];
                let c = dot(&vk, &e3);
                let mut e1 = [vk[0] - c * e3[0], vk[1] - c * e3[1], vk[2] - c * e3[2]];
                let n1 = dot(&e1, &e1).sqrt();
                e1 = [e1[0] / n1, e1[1] / n1, e1[2] / n1];
                let e2 = [
                    e3[1] * e1[2] - e3[2] * e1[1],
                    e3[2] * e1[0] - e3[0] * e1[2],
                    e3[0] * e1[1] - e3[1] * e1[0],
                ];
                let torque = [work[0][p], work[1][p], work[2][p]];
                let tc = dot(&torque, &e3);
                let t = [torque[0] - tc * e3[0], torque[1] - tc * e3[1], torque[2] - tc * e3[2]];
                let (m1, m2) = (dot(&t, &e1), dot(&t, &e2));
                let (mut maxres, mut minres) = (0.0f64, 0.0f64);
                let (mut jz, mut lz) = (None, None);
                for &m in &near {
                    if m == target {
                        continue;
                    }
                    let r = rel(m);
                    let res = dot(&e3, &r) / dot(&r, &r).sqrt();
                    if res > 0.0 {
                        if res > maxres {
                            jz = Some(m);
                            maxres = res;
                        }
                    } else if res < minres {
                        lz = Some(m);
                        minres = res;
                    }
                }
                let upper = jz.filter(|_| maxres.abs() >= AEPS);
                let lower = lz.filter(|_| minres.abs() >= AEPS);
                let (Some(jz), Some(lz)) = (upper, lower) else {
                    log.warn(CALLER, "A couple BC would need two nodes on opposite sides");
                    continue;
                };
                let rupper = dot(&e3, &rel(jz));
                let rlower = dot(&[-e3[0], -e3[1], -e3[2]], &rel(lz));
                let f2 = m1 / (rupper + rlower);
                let f1 = m2 / (rupper + rlower);
                let force = [f1 * e1[0] - f2 * e2[0], f1 * e1[1] - f2 * e2[1], f1 * e1[2] - f2 * e2[2]];
                let f = [dot(&force, &e[0]), dot(&force, &e[1]), dot(&force, &e[2])];
                let nc = if dofs > 2 { 3 } else { 2 };
                if let (Some(kj), Some(kl)) = (v.slot(jz), v.slot(lz)) {
                    for (c, &fc) in f.iter().enumerate().take(nc) {
                        add.push((kj * dofs + c, fc));
                    }
                    for (c, &fc) in f.iter().enumerate().take(nc) {
                        add.push((kl * dofs + c, -fc));
                    }
                }
            }
        }
        if let Some(sys) = self.slots[ix].system.as_mut() {
            for (k, x) in add {
                sys.rhs[k] += x;
            }
        }
        Ok(())
    }

    /// `CalculateLoads` into `<var> Loads` from the kept bulk system (after
    /// the solve, before the nonlinear change) or into `<var> Residual`
    /// from the solved system (after it), over the first `n` dofs.
    pub(crate) fn calculate_loads(
        &mut self,
        ix: usize,
        n: usize,
        from_bulk: bool,
        mesh: &dyn MeshView,
        log: &dyn Logger,
    ) -> Result<(), SolverError> {
        let slot = &self.slots[ix];
        let Some(id) = slot.variable else { return Ok(()) };
        let what = if from_bulk { "loads" } else { "residual" };
        let name = format!("{} {what}", self.vars.by_id(id).name);
        let Some(lid) = self.vars.id(&name) else { return Ok(()) };
        if slot.params.logical(&format!("calculate {what}")) == Some(false) {
            return Ok(());
        }
        for k in ["calculate energy inner product", "calculate energy norm"] {
            if slot.params.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(k.into()));
            }
        }
        let Some(sys) = slot.system.as_ref() else { return Ok(()) };
        let x = &self.vars.by_id(id).values[..n];
        let r = if from_bulk {
            let (values, rhs) = slot
                .asm
                .bulk
                .as_ref()
                .ok_or_else(|| SolverError::fatal("CalculateLoads", "\"BulkValues\" are not associated!"))?;
            residual(&sys.matrix, values, x, &rhs[..n])
        } else {
            residual(&sys.matrix, &sys.matrix.values, x, &sys.rhs[..n])
        };
        let loads = self.vars.by_id_mut(lid);
        let m = loads.values.len().min(r.len());
        loads.values[..m].copy_from_slice(&r[..m]);
        loads.values_changed = true;
        if from_bulk && self.slots[ix].params.logical("calculate boundary fluxes") == Some(true) {
            self.boundary_fluxes(id, lid, mesh, log)?;
        }
        Ok(())
    }

    /// The `Calculate Boundary Fluxes` part of `CalculateLoads`
    /// (SolveCore.F90:536-640): each boundary condition's sum of the nodal
    /// loads of its boundary nodes, a node counted for the first boundary
    /// that reaches it (`Calculate Boundary Flux` on some boundary limits the
    /// sums to those), stored as the Simulation entry `res: <var> Flux over
    /// BC i` (`<var> j Flux ...` per component) for SaveScalars.
    fn boundary_fluxes(&mut self, id: VarId, lid: VarId, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        let (name, dofs) = {
            let v = self.vars.by_id(id);
            (v.name.clone(), v.dofs)
        };
        log.info("CalculateLoads", "Computing boundary fluxes from nodal loads", 6);
        let nbc = self.model.boundary_conditions.len();
        let only_given = self.model.boundary_conditions.iter().any(|b| b.get("calculate boundary flux").is_some());
        let loads = self.vars.by_id(lid);
        let nslots = loads.values.len() / dofs.max(1);
        let mut sum = vec![0.0; nbc * dofs];
        let mut active = vec![false; nbc];
        let mut shared = vec![false; nbc];
        let mut summed = vec![0usize; nslots];
        let nb = mesh.n_bulk() as u32;
        for e in (nb..nb + mesh.n_boundary() as u32).map(ElemId) {
            let el = mesh.element(e);
            let Some(bc) = self.model.bc_of_tag(el.tag) else { continue };
            if only_given && self.model.boundary_conditions[bc - 1].logical("calculate boundary flux") != Some(true) {
                continue;
            }
            for &n in el.nodes {
                let Some(ind) = loads.slot(n) else { continue };
                if summed[ind] == 0 {
                    for j in 0..dofs {
                        sum[dofs * (bc - 1) + j] += loads.values[dofs * ind + j];
                    }
                    summed[ind] = bc;
                    active[bc - 1] = true;
                } else if summed[ind] != bc {
                    shared[bc - 1] = true;
                    shared[summed[ind] - 1] = true;
                }
            }
        }
        let mut n_active = 0;
        for i in 0..nbc {
            if !active[i] {
                continue;
            }
            if shared[i] {
                let msg = format!("Boundary {} includes inseparable dofs!", i + 1);
                if self.flux_shared_warned {
                    log.info("CalculateLoads", &msg, 8);
                } else {
                    log.warn("CalculateLoads", &msg);
                }
            }
            n_active += 1;
            for j in 0..dofs {
                let what = if dofs == 1 {
                    format!("{name} Flux over BC {}", i + 1)
                } else {
                    format!("{name} {} Flux over BC {}", j + 1, i + 1)
                };
                let v = sum[dofs * i + j];
                self.model.simulation.insert(&format!("res: {what}"), elmer_core::Value::real(v));
                log.info("CalculateLoads", &format!("{what}: {}", crate::fortran_fmt::es(v, 12, 5)), 6);
            }
        }
        self.flux_shared_warned = true;
        if n_active > 1 {
            for j in 0..dofs {
                let what = if dofs == 1 {
                    format!("{name} Flux over all BCs")
                } else {
                    format!("{name} {} Flux over all BCs", j + 1)
                };
                let total: f64 = sum.iter().skip(j).step_by(dofs).sum();
                log.info("CalculateLoads", &format!("{what}: {}", crate::fortran_fmt::es(total, 12, 5)), 6);
            }
        }
        Ok(())
    }
}

/// Whether a local system adds anything (a zero boundary contribution
/// leaves the bulk system as it was).
pub(crate) fn is_nonzero(local: &elmer_core::LocalSystem) -> bool {
    local.stiff.iter().chain(&local.force).any(|&v| v != 0.0)
}
