//! `SaveScalars` (`fem/src/modules/SaveData/SaveScalars.F90`) as a scheduled
//! [`Solver`].
//!
//! One visit walks the `Variable n` / `Operator n` keywords the way the
//! Fortran does (an operator or variable that is left out is inherited from
//! the previous number), collects `(name, value)` pairs, appends one row to
//! the `.dat` file through `elmer-saveio`, stores every value in the
//! `Simulation` section for later solvers and exports `Target Variable n`
//! as a global variable.

use std::path::PathBuf;

use elmer_core::{normalize_key, FieldKind, ParamMap, Params, Value};
use elmer_saveio::fortran::{fortran_es_e, fortran_i};
use elmer_saveio::savedata::{append_text, scalar_column_name, NamesMeta, ScalarScope, ScalarsWriter};
use elmer_solver::fortran_fmt::g;
use elmer_solver::{Solver, SolverCtx, SolverError, VarId, Variable};

use crate::integrals::{
    bc_is_active, boundary_integrals, boundary_statistics, bulk_integral, node_mask, vector_mean_deviation,
    vector_statistics, BulkRequest, Fe,
};
use crate::plan::{parse_operator, MaskKind, OperKind, PseudoNorm, ScalarPlan};
use crate::view::{component_name, resolve, Resolved};

const CALLER: &str = "SaveScalars";

/// Variables the Fortran saves as integers (SaveScalars.F90:503-506).
const INTEGER_GLOBALS: [&str; 4] = ["timestep", "nonlin iter", "coupled iter", "run"];

/// The `ValueNames`/`Values` lists and the host-associated `VariableName` of
/// the Fortran: `AddToSaveList` overwrites that variable with its `Target
/// Variable n` lookup, which is blank when the keyword is absent. The
/// column names of the second and later boundary conditions of a boundary
/// operator therefore come out without a variable name.
#[derive(Default)]
struct SaveList {
    names: Vec<String>,
    values: Vec<f64>,
    integer: Vec<bool>,
    variable_name: String,
}

impl SaveList {
    fn add(&mut self, plan: &ScalarPlan, ctx: &mut SolverCtx<'_>, name: &str, val: f64, integer: bool) -> Result<(), SolverError> {
        self.names.push(name.to_lowercase());
        self.values.push(val);
        self.integer.push(integer);
        let n = self.values.len();
        let target = plan.section.string(&format!("target variable {n}"));
        self.variable_name = target.unwrap_or("").to_string();
        if let Some(target) = target {
            match ctx.sim.vars.get_mut(target) {
                Some(v) if v.kind == FieldKind::Global => v.values[0] = val,
                Some(_) => {
                    return Err(SolverError::fatal(
                        CALLER,
                        format!("Target Variable {n}: '{target}' exists and is not a global scalar"),
                    ))
                }
                None => {
                    ctx.sim.vars.add(Variable::global(target, val))?;
                }
            }
            ctx.info(CALLER, &format!("Defining: {target} = {}", self.names[n - 1]), 8);
        }
        Ok(())
    }
}

/// The `SaveScalars` solver.
pub struct SaveScalars {
    plan: ScalarPlan,
    fe: Option<Fe>,
    meta: fn() -> NamesMeta,
    writer: Option<ScalarsWriter>,
    closest: Option<Vec<usize>>,
}

impl SaveScalars {
    pub(crate) fn new(plan: ScalarPlan, fe: Option<Fe>, meta: fn() -> NamesMeta) -> Self {
        SaveScalars {
            plan,
            fe,
            meta,
            writer: None,
            closest: None,
        }
    }

    pub fn plan(&self) -> &ScalarPlan {
        &self.plan
    }

    fn fe(&self) -> Result<&Fe, SolverError> {
        self.fe
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("element basis and quadrature (elmer-fe)".into()))
    }
}

fn slot_index(ctx: &SolverCtx<'_>) -> usize {
    let id = ctx.solver_id();
    ctx.sim.slots.iter().position(|s| s.id == id).expect("running solver has a slot")
}

/// `ListGetLogicalAnyBody` and friends: some section of the kind has
/// `name` true.
fn any_mask(ctx: &SolverCtx<'_>, kind: MaskKind, name: &str) -> bool {
    let m = &ctx.sim.model;
    let sections: &[ParamMap] = match kind {
        MaskKind::Body => &m.bodies,
        MaskKind::BodyForce => &m.body_forces,
        MaskKind::Material => &m.materials,
        MaskKind::None => return true,
    };
    sections.iter().any(|s| s.logical(name) == Some(true))
}

impl Solver for SaveScalars {
    /// `SaveScalars_init`: a solver that reports a norm gets the global
    /// variable `<Equation>_var`, unless the section names a `Variable`.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if !self.plan.reports_norm() || self.plan.section.get("variable").is_some() {
            return Ok(());
        }
        let eq = self.plan.equation.clone().unwrap_or_else(|| CALLER.to_string());
        let name = normalize_key(&format!("{eq}_var"));
        let id = match ctx.sim.vars.id(&name) {
            Some(id) => id,
            None => {
                let mut v = Variable::global(&name, 0.0);
                v.output = false;
                ctx.sim.vars.add(v)?
            }
        };
        let ix = slot_index(ctx);
        ctx.sim.slots[ix].variable = Some(id);
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.comm.size() > 1 {
            return Err(SolverError::Unsupported("SaveScalars: MPI reductions".into()));
        }
        ctx.info(CALLER, "-----------------------------------------", 4);
        ctx.info(CALLER, "Saving scalar values of various kinds", 4);

        let mut list = SaveList::default();
        self.collect(ctx, &mut list)?;
        if list.values.is_empty() {
            ctx.warn(CALLER, "Found no values to save");
            return Ok(());
        }
        ctx.info(CALLER, &format!("Found {} values to save in total", list.values.len()), 6);

        if self.plan.filename.is_some() {
            self.write_files(ctx, &list)?;
        }
        if self.plan.echo_values {
            ctx.info(CALLER, "Showing computed results:", 3);
            for (i, (n, v)) in list.names.iter().zip(&list.values).enumerate() {
                ctx.info(
                    CALLER,
                    &format!("{}: {n}{}", fortran_i(i as i64 + 1, 4), fortran_es_e(*v, 22, 12, Some(3))),
                    3,
                );
            }
        }

        // Data for later solvers in the Simulation section.
        for (n, v) in list.names.iter().zip(&list.values) {
            ctx.sim.model.simulation.insert(n, Value::real(*v));
        }
        self.set_pseudo_norm(ctx, &list)?;
        if let Some(i) = self.plan.cost_function_index.filter(|&i| i <= list.values.len()) {
            ctx.sim.model.simulation.insert("Cost Function", Value::real(list.values[i - 1]));
        }
        ctx.info(CALLER, "-----------------------------------------", 7);
        Ok(())
    }
}

impl SaveScalars {
    /// The `Variable n` / `Operator n` walk, then the point values and the
    /// scalars other solvers left in the Simulation section.
    fn collect(&mut self, ctx: &mut SolverCtx<'_>, list: &mut SaveList) -> Result<(), SolverError> {
        let plan = self.plan.clone();
        let p = &plan.section;
        let n_nodes = ctx.mesh.n_nodes();

        let mut got_var = true;
        let mut got_oper = false;
        let mut got_old_oper = false;
        let mut old_var: Option<String> = None;
        let mut old_variable_name = String::new();
        let mut old_oper0 = String::new();
        let mut old_mask = "save scalars".to_string();
        let mut nodal_oper = false;
        let mut got_nodal = false;
        let mut simulation_visited = false;
        let mut n = 0usize;

        while got_var || got_oper {
            got_oper = false;
            n += 1;
            let var_kw = p.string(&format!("variable {n}")).map(|s| normalize_key(s.trim_matches('"')));
            got_var = var_kw.is_some();
            list.variable_name = var_kw.clone().unwrap_or_default();

            let mut have_var = false;
            if got_var {
                let name = list.variable_name.clone();
                match resolve(&ctx.sim.vars, &name)? {
                    Some(Resolved::Global(id)) => {
                        old_var = Some(name.clone());
                        old_variable_name = name.clone();
                        let v = ctx.sim.vars.by_id(id);
                        if v.dofs == 1 {
                            let (val, int) = (v.values[0], INTEGER_GLOBALS.contains(&name.as_str()));
                            list.add(&plan, ctx, &format!("value: {name} scalar variable"), val, int)?;
                        } else {
                            for j in 1..=v.dofs {
                                let val = ctx.sim.vars.by_id(id).values[j - 1];
                                let label = format!("value: {} scalar variable", component_name(&name, j));
                                list.add(&plan, ctx, &label, val, false)?;
                            }
                        }
                        continue;
                    }
                    Some(_) => {
                        have_var = true;
                        old_var = Some(name.clone());
                        old_variable_name = name;
                    }
                    None => {
                        if let Some(val) = ctx.simulation_params().real(&name) {
                            ctx.info(CALLER, &format!("Adding keyword from simulation section as variable: {name}"), 10);
                            list.add(&plan, ctx, &name, val, false)?;
                            simulation_visited = true;
                        } else {
                            ctx.warn(CALLER, &format!("Requested variable does not exist: {name}"));
                        }
                        continue;
                    }
                }
                let kw = p.logical(&format!("nodal variable {n}"));
                got_nodal = kw.is_some();
                nodal_oper = kw.unwrap_or(false);
            } else if old_var.is_some() {
                list.variable_name = old_variable_name.clone();
                have_var = true;
            }
            if let Some(b) = p.logical(&format!("nodal variable {n}")) {
                got_nodal = true;
                nodal_oper = b;
            }

            let oper_kw = p.string(&format!("operator {n}"));
            got_oper = oper_kw.is_some();
            let mut oper0 = oper_kw.map(|s| normalize_key(s.trim_matches('"'))).unwrap_or_default();
            if !got_oper && got_old_oper {
                oper0 = old_oper0.clone();
            }
            let mask_kw = p.string(&format!("mask name {n}")).map(normalize_key);
            let got_mask_name = mask_kw.is_some();
            let mask_name = match mask_kw {
                Some(m) => {
                    old_mask = m.clone();
                    m
                }
                None => old_mask.clone(),
            };
            if !(got_oper || got_var || got_mask_name) {
                continue;
            }
            if got_oper {
                old_oper0 = oper0.clone();
                got_old_oper = true;
            } else if got_old_oper {
                oper0 = old_oper0.clone();
                got_oper = true;
            } else {
                ctx.info(CALLER, &format!("No operator given for variable: {}", list.variable_name), 3);
                continue;
            }

            let op = parse_operator(&oper0)?;
            let passive = p
                .integer("passive coordinate")
                .or_else(|| p.integer(&format!("passive coordinate {n}")))
                .filter(|&c| c > 0)
                .map(|c| c as usize);
            if op.mask != MaskKind::None && !any_mask(ctx, op.mask, &mask_name) {
                ctx.warn(CALLER, &format!("Masked operators require mask: {mask_name}"));
            }
            if op.kind != OperKind::Mesh && !have_var {
                ctx.warn(CALLER, &format!("Operator > {} < requires variable", op.oper));
                continue;
            }
            let variable = list.variable_name.clone();
            let mut save_name = String::new();
            if have_var {
                save_name = format!("{}: {variable}", op.oper0);
                if got_mask_name {
                    save_name = format!("{save_name} mask {mask_name}");
                }
                if got_nodal {
                    save_name.push_str(if nodal_oper { " nodal" } else { " non-nodal" });
                }
            }

            match op.kind {
                OperKind::Mesh => self.mesh_operator(ctx, list, &plan, &op.oper, &variable)?,
                OperKind::Meta => {
                    let (val, int) = meta_value(ctx, &variable, &op.oper)?;
                    list.add(&plan, ctx, &save_name, val, int)?;
                }
                OperKind::Statistic | OperKind::Deviation => {
                    let val = {
                        let view = match resolve(&ctx.sim.vars, &variable)? {
                            Some(Resolved::Field(v)) => v,
                            Some(Resolved::Components { first, .. }) => first,
                            _ => return Err(SolverError::Unsupported(format!("SaveScalars: {} of a global variable", op.oper))),
                        };
                        let mask = (op.mask != MaskKind::None)
                            .then(|| node_mask(ctx, op.mask, &mask_name, view.n_perm()));
                        if op.kind == OperKind::Deviation {
                            vector_mean_deviation(&view, mask.as_deref())
                        } else {
                            vector_statistics(&view, &op.oper, got_nodal.then_some(nodal_oper), mask.as_deref(), n_nodes)
                        }
                    };
                    list.add(&plan, ctx, &save_name, val, false)?;
                }
                OperKind::Bulk => {
                    let val = {
                        let view = match resolve(&ctx.sim.vars, &variable)? {
                            Some(Resolved::Field(v)) => v,
                            Some(Resolved::Components { first, .. }) => first,
                            _ => return Err(SolverError::Unsupported(format!("SaveScalars: {} of a global variable", op.oper))),
                        };
                        let req = BulkRequest {
                            oper: &op.oper,
                            mask: op.mask,
                            mask_name: &mask_name,
                            passive,
                        };
                        bulk_integral(ctx, self.fe()?, &view, &req)?
                    };
                    list.add(&plan, ctx, &save_name, val, false)?;
                }
                OperKind::BoundaryStatistic => {
                    if !ctx.sim.model.boundary_conditions.iter().any(|bc| bc_is_active(bc, &mask_name)) {
                        ctx.warn(CALLER, &format!("No flag > {mask_name} < active for operator: {}", op.oper));
                        continue;
                    }
                    let (mut val, count) = {
                        let view = match resolve(&ctx.sim.vars, &variable)? {
                            Some(Resolved::Field(v)) => v,
                            Some(Resolved::Components { first, .. }) => first,
                            _ => return Err(SolverError::Unsupported(format!("SaveScalars: {} of a global variable", op.oper))),
                        };
                        boundary_statistics(ctx, &view, &op.oper, &mask_name)
                    };
                    if op.oper == "boundary mean" && count > 0 {
                        val /= count as f64;
                    }
                    let name = scalar_column_name(&op.oper, &variable, &ScalarScope::BcMask(mask_name.clone()));
                    list.add(&plan, ctx, &name, val, false)?;
                }
                OperKind::BoundaryIntegral => {
                    let active: Vec<bool> = ctx
                        .sim
                        .model
                        .boundary_conditions
                        .iter()
                        .map(|bc| bc_is_active(bc, &mask_name))
                        .collect();
                    if !active.iter().any(|&a| a) {
                        return Err(SolverError::fatal(
                            CALLER,
                            format!("No flag > {mask_name}< active for operator: {}", op.oper),
                        ));
                    }
                    let totals = {
                        let view = match resolve(&ctx.sim.vars, &variable)? {
                            Some(Resolved::Field(v)) => Some(v),
                            Some(Resolved::Components { .. }) => {
                                return Err(SolverError::Unsupported(format!(
                                    "SaveScalars: {} of a component variable",
                                    op.oper
                                )))
                            }
                            _ => None,
                        };
                        boundary_integrals(ctx, self.fe()?, view.as_ref(), &op.oper, &mask_name, passive)?
                    };
                    for (j, _) in active.iter().enumerate().filter(|(_, &a)| a) {
                        let mut val = totals.fluxes[j];
                        if op.oper == "boundary int mean" && totals.areas[j] > 0.0 {
                            val /= totals.areas[j];
                        }
                        // The variable name here is whatever AddToSaveList left behind.
                        let name = scalar_column_name(&op.oper, &list.variable_name, &ScalarScope::Bc(j + 1));
                        list.add(&plan, ctx, &name, val, false)?;
                    }
                    // `Save Flux Range` (default true): the point flux range
                    // follows the per-BC values (SaveScalars.F90:888).
                    if op.oper == "diffusive flux" && plan.section.logical("save flux range").unwrap_or(true) {
                        let var = list.variable_name.clone();
                        list.add(&plan, ctx, &format!("min {}: {var}", op.oper), totals.range.0, false)?;
                        list.add(&plan, ctx, &format!("max {}: {var}", op.oper), totals.range.1, false)?;
                    }
                }
            }
        }
        ctx.info(CALLER, &format!("Performed {} reduction operations", n.saturating_sub(1)), 7);

        self.point_values(ctx, list, &plan)?;

        // Scalars other solvers left in the Simulation section.
        let prefix: String = plan.results_prefix.chars().take(4).collect();
        let found: Vec<(String, f64)> = ctx
            .simulation_params()
            .0
            .iter()
            .filter(|(k, _)| {
                let head: String = k.chars().take(4).collect();
                format!("{head:<4}") == format!("{prefix:<4}") && plan.results_prefix.len() <= 4
            })
            .filter_map(|(k, v)| match v {
                Value::Real { data, .. } => data.first().map(|x| (k.clone(), *x)),
                _ => None,
            })
            .collect();
        for (k, v) in found {
            if simulation_visited && list.names.contains(&k.to_lowercase()) {
                continue;
            }
            list.add(&plan, ctx, &k, v, false)?;
        }
        Ok(())
    }

    fn mesh_operator(
        &self,
        ctx: &mut SolverCtx<'_>,
        list: &mut SaveList,
        plan: &ScalarPlan,
        oper: &str,
        variable: &str,
    ) -> Result<(), SolverError> {
        match oper {
            "partitions" => list.add(plan, ctx, "value: number of partitions", ctx.comm.size() as f64, true),
            // Serial runs: the Fortran only fills these in for parallel meshes.
            "partition checksum" => list.add(plan, ctx, "value: partition checksum", 0.0, false),
            "partition neighbours checksum" => list.add(plan, ctx, "value: partition neighbours checksum", 0.0, false),
            "nodes" => list.add(plan, ctx, &format!("{oper}: {variable}"), ctx.mesh.n_nodes() as f64, true),
            "elements" => list.add(plan, ctx, &format!("{oper}: {variable}"), ctx.mesh.n_bulk() as f64, true),
            "bounding box" => {
                let (mut lo, mut hi) = ([f64::MAX; 3], [-f64::MAX; 3]);
                for n in 0..ctx.mesh.n_nodes() as u32 {
                    let x = ctx.mesh.coord(n);
                    for d in 0..3 {
                        lo[d] = lo[d].min(x[d]);
                        hi[d] = hi[d].max(x[d]);
                    }
                }
                let ndim = if ctx.mesh.dim() > 2 { 3 } else { 2 };
                for (d, axis) in ["x", "y", "z"].iter().enumerate().take(ndim) {
                    list.add(plan, ctx, &format!("{oper} min {axis}"), lo[d], false)?;
                    list.add(plan, ctx, &format!("{oper} max {axis}"), hi[d], false)?;
                }
                Ok(())
            }
            _ => unreachable!("mesh operator {oper} is classified in plan.rs"),
        }
    }

    /// `Save Points` and the nearest nodes of `Save Coordinates`: every
    /// scalar field with `Output` set, at each node (SaveScalars.F90:1000).
    fn point_values(&mut self, ctx: &mut SolverCtx<'_>, list: &mut SaveList, plan: &ScalarPlan) -> Result<(), SolverError> {
        if plan.save_points.is_empty() && plan.save_coordinates.is_empty() {
            return Ok(());
        }
        let mut nodes: Vec<usize> = Vec::new();
        for &l in &plan.save_points {
            if l < 0 || l as usize > ctx.mesh.n_nodes() {
                return Err(SolverError::fatal(CALLER, format!("Save Points: node {l} is not in the mesh")));
            }
            nodes.push(l as usize);
        }
        if self.closest.is_none() {
            let mesh = ctx.mesh;
            let found = plan
                .save_coordinates
                .iter()
                .map(|c| {
                    let mut best = (f64::MAX, 0usize);
                    for n in 0..mesh.n_nodes() {
                        let x = mesh.coord(n as u32);
                        let d: f64 = (0..3).map(|k| (x[k] - c[k]).powi(2)).sum();
                        if d < best.0 {
                            best = (d, n + 1);
                        }
                    }
                    best.1
                })
                .collect();
            self.closest = Some(found);
        }
        nodes.extend(self.closest.iter().flatten());

        for l in nodes {
            if l == 0 {
                continue;
            }
            // The mesh coordinates come first, then the variables in
            // registration order; vector fields appear as their components.
            let mut entries: Vec<(String, Option<f64>)> = Vec::new();
            let x = ctx.mesh.coord(l as u32 - 1);
            for (k, c) in x.iter().enumerate() {
                entries.push((format!("coordinate {}", k + 1), Some(*c)));
            }
            for (_, v) in ctx.sim.vars.iter() {
                if !v.output || v.kind == FieldKind::Global {
                    continue;
                }
                if v.kind != FieldKind::Nodal {
                    return Err(SolverError::Unsupported(format!(
                        "SaveScalars: point values of {:?} variable '{}'",
                        v.kind, v.name
                    )));
                }
                let slot = v.slot(l as u32 - 1);
                for k in 0..v.dofs {
                    let name = if v.dofs == 1 { v.name.clone() } else { component_name(&v.name, k + 1) };
                    entries.push((name, slot.map(|s| v.values[s * v.dofs + k])));
                }
            }
            for (name, val) in entries {
                if let Some(val) = val {
                    list.add(plan, ctx, &format!("value: {name} at node {l}"), val, false)?;
                }
            }
        }
        ctx.info(
            CALLER,
            &format!("Tabulated all field values at {} points", plan.save_points.len() + plan.save_coordinates.len()),
            3,
        );
        Ok(())
    }

    fn write_files(&mut self, ctx: &mut SolverCtx<'_>, list: &SaveList) -> Result<(), SolverError> {
        let visited = ctx.sim.slots[slot_index(ctx)].times_visited;
        let (display, path) = output_path(ctx, self.plan.filename.as_deref().unwrap_or_default(), &self.plan, visited == 0)?;
        if self.writer.is_none() {
            let mut meta = (self.meta)();
            meta.simulation_comment = ctx.simulation_params().string("comment").map(str::to_string);
            meta.comment = self.plan.section.string("comment").map(str::to_string);
            let mut w = ScalarsWriter::new(&path, &display, meta).with_file_append(self.plan.file_append);
            if self.plan.line_marker != 0 {
                w = w.with_line_marker(self.plan.line_marker);
            }
            self.writer = Some(w);
        }
        ctx.info(CALLER, &format!("Saving values to file: {display}"), 4);
        let row: Vec<(String, f64, bool)> = list
            .names
            .iter()
            .zip(&list.values)
            .zip(&list.integer)
            .map(|((n, v), i)| (n.clone(), *v, *i))
            .collect();
        self.writer
            .as_mut()
            .expect("created above")
            .row_flagged(&row)
            .map_err(|e| SolverError::fatal(CALLER, format!("Could not save {display}: {e}")))?;

        if visited == 0 && self.plan.file_append && self.plan.line_marker != 0 {
            if let Some(comment) = self.plan.section.string("comment") {
                let text = format!("{}: {comment}\n", fortran_i(self.plan.line_marker.into(), 6));
                append_text(&PathBuf::from(format!("{}.marker", path.display())), &text, true)
                    .map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
            }
        }
        if self.plan.live_graph {
            let csv = PathBuf::from(format!("{}.csv", path.display()));
            let append = visited > 0 || self.plan.file_append;
            let mut text = String::new();
            if !append {
                text.push_str(&list.names.join(","));
                text.push('\n');
            }
            let cells: Vec<String> = list.values.iter().map(|v| fortran_es_e(*v, 22, 12, Some(3))).collect();
            text.push_str(&cells.join(","));
            text.push('\n');
            append_text(&csv, &text, append).map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
        }
        Ok(())
    }

    /// `SetPseudoNorm` (SaveScalars.F90:3434).
    fn set_pseudo_norm(&self, ctx: &mut SolverCtx<'_>, list: &SaveList) -> Result<(), SolverError> {
        let Some(spec) = &self.plan.pseudo_norm else { return Ok(()) };
        let nvalues = list.values.len();
        let norm = match spec {
            PseudoNorm::Name(name) => match list.names.iter().position(|n| n == name) {
                Some(i) => {
                    ctx.info(CALLER, &format!("Associating scalar {} to norm!", i + 1), 8);
                    Some(list.values[i])
                }
                None => {
                    ctx.warn(CALLER, &format!("Could not find scalar for norm: {name}"));
                    None
                }
            },
            PseudoNorm::Index(i) => (*i <= nvalues).then(|| list.values[i - 1]),
            PseudoNorm::Reference(refs) => {
                if refs.len() > nvalues {
                    return Err(SolverError::fatal(
                        CALLER,
                        "Size of \"Reference Values\" bigger than computed values!",
                    ));
                }
                let sum: f64 = refs
                    .iter()
                    .zip(&list.values)
                    .map(|(r, v)| {
                        if r.abs() > f64::EPSILON {
                            let c = v / r;
                            c.max(1.0 / c)
                        } else {
                            1.0 + v.abs()
                        }
                    })
                    .sum();
                Some(sum / refs.len() as f64)
            }
        };
        let Some(norm) = norm else { return Ok(()) };
        if matches!(spec, PseudoNorm::Reference(_)) {
            // By construction the reference norm is 1 (ListAddNewConstReal).
            let ix = slot_index(ctx);
            if ctx.sim.slots[ix].params.get("reference norm").is_none() {
                ctx.sim.slots[ix].params.insert("Reference Norm", Value::real(1.0));
            }
        }
        let ix = slot_index(ctx);
        if let Some(id) = ctx.sim.slots[ix].variable {
            set_variable_norm(ctx, id, norm);
        }
        // `Name` is the equation name for the indexed and named norms. The
        // `Reference Values` branch never sets it, so the line ends empty.
        let name = match spec {
            PseudoNorm::Reference(_) => String::new(),
            _ => self.plan.equation.clone().unwrap_or_else(|| CALLER.to_string()),
        };
        let msg = format!("SS (ITER=1) (NRM,RELC): ({}{} ) :: {name}", g(norm, 15, 8), g(0.0, 15, 8));
        ctx.info("ComputeChange", &msg, 3);
        Ok(())
    }
}

fn set_variable_norm(ctx: &mut SolverCtx<'_>, id: VarId, norm: f64) {
    let v = ctx.sim.vars.by_id_mut(id);
    v.values.iter_mut().for_each(|x| *x = norm);
    v.conv.norm = norm.abs();
}

/// `norm`, `nonlin change`, `dofs`, ... of a primary variable.
fn meta_value(ctx: &SolverCtx<'_>, variable: &str, oper: &str) -> Result<(f64, bool), SolverError> {
    let id = ctx.sim.vars.id(variable).ok_or_else(|| {
        SolverError::Unsupported(format!("SaveScalars: '{oper}' of '{variable}', which is not a primary variable"))
    })?;
    let v = ctx.sim.vars.by_id(id);
    Ok(match oper {
        "norm" => (v.conv.norm, false),
        "nonlin change" => (v.conv.nonlin_change, false),
        "steady state change" => (v.conv.steady_change, false),
        "nonlin iter" => (v.conv.nonlin_iter as f64, true),
        "nonlin converged" => (v.conv.nonlin_converged as f64, true),
        "steady converged" => (v.conv.steady_converged as f64, true),
        "dofs" => (v.values.len() as f64, true),
        "nans" => (v.values.iter().filter(|x| x.is_nan()).count() as f64, true),
        _ => unreachable!("meta operator {oper} is classified in plan.rs"),
    })
}

/// `SolverOutputDirectory` (SaveUtils.F90:788) plus the `Partition
/// Numbering` suffix: the file name as it appears in the metadata and the
/// path to write. Creates a one-level directory on the first visit.
pub(crate) fn output_path(
    ctx: &SolverCtx<'_>,
    filename: &str,
    plan: &ScalarPlan,
    first_visit: bool,
) -> Result<(String, PathBuf), SolverError> {
    let mut file = filename.to_string();
    if plan.partition_numbering {
        let np = ctx.comm.size();
        file = match file.rfind('.') {
            Some(i) => format!("{}_np{np}{}", &file[..i], &file[i..]),
            None => format!("{file}_np{np}"),
        };
    }
    split_output_path(ctx, plan.section.string("output directory"), &file, first_visit)
}

pub(crate) fn split_output_path(
    ctx: &SolverCtx<'_>,
    solver_dir: Option<&str>,
    filename: &str,
    first_visit: bool,
) -> Result<(String, PathBuf), SolverError> {
    let found = solver_dir.or_else(|| ctx.simulation_params().string("output directory"));
    let mut dir = found.unwrap_or(".").trim().to_string();
    let mut nd = if found.is_none() && dir == "." { 0 } else { dir.len() };
    let abs = filename.contains(':') || filename.starts_with('/') || filename.starts_with('\\');
    if nd > 0 && !abs {
        if let Some(rest) = dir.strip_prefix("~/") {
            if let Ok(home) = std::env::var("HOME") {
                dir = format!("{home}/{rest}");
                nd = dir.len();
            }
        }
        if first_visit {
            let _ = std::fs::create_dir(&dir);
        }
    }
    let mut file = filename.to_string();
    if let Some(n) = file.find('/').map(|i| i + 1).filter(|&n| n > 2) {
        if abs || nd == 0 {
            dir = file[..n - 1].to_string();
        } else {
            dir = format!("{dir}/{}", &file[..n - 1]);
        }
        file = file[n..].to_string();
        if first_visit {
            let _ = std::fs::create_dir(&dir);
        }
    }
    if ctx.simulation_params().logical("output partitioning directory") == Some(true) {
        return Err(SolverError::Unsupported("Output Partitioning Directory".into()));
    }
    let display = format!("{dir}/{file}");
    let path = PathBuf::from(&display);
    Ok((display, path))
}
