//! `SaveLine` (`fem/src/modules/SaveData/SaveLine.F90`) for the existing-
//! boundary mode as a scheduled [`Solver`].
//!
//! Every boundary element whose boundary condition has the `Save Mask`
//! keyword (default `Save Line`) true contributes its nodes. The nodes are
//! numbered by first touch and, unless `Optimize Node Ordering = False`, in
//! Elmer's bandwidth-optimal order for the line (MeshBasics.F90
//! `MakePermUsingMask`), so a boundary curve comes out as an ordered
//! polyline. One row per node: the label columns, the three coordinates and
//! the requested variables. Polylines, circles, isocurves, axis lines, heat
//! flux columns and intersections with element edges are refused when the
//! solver is built.

use elmer_assembly::{rcm_order, ListMatrix};
use elmer_core::{ElemId, Params};
use elmer_saveio::savedata::{append_text, line_names_text, line_row_text, LineColumns, LineNodeCounts, NamesMeta};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::plan::LinePlan;
use crate::scalars::split_output_path;
use crate::view::{component_name, resolve, Resolved, VarView};

const CALLER: &str = "SaveLine";

/// One value column of the file.
struct Column<'a> {
    name: String,
    view: VarView<'a>,
}

/// Nodes of the saved line in file order, each with the boundary condition
/// that last touched it.
#[derive(Clone, Debug, Default)]
struct LineNodes {
    /// 0-based node indices.
    nodes: Vec<usize>,
    /// 1-based BC id per entry of `nodes`.
    bc: Vec<usize>,
}

/// The `SaveLine` solver.
pub struct SaveLine {
    plan: LinePlan,
    meta: fn() -> NamesMeta,
    line: Option<LineNodes>,
}

impl SaveLine {
    pub(crate) fn new(plan: LinePlan, meta: fn() -> NamesMeta) -> Self {
        SaveLine {
            plan,
            meta,
            line: None,
        }
    }

    pub fn plan(&self) -> &LinePlan {
        &self.plan
    }
}

/// `VariableGetN` for every `Variable n` up to the first one that does not
/// resolve, then the columns they add.
fn columns<'a>(ctx: &'a SolverCtx<'_>, plan: &LinePlan) -> Result<Vec<Column<'a>>, SolverError> {
    let mut cols = Vec::new();
    for (i, name) in plan.variables.iter().enumerate() {
        match resolve(&ctx.sim.vars, name)? {
            None => {
                ctx.warn("VariableGetN", &format!("Variable given but not found: {name}"));
                break;
            }
            Some(Resolved::Global(_)) => {
                return Err(SolverError::Unsupported(format!(
                    "{CALLER}: Variable {} is a global variable",
                    i + 1
                )))
            }
            Some(Resolved::Field(v)) if v.ncomp == 1 => cols.push(Column {
                name: name.clone(),
                view: v,
            }),
            Some(Resolved::Field(v)) => {
                for k in 0..v.ncomp {
                    cols.push(Column {
                        name: component_name(name, k + 1),
                        view: VarView {
                            offset: v.offset + k,
                            ncomp: 1,
                            ..v
                        },
                    });
                }
            }
            Some(Resolved::Components { count, .. }) => {
                for k in 1..=count {
                    let alias = component_name(name, k);
                    let Some(Resolved::Field(v)) = resolve(&ctx.sim.vars, &alias)? else {
                        return Err(SolverError::fatal(CALLER, format!("component {alias} vanished")));
                    };
                    cols.push(Column { name: alias, view: v });
                }
            }
        }
    }
    Ok(cols)
}

/// `MakePermUsingMask` over the boundary elements and the `BoundaryIndex`
/// bookkeeping of `SaveExistingLines`.
fn find_line(ctx: &SolverCtx<'_>, plan: &LinePlan) -> Result<LineNodes, SolverError> {
    let m = &ctx.sim.model;
    if m.body_forces.iter().any(|b| b.get(&plan.mask).is_some()) {
        return Err(SolverError::Unsupported(format!(
            "{CALLER}: Save Mask '{}' on a Body Force",
            plan.mask
        )));
    }
    let nb = ctx.mesh.n_bulk() as u32;
    let nbound = ctx.mesh.n_boundary() as u32;
    // Elements whose BC has the mask true, in mesh order.
    let hit: Vec<ElemId> = (nb..nb + nbound)
        .map(ElemId)
        .filter(|&e| ctx.bc_params(e).is_some_and(|bc| bc.logical(&plan.mask) == Some(true)))
        .collect();
    if hit.is_empty() {
        return Ok(LineNodes::default());
    }

    // First-touch numbering, 1-based as in the Fortran.
    let n_all = ctx.mesh.n_nodes();
    let mut perm = vec![0usize; n_all];
    let mut k = 0;
    for &e in &hit {
        for &n in ctx.mesh.element(e).nodes {
            if perm[n as usize] == 0 {
                k += 1;
                perm[n as usize] = k;
            }
        }
    }

    if plan.optimize_order {
        // The first node of the numbering becomes the corner with the least
        // x + y + z.
        let mut best = (f64::MAX, 0usize);
        for (i, &p) in perm.iter().enumerate() {
            if p > 0 {
                let c = ctx.mesh.coord(i as u32);
                let d = c[0] + c[1] + c[2];
                if d < best.0 {
                    best = (d, i);
                }
            }
        }
        let first = perm.iter().position(|&p| p == 1).expect("node numbered 1");
        perm[first] = perm[best.1];
        perm[best.1] = 1;

        let mut list = ListMatrix::new(k);
        for &e in &hit {
            let nodes = ctx.mesh.element(e).nodes;
            for &a in nodes {
                for &b in nodes {
                    list.add(perm[a as usize] - 1, (perm[b as usize] - 1) as u32);
                }
            }
        }
        let new_index = rcm_order(&list);
        for p in perm.iter_mut().filter(|p| **p > 0) {
            *p = new_index[*p - 1] as usize + 1;
        }
    }

    let mut inv = vec![0usize; k];
    for (i, &p) in perm.iter().enumerate() {
        if p > 0 {
            inv[p - 1] = i;
        }
    }
    // The last boundary element that touches a node names its BC. Like the
    // Fortran, only elements whose nodes are all on the line count, and the
    // keyword has to be present on the BC (its value is not looked at).
    let mut bc_of = vec![0usize; k];
    for e in (nb..nb + nbound).map(ElemId) {
        let nodes = ctx.mesh.element(e).nodes;
        if nodes.iter().any(|&n| perm[n as usize] == 0) {
            continue;
        }
        let Some(id) = ctx.bc_id(e) else { continue };
        if ctx.sim.model.boundary_conditions[id - 1].get(&plan.mask).is_none() {
            continue;
        }
        for &n in nodes {
            bc_of[perm[n as usize] - 1] = id;
        }
    }
    Ok(LineNodes { nodes: inv, bc: bc_of })
}

impl Solver for SaveLine {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if ctx.comm.size() > 1 {
            return Err(SolverError::Unsupported("SaveLine: MPI reductions".into()));
        }
        let plan = self.plan.clone();
        ctx.info(CALLER, "-----------------------------------------", 4);
        ctx.info(CALLER, &format!("Saving data on lines to file: {}", plan.filename), 4);

        let slot = ctx.sim.slots.iter().position(|s| s.id == ctx.solver_id()).expect("running solver has a slot");
        let visited = ctx.sim.slots[slot].times_visited;
        let done_time = ctx.sim.slots[slot].done_time;

        if self.line.is_none() {
            self.line = Some(find_line(ctx, &plan)?);
        }
        let line = self.line.as_ref().expect("set above");
        let cols = columns(ctx, &plan)?;
        // Coordinates plus the variables.
        let n_results = 3 + cols.len();
        if cols.is_empty() {
            ctx.info(CALLER, "All done", 3);
            return Ok(());
        }

        let (display, path) = split_output_path(ctx, plan.section.string("output directory"), &plan.filename, visited == 0)?;

        let mut text = String::new();
        let mut found_nan = 0usize;
        for (idx, &node) in line.nodes.iter().enumerate() {
            let mut labels: Vec<i64> = Vec::new();
            if !plan.skip_boundary_info {
                if plan.line_marker != 0 {
                    labels.push(plan.line_marker.into());
                }
                if transient {
                    labels.push(done_time as i64);
                }
                labels.push(visited as i64 + 1);
                labels.push(line.bc[idx] as i64);
                labels.push(node as i64 + 1);
            }
            let mut values = Vec::with_capacity(n_results);
            values.extend(ctx.mesh.coord(node as u32));
            for c in &cols {
                values.push(c.view.slot(node).map_or(0.0, |s| c.view.at(s, 0)));
            }
            for v in &mut values {
                if v.is_nan() {
                    found_nan += 1;
                    *v = -1.0;
                }
            }
            text.push_str(&line_row_text(&labels, &values));
        }
        if !text.is_empty() {
            append_text(&path, &text, visited > 0 || plan.file_append)
                .map_err(|e| SolverError::fatal(CALLER, format!("Problems saving {display}: {e}")))?;
        }

        if visited == 0 {
            let mut meta = (self.meta)();
            meta.simulation_comment = ctx.simulation_params().string("comment").map(str::to_string);
            meta.comment = plan.section.string("comment").map(str::to_string);
            let mut names = vec!["coordinate 1".to_string(), "coordinate 2".into(), "coordinate 3".into()];
            names.extend(cols.iter().map(|c| c.name.clone()));
            let counts = LineNodeCounts {
                bc: line.nodes.len(),
                ..Default::default()
            };
            let lc = LineColumns {
                line_marker: plan.line_marker != 0,
                timestep: transient,
                skip_boundary_info: plan.skip_boundary_info,
            };
            let names_path = format!("{display}.names");
            std::fs::write(&names_path, line_names_text(&meta, &display, counts, lc, &names))
                .map_err(|e| SolverError::fatal(CALLER, format!("Could not save {names_path}: {e}")))?;
        }
        if found_nan > 0 {
            ctx.warn(CALLER, &format!("Replaced {found_nan} NaN entries with -1"));
        }
        ctx.info(CALLER, "All done", 3);
        Ok(())
    }
}
