//! `SaveLine` (`fem/src/modules/SaveData/SaveLine.F90`) for the existing-
//! boundary mode as a scheduled [`Solver`].
//!
//! Every boundary element whose boundary condition has the `Save Mask`
//! keyword (default `Save Line`) true contributes its nodes. The nodes are
//! numbered by first touch and, unless `Optimize Node Ordering = False`, in
//! Elmer's bandwidth-optimal order for the line (MeshBasics.F90
//! `MakePermUsingMask`), so a boundary curve comes out as an ordered
//! polyline. One row per node: the label columns, the three coordinates and
//! the requested variables.
//!
//! `Polyline Coordinates` with `Polyline Divisions` add points along each
//! line (`SavePolyLines`): a point is located in a bulk element with
//! `PointInElement` and the fields are interpolated there. Those rows are
//! tabulated and written sorted by their position along the lines, after
//! the others. `Show Norm Index` sums one column into the solver's norm.
//! Polylines without divisions (edge intersections), circles, isocurves,
//! axis lines, heat flux columns and intersections with element edges are
//! refused when the solver is built.

use elmer_assembly::{rcm_order, ListMatrix};
use elmer_core::{normalize_key, ElemId, Params};
use elmer_interp::geom::{point_in_element, PointTol};
use elmer_interp::shape::Shape;
use elmer_saveio::savedata::{append_text, line_names_text, line_row_text, LineColumns, LineNodeCounts, NamesMeta};
use elmer_solver::{tangent_directions, Solver, SolverCtx, SolverError, Variable};

use crate::plan::LinePlan;
use crate::scalars::{set_variable_norm, split_output_path};
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

/// `GetElementCorners`.
fn corners(code: u16) -> usize {
    match code / 100 {
        1 => 1,
        2 => 2,
        3 => 3,
        4 | 5 => 4,
        6 => 5,
        7 => 6,
        _ => 8,
    }
}

/// One `SavePolyLines` point: its position along the lines, labels, values.
struct PolyRow {
    pos: f64,
    labels: Vec<i64>,
    values: Vec<f64>,
}

/// `SavePolyLines` with `Polyline Divisions` (SaveLine.F90:1608): rows in
/// visiting order, and the number of points (`SaveNodes(2)`). Index loops as
/// in the Fortran.
#[allow(clippy::needless_range_loop)]
fn polyline_rows(
    ctx: &SolverCtx<'_>,
    plan: &LinePlan,
    cols: &[Column<'_>],
    label: &dyn Fn(usize, usize) -> Vec<i64>,
) -> Result<(Vec<PolyRow>, usize), SolverError> {
    let dim = ctx.mesh.dim();
    let eps = 1.0e-5;
    let mask_perm = ctx.primary().and_then(|f| f.perm.map(<[u32]>::to_vec));
    for c in cols {
        if c.view.elemental || c.view.dg || c.view.n_perm() > ctx.mesh.n_nodes() {
            return Err(SolverError::Unsupported(format!(
                "{CALLER}: polyline values of the non-nodal field '{}'",
                c.name
            )));
        }
    }
    let mut rows = Vec::new();
    let mut count = 0usize;
    let mut basis = vec![0.0; 27];
    for (li, &(a, b)) in plan.polylines.iter().enumerate() {
        let line = li + 1;
        let mut r0 = a;
        let mut r1 = b;
        if dim != 3 || plan.polyline_cols < 3 {
            r0[2] = 0.0;
            r1[2] = 0.0;
        }
        let dr = [r1[0] - r0[0], r1[1] - r0[1], r1[2] - r0[2]];
        let s = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
        let ln = dr.map(|x| x / s);
        let (t1, t2) = tangent_directions(ln);
        let nsize = plan.divisions[li].max(0) as usize;
        let mut tag = vec![false; nsize + 1];
        for e in (0..ctx.mesh.n_bulk() as u32).map(ElemId) {
            let el = ctx.mesh.element(e);
            if let Some(p) = &mask_perm {
                if el.nodes.iter().any(|&k| p.get(k as usize).is_none_or(|&q| q == elmer_core::NO_DOF)) {
                    continue;
                }
            }
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let (mut mn, mut mx) = ([0.0; 3], [0.0; 3]);
            for (i, x) in xyz.iter().take(corners(el.etype.0)).enumerate() {
                let ds = [(x[0] - r0[0]) / s, (x[1] - r0[1]) / s, (x[2] - r0[2]) / s];
                let dot = |v: &[f64; 3]| ds[0] * v[0] + ds[1] * v[1] + ds[2] * v[2];
                let lc = [dot(&ln), dot(&t1), dot(&t2)];
                if i == 0 {
                    (mn, mx) = (lc, lc);
                } else {
                    for j in 0..3 {
                        mn[j] = mn[j].min(lc[j]);
                        mx[j] = mx[j].max(lc[j]);
                    }
                }
            }
            if mn[1] * mx[1] > eps || (dim == 3 && mn[2] * mx[2] > eps) {
                continue;
            }
            let imin = ((nsize as f64 * mn[0]).ceil()).max(0.0) as usize;
            let fmax = (nsize as f64 * mx[0]).floor();
            if fmax < 0.0 {
                continue;
            }
            let imax = (fmax as usize).min(nsize);
            let shape = Shape::get(el.etype.0)
                .ok_or_else(|| SolverError::Unsupported(format!("{CALLER}: polylines through element type {}", el.etype.0)))?;
            for i in imin..=imax {
                if tag[i] {
                    continue;
                }
                let g: [f64; 3] = std::array::from_fn(|d| r0[d] + i as f64 * dr[d] / nsize as f64);
                let tol = PointTol { local: Some(eps), ..PointTol::default() };
                let probe = point_in_element(shape, &xyz, g, tol, false, dim);
                if !probe.inside {
                    continue;
                }
                let n = xyz.len();
                basis.resize(n, 0.0);
                shape.basis(probe.uvw, &mut basis[..n]);
                tag[i] = true;
                count += 1;
                let pos = i as f64 / nsize as f64 + 2.0 * (line - 1) as f64;
                // EvaluateVariableAtGivenPoint: SUM( Basis * Values ) over the
                // element's nodes, 0 where some node has no value.
                let mut values = Vec::with_capacity(3 + cols.len());
                for d in 0..3 {
                    let mut acc = 0.0;
                    for k in 0..n {
                        acc += basis[k] * xyz[k][d];
                    }
                    values.push(acc);
                }
                for c in cols {
                    let slots: Option<Vec<usize>> = el.nodes.iter().map(|&k| c.view.slot(k as usize)).collect();
                    values.push(match slots {
                        Some(sl) => {
                            let mut acc = 0.0;
                            for k in 0..n {
                                acc += basis[k] * c.view.at(sl[k], 0);
                            }
                            acc
                        }
                        None => 0.0,
                    });
                }
                rows.push(PolyRow { pos, labels: label(line, i), values });
            }
        }
    }
    Ok((rows, count))
}

impl Solver for SaveLine {
    /// `SaveLine_init`: `Show Norm Index` gets the norm variable
    /// `<Equation>_var` unless the section names a `Variable`.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if self.plan.norm_index <= 0 || self.plan.section.get("variable").is_some() {
            return Ok(());
        }
        let eq = self.plan.equation.clone().unwrap_or_default();
        let name = normalize_key(&format!("{eq}_var"));
        let id = match ctx.sim.vars.id(&name) {
            Some(id) => id,
            None => {
                let mut v = Variable::global(&name, 0.0);
                v.output = false;
                ctx.sim.vars.add(v)?
            }
        };
        let ix = ctx.sim.slots.iter().position(|s| s.id == ctx.solver_id()).expect("running solver has a slot");
        ctx.sim.slots[ix].variable = Some(id);
        Ok(())
    }

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
        let line = self.line.as_ref().expect("set above").clone();
        let cols = columns(ctx, &plan)?;
        // Coordinates plus the variables.
        let n_results = 3 + cols.len();
        let set_norm = |ctx: &mut SolverCtx<'_>, norm: f64| {
            if plan.norm_index > 0 {
                if let Some(id) = ctx.sim.slots.iter().find(|s| s.id == ctx.solver_id()).and_then(|s| s.variable) {
                    set_variable_norm(ctx, id, norm);
                }
            }
        };
        if cols.is_empty() {
            set_norm(ctx, 0.0);
            ctx.info(CALLER, "All done", 3);
            return Ok(());
        }
        let label = |bc: usize, node: usize| -> Vec<i64> {
            let mut labels: Vec<i64> = Vec::new();
            if !plan.skip_boundary_info {
                if plan.line_marker != 0 {
                    labels.push(plan.line_marker.into());
                }
                if transient {
                    labels.push(done_time as i64);
                }
                labels.push(visited as i64 + 1);
                labels.push(bc as i64);
                labels.push(node as i64);
            }
            labels
        };
        // WriteFieldsAtElement's norm: a label column, or a value column.
        let mut norm = 0.0;
        let ni = plan.norm_index.max(0) as usize;
        let mut add_norm = |labels: &[i64], values: &[f64]| {
            if ni == 0 {
                return;
            }
            if ni <= labels.len() {
                norm += labels[ni - 1] as f64;
            } else if let Some(v) = values.get(ni - labels.len() - 1) {
                norm += v;
            }
        };

        let (display, path) = split_output_path(ctx, plan.section.string("output directory"), &plan.filename, visited == 0)?;

        let mut text = String::new();
        let mut found_nan = 0usize;
        for (idx, &node) in line.nodes.iter().enumerate() {
            let labels = label(line.bc[idx], node + 1);
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
            add_norm(&labels, &values);
            text.push_str(&line_row_text(&labels, &values));
        }
        // SavePolyLines. A point at position 0 is written at once; the rest
        // are tabulated and written sorted by position (CloseLineFile).
        let (mut rows, poly_count) = if plan.polylines.is_empty() {
            (Vec::new(), 0)
        } else {
            polyline_rows(ctx, &plan, &cols, &label)?
        };
        let mut table = Vec::new();
        for mut r in rows.drain(..) {
            for v in &mut r.values {
                if v.is_nan() {
                    found_nan += 1;
                    *v = -1.0;
                }
            }
            add_norm(&r.labels, &r.values);
            if r.pos > 0.0 {
                table.push(r);
            } else {
                text.push_str(&line_row_text(&r.labels, &r.values));
            }
        }
        table.sort_by(|a, b| a.pos.total_cmp(&b.pos));
        let mut prev: Option<f64> = None;
        for r in &table {
            if let Some(p) = prev {
                if r.pos > 0.0 && (r.pos - p).abs() < f64::EPSILON {
                    continue;
                }
            }
            prev = Some(r.pos);
            text.push_str(&line_row_text(&r.labels, &r.values));
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
                polyline: poly_count,
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
        if plan.norm_index > 0 {
            set_norm(ctx, norm);
        }
        if found_nan > 0 {
            ctx.warn(CALLER, &format!("Replaced {found_nan} NaN entries with -1"));
        }
        ctx.info(CALLER, "All done", 3);
        Ok(())
    }
}
