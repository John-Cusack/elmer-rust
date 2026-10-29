//! `StructuredProjectToPlane` (StructuredProjectToPlane.F90:84): reduces
//! fields of a mesh that could have been extruded along a coordinate axis to
//! its top (or bottom, or middle) plane with an operator (`sum`, `min`,
//! `max`, `int`, `int mean`, `top`, `bottom`, `middle`, `isosurface`,
//! `layer below top`, `layer above bottom`, `thickness`), or computes full
//! dimensional column fields (`index`, `depth`, `height`, `distance`).
//! `upper <op>` and `lower <op>` restrict the operator to one side of the
//! `Mid Surface`.
//!
//! Not ported, and refused: the mask keywords (`Mapping Mask Variable`,
//! `Mapping Mask Name`), `Mapping Coordinate Transformation`, the `error ...`
//! operators and discontinuous Galerkin input variables.

use elmer_core::{FieldKind, ParamMap, Params, NO_DOF};
use elmer_solver::norm::{compute_norm, NormSpec};
use elmer_solver::{Solver, SolverCtx, SolverError, Variable};

use crate::structure::{detect_extruded_structure, Structure, Want};
use crate::util::{coords_of, slot_index};

const CALLER: &str = "StructuredProjectToPlane";

/// `-nooutput -global savescalars_var`, added by `_init` for `Show Norm Index`.
const NORM_VARIABLE: &str = "savescalars_var";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(StructuredProjectToPlane::default()))
}

fn fatal(msg: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, msg)
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

/// Layer bookkeeping built on the first visit.
struct Init {
    s: Structure,
    /// Slot of a top / bottom / middle node in a projected field, `NO_DOF` elsewhere.
    top_perm: Vec<u32>,
    bot_perm: Vec<u32>,
    mid_perm: Option<Vec<u32>>,
    top_nodes: usize,
}

#[derive(Default)]
pub struct StructuredProjectToPlane {
    init: Option<Init>,
}

/// A field read by an operator (`FieldIn`, `DofsIn`, `PermIn`).
#[derive(Clone)]
struct Input {
    values: Vec<f64>,
    dofs: usize,
    perm: Option<Vec<u32>>,
}

impl Input {
    /// `FieldIn(Dofs*(PermIn(i)-1)+dof)`; `None` where the field has no dof.
    fn at(&self, dofs: usize, i: usize, dof: usize) -> Option<f64> {
        let k = match &self.perm {
            Some(p) => {
                let k = p[i];
                if k == NO_DOF {
                    return None;
                }
                k as usize
            }
            None => i,
        };
        self.values.get(dofs * k + dof).copied()
    }
}

/// `VariableGet(Model % Variables, name)`; `Coordinate k` is the mesh.
fn input_variable(ctx: &SolverCtx<'_>, name: &str) -> Result<Input, SolverError> {
    let key = elmer_core::normalize_key(name);
    if let Some(k) = key
        .strip_prefix("coordinate ")
        .and_then(|k| k.parse::<usize>().ok())
        .filter(|k| (1..=3).contains(k))
    {
        let values = coords_of(ctx.mesh).iter().map(|c| c[k - 1]).collect();
        return Ok(Input {
            values,
            dofs: 1,
            perm: None,
        });
    }
    let v = ctx
        .sim
        .vars
        .get(name)
        .ok_or_else(|| fatal(format!("Variable does not exist: {name}")))?;
    if v.kind != FieldKind::Nodal {
        return Err(unsupported("input variables that are not nodal"));
    }
    Ok(Input {
        values: v.values.clone(),
        dofs: v.dofs,
        perm: v.perm.as_ref().map(|p| p.to_vec()),
    })
}

impl Solver for StructuredProjectToPlane {
    /// `StructuredProjectToPlane_init`: a variable for `Show Norm Index`.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let ix = slot_index(ctx);
        let wanted = ctx
            .params()
            .integer("show norm index")
            .is_some_and(|n| n > 0);
        if wanted && ctx.params().get("variable").is_none() && ctx.sim.slots[ix].variable.is_none()
        {
            let mut v = Variable::new(NORM_VARIABLE, FieldKind::Global, 1, None, vec![0.0]);
            v.output = false;
            let id = ctx.sim.vars.add(v)?;
            ctx.sim.slots[ix].variable = Some(id);
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "------------------------------------------", 4);
        ctx.info(CALLER, "Performing projection on a structured mesh ", 4);
        ctx.info(CALLER, "------------------------------------------", 4);
        let params = ctx.params().clone();
        for k in [
            "mapping mask variable",
            "mapping mask name",
            "mapping coordinate transformation",
        ] {
            if params.get(k).is_some() {
                return Err(unsupported(k));
            }
        }
        let nnodes = ctx.mesh.n_nodes();
        if self.init.is_none() {
            self.init = Some(initialize(ctx, &params)?);
        }
        let init = self.init.take().expect("initialised above");
        let r = project(ctx, &params, &init, nnodes);
        self.init = Some(init);
        r
    }
}

/// The first visit: structure and the top / bottom / middle numbering.
fn initialize(ctx: &SolverCtx<'_>, params: &ParamMap) -> Result<Init, SolverError> {
    let active = params
        .integer("active coordinate")
        .ok_or_else(|| fatal("Active Coordinate is missing"))?;
    if !(1..=3).contains(&active) {
        return Err(fatal("Invalid value for Active Coordinate"));
    }
    let mut unit = [0.0; 3];
    unit[active as usize - 1] = if params.logical("project to bottom").unwrap_or(false) {
        -1.0
    } else {
        1.0
    };
    let eps = ctx
        .creal(params, "dot product tolerance")?
        .unwrap_or(1.0e-4);
    let s = detect_extruded_structure(
        ctx.mesh,
        &ctx.sim.model,
        unit,
        eps,
        Want {
            layers: false,
            mid: true,
        },
    )?;
    let nnodes = ctx.mesh.n_nodes();

    let number = |pointer: &[usize]| -> (Vec<u32>, usize) {
        let mut perm = vec![NO_DOF; nnodes];
        let mut n = 0;
        for i in 0..nnodes {
            if pointer[i] == i {
                perm[i] = n as u32;
                n += 1;
            }
        }
        (perm, n)
    };
    let (top_perm, top_nodes) = number(&s.top);
    ctx.info(CALLER, &format!("Number of top nodes: {top_nodes}"), 10);
    let (bot_perm, bot_nodes) = number(&s.bot);
    ctx.info(CALLER, &format!("Number of bot nodes: {bot_nodes}"), 10);
    if bot_nodes == 0 {
        return Err(fatal("Cannot continue with zero BotNodes!"));
    }
    let no_layers = nnodes / bot_nodes;
    if no_layers < 2 {
        return Err(fatal(format!(
            "Solver does not makse sense with {no_layers} layers!"
        )));
    }
    if bot_nodes != top_nodes {
        ctx.warn(
            CALLER,
            &format!("Conflicting BotNodes vs. TopNodes: {bot_nodes} - {top_nodes}"),
        );
    }
    let mid_perm = s.mid.as_ref().map(|m| {
        let (perm, n) = number(m);
        if n != bot_nodes {
            ctx.warn(
                CALLER,
                &format!("Conflicting BotNodes vs. MidNodes: {bot_nodes} - {n}"),
            );
        }
        perm
    });
    Ok(Init {
        s,
        top_perm,
        bot_perm,
        mid_perm,
        top_nodes,
    })
}

fn project(
    ctx: &mut SolverCtx<'_>,
    params: &ParamMap,
    init: &Init,
    nnodes: usize,
) -> Result<(), SolverError> {
    let s = &init.s;
    let coord_dir = params
        .integer("active coordinate")
        .expect("checked in initialize") as usize
        - 1;
    let coord: Vec<f64> = coords_of(ctx.mesh).iter().map(|c| c[coord_dir]).collect();
    let norm_ind = params.integer("show norm index").unwrap_or(0);
    let mid_exists = s.mid.is_some();

    let mut got_old_var = false;
    let mut old_var_name = String::new();
    let mut old_oper = String::new();
    let mut field_in: Option<Input> = None;

    let mut no_var = 0usize;
    loop {
        no_var += 1;
        let var_name = params
            .string(&format!("variable {no_var}"))
            .map(str::to_string);
        let mut got_var = var_name.is_some();
        let mut var_name = var_name.unwrap_or_default();
        let mut dofs_in = 1;
        let mut dofs = 1;
        if got_var {
            let inp = input_variable(ctx, &var_name)?;
            dofs_in = inp.dofs;
            dofs = inp.dofs;
            field_in = Some(inp);
            got_old_var = true;
            old_var_name = var_name.clone();
        }

        let oper_key = params
            .string(&format!("operator {no_var}"))
            .map(|s| s.trim().to_ascii_lowercase());
        let got_oper = oper_key.is_some();
        let oper = match oper_key {
            Some(o) => o,
            None if no_var > 1 => old_oper.clone(),
            None => String::new(),
        };
        // Either a new field or a new operator is needed.
        if !(got_var || got_oper) {
            if no_var == 1 {
                ctx.warn(CALLER, "Not even one field and operator to treat?");
            }
            break;
        }
        old_oper = oper.clone();
        if got_oper && !got_var {
            got_var = got_old_var;
            if got_var {
                var_name = old_var_name.clone();
            }
        }

        let oper0 = oper.clone();
        let (mut upper, mut lower) = (false, false);
        let mut oper = oper;
        if let Some(rest) = oper0.strip_prefix("upper ") {
            if !mid_exists {
                return Err(fatal("Upper operator cannot exist without midlayer"));
            }
            upper = true;
            oper = rest.trim().to_string();
        } else if let Some(rest) = oper0.strip_prefix("lower ") {
            if !mid_exists {
                return Err(fatal("Lower operator cannot exist without midlayer"));
            }
            lower = true;
            oper = rest.trim().to_string();
        }
        let tmp_top: &[usize] = if lower {
            s.mid.as_deref().expect("checked")
        } else {
            &s.top
        };
        let tmp_bot: &[usize] = if upper {
            s.mid.as_deref().expect("checked")
        } else {
            &s.bot
        };

        let column_op = matches!(
            oper.as_str(),
            "height" | "depth" | "index" | "thickness" | "distance"
        );
        if !column_op && !got_var {
            return Err(fatal(format!(
                "Variable required for this operator: {oper}"
            )));
        }
        if oper.starts_with("error") {
            return Err(unsupported("error operators"));
        }

        let target_name = params
            .string(&format!("target variable {no_var}"))
            .map(str::to_string)
            .unwrap_or_else(|| format!("{oper0} {var_name}"));
        let rdofs = dofs;
        let reduced = !matches!(oper.as_str(), "height" | "depth" | "index" | "distance");
        let everywhere = params
            .logical("project to everywhere")
            .or_else(|| params.logical(&format!("target variable {no_var} everywhere")))
            .unwrap_or(false);

        // The projected variable, created if needed.
        if ctx.sim.vars.get(&target_name).is_none() {
            let perm: Vec<u32> = if reduced && !everywhere {
                if params
                    .logical(&format!("target variable {no_var} at bottom"))
                    .unwrap_or(false)
                {
                    init.bot_perm.clone()
                } else if mid_exists
                    && params
                        .logical(&format!("target variable {no_var} at middle"))
                        .unwrap_or(false)
                {
                    init.mid_perm.clone().expect("mid exists")
                } else {
                    init.top_perm.clone()
                }
            } else {
                (0..nnodes as u32).collect()
            };
            let slots = perm.iter().filter(|&&p| p != NO_DOF).count();
            let v = Variable::new(
                &target_name,
                FieldKind::Nodal,
                rdofs,
                Some(perm.into()),
                vec![0.0; slots * rdofs],
            );
            ctx.sim.vars.add(v)?;
            ctx.info(CALLER, &format!("Created variable: {target_name}"), 9);
        }
        let (out_dofs, out_perm) = {
            let v = ctx.sim.vars.get(&target_name).expect("exists");
            (v.dofs, v.perm.as_ref().map(|p| p.to_vec()))
        };
        if out_dofs != rdofs {
            return Err(fatal("Mismatch in the dofs in fields!"));
        }
        let mut field_out = vec![0.0; ctx.sim.vars.get(&target_name).expect("exists").values.len()];

        let mut levelset: Option<(Input, f64)> = None;
        let mut layer = 0usize;
        if oper == "isosurface" {
            let inp = match params.string(&format!("isosurface variable {no_var}")) {
                Some(n) => input_variable(ctx, n)?,
                None => Input {
                    values: coord.clone(),
                    dofs: 1,
                    perm: None,
                },
            };
            let level = params
                .real(&format!("isosurface value {no_var}"))
                .unwrap_or(0.0);
            levelset = Some((inp, level));
        } else if oper.starts_with("layer") {
            layer = params
                .integer(&format!("layer index {no_var}"))
                .ok_or_else(|| fatal(format!("no \"Layer Index {no_var}\" indicated for output")))?
                as usize;
        }

        let skip = |j: usize| -> bool {
            let m = s.mid.as_ref();
            if upper {
                coord[j] < coord[m.expect("checked")[j]]
            } else if lower {
                coord[j] > coord[m.expect("checked")[j]]
            } else {
                false
            }
        };
        let fin = field_in.clone();
        let inp =
            |i: usize, dof: usize| -> Option<f64> { fin.as_ref().and_then(|f| f.at(dofs, i, dof)) };
        let mut top_field = vec![0.0; init.top_nodes];
        let tp = |i: usize| init.top_perm[i] as usize;

        for dof in 0..rdofs {
            match oper.as_str() {
                "sum" | "min" | "max" => {
                    top_field.iter_mut().for_each(|v| {
                        *v = match oper.as_str() {
                            "sum" => 0.0,
                            "min" => f64::MAX,
                            _ => -f64::MAX,
                        }
                    });
                    for i in 0..nnodes {
                        if skip(i) {
                            continue;
                        }
                        let itop = s.top[i];
                        let Some(x) = inp(i, dof) else { continue };
                        let t = &mut top_field[tp(itop)];
                        *t = match oper.as_str() {
                            "sum" => *t + x,
                            "min" => x.min(*t),
                            _ => x.max(*t),
                        };
                    }
                }
                "bottom" => {
                    top_field.iter_mut().for_each(|v| *v = 0.0);
                    for i in 0..nnodes {
                        if skip(i) {
                            continue;
                        }
                        if i == tmp_bot[i] {
                            if let Some(x) = inp(i, dof) {
                                top_field[tp(s.top[i])] = x;
                            }
                        }
                    }
                }
                "top" => {
                    top_field.iter_mut().for_each(|v| *v = 0.0);
                    for i in 0..nnodes {
                        if i == tmp_top[i] {
                            if let Some(x) = inp(i, dof) {
                                top_field[tp(s.top[i])] = x;
                            }
                        }
                    }
                }
                "middle" => {
                    top_field.iter_mut().for_each(|v| *v = 0.0);
                    let mid = s
                        .mid
                        .as_ref()
                        .ok_or_else(|| fatal("Middle operator cannot exist without midlayer"))?;
                    for i in 0..nnodes {
                        if i == mid[i] {
                            if let Some(x) = inp(i, dof) {
                                top_field[tp(s.top[i])] = x;
                            }
                        }
                    }
                }
                "layer below top" => {
                    top_field.iter_mut().for_each(|v| *v = 0.0);
                    for i in 0..nnodes {
                        if i == tmp_top[i] {
                            let mut l = i;
                            for _ in 0..layer {
                                l = s.down[l];
                            }
                            if let Some(x) = inp(l, dof) {
                                top_field[tp(i)] = x;
                            }
                        }
                    }
                }
                "layer above bottom" => {
                    top_field.iter_mut().for_each(|v| *v = 0.0);
                    for i in 0..nnodes {
                        if skip(i) {
                            continue;
                        }
                        if i == tmp_bot[i] {
                            let mut l = i;
                            for _ in 0..layer {
                                l = s.up[l];
                            }
                            if let Some(x) = inp(l, dof) {
                                top_field[tp(s.top[i])] = x;
                            }
                        }
                    }
                }
                "isosurface" => {
                    top_field.iter_mut().for_each(|v| *v = 0.0);
                    let (ls, level) = levelset.as_ref().expect("set above");
                    for i in 0..nnodes {
                        let iup = s.up[i];
                        let (j, jup) = match &ls.perm {
                            Some(p) => {
                                if p[i] == NO_DOF || p[iup] == NO_DOF {
                                    continue;
                                }
                                (p[i] as usize, p[iup] as usize)
                            }
                            None => (i, iup),
                        };
                        if (ls.values[jup] - level) * (ls.values[j] - level) <= 0.0 {
                            let itop = s.top[i];
                            let dx = (ls.values[jup] - ls.values[j]).abs();
                            if dx.abs() < f64::EPSILON {
                                return Err(fatal("dx smaller than machine Epsilon!"));
                            }
                            let q = (ls.values[jup] - level).abs() / dx;
                            if let (Some(a), Some(b)) = (inp(i, dof), inp(iup, dof)) {
                                top_field[tp(itop)] = q * a + (1.0 - q) * b;
                            }
                        }
                    }
                }
                "int" | "int mean" => {
                    top_field.iter_mut().for_each(|v| *v = 0.0);
                    for i in 0..nnodes {
                        if skip(i) {
                            continue;
                        }
                        let height = if oper == "int mean" {
                            (coord[s.top[i]] - coord[s.bot[i]]).abs()
                        } else {
                            1.0
                        };
                        // For top and bottom this reduces the distance to half.
                        let iup = if i == tmp_top[i] { i } else { s.up[i] };
                        let idown = if i == tmp_bot[i] { i } else { s.down[i] };
                        let mut dx = 0.5 * (coord[iup] - coord[idown]);
                        if height.abs() < f64::EPSILON {
                            return Err(fatal("height smaller than machine Epsilon!"));
                        }
                        dx = dx.abs() / height;
                        let Some(x) = inp(i, dof) else { continue };
                        top_field[tp(s.top[i])] += dx * x;
                    }
                }
                "thickness" => {
                    top_field.iter_mut().for_each(|v| *v = 0.0);
                    for i in 0..nnodes {
                        if skip(i) {
                            continue;
                        }
                        let iup = if i == tmp_top[i] { i } else { s.up[i] };
                        let idown = if i == tmp_bot[i] { i } else { s.down[i] };
                        let dx = (0.5 * (coord[iup] - coord[idown])).abs();
                        top_field[tp(s.top[i])] += dx;
                    }
                }
                "index" | "depth" | "height" | "distance" => {
                    let slot = |l: usize| -> Option<usize> {
                        match &out_perm {
                            Some(p) => (p[l] != NO_DOF).then_some(p[l] as usize),
                            None => Some(l),
                        }
                    };
                    let nsize = nnodes;
                    // Walks a column from its top downwards: `index` and `depth`.
                    let from_top = |field_out: &mut Vec<f64>, index: bool| {
                        for i in 0..nnodes {
                            if skip(i) || i != tmp_top[i] {
                                continue;
                            }
                            let mut l = i;
                            let mut depth = 0.0;
                            for k in 1..=nsize {
                                if k > 1 {
                                    depth += (coord[s.up[l]] - coord[l]).abs();
                                }
                                if let Some(p) = slot(l) {
                                    field_out[p] = if index { k as f64 } else { depth };
                                }
                                if l == tmp_bot[l] {
                                    break;
                                }
                                l = s.down[l];
                            }
                        }
                    };
                    // Walks a column from its bottom upwards: `height`.
                    let from_bottom = |field_out: &mut Vec<f64>, min: bool| {
                        for i in 0..nsize {
                            if skip(i) || i != tmp_bot[i] {
                                continue;
                            }
                            let mut l = i;
                            let mut height = 0.0;
                            for k in 1..=nsize {
                                if k > 1 {
                                    height += (coord[l] - coord[s.down[l]]).abs();
                                }
                                if let Some(p) = slot(l) {
                                    field_out[p] = if min {
                                        height.min(field_out[p])
                                    } else {
                                        height
                                    };
                                }
                                if l == tmp_top[l] {
                                    break;
                                }
                                l = s.up[l];
                            }
                        }
                    };
                    field_out.iter_mut().for_each(|v| *v = 0.0);
                    match oper.as_str() {
                        "index" => from_top(&mut field_out, true),
                        "depth" => from_top(&mut field_out, false),
                        "height" => from_bottom(&mut field_out, false),
                        _ => {
                            from_top(&mut field_out, false);
                            from_bottom(&mut field_out, true);
                        }
                    }
                }
                other => return Err(fatal(format!("Unknown operator: {other}"))),
            }

            // Copy the projected values to the target variable.
            if reduced {
                for i in 0..nnodes {
                    let j = match &out_perm {
                        Some(p) => {
                            if p[i] == NO_DOF {
                                continue;
                            }
                            p[i] as usize
                        }
                        None => i,
                    };
                    field_out[rdofs * j + dof] = top_field[tp(s.top[i])];
                }
            }
        }

        let _ = dofs_in;
        let v = ctx.sim.vars.get_mut(&target_name).expect("exists");
        v.values.copy_from_slice(&field_out);
        v.values_changed = true;

        if norm_ind as usize == no_var {
            let spec = NormSpec::from_params(ctx.params(), 1);
            let norm = compute_norm(&spec, &field_out, field_out.len(), ctx.comm);
            let ix = slot_index(ctx);
            if let Some(id) = ctx.sim.slots[ix].variable {
                ctx.sim
                    .vars
                    .by_id_mut(id)
                    .values
                    .iter_mut()
                    .for_each(|v| *v = norm);
            }
        }
    }
    ctx.info(CALLER, "------------------------------------------", 4);
    Ok(())
}
