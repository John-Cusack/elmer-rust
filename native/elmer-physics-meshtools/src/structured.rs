//! `StructuredMeshMapper` (StructuredMeshMapper.F90:76): maps a mesh that
//! could have been extruded along one coordinate axis between a bottom and a
//! top surface (and optionally a mid surface), or between several fixed
//! layers, by linear interpolation along each column of nodes.
//!
//! Not ported, and refused: `Mapping Mask Variable` / `Mapping Mask Name`,
//! `Mapping Coordinate Transformation`, `Mapping Original Coordinates`,
//! `Base Displacement Variable`, `Mapped Mesh Name`, `Mesh Height Map` and
//! `Correct Surface Mask`.

use elmer_core::{ElemId, FieldKind, MeshView, ParamMap, Params, Value};
use elmer_solver::{Solver, SolverCtx, SolverError, Variable};

use crate::structure::{detect_extruded_structure, Structure, Want};
use crate::util::{apply_coords, coords_of};

const CALLER: &str = "StructuredMeshMapper";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(StructuredMeshMapper::default()))
}

#[derive(Default)]
pub struct StructuredMeshMapper {
    visited: bool,
    structure: Option<Structure>,
    /// `FixedLayers`, with the first and last entry forced to the end layers.
    fixed_layers: Vec<i32>,
    /// `Field`: top, bottom and mid surface values from the boundary
    /// conditions, kept between visits.
    field: Option<Vec<f64>>,
    /// `MinHeight` of the current visit.
    min_height: f64,
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

fn fatal(msg: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, msg)
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
enum Mode {
    /// Keep the coordinate (or displacement 0).
    None,
    Level,
    Variable,
    Bc,
}

/// A one-dof nodal variable the mapper writes, resolved by name, with its
/// node -> slot permutation (`VeloVar`, `UpdateVar`).
struct Target {
    name: String,
    perm: Option<Vec<u32>>,
}

impl Target {
    fn set(&self, ctx: &mut SolverCtx<'_>, node: usize, v: f64) {
        let slot = match &self.perm {
            Some(p) => {
                let s = p[node];
                if s == elmer_core::NO_DOF {
                    return;
                }
                s as usize
            }
            None => node,
        };
        if let Some(mut c) = ctx.component_mut(&self.name) {
            c.set(slot, v);
        }
    }
}

/// `VariableGet`, or `DefaultVariableAdd` when missing: a name ending in a
/// component number creates the whole vector variable.
fn velocity_or_update(
    ctx: &mut SolverCtx<'_>,
    key: &str,
    nnodes: usize,
    dim: usize,
) -> Result<Option<Target>, SolverError> {
    let Some(name) = ctx.params().string(key).map(str::to_string) else {
        return Ok(None);
    };
    let existing = ctx
        .component(&name)
        .map(|c| (c.len(), c.perm.map(<[u32]>::to_vec)));
    if let Some((len, perm)) = existing {
        if len != nnodes {
            return Err(fatal(format!("Sizes must agree: {name}")));
        }
        if ctx.sim.vars.get(&name).is_some_and(|v| v.dofs != 1) {
            return Err(fatal(format!("The size of {key} must be one")));
        }
        return Ok(Some(Target { name, perm }));
    }
    let perm: std::sync::Arc<[u32]> = (0..nnodes as u32).collect::<Vec<_>>().into();
    let t = name.trim_end();
    if t.ends_with(['1', '2', '3']) {
        let base = t[..t.len() - 1].trim();
        let v = Variable::new(
            base,
            FieldKind::Nodal,
            dim,
            Some(perm.clone()),
            vec![0.0; nnodes * dim],
        );
        ctx.sim.vars.add(v)?;
    } else {
        ctx.sim.vars.add(Variable::new(
            t,
            FieldKind::Nodal,
            1,
            Some(perm.clone()),
            vec![0.0; nnodes],
        ))?;
    }
    Ok(Some(Target {
        name,
        perm: Some(perm.to_vec()),
    }))
}

impl Solver for StructuredMeshMapper {
    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "---------------------------------------", 4);
        ctx.info(CALLER, "Performing mapping on a structured mesh ", 4);
        ctx.info(CALLER, "---------------------------------------", 4);

        let params = ctx.params().clone();
        for k in [
            "mapping mask variable",
            "mapping mask name",
            "mapping coordinate transformation",
            "mapping original coordinates",
            "base displacement variable",
            "mapped mesh name",
            "mesh height map",
            "correct surface mask",
        ] {
            if params.get(k).is_some() {
                return Err(unsupported(k));
            }
        }
        let dim = ctx.mesh.dim();
        let nnodes = ctx.mesh.n_nodes();
        let reinitialize = params.logical("always detect structure").unwrap_or(false);
        if reinitialize {
            self.field = None;
        }
        let fixed = match params.get("fixed layer indexes") {
            Some(Value::Integer(v)) => Some(v.clone()),
            _ => None,
        };
        let multi_layer = fixed.is_some();
        let bot_proj = params.logical("project to bottom").unwrap_or(false);
        let active = params
            .integer("active coordinate")
            .ok_or_else(|| fatal("Active Coordinate is missing"))?;
        if !(1..=3).contains(&active) {
            return Err(fatal("Invalid value for Active Coordinate"));
        }
        let dir = active as usize - 1;

        if self.structure.is_none() || reinitialize {
            let mut unit = [0.0; 3];
            unit[dir] = if bot_proj { -1.0 } else { 1.0 };
            let eps = ctx
                .creal(&params, "dot product tolerance")?
                .unwrap_or(1.0e-4);
            let want = Want {
                layers: multi_layer,
                mid: !multi_layer,
            };
            ctx.info(CALLER, "Determining extruded structure", 6);
            let mut s = detect_extruded_structure(ctx.mesh, &ctx.sim.model, unit, eps, want)?;
            if let Some(mut f) = fixed {
                s.layers += 1;
                if f[0] != 1 {
                    ctx.warn(
                        CALLER,
                        &format!("Enforcing first fixed layer to: 1 (was {})", f[0]),
                    );
                    f[0] = 1;
                }
                let last = f.len() - 1;
                if f[last] as usize != s.layers {
                    ctx.warn(
                        CALLER,
                        &format!(
                            "Enforcing last fixed layer to: {} (was {})",
                            s.layers, f[last]
                        ),
                    );
                    f[last] = s.layers as i32;
                }
                self.fixed_layers = f;
            }
            self.structure = Some(s);
        }
        let s = self.structure.take().expect("structure detected above");
        let r = self.map(
            ctx,
            &params,
            &s,
            dt,
            dir,
            dim,
            nnodes,
            multi_layer,
            bot_proj,
        );
        self.structure = Some(s);
        r
    }
}

impl StructuredMeshMapper {
    #[allow(clippy::too_many_arguments)]
    fn map(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        params: &ParamMap,
        s: &Structure,
        dt: f64,
        dir: usize,
        dim: usize,
        nnodes: usize,
        multi_layer: bool,
        bot_proj: bool,
    ) -> Result<(), SolverError> {
        let all = coords_of(ctx.mesh);
        let orig: Vec<f64> = all.iter().map(|c| c[dir]).collect();
        let mut coord = orig.clone();

        if params.get("base displacement first").is_some() {
            return Err(unsupported("Base Displacement First"));
        }
        let velo = velocity_or_update(ctx, "mesh velocity variable", nnodes, dim)?;
        let update = velocity_or_update(ctx, "mesh update variable", nnodes, dim)?;
        let displacement_mode = params.logical("displacement mode").unwrap_or(false);
        let min_height = match ctx.creal(params, "minimum mesh height")? {
            Some(h) => h,
            None => ctx.creal(params, "minimum height")?.unwrap_or(f64::EPSILON),
        };
        self.min_height = min_height;
        ctx.info(
            CALLER,
            &format!("Adjusting upper surface to maintain minimum height to:{min_height:11.4E}"),
            6,
        );

        let (mut tangled, mut limited) = (0usize, 0usize);
        if multi_layer {
            self.multi_layer_mapper(
                ctx,
                params,
                s,
                dt,
                &orig,
                &mut coord,
                displacement_mode,
                bot_proj,
                &velo,
                &update,
                &mut tangled,
                &mut limited,
            )?;
        } else {
            self.binary_layer_mapper(
                ctx,
                params,
                s,
                dt,
                &orig,
                &mut coord,
                displacement_mode,
                &velo,
                &update,
                &mut tangled,
                &mut limited,
            )?;
        }
        if limited > 0 {
            ctx.info(
                CALLER,
                &format!("There seems to be {limited} (out of {nnodes}) limited heights!"),
                6,
            );
        }
        if tangled > 0 {
            ctx.info(
                CALLER,
                &format!("There seems to be {tangled} (out of {nnodes}) tangled nodes!"),
                5,
            );
        }
        if params.logical("mesh mapping passive").unwrap_or(false) {
            ctx.info(CALLER, "Taking back the suggested mapping!", 5);
            coord.clone_from(&orig);
        }
        if let Some(v) = velo.as_ref().filter(|_| {
            !self.visited && params.logical("mesh velocity first zero").unwrap_or(false)
        }) {
            for n in 0..nnodes {
                v.set(ctx, n, 0.0);
            }
        }
        let mut new = all.clone();
        for (c, x) in new.iter_mut().zip(&coord) {
            c[dir] = *x;
        }
        apply_coords(CALLER, &all, &new)?;
        self.visited = true;
        Ok(())
    }

    /// The surface values of `key` from the boundary conditions, written
    /// into `Field` at the element nodes (later elements win).
    fn bc_surface(
        &mut self,
        ctx: &SolverCtx<'_>,
        key: &str,
        nnodes: usize,
    ) -> Result<(), SolverError> {
        let field = self.field.get_or_insert_with(|| vec![0.0; nnodes]);
        let elems: Vec<ElemId> = ctx.boundary_elements().collect();
        for e in elems {
            let bc = ctx.bc_params(e).expect("boundary element with a BC");
            let Some(vals) = ctx.real_at_nodes(bc, key, e)? else {
                continue;
            };
            for (&n, &x) in ctx.mesh.element(e).nodes.iter().zip(&vals) {
                field[n as usize] = x;
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn binary_layer_mapper(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        params: &ParamMap,
        s: &Structure,
        dt: f64,
        orig: &[f64],
        coord: &mut [f64],
        displacement_mode: bool,
        velo: &Option<Target>,
        update: &Option<Target>,
        tangled_count: &mut usize,
        limited_count: &mut usize,
    ) -> Result<(), SolverError> {
        let nnodes = orig.len();
        let detangle = params.logical("correct surface").unwrap_or(false);
        let min_height = self.min_height;
        if detangle {
            ctx.info(
                CALLER,
                "> Correct Surface < in case of intersecting upper and lower surface",
                4,
            );
        }

        let (top_mode, top_val0, top_var) = self.surface_mode(
            ctx,
            params,
            "top surface level",
            "top surface variable name",
            "top surface",
            nnodes,
        )?;
        if top_mode == Mode::None {
            return Err(fatal("Top surface BC entry is missing"));
        }
        let (bot_mode, bot_val0, bot_var) = self.surface_mode(
            ctx,
            params,
            "bottom surface level",
            "bottom surface variable name",
            "bottom surface",
            nnodes,
        )?;
        if let Some(mid) = &s.mid {
            let _ = mid;
            self.bc_surface(ctx, "mid surface", nnodes)?;
        }
        let field = self.field.clone();

        let var_value = |v: &Option<(Vec<f64>, Option<Vec<u32>>)>,
                         node: usize,
                         what: &str|
         -> Result<f64, SolverError> {
            let (vals, perm) = v.as_ref().expect("variable mode has a variable");
            match perm {
                Some(p) => {
                    let sl = p[node];
                    if sl == elmer_core::NO_DOF {
                        return Err(fatal(format!("{what} surface variable perm is zero!")));
                    }
                    Ok(vals[sl as usize])
                }
                None => Ok(vals[node]),
            }
        };

        for i in 0..nnodes {
            let (itop, ibot) = (s.top[i], s.bot[i]);
            let (x0top, x0bot, x0loc) = (orig[itop], orig[ibot], orig[i]);
            let fld = |n: usize| field.as_ref().map_or(0.0, |f| f[n]);

            let mut top_val = match top_mode {
                Mode::Level => top_val0,
                Mode::Variable => var_value(&top_var, itop, "Top")?,
                Mode::Bc => fld(itop),
                Mode::None => unreachable!("refused above"),
            };
            let bot_val = match bot_mode {
                Mode::Level => bot_val0,
                Mode::Variable => var_value(&bot_var, ibot, "Bottom")?,
                Mode::Bc => fld(ibot),
                Mode::None => {
                    if displacement_mode {
                        0.0
                    } else {
                        x0bot
                    }
                }
            };
            let (mid_val, x0mid) = match &s.mid {
                Some(m) => (fld(m[i]), orig[m[i]]),
                None => (bot_val, x0bot),
            };

            // Check whether the mesh gets tangled (only the upper part with a mid layer).
            let mut dx = if displacement_mode {
                top_val + x0top - mid_val - x0mid
            } else {
                top_val - mid_val
            };
            if dx < min_height {
                if dx < f64::MIN_POSITIVE {
                    *tangled_count += 1;
                }
                *limited_count += 1;
                if detangle {
                    dx = min_height;
                    top_val = if displacement_mode {
                        mid_val + x0mid + x0top + dx
                    } else {
                        mid_val + dx
                    };
                }
            }

            let xloc = if s.mid.is_some() {
                if (x0top - x0mid) * (x0loc - x0mid) > 0.0 {
                    let w = (x0loc - x0mid) / (x0top - x0mid);
                    w * top_val + (1.0 - w) * mid_val
                } else {
                    let w = (x0loc - x0bot) / (x0mid - x0bot);
                    w * mid_val + (1.0 - w) * bot_val
                }
            } else {
                let w = (x0loc - x0bot) / (x0top - x0bot);
                w * top_val + (1.0 - w) * bot_val
            };

            if displacement_mode {
                if let Some(v) = velo {
                    v.set(ctx, i, xloc / dt);
                }
                coord[i] = orig[i] + xloc;
            } else {
                if let Some(v) = velo {
                    v.set(ctx, i, (xloc - orig[i]) / dt);
                }
                coord[i] = xloc;
            }
            if let Some(u) = update {
                u.set(ctx, i, coord[i] - orig[i]);
            }
        }
        Ok(())
    }

    /// Level, variable or boundary condition source of one surface
    /// (StructuredMeshMapper.F90:395-433 for the top, 435-477 for the bottom).
    #[allow(clippy::type_complexity)]
    fn surface_mode(
        &mut self,
        ctx: &SolverCtx<'_>,
        params: &ParamMap,
        level: &str,
        var_key: &str,
        bc_key: &str,
        nnodes: usize,
    ) -> Result<(Mode, f64, Option<(Vec<f64>, Option<Vec<u32>>)>), SolverError> {
        if let Some(v) = ctx.creal(params, level)? {
            return Ok((Mode::Level, v, None));
        }
        if let Some(name) = params.string(var_key) {
            let Some(c) = ctx.component(name) else {
                return Err(fatal(format!(
                    "{} surface variable is missing: {name}",
                    if bc_key.starts_with("top") {
                        "Top"
                    } else {
                        "Bottom"
                    }
                )));
            };
            if ctx.sim.vars.get(name).is_some_and(|v| v.dofs != 1) {
                return Err(fatal(format!(
                    "{} surface variable should have only 1 dof",
                    if bc_key.starts_with("top") {
                        "Top"
                    } else {
                        "Bottom"
                    }
                )));
            }
            let vals: Vec<f64> = c.iter().collect();
            return Ok((
                Mode::Variable,
                0.0,
                Some((vals, c.perm.map(<[u32]>::to_vec))),
            ));
        }
        if ctx
            .sim
            .model
            .boundary_conditions
            .iter()
            .any(|b| b.get(bc_key).is_some())
        {
            self.bc_surface(ctx, bc_key, nnodes)?;
            return Ok((Mode::Bc, 0.0, None));
        }
        Ok((Mode::None, 0.0, None))
    }

    #[allow(clippy::too_many_arguments)]
    fn multi_layer_mapper(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        params: &ParamMap,
        s: &Structure,
        dt: f64,
        orig: &[f64],
        coord: &mut [f64],
        displacement_mode: bool,
        bot_proj: bool,
        velo: &Option<Target>,
        update: &Option<Target>,
        tangled_count: &mut usize,
        limited_count: &mut usize,
    ) -> Result<(), SolverError> {
        let nnodes = orig.len();
        let nfixed = self.fixed_layers.len();
        let nlayers = s.layers;
        ctx.info(CALLER, &format!("Mapping using {nfixed} fixed layers"), 6);
        let detangle = params.logical("correct surface").unwrap_or(false);
        let name = params
            .string("fixed layer variable")
            .ok_or_else(|| fatal("No > Fixed Layer Variable < found"))?
            .to_string();
        let (fixed_dofs, fixed_vals, fixed_perm) = {
            let v = ctx
                .sim
                .vars
                .get(&name)
                .ok_or_else(|| fatal(format!("Could not find variable: {name}")))?;
            (
                v.dofs,
                v.values.clone(),
                v.perm.as_ref().map(|p| p.to_vec()),
            )
        };
        if fixed_dofs != nfixed {
            return Err(fatal(format!(
                "Invalid number of components in fixed layer variable:{fixed_dofs}"
            )));
        }
        let min_height = self.min_height;
        let fl = &self.fixed_layers;

        for i in 0..nnodes {
            let ibot = s.bot[i];
            // Start mapping from the bottom.
            if ibot != i {
                continue;
            }
            let mut stride_ind = vec![i; nlayers];
            let mut stride_coord = vec![orig[i]; nlayers];
            let mut j = i;
            let mut itop = i;
            for k in 1..nlayers {
                j = s.up[j];
                stride_coord[k] = orig[j];
                stride_ind[k] = j;
                itop = j;
            }

            // The projection matrix of this column.
            let mut proj = vec![vec![0.0f64; nfixed]; nlayers];
            for k in 1..=nlayers {
                let mut hit = false;
                for j in 1..=nfixed {
                    let fj = fl[j - 1] as usize;
                    if fj == k {
                        proj[k - 1][j - 1] = 1.0;
                        hit = true;
                    } else if j < nfixed && fj < k && (fl[j] as usize) > k {
                        let q = (stride_coord[k - 1] - stride_coord[fj - 1])
                            / (stride_coord[fl[j] as usize - 1] - stride_coord[fj - 1]);
                        proj[k - 1][j] = q;
                        proj[k - 1][j - 1] = 1.0 - q;
                        hit = true;
                    }
                    if hit {
                        break;
                    }
                }
                if !hit {
                    return Err(fatal(format!("Could not find mapping for layer: {k}")));
                }
            }

            // The fixed layer variable is at the bottom or the top of the column.
            let slot = |n: usize| match &fixed_perm {
                Some(p) => (p[n] != elmer_core::NO_DOF).then_some(p[n] as usize),
                None => Some(n),
            };
            let k = slot(ibot)
                .or_else(|| slot(itop))
                .ok_or_else(|| fatal("Could not find fixed fields at top or bottom!"))?;
            let fixed_coord = &fixed_vals[nfixed * k..nfixed * (k + 1)];

            for (l, row) in proj.iter().enumerate() {
                let p = row.iter().zip(fixed_coord).fold(0.0, |a, (x, y)| a + x * y);
                stride_coord[l] = if displacement_mode {
                    stride_coord[l] + p
                } else {
                    p
                };
            }

            if detangle {
                let n = nlayers;
                if bot_proj {
                    let cnt = (0..n - 1)
                        .filter(|&l| stride_coord[l] - stride_coord[l + 1] < min_height)
                        .count();
                    if cnt > 0 {
                        *limited_count += cnt;
                        *tangled_count += (0..n - 1)
                            .filter(|&l| stride_coord[l] - stride_coord[l + 1] < f64::MIN_POSITIVE)
                            .count();
                        for k in 1..n {
                            stride_coord[k] = stride_coord[k].min(stride_coord[k - 1] - min_height);
                        }
                    }
                } else {
                    let cnt = (0..n - 1)
                        .filter(|&l| stride_coord[l + 1] - stride_coord[l] < min_height)
                        .count();
                    if cnt > 0 {
                        *limited_count += cnt;
                        *tangled_count += (0..n - 1)
                            .filter(|&l| stride_coord[l + 1] - stride_coord[l] < f64::MIN_POSITIVE)
                            .count();
                        for k in 1..n {
                            stride_coord[k] = stride_coord[k].max(stride_coord[k - 1] + min_height);
                        }
                    }
                }
            }

            for (l, &n) in stride_ind.iter().enumerate() {
                coord[n] = stride_coord[l];
            }
            for &n in &stride_ind {
                if let Some(v) = velo {
                    v.set(ctx, n, (coord[n] - orig[n]) / dt);
                }
                if let Some(u) = update {
                    u.set(ctx, n, coord[n] - orig[n]);
                }
            }
        }
        ctx.info(CALLER, "Finished multilayer mapping", 8);
        Ok(())
    }
}

#[allow(dead_code)]
fn _unused(_: &dyn MeshView) {}
