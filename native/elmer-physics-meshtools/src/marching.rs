//! `MarchingODESolver` (MarchingODESolver.F90): the ODE `dx/dt = f`
//! integrated along the columns of an extruded mesh, the march coordinate
//! being the `Active Coordinate` travelled at the `Draw Velocity`, by the
//! Newmark (`Newmark Beta`, or the `Timestepping Method`) rule. A transient
//! run advances the solution by as many layers as the timestep covers.
//!
//! Ported: one-component variables, `<var>: Source` and `<var>: Time
//! Derivative Coefficient` material keywords, `Apply BCs Only`, `Parabolic
//! Model`, `Minimum Cons`, the nonlinear iterations of a layer and the
//! transient marching. Refused: several components, `Use Internal Values`
//! (with `Additional Internal Variable` and `Timestep Divisions`),
//! `Enforce Unity Sum` and discontinuous Galerkin solvers.

use elmer_core::{normalize_key, ElemId, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::structure::{detect_extruded_structure, extruded_coordinate, Mask, Want};

const CALLER: &str = "MarchingODESolver";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(MarchingOde::default()))
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// `MarchingODESolver_init`: `Time Derivative Order = 1`, `Use Global Mass
/// Matrix` and `No Matrix` (each unless given). Call it before the solver's
/// variable and matrix are created.
pub fn marching_defaults(p: &mut ParamMap) {
    for (k, v) in [
        ("Time derivative order", Value::Integer(vec![1])),
        ("Use Global Mass Matrix", Value::Logical(true)),
        ("No Matrix", Value::Logical(true)),
    ] {
        if p.get(k).is_none() {
            p.insert(k, v);
        }
    }
}

/// The structure found on the first visit (the Fortran's `SAVE`s).
struct Columns {
    /// `UpPointer`, by node.
    up: Vec<usize>,
    /// The bottom nodes of the marched columns, in node order (`InvPerm`).
    bottom: Vec<usize>,
    layers: usize,
    /// `Coord`: the extruded coordinate by node.
    coord: Vec<f64>,
    hparam: f64,
}

#[derive(Default)]
pub struct MarchingOde {
    columns: Option<Columns>,
}

impl Solver for MarchingOde {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        for k in ["use internal values", "enforce unity sum"] {
            if flag(p, k) {
                return Err(unsupported(k));
            }
        }
        if p.get("additional internal variable").is_some() {
            return Err(unsupported("Additional Internal Variable"));
        }
        if p.integer("timestep divisions").unwrap_or(1) > 1 {
            return Err(unsupported("Timestep Divisions"));
        }
        if ctx.sim.model.solvers.iter().any(|s| flag(s, "discontinuous galerkin")) {
            return Err(unsupported("discontinuous Galerkin solvers"));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "-----------------------------------------------------", 6);
        ctx.info(CALLER, "Solving ODE on moving coordinates in structured mesh", 4);
        let p = ctx.params().clone();
        let (var_name, dofs, perm) = {
            let v = ctx
                .primary()
                .ok_or_else(|| SolverError::fatal(CALLER, "Solver variable not associated!"))?;
            let perm = v
                .perm
                .ok_or_else(|| SolverError::fatal(CALLER, "the variable has no Perm"))?
                .to_vec();
            (v.name.to_string(), v.dofs, perm)
        };
        if dofs != 1 {
            return Err(unsupported("several components"));
        }
        if self.columns.is_none() {
            self.columns = Some(self.initialize(ctx, &p, &perm)?);
        }
        let cols = self.columns.as_ref().expect("initialized");
        let n = cols.bottom.len();

        let beta = match ctx.creal(&p, "newmark beta")? {
            Some(b) => b,
            None => match p.string("timestepping method").map(normalize_key).as_deref() {
                Some("explicit euler") => 0.0,
                Some("crank-nicolson") => 0.5,
                Some("newmark") => {
                    return Err(SolverError::fatal(CALLER, "> Newmark Beta < not found"));
                }
                _ => 1.0,
            },
        };
        let parabolic = flag(&p, "parabolic model");
        let velo = match ctx.creal(ctx.simulation_params(), "draw velocity")? {
            Some(v) => v,
            None => ctx
                .creal(&p, "draw velocity")?
                .ok_or_else(|| SolverError::fatal(CALLER, ">Draw Velocity< is needed for marching solver!"))?,
        };
        let max_iter = p.integer("nonlinear system max iterations").unwrap_or(1).max(0) as usize;
        let tol = ctx.creal(&p, "nonlinear system convergence tolerance")?.unwrap_or(0.0);
        let min_cons = ctx.creal(&p, "minimum cons")?;

        // FirstExtrudedMaterial: the material of the first active element.
        let first = ctx
            .active_elements()
            .next()
            .ok_or_else(|| SolverError::fatal(CALLER, "Could not set material for extruded domain!"))?;
        let material = ctx
            .body_params(first)
            .and_then(|b| b.material)
            .cloned()
            .ok_or_else(|| SolverError::fatal(CALLER, "Could not set material for extruded domain!"))?;
        let f_key = format!("{var_name}: source");
        let c_key = format!("{var_name}: time derivative coefficient");
        if material.get(&format!("{var_name}: reaction coefficient")).is_some() {
            return Err(SolverError::fatal(CALLER, "Code some more to account for reaction term!"));
        }
        let have_f = material.get(&f_key).is_some();
        let have_c = material.get(&c_key).is_some();

        let mut dtn = 0usize;
        if transient {
            let dth = cols.hparam / velo;
            let steps = (dt / dth).round();
            if steps >= cols.layers as f64 {
                ctx.info(CALLER, "Timestep so large than we can use steady algo!", 4);
            } else if (dt / dth - steps).abs() > 0.01 {
                return Err(SolverError::fatal(CALLER, "Timesteps are not matching"));
            } else {
                dtn = steps as usize;
                ctx.info(CALLER, &format!("Number of marching steps for each timestep: {dtn}"), 5);
            }
        }
        let do_transient = dtn > 0;

        // The layer's nodes and the field there.
        let get = |ctx: &SolverCtx<'_>, nodes: &[usize]| -> Vec<f64> {
            let v = ctx.primary().expect("checked");
            nodes
                .iter()
                .map(|&k| {
                    let x = v.values[perm[k] as usize];
                    if parabolic {
                        0.5 * x.powi(2)
                    } else {
                        x
                    }
                })
                .collect()
        };
        let set = |ctx: &mut SolverCtx<'_>, nodes: &[usize], x: &[f64]| {
            let v = ctx.primary_mut().expect("checked");
            for (&k, &xi) in nodes.iter().zip(x) {
                v.values[perm[k] as usize] = if parabolic { (2.0 * xi).sqrt() } else { xi };
            }
        };
        // GetCoefficients: the source (and coefficient) at the layer.
        let coefficients = |ctx: &SolverCtx<'_>, nodes: &[usize]| -> Result<(Vec<f64>, Vec<f64>), SolverError> {
            let ids: Vec<u32> = nodes.iter().map(|&k| k as u32).collect();
            let eval = |key: &str| -> Result<Vec<f64>, SolverError> {
                match material.get(key) {
                    Some(v) => elmer_solver::expr::eval_at_nodes_with(
                        v,
                        &ids,
                        &ctx.fields(),
                        ctx.sim.backends.matc.as_deref(),
                        &ctx.sim.backends.udfs,
                    ),
                    None => Ok(vec![0.0; ids.len()]),
                }
            };
            let f = if have_f { eval(&f_key)? } else { vec![0.0; ids.len()] };
            let c = if have_c { eval(&c_key)? } else { vec![0.0; ids.len()] };
            Ok((f, c))
        };

        for _dti in 0..dtn.max(1) {
            let mut inv: Vec<usize> = cols.bottom.clone();
            let mut xvec = get(ctx, &inv);
            let (mut f0, mut c0) = coefficients(ctx, &inv)?;
            let mut x0 = xvec.clone();
            for layer in 1..=cols.layers {
                let prev = inv.clone();
                inv = prev.iter().map(|&k| cols.up[k]).collect();
                if n > 0 && inv[0] == prev[0] {
                    return Err(SolverError::fatal(
                        CALLER,
                        format!("InvPerm is the same on different layers: {}", inv[0] + 1),
                    ));
                }
                let dz = if n > 0 { cols.coord[inv[0]] - cols.coord[prev[0]] } else { 0.0 };
                let dtime = dz / velo;
                if dtime < f64::EPSILON {
                    return Err(SolverError::fatal(CALLER, "Cannot march if timestep is ~zero"));
                }
                let mut xivec = get(ctx, &inv);
                for _iter in 0..max_iter {
                    if max_iter > 1 {
                        xivec.copy_from_slice(&xvec);
                    }
                    let (mut f, mut c) = coefficients(ctx, &inv)?;
                    if (beta - 1.0).abs() >= f64::EPSILON {
                        if have_f {
                            for i in 0..n {
                                f[i] = beta * f[i] + (1.0 - beta) * f0[i];
                            }
                        }
                        if have_c {
                            for i in 0..n {
                                c[i] = beta * c[i] + (1.0 - beta) * c0[i];
                            }
                        }
                    }
                    for i in 0..n {
                        xvec[i] = if have_c {
                            x0[i] + (dtime / 1.0) * f[i] / c[i]
                        } else {
                            x0[i] + (dtime / 1.0) * f[i]
                        };
                    }
                    if let Some(mc) = min_cons {
                        for x in &mut xvec {
                            *x = x.max(mc);
                        }
                    }
                    let mut nrm = 0.0;
                    let mut chg = 0.0;
                    for i in 0..n {
                        nrm += xvec[i] * xvec[i];
                        let d = xvec[i] - xivec[i];
                        chg += d * d;
                    }
                    let change = chg.sqrt() / nrm.sqrt();
                    set(ctx, &inv, &xvec);
                    if change < tol {
                        break;
                    }
                }
                if layer < cols.layers {
                    if do_transient {
                        set(ctx, &inv, &xivec);
                        (f0, c0) = coefficients(ctx, &inv)?;
                        x0 = xivec.clone();
                        set(ctx, &inv, &xvec);
                    } else {
                        (f0, c0) = coefficients(ctx, &inv)?;
                        x0 = xvec.clone();
                    }
                }
            }
        }
        ctx.info(CALLER, "All done", 5);
        Ok(())
    }
}

impl MarchingOde {
    /// The first visit: `DetectExtrudedStructure` with the variable as the
    /// mask (its `Extruded Coordinate`), the constant mesh parameter, and
    /// the bottom nodes of the columns that stay in the variable (and on
    /// the boundary with `Apply BCs Only`).
    fn initialize(&self, ctx: &mut SolverCtx<'_>, p: &ParamMap, perm: &[u32]) -> Result<Columns, SolverError> {
        ctx.info(CALLER, "Initializing structured mesh and ODE structures", 6);
        let dir = p.integer("active coordinate").unwrap_or(0);
        if !(1..=3).contains(&dir) {
            return Err(SolverError::fatal("StructuredMeshMapper", "Invalid value for Active Coordinate"));
        }
        let dir = dir as usize - 1;
        if flag(p, "mapping original coordinates") || p.get("mapping coordinate transformation").is_some() {
            return Err(unsupported("Mapping Original Coordinates / Coordinate Transformation"));
        }
        let mut unit = [0.0; 3];
        unit[dir] = if flag(p, "project to bottom") { -1.0 } else { 1.0 };
        let eps = ctx.creal(p, "dot product tolerance")?.unwrap_or(1.0e-4);
        let size = perm.iter().filter(|&&s| s != NO_DOF).map(|&s| s as usize + 1).max().unwrap_or(0);
        let mask = Mask {
            perm: perm.to_vec(),
            size,
        };
        extruded_coordinate(ctx, &mask, dir)?;
        let want = Want { layers: true, mid: false };
        let s = detect_extruded_structure(ctx.mesh, &ctx.sim.model, unit, eps, want, Some(&mask))?;
        ctx.info(CALLER, &format!("Number of element layers: {}", s.layers), 7);
        let ext = ctx.variable("extruded coordinate").expect("added above");
        let ep = ext.perm.expect("on the mask");
        let nn = ctx.mesh.n_nodes();
        let mut coord = vec![0.0; nn];
        let (mut cmin, mut cmax) = (f64::INFINITY, f64::NEG_INFINITY);
        for k in 0..nn {
            if ep[k] != NO_DOF {
                coord[k] = ext.values[ep[k] as usize];
            }
        }
        for &v in ext.values.iter() {
            cmin = cmin.min(v);
            cmax = cmax.max(v);
        }
        if cmax - cmin < f64::EPSILON {
            return Err(SolverError::fatal(CALLER, "We cannot march when dz is ~zero!"));
        }
        let hparam = (cmax - cmin) / s.layers as f64;
        let mut bc_node: Option<Vec<bool>> = None;
        if flag(p, "apply bcs only") {
            let mut b = vec![false; nn];
            let nb = ctx.mesh.n_bulk();
            for t in nb..nb + ctx.mesh.n_boundary() {
                for &k in ctx.mesh.element(ElemId(t as u32)).nodes {
                    b[k as usize] = true;
                }
            }
            if b.iter().any(|&x| x) {
                bc_node = Some(b);
            }
        }
        let on_bc = |k: usize| bc_node.as_ref().is_none_or(|b| b[k]);
        let mut bottom = Vec::new();
        for i in 0..nn {
            if perm[i] == NO_DOF || !on_bc(i) {
                continue;
            }
            if s.bot[i] != i {
                continue;
            }
            let up = s.up[i];
            if perm[up] == NO_DOF || !on_bc(up) {
                continue;
            }
            bottom.push(i);
        }
        ctx.info(CALLER, &format!("Number of bottom nodes: {}", bottom.len()), 7);
        Ok(Columns {
            up: s.up,
            bottom,
            layers: s.layers,
            coord,
            hparam,
        })
    }
}
