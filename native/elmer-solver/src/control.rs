//! `Apply Source Control` (`ControlLinearSystem`, SolveCore.F90:6478): the
//! amplitudes of the `<var> Control [i]` source terms are chosen so that the
//! solution takes `Control Target Value [i]` at the control nodes (or, for
//! one control without a node, at the extremum easiest reached).
//!
//! Before the solve the sources, times their current amplitudes, are added
//! to the right-hand side. After it each source is solved for on its own
//! (`A dx = f_i`), the sensitivities give the amplitude changes `dc`
//! (`LuSolve`), and the field moves by `dc_i dx`. As in the Fortran that is
//! the `dx` of the last control for every `i`, and the nonlinear change of
//! the solve is the one before this correction.

use elmer_core::{BasisEval, ElemId, ParamMap, Params, Value};

use crate::error::SolverError;
use crate::solver::{solve_linear, SolverCtx};

const CALLER: &str = "ControlLinearSystem";

fn numbered(base: &str, i: usize, n: usize) -> String {
    if n > 1 {
        format!("{base} {i}")
    } else {
        base.to_string()
    }
}

/// `SetNodalSources`: the load vector of `name` from the body forces and
/// boundary conditions whose keywords start with it, integrated with the
/// framework's basis (`LocalSourceAssembly`), over the variable's dofs.
fn nodal_sources(ctx: &SolverCtx<'_>, name: &str) -> Result<Vec<f64>, SolverError> {
    let id = ctx.sim.slots[ctx.slot].variable.expect("checked by the caller");
    let v = ctx.sim.vars.by_id(id);
    let (dofs, nsize) = (v.dofs, v.values.len());
    let mut f = vec![0.0; nsize];
    let model = &ctx.sim.model;
    let prefixed = |p: &ParamMap| p.0.keys().any(|k| k.starts_with(name));
    let active_bc: Vec<bool> = model
        .boundary_conditions
        .iter()
        .map(|b| b.get("target boundaries").is_some() && prefixed(b))
        .collect();
    let active_bf: Vec<bool> = model.body_forces.iter().map(prefixed).collect();
    let (any_bc, any_bf) = (active_bc.contains(&true), active_bf.contains(&true));
    if !any_bc && !any_bf {
        return Ok(f);
    }
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature for source control (elmer-fe)".into()))?;
    let axi = model
        .simulation
        .string("coordinate system")
        .is_some_and(|c| elmer_core::normalize_key(c).starts_with("axi"));
    let nb = ctx.mesh.n_bulk() as u32;
    let first = if any_bf { 0 } else { nb };
    let last = if any_bc { nb + ctx.mesh.n_boundary() as u32 } else { nb };
    let mut ev = BasisEval::default();
    for e in (first..last).map(ElemId) {
        let section = if e.0 >= nb {
            match model.bc_of_tag(ctx.mesh.element(e).tag) {
                Some(bc) if active_bc[bc - 1] => &model.boundary_conditions[bc - 1],
                _ => continue,
            }
        } else {
            let Some(b) = ctx.body_params(e) else { continue };
            let Some(bf) = b.body.integer("body force").filter(|&k| k >= 1) else { continue };
            if !active_bf.get(bf as usize - 1).copied().unwrap_or(false) {
                continue;
            }
            &model.body_forces[bf as usize - 1]
        };
        let el = ctx.mesh.element(e);
        let n = el.nodes.len();
        let mut src = vec![vec![0.0; n]; dofs];
        if dofs == 1 {
            src[0] = ctx.real_at_nodes(section, name, e)?.unwrap_or_else(|| vec![0.0; n]);
        } else {
            let mut found = 0;
            for (i, s) in src.iter_mut().enumerate() {
                if let Some(x) = ctx.real_at_nodes(section, &format!("{name} {}", i + 1), e)? {
                    *s = x;
                    found += 1;
                }
            }
            if found == 0 {
                return Err(SolverError::fatal("SetNodalSources", format!("Could not find for any component: {name}")));
            }
        }
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
        let mut force = vec![vec![0.0; n]; dofs];
        for qp in fe.quadrature.gauss_points(el.etype, None) {
            fe.basis
                .eval(el.etype, &xyz, qp.uvw, &mut ev)
                .map_err(|m| SolverError::fatal("SetNodalSources", m))?;
            let mut w = qp.weight * ev.det_j;
            if axi {
                w *= xyz.iter().zip(&ev.basis).map(|(c, b)| c[0] * b).sum::<f64>();
            }
            for (i, s) in src.iter().enumerate() {
                let at_ip: f64 = s.iter().zip(&ev.basis).map(|(a, b)| a * b).sum();
                for p in 0..n {
                    force[i][p] += w * ev.basis[p] * at_ip;
                }
            }
        }
        for (p, &node) in el.nodes.iter().enumerate() {
            let Some(slot) = v.slot(node) else { continue };
            for (i, fi) in force.iter().enumerate() {
                f[dofs * slot + i] += fi[p];
            }
        }
    }
    Ok(f)
}

/// `GetControlNode`: `Control Node Index [i]` (1-based), else the node
/// with a dof nearest to `Control Node Coordinates` (row `i` of a
/// multi-row array, or `Control Node Coordinates i`), remembered as
/// `Control Node Index i`. `None` when neither is given.
pub(crate) fn control_node(ctx: &mut SolverCtx<'_>, i: usize) -> Option<usize> {
    let p = ctx.params();
    let idx = p
        .integer(&format!("control node index {i}"))
        .or_else(|| (i == 1).then(|| p.integer("control node index")).flatten());
    if let Some(k) = idx {
        return usize::try_from(k).ok().filter(|&k| k > 0);
    }
    let (arr, row) = match p.get("control node coordinates") {
        Some(Value::Real { rows, cols, data }) => ((*rows, *cols, data.clone()), i),
        _ => match p.get(&format!("control node coordinates {i}")) {
            Some(Value::Real { rows, cols, data }) => ((*rows, *cols, data.clone()), 1),
            _ => return None,
        },
    };
    let (rows, cols, data) = arr;
    let mut c = [0.0; 3];
    if cols == 1 {
        for (k, ck) in c.iter_mut().enumerate().take(rows.min(3)) {
            *ck = data[k];
        }
    } else {
        for (k, ck) in c.iter_mut().enumerate().take(cols.min(3)) {
            *ck = data[(row - 1) * cols + k];
        }
    }
    let id = ctx.sim.slots[ctx.slot].variable?;
    let v = ctx.sim.vars.by_id(id);
    let mut best = (f64::MAX, 0usize);
    for node in 0..ctx.mesh.n_nodes() as u32 {
        if v.slot(node).is_none() {
            continue;
        }
        let x = ctx.mesh.coord(node);
        let d = (x[0] - c[0]).powi(2) + (x[1] - c[1]).powi(2) + (x[2] - c[2]).powi(2);
        if d < best.0 {
            best = (d, node as usize + 1);
        }
    }
    if best.1 > 0 {
        ctx.info("GetControlNode", &format!("Control Node located to index: {}", best.1), 6);
        ctx.sim.slots[ctx.slot]
            .params
            .insert(&format!("Control Node Index {i}"), Value::Integer(vec![best.1 as i32]));
        Some(best.1)
    } else {
        None
    }
}

/// `LuSolve(n, A, x)`: Gaussian elimination with partial pivoting; `a` is
/// row-major `n x n`.
fn lu_solve(n: usize, a: &mut [f64], x: &mut [f64]) -> Result<(), SolverError> {
    for k in 0..n {
        let piv = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .expect("k < n");
        if a[piv * n + k] == 0.0 {
            return Err(SolverError::fatal("LuSolve", "singular control sensitivity matrix"));
        }
        if piv != k {
            for j in 0..n {
                a.swap(k * n + j, piv * n + j);
            }
            x.swap(k, piv);
        }
        for i in k + 1..n {
            let m = a[i * n + k] / a[k * n + k];
            for j in k..n {
                a[i * n + j] -= m * a[k * n + j];
            }
            x[i] -= m * x[k];
        }
    }
    for k in (0..n).rev() {
        let s: f64 = (k + 1..n).map(|j| a[k * n + j] * x[j]).sum();
        x[k] = (x[k] - s) / a[k * n + k];
    }
    Ok(())
}

/// `DefaultSolve` with `Apply Source Control`.
pub(crate) fn solve_controlled(ctx: &mut SolverCtx<'_>) -> Result<f64, SolverError> {
    if ctx.comm.size() > 1 {
        return Err(SolverError::fatal(CALLER, "Controlling of source terms implemented only in serial!"));
    }
    let p = ctx.params().clone();
    // `Control Use Loads`: the target and the sensitivities are the nodal
    // loads of the bulk system (CalculateLoads with NodalValues) rather
    // than the solution; the amplitude still moves the solution.
    let use_loads = p.logical("control use loads") == Some(true);
    if use_loads && p.logical("calculate boundary fluxes") == Some(true) {
        return Err(SolverError::Unsupported("Control Use Loads with Calculate Boundary Fluxes".into()));
    }
    let id = ctx.sim.slots[ctx.slot]
        .variable
        .ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    let (name, dofs) = {
        let v = ctx.sim.vars.by_id(id);
        (v.name.clone(), v.dofs)
    };
    let nc = p.integer("number of controls").unwrap_or(1).max(1) as usize;

    // PreSolve: the sources and their current amplitudes on the rhs.
    ctx.info(CALLER, "Applying controlled sources", 7);
    let mut f = Vec::with_capacity(nc);
    for i in 1..=nc {
        let mut amp = ctx.creal(&p, &numbered("control amplitude", i, nc))?;
        if amp.is_none() {
            amp = ctx.creal(&p, &numbered("initial control amplitude", i, nc))?;
        }
        let amp = amp.unwrap_or(0.0);
        let src = numbered(&format!("{name} control"), i, nc);
        ctx.info(CALLER, &format!("Computing source term for: {src}"), 7);
        let mut fi = nodal_sources(ctx, &src)?;
        for &(row, _) in &ctx.sim.slots[ctx.slot].asm.constrained {
            if let Some(x) = fi.get_mut(row as usize) {
                *x = 0.0;
            }
        }
        if amp.abs() > 1.0e-20 {
            let sys = ctx.sim.slots[ctx.slot].system.as_mut().ok_or_else(|| SolverError::fatal(CALLER, "no assembled system"))?;
            for (b, x) in sys.rhs.iter_mut().zip(&fi) {
                *b += amp * x;
            }
        }
        f.push(fi);
    }

    ctx.solve_system()?;

    // PostSolve: the amplitude changes.
    ctx.info(CALLER, "Dertermining source term amplitude", 7);
    let diag = p.logical("control diagonal") == Some(true);
    let dof0 = if dofs > 1 {
        p.integer("control target component")
            .ok_or_else(|| SolverError::fatal("ListGetInteger", "Requested property: [control target component] not found"))?
            as usize
    } else {
        1
    };
    let mut target = vec![0.0; nc];
    let mut cdof: Vec<Option<usize>> = vec![None; nc];
    for i in 1..=nc {
        let key = numbered("control target value", i, nc);
        target[i - 1] = ctx
            .creal(&p, &key)?
            .ok_or_else(|| SolverError::fatal("ListGetCReal", format!("Requested property: [{key}] not found")))?;
        if let Some(node) = control_node(ctx, i) {
            let v = ctx.sim.vars.by_id(id);
            cdof[i - 1] = v.slot(node as u32 - 1).map(|s| dofs * s + dof0 - 1);
        }
    }
    let extremum = cdof.iter().any(Option::is_none);
    if extremum && nc > 1 {
        return Err(SolverError::fatal(CALLER, format!("Extremum control cannot be used with {nc} controls!")));
    }
    if extremum && dofs > 1 {
        return Err(SolverError::Unsupported("extremum source control of a multi-dof variable".into()));
    }
    let mut sens = vec![0.0; nc * nc];
    let mut cval = vec![0.0; nc];
    let mut dx = Vec::new();
    let mut r0 = Vec::new();
    let schur = crate::solver::schur_matrix(ctx.sim, ctx.params());
    for i in 0..nc {
        let (log, comm) = (ctx.log, ctx.comm);
        let sim = &mut *ctx.sim;
        let slot = &mut sim.slots[ctx.slot];
        let sys = slot.system.as_mut().expect("solved above");
        let mut b = f[i].clone();
        b.truncate(sys.matrix.nrows);
        dx = vec![0.0; sys.matrix.nrows];
        let bnorm = b.iter().map(|x| x * x).sum::<f64>().sqrt();
        if bnorm > f64::MIN_POSITIVE {
            let env = crate::solver::BlockEnv {
                dofs: sim.vars.by_id(id).dofs,
                simulation: &sim.model.simulation,
                schur: schur.as_ref(),
                multigrid: None,
            };
            solve_linear(&sim.backends.linear, &env, &slot.params, &mut sys.matrix, &mut b, &mut dx, bnorm, false, log, comm)?;
        }
        let x0 = &sim.vars.by_id(id).values;
        let dr;
        let (y0, dy): (&[f64], &[f64]) = if use_loads {
            let n = sys.matrix.nrows;
            let (values, rhs) = slot
                .asm
                .bulk
                .as_ref()
                .ok_or_else(|| SolverError::fatal("CalculateLoads", "\"BulkValues\" are not associated!"))?;
            if i == 0 {
                r0 = crate::loads::residual(&sys.matrix, values, &x0[..n], &rhs[..n]);
            }
            // The control source stands in for BulkRhs.
            let mut fi = f[i].clone();
            fi.resize(n, 0.0);
            dr = crate::loads::residual(&sys.matrix, values, &dx, &fi);
            (&r0, &dr)
        } else {
            (x0, &dx)
        };
        if extremum {
            let val = target[0];
            let mut mincand = f64::MAX;
            for j in 0..dy.len() {
                if dy[j].abs() < f64::MIN_POSITIVE {
                    continue;
                }
                let cand = (val - y0[j]) / dy[j];
                if cand.abs() < mincand.abs() {
                    mincand = cand;
                    sens[0] = dy[j];
                    cval[i] = y0[j];
                }
            }
            ctx.info(CALLER, "Extremum value is easiest found in dof: 0", 7);
        } else {
            for (j, c) in cdof.iter().enumerate() {
                if diag && j != i {
                    continue;
                }
                sens[j * nc + i] = dy[c.expect("not extremum")];
            }
            cval[i] = y0[cdof[i].expect("not extremum")];
        }
    }
    let mut dc: Vec<f64> = target.iter().zip(&cval).map(|(t, v)| t - v).collect();
    lu_solve(nc, &mut sens, &mut dc)?;
    if let Some(r) = ctx.creal(&p, "control relaxation factor")? {
        dc.iter_mut().for_each(|d| *d *= r);
    }
    for i in 1..=nc {
        let key = numbered("Control Amplitude", i, nc);
        let amp = ctx.creal(ctx.params(), &key)?.unwrap_or(0.0) + dc[i - 1];
        ctx.sim.slots[ctx.slot].params.insert(&key, Value::real(amp));
        let v = ctx.sim.vars.by_id_mut(id);
        for (x, d) in v.values.iter_mut().zip(&dx) {
            *x += dc[i - 1] * d;
        }
        v.values_changed = true;
        ctx.info(CALLER, &format!("Applied {key}: {}", crate::fortran_fmt::es(amp, 15, 6)), 5);
    }
    ctx.sync_model();
    Ok(ctx.convergence().map_or(0.0, |c| c.norm))
}


/// `ApplyExplicitControl` (SolveCore.F90:6441): for each of `Number Of
/// Controls`, the solver's field at its control node minus `Control Target
/// Value [i]`, kept in the global vector `cpar` that keywords may depend on.
pub(crate) fn apply_explicit_control(ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
    const CALLER: &str = "ApplyExplicitControl";
    let n = ctx.params().integer("number of controls").unwrap_or(0).max(0) as usize;
    if n == 0 {
        ctx.warn(CALLER, "Explicit control points requested but no controls available!");
        return Ok(());
    }
    let p = ctx.params();
    if p.get("control variable").is_some() || (1..=n).any(|i| p.get(&format!("control variable {i}")).is_some()) {
        return Err(SolverError::Unsupported("Control Variable other than the solver's (explicit control)".into()));
    }
    let id = ctx.sim.slots[ctx.slot]
        .variable
        .ok_or_else(|| SolverError::fatal("GetControlValue", "Could not find keyword for: Control Variable"))?;
    let mut vals = Vec::with_capacity(n);
    for i in 1..=n {
        let node = control_node(ctx, i).ok_or_else(|| SolverError::fatal("GetControlValue", "Could not find control node!"))?;
        let v = ctx.sim.vars.by_id(id);
        let dof = if v.dofs > 1 {
            ctx.params().integer("control target component").ok_or_else(|| {
                SolverError::fatal("ListGetInteger", "Requested property: [control target component] not found")
            })? as usize
        } else {
            1
        };
        let s = v
            .slot(node as u32 - 1)
            .ok_or_else(|| SolverError::fatal("GetControlValue", "Could not find control node!"))?;
        let val = v.values[v.dofs * s + dof - 1];
        let p = ctx.params().clone();
        let target = match ctx.creal(&p, "control target value")? {
            Some(t) => t,
            None => ctx.creal(&p, &format!("control target value {i}"))?.unwrap_or(0.0),
        };
        vals.push(val - target);
    }
    let cpar = match ctx.sim.vars.id("cpar") {
        Some(c) => c,
        None => ctx
            .sim
            .vars
            .add(crate::Variable::new("cpar", elmer_core::FieldKind::Global, n, None, vec![0.0; n]))?,
    };
    let c = ctx.sim.vars.by_id_mut(cpar);
    c.values[..n].copy_from_slice(&vals);
    c.values_changed = true;
    ctx.sim.publish_live();
    Ok(())
}
