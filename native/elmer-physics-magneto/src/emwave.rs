//! `EMWaveSolver.F90` (fem/src/modules): electromagnetic waves in time,
//! `eps d2E/dt2 + sigma dE/dt + curl (1/mu) curl E = -dJ/dt`, on edge
//! elements (`Procedure = "EMWaveSolver" "EMWaveSolver"`), and its nodal
//! post-processor `EMWaveCalcFields`.
//!
//! [`init_section`] is `EMWaveSolver_Init0`: the `Element` `EdgeElementStyle`
//! picks, `Hcurl Basis`, `Use Global Mass Matrix` with `Constant Bulk
//! Matrix`, the variable `E` (not saved), `Time Derivative Order = 2`.
//! [`calc_init_section`] is `EMWaveCalcFields_Init0` + `_Init`: the
//! primary solver index, the `-nooutput hr_dummy` variable and the exported
//! nodal `Elfield[Elfield:3]`.
//!
//! The solver: `LocalMatrix` (`STIFF = 1/mu (curl, curl)`, `MASS = eps (.,.)`,
//! `DAMP = sigma (.,.)`, `FORCE = -(dJ/dt, .)`) on `GaussPoints(EdgeBasis,
//! PReferenceElement=Piola)` with the Piola basis, through
//! `Default2ndOrderTime` (with `Use Global Mass Matrix`: the global mass and
//! damping matrices and `Time2ndOrder_CRS` at `DefaultFinishAssembly`); on
//! the boundary edges `Electric Damping Coefficient` (`DAMP = B/mu (w, w)`)
//! and `Magnetic Boundary Load` with the trace of the edge basis
//! (`EdgeElementInfo` on the line), and the tangential `E {e}` Dirichlet
//! values (`VectorElementEdgeDOFs`). `DefaultSolve`.
//!
//! The constant bulk matrix is reassembled at each step: it is the same
//! matrix `RestoreBulkMatrix` gives back.
//!
//! Ported: 2D, first-order Piola elements (triangles, quadrilaterals).
//! Refused: 3D, the second-order and second-kind bases, `TEM Potential`,
//! eigen analysis, the elemental (DG) fields and the time derivatives of the
//! post-processor.

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementType, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::case::has_prefix;
use crate::topology::edge_map;

const CALLER: &str = "EMWaveSolver";

pub const PROC: (&str, &str) = ("EMWaveSolver", "EMWaveSolver");
pub const CALC: (&str, &str) = ("EMWaveSolver", "EMWaveCalcFields");

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROC.0, PROC.1, |_| Ok(Box::new(EmWave)));
    registry.register(CALC.0, CALC.1, |_| Ok(Box::new(EmWaveCalc)));
}

fn is_proc(p: &ParamMap, proc: (&str, &str)) -> bool {
    elmer_solver::solver::procedure_of(p).is_some_and(|(l, e)| {
        l.trim_start_matches("./").eq_ignore_ascii_case(proc.0) && e.eq_ignore_ascii_case(proc.1)
    })
}

pub fn is_section(p: &ParamMap) -> bool {
    is_proc(p, PROC)
}

pub fn is_calc_section(p: &ParamMap) -> bool {
    is_proc(p, CALC)
}

fn unsupported(m: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {}", m.into()))
}

/// `EdgeElementStyle`: (Piola, second kind, degree).
fn style(p: &ParamMap) -> (bool, bool, u8) {
    let quadratic = p.logical("quadratic approximation") == Some(true);
    let second = p.logical("second kind basis") == Some(true);
    let piola = p.logical("use piola transform") == Some(true) || quadratic || second;
    (piola, second, if quadratic { 2 } else { 1 })
}

fn add_new(p: &mut ParamMap, key: &str, v: Value) {
    if p.get(key).is_none() {
        p.insert(key, v);
    }
}

/// `EMWaveSolver_Init0`.
pub fn init_section(p: &mut ParamMap) {
    if p.get("element").is_none() {
        let (piola, second, k) = style(p);
        let e = if k == 2 {
            "n:0 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2"
        } else if second {
            "n:0 e:2"
        } else if piola {
            "n:0 e:1 -quad b:2 -brick b:3 -quad_face b:2"
        } else {
            "n:0 e:1"
        };
        p.insert("Element", Value::Str(e.into()));
    }
    add_new(p, "Hcurl Basis", Value::Logical(true));
    if p.logical("constant bulk matrix") == Some(true) || p.logical("eigen analysis") == Some(true) {
        add_new(p, "Use Global Mass Matrix", Value::Logical(true));
    }
    add_new(p, "Variable Output", Value::Logical(false));
    add_new(p, "Variable", Value::Str("E".into()));
    add_new(p, "Linear System Complex", Value::Logical(false));
    add_new(p, "Time derivative order", Value::Integer(vec![2]));
}

/// `EMWaveCalcFields_Init0` and `EMWaveCalcFields_Init` for solver
/// section `i` of `solvers`.
pub fn calc_init_section(solvers: &mut [ParamMap], i: usize) {
    let soln = solvers.iter().position(|s| {
        s.string("procedure").is_some_and(|p| p.contains("EMWaveSolver") && !p.contains("EMWaveCalcFields"))
    });
    let p = &mut solvers[i];
    if let Some(s) = soln {
        p.insert("Primary Solver Index", Value::Integer(vec![s as i32 + 1]));
    }
    p.insert("Variable", Value::Str("-nooutput hr_dummy".into()));
    p.insert("Linear System Refactorize", Value::Logical(false));
    add_new(p, "Skip Compute Nonlinear Change", Value::Logical(true));
    if p.logical("calculate nodal fields") == Some(false) {
        return;
    }
    let mut k = 1;
    while p.get(&format!("exported variable {k}")).is_some() {
        k += 1;
    }
    p.insert(&format!("Exported Variable {k}"), Value::Str("Elfield[Elfield:3]".into()));
}

/// `ListGetElementReal(h, Basis, Element)`: a constant as it is, else
/// `SUM(Basis * nodal values)`; `None` when the keyword is not given.
fn elem_real(ctx: &SolverCtx<'_>, p: Option<&ParamMap>, key: &str, e: ElemId, basis: &[f64]) -> Result<Option<f64>, SolverError> {
    let Some(p) = p else { return Ok(None) };
    match p.get(key) {
        None => Ok(None),
        Some(Value::Real { data, .. }) if data.len() == 1 => Ok(Some(data[0])),
        Some(_) => {
            let v = ctx.real_at_nodes(p, key, e)?.unwrap_or_default();
            let mut s = 0.0;
            for (b, x) in basis.iter().zip(&v) {
                s += b * x;
            }
            Ok(Some(s))
        }
    }
}

/// `GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=.TRUE.)`:
/// `EdgeElementGaussPoints` of a first-order basis.
fn bulk_rule(code: u16) -> Result<elmer_fe::quadrature::GaussRule, SolverError> {
    let fatal = |m: String| SolverError::fatal(CALLER, m);
    match code {
        404 => elmer_fe::quadrature::gauss_points(ElementType(404), Some(9)).map_err(|e| fatal(e.to_string())),
        c => Err(unsupported(format!("element type {c}"))),
    }
}

fn consts(ctx: &SolverCtx<'_>) -> Result<(f64, f64), SolverError> {
    let c = ctx.constants();
    let eps0 = c
        .real("permittivity of vacuum")
        .ok_or_else(|| SolverError::fatal("EMWaveSolver_Init0", "> Permittivity of Vacuum < is required"))?;
    let mu0 = c
        .real("permeability of vacuum")
        .ok_or_else(|| SolverError::fatal("EMWaveSolver_Init0", "> Permeability of Vacuum < is required"))?;
    Ok((eps0, mu0))
}

struct EmWave;

fn check(ctx: &SolverCtx<'_>, sp: &ParamMap) -> Result<u8, SolverError> {
    let (piola, second, k) = style(sp);
    if ctx.mesh.dim() == 2 && !piola {
        return Err(SolverError::fatal(CALLER, "A 2D model needs Use Piola Transform = True"));
    }
    if ctx.mesh.dim() != 2 {
        return Err(unsupported("3D"));
    }
    if second || k != 1 {
        return Err(unsupported("the second-order and second-kind bases"));
    }
    if sp.logical("eigen analysis") == Some(true) {
        return Err(unsupported("Eigen Analysis"));
    }
    for bc in &ctx.sim.model.boundary_conditions {
        if has_prefix(bc, "tem potential") {
            return Err(unsupported("'TEM Potential'"));
        }
    }
    Ok(k)
}

/// The parent's local edge of a boundary line and the slot of its DOF.
fn boundary_edge(ctx: &SolverCtx<'_>, b: ElemId, slots: &mut Vec<u32>) -> Result<Option<(ElemId, usize, u32)>, SolverError> {
    let el = ctx.mesh.element(b);
    let nb = ctx.mesh.n_bulk();
    let Some(parent) = el.parents.iter().flatten().copied().find(|p| (p.0 as usize) < nb) else {
        return Ok(None);
    };
    let pel = ctx.mesh.element(parent);
    let emap = edge_map(pel.etype.0 / 100);
    let Some(l) = emap
        .iter()
        .position(|&[i, j]| el.nodes[..2].contains(&pel.nodes[i]) && el.nodes[..2].contains(&pel.nodes[j]))
    else {
        return Ok(None);
    };
    ctx.element_slots(parent, slots);
    Ok(Some((parent, l, slots[l])))
}

/// `VectorElementEdgeDOFs` of a first-order basis on 2D boundary edges:
/// `Integral(1) = int (L + VL.e)` along the parent's edge `(j, k)`,
/// negated when `k < j` globally.
fn edge_dirichlet(ctx: &SolverCtx<'_>, var: &str) -> Result<Vec<(u32, f64)>, SolverError> {
    let prefix = format!("{var} {{e}}");
    let rule = elmer_fe::quadrature::gauss_points(ElementType(202), None).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
    let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
    let mut slots = Vec::new();
    let mut fixed = Vec::new();
    let mut ev = BasisEval::default();
    for &b in &boundary {
        let Some(bc) = ctx.bc_params(b) else { continue };
        if !has_prefix(bc, &prefix) {
            continue;
        }
        let el = ctx.mesh.element(b);
        if el.etype.0 != 202 {
            return Err(unsupported(format!("'{prefix}' on element type {}", el.etype.0)));
        }
        let Some((parent, l, slot)) = boundary_edge(ctx, b, &mut slots)? else { continue };
        let pel = ctx.mesh.element(parent);
        let [j, k] = edge_map(pel.etype.0 / 100)[l];
        let (nj, nk) = (pel.nodes[j], pel.nodes[k]);
        let order: [usize; 2] = if el.nodes[0] == nj { [0, 1] } else { [1, 0] };
        let at = |key: &str| -> Result<[f64; 2], SolverError> {
            Ok(match ctx.real_at_nodes(bc, key, b)? {
                Some(v) => [v[order[0]], v[order[1]]],
                None => [0.0; 2],
            })
        };
        let load = at(&prefix)?;
        let vload = [at(&format!("{prefix} 1"))?, at(&format!("{prefix} 2"))?, at(&format!("{prefix} 3"))?];
        let (pj, pk) = (ctx.mesh.coord(nj), ctx.mesh.coord(nk));
        let d = [pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let dir = d.map(|x| x / len);
        let xyz = [pj, pk];
        let mut integral = 0.0f64;
        for t in 0..rule.len() {
            elmer_fe::LagrangeBasis::new(3)
                .eval(ElementType(202), &xyz, [rule.u[t], 0.0, 0.0], &mut ev)
                .map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let s = rule.s[t] * ev.det_j;
            let (b0, b1) = (ev.basis[0], ev.basis[1]);
            let lv = load[0] * b0 + load[1] * b1;
            let vl: [f64; 3] = std::array::from_fn(|c| vload[c][0] * b0 + vload[c][1] * b1);
            let ve = (vl[0] * dir[0] + vl[1] * dir[1]) + vl[2] * dir[2];
            integral += s * (lv + ve);
        }
        if nk < nj {
            integral = -integral;
        }
        if slot != NO_DOF {
            fixed.push((slot, integral));
        }
    }
    Ok(fixed)
}

impl Solver for EmWave {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        if !transient {
            return Err(unsupported("a steady-state run"));
        }
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::fatal(CALLER, "the edge element DOFs (no DofsBackend)"));
        }
        let sp = ctx.params().clone();
        let degree = check(ctx, &sp)?;
        if ctx.primary().map(|v| v.dofs) != Some(1) {
            return Err(SolverError::fatal(CALLER, "Invalid variable size"));
        }
        let (eps0, mu0) = consts(ctx)?;
        let var = ctx.variable_name().unwrap_or("e").to_ascii_lowercase();
        let max_iter = sp.integer("nonlinear system max iterations").unwrap_or(1).max(1);
        for _ in 0..max_iter {
            // DoBulkAssembly.
            ctx.default_initialize()?;
            let nbulk = ctx.mesh.n_bulk();
            let active: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nbulk).collect();
            let mut slots = Vec::new();
            for &e in &active {
                ctx.element_slots(e, &mut slots);
                if slots.contains(&NO_DOF) {
                    return Err(SolverError::fatal(CALLER, format!("element {} lacks DOFs", e.0)));
                }
                let nd = slots.len();
                let el = ctx.mesh.element(e);
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let gindex: Vec<i64> = el.nodes.iter().map(|&k| k as i64 + 1).collect();
                let etype = el.etype;
                let bp = ctx.body_params(e);
                let (material, body_force) = (bp.as_ref().and_then(|b| b.material), bp.as_ref().and_then(|b| b.body_force));
                let rule = bulk_rule(etype.0)?;
                let mut ls = LocalSystem {
                    n: nd,
                    stiff: vec![0.0; nd * nd],
                    mass: Some(vec![0.0; nd * nd]),
                    damp: Some(vec![0.0; nd * nd]),
                    force: vec![0.0; nd],
                };
                for t in 0..rule.len() {
                    let (ev, basis, _) = crate::piola::edge_info_cdim(etype, &xyz, &gindex, [rule.u[t], rule.v[t], rule.w[t]], degree, ctx.mesh.dim())?;
                    if ev.dofs != nd {
                        return Err(SolverError::fatal(CALLER, format!("{} edge basis functions for {nd} DOFs", ev.dofs)));
                    }
                    let weight = ev.det * rule.s[t];
                    let eps = eps0 * elem_real(ctx, material, "relative permittivity", e, &basis)?.unwrap_or(0.0);
                    let mu = mu0 * elem_real(ctx, material, "relative permeability", e, &basis)?.unwrap_or(0.0);
                    let muinv = 1.0 / mu;
                    let cond = elem_real(ctx, material, "electric conductivity", e, &basis)?.unwrap_or(0.0);
                    let l: [f64; 3] = [
                        elem_real(ctx, body_force, "current density rate 1", e, &basis)?.unwrap_or(0.0),
                        elem_real(ctx, body_force, "current density rate 2", e, &basis)?.unwrap_or(0.0),
                        elem_real(ctx, body_force, "current density rate 3", e, &basis)?.unwrap_or(0.0),
                    ];
                    let (w, r) = (&ev.edge, &ev.curl);
                    let dot3 = |a: &[f64; 3], b: &[f64; 3]| (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2];
                    let (mass, damp) = (ls.mass.as_mut().unwrap(), ls.damp.as_mut().unwrap());
                    for i in 0..nd {
                        ls.force[i] -= dot3(&l, &w[i]) * weight;
                        for j in 0..nd {
                            ls.stiff[i * nd + j] += muinv * dot3(&r[i], &r[j]) * weight;
                            mass[i * nd + j] += eps * dot3(&w[j], &w[i]) * weight;
                            damp[i * nd + j] += cond * dot3(&w[j], &w[i]) * weight;
                        }
                    }
                }
                ctx.default_second_order_time(&mut ls, e, dt)?;
                ctx.default_update(&ls, e)?;
            }
            ctx.default_finish_bulk_assembly()?;

            // DoBoundaryAssembly: the boundary edges.
            let boundary: Vec<ElemId> = ctx.boundary_elements().collect();
            let mut bslots = Vec::new();
            for &b in &boundary {
                let Some(bc) = ctx.bc_params(b) else { continue };
                let el = ctx.mesh.element(b);
                if el.etype.0 / 100 == 1 {
                    continue;
                }
                if el.etype.0 != 202 {
                    return Err(unsupported(format!("boundary elements of type {}", el.etype.0)));
                }
                let Some((parent, l, slot)) = boundary_edge(ctx, b, &mut slots)? else { continue };
                ctx.element_slots(b, &mut bslots);
                if bslots.as_slice() != [slot] {
                    return Err(SolverError::fatal(CALLER, "a boundary edge without its edge DOF"));
                }
                let pel = ctx.mesh.element(parent);
                let [j, k] = edge_map(pel.etype.0 / 100)[l];
                let (nj, nk) = (pel.nodes[j], pel.nodes[k]);
                // `Mesh % Edges(k)`: its only element's node order.
                let xs = [ctx.mesh.coord(nj), ctx.mesh.coord(nk)];
                let reverse = nk < nj;
                let pmat = ctx.body_params(parent).and_then(|bp| bp.material);
                let rule = elmer_fe::quadrature::gauss_points(ElementType(202), Some(2)).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
                let mut damp = 0.0f64;
                let mut force = 0.0f64;
                let mut any_damp = false;
                for t in 0..rule.len() {
                    let u = rule.u[t];
                    // LineNodalPBasis and its derivative.
                    let basis = [0.5 * (1.0 - u), 0.5 * (1.0 + u)];
                    let dl = [-0.5, 0.5];
                    // ElementMetric of the line in 2D.
                    let dx: [f64; 3] = std::array::from_fn(|c| xs[0][c] * dl[0] + xs[1][c] * dl[1]);
                    let g = dx[0] * dx[0] + dx[1] * dx[1];
                    let metric = 1.0 / g;
                    let det_j = g.sqrt();
                    let eb = if reverse { -0.5 } else { 0.5 };
                    let w: [f64; 3] = std::array::from_fn(|c| if c < 2 { (dx[c] * metric) * eb } else { 0.0 });
                    let weight = det_j * rule.s[t];
                    let mu = mu0 * elem_real(ctx, pmat, "relative permeability", parent, &basis)?.unwrap_or(0.0);
                    let muinv = 1.0 / mu;
                    let lb: [f64; 3] = [
                        elem_real(ctx, Some(bc), "magnetic boundary load 1", b, &basis)?.unwrap_or(0.0),
                        elem_real(ctx, Some(bc), "magnetic boundary load 2", b, &basis)?.unwrap_or(0.0),
                        elem_real(ctx, Some(bc), "magnetic boundary load 3", b, &basis)?.unwrap_or(0.0),
                    ];
                    let lw = (lb[0] * w[0] + lb[1] * w[1]) + lb[2] * w[2];
                    force -= muinv * lw * weight;
                    if let Some(bv) = elem_real(ctx, Some(bc), "electric damping coefficient", b, &basis)? {
                        any_damp = true;
                        let ww = (w[0] * w[0] + w[1] * w[1]) + w[2] * w[2];
                        damp += muinv * bv * ww * weight;
                    }
                }
                let mut ls = LocalSystem { n: 1, stiff: vec![0.0], mass: Some(vec![0.0]), damp: Some(vec![if any_damp { damp } else { 0.0 }]), force: vec![force] };
                ctx.default_second_order_time(&mut ls, b, dt)?;
                ctx.default_update(&ls, b)?;
            }
            ctx.default_finish_boundary_assembly()?;
            let fixed = edge_dirichlet(ctx, &var)?;
            ctx.default_dirichlet_after(&fixed)?;
            ctx.default_solve()?;
            if ctx.convergence().is_some_and(|c| c.nonlin_converged == 1) {
                break;
            }
        }
        Ok(())
    }
}

struct EmWaveCalc;

impl Solver for EmWaveCalc {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        const CALC_CALLER: &str = "EMWaveCalcFields";
        let sp = ctx.params().clone();
        let soln = sp
            .integer("primary solver index")
            .ok_or_else(|| SolverError::fatal(CALC_CALLER, "We should know > Primary Solver Index <"))? as usize;
        let psp = ctx.sim.model.solvers.get(soln - 1).cloned().ok_or_else(|| SolverError::fatal(CALC_CALLER, "no primary solver"))?;
        let (piola, _, degree) = style(&psp);
        if !piola || degree != 1 || ctx.mesh.dim() != 2 {
            return Err(SolverError::Unsupported(format!("{CALC_CALLER}: other than the 2D first-order Piola basis")));
        }
        if psp.logical("eigen analysis") == Some(true) {
            return Err(SolverError::Unsupported(format!("{CALC_CALLER}: Eigen Analysis")));
        }
        for k in ["calculate electric field derivatives", "average within materials"] {
            if sp.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(format!("{CALC_CALLER}: '{k}'")));
            }
        }
        let pix = ctx
            .sim
            .slots
            .iter()
            .position(|s| s.id == soln)
            .ok_or_else(|| SolverError::fatal(CALC_CALLER, "primary solver slot not found"))?;
        let pvar = ctx.sim.slots[pix].variable.ok_or_else(|| SolverError::fatal(CALC_CALLER, "no potential variable"))?;
        if ctx.sim.vars.by_id(pvar).dofs != 1 {
            return Err(SolverError::fatal(CALC_CALLER, "Primary variable should have 1 dofs"));
        }
        let Some(ef_id) = ctx.sim.vars.id("elfield") else { return Ok(()) };
        if ctx.sim.vars.id("elfield e").is_some() {
            return Err(SolverError::Unsupported(format!("{CALC_CALLER}: the elemental fields")));
        }
        ctx.default_initialize()?;
        let n_rows = ctx.system().map_or(0, |s| s.rhs.len());
        let mut gforce = vec![vec![0.0; n_rows]; 3];
        let nbulk = ctx.mesh.n_bulk();
        let active: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < nbulk).collect();
        let mut pslots = Vec::new();
        let mut slots = Vec::new();
        for &e in &active {
            let el = ctx.mesh.element(e);
            let etype = el.etype;
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let gindex: Vec<i64> = el.nodes.iter().map(|&k| k as i64 + 1).collect();
            {
                let d = ctx.sim.slot_dofs(pix).ok_or_else(|| SolverError::fatal(CALC_CALLER, "the potential has no DOF layout"))?;
                d.element_slots(e, &mut pslots);
            }
            let pv = ctx.sim.vars.by_id(pvar);
            let sol: Vec<f64> = pslots.iter().map(|&s| if s == NO_DOF { 0.0 } else { pv.values[s as usize] }).collect();
            ctx.element_slots(e, &mut slots);
            let n = slots.len();
            let rule = bulk_rule(etype.0)?;
            let mut mass = vec![0.0; n * n];
            let mut force = vec![[0.0f64; 3]; n];
            for t in 0..rule.len() {
                let (ev, basis, _) = crate::piola::edge_info_cdim(etype, &xyz, &gindex, [rule.u[t], rule.v[t], rule.w[t]], degree, ctx.mesh.dim())?;
                let nd = ev.dofs.min(sol.len());
                let mut ef = [0.0f64; 3];
                for (c, efc) in ef.iter_mut().enumerate() {
                    let mut acc = 0.0;
                    for i in 0..nd {
                        acc += sol[i] * ev.edge[i][c];
                    }
                    *efc = acc;
                }
                let s = rule.s[t] * ev.det;
                for p in 0..n {
                    for q in 0..n {
                        mass[p * n + q] += s * basis[p] * basis[q];
                    }
                    for c in 0..3 {
                        force[p][c] += s * ef[c] * basis[p];
                    }
                }
            }
            let ls = LocalSystem { n, stiff: mass, mass: None, damp: None, force: force.iter().map(|f| f[0]).collect() };
            ctx.default_update(&ls, e)?;
            for (c, g) in gforce.iter_mut().enumerate() {
                for (p, &s) in slots.iter().enumerate() {
                    if s != NO_DOF {
                        g[s as usize] += force[p][c];
                    }
                }
            }
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        // GlobalSol.
        for (c, g) in gforce.into_iter().enumerate() {
            if let Some(sys) = ctx.system_mut() {
                sys.rhs = g;
            }
            if let Some(v) = ctx.primary_mut() {
                v.values.fill(0.0);
            }
            ctx.default_solve()?;
            let x = ctx.primary().map(|v| v.values.to_vec()).unwrap_or_default();
            let ef = ctx.sim.vars.by_id_mut(ef_id);
            for (i, xi) in x.iter().enumerate() {
                if 3 * i + c < ef.values.len() {
                    ef.values[3 * i + c] = *xi;
                }
            }
            ef.values_changed = true;
        }
        Ok(())
    }
}
