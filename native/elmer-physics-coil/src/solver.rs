//! The `CoilSolver` procedure: `CoilSolver_init` and the solver of
//! CoilSolver.F90:170, for open coils.

use elmer_core::{normalize_key, BasisEval, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

/// `Procedure = "CoilSolver" "CoilSolver"`.
pub const PROCEDURE: (&str, &str) = ("CoilSolver", "CoilSolver");

const CALLER: &str = "CoilSolver";

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(CoilSolver))
}

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("CoilSolver: {what}"))
}

/// `NextFreeKeyword('Exported Variable', Params)`.
fn add_exported(p: &mut ParamMap, spec: &str) {
    let mut k = 1;
    while p.get(&format!("exported variable {k}")).is_some() {
        k += 1;
    }
    p.insert(&format!("Exported Variable {k}"), Value::Str(spec.into()));
}

/// The solver's `Coil Closed` (`Component` sections, which may also set
/// it, are refused before a solver is built).
fn coil_closed(p: &ParamMap) -> bool {
    flag(p, "coil closed")
}

/// The keywords `CoilSolver_init` adds to the solver section: the
/// no-output `CoilTmp` variable, the exported `CoilPot`, `CoilFix`,
/// `CoilPotB`, `PotSelect`, `CoilSet(B)` and the `dim`-component nodal
/// `CoilCurrent`, and `Calculate Loads`. Call it before the variables are
/// created. `dim` is `CoordinateSystemDimension()`. The elemental
/// `CoilCurrent e` and `CoilIndex e` fields are left to the solver to
/// refuse.
pub fn solver_defaults(p: &mut ParamMap, dim: usize) {
    let nodal = p.logical("calculate nodal fields").unwrap_or(true);
    if p.get("variable").is_none() {
        p.insert("Variable", Value::Str("-nooutput CoilTmp".into()));
    }
    add_exported(p, "CoilPot");
    if flag(p, "coil conductivity fix") {
        add_exported(p, "CoilFix");
    }
    if coil_closed(p) {
        add_exported(p, "CoilPotB");
        add_exported(p, "PotSelect");
        if flag(p, "save coil set") {
            add_exported(p, "CoilSetB");
        }
    }
    if flag(p, "save coil set") {
        add_exported(p, "CoilSet");
    }
    if nodal {
        add_exported(p, &format!("-dofs {dim} CoilCurrent"));
    }
    p.insert("Calculate Loads", Value::Logical(true));
}

pub struct CoilSolver;

/// The element's basis and Gauss points, as `ElementInfo` gives them.
pub(crate) struct Point {
    pub weight: f64,
    pub basis: Vec<f64>,
    pub dbasis: Vec<[f64; 3]>,
}

pub(crate) fn element_points(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Vec<Point>, SolverError> {
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
    let el = ctx.mesh.element(e);
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let mut ev = BasisEval::default();
    let mut out = Vec::new();
    for qp in fe.quadrature.gauss_points(el.etype, None) {
        fe.basis
            .eval(el.etype, &xyz, qp.uvw, &mut ev)
            .map_err(|m| SolverError::fatal("ElementInfo", m))?;
        out.push(Point {
            weight: qp.weight * ev.det_j,
            basis: ev.basis.clone(),
            dbasis: ev.dbasis_dx.clone(),
        });
    }
    Ok(out)
}

/// `SUM( a(1:n) * b(1:n) )`, in order.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// `SUM( dBasisdx(1:n, i) * v(1:n) )`.
fn grad_i(db: &[[f64; 3]], i: usize, v: &[f64]) -> f64 {
    let mut s = 0.0;
    for (d, x) in db.iter().zip(v) {
        s += d[i] * x;
    }
    s
}

/// `SQRT( SUM( v**2 ) )`.
fn norm(v: &[f64]) -> f64 {
    let mut s = 0.0;
    for x in v {
        s += x * x;
    }
    s.sqrt()
}

/// The field values at the element's nodes through the solver's Perm.
fn nodal(values: &[f64], perm: &[u32], nodes: &[u32]) -> Vec<f64> {
    nodes
        .iter()
        .map(|&k| values[perm[k as usize] as usize])
        .collect()
}

/// What the local matrices read of the solver's options.
struct Opts {
    fix_conductivity: bool,
    anisotropic: bool,
    cond_name: String,
    unity: bool,
    dim: usize,
}

impl Opts {
    /// The conductivity vector at a Gauss point: `CondAtIp` of
    /// LocalPotMatrix (`fix` the `CoilFix` values when they apply).
    fn cond(
        &self,
        pt: &Point,
        el_cond: Option<&[f64]>,
        fix: Option<&[f64]>,
        pot: &[f64],
    ) -> [f64; 3] {
        let mut cond = [1.0f64; 3];
        let fix_at = match fix {
            Some(f) => {
                let fix_at = dot(&pt.basis, f);
                if self.anisotropic {
                    let grad = [
                        grad_i(&pt.dbasis, 0, pot),
                        grad_i(&pt.dbasis, 1, pot),
                        grad_i(&pt.dbasis, 2, pot),
                    ];
                    let abs = norm(&grad);
                    for i in 0..3 {
                        cond[i] = grad[i].abs() / abs;
                    }
                }
                fix_at
            }
            None => 1.0,
        };
        match el_cond {
            Some(c) => {
                let s = dot(&pt.basis, c) * fix_at;
                for c in &mut cond {
                    *c *= s;
                }
            }
            None => {
                for c in &mut cond {
                    *c *= fix_at;
                }
            }
        }
        cond
    }
}

/// `GetReal( Material, CondName, GotElCond )` unless `Use Unity
/// Conductivity`.
fn el_cond(ctx: &SolverCtx<'_>, o: &Opts, e: ElemId) -> Result<Option<Vec<f64>>, SolverError> {
    if o.unity {
        return Ok(None);
    }
    let Some(m) = ctx.body_params(e).and_then(|b| b.material) else {
        return Ok(None);
    };
    ctx.real_at_nodes(m, &o.cond_name, e)
}

/// `LocalPotMatrix`: `div( sigma grad u ) = 0` with the conductivity of
/// [`Opts::cond`] (the `CoilFix` values from the second iteration on).
fn local_pot(
    ctx: &SolverCtx<'_>,
    o: &Opts,
    e: ElemId,
    pot: &[f64],
    fix: Option<&[f64]>,
) -> Result<LocalSystem, SolverError> {
    let n = ctx.mesh.element(e).nodes.len();
    let elc = el_cond(ctx, o, e)?;
    let mut ls = LocalSystem::zeros(n);
    for pt in element_points(ctx, e)? {
        let cond = o.cond(&pt, elc.as_deref(), fix, pot);
        let w = pt.weight;
        for p in 0..n {
            for q in 0..n {
                let mut s = ls.stiff[p * n + q];
                for i in 0..o.dim {
                    s += w * cond[i] * pt.dbasis[p][i] * pt.dbasis[q][i];
                }
                ls.stiff[p * n + q] = s;
            }
        }
    }
    Ok(ls)
}

/// `LocalFixMatrix`: `-div( D grad c ) + |grad u| c = 1` for the
/// conductivity fix `c` (`CFix Diffusion` D).
fn local_fix(
    ctx: &SolverCtx<'_>,
    e: ElemId,
    pot: &[f64],
    diff: &[f64],
) -> Result<LocalSystem, SolverError> {
    let n = ctx.mesh.element(e).nodes.len();
    let mut ls = LocalSystem::zeros(n);
    for pt in element_points(ctx, e)? {
        let grad = [
            grad_i(&pt.dbasis, 0, pot),
            grad_i(&pt.dbasis, 1, pot),
            grad_i(&pt.dbasis, 2, pot),
        ];
        let abs = norm(&grad);
        let dreg = dot(&pt.basis, diff);
        let w = pt.weight;
        let wd = w * dreg;
        for p in 0..n {
            for q in 0..n {
                // MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
                let mut m = 0.0;
                for k in 0..3 {
                    m += pt.dbasis[p][k] * pt.dbasis[q][k];
                }
                ls.stiff[p * n + q] += wd * m;
            }
        }
        for p in 0..n {
            for q in 0..n {
                ls.stiff[p * n + q] += w * abs * pt.basis[q] * pt.basis[p];
            }
        }
        for p in 0..n {
            ls.force[p] += w * 1.0 * pt.basis[p];
        }
    }
    Ok(ls)
}

/// `LocalFluxMatrix`: the mass matrix and the component `dimi` of `sigma
/// grad u`.
fn local_flux(
    ctx: &SolverCtx<'_>,
    o: &Opts,
    e: ElemId,
    pot: &[f64],
    fix: Option<&[f64]>,
    dimi: usize,
) -> Result<LocalSystem, SolverError> {
    let n = ctx.mesh.element(e).nodes.len();
    let elc = el_cond(ctx, o, e)?;
    let mut ls = LocalSystem::zeros(n);
    for pt in element_points(ctx, e)? {
        let grad = [
            grad_i(&pt.dbasis, 0, pot),
            grad_i(&pt.dbasis, 1, pot),
            grad_i(&pt.dbasis, 2, pot),
        ];
        let cond = o.cond(&pt, elc.as_deref(), fix, pot);
        let load = cond[dimi] * grad[dimi];
        let w = pt.weight;
        for p in 0..n {
            for q in 0..n {
                ls.stiff[p * n + q] += w * pt.basis[q] * pt.basis[p];
            }
        }
        for p in 0..n {
            ls.force[p] += w * load * pt.basis[p];
        }
    }
    Ok(ls)
}

impl Solver for CoilSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        for k in [
            "single coil cut",
            "test coil cut",
            "calculate elemental fields",
            "fix input current density",
            "use wall distance",
            "save coil index",
        ] {
            if flag(p, k) {
                return Err(unsupported(k));
            }
        }
        if p.get("coil helicity").is_some() {
            return Err(unsupported("Coil Helicity"));
        }
        if p.logical("calculate nodal fields") == Some(false) {
            return Err(unsupported("Calculate Nodal Fields = False"));
        }
        let cs = ctx
            .simulation_params()
            .string("coordinate system")
            .map(normalize_key)
            .unwrap_or_else(|| "cartesian".into());
        if !cs.starts_with("cartesian") {
            return Err(unsupported("non-Cartesian coordinates"));
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "--------------------------------------", 4);
        ctx.info(CALLER, "Solving current distribution in a coil", 4);
        ctx.info(CALLER, "--------------------------------------", 4);
        let p = ctx.params().clone();
        let Some(prim) = ctx.primary() else {
            return Err(SolverError::fatal(CALLER, "no variable"));
        };
        let sol_name = prim.name.to_string();
        let nsize = prim.values.len();
        let perm: Vec<u32> = prim
            .perm
            .ok_or_else(|| SolverError::fatal(CALLER, "the variable has no Perm"))?
            .to_vec();
        let dim = ctx.mesh.dim();
        let o = Opts {
            fix_conductivity: flag(&p, "coil conductivity fix"),
            anisotropic: flag(&p, "coil anisotropic"),
            cond_name: p
                .string("electric conductivity name")
                .unwrap_or("Electric Conductivity")
                .to_string(),
            unity: flag(&p, "use unity conductivity"),
            dim,
        };
        let closed = coil_closed(&p);
        let parts = if closed { 2 } else { 1 };
        let bcs = &ctx.sim.model.boundary_conditions;
        if closed {
            ctx.info(CALLER, "Assuming that all coils are closed!", 5);
        } else if !bcs.iter().any(|b| flag(b, "coil start")) {
            return Err(SolverError::fatal(
                CALLER,
                "> Coil Start < must be defined on some BC",
            ));
        }
        if !closed && !bcs.iter().any(|b| flag(b, "coil end")) {
            return Err(SolverError::fatal(
                CALLER,
                "> Coil End < must be defined on some BC",
            ));
        }
        let calc_curr = p.logical("calculate coil current").unwrap_or(true);
        let max_iter = match (o.fix_conductivity, o.anisotropic) {
            (true, true) => 3,
            (true, false) => 2,
            _ => 1,
        };
        if ctx.variable("coilpot").is_none() {
            return Err(SolverError::fatal(CALLER, "CoilPot not associated!"));
        }
        if closed {
            for v in ["CoilPotB", "PotSelect"] {
                if ctx.variable(v).is_none() {
                    return Err(SolverError::fatal(CALLER, format!("{v} not associated!")));
                }
            }
        }
        if o.fix_conductivity {
            let f = ctx
                .variable_mut("coilfix")
                .ok_or_else(|| SolverError::fatal(CALLER, "CoilFix not associated!"))?;
            f.values.fill(1.0);
        }
        let loads_name = format!("{sol_name} loads");
        if ctx.variable(&loads_name).is_none() {
            return Err(SolverError::fatal(
                CALLER,
                format!("> {sol_name} < Loads not associated!"),
            ));
        }

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut set_a = vec![0i32; nsize];
        let mut set_b = vec![0i32; if closed { nsize } else { 0 }];
        if closed {
            let center = crate::fit::coil_center(ctx, &p, &elems)?;
            let axes = crate::fit::coil_axes(ctx, &p, &elems, center)?;
            choose_fixed_bulk_nodes(ctx, &perm, &mut set_a, 1, center, &axes)?;
            choose_fixed_bulk_nodes(ctx, &perm, &mut set_b, 2, center, &axes)?;
        } else {
            choose_fixed_end_nodes(ctx, &perm, &mut set_a);
        }
        let desired_curr = ctx.creal(&p, "desired coil current")?;
        let got_curr = desired_curr.is_some();
        let mut desired_curr = desired_curr.unwrap_or(1.0);
        let desired_dens = ctx.creal(&p, "desired current density")?;
        let mut got_dens = desired_dens.is_some();
        let mut desired_dens = desired_dens.unwrap_or(1.0);
        let mut got_curr = got_curr;
        let normalize = flag(&p, "normalize coil current");
        if let Some(cross) = ctx.creal(&p, "coil cross section")? {
            if got_curr {
                desired_dens = desired_curr / cross;
                got_dens = true;
            } else if got_dens {
                desired_curr = desired_dens * cross;
                got_curr = true;
            }
        }
        ctx.info(CALLER, "Coil system consists of 1 coils", 7);
        for (k, set) in [&set_a, &set_b].into_iter().take(parts).enumerate() {
            // CountFixingNodes: in parallel the owned nodes, summed.
            let mut n = [0.0f64; 2];
            for (i, &j) in perm.iter().enumerate() {
                if j == NO_DOF || ctx.node_owned().is_some_and(|o| !o[i]) {
                    continue;
                }
                match set[j as usize] {
                    s if s > 0 => n[0] += 1.0,
                    s if s < 0 => n[1] += 1.0,
                    _ => {}
                }
            }
            ctx.comm.sum_f64(&mut n);
            let (nplus, nminus) = (n[0] as usize, n[1] as usize);
            ctx.info(
                CALLER,
                &format!("Set{} : {nplus} +nodes and {nminus} -nodes", k + 1),
                4,
            );
            if nplus == 0 || nminus == 0 {
                ctx.warn(CALLER, "Cannot set Dirichlet conditions with this set");
            }
        }

        let pot_of = |ctx: &SolverCtx<'_>, name: &str| -> Vec<f64> {
            ctx.variable(name)
                .map(|v| v.values.to_vec())
                .unwrap_or_default()
        };
        // NodalPot: CoilPotB on the elements PotSelect puts wholly on the
        // positive side of a closed coil, CoilPot elsewhere.
        let element_pot = |ctx: &SolverCtx<'_>, e: ElemId| -> Vec<f64> {
            let nodes = ctx.mesh.element(e).nodes;
            let name = if closed {
                let sel = ctx.variable("potselect").expect("checked");
                let sp = sel.perm.expect("exported on the solver's Perm");
                let min = nodes
                    .iter()
                    .map(|&k| sel.values[sp[k as usize] as usize])
                    .fold(f64::INFINITY, f64::min);
                if min > 0.0 {
                    "coilpotb"
                } else {
                    "coilpot"
                }
            } else {
                "coilpot"
            };
            let v = ctx.variable(name).expect("checked");
            nodal(v.values, &perm, nodes)
        };
        for iter in 1..=max_iter {
            if iter > 1 {
                ctx.info(CALLER, "Fixing the conductivity field", 4);
                ctx.default_initialize()?;
                for &e in &elems {
                    let n = ctx.mesh.element(e).nodes.len();
                    let diff = ctx
                        .real_at_nodes(&p, "cfix diffusion", e)?
                        .unwrap_or_else(|| vec![0.0; n]);
                    let ls = local_fix(ctx, e, &element_pot(ctx, e), &diff)?;
                    ctx.default_update(&ls, e)?;
                }
                ctx.set_param("Skip Compute Nonlinear Change", Value::Logical(true));
                // The loads of this solve are those of the stale bulk system
                // in the Fortran and are overwritten by the next one.
                ctx.set_param("Calculate Loads", Value::Logical(false));
                ctx.default_solve()?;
                ctx.set_param("Calculate Loads", Value::Logical(true));
                let sol = pot_of(ctx, &sol_name);
                ctx.variable_mut("coilfix")
                    .expect("checked")
                    .values
                    .copy_from_slice(&sol);
            }
            ctx.info(CALLER, "Computing the dummy potential field", 4);
            for part in 0..parts {
                let (set, pot_name) = if part == 0 {
                    (&set_a, "coilpot")
                } else {
                    (&set_b, "coilpotb")
                };
                ctx.default_initialize()?;
                let fix = if iter > 1 {
                    Some(pot_of(ctx, "coilfix"))
                } else {
                    None
                };
                for &e in &elems {
                    let nodes = ctx.mesh.element(e).nodes.to_vec();
                    let fx = fix.as_ref().map(|f| nodal(f, &perm, &nodes));
                    let ls = local_pot(ctx, &o, e, &element_pot(ctx, e), fx.as_deref())?;
                    ctx.default_update(&ls, e)?;
                }
                ctx.default_finish_bulk_assembly()?;
                let fixed: Vec<(u32, f64)> = set
                    .iter()
                    .enumerate()
                    .filter(|(_, &s)| s != 0)
                    .map(|(i, &s)| (i as u32, if s > 0 { 1.0 } else { 0.0 }))
                    .collect();
                cut_interface_connections(ctx, set)?;
                ctx.default_dirichlet_after(&fixed)?;
                ctx.set_param("Skip Compute Nonlinear Change", Value::Logical(false));
                ctx.default_solve()?;
                let sol = pot_of(ctx, &sol_name);
                ctx.variable_mut(pot_name)
                    .expect("checked")
                    .values
                    .copy_from_slice(&sol);
                scale_potential(ctx, set, desired_curr, &loads_name, pot_name)?;
            }
        }

        if calc_curr {
            ctx.info(CALLER, "Computing current densities", 20);
            ctx.set_param("Calculate Loads", Value::Logical(false));
            let fix = if o.fix_conductivity {
                Some(pot_of(ctx, "coilfix"))
            } else {
                None
            };
            let mut min_curr = [0.0f64; 3];
            let mut max_curr = [0.0f64; 3];
            // The Fortran repeats the components for each part; the
            // repeat restarts the iterative solves from the last solution.
            for (part, dimi) in (0..parts).flat_map(|p| (0..dim).map(move |d| (p, d))) {
                let _ = part;
                ctx.info(
                    CALLER,
                    &format!("Computing current component: {}", dimi + 1),
                    6,
                );
                if ctx.variable("coilcurrent").is_none() {
                    return Err(SolverError::fatal(CALLER, "CoilCurrent not associated!"));
                }
                ctx.default_initialize()?;
                for &e in &elems {
                    let nodes = ctx.mesh.element(e).nodes.to_vec();
                    let fx = fix.as_ref().map(|f| nodal(f, &perm, &nodes));
                    let ls = local_flux(ctx, &o, e, &element_pot(ctx, e), fx.as_deref(), dimi)?;
                    ctx.default_update(&ls, e)?;
                }
                ctx.set_param("Skip Compute Nonlinear Change", Value::Logical(true));
                ctx.default_solve()?;
                let sol = pot_of(ctx, &sol_name);
                min_curr[dimi] = sol.iter().copied().fold(f64::INFINITY, f64::min);
                max_curr[dimi] = sol.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let cur = ctx.variable_mut("coilcurrent").expect("checked");
                for (slot, v) in sol.iter().enumerate() {
                    cur.values[slot * dim + dimi] = *v;
                }
            }
            if ctx.comm.size() > 1 {
                let mut neg = min_curr.map(|x| -x);
                ctx.comm.max_f64(&mut neg);
                min_curr = neg.map(|x| -x);
                ctx.comm.max_f64(&mut max_curr);
            }
            ctx.info(
                CALLER,
                &format!("Minimum current components: {min_curr:?}"),
                7,
            );
            ctx.info(
                CALLER,
                &format!("Maximum current components: {max_curr:?}"),
                7,
            );
            if normalize {
                normalize_current_density(ctx, &elems, dim, got_dens, got_curr, desired_dens)?;
            } else {
                ctx.info(CALLER, "No normalization of current requested!", 20);
            }
            ctx.set_param("Calculate Loads", Value::Logical(true));
        }

        // `CoilSet` (and `CoilSetB`): the fixing sets by Perm entry.
        for (name, set) in [("coilset", &set_a), ("coilsetb", &set_b)]
            .into_iter()
            .take(parts)
        {
            if let Some(cs) = ctx.variable_mut(name) {
                if let Some(cperm) = cs.perm {
                    let cperm = cperm.to_vec();
                    for &j in &cperm {
                        if j != NO_DOF {
                            cs.values[j as usize] = f64::from(set[j as usize]);
                        }
                    }
                }
            }
        }
        let pot = pot_of(ctx, "coilpot");
        ctx.primary_mut()
            .expect("checked")
            .values
            .copy_from_slice(&pot);
        ctx.info(CALLER, "All done", 7);
        Ok(())
    }
}

/// `ChooseFixedEndNodes`: `+12` on `Coil Start`, `-12` on `Coil End`
/// (`2 + ioffset`, `ioffset = 10 * NoCoils`).
fn choose_fixed_end_nodes(ctx: &SolverCtx<'_>, perm: &[u32], set: &mut [i32]) {
    let bnd: Vec<ElemId> = ctx.boundary_elements().collect();
    for &e in &bnd {
        let nodes = ctx.mesh.element(e).nodes;
        if nodes.iter().any(|&k| perm[k as usize] == NO_DOF) {
            continue;
        }
        let Some(bc) = ctx.bc_params(e) else { continue };
        let v = if flag(bc, "coil start") {
            12
        } else if flag(bc, "coil end") {
            -12
        } else {
            continue;
        };
        for &k in nodes {
            set[perm[k as usize] as usize] = v;
        }
    }
}

/// `ChooseFixedBulkNodesNarrow( Set, SetNo, .FALSE. )` of a closed coil:
/// the bulk elements on the `tangent1 >= 0` half (mirrored for set 2) that
/// the plane `tangent2 = 0` cuts are fixed, `+12` on its positive side and
/// `-12` on the other. Set 1 also stores `tangent1 . r` in `PotSelect`.
fn choose_fixed_bulk_nodes(
    ctx: &mut SolverCtx<'_>,
    perm: &[u32],
    set: &mut [i32],
    set_no: usize,
    center: [f64; 3],
    axes: &crate::fit::Axes,
) -> Result<(), SolverError> {
    ctx.info(
        CALLER,
        &format!("Choosing fixing nodes for set: {set_no}"),
        4,
    );
    let mirror = set_no == 2;
    let sp: Vec<u32> = ctx
        .variable("potselect")
        .and_then(|v| v.perm.map(<[u32]>::to_vec))
        .ok_or_else(|| SolverError::fatal(CALLER, "PotSelect not associated!"))?;
    let dot3 = |a: &[f64; 3], b: &[f64; 3]| {
        let mut s = 0.0;
        for k in 0..3 {
            s += a[k] * b[k];
        }
        s
    };
    for t in 0..ctx.mesh.n_bulk() {
        let nodes = ctx.mesh.element(ElemId(t as u32)).nodes.to_vec();
        let mut wrong_side = false;
        let mut cut: Vec<f64> = Vec::with_capacity(nodes.len());
        for &i in &nodes {
            if perm[i as usize] == NO_DOF {
                continue;
            }
            let x = ctx.mesh.coord(i);
            let mut r = [x[0] - center[0], x[1] - center[1], x[2] - center[2]];
            if mirror {
                r = [-r[0], -r[1], -r[2]];
            }
            let rp1 = dot3(&axes.tangent1, &r);
            let rp2 = dot3(&axes.tangent2, &r);
            if set_no == 1 {
                let sel = ctx.variable_mut("potselect").expect("checked");
                sel.values[sp[i as usize] as usize] = rp1;
            }
            if rp1 < 0.0 {
                wrong_side = true;
                break;
            }
            cut.push(rp2);
        }
        if wrong_side || cut.len() < 2 {
            continue;
        }
        let max = cut.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let min = cut.iter().copied().fold(f64::INFINITY, f64::min);
        if max >= -f64::EPSILON && min <= f64::EPSILON {
            let mut l = 0;
            for &i in &nodes {
                let j = perm[i as usize];
                if j == NO_DOF {
                    continue;
                }
                set[j as usize] = if cut[l] > 0.0 { 12 } else { -12 };
                l += 1;
            }
        }
    }
    if ctx.comm.size() > 1 {
        // The holders of a node agree on the value of greatest magnitude,
        // the positive one on a tie (`ParallelSumVectorInt` max and min).
        let of_nodes = |set: &[i32]| -> Vec<f64> {
            perm.iter().map(|&j| if j == NO_DOF { 0.0 } else { f64::from(set[j as usize]) }).collect()
        };
        let mut max = of_nodes(set);
        let mut min = max.clone();
        ctx.nodes_max(&mut max)?;
        ctx.nodes_min(&mut min)?;
        for (i, &j) in perm.iter().enumerate() {
            if j != NO_DOF {
                set[j as usize] = if min[i].abs() > max[i] { min[i] } else { max[i] } as i32;
            }
        }
    }
    Ok(())
}

/// `CutInterfaceConnections( A, Set )`: the bulk couplings between the two
/// fixing sets move onto the diagonal, so that the nodal loads of each end
/// see only its own side.
fn cut_interface_connections(ctx: &mut SolverCtx<'_>, set: &[i32]) -> Result<(), SolverError> {
    let (a, bulk) = ctx
        .bulk_values_mut()
        .ok_or_else(|| SolverError::fatal("CutInterfaceConnections", "Do not have bulk values!"))?;
    let mut cnt = 0usize;
    for i in 0..a.nrows {
        let s1 = set[i];
        if s1 == 0 {
            continue;
        }
        let r = a.row_ptr[i]..a.row_ptr[i + 1];
        let dia = match &a.diag {
            Some(d) => d[i],
            None => r
                .clone()
                .find(|&j| a.col_idx[j] as usize == i)
                .ok_or_else(|| {
                    SolverError::fatal("CutInterfaceConnections", "no diagonal entry")
                })?,
        };
        for j in r {
            let s2 = set[a.col_idx[j] as usize];
            if s1 * s2 < 0 {
                bulk[dia] += bulk[j];
                bulk[j] = 0.0;
                cnt += 1;
            }
        }
    }
    let mut c = [cnt as f64];
    ctx.comm.sum_f64(&mut c);
    let cnt = c[0] as usize;
    ctx.info(
        "CutInterfaceConnections",
        &format!("Number of connections cut: {cnt}"),
        7,
    );
    Ok(())
}

/// `ScalePotential`: the coil current from the nodal loads of the two ends
/// (in node order), and `CoilPot` scaled to the desired current.
fn scale_potential(
    ctx: &mut SolverCtx<'_>,
    set: &[i32],
    desired: f64,
    loads_name: &str,
    pot_name: &str,
) -> Result<(), SolverError> {
    let (possum, negsum) = {
        let lv = ctx.variable(loads_name).expect("checked");
        let lperm = lv
            .perm
            .ok_or_else(|| SolverError::fatal(CALLER, "loads without Perm"))?;
        let (mut possum, mut negsum) = (0.0, 0.0);
        let owned = ctx.node_owned();
        for (i, &j) in lperm.iter().enumerate() {
            if j == NO_DOF || owned.is_some_and(|o| !o[i]) {
                continue;
            }
            let sgn = set[j as usize];
            if sgn == 0 {
                continue;
            }
            if sgn.rem_euclid(10) == 2 {
                possum += lv.values[j as usize];
            } else if sgn.rem_euclid(10) - 10 == -2 {
                negsum += lv.values[j as usize];
            }
        }
        let mut s = [possum, negsum];
        if ctx.comm.size() > 1 {
            ctx.comm.sum_f64(&mut s);
        }
        (s[0], s[1])
    };
    let mut ok = true;
    if possum.abs() < f64::EPSILON {
        ctx.warn(CALLER, "No positive current sources on coil 1 end!");
        ok = false;
    }
    if negsum.abs() < f64::EPSILON {
        ctx.warn(CALLER, "No negative current sources on coil 1 end!");
        ok = false;
    }
    if !ok {
        ctx.warn(CALLER, "Crappy potentials in coil 1");
        return Err(SolverError::fatal(CALLER, "Scaling of potential failed!"));
    }
    let sumerr = 2.0 * (possum.abs() - negsum.abs()).abs() / (possum.abs() + negsum.abs());
    if sumerr > 0.5 {
        ctx.warn(CALLER, "Positive and negative sums differ quite a bit!");
    }
    let initial = (possum - negsum) / 2.0;
    let coeff = desired / initial;
    ctx.info(CALLER, &format!("Coil potential multiplier: {coeff:e}"), 5);
    let pv = ctx.variable_mut(pot_name).expect("checked");
    for v in pv.values.iter_mut() {
        *v *= coeff;
    }
    Ok(())
}

/// `NormalizeCurrentDensity` for the nodal `CoilCurrent` of one coil: each
/// node's current scaled to the desired density, or to the average density
/// when only the current is given.
fn normalize_current_density(
    ctx: &mut SolverCtx<'_>,
    elems: &[ElemId],
    dim: usize,
    got_dens: bool,
    got_curr: bool,
    desired_dens: f64,
) -> Result<(), SolverError> {
    ctx.info(
        CALLER,
        "Normalizing current density (as for stranded coil)",
        5,
    );
    let (values, cperm) = {
        let f = ctx.variable("coilcurrent").ok_or_else(|| {
            SolverError::fatal(CALLER, "CoilCurrent (nodal and elemental) not associated!")
        })?;
        let cperm = f
            .perm
            .ok_or_else(|| SolverError::fatal(CALLER, "CoilCurrent without Perm"))?;
        (f.values.to_vec(), cperm.to_vec())
    };
    let target = if !got_dens && got_curr {
        let (mut tot_vol, mut tot_curr) = (0.0, 0.0);
        for &e in elems {
            let nodes = ctx.mesh.element(e).nodes.to_vec();
            let n = nodes.len();
            let mut nc = vec![vec![0.0; n]; dim];
            for d in 0..dim {
                for (i, &k) in nodes.iter().enumerate() {
                    nc[d][i] = values[dim * cperm[k as usize] as usize + d];
                }
            }
            for pt in element_points(ctx, e)? {
                let mut local = [0.0f64; 3];
                for d in 0..dim {
                    local[d] = dot(&pt.basis, &nc[d]);
                }
                let abs = norm(&local[..dim]);
                tot_vol += pt.weight;
                tot_curr += pt.weight * abs;
            }
        }
        let t = tot_curr / tot_vol;
        ctx.info(CALLER, &format!("Average current density: {t:e}"), 4);
        t
    } else {
        desired_dens
    };
    let n_nodes = ctx.mesh.n_nodes();
    let f = ctx.variable_mut("coilcurrent").expect("checked");
    for i in 0..n_nodes {
        let j = cperm[i];
        if j == NO_DOF {
            continue;
        }
        let j = j as usize;
        let local: Vec<f64> = f.values[dim * j..dim * j + dim].to_vec();
        let abs = norm(&local);
        if abs > f64::MIN_POSITIVE {
            let scale = target / abs;
            for d in 0..dim {
                f.values[dim * j + d] = scale * local[d];
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(p: &ParamMap, k: &str) -> Option<String> {
        p.string(k).map(str::to_string)
    }

    #[test]
    fn init_adds_the_coil_fields() {
        let mut p = ParamMap::new();
        p.insert("Coil Conductivity Fix", Value::Logical(true));
        p.insert("Save Coil Set", Value::Logical(true));
        solver_defaults(&mut p, 3);
        assert_eq!(s(&p, "variable").as_deref(), Some("-nooutput CoilTmp"));
        let exported: Vec<String> = (1..=4)
            .filter_map(|k| s(&p, &format!("exported variable {k}")))
            .collect();
        assert_eq!(
            exported,
            ["CoilPot", "CoilFix", "CoilSet", "-dofs 3 CoilCurrent"]
        );
        assert_eq!(p.logical("calculate loads"), Some(true));
    }
}
