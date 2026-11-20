//! The `TopoOpt` procedure: `TopoOpt_init0`, `TopoOpt_init` and the solver
//! of TopoOpt.F90:107.

use elmer_core::{normalize_key, BasisEval, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::filter::{apply_filter, distance_filter, simple_filter, Crs};

/// `Procedure = "TopoOpt" "TopoOpt"`.
pub const PROCEDURE: (&str, &str) = ("TopoOpt", "TopoOpt");

const CALLER: &str = "TopoOpt";

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(TopoOpt::default()))
}

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("TopoOpt: {what}"))
}

fn add_exported(p: &mut ParamMap, spec: &str) {
    let mut k = 1;
    while p.get(&format!("exported variable {k}")).is_some() {
        k += 1;
    }
    p.insert(&format!("Exported Variable {k}"), Value::Str(spec.into()));
}

fn is_pde(p: &ParamMap) -> bool {
    p.string("filter type")
        .is_some_and(|s| normalize_key(s) == "pde")
}

/// `TopoOpt_init0` and `TopoOpt_init`: `Primary Solver Index`, the
/// elemental `topo rho`, `topo mult`, `topo ce`, `topo dc` (and `topo bw`)
/// fields and the variable: the nodal `xNodal` of the PDE filter, else the
/// global `topoopt_nrm` that carries the objective. Call it before the
/// variables are created.
pub fn solver_defaults(p: &mut ParamMap) {
    if p.get("primary solver index").is_none() {
        p.insert("Primary Solver Index", Value::Integer(vec![1]));
    }
    for f in ["topo rho", "topo mult", "topo ce", "topo dc"] {
        add_exported(p, &format!("-elem {f}"));
    }
    if flag(p, "create bw topology") {
        add_exported(p, "-elem topo bw");
    }
    if p.get("variable").is_none() {
        let v = if is_pde(p) {
            "xNodal"
        } else {
            "-nooutput -global topoopt_nrm"
        };
        p.insert("Variable", Value::Str(v.into()));
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Method {
    Sensitivity,
    Density,
    None,
}

/// What `TopoOpt` keeps between visits (its `SAVE` variables).
#[derive(Default)]
pub struct TopoOpt {
    times_visited: usize,
    state: Option<State>,
    /// `DoMatrix` of `ApplyPDEFilter`: the PDE filter's matrix is
    /// assembled on the first call only.
    pde_matrix_done: bool,
}

struct State {
    method: Method,
    fmat: Option<Crs>,
    niter: usize,
    pde_diff: Option<[f64; 3]>,
    interface_node: Option<Vec<bool>>,
    phys: usize,
    dofs: usize,
    x_topo: Vec<f64>,
    dv0: Vec<f64>,
    dv: Vec<f64>,
    wmin: f64,
}

/// `ElementSize( Element, Nodes )`: `detJ` at the family's fixed point.
fn element_size(ctx: &SolverCtx<'_>, e: ElemId) -> Result<f64, SolverError> {
    let el = ctx.mesh.element(e);
    let uvw = match el.etype.0 / 100 {
        1 => return Ok(1.0),
        2 | 4 | 6 | 8 => [0.0, 0.0, 0.0],
        3 | 7 => [0.5, 0.5, 0.0],
        5 => [0.5, 0.5, 0.5],
        _ => {
            return Err(SolverError::fatal(
                "ElementSize",
                "Not implemented for elementtype",
            ))
        }
    };
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let mut ev = BasisEval::default();
    fe.basis
        .eval(el.etype, &xyz, uvw, &mut ev)
        .map_err(|m| SolverError::fatal("ElementInfo", m))?;
    Ok(ev.det_j)
}

/// The values and permutation of an elemental field.
fn elemental(ctx: &SolverCtx<'_>, name: &str) -> Result<(Vec<f64>, Vec<u32>), SolverError> {
    let v = ctx
        .variable(name)
        .ok_or_else(|| SolverError::fatal(CALLER, format!("could not find \"{name}\"")))?;
    let perm = v
        .perm
        .map(<[u32]>::to_vec)
        .unwrap_or_else(|| (0..v.values.len() as u32).collect());
    Ok((v.values.to_vec(), perm))
}

fn set_values(ctx: &mut SolverCtx<'_>, name: &str, values: &[f64]) {
    if let Some(v) = ctx.variable_mut(name) {
        v.values.copy_from_slice(values);
    }
}

impl TopoOpt {
    /// `ApplyPDEFilter( x, xf, Diff, w, wmin )`: the Helmholtz-type filter
    /// `(N, u) + (D grad N, grad u) = (N, x)` on the `xNodal` field, `xf`
    /// the element means of `u` by active element.
    fn pde_filter(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        st: &State,
        x: &[f64],
        xf: &mut [f64],
        w: Option<(&[f64], f64)>,
    ) -> Result<(), SolverError> {
        let diff = st.pde_diff.expect("PDE filter");
        let dim = ctx.mesh.dim();
        let xtmp: Vec<f64> = match w {
            Some((w, _)) => x.iter().zip(w).map(|(a, b)| a * b).collect(),
            None => x.to_vec(),
        };
        let do_matrix = !self.pde_matrix_done;
        if do_matrix {
            ctx.default_initialize()?;
        } else if let Some(s) = ctx.system_mut() {
            s.rhs.iter_mut().for_each(|v| *v = 0.0);
        }
        let (_, eperm) = elemental(ctx, "topo rho")?;
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        for &e in &elems {
            let (etype, nodes) = {
                let el = ctx.mesh.element(e);
                (el.etype, el.nodes.to_vec())
            };
            let n = nodes.len();
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let fe = ctx.sim.backends.fe.as_ref().expect("checked in init");
            let mut ls = LocalSystem::zeros(n);
            let mut ev = BasisEval::default();
            if do_matrix {
                for qp in fe.quadrature.gauss_points(etype, None) {
                    fe.basis
                        .eval(etype, &xyz, qp.uvw, &mut ev)
                        .map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let weight = qp.weight * ev.det_j;
                    for p in 0..n {
                        for q in 0..n {
                            ls.stiff[p * n + q] += weight * ev.basis[p] * ev.basis[q];
                        }
                    }
                    for p in 0..n {
                        if st
                            .interface_node
                            .as_ref()
                            .is_some_and(|f| f[nodes[p] as usize])
                        {
                            continue;
                        }
                        for q in 0..n {
                            let mut s = 0.0;
                            for k in 0..dim {
                                s += diff[k] * ev.dbasis_dx[p][k] * ev.dbasis_dx[q][k];
                            }
                            ls.stiff[p * n + q] += weight * s;
                        }
                    }
                }
            }
            let xi = xtmp[eperm[e.0 as usize] as usize];
            for qp in fe.quadrature.gauss_points(etype, Some(1)) {
                fe.basis
                    .eval(etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let weight = qp.weight * ev.det_j;
                for p in 0..n {
                    ls.force[p] += weight * ev.basis[p] * xi;
                }
            }
            ctx.default_update(&ls, e)?;
        }
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        let v = ctx.primary().expect("checked in init");
        let perm = v.perm.expect("nodal");
        for (t, &e) in elems.iter().enumerate() {
            let nodes = ctx.mesh.element(e).nodes;
            let mut s = 0.0;
            for &k in nodes {
                s += v.values[perm[k as usize] as usize];
            }
            xf[t] = s / nodes.len() as f64;
        }
        if let Some((w, wmin)) = w {
            for (v, wi) in xf.iter_mut().zip(w) {
                *v /= wi.max(wmin);
            }
        }
        self.pde_matrix_done = true;
        Ok(())
    }

    /// The filter of `st`: PDE or matrix, plain or transposed or weighted.
    fn filter(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        st: &State,
        x: &[f64],
        xf: &mut [f64],
        w: Option<(&[f64], f64)>,
        trans: bool,
    ) -> Result<(), SolverError> {
        match &st.fmat {
            Some(f) => {
                apply_filter(f, x, xf, st.niter, w, trans);
                Ok(())
            }
            None => self.pde_filter(ctx, st, x, xf, w),
        }
    }

    fn first_visit(&mut self, ctx: &mut SolverCtx<'_>) -> Result<State, SolverError> {
        let p = ctx.params().clone();
        let mesh = ctx.mesh;
        ctx.info(CALLER, &format!("Number of nodes: {}", mesh.n_nodes()), 4);
        ctx.info(
            CALLER,
            &format!("Number of bulk elements: {}", mesh.n_bulk()),
            4,
        );
        let method = match p.string("filter method").map(normalize_key).as_deref() {
            Some("sensitivity") => Method::Sensitivity,
            Some("density") => Method::Density,
            Some("none") => Method::None,
            Some(m) => {
                return Err(SolverError::fatal(
                    CALLER,
                    format!("Uknown filtering method: {m}"),
                ))
            }
            None => return Err(SolverError::fatal(CALLER, "> Filter Method < not found")),
        };
        let ftype = p
            .string("filter type")
            .map(normalize_key)
            .ok_or_else(|| SolverError::fatal(CALLER, "> Filter Type < not found"))?;
        if !["simple", "pde", "distance"].contains(&ftype.as_str()) {
            return Err(SolverError::fatal(
                CALLER,
                format!("Unknown filter type: {ftype}"),
            ));
        }
        let phys = p.integer("primary solver index").unwrap_or(1) as usize;
        let slot = ctx
            .sim
            .slots
            .iter()
            .find(|s| s.id == phys)
            .ok_or_else(|| SolverError::fatal(CALLER, "no primary solver"))?;
        if slot.params.logical("local matrix storage") != Some(true) {
            return Err(SolverError::fatal(
                CALLER,
                "Primary solver should have active \"Local Matrix Storage\"",
            ));
        }
        let uvar = slot
            .variable
            .ok_or_else(|| SolverError::fatal(CALLER, "primary solver has no variable"))?;
        let dofs = ctx.sim.vars.by_id(uvar).dofs;
        let (x_phys, eperm) = elemental(ctx, "topo rho")?;
        let nsize = x_phys.len();
        let mut dv0 = vec![0.0; nsize];
        for i in 0..mesh.n_bulk() {
            let j = eperm.get(i).copied().unwrap_or(NO_DOF);
            if j == NO_DOF {
                continue;
            }
            dv0[j as usize] = element_size(ctx, ElemId(i as u32))?;
        }
        // The filters index the bulk elements; the fields their Perm.
        let identity = eperm.len() >= mesh.n_bulk()
            && (0..mesh.n_bulk()).all(|i| eperm[i] == i as u32)
            && nsize == mesh.n_bulk();
        let wmin = ctx
            .creal(&p, "sensitivity filter threshold")?
            .unwrap_or(1.0e-3f32 as f64);
        let (mut fmat, mut niter, mut pde_diff, mut interface_node) = (None, 1, None, None);
        if ftype == "pde" {
            let dim = mesh.dim();
            let d = match p.get("pde filter diffusion constant") {
                Some(Value::Real { rows, cols, data }) => {
                    let mut d = [0.0f64; 3];
                    if *rows == 1 && *cols == 1 {
                        d = [data[0]; 3];
                    } else if (*rows == 1 && *cols >= dim) || (*rows >= dim && *cols == 1) {
                        d[..dim].copy_from_slice(&data[..dim]);
                    } else {
                        return Err(SolverError::fatal(
                            CALLER,
                            format!("Invalid size for \"PDE Filter Diffusion Constant\": {rows} x {cols}"),
                        ));
                    }
                    d
                }
                _ => {
                    return Err(SolverError::fatal(
                        CALLER,
                        "> PDE Filter Diffusion Constant < not found",
                    ))
                }
            };
            pde_diff = Some(d);
            if flag(&p, "pde filter skip interface") {
                let mut count = vec![0usize; mesh.n_nodes()];
                for t in 0..mesh.n_bulk() {
                    for &k in mesh.element(ElemId(t as u32)).nodes {
                        count[k as usize] += 1;
                    }
                }
                interface_node = Some(count.iter().map(|&c| c < 4).collect());
            }
        } else {
            if ctx.comm.size() > 1 {
                return Err(SolverError::fatal(
                    CALLER,
                    "Only PDE Filter is implemented in parallel!",
                ));
            }
            if !identity {
                return Err(unsupported(
                    "element filters on a solver not active on every bulk element",
                ));
            }
            let mut f = if ftype == "simple" {
                niter = p.integer("simple filter iterations").unwrap_or(0).max(1) as usize;
                simple_filter(mesh)?
            } else {
                let r = ctx.creal(&p, "distance filter radius")?.ok_or_else(|| {
                    SolverError::fatal(CALLER, "> Distance Filter Radius < not found")
                })?;
                distance_filter(mesh, r)?
            };
            f.normalize();
            fmat = Some(f);
        }
        let mut st = State {
            method,
            fmat,
            niter,
            pde_diff,
            interface_node,
            phys,
            dofs,
            x_topo: x_phys,
            dv0: dv0.clone(),
            dv: Vec::new(),
            wmin,
        };
        st.dv = match method {
            Method::Density => {
                let mut dv = vec![0.0; nsize];
                self.filter(ctx, &st, &dv0, &mut dv, None, true)?;
                dv
            }
            _ => dv0,
        };
        Ok(st)
    }
}

impl Solver for TopoOpt {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        for k in [
            "solve adjoint problem",
            "use adjoint solution",
            "create zero levelset",
        ] {
            if flag(p, k) {
                return Err(unsupported(k));
            }
        }
        if p.integer("periodic physsolver").unwrap_or(0) > 1 {
            return Err(unsupported("Periodic PhysSolver"));
        }
        if p.integer("no modes").unwrap_or(0) > 0 {
            return Err(unsupported("constraint modes (No Modes)"));
        }
        if ctx
            .simulation_params()
            .integer("periodic timesteps")
            .unwrap_or(0)
            != 0
        {
            return Err(unsupported("Periodic Timesteps cycles"));
        }
        let cs = ctx
            .simulation_params()
            .string("coordinate system")
            .map(normalize_key)
            .unwrap_or_else(|| "cartesian".into());
        if !cs.starts_with("cartesian") {
            return Err(unsupported("axisymmetric weighting"));
        }
        if ctx.sim.backends.fe.is_none() {
            return Err(SolverError::NotWired(
                "basis and quadrature (elmer-fe)".into(),
            ));
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "-------------------------------------------", 4);
        ctx.info(CALLER, "Updating density for topology optimization", 4);
        let st = match self.state.take() {
            Some(s) => s,
            None => self.first_visit(ctx)?,
        };
        self.visit(ctx, st)
    }
}

impl TopoOpt {
    fn visit(&mut self, ctx: &mut SolverCtx<'_>, mut st: State) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        let need = |ctx: &SolverCtx<'_>, k: &str| -> Result<f64, SolverError> {
            ctx.creal(&p, k)?
                .ok_or_else(|| SolverError::fatal(CALLER, format!("> {k} < not found")))
        };
        let penal = need(ctx, "penalty exponent")?;
        let vol_frac = need(ctx, "volume fraction")?;
        let emin = need(ctx, "minimum relative density")?;
        let efrac = 1.0 - emin;
        let (mut x_phys, eperm) = elemental(ctx, "topo rho")?;
        let nsize = x_phys.len();
        // ObjectiveGradients (cycle mode 0: obj and dc start from zero).
        let (mut ce, _) = elemental(ctx, "topo ce")?;
        let mut dc = vec![0.0; nsize];
        let mut obj = 0.0;
        {
            let slot = ctx
                .sim
                .slots
                .iter()
                .find(|s| s.id == st.phys)
                .expect("checked on the first visit");
            let store = slot.local_store.as_ref().ok_or_else(|| {
                SolverError::fatal(CALLER, "Could not find local stiffness matrix!")
            })?;
            let u = ctx.sim.vars.by_id(slot.variable.expect("checked"));
            let uperm = u.perm.as_deref().expect("nodal");
            let dofs = st.dofs;
            let elems: Vec<ElemId> = ctx.active_elements().collect();
            for &e in &elems {
                let nodes = ctx.mesh.element(e).nodes;
                let m = dofs * nodes.len();
                let (k, _) = store
                    .get(e)
                    .filter(|(k, _)| k.len() == m * m)
                    .ok_or_else(|| {
                        SolverError::fatal(CALLER, "Could not find local stiffness matrix!")
                    })?;
                let mut sol = vec![0.0; m];
                for (i, &nd) in nodes.iter().enumerate() {
                    let s = uperm[nd as usize];
                    if s == NO_DOF {
                        continue;
                    }
                    for d in 0..dofs {
                        sol[i * dofs + d] = u.values[s as usize * dofs + d];
                    }
                }
                let mut sum = 0.0;
                for i in 0..m {
                    let mut act = 0.0;
                    for j in 0..m {
                        act += k[i * m + j] * sol[j];
                    }
                    sum += sol[i] * act;
                }
                ce[eperm[e.0 as usize] as usize] = sum;
            }
        }
        for i in 0..nsize {
            dc[i] -= penal * x_phys[i].powf(penal - 1.0) * efrac * ce[i];
        }
        for i in 0..nsize {
            obj += (emin + efrac * x_phys[i].powf(penal)) * ce[i];
        }
        ctx.info(CALLER, &format!(" Objective function: {obj:e}"), 4);
        // ApplyFilter( Fmat, dc, dc, ... ): the output starts as dc.
        let mut dcf = dc.clone();
        match st.method {
            Method::Sensitivity => {
                let xt = st.x_topo.clone();
                self.filter(ctx, &st, &dc, &mut dcf, Some((&xt, st.wmin)), false)?;
                dc = dcf;
            }
            Method::Density => {
                self.filter(ctx, &st, &dc, &mut dcf, None, true)?;
                dc = dcf;
            }
            Method::None => ctx.info(CALLER, "Applying no filtering", 4),
        }
        update_densities(ctx, &p, &mut st.x_topo, &dc, &st.dv, vol_frac)?;
        match st.method {
            Method::Density => {
                let xt = st.x_topo.clone();
                self.filter(ctx, &st, &xt, &mut x_phys, None, false)?;
            }
            _ => x_phys.copy_from_slice(&st.x_topo),
        }
        if let Some(v) = ctx.primary_mut() {
            if v.values.len() == 1 {
                v.values[0] = obj;
            }
        }
        if ctx.variable("topo bw").is_some() {
            let bw = topology_bw(ctx, &st.dv0, &x_phys, vol_frac);
            set_values(ctx, "topo bw", &bw);
        }
        let mult: Vec<f64> = x_phys
            .iter()
            .map(|x| emin + efrac * x.powf(penal))
            .collect();
        set_values(ctx, "topo rho", &x_phys);
        set_values(ctx, "topo ce", &ce);
        set_values(ctx, "topo dc", &dc);
        set_values(ctx, "topo mult", &mult);
        self.times_visited += 1;
        self.state = Some(st);
        Ok(())
    }
}

/// `UpdateDensities( x, dc, dv, g )`: the optimality-criteria update with
/// the Lagrange multiplier of the volume constraint found by bisection.
fn update_densities(
    ctx: &SolverCtx<'_>,
    p: &ParamMap,
    x: &mut [f64],
    dc: &[f64],
    dv: &[f64],
    vol_frac: f64,
) -> Result<(), SolverError> {
    let mv = ctx.creal(p, "bisection search max change")?.unwrap_or(0.2);
    let tol = ctx
        .creal(p, "bisection search tolerance")?
        .unwrap_or(1.0e-6);
    let damp = ctx
        .creal(p, "bisection search damping exponent")?
        .unwrap_or(0.5);
    let (mut l1, mut l2) = (0.0f64, 1.0e9f64);
    let mut sdv = 0.0;
    for v in dv {
        sdv += v;
    }
    let v0 = vol_frac * sdv;
    let floor = 1.0e-10f32 as f64;
    let mut xnew = vec![0.0; x.len()];
    let mut vi = 0.0;
    let mut k = 0;
    for it in 1..=1000 {
        k = it;
        let lmid = 0.5 * (l2 + l1);
        for i in 0..x.len() {
            let r = -dc[i] / (dv[i] * lmid);
            let t = x[i] * floor.max(r).powf(damp);
            xnew[i] = 0.0f64.max((x[i] - mv).max(1.0f64.min((x[i] + mv).min(t))));
        }
        vi = 0.0;
        for i in 0..x.len() {
            vi += dv[i] * xnew[i];
        }
        if vi > v0 {
            l1 = lmid;
        } else {
            l2 = lmid;
        }
        let err = (l2 - l1) / (l1 + l2);
        if err < tol {
            break;
        }
    }
    x.copy_from_slice(&xnew);
    ctx.info(CALLER, &format!("Number of bisection iterations: {k}"), 7);
    ctx.info(
        CALLER,
        &format!("Volume target and accuracy: {v0:e} {:e}", vi - v0),
        7,
    );
    Ok(())
}

/// `DefineTopologyBW`: 1 where the density exceeds the volume-conserving
/// threshold, found by ten refinements of a 100-bin histogram.
fn topology_bw(ctx: &SolverCtx<'_>, dv0: &[f64], x: &[f64], vol_frac: f64) -> Vec<f64> {
    let m = 100usize;
    let (mut xlow, mut xup) = (0.0f64, 1.0f64);
    let mut vtot = 0.0;
    for v in dv0 {
        vtot += v;
    }
    let vtarget = vol_frac * vtot;
    let mut q = 1.0f64;
    let mut iters = 0;
    for iter in 1..=10 {
        iters = iter;
        let h = (xup - xlow) / m as f64;
        let mut histv = vec![0.0f64; m + 1];
        for i in 0..x.len() {
            let j = ((x[i] - xlow) / h).ceil();
            let j = (j.min(m as f64)).max(0.0) as usize;
            histv[j] += dv0[i];
        }
        let mut cumv = vec![0.0f64; m + 1];
        cumv[0] = histv[0];
        for i in 1..=m {
            cumv[i] = cumv[i - 1] + histv[i];
        }
        let mut hit = false;
        q = 1.0;
        for i in 1..=m {
            if cumv[i - 1] < vtarget && cumv[i] > vtarget {
                xlow += (i - 1) as f64 * h;
                xup = xlow + h;
                q = (vtarget - cumv[i - 1]) / (cumv[i] - cumv[i - 1]);
                hit = true;
                break;
            } else if (cumv[i - 1] - vtarget).abs() < f64::EPSILON {
                xlow += (i - 1) as f64 * h;
                xup = xlow;
                break;
            } else if (cumv[i] - vtarget).abs() < f64::EPSILON {
                xlow += i as f64 * h;
                xup = xlow;
                break;
            }
        }
        if !hit {
            break;
        }
    }
    let xmid = (1.0 - q) * xlow + q * xup;
    ctx.info(
        CALLER,
        &format!("Mass conserving B&W limit after {iters} iters: {xmid:e}"),
        7,
    );
    x.iter()
        .map(|&v| if v > xmid { 1.0 } else { 0.0 })
        .collect()
}
