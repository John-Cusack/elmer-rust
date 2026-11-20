//! `SteadyPhaseChange` and its user function `MeltingHeat`
//! (fem/src/modules/SteadyPhaseChange.F90): the steady position of a
//! solid-liquid interface in 2D, moved each visit onto the melting-point
//! isotherm of the temperature field.
//!
//! The interface is the boundary elements of the BCs with `Phase Change =
//! True`; the solver's variable (`PhaseSurface`) is the height of the
//! interface in its normal direction (the coordinate in which the interface
//! nodes span the least). On each visit (`isoterm formulation`):
//!
//! - `CreateIsotherm`: the segments where `T - MeltPoint` changes sign in
//!   the bulk elements of the solid and the liquid body;
//! - each interface node moves by its normal distance to the nearest
//!   isotherm segment (`SurfaceMove`, saved as `<var>Diff`), relaxed by
//!   `Nonlinear System Relaxation Factor`;
//! - the norm is `sqrt(sum s^2 / NoBNodes)` over the interface nodes.
//!
//! The Newton formulation (after `Nonlinear System Newton After
//! Iterations` visits, or once the change is below `... Newton After
//! Tolerance`) is ported as written, including its use of the boundary
//! element counter as a node index (`Temperature(TempPerm(t))`). The mesh
//! itself is moved by the case's mesh update solver (`Mesh Update 2 = Equals
//! PhaseSurface`); `NewY` holds the mesh's current y coordinates.
//!
//! `MeltingHeat` (`Heat Flux = Variable Coordinate 1; Real Procedure
//! "SteadyPhaseChange" "MeltingHeat"`) is the latent heat released by the
//! pull velocity through the interface: `L rho_solid (n . v_pull)`, the
//! normal of the current boundary element at its centre (`NormalVector`
//! with `Check`, `Normal Target Body` honoured), from the moving mesh.
//!
//! Refused: `Passive Steps`, `Lumped Acceleration After Iterations`,
//! `Internal Mesh Movement`, `Save Isotherm`, a `Normals` field, and other
//! than 2D meshes.

use std::f64::consts::PI;
use std::sync::Arc;

use elmer_core::{BasisEval, ElemId, ParamMap, Params, Value, NO_DOF};
use elmer_solver::elmer_expr::{self, Udf};
use elmer_solver::variable::Variable;
use elmer_solver::{Model, Solver, SolverCtx, SolverError};

const CALLER: &str = "SteadyPhaseChange";
const AEPS: f64 = 10.0 * f64::EPSILON;

pub fn factory(p: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    for key in ["passive steps", "lumped acceleration after iterations"] {
        if p.get(key).is_some() {
            return Err(SolverError::Unsupported(format!("SteadyPhaseChange: {key}")));
        }
    }
    for key in ["internal mesh movement", "save isotherm"] {
        if p.logical(key) == Some(true) {
            return Err(SolverError::Unsupported(format!("SteadyPhaseChange: {key}")));
        }
    }
    Ok(Box::new(SteadyPhaseChange::default()))
}

/// The SAVEd state of the Fortran routine.
#[derive(Default)]
pub struct SteadyPhaseChange {
    visited: usize,
    ready: bool,
    newton: bool,
    /// Interface nodes (`BoundaryMarker`) and boundary elements.
    marker: Vec<bool>,
    phase_elements: Vec<(usize, ElemId)>,
    nob_nodes: usize,
    normal_dir: usize,
    tangent_dir: usize,
    width: f64,
    trip_node: usize,
    solid_body: usize,
    liquid_body: usize,
    prev_temp: Vec<f64>,
}

fn info_i(ctx: &SolverCtx<'_>, what: &str, v: usize) {
    ctx.info(CALLER, &format!("{what:<34}{v:>12}"), 4);
}

fn info_r(ctx: &SolverCtx<'_>, what: &str, v: f64) {
    ctx.info(CALLER, &format!("{what:<34}{}", elmer_solver::elmer_saveio::fortran::fortran_es(v, 12, 4)), 4);
}

impl Solver for SteadyPhaseChange {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        register_melting_heat(ctx);
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "--------------------------------------------", 4);
        ctx.info(CALLER, "Using steady algorithm to find the isotherm", 4);
        ctx.info(CALLER, "--------------------------------------------", 4);
        self.visited += 1;
        if ctx.mesh.dim() != 2 {
            return Err(SolverError::fatal(CALLER, "Implemented only in 2D"));
        }
        let p = ctx.params().clone();
        let (surf_name, surf_perm) = {
            let v = ctx.primary().ok_or_else(|| SolverError::Unsupported("SteadyPhaseChange: no Variable (Surface Variable)".into()))?;
            let perm = v.perm.map(|q| q.to_vec()).ok_or_else(|| SolverError::fatal(CALLER, "Surface field needed for Phase Change"))?;
            (v.name.to_string(), perm)
        };
        if surf_perm.iter().all(|&s| s == NO_DOF) {
            return Err(SolverError::fatal(CALLER, "Surface field needed for Phase Change"));
        }
        let nslots = ctx.primary().expect("checked").values.len();
        // `<var>Diff` on the surface's permutation.
        let diff_name = format!("{surf_name}diff");
        if ctx.variable(&diff_name).is_none() {
            let perm = ctx.sim.vars.get(&surf_name).and_then(|v| v.perm.clone());
            ctx.sim.vars.add(Variable::new(&diff_name, elmer_core::FieldKind::Nodal, 1, perm, vec![0.0; nslots]))?;
        }
        let temp_name = p.string("phase change variable").map_or_else(|| "temperature".to_string(), elmer_core::normalize_key);
        let (temp, temp_perm) = {
            let v = ctx.variable(&temp_name).ok_or_else(|| SolverError::fatal(CALLER, "Temperature field needed for Phase Change"))?;
            (v.values.to_vec(), v.perm.map(|q| q.to_vec()))
        };
        let tslot = |k: usize| -> Option<usize> {
            match &temp_perm {
                None => Some(k),
                Some(q) => (q[k] != NO_DOF).then(|| q[k] as usize),
            }
        };
        if temp_perm.as_ref().is_some_and(|q| q.iter().all(|&s| s == NO_DOF)) {
            return Err(SolverError::fatal(CALLER, "Temperature field has zero Perm vector!"));
        }
        let relax = p.real("nonlinear system relaxation factor").unwrap_or(1.0);
        let mesh = ctx.mesh;
        let model: &Model = &ctx.sim.model;
        let nb = mesh.n_bulk();
        let phase_bc = |e: ElemId| -> bool {
            let el = mesh.element(e);
            el.etype.family() != 1
                && model.bc_of_tag(el.tag).is_some_and(|b| model.boundary_conditions[b - 1].logical("phase change") == Some(true))
        };

        if !self.ready {
            let n_nodes = mesh.n_nodes();
            self.marker = vec![false; n_nodes];
            for t in 0..mesh.n_boundary() {
                let e = ElemId((nb + t) as u32);
                if !phase_bc(e) {
                    continue;
                }
                for &k in mesh.element(e).nodes {
                    self.marker[k as usize] = true;
                }
                self.phase_elements.push((t + 1, e));
            }
            self.nob_nodes = self.marker.iter().filter(|&&m| m).count();
            info_i(ctx, "Number of interface nodes:", self.nob_nodes);
            if self.nob_nodes == 0 {
                ctx.warn(CALLER, "There is no Phase Change boundary");
            }
            info_i(ctx, "Number of interface elements:", self.phase_elements.len());
            let dim = 2;
            let (mut cmax, mut cmin) = ([-f64::MAX; 3], [f64::MAX; 3]);
            let (mut imax, mut imin) = ([0usize; 3], [0usize; 3]);
            for k in 0..n_nodes {
                if !self.marker[k] {
                    continue;
                }
                let c = mesh.coord(k as u32);
                for j in 0..dim {
                    if c[j] > cmax[j] {
                        cmax[j] = c[j];
                        imax[j] = k;
                    }
                    if c[j] < cmin[j] {
                        cmin[j] = c[j];
                        imin[j] = k;
                    }
                }
            }
            let mut j = 0;
            for i in 0..dim {
                if cmax[i] - cmin[i] < cmax[j] - cmin[j] {
                    j = i;
                }
            }
            self.normal_dir = j;
            let mut j = 0;
            for i in 0..dim {
                if cmax[i] - cmin[i] > cmax[j] - cmin[j] {
                    j = i;
                }
            }
            self.width = cmax[j] - cmin[j];
            self.tangent_dir = j;
            info_r(ctx, "Width of the interface", self.width);
            self.trip_node = imax[self.tangent_dir];
            info_i(ctx, "Index of the triple point: ", self.trip_node + 1);
            info_i(ctx, "Index of the axis point: ", imin[self.tangent_dir] + 1);
            // `GotLiquid` and `GotSolid` start false and are ANDed, so the
            // bodies always come from the first interface element's parents:
            // the solid is the one reaching higher.
            if let Some(&(_, e)) = self.phase_elements.first() {
                let [l, r] = mesh.element(e).parents;
                let (Some(l), Some(r)) = (l, r) else {
                    return Err(SolverError::fatal(CALLER, "a phase change boundary without two parents"));
                };
                let ymax = |pe: ElemId| mesh.element(pe).nodes.iter().map(|&k| mesh.coord(k)[1]).fold(f64::MIN, f64::max);
                let body = |pe: ElemId| model.body_of_tag(mesh.element(pe).tag).unwrap_or(0);
                if ymax(l) > ymax(r) {
                    (self.solid_body, self.liquid_body) = (body(l), body(r));
                } else {
                    (self.solid_body, self.liquid_body) = (body(r), body(l));
                }
            }
            info_i(ctx, "Body index for solid: ", self.solid_body);
            info_i(ctx, "Body index for liquid: ", self.liquid_body);
            self.prev_temp = vec![0.0; nslots];
            self.ready = true;
        }

        let trip_temp = tslot(self.trip_node).map_or(0.0, |s| temp[s]);
        if let Some(n) = p.integer("nonlinear system newton after iterations") {
            if self.visited > n.max(0) as usize {
                self.newton = true;
            }
        }
        let newton_after_tol = p.real("nonlinear system newton after tolerance").unwrap_or(0.0);
        ctx.info(CALLER, if self.newton { "Steady state newton formulation" } else { "Steady state isoterm formulation" }, 4);
        let triple_fixed = p.logical("triple point fixed") == Some(true);
        let melt = if triple_fixed {
            info_r(ctx, "Melting point set: ", trip_temp);
            trip_temp
        } else {
            match model.materials.iter().find_map(|m| m.real("melting point")) {
                Some(m) => {
                    info_r(ctx, "Melting point found: ", m);
                    m
                }
                None => {
                    ctx.info(CALLER, "Could not find melting point in any material!", 4);
                    0.0
                }
            }
        };
        let (iso, xmin, xmax) = if self.newton {
            (Vec::new(), 0.0, 0.0)
        } else {
            let (iso, xmin, xmax) = self.isotherm(ctx, &temp, temp_perm.as_deref(), melt)?;
            info_i(ctx, "Number of isotherm segments:", iso.len());
            (iso, xmin, xmax)
        };

        // The surface moves, node by node over the interface elements.
        let mut moves: Vec<f64> = ctx.variable(&diff_name).expect("added").values.to_vec();
        let mut done = vec![false; nslots];
        let (mut area, mut volume, mut tave, mut tabs, mut volabs, mut max_tdiff) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0f64);
        let cartesian = model
            .simulation
            .string("coordinate system")
            .is_none_or(|c| elmer_core::normalize_key(c).starts_with("cartesian"));
        let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no FE backend"))?;
        let mut be = BasisEval::default();
        let (td, nd) = (self.tangent_dir, self.normal_dir);
        for &(t, e) in &self.phase_elements {
            let el = mesh.element(e);
            if !(200..=203).contains(&el.etype.0) {
                return Err(SolverError::fatal("PhaseChangeSolve", "Implemented only for elements 202 and 203!"));
            }
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
            let tdiff: Vec<f64> = el
                .nodes
                .iter()
                .map(|&k| tslot(k as usize).map_or(f64::NAN, |s| temp[s]) - melt)
                .collect();
            max_tdiff = max_tdiff.max(tdiff.iter().fold(0.0f64, |m, x| m.max(x.abs())));
            for nn in 0..n {
                let node = el.nodes[nn] as usize;
                let k = surf_perm[node];
                if k == NO_DOF || done[k as usize] {
                    continue;
                }
                let k = k as usize;
                done[k] = true;
                if triple_fixed && node == self.trip_node {
                    moves[k] = 0.0;
                    ctx.info(CALLER, "Triple point position fixed", 4);
                    continue;
                }
                if nn == 2 {
                    let a = surf_perm[el.nodes[0] as usize] as usize;
                    let b = surf_perm[el.nodes[1] as usize] as usize;
                    moves[k] = 0.5 * (moves[a] + moves[b]);
                }
                // `Temperature(TempPerm(t))`: the element counter as a node.
                let ttemp = tslot(t - 1).map_or(0.0, |s| temp[s]);
                let update = if self.newton {
                    let mut dtdz = ttemp - self.prev_temp[k];
                    if dtdz.abs() < AEPS {
                        ctx.warn(CALLER, "Very small temperature update.");
                        dtdz = 1.0;
                    }
                    moves[k] * (melt - ttemp) / dtdz
                } else {
                    let eps = 1.0e-6 * (xmax - xmin);
                    let (xx, yy) = if td == 0 { (xyz[nn][0], xyz[nn][1]) } else { (xyz[nn][1], xyz[nn][0]) };
                    let (mut dxmin, mut dymin) = (f64::MAX, f64::MAX);
                    let mut stat = false;
                    let mut imin = 0usize;
                    let mut i = 0;
                    while i + 1 < iso.len() {
                        let (x1, x2) = (iso[i][td], iso[i + 1][td]);
                        let (y1, y2) = (iso[i][nd], iso[i + 1][nd]);
                        if xx > x1 - eps && xx < x2 + eps {
                            dxmin = 0.0;
                            let mut d = (yy - y1).abs().min((yy - y2).abs());
                            d += (x1 - xx).max(0.0);
                            d += (xx - x2).max(0.0);
                            if d <= dymin {
                                stat = true;
                                dymin = d;
                                imin = i;
                            }
                        }
                        if !stat {
                            let d = (xx - x1).abs().min((xx - x2).abs());
                            if d <= dxmin {
                                dxmin = d;
                                imin = i;
                            }
                        }
                        i += 2;
                    }
                    let (x1, x2) = (iso[imin][td], iso[imin + 1][td]);
                    let (y1, y2) = (iso[imin][nd], iso[imin + 1][nd]);
                    if !stat && dxmin > 1.0e-2 * (x1 - x2).abs() {
                        ctx.warn(CALLER, "Isotherm error?");
                    }
                    if (x2 - x1).abs() > AEPS {
                        (y1 - yy) + (xx - x1) * (y2 - y1) / (x2 - x1)
                    } else {
                        0.5 * (y1 + y2) - yy
                    }
                };
                self.prev_temp[k] = ttemp;
                moves[k] = update;
            }
            for qp in fe.quadrature.gauss_points(el.etype, None) {
                fe.basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let b = &be.basis[..n];
                let mut s = qp.weight * be.det_j;
                if !cartesian {
                    s = s * (0..n).fold(0.0, |acc, j| acc + b[j] * xyz[j][0]) * 2.0 * PI;
                }
                let sm = |j: usize| moves[surf_perm[el.nodes[j] as usize] as usize];
                area += s;
                volume += s * (0..n).fold(0.0, |acc, j| acc + b[j] * sm(j));
                tave += s * (0..n).fold(0.0, |acc, j| acc + b[j] * tdiff[j]);
                volabs += s * (0..n).fold(0.0, |acc, j| acc + b[j] * sm(j).abs());
                tabs += s * (0..n).fold(0.0, |acc, j| acc + b[j] * tdiff[j].abs());
            }
        }
        let _ = (tave / area, tabs / area, volume / area, volabs / area);

        for m in moves.iter_mut() {
            *m *= relax;
        }
        let dpos;
        let (norm, change, max_move, max_surf) = {
            let v = ctx.primary_mut().expect("checked");
            for (s, m) in v.values.iter_mut().zip(&moves) {
                *s += m;
            }
            dpos = moves[surf_perm[self.trip_node] as usize];
            let (mut norm, mut change) = (0.0, 0.0);
            for k in 0..self.marker.len() {
                if !self.marker[k] {
                    continue;
                }
                let s = surf_perm[k] as usize;
                norm += v.values[s] * v.values[s];
                change += moves[s] * moves[s];
            }
            let change = (change / norm).sqrt();
            let norm = (norm / self.nob_nodes as f64).sqrt();
            let max_move = moves.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            let max_surf = v.values.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            (norm, change, max_move, max_surf)
        };
        ctx.variable_mut(&diff_name).expect("added").values.copy_from_slice(&moves);
        if change.abs() < newton_after_tol {
            self.newton = true;
        }
        info_r(ctx, "Result Norm:", norm);
        info_r(ctx, "Relative Change: ", change);
        info_r(ctx, "Max Height:", max_surf);
        info_r(ctx, "Max Change in Height:", max_move);
        info_r(ctx, "Max Angle:       ", max_surf / self.width);
        info_r(ctx, "Maximum temperature difference: ", max_tdiff);
        if let Some(id) = ctx.sim.vars.id(&surf_name) {
            ctx.sim.vars.by_id_mut(id).conv.norm = norm;
        }

        let sim = &mut ctx.sim.model.simulation;
        sim.insert("res: Triple point temperature", Value::real(trip_temp));
        sim.insert("res: triple point movement", Value::real(dpos));
        // `NewY`: the mesh's y coordinates (a pointer to them in the Fortran).
        let y: Vec<f64> = (0..ctx.mesh.n_nodes() as u32).map(|k| ctx.mesh.coord(k)[1]).collect();
        match ctx.variable_mut("newy") {
            Some(v) => v.values.copy_from_slice(&y),
            None => {
                ctx.sim.vars.add(Variable::new("NewY", elmer_core::FieldKind::Nodal, 1, None, y))?;
            }
        }
        Ok(())
    }
}

impl SteadyPhaseChange {
    /// `CreateIsotherm`: the end points of the melting-point isotherm in
    /// the bulk elements of the solid and liquid bodies, two per element,
    /// ordered along the tangent direction; and the x range of the points.
    #[allow(clippy::type_complexity)]
    fn isotherm(
        &self,
        ctx: &SolverCtx<'_>,
        temp: &[f64],
        temp_perm: Option<&[u32]>,
        melt: f64,
    ) -> Result<(Vec<[f64; 2]>, f64, f64), SolverError> {
        let mesh = ctx.mesh;
        let model = &ctx.sim.model;
        let (mut xmin, mut xmax) = (f64::MAX, -f64::MAX);
        let mut iso: Vec<[f64; 2]> = Vec::new();
        for t in 0..mesh.n_bulk() {
            let el = mesh.element(ElemId(t as u32));
            if el.etype.0 < 300 {
                continue;
            }
            let body = model.body_of_tag(el.tag).unwrap_or(0);
            if body != self.liquid_body && body != self.solid_body {
                continue;
            }
            let slot = |k: u32| temp_perm.map_or(Some(k as usize), |q| (q[k as usize] != NO_DOF).then(|| q[k as usize] as usize));
            let Some(slots) = el.nodes.iter().map(|&k| slot(k)).collect::<Option<Vec<usize>>>() else { continue };
            let tdiff: Vec<f64> = slots.iter().map(|&s| temp[s] - melt).collect();
            if tdiff.iter().all(|&x| x < 0.0) || tdiff.iter().all(|&x| x > 0.0) {
                continue;
            }
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
            let vertex = (el.etype.0 / 100) as usize;
            let mut n = 0;
            for nn in 0..vertex {
                let next = (nn + 1) % vertex;
                let (t1, t2) = (tdiff[nn], tdiff[next]);
                if (t1 < 0.0 && 0.0 <= t2) || (t2 <= 0.0 && 0.0 < t1) {
                    n += 1;
                    if n <= 2 {
                        let px = xyz[nn][0] + t1 * ((xyz[next][0] - xyz[nn][0]) / (t1 - t2));
                        let py = xyz[nn][1] + t1 * ((xyz[next][1] - xyz[nn][1]) / (t1 - t2));
                        xmin = px.min(xmin);
                        xmax = px.max(xmax);
                        iso.push([px, py]);
                    } else {
                        ctx.warn(CALLER, "Wiggly Isotherm");
                    }
                }
            }
            if n == 1 {
                iso.pop();
                continue;
            }
            if n >= 2 {
                let l = iso.len();
                if iso[l - 2][self.tangent_dir] > iso[l - 1][self.tangent_dir] {
                    iso.swap(l - 2, l - 1);
                }
            }
        }
        if iso.is_empty() {
            return Err(SolverError::fatal(CALLER, "Isotherm is empty thus cannot map phase change surface"));
        }
        Ok((iso, xmin, xmax))
    }
}

/// The facts `MeltingHeat` reads about a boundary element: its body (the
/// BC's `Body Id`), its parents' bodies and corner nodes, and the BC's
/// `Normal Target Body`.
struct BoundaryFacts {
    body: Option<usize>,
    parents: [Option<(usize, Vec<u32>)>; 2],
    out_body: Option<usize>,
}

/// Registers `("SteadyPhaseChange", "MeltingHeat")` for the run.
fn register_melting_heat(ctx: &SolverCtx<'_>) {
    let mesh = ctx.mesh;
    let model: Arc<Model> = Arc::new(ctx.sim.model.clone());
    let nb = mesh.n_bulk();
    let mut facts: Vec<Option<BoundaryFacts>> = Vec::new();
    for t in 0..mesh.n_boundary() {
        let el = mesh.element(ElemId((nb + t) as u32));
        let bc = model.bc_of_tag(el.tag);
        let parent = |p: Option<ElemId>| {
            p.and_then(|p| {
                let pe = mesh.element(p);
                Some((model.body_of_tag(pe.tag)?, pe.nodes.to_vec()))
            })
        };
        facts.push(Some(BoundaryFacts {
            body: model.boundary_body(el.tag),
            parents: [parent(el.parents[0]), parent(el.parents[1])],
            out_body: bc.and_then(|b| model.boundary_conditions[b - 1].integer("normal target body")).filter(|&b| b > 0).map(|b| b as usize),
        }));
    }
    // The parents' centres from the mesh as it is at setup: only the sign
    // of the normal depends on them.
    let coords: Vec<[f64; 3]> = (0..mesh.n_nodes() as u32).map(|k| mesh.coord(k)).collect();
    let facts = Arc::new(facts);
    let coords = Arc::new(coords);
    let udf = move |node: u32, _deps: &[f64]| -> f64 {
        melting_heat(&model, &facts, &coords, nb, node).unwrap_or(f64::NAN)
    };
    ctx.sim.backends.udfs.register_fortran("SteadyPhaseChange", "MeltingHeat", Udf::Node(Arc::new(udf)));
}

fn material(model: &Model, body: usize) -> Option<&ParamMap> {
    let m = model.bodies.get(body.checked_sub(1)?)?.integer("material")?;
    model.materials.get((m as usize).checked_sub(1)?)
}

/// `MeltingHeat(Model, Node, t)` for the current boundary element.
fn melting_heat(model: &Model, facts: &[Option<BoundaryFacts>], coords: &[[f64; 3]], nb: usize, node: u32) -> Option<f64> {
    let cur = elmer_expr::current_element()?;
    let f = facts.get((cur.index as usize).checked_sub(nb)?)?.as_ref()?;
    if !cur.nodes.contains(&node) {
        return None;
    }
    let mat = material(model, f.body?)?;
    let latent = mat.real("latent heat")?;
    let sim = &model.simulation;
    let pull = [
        sim.real("res: pull velocity 1").or_else(|| mat.real("convection velocity 1")).unwrap_or(0.0),
        sim.real("res: pull velocity 2").or_else(|| mat.real("convection velocity 2")).unwrap_or(0.0),
        0.0,
    ];
    // The solid side's density: the left parent if solid, else the right.
    let left = f.parents[0].as_ref()?;
    let density = match material(model, left.0) {
        Some(m) if m.logical("solid") == Some(true) => m.real("density")?,
        _ => material(model, f.parents[1].as_ref()?.0)?.real("density")?,
    };
    // NormalVector(Element, Nodes, 0, 0, .TRUE.) of the line at its centre.
    let x = elmer_expr::current_coords();
    if x.len() != 2 {
        return None;
    }
    let d = [-0.5, 0.5];
    let (dxdu, dydu) = (d[0] * x[0][0] + d[1] * x[1][0], d[0] * x[0][1] + d[1] * x[1][1]);
    let det = 1.0 / (dxdu * dxdu + dydu * dydu).sqrt();
    let mut normal = [-dydu * det, dxdu * det, 0.0];
    // CheckNormalDirection: the parent the normal must point away from.
    let pick = match (&f.parents[0], &f.parents[1]) {
        (Some(l), Some(r)) => match f.out_body {
            Some(k) if l.0 == k => Some(r),
            Some(_) => Some(l),
            None if l.0 > r.0 => Some(l),
            None if l.0 < r.0 => Some(r),
            None => None,
        },
        (Some(l), None) => Some(l),
        (None, r) => r.as_ref(),
    };
    if let Some((_, pn)) = pick {
        let corners = pn.len().min(4);
        let c: [f64; 3] = std::array::from_fn(|a| pn[..corners].iter().map(|&k| coords[k as usize][a]).sum::<f64>() / corners as f64);
        let point: [f64; 3] = std::array::from_fn(|a| 0.5 * (x[0][a] + x[1][a]));
        let dot: f64 = (0..3).map(|a| (c[a] - point[a]) * normal[a]).sum();
        if dot > 0.0 {
            normal = normal.map(|c| -c);
        }
    }
    let normal_pull = normal.iter().zip(&pull).fold(0.0, |s, (n, v)| s + n * v);
    Some(latent * density * normal_pull)
}
