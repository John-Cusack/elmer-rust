//! Adaptive mesh refinement (Adaptive.F90): `RefineMesh` with the a
//! posteriori residual estimate of `ComputeError` and the internal
//! splitting of `SplitMesh` (`Adaptive Remesh = False`) or the external
//! mesher of `External_ReMesh` (`Adaptive Remesh = True`, [`remesh`]),
//! moving the run to the new mesh through `elmer_solver::remesh`.
//!
//! A solver supplies its residuals through [`Estimator`] (HeatSolve's
//! `HeatSolver_Inside_Residual`, `_Edge_Residual`, `_Boundary_Residual`).
//!
//! Refused: remeshing with MMG or Gmsh, `Flux
//! Recovery`, coarsening, `Adaptive Pre Smoothing`, outlier fractions and
//! non-triangular meshes (SplitMesh's own Fatal).
#![forbid(unsafe_code)]
#![allow(clippy::needless_range_loop)]

pub mod diffusion;
pub mod elasticity;
pub mod mesh;
pub mod remesh;
mod split;

use std::sync::Arc;

use elmer_core::{ElemId, FieldKind, MeshView, ParamMap, Params};
use elmer_solver::remesh::MeshChange;
use elmer_solver::{SolverCtx, SolverError, Variable};

pub use mesh::{Edge, Edges, OwnedMesh};

/// The keyword elmer-run gives refining solvers: the case directory, where
/// `Adaptive Save Mesh` writes `<mesh name>/mesh.*`.
pub const OUTPUT_DIRECTORY_KEY: &str = "Adaptive Output Directory";

/// `AEPS = 10 * EPSILON(1.0_dp)`.
const AEPS: f64 = 10.0 * f64::EPSILON;
const CALLER: &str = "RefineMesh";

#[derive(Debug)]
pub enum AdaptiveError {
    Fatal(String),
}

impl From<AdaptiveError> for SolverError {
    fn from(e: AdaptiveError) -> Self {
        match e {
            AdaptiveError::Fatal(m) => SolverError::fatal(CALLER, m),
        }
    }
}

/// The residual indicators of one solver, each `(norm, Indicator(1:2))`.
pub trait Estimator {
    /// `InsideResidual` of bulk element `e`: `(Fnorm, Indicator)`.
    fn inside(&mut self, ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, [f64; 2]), SolverError>;
    /// `EdgeResidual` of an interior edge.
    fn edge(&mut self, ctx: &SolverCtx<'_>, edge: &Edge) -> Result<[f64; 2], SolverError>;
    /// `BoundaryResidual` of boundary element `e`: `(Gnorm, Indicator)`.
    fn boundary(&mut self, ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, [f64; 2]), SolverError>;
}

/// `Element % hK` as `MeshStabParams` sets it: `StabParam` when a solver of
/// the mesh stabilises, else `ElementDiameter`.
pub fn element_hk(ctx: &SolverCtx<'_>, e: ElemId) -> Result<f64, SolverError> {
    let stabilize = ctx.sim.model.solvers.iter().any(|p| {
        p.logical("stabilize") == Some(true)
            || p.string("stabilization method").is_some_and(|s| {
                let s = elmer_core::normalize_key(s);
                s == "vms" || s == "stabilized"
            })
            || p.logical("pressure stabilization") == Some(true)
    });
    let long = ctx
        .simulation_params()
        .logical("stabilization use longest element edge")
        .unwrap_or(false);
    let el = ctx.mesh.element(e);
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    if stabilize {
        let (_, hk) = elmer_fe::stab_param(el.etype, &xyz, ctx.mesh.dim(), long)
            .map_err(|m| SolverError::fatal("MeshStabParams", m.to_string()))?;
        Ok(hk)
    } else {
        Ok(elmer_fe::element_diameter(el.etype, &xyz, long))
    }
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("RefineMesh: {what}"))
}

/// The keywords of RefineMesh this port does not follow.
pub fn check_params(p: &ParamMap) -> Result<(), SolverError> {
    for k in [
        "flux recovery",
        "adaptive error histogram",
    ] {
        if p.logical(k) == Some(true) {
            return Err(unsupported(k));
        }
    }
    for k in [
        "adaptive pre smoothing",
        "adaptive max outlier fraction",
        "adaptive use nodal error as limit",
    ] {
        if p.get(k).is_some() {
            return Err(unsupported(k));
        }
    }
    // `Adaptive Min H`, `Max H` and `Coarsening` only steer the remeshers.
    Ok(())
}

/// The SAVEd state of RefineMesh.
#[derive(Default)]
pub struct RefineState {
    do_final_ref: bool,
}

/// A nodal field of the current mesh with the identity permutation, created
/// with `init` when missing (`VariableAdd` without a Perm). Returns whether
/// it was created.
fn field(
    ctx: &mut SolverCtx<'_>,
    name: &str,
    output: bool,
    init: f64,
) -> Result<bool, SolverError> {
    if ctx.sim.vars.get(name).is_some() {
        return Ok(false);
    }
    let n = ctx.mesh.n_nodes();
    let mut v = Variable::new(name, FieldKind::Nodal, 1, None, vec![init; n]);
    v.output = output;
    ctx.sim.vars.add(v)?;
    Ok(true)
}

fn values(ctx: &SolverCtx<'_>, name: &str) -> Vec<f64> {
    ctx.sim.vars.get(name).expect("created").values.clone()
}

fn store(ctx: &mut SolverCtx<'_>, name: &str, v: Vec<f64>) {
    ctx.sim.vars.get_mut(name).expect("created").values = v;
}

/// `ComputeError`: the element error indicators and their maximum.
fn compute_error(
    ctx: &SolverCtx<'_>,
    est: &mut dyn Estimator,
) -> Result<(Vec<f64>, f64), SolverError> {
    let mesh = ctx.mesh;
    let nb = mesh.n_bulk();
    let edges = Edges::find(mesh);
    let mut fnorm = 0.0;
    let mut t = vec![[0.0f64; 2]; nb];
    for i in 0..nb {
        let (f, ind) = est.inside(ctx, ElemId(i as u32))?;
        fnorm += f;
        t[i][0] += ind[0];
        t[i][1] += ind[1];
    }
    if mesh.dim() != 2 {
        return Err(unsupported("three-dimensional meshes (face residuals)"));
    }
    for e in &edges.edges {
        let (Some(l), Some(r)) = (e.left, e.right) else {
            continue;
        };
        let ind = est.edge(ctx, e)?;
        t[l][0] += ind[0];
        t[l][1] += ind[1];
        t[r][0] += ind[0];
        t[r][1] += ind[1];
    }
    for b in nb..nb + mesh.n_boundary() {
        let el = mesh.element(ElemId(b as u32));
        if el.etype.0 == 101 {
            continue;
        }
        let parents = el.parents;
        let (g, ind) = est.boundary(ctx, ElemId(b as u32))?;
        fnorm += g;
        for p in parents.into_iter().flatten() {
            let p = p.0 as usize;
            t[p][0] += ind[0];
            t[p][1] += ind[1];
        }
    }
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    for x in &t {
        s1 += x[1];
    }
    for x in &t {
        s2 += x[0];
    }
    let s = s1.sqrt() / s2.sqrt();
    let mut ei: Vec<f64> = t
        .iter()
        .map(|x| (x[0] / (2.0 * s) + s * x[1] / 2.0).sqrt())
        .collect();
    if fnorm > AEPS {
        let sf = fnorm.sqrt();
        for x in &mut ei {
            *x /= sf;
        }
    }
    let max = ei.iter().copied().fold(f64::MIN, f64::max);
    Ok((ei, max))
}

/// `RefineMesh` with the internal splitting: estimates the error of the
/// solver's solution, decides whether the mesh has converged and, if not,
/// asks the run to go on with the split mesh.
pub fn refine_mesh(
    ctx: &mut SolverCtx<'_>,
    state: &mut RefineState,
    est: &mut dyn Estimator,
) -> Result<(), SolverError> {
    ctx.info(CALLER, " ", 5);
    ctx.info(
        CALLER,
        "----------- M E S H   R E F I N E M E N T --------------",
        5,
    );
    let p = ctx.params().clone();
    check_params(&p)?;
    let depth = ctx.adaptive_depth();
    let int = |k: &str| match p.get(k) {
        Some(elmer_core::Value::Integer(v)) if !v.is_empty() => Some(v[0] as i64),
        _ => None,
    };
    let min_depth = int("adaptive min depth").unwrap_or(0);
    if let Some(max_depth) = int("adaptive max depth") {
        if depth as i64 > max_depth {
            ctx.set_adaptive_finished(true);
            ctx.info(CALLER, "Max adaptive depth reached! Doing nothing!", 5);
            return Ok(());
        }
    }
    let init = depth == 0;
    if init {
        ctx.info(CALLER, "Initializing stuff on the coarsest level!", 4);
        state.do_final_ref = false;
    }
    if state.do_final_ref {
        ctx.info(CALLER, "Final refinement done. Nothing to do!", 6);
        ctx.set_adaptive_finished(true);
        return Ok(());
    }
    ctx.set_adaptive_finished(false);
    let adaptive_output = p.logical("adaptive output").unwrap_or(false);
    let (ei, max_error) = compute_error(ctx, est)?;
    let n = ei.len();
    let mut s = 0.0;
    for x in &ei {
        s += x.powi(2);
    }
    let error_estimate = (s / n as f64).sqrt();
    let nn = ctx.mesh.n_nodes();
    let nb = ctx.mesh.n_bulk();

    // Hvalue: the mean ElementDiameter at each node.
    field(ctx, "hvalue", adaptive_output, 0.0)?;
    let prev_h0 = values(ctx, "hvalue");
    let mut prev_h = if init { vec![0.0; nn] } else { prev_h0 };
    let mut hvalue = vec![0.0f64; nn];
    let mut referenced = vec![0usize; nn];
    for i in 0..nb {
        let el = ctx.mesh.element(ElemId(i as u32));
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
        let d = elmer_fe::element_diameter(el.etype, &xyz, false);
        for &k in el.nodes {
            hvalue[k as usize] += d;
            referenced[k as usize] += 1;
        }
    }
    for k in 0..nn {
        if referenced[k] > 0 {
            hvalue[k] /= referenced[k] as f64;
        }
    }
    let created = field(ctx, "hconvergence", adaptive_output, 1.0)?;
    let mut hconv = values(ctx, "hconvergence");
    if init && !created {
        hconv = vec![1.0; nn];
    }
    let name = ctx
        .variable_name()
        .ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?
        .to_string();
    let err_name = format!("{name}.error");
    let perr_name = format!("{name}.perror");
    let eref_name = format!("{name}.eref");
    field(ctx, &err_name, true, 0.0)?;
    field(ctx, &perr_name, adaptive_output, 0.0)?;
    let mut prev_nodal = values(ctx, &perr_name);
    if init {
        prev_nodal = vec![0.0; nn];
    }
    let mut nodal = vec![0.0f64; nn];
    let mut referenced = vec![0usize; nn];
    for i in 0..nb {
        for &k in ctx.mesh.element(ElemId(i as u32)).nodes {
            referenced[k as usize] += 1;
            nodal[k as usize] += ei[i];
        }
    }
    for k in 0..nn {
        if referenced[k] > 0 {
            nodal[k] /= referenced[k] as f64;
        }
    }
    let created = field(ctx, &eref_name, adaptive_output, 0.0)?;
    let mut eref = values(ctx, &eref_name);
    if init || created {
        eref = nodal.clone();
    }
    let created = field(ctx, "href", adaptive_output, 0.0)?;
    let mut href = values(ctx, "href");
    if init || created {
        href = hvalue.clone();
    }
    // `eRef = MAX( eRef, 1.0d-12 )`, `hRef = MAX( hRef, 1.0d-12 )` on the
    // stored arrays themselves.
    for x in &mut eref {
        *x = x.max(1.0e-12);
    }
    for x in &mut href {
        *x = x.max(1.0e-12);
    }
    store(ctx, &eref_name, eref.clone());
    store(ctx, "href", href.clone());
    store(ctx, &err_name, nodal.clone());

    let error_limit = p.real("adaptive error limit").unwrap_or(0.5);
    let max_scale = p.real("adaptive max error scale").unwrap_or(1.0);
    let ave_scale = p.real("adaptive average error scale").unwrap_or(1.0);
    let (max_frac, out_frac) = (1.0, 0.0);
    let mut limit = error_limit;
    let mut converged = false;
    if depth as i64 > min_depth {
        let conv = error_estimate < ave_scale * error_limit;
        if max_error < max_scale * error_limit && conv && out_frac < max_frac {
            match p.real("adaptive final refinement") {
                Some(f) => {
                    state.do_final_ref = true;
                    ctx.info(CALLER, "Performing one final refinement", 6);
                    limit = f * error_limit;
                }
                None => {
                    ctx.info(CALLER, "Mesh convergence limit reached. Nothing to do!", 6);
                    ctx.set_adaptive_finished(true);
                    converged = true;
                }
            }
        }
    }
    let mut new_mesh = None;
    if !converged {
        let mc = p
            .real("adaptive max change")
            .filter(|&m| m > AEPS)
            .unwrap_or(3.0);
        for k in 0..nn {
            if eref[k] > 0.0 {
                prev_nodal[k] += (hvalue[k] / href[k]).ln() * (nodal[k] / eref[k]).ln();
            }
        }
        for k in 0..nn {
            prev_h[k] += (hvalue[k] / href[k]).ln().powi(2);
        }
        if depth > 0 {
            for k in 0..nn {
                hconv[k] = if prev_h[k] > 0.0 {
                    (prev_nodal[k] / prev_h[k]).max(0.25)
                } else {
                    0.25
                };
            }
        }
        store(ctx, &perr_name, prev_nodal.clone());
        store(ctx, "hconvergence", hconv.clone());
        if p.logical("adaptive remesh") == Some(true) {
            let temperature = ctx.primary().is_some_and(|f| f.name.eq_ignore_ascii_case("temperature"));
            remesh::check(&p, ctx.simulation_params(), temperature)?;
            // The case directory of a bare `case.sif` is "": the working one.
            let case_dir = match p.string(OUTPUT_DIRECTORY_KEY).filter(|d| !d.is_empty()) {
                Some(d) => std::path::PathBuf::from(d),
                None => std::path::PathBuf::from("."),
            };
            ctx.info(CALLER, "Using file I/O for mesh refinement", 5);
            let m = remesh::external_remesh(
                ctx,
                &remesh::RemeshInput {
                    hvalue: &hvalue,
                    nodal_error: &nodal,
                    h_convergence: &hconv,
                    error_limit: limit / 3.0,
                    min_h: p.real("adaptive min h").unwrap_or(0.0),
                    max_h: p.real("adaptive max h").unwrap_or(0.0),
                    max_change: mc,
                    coarsening: p.logical("adaptive coarsening").unwrap_or(true),
                    depth,
                    case_dir: &case_dir,
                },
            )?;
            new_mesh = Some(m);
        } else {
            new_mesh = split::split_mesh(ctx.mesh, depth, &ei, limit, &hconv, mc)?;
        }
        // `Hvalue(1:nn) = PrevHValue(1:nn)`.
        store(ctx, "hvalue", prev_h.clone());
        if new_mesh.is_none() {
            ctx.info(CALLER, "Current mesh seems fine. Nothing to do.", 6);
        }
    }
    if let Some(m) = new_mesh {
        ctx.info(CALLER, "The new mesh consists of: ", 5);
        ctx.info(CALLER, &format!("Nodal points: {}", m.n_nodes()), 5);
        ctx.info(CALLER, &format!("Bulk elements: {}", m.n_bulk()), 5);
        ctx.info(CALLER, &format!("Boundary elements: {}", m.n_boundary()), 5);
        let base = p
            .string("adaptive mesh name")
            .unwrap_or("RefinedMesh")
            .to_string();
        let numbering = p.logical("adaptive mesh numbering").unwrap_or(true);
        let name = if numbering {
            format!("{base}{}", depth + 1)
        } else {
            base
        };
        if p.logical("adaptive save mesh") == Some(true) {
            // `OutputPath // '/' // NewMesh % Name`; elmer-run gives the
            // directory of the case.
            let base_dir = p.string(OUTPUT_DIRECTORY_KEY).unwrap_or(".").to_string();
            mesh::write_mesh_to_disk(&m, &std::path::Path::new(&base_dir).join(&name), &|t| t)
                .map_err(|e| SolverError::fatal("WriteMeshToDisk", e.to_string()))?;
        }
        ctx.request_mesh_change(MeshChange {
            mesh: Arc::new(m),
            name,
        });
    }
    ctx.info(
        CALLER,
        "----------- E N D   M E S H   R E F I N E M E N T --------------",
        5,
    );
    Ok(())
}
