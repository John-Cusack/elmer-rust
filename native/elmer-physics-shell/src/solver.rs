//! The `ShellSolver` procedure: `ShellSolver_Init0` and the solver of
//! ShellSolver.F90:154 for quadrilateral shells described by the mesh
//! (the blending surface of the element normals and the curved edges).

use std::collections::HashMap;

use elmer_core::{normalize_key, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_fe::basis::{nodal_basis, nodal_first_derivatives};
use elmer_fe::element_type;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::geom::{cross, dot, element_frame, hermite_form, norm, Frame, V3};
use crate::kernel::{shell_local_matrix, KernelInput, KernelOpts, M};
use crate::ShellError;

/// `Procedure = "ShellSolver" "ShellSolver"`.
pub const PROCEDURE: (&str, &str) = ("ShellSolver", "ShellSolver");

const CALLER: &str = "ShellSolver";
/// `AEPS = 10 * EPSILON(1.0_dp)`.
const AEPS: f64 = 10.0 * f64::EPSILON;

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(ShellSolver))
}

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("ShellSolver: {what}"))
}

impl From<ShellError> for SolverError {
    fn from(e: ShellError) -> Self {
        match e {
            ShellError::Unsupported(m) => unsupported(&m),
            ShellError::Fatal(m) => SolverError::fatal(CALLER, m),
        }
    }
}

fn add_new(p: &mut ParamMap, key: &str, v: Value) {
    if p.get(key).is_none() {
        p.insert(key, v);
    }
}

fn real(x: f64) -> Value {
    Value::Real {
        rows: 1,
        cols: 1,
        data: vec![x],
    }
}

/// `ShellSolver_Init0`: the six-component `Deflection[U:3 DNU:3]`
/// variable, `Large Deflection` (default on, with 50 nonlinear iterations
/// to 1e-5), `Time Derivative Order = 2`, `Skip Compute Nonlinear Change`
/// and the Simulation's `Initialize Dirichlet Conditions = False`. Call it
/// before the variables are created.
pub fn solver_defaults(p: &mut ParamMap, simulation: &mut ParamMap) {
    add_new(p, "Variable DOFs", Value::Integer(vec![6]));
    p.insert("Bubbles in Global System", Value::Logical(true));
    simulation.insert("Initialize Dirichlet Conditions", Value::Logical(false));
    if flag(p, "rotate dofs") {
        add_new(p, "Variable", Value::Str("Deflection[U:3 Theta:3]".into()));
    } else {
        add_new(p, "Variable", Value::Str("Deflection[U:3 DNU:3]".into()));
    }
    if flag(p, "eigen analysis") {
        p.insert("Large Deflection", Value::Logical(false));
    } else {
        add_new(p, "Large Deflection", Value::Logical(true));
        add_new(p, "Nonlinear System Max Iterations", Value::Integer(vec![50]));
        add_new(p, "Nonlinear System Convergence Tolerance", real(1.0e-5));
    }
    add_new(p, "Time derivative order", Value::Integer(vec![2]));
    add_new(p, "Skip Compute Nonlinear Change", Value::Logical(true));
    p.insert("Shell Solver", Value::Logical(true));
}

pub struct ShellSolver;

/// `ReadSurfaceDirector`: the unit nodal directors of `mesh.director`
/// (`k d1 d2 d3` per line, `k` 1-based), by node.
fn read_directors(path: &str) -> Result<HashMap<u32, V3>, SolverError> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| SolverError::fatal("ReadSurfaceDirector", format!("Opening mesh.director file failed: {e}")))?;
    let mut out = HashMap::new();
    for line in text.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        if w.len() < 4 {
            continue;
        }
        let bad = || SolverError::fatal("ReadSurfaceDirector", format!("bad line in mesh.director: {line}"));
        let k: u32 = w[0].parse().map_err(|_| bad())?;
        let mut d = [0.0f64; 3];
        for i in 0..3 {
            d[i] = w[1 + i].replace(['d', 'D'], "E").parse().map_err(|_| bad())?;
        }
        let nr = norm(&d);
        out.insert(k.saturating_sub(1), [d[0] / nr, d[1] / nr, d[2] / nr]);
    }
    Ok(out)
}

/// The geometry of one background element.
struct ElementGeometry {
    frame: Frame,
}

/// `NormalVector( Element, Nodes, Check=.TRUE. )` of a four-node face: the
/// unit normal at the centre, turned away from the parent body
/// (`CheckNormalDirection`).
fn checked_normal(ctx: &SolverCtx<'_>, e: ElemId) -> Result<V3, SolverError> {
    let mesh = ctx.mesh;
    let el = mesh.element(e);
    let def = element_type(el.etype.0).map_err(|m| SolverError::fatal("NormalVector", m.to_string()))?;
    let n = el.nodes.len();
    let xyz: Vec<V3> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
    let mut d = vec![[0.0f64; 3]; n];
    nodal_first_derivatives(def, 0.0, 0.0, 0.0, &mut d);
    let col = |k: usize, c: usize| -> f64 {
        let mut s = 0.0;
        for i in 0..n {
            s += d[i][k] * xyz[i][c];
        }
        s
    };
    let (dxdu, dydu, dzdu) = (col(0, 0), col(0, 1), col(0, 2));
    let (dxdv, dydv, dzdv) = (col(1, 0), col(1, 1), col(1, 2));
    let auu = dxdu * dxdu + dydu * dydu + dzdu * dzdu;
    let auv = dxdu * dxdv + dydu * dydv + dzdu * dzdv;
    let avv = dxdv * dxdv + dydv * dydv + dzdv * dzdv;
    let deta = 1.0 / (auu * avv - auv * auv).sqrt();
    let mut normal = [
        (dydu * dzdv - dydv * dzdu) * deta,
        (dxdv * dzdu - dxdu * dzdv) * deta,
        (dxdu * dydv - dxdv * dydu) * deta,
    ];
    // The boundary's centre.
    let mut basis = vec![0.0f64; n];
    nodal_basis(def, 0.0, 0.0, 0.0, &mut basis);
    let centre = |pts: &[V3], b: &[f64]| -> V3 {
        let mut c = [0.0f64; 3];
        for k in 0..3 {
            let mut s = 0.0;
            for i in 0..pts.len() {
                s += b[i] * pts[i][k];
            }
            c[k] = s;
        }
        c
    };
    let x = centre(&xyz, &basis);
    // CheckNormalDirection: the parent the normal points away from.
    let model = &ctx.sim.model;
    let out_body = model
        .bc_of_tag(el.tag)
        .and_then(|b| model.boundary_conditions[b - 1].integer("normal target body"))
        .unwrap_or(0) as usize;
    let body = |p: ElemId| model.body_of_tag(mesh.element(p).tag).unwrap_or(0);
    let parent = match el.parents {
        [Some(l), Some(r)] => {
            let (bl, br) = (body(l), body(r));
            if out_body > 0 {
                Some(if bl == out_body { r } else { l })
            } else if bl > br {
                Some(l)
            } else if bl < br {
                Some(r)
            } else {
                return Err(unsupported("director of a face between elements of one body"));
            }
        }
        [Some(l), None] => Some(l),
        [None, r] => r,
    };
    let Some(parent) = parent else { return Ok(normal) };
    let pel = mesh.element(parent);
    let pdef = element_type(pel.etype.0).map_err(|m| SolverError::fatal("CheckNormalDirection", m.to_string()))?;
    let uvw = match pel.etype.0 / 100 {
        2 | 4 | 8 => [0.0, 0.0, 0.0],
        3 | 7 => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        5 => [1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0],
        6 => [0.0, 0.0, 1.0 / 3.0],
        _ => return Err(SolverError::fatal("CheckNormalDirection", "Invalid elementcode for parent element!")),
    };
    let pxyz: Vec<V3> = pel.nodes.iter().map(|&k| mesh.coord(k)).collect();
    let mut pb = vec![0.0f64; pel.nodes.len()];
    nodal_basis(pdef, uvw[0], uvw[1], uvw[2], &mut pb);
    let pc = centre(&pxyz, &pb);
    let dc = [pc[0] - x[0], pc[1] - x[1], pc[2] - x[2]];
    let mut s = 0.0;
    for k in 0..3 {
        s += dc[k] * normal[k];
    }
    let pbody = body(parent);
    if (s > 0.0 && pbody != out_body) || (s <= 0.0 && pbody == out_body) {
        normal = [-normal[0], -normal[1], -normal[2]];
    }
    Ok(normal)
}

/// The geometry of every active quadrilateral: `GetElementalDirector`
/// (the face normal at every node), `CheckSurfaceOrientation`,
/// `CreateCurvedEdges`, `AverageDirector` and the lines-of-curvature frame.
fn element_geometry(
    ctx: &SolverCtx<'_>,
    elems: &[ElemId],
    directors: Option<&HashMap<u32, V3>>,
) -> Result<HashMap<u32, ElementGeometry>, SolverError> {
    let mut out = HashMap::new();
    for &e in elems {
        let el = ctx.mesh.element(e);
        if el.etype.0 != 404 {
            return Err(unsupported(&format!("background element type {}", el.etype.0)));
        }
        let nodes: [V3; 4] = std::array::from_fn(|i| ctx.mesh.coord(el.nodes[i]));
        let dirs: [V3; 4] = match directors {
            Some(map) => {
                let mut d = [[0.0f64; 3]; 4];
                for i in 0..4 {
                    d[i] = *map.get(&el.nodes[i]).ok_or_else(|| {
                        SolverError::fatal("ReadSurfaceDirector", format!("no director for node {}", el.nodes[i] + 1))
                    })?;
                }
                d
            }
            None => [checked_normal(ctx, e)?; 4],
        };
        // CheckSurfaceOrientation: the directors on one side of the face.
        let e1 = {
            let v = [nodes[1][0] - nodes[0][0], nodes[1][1] - nodes[0][1], nodes[1][2] - nodes[0][2]];
            let nr = norm(&v);
            [v[0] / nr, v[1] / nr, v[2] / nr]
        };
        let e2 = [nodes[2][0] - nodes[1][0], nodes[2][1] - nodes[1][1], nodes[2][2] - nodes[1][2]];
        let e3 = cross(e1, e2);
        let nr = norm(&e3);
        let e3 = [e3[0] / nr, e3[1] / nr, e3[2] / nr];
        let c = dot(&dirs[0], &e3);
        if dirs.iter().skip(1).any(|d| c * dot(d, &e3) < 0.0) {
            return Err(SolverError::fatal(
                "CheckSurfaceOrientation",
                "Director data does not define a unique upper/lower surface.",
            ));
        }
        // CreateCurvedEdges: the Hermite curves of the four edges.
        let pairs = [(0usize, 1usize), (1, 2), (3, 2), (0, 3)];
        let edge: [[f64; 6]; 4] = std::array::from_fn(|k| {
            let (i, j) = pairs[k];
            hermite_form(nodes[i], nodes[j], dirs[i], dirs[j])
        });
        // AverageDirector with the plate check.
        let mut d = [0.0f64; 3];
        for dir in &dirs {
            for k in 0..3 {
                d[k] += dir[k];
            }
        }
        let d = [1.0 / 4.0 * d[0], 1.0 / 4.0 * d[1], 1.0 / 4.0 * d[2]];
        let eps = f32::EPSILON as f64;
        let mut plate = true;
        for dir in dirs.iter().skip(1) {
            let m = (0..3).map(|k| (dirs[0][k] - dir[k]).abs()).fold(0.0f64, f64::max);
            plate = !(m > eps);
            if !plate {
                break;
            }
        }
        let spherical = ctx
            .body_params(e)
            .is_some_and(|b| b.body.logical("spherical body") == Some(true));
        let frame = element_frame(&nodes, &edge, d, plate, spherical)?;
        out.insert(e.0, ElementGeometry { frame });
    }
    Ok(out)
}

/// The keyword at the element's nodes, `0` when the section lacks it.
fn at_nodes(ctx: &SolverCtx<'_>, p: Option<&ParamMap>, key: &str, e: ElemId) -> Result<Option<Vec<f64>>, SolverError> {
    match p {
        Some(p) => ctx.real_at_nodes(p, key, e),
        None => Ok(None),
    }
}

impl Solver for ShellSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        if p.string("element").is_some_and(|s| s.contains("p:")) {
            return Err(unsupported("p elements (Cartesian Formulation)"));
        }
        for k in [
            "cartesian formulation",
            "skip surface reconstruction",
            "mesh reparametrization",
            "drilling dofs",
            "rotate dofs",
            "bubbles",
            "principal axes output",
            "eigen analysis",
            "harmonic mode",
            "harmonic analysis",
            "displace mesh",
            "stability analysis",
            "geometric stiffness",
        ] {
            if flag(p, k) {
                return Err(unsupported(k));
            }
        }
        for k in ["strain reduction operator", "membrane strain reduction operator"] {
            if p.get(k).is_some() {
                return Err(unsupported(k));
            }
        }
        if p.integer("variable dofs").unwrap_or(6) != 6 {
            return Err(unsupported("Variable DOFs other than 6"));
        }
        if ctx.sim.vars.get("director").is_some() {
            return Err(unsupported("a Director field"));
        }
        if ctx.sim.vars.get("displacement").is_some() {
            return Err(unsupported("shell-solid coupling (a Displacement field)"));
        }
        let simtype = ctx.simulation_params().string("simulation type").map(normalize_key);
        if simtype.as_deref() == Some("transient") {
            return Err(unsupported("transient shells (mass matrix)"));
        }
        if ctx.comm.size() > 1 {
            return Err(unsupported("more than one partition"));
        }
        // ShellBoundaryMatrix: loads on edges.
        let nb = ctx.mesh.n_bulk();
        for t in nb..nb + ctx.mesh.n_boundary() {
            let el = ctx.mesh.element(ElemId(t as u32));
            if el.etype.0 / 100 > 2 {
                continue;
            }
            if let Some(bc) = ctx.bc_params(ElemId(t as u32)) {
                if bc.0.keys().any(|k| {
                    ["force", "load", "couple", "moment", "traction", "pressure"]
                        .iter()
                        .any(|w| k.contains(w))
                }) {
                    return Err(unsupported("edge loads (ShellBoundaryMatrix)"));
                }
            }
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        let ld = flag(&p, "large deflection");
        let max_iters = if ld {
            p.integer("nonlinear system max iterations").unwrap_or(0).max(0) as usize
        } else {
            0
        };
        let tol = ctx.creal(&p, "nonlinear system convergence tolerance")?.unwrap_or(0.0);
        let relaxation = ctx.creal(&p, "nonlinear system relaxation factor")?;
        let relax = relaxation.filter(|r| (r - 1.0).abs() > f64::EPSILON);
        let opts = KernelOpts {
            large_deflection: ld,
            nonlinear_bending: p.logical("nonlinear bending strains").unwrap_or(true),
            shear_alpha: ctx.creal(&p, "shear relaxation alpha")?.unwrap_or(1.0),
            membrane_alpha: ctx.creal(&p, "membrane relaxation alpha")?.unwrap_or(0.0),
            stretch_alpha: ctx.creal(&p, "stretch relaxation alpha")?.unwrap_or(0.0),
            benchmark: !ld && flag(&p, "benchmark problem"),
        };
        let (perm, nsize) = {
            let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
            (v.perm.ok_or_else(|| SolverError::fatal(CALLER, "no Perm"))?.to_vec(), v.values.len())
        };
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        for &e in &elems {
            let fam = ctx.mesh.element(e).etype.0 / 100;
            if fam == 2 {
                return Err(unsupported("beam elements (BeamStiffnessMatrix)"));
            }
        }
        let shells: Vec<ElemId> = elems
            .iter()
            .copied()
            .filter(|&e| matches!(ctx.mesh.element(e).etype.0 / 100, 3 | 4))
            .collect();
        let directors = match p.string("mesh director file") {
            Some(f) => Some(read_directors(f)?),
            None => None,
        };
        let geometry = element_geometry(ctx, &shells, directors.as_ref())?;
        let fe = ctx
            .sim
            .backends
            .fe
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
        let points: Vec<(f64, f64, f64)> = fe
            .quadrature
            .gauss_points(elmer_core::ElementType(404), None)
            .iter()
            .map(|q| (q.uvw[0], q.uvw[1], q.weight))
            .collect();

        let mut bulk_rhs = vec![0.0f64; nsize];
        let mut res0 = 1.0;
        for iter in 1..=max_iters + 1 {
            if iter < max_iters + 1 {
                ctx.info(CALLER, &format!("Nonlinear iteration: {iter:4}"), 4);
            }
            let total: Vec<f64> = ctx.primary().expect("checked").values.to_vec();
            ctx.default_initialize()?;
            for &e in &shells {
                let nodes = ctx.mesh.element(e).nodes.to_vec();
                let mut sol = vec![0.0f64; M * nodes.len()];
                if ld {
                    for (i, &k) in nodes.iter().enumerate() {
                        let s = perm[k as usize];
                        if s == NO_DOF {
                            continue;
                        }
                        for c in 0..M {
                            sol[i * M + c] = total[s as usize * M + c];
                        }
                    }
                }
                let bp = ctx.body_params(e);
                let material = bp.and_then(|b| b.material);
                let need = |key: &str| -> Result<Vec<f64>, SolverError> {
                    at_nodes(ctx, material, key, e)?
                        .ok_or_else(|| SolverError::fatal("ShellLocalMatrix", format!("{key} not found")))
                };
                let poisson = need("poisson ratio")?;
                let youngs = need("youngs modulus")?;
                let thickness = match at_nodes(ctx, material, "shell thickness", e)? {
                    Some(t) => t,
                    None => at_nodes(ctx, material, "thickness", e)?.ok_or_else(|| {
                        SolverError::fatal("ShellLocalMatrix", "Shell Thickness/Thickness undefined")
                    })?,
                };
                let pressure = at_nodes(ctx, bp.and_then(|b| b.body_force), "normal pressure", e)?
                    .unwrap_or_else(|| vec![0.0; nodes.len()]);
                let g = &geometry[&e.0];
                let out = shell_local_matrix(
                    &g.frame,
                    &KernelInput {
                        thickness: &thickness,
                        poisson: &poisson,
                        youngs: &youngs,
                        pressure: &pressure,
                        local_sol: &sol,
                        points: &points,
                    },
                    &opts,
                )?;
                let n = out.force.len();
                let ls = LocalSystem {
                    n,
                    stiff: out.stiff,
                    mass: None,
                    damp: None,
                    force: out.force,
                };
                ctx.default_update(&ls, e)?;
                if ld && iter == 1 {
                    for (i, &k) in nodes.iter().enumerate() {
                        let s = perm[k as usize];
                        if s == NO_DOF {
                            continue;
                        }
                        for c in 0..M {
                            bulk_rhs[s as usize * M + c] += out.rhs_force[i * M + c];
                        }
                    }
                }
            }
            ctx.default_finish_bulk_assembly()?;
            ctx.default_finish_boundary_assembly()?;
            // DefaultDirichletBCs, then DValues - Values: the increment.
            let (values, rows) = ctx.dirichlet_values_with_rows()?;
            let cur = ctx.primary().expect("checked").values.to_vec();
            let fixed: Vec<(u32, f64)> = values.iter().map(|&(k, v)| (k, v - cur[k as usize])).collect();
            ctx.apply_dirichlet_values_with_rows(&fixed, &rows)?;
            if ld {
                if iter == 1 {
                    let mx = bulk_rhs.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                    if mx < AEPS {
                        ctx.info(CALLER, "No distributed loads ... ", 4);
                        res0 = 1.0;
                    } else {
                        let mut s = 0.0;
                        for x in &bulk_rhs {
                            s += x.powi(2);
                        }
                        res0 = s.sqrt();
                    }
                }
                let rhs = ctx.system().map(|s| s.rhs.clone()).unwrap_or_default();
                let mut s = 0.0;
                for x in &rhs {
                    s += x.powi(2);
                }
                let res = s.sqrt() / res0;
                ctx.info(CALLER, &format!("Residual for nonlinear iterate{:4} {res:e}", iter - 1), 3);
                if res < tol {
                    ctx.info(CALLER, "Nonlinear iteration is terminated succesfully", 3);
                    break;
                } else if iter > max_iters {
                    if flag(&p, "nonlinear system abort not converged") {
                        return Err(SolverError::fatal(CALLER, "Nonlinear iteration did not converge to tolerance"));
                    }
                    ctx.info(CALLER, "Nonlinear iteration did not converge to tolerance", 6);
                    break;
                }
                ctx.primary_mut().expect("checked").values.fill(0.0);
            }
            ctx.default_solve()?;
            if ld {
                let v = ctx.primary_mut().expect("checked");
                for (x, t) in v.values.iter_mut().zip(&total) {
                    *x = match relax {
                        Some(r) => t + r * *x,
                        None => t + *x,
                    };
                }
            }
            if !ld {
                break;
            }
        }
        Ok(())
    }
}
