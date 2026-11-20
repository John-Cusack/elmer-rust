//! `ReynoldsSolver` and `ReynoldsPostprocess` of the `ReynoldsSolver`
//! module (fem/src/modules/ReynoldsSolver.F90): the Reynolds equation of a
//! thin film on 1D or 2D elements, and the fields that follow from its
//! pressure.
//!
//! **ReynoldsSolver** (`LocalBulkMatrix`). With density `rho(p)`
//! (`Compressibility Model`: incompressible, weakly compressible
//! `exp(p/K)`, isothermal or adiabatic ideal gas), gap `h` and viscosity
//! `mu` (`Viscosity Model` newtonian or rarefied):
//! `div(-rho h^3/(12 mu) grad p) - 1/2 rho' h u . grad p - Y rho p
//!  = rho v_n + 1/2 rho div(h u) - Y rho p_ext`,
//! the surface velocity `u` split into its normal and tangential parts,
//! `Y` the `Flow Admittance`. Boundaries (`LocalBoundaryMatrix`): `Open
//! Side`, `Filmpressure Flux`, `Filmpressure Velocity` and the Robin pair
//! `External FilmPressure` / `Filmpressure Transfer Coefficient`. The
//! nonlinear iterations run until `NonlinConverged`.
//!
//! As in the Fortran: the body force is looked up through `Body(<equation
//! id>)`, and the boundary terms use the compressibility, viscosity model
//! and element pressure left by the last bulk element.
//!
//! **ReynoldsPostprocess**: for `Calculate Force`, `Calculate Flux`,
//! `Calculate Mean Velocity` and `Calculate Heating`, a mass-matrix
//! projection of each component (`LocalMatrix`), solved with the solver's
//! own variable (`<p> Heating`), copied to `<p> Force` / `Flux` / `Mean
//! Velocity`; the integrated pressure and sliding forces and the heating go
//! to the Simulation as `res:` values.
//!
//! Refused: bubbles, `Manning Model`, `Gap Sensitivity`, `Calculate
//! Sensitivity`, artificial compressibility, line searches, non-Cartesian
//! coordinates, p elements, and a surface velocity with a normal component
//! on elements outside the `xy` plane or on line elements.

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "ReynoldsSolver";
const POST: &str = "ReynoldsPostprocess";

pub fn solver_factory(p: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    for key in ["bubbles", "manning model", "gap sensitivity", "calculate sensitivity", "nonlinear system linesearch"] {
        if p.logical(key) == Some(true) {
            return Err(SolverError::Unsupported(format!("ReynoldsSolver: {key}")));
        }
    }
    Ok(Box::new(Reynolds))
}

pub fn post_factory(p: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    if p.logical("manning model") == Some(true) {
        return Err(SolverError::Unsupported("ReynoldsPostprocess: Manning Model".into()));
    }
    Ok(Box::new(ReynoldsPost))
}

pub struct Reynolds;
pub struct ReynoldsPost;

#[derive(Clone, Copy, PartialEq)]
enum Compress {
    None,
    Weak(f64),
    GasIsothermal,
    GasAdiabatic(f64),
}

#[derive(Clone, Copy)]
enum Visc {
    Newtonian,
    Rarefied(f64),
}

fn fatal(m: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, m)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |s, (x, y)| s + x * y)
}

fn has_prefix(p: Option<&ParamMap>, prefix: &str) -> bool {
    p.is_some_and(|p| p.0.keys().any(|k| k.starts_with(prefix)))
}

/// The sections of a bulk element: material and the body force of
/// `Body(<equation id>)` (the Fortran's `Model % Bodies(ent_id)`).
fn sections(ctx: &SolverCtx<'_>, e: ElemId) -> Result<(usize, ParamMap, Option<ParamMap>), SolverError> {
    let model = &ctx.sim.model;
    let body_id = model.body_of_tag(ctx.mesh.element(e).tag).ok_or_else(|| fatal("an element without a body"))?;
    let body = &model.bodies[body_id - 1];
    let eq_id = body.integer("equation").ok_or_else(|| fatal("a body without an Equation"))? as usize;
    let mat_id = body.integer("material").ok_or_else(|| fatal("a body without a Material"))? as usize;
    let mat = model.materials.get(mat_id - 1).cloned().ok_or_else(|| fatal("no such Material"))?;
    let bf = model
        .bodies
        .get(eq_id.wrapping_sub(1))
        .and_then(|b| b.integer("body force"))
        .filter(|&i| i > 0)
        .and_then(|i| model.body_forces.get(i as usize - 1))
        .cloned();
    Ok((mat_id, mat, bf))
}

fn nodal(ctx: &SolverCtx<'_>, p: Option<&ParamMap>, key: &str, e: ElemId, n: usize) -> Result<(Vec<f64>, bool), SolverError> {
    match p {
        Some(p) => Ok(match ctx.real_at_nodes(p, key, e)? {
            Some(v) => (v, true),
            None => (vec![0.0; n], false),
        }),
        None => Ok((vec![0.0; n], false)),
    }
}

/// The nodal surface or tangent velocity and `UseVelocity`.
fn velocity(
    ctx: &SolverCtx<'_>,
    bf: Option<&ParamMap>,
    mat: Option<&ParamMap>,
    e: ElemId,
    n: usize,
) -> Result<(Vec<[f64; 3]>, bool), SolverError> {
    let mut velo = vec![[0.0; 3]; n];
    let mut got = [false; 3];
    let read = |p: Option<&ParamMap>, base: &str, velo: &mut Vec<[f64; 3]>, got: &mut [bool; 3]| -> Result<(), SolverError> {
        for i in 0..3 {
            let (v, g) = nodal(ctx, p, &format!("{base} {}", i + 1), e, n)?;
            for (o, x) in velo.iter_mut().zip(v) {
                o[i] = x;
            }
            got[i] = g;
        }
        Ok(())
    };
    let mut use_velocity = false;
    if has_prefix(bf, "surface velocity") {
        read(bf, "surface velocity", &mut velo, &mut got)?;
        use_velocity = got.iter().any(|&g| g);
    }
    if !use_velocity && has_prefix(mat, "surface velocity") {
        read(mat, "surface velocity", &mut velo, &mut got)?;
        use_velocity = got.iter().any(|&g| g);
    }
    if !use_velocity {
        if has_prefix(bf, "tangent velocity") {
            read(bf, "tangent velocity", &mut velo, &mut got)?;
        }
        if !got.iter().any(|&g| g) && has_prefix(mat, "tangent velocity") {
            read(mat, "tangent velocity", &mut velo, &mut got)?;
        }
    }
    Ok((velo, use_velocity))
}

/// `NormalVector(Element, Nodes, u, v, .FALSE.)` of a bulk element in the
/// `xy` plane: `(0, 0, ±1)` by the element's orientation. `None` for other
/// elements.
/// The sign is that of the corner polygon's area, which a valid element
/// shares with its Jacobian `dx/du dy/dv - dx/dv dy/du`.
fn plane_normal(xyz: &[[f64; 3]], code: u16) -> Option<[f64; 3]> {
    if element_dim(code) != 2 || xyz.iter().any(|c| c[2] != 0.0) {
        return None;
    }
    let k = if code / 100 == 3 { 3 } else { 4 };
    let area: f64 = (0..k).map(|i| {
        let (a, b) = (xyz[i], xyz[(i + 1) % k]);
        a[0] * b[1] - b[0] * a[1]
    }).sum();
    Some([0.0, 0.0, if area < 0.0 { -1.0 } else { 1.0 }])
}

fn element_dim(code: u16) -> usize {
    match code / 100 {
        1 => 0,
        2 => 1,
        3 | 4 => 2,
        _ => 3,
    }
}

/// The material's models, read when the material changes.
fn models(mat: &ParamMap, linear: bool) -> Result<(f64, Visc, Compress), SolverError> {
    let reference = mat.real("reference pressure").unwrap_or(0.0);
    let visc = match mat.string("viscosity model").map(elmer_core::normalize_key).as_deref() {
        None | Some("newtonian") => Visc::Newtonian,
        Some("rarefied") => Visc::Rarefied(mat.real("mean free path").ok_or_else(|| fatal("Mean Free Path not found"))?),
        Some(m) => return Err(SolverError::Unsupported(format!("ReynoldsSolver: Viscosity Model {m}"))),
    };
    let mut compress = Compress::None;
    if !linear {
        if let Some(m) = mat.string("compressibility model") {
            compress = match elmer_core::normalize_key(m).trim_end_matches('"') {
                "incompressible" => Compress::None,
                "weakly compressible" => Compress::Weak(mat.real("bulk modulus").ok_or_else(|| fatal("Bulk Modulus not found"))?),
                "isothermal ideal gas" => Compress::GasIsothermal,
                "adiabatic ideal gas" => {
                    Compress::GasAdiabatic(mat.real("specific heat ratio").ok_or_else(|| fatal("Specific Heat Ratio not found"))?)
                }
                "artificial compressible" => {
                    return Err(SolverError::Unsupported("ReynoldsSolver: artificial compressibility".into()))
                }
                other => return Err(fatal(format!("Unknown compressibility model: {other}"))),
            };
        }
    }
    Ok((reference, visc, compress))
}

fn rarefied(visc: f64, mfp0: f64, reference: f64, tot: f64, gap: f64) -> (f64, f64) {
    let mfp = mfp0 * reference / tot;
    let kn = mfp / gap.abs();
    (visc / (1.0 + 9.638 * kn.powf(1.159)), kn)
}

/// The scalar values of `var` at the element's nodes.
fn local_values(values: &[f64], perm: Option<&[u32]>, nodes: &[u32]) -> Vec<f64> {
    nodes
        .iter()
        .map(|&k| {
            let s = perm.map_or(k, |p| p[k as usize]);
            if s == NO_DOF {
                0.0
            } else {
                values[s as usize]
            }
        })
        .collect()
}

impl Solver for Reynolds {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "---------------------------------------", 5);
        ctx.info(CALLER, if transient { "Solving the transient Reynolds equation" } else { "Solving the steady-state Reynolds equation" }, 5);
        ctx.info(CALLER, "---------------------------------------", 5);
        let p = ctx.params().clone();
        if ctx.simulation_params().string("coordinate system").is_some_and(|c| !elmer_core::normalize_key(c).starts_with("cartesian")) {
            return Err(SolverError::Unsupported("ReynoldsSolver: non-Cartesian coordinates".into()));
        }
        let (dofs, perm) = {
            let v = ctx.primary().ok_or_else(|| fatal("no variable"))?;
            (v.dofs, v.perm.map(|p| p.to_vec()))
        };
        if dofs != 1 {
            return Err(fatal("Impossible number of dofs! (should be 1)"));
        }
        if perm.as_ref().is_some_and(|p| p.iter().all(|&s| s == NO_DOF)) {
            return Ok(());
        }
        if p.logical("skip first solution") == Some(true) {
            return Err(SolverError::Unsupported("ReynoldsSolver: Skip First Solution".into()));
        }
        let model = &ctx.sim.model;
        let any_bc = model.boundary_conditions.iter().any(|b| {
            b.logical("open side") == Some(true)
                || ["filmpressure flux", "filmpressure velocity", "filmpressure transfer coefficient"].iter().any(|k| b.get(k).is_some())
        });
        let min_gap = p.real("min gap height");
        let fsi_mult = match p.real("fsi velocity multiplier") {
            _ if ctx.scalar("timestep").unwrap_or(1.0).round() as i32 == 1 && p.logical("steady state fsi start") == Some(true) => 0.0,
            Some(m) => m,
            None => 1.0,
        };
        let iterations = p.integer("nonlinear system max iterations").unwrap_or(1).max(1);
        let cdim = ctx.mesh.dim();
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut be = BasisEval::default();
        ctx.info(CALLER, "-------------------------------------------------", 5);
        for iter in 1..=iterations {
            let linear = iter == 1 && p.logical("linear first iteration") == Some(true);
            ctx.info(CALLER, &format!("{:<34}{iter:>5}", "Reynolds iteration:"), 5);
            ctx.default_initialize()?;
            let pres_all: Vec<f64> = ctx.primary().expect("checked").values.to_vec();
            let mut mat_old = None;
            let (mut reference, mut visc_model, mut compress) = (0.0, Visc::Newtonian, Compress::None);
            let mut last_pres: Vec<f64> = Vec::new();
            for &e in &elems {
                let el = ctx.mesh.element(e);
                if el.etype.0 > 500 {
                    return Err(fatal("This is a reduced dimensional solver for 1D and 2D only!"));
                }
                let nodes = el.nodes.to_vec();
                let etype = el.etype;
                let n = nodes.len();
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let pres = local_values(&pres_all, perm.as_deref(), &nodes);
                let (mat_id, mat, bf) = sections(ctx, e)?;
                let (velo, use_velocity) = velocity(ctx, bf.as_ref(), Some(&mat), e, n)?;
                let (mut vn, got) = nodal(ctx, bf.as_ref(), "normal velocity", e, n)?;
                if !got {
                    vn = nodal(ctx, Some(&mat), "normal velocity", e, n)?.0;
                }
                let (fsi, got) = nodal(ctx, bf.as_ref(), "fsi velocity", e, n)?;
                vn.iter_mut().zip(&fsi).for_each(|(a, b)| *a += fsi_mult * b);
                if !got {
                    let fsi = nodal(ctx, Some(&mat), "fsi velocity", e, n)?.0;
                    vn.iter_mut().zip(&fsi).for_each(|(a, b)| *a += fsi_mult * b);
                }
                let mut gap = nodal(ctx, Some(&mat), "gap height", e, n)?.0;
                if let Some(m) = min_gap {
                    gap.iter_mut().for_each(|g| *g = g.max(m));
                }
                let adm = nodal(ctx, Some(&mat), "flow admittance", e, n)?.0;
                let ext = nodal(ctx, bf.as_ref(), "external filmpressure", e, n)?.0;
                let visc_n = nodal(ctx, Some(&mat), "viscosity", e, n)?.0;
                if mat_old != Some(mat_id) {
                    mat_old = Some(mat_id);
                    (reference, visc_model, compress) = models(&mat, linear)?;
                }
                let normal = if use_velocity {
                    Some(plane_normal(&xyz, etype.0).ok_or_else(|| {
                        SolverError::Unsupported("ReynoldsSolver: a surface velocity on line elements or elements outside the xy plane".into())
                    })?)
                } else {
                    None
                };
                let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| fatal("no FE backend"))?;
                let dim = cdim;
                let mut local = LocalSystem::zeros(n);
                let mut mass = vec![0.0; n * n];
                for qp in fe.quadrature.gauss_points(etype, None) {
                    fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let (b, g) = (&be.basis[..n], &be.dbasis_dx[..n]);
                    let s = qp.weight * be.det_j;
                    let mut normal_velo = dot(b, &vn);
                    let tangent: [f64; 3] = match normal {
                        Some(nv) => {
                            let v: [f64; 3] = std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + b[k] * velo[k][i]));
                            normal_velo += nv.iter().zip(&v).fold(0.0, |acc, (x, y)| acc + x * y);
                            std::array::from_fn(|i| v[i] - normal_velo * nv[i])
                        }
                        None => std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + b[k] * velo[k][i])),
                    };
                    let damp = dot(b, &adm);
                    let ext_p = dot(b, &ext);
                    let pr = dot(b, &pres);
                    let gp = dot(b, &gap);
                    let tot = reference + pr;
                    let mut visc = dot(b, &visc_n);
                    if let Visc::Rarefied(mfp0) = visc_model {
                        visc = rarefied(visc, mfp0, reference, tot, gp).0;
                    }
                    let (density, density_der) = match compress {
                        Compress::None => (1.0, 0.0),
                        Compress::Weak(k) => {
                            let d = (tot / k).exp();
                            (d, d / k)
                        }
                        Compress::GasIsothermal => (tot, 1.0),
                        Compress::GasAdiabatic(r) => (tot.powf(1.0 / r), (1.0 / r) * tot.powf(1.0 / r - 1.0)),
                    };
                    let ms = -density * (gp * gp * gp) / (12.0 * visc);
                    let hr = match compress {
                        Compress::GasIsothermal => -damp * (density + ext_p) / 2.0,
                        Compress::GasAdiabatic(r) => -damp * (density + ext_p.powf(1.0 / r)) / 2.0,
                        _ => -damp * density,
                    };
                    let mm = -density_der * gp;
                    let ma = 0.0;
                    let l = density * normal_velo;
                    let mut slr = 0.0;
                    for i in 0..dim {
                        let sum = (0..n).fold(0.0, |acc, k| acc + g[k][i] * velo[k][i] * gap[k]);
                        slr += 0.5 * density * sum;
                    }
                    let sll: [f64; 3] = std::array::from_fn(|i| -0.5 * density_der * gp * tangent[i]);
                    for pp in 0..n {
                        for q in 0..n {
                            let mut a = (ma + hr) * b[q] * b[pp];
                            for i in 0..dim {
                                for j in 0..dim {
                                    let metric = if i == j { 1.0 } else { 0.0 };
                                    a += ms * metric * g[q][i] * g[pp][j];
                                    a += sll[j] * metric * g[q][i] * b[pp];
                                }
                            }
                            local.stiff[pp * n + q] += s * a;
                            if transient {
                                mass[pp * n + q] += s * (mm * b[q] * b[pp]);
                            }
                        }
                        let f = l + slr + hr * ext_p;
                        local.force[pp] += s * b[pp] * f;
                    }
                }
                if transient {
                    local.mass = Some(mass);
                    ctx.default_first_order_time(&mut local, e, dt)?;
                    local.mass = None;
                }
                ctx.default_update(&local, e)?;
                last_pres = pres;
            }
            ctx.default_finish_bulk_assembly()?;
            if any_bc {
                self.boundary(ctx, perm.as_deref(), min_gap, (reference, visc_model, compress), &last_pres, dt, transient)?;
            }
            ctx.default_finish_boundary_assembly()?;
            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.convergence().is_some_and(|c| c.nonlin_converged > 0) {
                break;
            }
        }
        ctx.info(CALLER, "-------------------------------------------------", 5);
        Ok(())
    }
}

impl Reynolds {
    /// `GlobalBoundaryAssemby`: `models` and `last_pres` are what the last
    /// bulk element left (`ElemPressure` is not re-read here).
    #[allow(clippy::too_many_arguments)]
    fn boundary(
        &self,
        ctx: &mut SolverCtx<'_>,
        perm: Option<&[u32]>,
        min_gap: Option<f64>,
        models: (f64, Visc, Compress),
        last_pres: &[f64],
        dt: f64,
        transient: bool,
    ) -> Result<(), SolverError> {
        let (reference, visc_model, compress) = models;
        let in_perm = |k: u32| perm.is_none_or(|p| p[k as usize] != NO_DOF);
        let bes: Vec<ElemId> = ctx.boundary_elements().collect();
        let mut be = BasisEval::default();
        for e in bes {
            let el = ctx.mesh.element(e);
            let nodes = el.nodes.to_vec();
            let etype = el.etype;
            let parents = el.parents;
            let n = nodes.len();
            let Some(bc) = ctx.bc_params(e).cloned() else { continue };
            let open = bc.logical("open side") == Some(true);
            let (flux, got_flux) = nodal(ctx, Some(&bc), "filmpressure flux", e, n)?;
            let (vpres, got_velo) = nodal(ctx, Some(&bc), "filmpressure velocity", e, n)?;
            let (coeff, got_coeff) = nodal(ctx, Some(&bc), "filmpressure transfer coefficient", e, n)?;
            let (extp, got_ext) = nodal(ctx, Some(&bc), "external filmpressure", e, n)?;
            if got_ext != got_coeff {
                return Err(fatal("Give neither or both keywords for Robin BC!"));
            }
            if !(open || got_ext || got_flux || got_velo) || etype.family() == 1 {
                continue;
            }
            if !nodes.iter().all(|&k| in_perm(k)) {
                continue;
            }
            let full = |pe: ElemId| ctx.mesh.element(pe).nodes.iter().all(|&k| in_perm(k));
            let parent = match parents {
                [Some(l), _] if full(l) => l,
                [_, Some(r)] if full(r) => r,
                _ => return Err(fatal("No proper parent element available for specified boundary")),
            };
            let (_, mat, _) = sections(ctx, parent)?;
            let mut gap = nodal(ctx, Some(&mat), "gap height", e, n)?.0;
            if let Some(m) = min_gap {
                gap.iter_mut().for_each(|g| *g = g.max(m));
            }
            let visc_n = nodal(ctx, Some(&mat), "viscosity", e, n)?.0;
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| fatal("no FE backend"))?;
            let mut local = LocalSystem::zeros(n);
            for qp in fe.quadrature.gauss_points(etype, None) {
                fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let b = &be.basis[..n];
                let s = qp.weight * be.det_j;
                let gp = gap.iter().zip(b).fold(0.0, |acc, (x, y)| acc + x * y);
                let pr = last_pres.iter().zip(b).fold(0.0, |acc, (x, y)| acc + x * y);
                let tot = reference + pr;
                let (visc, dl) = match visc_model {
                    Visc::Newtonian => (dot(b, &visc_n), 0.8488),
                    Visc::Rarefied(mfp0) => {
                        let (v, kn) = rarefied(dot(b, &visc_n), mfp0, reference, tot, gp);
                        (v, 0.8488 * (1.0 + 2.676 * kn.powf(0.659)))
                    }
                };
                let damp = gp / (12.0 * dl * visc);
                let density = match compress {
                    Compress::None => 1.0,
                    Compress::Weak(k) => (tot / k).exp(),
                    Compress::GasIsothermal => tot,
                    Compress::GasAdiabatic(r) => tot.powf(1.0 / r),
                };
                let (mut a, mut bb) = (0.0, 0.0);
                if got_ext {
                    a = -dot(b, &coeff);
                    bb = a * dot(b, &extp);
                }
                if got_flux {
                    bb -= dot(b, &flux);
                }
                if got_velo {
                    bb -= gp * dot(b, &vpres);
                }
                if open {
                    a -= damp * gp * density;
                }
                for pp in 0..n {
                    local.force[pp] += s * b[pp] * bb;
                    for q in 0..n {
                        local.stiff[pp * n + q] += s * b[q] * b[pp] * a;
                    }
                }
            }
            if transient {
                local.mass = Some(vec![0.0; n * n]);
                ctx.default_first_order_time(&mut local, e, dt)?;
                local.mass = None;
            }
            ctx.default_update(&local, e)?;
        }
        Ok(())
    }
}

impl Solver for ReynoldsPost {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(POST, "--------------------------------------------------", 5);
        ctx.info(POST, "Computing postprocessing fields from film pressure", 5);
        ctx.info(POST, "--------------------------------------------------", 5);
        let p = ctx.params().clone();
        if ctx.primary().is_none_or(|v| v.dofs != 1) {
            return Err(SolverError::fatal(POST, "Impossible number of dofs! (should be 1)"));
        }
        if ctx.simulation_params().string("coordinate system").is_some_and(|c| !elmer_core::normalize_key(c).starts_with("cartesian")) {
            return Err(SolverError::Unsupported("ReynoldsPostprocess: non-Cartesian coordinates".into()));
        }
        let pname = p.string("reynolds pressure variable name").map_or_else(|| "filmpressure".to_string(), elmer_core::normalize_key);
        let Some((pres_all, pperm)) = ctx.variable(&pname).map(|v| (v.values.to_vec(), v.perm.map(|q| q.to_vec()))) else {
            return Err(SolverError::fatal(POST, format!("Could not find primary variable: {pname}")));
        };
        if ctx.variable(&format!("{pname} gap sensitivity")).is_some() {
            return Err(SolverError::Unsupported("ReynoldsPostprocess: a gap sensitivity".into()));
        }
        if pperm.as_ref().is_some_and(|q| q.iter().all(|&s| s == NO_DOF)) {
            return Ok(());
        }
        let min_gap = p.real("min gap height");
        let ambient = p.real("ambient pressure").unwrap_or(0.0);
        let moment_about: Option<[f64; 3]> = match p.get("moment about") {
            Some(Value::Real { data, .. }) if data.len() > 1 => Some(std::array::from_fn(|i| data.get(i).copied().unwrap_or(0.0))),
            _ => {
                let m: [Option<f64>; 3] = std::array::from_fn(|i| p.real(&format!("moment about {}", i + 1)));
                m.iter().any(Option::is_some).then(|| std::array::from_fn(|i| m[i].unwrap_or(0.0)))
            }
        };
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut be = BasisEval::default();
        let (mut heat_slide, mut heat_pres) = (0.0, 0.0);
        let (mut pforce, mut vforce, mut moment) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        let own = ctx.variable_name().map(str::to_string).unwrap_or_default();
        ctx.info(POST, &format!("Primary variable name: {own}"), 4);
        for mode in 1..=4 {
            let name = match mode {
                1 => format!("{pname} force"),
                2 => format!("{pname} flux"),
                3 => format!("{pname} mean velocity"),
                _ => format!("{pname} heating"),
            };
            let Some(components) = ctx.variable(&name).map(|v| v.dofs) else { continue };
            for component in 1..=components {
                ctx.default_initialize()?;
                let mut mat_old = None;
                let (mut reference, mut visc_model, mut opposing) = (0.0, Visc::Newtonian, false);
                for &e in &elems {
                    let el = ctx.mesh.element(e);
                    let nodes = el.nodes.to_vec();
                    let etype = el.etype;
                    let n = nodes.len();
                    let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                    let pres = local_values(&pres_all, pperm.as_deref(), &nodes);
                    let (mat_id, mat, bf) = sections(ctx, e)?;
                    let (velo, use_velocity) = velocity(ctx, bf.as_ref(), Some(&mat), e, n)?;
                    let mut gap = nodal(ctx, Some(&mat), "gap height", e, n)?.0;
                    if let Some(m) = min_gap {
                        gap.iter_mut().for_each(|g| *g = g.max(m));
                    }
                    let visc_n = nodal(ctx, Some(&mat), "viscosity", e, n)?.0;
                    if mat_old != Some(mat_id) {
                        mat_old = Some(mat_id);
                        reference = mat.real("reference pressure").unwrap_or(0.0);
                        visc_model = models(&mat, true)?.1;
                        opposing = mat.logical("opposing wall") == Some(true);
                    }
                    let normal = match plane_normal(&xyz, etype.0) {
                        Some(nv) => nv,
                        None => {
                            return Err(SolverError::Unsupported(
                                "ReynoldsPostprocess: elements outside the xy plane or line elements".into(),
                            ))
                        }
                    };
                    let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| SolverError::fatal(POST, "no FE backend"))?;
                    let mut local = LocalSystem::zeros(n);
                    for qp in fe.quadrature.gauss_points(etype, None) {
                        fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                        let (b, g) = (&be.basis[..n], &be.dbasis_dx[..n]);
                        let s = qp.weight * be.det_j;
                        let radius = moment_about.map(|m| {
                            let x: [f64; 3] = std::array::from_fn(|a| (0..n).fold(0.0, |acc, k| acc + xyz[k][a] * b[k]));
                            [x[0] - m[0], x[1] - m[1], x[2] - m[2]]
                        });
                        let v: [f64; 3] = std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + b[k] * velo[k][i]));
                        let tangent: [f64; 3] = if use_velocity {
                            let vn = normal.iter().zip(&v).fold(0.0, |acc, (x, y)| acc + x * y);
                            std::array::from_fn(|i| v[i] - vn * normal[i])
                        } else {
                            v
                        };
                        let pr = dot(b, &pres);
                        let mut tot = pr + reference;
                        let gp = dot(b, &gap);
                        let grad: [f64; 3] = std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + g[k][i] * pres[k]));
                        let mut visc = dot(b, &visc_n);
                        if let Visc::Rarefied(mfp0) = visc_model {
                            visc = rarefied(visc, mfp0, reference, tot, gp).0;
                        }
                        tot -= ambient;
                        let c = component - 1;
                        let source = match mode {
                            1 => {
                                let mut spres = -tot * normal[c];
                                if opposing {
                                    spres = -spres;
                                }
                                spres -= gp * grad[c] / 2.0;
                                let sslide = visc * tangent[c] / gp;
                                let source = spres + sslide;
                                pforce[c] += s * spres;
                                vforce[c] += s * sslide;
                                if let Some(r) = radius {
                                    match component {
                                        1 => {
                                            moment[1] += r[2] * s * source;
                                            moment[2] -= r[1] * s * source;
                                        }
                                        2 => {
                                            moment[0] -= r[2] * s * source;
                                            moment[2] += r[0] * s * source;
                                        }
                                        _ => {
                                            moment[0] += r[1] * s * source;
                                            moment[1] -= r[0] * s * source;
                                        }
                                    }
                                }
                                source
                            }
                            2 => -(gp * gp * gp / (12.0 * visc)) * grad[c] + gp * tangent[c] / 2.0,
                            3 => -(gp * gp / (12.0 * visc)) * grad[c] + tangent[c] / 2.0,
                            _ => {
                                let spres = (gp * gp * gp / (12.0 * visc)) * (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
                                let sslide = (visc / gp) * (tangent[0] * tangent[0] + tangent[1] * tangent[1] + tangent[2] * tangent[2]);
                                heat_pres += s * spres;
                                heat_slide += s * sslide;
                                spres + sslide
                            }
                        };
                        for pp in 0..n {
                            for q in 0..n {
                                local.stiff[pp * n + q] += s * (b[q] * b[pp]);
                            }
                            local.force[pp] += s * b[pp] * source;
                        }
                    }
                    ctx.default_update(&local, e)?;
                }
                ctx.default_dirichlet()?;
                ctx.info(POST, "Dirichlet conditions done", 4);
                ctx.default_solve()?;
                if mode != 4 {
                    let sol = ctx.primary().expect("checked").values.to_vec();
                    let out = ctx.variable_mut(&name).expect("checked");
                    for (k, x) in sol.iter().enumerate() {
                        if let Some(o) = out.values.get_mut(component - 1 + k * components) {
                            *o = *x;
                        }
                    }
                }
            }
            let es = |x: f64| elmer_solver::elmer_saveio::fortran::fortran_es(x, 15, 4);
            let sim = &mut ctx.sim.model.simulation;
            if mode == 1 {
                for (i, f) in pforce.iter().enumerate() {
                    ctx.log.info(POST, &format!("{:<34}{}", format!("Pressure force {} (N):", i + 1), es(*f)), 5);
                }
                for (i, f) in pforce.iter().enumerate() {
                    sim.insert(&format!("res: Pressure force {}", i + 1), Value::real(*f));
                }
                for (i, f) in vforce.iter().enumerate() {
                    ctx.log.info(POST, &format!("{:<34}{}", format!("Sliding force {} (N):", i + 1), es(*f)), 5);
                    sim.insert(&format!("res: Sliding force {}", i + 1), Value::real(*f));
                }
                let tot: f64 = (0..3).map(|i| (pforce[i] + vforce[i]).powi(2)).fold(0.0, |a, x| a + x).sqrt();
                ctx.log.info(POST, &format!("{:<34}{}", "Reynolds force (N): ", es(tot)), 5);
                sim.insert("res: Reynolds force", Value::real(tot));
                if moment_about.is_some() {
                    for (i, m) in moment.iter().enumerate() {
                        ctx.log.info(POST, &format!("{:<34}{}", format!("Reynolds moment {} (Nm):", i + 1), es(*m)), 5);
                        sim.insert(&format!("res: Reynolds moment {}", i + 1), Value::real(*m));
                    }
                }
            } else if mode == 4 {
                let total = heat_pres + heat_slide;
                ctx.log.info(POST, &format!("{:<34}{}", "Pressure heating (W): ", es(heat_pres)), 5);
                ctx.log.info(POST, &format!("{:<34}{}", "Sliding heating (W): ", es(heat_slide)), 5);
                ctx.log.info(POST, &format!("{:<34}{}", "Reynolds heating (W): ", es(total)), 5);
                sim.insert("res: Reynolds heating", Value::real(total));
            }
        }
        Ok(())
    }
}
