//! `FreeSurfaceSolver` (fem/src/modules/FreeSurfaceSolver.F90): the
//! kinematic free-surface equation `dz/dt + u dz/dx (+ v dz/dy) = w + a`
//! on the surface elements (a `Body Id` boundary in 2D), with the velocity
//! of `Flow Solution Name` (`Convection = Computed`), constant
//! (`Convection Velocity i`) or none, the accumulation `<var>
//! Accumulation` (times `sqrt(1 + |grad z|^2)`) or the flux `<var>
//! Accumulation Flux i`, and residual-squared (SUPG) stabilization. The
//! element is projected to the horizontal plane (the vertical coordinate
//! set to zero) before the integration.
//!
//! Ported: the `stabilized` method, the Eulerian and the horizontal ALE
//! formulation (`ALE Formulation`, the mesh velocity subtracted),
//! `Velocity Implicity`, the nonlinear iterations, transient runs through
//! `Default1stOrderTime` (the predictor-corrector Adams schemes
//! included). Refused: `bubbles`, `Apply Dirichlet` (the limiters),
//! relaxation and displacement limits, rotated free surfaces, `Use linear
//! elements` on quadratic elements.

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "FreeSurfaceSolver";

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(FreeSurfaceSolver))
}

pub struct FreeSurfaceSolver;

impl Solver for FreeSurfaceSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        let name = ctx.variable_name().unwrap_or_default().to_string();
        let caller = format!("{CALLER} ({name})");
        ctx.info(&caller, "Solving for free surface", 5);
        let (perm, nvals) = {
            let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "Variable values not associated"))?;
            (v.perm.map(|p| p.to_vec()), v.values.len())
        };
        let perm = perm.ok_or_else(|| SolverError::fatal(CALLER, "no permutation"))?;
        if !perm.iter().any(|&s| s != NO_DOF) {
            ctx.warn(&caller, "A serial run without a free surface, but the solver switched in - weird!");
            return Ok(());
        }
        let dim = ctx.mesh.dim();
        if !(2..=3).contains(&dim) {
            return Err(SolverError::fatal(CALLER, format!("It is not possible to compute free-surface problems in DIM={dim} dimensions. Aborting")));
        }
        let cv = ctx.creal(&p, "velocity implicity")?.unwrap_or(1.0);
        if p.real("linear system convergence tolerance").is_none() {
            return Err(SolverError::fatal(CALLER, "No > Linear System Convergence Tolerance < found"));
        }
        let nonlinear_tol = p.real("nonlinear system convergence tolerance").unwrap_or(0.0);
        let nonlinear_iter = p.integer("nonlinear system max iterations").unwrap_or(1).max(1);
        for (key, what) in [
            ("maximum displacement", "Maximum Displacement"),
            ("compute local maximum displacement", "Compute Local Maximum Displacement"),
        ] {
            if p.get(key).is_some() {
                return Err(unsupported(what));
            }
        }
        if ctx.creal(&p, "relaxation factor")?.is_some_and(|r| r < 1.0) {
            return Err(unsupported("Relaxation Factor"));
        }
        if p.logical("apply dirichlet") == Some(true) {
            return Err(unsupported("Apply Dirichlet (the limiters)"));
        }
        if p.logical("use linear elements") == Some(true) {
            return Err(unsupported("Use linear elements"));
        }
        let ale = p.logical("ale formulation") == Some(true);
        let bubbles = match p.string("stabilization method").map(elmer_core::normalize_key).as_deref() {
            Some("bubbles") => return Err(unsupported("Stabilization Method = bubbles")),
            _ => false,
        };
        let _ = bubbles;
        if ctx.variable(&format!("{name} orientation")).is_some() {
            return Err(unsupported("rotated free surfaces (Orientation)"));
        }
        ctx.info(&caller, &format!("Size of problem: {nvals}"), 8);

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut be = BasisEval::default();
        let mut prev_norm;
        for iter in 1..=nonlinear_iter {
            ctx.info(&caller, &format!("Non-linear Iteration {iter} out of max. {nonlinear_iter}"), 4);
            ctx.default_initialize()?;
            for &e in &elems {
                let el = ctx.mesh.element(e);
                if el.etype.family() == 1 {
                    continue;
                }
                let nodes = el.nodes.to_vec();
                let etype = el.etype;
                let n = nodes.len();
                // The element in the horizontal plane.
                let xyz: Vec<[f64; 3]> = nodes
                    .iter()
                    .map(|&k| {
                        let c = ctx.mesh.coord(k);
                        if dim == 2 {
                            [c[0], 0.0, 0.0]
                        } else {
                            [c[0], c[1], 0.0]
                        }
                    })
                    .collect();
                let bp = ctx.body_params(e).ok_or_else(|| SolverError::fatal(CALLER, "a surface element without a body"))?;
                let (mat, bf, eq) = (bp.material.cloned(), bp.body_force.cloned(), bp.equation.cloned());
                let flag = eq
                    .as_ref()
                    .and_then(|q| q.string("convection"))
                    .map(elmer_core::normalize_key)
                    .ok_or_else(|| SolverError::fatal(CALLER, "No string for keyword > Convection < found in Equation"))?;
                let mut velo = vec![[0.0; 3]; n];
                match flag.as_str() {
                    "constant" => {
                        for i in 0..3 {
                            let key = format!("convection velocity {}", i + 1);
                            let v = match mat.as_ref().and_then(|m| m.get(&key)).is_some() {
                                true => ctx.real_at_nodes(mat.as_ref().expect("checked"), &key, e)?,
                                false => match &eq {
                                    Some(q) => ctx.real_at_nodes(q, &key, e)?,
                                    None => None,
                                },
                            };
                            if let Some(v) = v {
                                for (o, x) in velo.iter_mut().zip(v) {
                                    o[i] = x;
                                }
                            }
                        }
                    }
                    "computed" => {
                        let fname = eq
                            .as_ref()
                            .and_then(|q| q.string("flow solution name"))
                            .map_or_else(|| "flow solution".to_string(), elmer_core::normalize_key);
                        let f = ctx.variable(&fname).ok_or_else(|| {
                            SolverError::fatal(CALLER, format!("Convection flag set to > computed <, but no variable >{fname}< found"))
                        })?;
                        let ns = f.dofs;
                        if !(dim == 2 && ns == 3 || dim == 3 && ns == 4) {
                            return Err(SolverError::fatal(CALLER, format!("DIM={dim} NSDOFs={ns} does not combine. Aborting")));
                        }
                        let fperm = f.perm.map(|p| p.to_vec());
                        let prev: Option<Vec<f64>> = ctx
                            .sim
                            .vars
                            .get(&fname)
                            .and_then(|v| v.prev_values.first().cloned());
                        let fv = ctx.variable(&fname).expect("checked").values.to_vec();
                        for (o, &k) in velo.iter_mut().zip(&nodes) {
                            let s = fperm.as_ref().map_or(k, |p| p[k as usize]);
                            if s == NO_DOF {
                                return Err(SolverError::fatal(CALLER, "a surface node outside the flow solution"));
                            }
                            let j = s as usize * ns;
                            for i in 0..dim {
                                o[i] = if transient && (cv - 1.0).abs() > 0.001 {
                                    let pv = prev.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no previous flow solution"))?;
                                    cv * fv[j + i] + (1.0 - cv) * pv[j + i]
                                } else {
                                    fv[j + i]
                                };
                            }
                        }
                    }
                    _ => {}
                }
                // GetVectorLocalSolution(MeshVelocity, 'Mesh Velocity').
                let mut mesh_velo = vec![[0.0; 3]; n];
                if ale {
                    if let Some(mv) = ctx.variable("mesh velocity") {
                        let mperm = mv.perm.map(|p| p.to_vec());
                        for (o, &k) in mesh_velo.iter_mut().zip(&nodes) {
                            let s = mperm.as_ref().map_or(k, |p| p[k as usize]);
                            if s != NO_DOF {
                                for i in 0..mv.dofs.min(3) {
                                    o[i] = mv.values[s as usize * mv.dofs + i];
                                }
                            }
                        }
                    }
                }
                // The accumulation (a normal flux) or a flux vector.
                let mut source = vec![0.0; n];
                let mut flux = vec![[0.0; 3]; n];
                let mut normal_flux = false;
                let substantial = bf.is_none();
                if let Some(b) = &bf {
                    match ctx.real_at_nodes(b, &format!("{name} accumulation"), e)? {
                        Some(v) => {
                            source = v;
                            normal_flux = true;
                        }
                        None => {
                            for i in 0..dim {
                                if let Some(v) = ctx.real_at_nodes(b, &format!("{name} accumulation flux {}", i + 1), e)? {
                                    for (o, x) in flux.iter_mut().zip(v) {
                                        o[i] = x;
                                    }
                                }
                            }
                        }
                    }
                }
                // ElemFreeSurf: the previous step's surface.
                let old: Vec<f64> = if transient {
                    let v = ctx.sim.vars.get(&name).ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
                    let h = v.prev_values.first().ok_or_else(|| SolverError::fatal(CALLER, "no history"))?;
                    nodes.iter().map(|&k| h[perm[k as usize] as usize]).collect()
                } else {
                    vec![0.0; n]
                };

                // LocalMatrix.
                let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no FE backend"))?;
                let hk = elmer_fe::element_diameter(etype, &xyz, false);
                let mut local = LocalSystem::zeros(n);
                let mut mass = vec![0.0; n * n];
                for qp in fe.quadrature.gauss_points(etype, None) {
                    fe.basis.eval(etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let (b, g) = (&be.basis[..n], &be.dbasis_dx[..n]);
                    let s = qp.weight * be.det_j;
                    let mut grad_fs = [0.0; 3];
                    for (i, gi) in grad_fs.iter_mut().enumerate().take(dim - 1) {
                        *gi = (0..n).fold(0.0, |acc, k| acc + g[k][i] * old[k]);
                    }
                    grad_fs[dim - 1] = 1.0;
                    let mut vg = [0.0; 3];
                    for (i, v) in vg.iter_mut().enumerate().take(dim) {
                        *v = if ale && i < dim - 1 {
                            (0..n).fold(0.0, |acc, k| acc + b[k] * (velo[k][i] - mesh_velo[k][i]))
                        } else {
                            (0..n).fold(0.0, |acc, k| acc + b[k] * velo[k][i])
                        };
                    }
                    let norm_grad = if dim == 3 {
                        (1.0 + grad_fs[0].powi(2) + grad_fs[1].powi(2)).sqrt()
                    } else {
                        (1.0 + grad_fs[0].powi(2)).sqrt()
                    };
                    let unorm = vg[..dim - 1].iter().map(|v| v * v).sum::<f64>().sqrt();
                    let tau = if unorm != 0.0 { hk / (2.0 * unorm) } else { 0.0 };
                    let su: Vec<f64> = (0..n).map(|p| (0..dim - 1).fold(0.0, |acc, i| acc + vg[i] * g[p][i])).collect();
                    let sw = su.clone();
                    for pp in 0..n {
                        for q in 0..n {
                            let a = &mut local.stiff[pp * n + q];
                            for i in 0..dim - 1 {
                                *a += s * vg[i] * g[q][i] * b[pp];
                            }
                            *a += s * tau * su[q] * sw[pp];
                        }
                    }
                    if transient {
                        for pp in 0..n {
                            for q in 0..n {
                                mass[pp * n + q] += s * b[q] * (b[pp] + tau * sw[pp]);
                            }
                        }
                    }
                    let src = if substantial {
                        0.0
                    } else if normal_flux {
                        norm_grad * (0..n).fold(0.0, |acc, k| acc + source[k] * b[k])
                    } else {
                        (0..dim).fold(0.0, |acc, i| acc + (0..n).fold(0.0, |a2, k| a2 + b[k] * flux[k][i]) * grad_fs[i])
                    };
                    for pp in 0..n {
                        local.force[pp] += (vg[dim - 1] + src) * (b[pp] + tau * sw[pp]) * s;
                    }
                }
                if transient {
                    local.mass = Some(mass);
                    ctx.default_first_order_time(&mut local, e, dt)?;
                }
                local.mass = None;
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;
            ctx.default_finish_boundary_assembly()?;
            ctx.default_dirichlet()?;
            prev_norm = ctx.primary().map_or(0.0, |_| ctx.convergence().map_or(0.0, |c| c.norm));
            let norm = ctx.default_solve()?;
            let change = if prev_norm + norm != 0.0 { 2.0 * (prev_norm - norm).abs() / (prev_norm + norm) } else { 0.0 };
            ctx.info(&caller, &format!(" Result Norm   : {norm:>25.17E}"), 4);
            ctx.info(&caller, &format!(" Relative Change : {change:>25.17E}"), 4);
            if change < nonlinear_tol {
                ctx.info(&caller, &format!("Converged after{iter} iterations"), 4);
                break;
            }
        }
        Ok(())
    }
}
