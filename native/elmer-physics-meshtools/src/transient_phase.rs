//! `TransientPhaseChange` (fem/src/modules/TransientPhaseChange.F90): the
//! transient motion of a solid-liquid interface from the jump of the heat
//! flux across it (Stefan condition), on the interface's own elements.
//!
//! Each visit solves for the interface velocity `<var> Velo` (an exported
//! field), `rho L (n . e_N) v = -[k grad T] . n` integrated over the
//! interface elements, with `Velocity Smoothing Factor` smoothing the
//! velocity along the interface away from its end nodes (`Phase Change Side`
//! BCs, and the triple and axis points). With `Use Nodal Loads` the
//! right-hand side is the temperature solver's `<T> Loads` instead
//! (negated by `Nodal Loads Negative`, averaged with the previous step's by
//! `Loads Relaxation Factor`). The velocity is relaxed with the previous
//! one (`Velocity Relaxation Factor`), and the surface moves by the simple
//! update `Surface = PrevSurface + SpeedUp dt (Velo + UPull)`, the pull
//! velocity being the material's `Convection Velocity <normal direction>`.
//!
//! The velocity solve's norm is the solver's (`SolveSystem` into the
//! velocity, with the solver's variable keeping the norms).
//!
//! Refused: `Pull Rate Control`, `Triple Point Fixed`, a `Normals` /
//! `Normal Variable` field, `Passive Steps`, the PDE surface update
//! (`Surface Smoothing Factor` or a time order above 1), and 3D.

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Model, Solver, SolverCtx, SolverError};

const CALLER: &str = "TransientPhaseChange";

/// `TransientPhaseChange_Init`: the `-nooutput <var> Velo` export and a
/// first-order time derivative.
pub fn transient_defaults(params: &mut ParamMap) {
    let Some(var) = params.string("variable").map(|v| v.trim().to_string()) else { return };
    let name = var.split('[').next().unwrap_or(&var).trim().to_string();
    let k = (1..).find(|k| params.get(&format!("exported variable {k}")).is_none()).expect("a free index");
    params.insert(&format!("Exported Variable {k}"), Value::Str(format!("-nooutput {name} Velo")));
    if params.get("time derivative order").is_none() {
        params.insert("Time Derivative Order", Value::Integer(vec![1]));
    }
}

pub fn factory(p: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    for key in ["pull rate control", "triple point fixed"] {
        if p.logical(key) == Some(true) {
            return Err(SolverError::Unsupported(format!("TransientPhaseChange: {key}")));
        }
    }
    for key in ["passive steps", "normal variable"] {
        if p.get(key).is_some() {
            return Err(SolverError::Unsupported(format!("TransientPhaseChange: {key}")));
        }
    }
    if p.real("surface smoothing factor").is_some_and(|f| f.abs() >= f64::EPSILON) {
        return Err(SolverError::Unsupported("TransientPhaseChange: the surface update by a PDE (Surface Smoothing Factor)".into()));
    }
    Ok(Box::new(TransientPhaseChange::default()))
}

#[derive(Default)]
pub struct TransientPhaseChange {
    ready: bool,
    /// 0-based normal coordinate.
    normal_dir: usize,
    boundary_node: Vec<bool>,
    upull: Option<f64>,
    prev_loads: Vec<f64>,
    current_loads: Vec<f64>,
}

fn material(model: &Model, body: usize) -> Option<&ParamMap> {
    let m = model.bodies.get(body.checked_sub(1)?)?.integer("material")?;
    model.materials.get((m as usize).checked_sub(1)?)
}

/// `NodeU/V/W` of a parent element's nodes.
fn node_uvw(code: u16, k: usize) -> Option<[f64; 3]> {
    match code {
        303 => [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]].get(k).copied(),
        404 => [[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]].get(k).copied(),
        _ => None,
    }
}

impl Solver for TransientPhaseChange {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "--------------------------------------------", 4);
        ctx.info(CALLER, "Using transient algorithm for surface update", 4);
        ctx.info(CALLER, "--------------------------------------------", 4);
        if !transient {
            return Err(SolverError::fatal(CALLER, "This only makes sense in a transient setting"));
        }
        if ctx.mesh.dim() != 2 {
            return Err(SolverError::Unsupported("TransientPhaseChange: 3D".into()));
        }
        let method = ctx.simulation_params().string("timestepping method").map(elmer_core::normalize_key).unwrap_or_default();
        let bdf1 = ctx.simulation_params().integer("bdf order") == Some(1);
        if !(method == "implicit euler" || (method == "bdf" && bdf1)) {
            return Err(SolverError::Unsupported("TransientPhaseChange: the surface update by a PDE (time order above 1)".into()));
        }
        let p = ctx.params().clone();
        let dim = 2;
        let (surf_name, surf_perm) = {
            let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "Surface field needed for Phase Change"))?;
            (v.name.to_string(), v.perm.map(|q| q.to_vec()).ok_or_else(|| SolverError::fatal(CALLER, "Surface field needed for Phase Change"))?)
        };
        if surf_perm.iter().all(|&s| s == NO_DOF) {
            return Err(SolverError::fatal(CALLER, "Surface field needed for Phase Change"));
        }
        let nslots = ctx.primary().expect("checked").values.len();
        let temp_name = p.string("phase change variable").map_or_else(|| "temperature".to_string(), elmer_core::normalize_key);
        let (temp, temp_perm) = {
            let v = ctx.variable(&temp_name).ok_or_else(|| SolverError::fatal(CALLER, "Temperature field needed for Phase Change"))?;
            (v.values.to_vec(), v.perm.map(|q| q.to_vec()))
        };
        let tslot = |k: u32| temp_perm.as_ref().map_or(Some(k as usize), |q| (q[k as usize] != NO_DOF).then(|| q[k as usize] as usize));
        let coupled_iter = ctx.scalar("coupled iter").unwrap_or(1.0).round() as i32;
        let timestep = ctx.scalar("timestep").unwrap_or(1.0).round() as i32;
        let velo_name = format!("{surf_name} velo");
        if ctx.variable(&velo_name).is_none() {
            return Err(SolverError::fatal(CALLER, "Surface Velo field should exist!"));
        }

        if !self.ready {
            ctx.info(CALLER, "Doing some first time initializations", 7);
            let (mut cmax, mut cmin) = ([-f64::MAX; 3], [f64::MAX; 3]);
            let (mut imax, mut imin) = ([0usize; 3], [0usize; 3]);
            for k in 0..ctx.mesh.n_nodes() {
                if surf_perm[k] == NO_DOF {
                    continue;
                }
                let c = ctx.mesh.coord(k as u32);
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
            self.normal_dir = match p.integer("normal direction") {
                Some(d) if (1..=dim as i32).contains(&d) => d as usize - 1,
                Some(d) => return Err(SolverError::fatal(CALLER, format!("Normal Direction {d}"))),
                None => {
                    let mut j = 0;
                    for i in 0..dim {
                        if cmax[i] - cmin[i] < cmax[j] - cmin[j] {
                            j = i;
                        }
                    }
                    ctx.info(CALLER, &format!("Normal coordinate set to: {}", j + 1), 7);
                    j
                }
            };
            self.boundary_node = vec![false; nslots];
            let tangent = 1 - self.normal_dir;
            for k in [imax[tangent], imin[tangent]] {
                self.boundary_node[surf_perm[k] as usize] = true;
            }
            let model = &ctx.sim.model;
            let nb = ctx.mesh.n_bulk();
            for t in 0..ctx.mesh.n_boundary() {
                let el = ctx.mesh.element(ElemId((nb + t) as u32));
                let side = model.bc_of_tag(el.tag).is_some_and(|b| model.boundary_conditions[b - 1].logical("phase change side") == Some(true));
                if !side {
                    continue;
                }
                for &k in el.nodes {
                    if surf_perm[k as usize] != NO_DOF {
                        self.boundary_node[surf_perm[k as usize] as usize] = true;
                    }
                }
            }
            let n = self.boundary_node.iter().filter(|&&b| b).count();
            ctx.info(CALLER, &format!("Number of boundary nodes: {n}"), 7);
            if ctx.variable("normals").is_some() {
                return Err(SolverError::Unsupported("TransientPhaseChange: a Normals field".into()));
            }
            self.ready = true;
        }
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        if self.upull.is_none() {
            let first = *elems.first().ok_or_else(|| SolverError::fatal(CALLER, "no active elements"))?;
            let model = &ctx.sim.model;
            let body = model.boundary_body(ctx.mesh.element(first).tag).or_else(|| model.body_of_tag(ctx.mesh.element(first).tag));
            let mat = body.and_then(|b| material(model, b));
            self.upull = Some(mat.and_then(|m| m.real(&format!("convection velocity {}", self.normal_dir + 1))).unwrap_or(0.0));
        }
        let upull = self.upull.expect("set");

        let use_loads = p.logical("use nodal loads") == Some(true);
        let mut use_first_loads = false;
        if use_loads {
            let loads_name = format!("{temp_name} loads");
            let (lv, lp) = {
                let v = ctx.variable(&loads_name).ok_or_else(|| SolverError::fatal(CALLER, "Loads are requested to be used but missing!"))?;
                (v.values.to_vec(), v.perm.map(|q| q.to_vec()))
            };
            ctx.info(CALLER, "Using nodal loads to determine phase change", 6);
            let sign = if p.logical("nodal loads negative") == Some(true) {
                ctx.info(CALLER, "Assigning negative sign to nodal loads!", 6);
                -1.0
            } else {
                1.0
            };
            if self.prev_loads.is_empty() {
                self.prev_loads = vec![0.0; nslots];
                self.current_loads = vec![0.0; nslots];
            }
            use_first_loads = p.logical("use first loads") == Some(true);
            if !use_first_loads && timestep > 1 && coupled_iter == 1 {
                self.prev_loads.clone_from(&self.current_loads);
            }
            for t in 0..ctx.mesh.n_nodes() {
                let i = surf_perm[t];
                if i == NO_DOF {
                    continue;
                }
                let j = lp.as_ref().map_or(t as u32, |q| q[t]);
                if j == NO_DOF {
                    continue;
                }
                self.current_loads[i as usize] = sign * lv[j as usize];
            }
            if use_first_loads && coupled_iter == 1 {
                self.prev_loads.clone_from(&self.current_loads);
            }
            let es = |x: f64| elmer_solver::elmer_saveio::fortran::fortran_es(x, 12, 5);
            let mn = self.current_loads.iter().copied().fold(f64::MAX, f64::min);
            let mx = self.current_loads.iter().copied().fold(-f64::MAX, f64::max);
            ctx.info(CALLER, &format!("Minimum nodal load at interface: {}", es(mn)), 4);
            ctx.info(CALLER, &format!("Maximum nodal load at interface: {}", es(mx)), 4);
        }
        let velo_relax = p.real("velocity relaxation factor");
        let prev_velo: Vec<f64> = ctx.variable(&velo_name).expect("checked").values.to_vec();
        let loads_relax = p.real("loads relaxation factor");
        let loads_order = match loads_relax {
            Some(_) if !use_first_loads || coupled_iter > 1 => (timestep - 1).min(2),
            _ => 0,
        };
        let stab = p.real("velocity smoothing factor").unwrap_or(0.0);
        let iterations = p.integer("nonlinear system max iterations").unwrap_or(1).max(1);
        let cartesian = ctx
            .simulation_params()
            .string("coordinate system")
            .is_none_or(|c| elmer_core::normalize_key(c).starts_with("cartesian"));
        let nd = self.normal_dir;
        let mut be = BasisEval::default();
        let mut pb = BasisEval::default();
        for _ in 0..iterations {
            ctx.default_initialize()?;
            if use_loads {
                let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no system"))?;
                for t in 0..surf_perm.len() {
                    let i = surf_perm[t];
                    if i == NO_DOF {
                        continue;
                    }
                    let i = i as usize;
                    sys.rhs[i] = match (loads_order, loads_relax) {
                        (1 | 2, Some(r)) => r * self.current_loads[i] + (1.0 - r) * self.prev_loads[i],
                        _ => self.current_loads[i],
                    };
                }
            }
            for &e in &elems {
                let el = ctx.mesh.element(e);
                if el.etype.0 != 202 {
                    return Err(SolverError::Unsupported(format!("TransientPhaseChange: interface elements of type {}", el.etype.0)));
                }
                let nodes = el.nodes.to_vec();
                let parents = el.parents;
                let tag = el.tag;
                let n = nodes.len();
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let model = &ctx.sim.model;
                let body = model.boundary_body(tag).ok_or_else(|| SolverError::fatal(CALLER, "an interface element without a body"))?;
                let mat = material(model, body).cloned().ok_or_else(|| SolverError::fatal(CALLER, "no material"))?;
                let out_body = model
                    .bc_of_tag(tag)
                    .and_then(|b| model.boundary_conditions[b - 1].integer("normal target body"))
                    .filter(|&b| b > 0)
                    .map(|b| b as usize);
                let latent = ctx.real_at_nodes(&mat, "latent heat", e)?.ok_or_else(|| SolverError::fatal(CALLER, "Latent Heat not found"))?;
                let density = ctx.real_at_nodes(&mat, "density", e)?.ok_or_else(|| SolverError::fatal(CALLER, "Density not found"))?;
                // The parents: body, material, nodes, coordinates.
                let mut par = Vec::new();
                for pe in parents.iter().flatten() {
                    let pel = ctx.mesh.element(*pe);
                    let pbody = ctx.sim.model.body_of_tag(pel.tag).unwrap_or(0);
                    par.push((*pe, pbody, pel.etype, pel.nodes.to_vec()));
                }
                let mut local = LocalSystem::zeros(n);
                let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no FE backend"))?;
                let rule = fe.quadrature.gauss_points(el.etype, None);
                for qp in rule {
                    fe.basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                    let b: Vec<f64> = be.basis[..n].to_vec();
                    let g: Vec<[f64; 3]> = be.dbasis_dx[..n].to_vec();
                    let mut s = qp.weight * be.det_j;
                    let xcoord = (0..n).fold(0.0, |acc, j| acc + xyz[j][0] * b[j]);
                    if use_loads && !cartesian {
                        s *= xcoord;
                    }
                    let normal = line_normal(&xyz, &b, &par, out_body, ctx.mesh);
                    let local_heat = (0..n).fold(0.0, |acc, j| acc + b[j] * latent[j]);
                    let local_dens = (0..n).fold(0.0, |acc, j| acc + b[j] * density[j]);
                    let stab_coeff = stab * local_heat * local_dens;
                    let mut flux = 0.0;
                    if !use_loads {
                        let mut tgrad = [[0.0; 3]; 2];
                        let mut l = 0usize;
                        for (_, pbody, ptype, pnodes) in &par {
                            let pm = material(&ctx.sim.model, *pbody).cloned().ok_or_else(|| SolverError::fatal(CALLER, "no parent material"))?;
                            if pm.logical("solid") == Some(true) {
                                if l == 2 {
                                    return Err(SolverError::fatal(CALLER, "Both materials cannot be solid!"));
                                }
                                l = 2;
                            } else {
                                if l == 1 {
                                    return Err(SolverError::fatal(CALLER, "Both materials cannot be liquid!"));
                                }
                                l = 1;
                            }
                            let pe = par.iter().find(|q| q.3 == *pnodes).map(|q| q.0).expect("listed");
                            let cond = ctx.real_at_nodes(&pm, "heat conductivity", pe)?.ok_or_else(|| SolverError::fatal(CALLER, "Heat Conductivity not found"))?;
                            let ptemp: Vec<f64> = pnodes
                                .iter()
                                .map(|&k| tslot(k).map(|s| temp[s]).ok_or_else(|| SolverError::fatal(CALLER, "a parent node without temperature")))
                                .collect::<Result<_, _>>()?;
                            let mut ref_xyz = [[0.0; 3]; 2];
                            for j in 0..n {
                                let k = pnodes.iter().position(|&q| q == nodes[j]).ok_or_else(|| SolverError::fatal(CALLER, "a node not in its parent"))?;
                                ref_xyz[j] = node_uvw(ptype.0, k).ok_or_else(|| SolverError::Unsupported(format!("TransientPhaseChange: parent elements of type {}", ptype.0)))?;
                            }
                            let puvw: [f64; 3] = std::array::from_fn(|a| (0..n).fold(0.0, |acc, j| acc + b[j] * ref_xyz[j][a]));
                            let pxyz: Vec<[f64; 3]> = pnodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                            fe.basis.eval(*ptype, &pxyz, puvw, &mut pb).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                            let pn = pnodes.len();
                            let kc = (0..pn).fold(0.0, |acc, j| acc + cond[j] * pb.basis[j]);
                            for jd in 0..dim {
                                tgrad[l - 1][jd] = kc * (0..pn).fold(0.0, |acc, j| acc + pb.dbasis_dx[j][jd] * ptemp[j]);
                            }
                        }
                        flux = (0..3).fold(0.0, |acc, a| acc + (tgrad[0][a] - tgrad[1][a]) * normal[a]);
                    }
                    for pp in 0..n {
                        for q in 0..n {
                            local.stiff[pp * n + q] += s * b[pp] * b[q] * normal[nd] * local_dens * local_heat;
                            if self.boundary_node[surf_perm[nodes[pp] as usize] as usize] {
                                continue;
                            }
                            for i in 0..dim {
                                if i == nd {
                                    continue;
                                }
                                local.stiff[pp * n + q] += s * stab_coeff * g[q][i] * g[pp][i];
                            }
                        }
                        if !use_loads {
                            local.force[pp] -= s * b[pp] * flux;
                        }
                    }
                }
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;
            ctx.default_finish_boundary_assembly()?;
            ctx.default_dirichlet()?;
            // SolveSystem into the velocity: the solver's variable keeps the norms.
            let surface: Vec<f64> = ctx.primary().expect("checked").values.to_vec();
            let velo: Vec<f64> = ctx.variable(&velo_name).expect("checked").values.to_vec();
            ctx.primary_mut().expect("checked").values.copy_from_slice(&velo);
            let solved = ctx.default_solve();
            let new_velo: Vec<f64> = ctx.primary().expect("checked").values.to_vec();
            ctx.primary_mut().expect("checked").values.copy_from_slice(&surface);
            solved?;
            ctx.variable_mut(&velo_name).expect("checked").values.copy_from_slice(&new_velo);
            if ctx.convergence().is_some_and(|c| c.nonlin_converged == 1) {
                break;
            }
        }
        if let Some(r) = velo_relax {
            let v = ctx.variable_mut(&velo_name).expect("checked");
            for (x, pv) in v.values.iter_mut().zip(&prev_velo) {
                *x = r * *x + (1.0 - r) * pv;
            }
        }
        let speedup = p.real("transient speedup").unwrap_or(1.0);
        ctx.info(CALLER, "Updating surface by simple update", 7);
        let velo: Vec<f64> = ctx.variable(&velo_name).expect("checked").values.to_vec();
        let prev: Vec<f64> = {
            let v = ctx.sim.vars.get(&surf_name).expect("the surface");
            v.prev_values.first().cloned().ok_or_else(|| SolverError::fatal(CALLER, "no previous surface"))?
        };
        let v = ctx.primary_mut().expect("checked");
        for ((s, pv), w) in v.values.iter_mut().zip(&prev).zip(&velo) {
            *s = pv + speedup * dt * (w + upull);
        }
        Ok(())
    }
}

/// `NormalVector(Element, Nodes, u, v, .TRUE.)` of a 202 line at the basis
/// `b`: the line's normal, turned by `CheckNormalDirection` (`Normal Target
/// Body`, else the higher body id) to point away from the chosen parent.
fn line_normal(
    x: &[[f64; 3]],
    b: &[f64],
    par: &[(ElemId, usize, elmer_core::ElementType, Vec<u32>)],
    out_body: Option<usize>,
    mesh: &dyn elmer_core::MeshView,
) -> [f64; 3] {
    let (dxdu, dydu) = (-0.5 * x[0][0] + 0.5 * x[1][0], -0.5 * x[0][1] + 0.5 * x[1][1]);
    let det = 1.0 / (dxdu * dxdu + dydu * dydu).sqrt();
    let mut normal = [-dydu * det, dxdu * det, 0.0];
    let pick = match par {
        [l, r] => match out_body {
            Some(k) if l.1 == k => Some(r),
            Some(_) => Some(l),
            None if l.1 > r.1 => Some(l),
            None if l.1 < r.1 => Some(r),
            None => None,
        },
        [l] => Some(l),
        _ => None,
    };
    if let Some((_, _, _, pn)) = pick {
        let corners = pn.len().min(4);
        let c: [f64; 3] = std::array::from_fn(|a| pn[..corners].iter().map(|&k| mesh.coord(k)[a]).sum::<f64>() / corners as f64);
        let point: [f64; 3] = std::array::from_fn(|a| b[0] * x[0][a] + b[1] * x[1][a]);
        let dot: f64 = (0..3).map(|a| (c[a] - point[a]) * normal[a]).sum();
        if dot > 0.0 {
            normal = normal.map(|c| -c);
        }
    }
    normal
}
