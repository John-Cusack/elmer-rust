//! `StreamSolver` (fem/src/modules/StreamSolver.F90): the stream function
//! of a 2D velocity field, `-lap psi = d v/dx - d u/dy` with the natural
//! boundary term `(-v, u) . n`, optionally the Stokes stream function of
//! an axisymmetric flow (`Stokes Stream Function`, default in axisymmetric
//! coordinates). `Stream Function First Node` fixes one node to zero; after
//! the solve the values are shifted to a zero minimum (`Stream Function
//! Shifting`, default on) and scaled to a unit range (`Stream Function
//! Scaling`). The solver's norm is that of the solve, before the shift.

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "StreamSolver";
/// `AEPS` (Types.F90).
const AEPS: f64 = 10.0 * f64::EPSILON;

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(StreamSolver))
}

pub struct StreamSolver;

/// The (shifted) velocity components at the element's nodes, as LOAD
/// holds them: bulk `(v, -u)`, boundary `(-v, u)`.
fn nodal_load(
    nodes: &[u32],
    perm: Option<&[u32]>,
    values: &[f64],
    dofs: usize,
    sign: f64,
) -> Result<Vec<[f64; 2]>, SolverError> {
    nodes
        .iter()
        .map(|&k| {
            let s = perm.map_or(k, |p| p[k as usize]);
            if s == NO_DOF {
                return Err(SolverError::fatal(CALLER, format!("node {} has no velocity", k + 1)));
            }
            let b = s as usize * dofs;
            Ok([sign * values[b + 1], -sign * values[b]])
        })
        .collect()
}

impl Solver for StreamSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(SolverError::Unsupported("StreamSolver: transient runs".into()));
        }
        ctx.info(CALLER, "STREAMSOLVER", 4);
        let p = ctx.params().clone();
        let flow_name = p.string("stream function velocity variable").map_or_else(|| "flow solution".to_string(), elmer_core::normalize_key);
        let Some((fdofs, fperm, fvals)) = ctx.variable(&flow_name).map(|f| (f.dofs, f.perm.map(|p| p.to_vec()), f.values.to_vec())) else {
            ctx.warn(CALLER, "No variable for velocity associated.");
            ctx.warn(CALLER, "Quitting execution of StreamSolver.");
            return Ok(());
        };
        if fdofs < 2 {
            return Err(SolverError::fatal(CALLER, "the velocity variable has fewer than 2 components"));
        }
        let nn = ctx.mesh.n_nodes() as i64;
        let first = p.integer("stream function first node").map(|f| i64::from(f).clamp(1, nn));
        let shifting = p.logical("stream function shifting").unwrap_or(true);
        let scaling = p.logical("stream function scaling").unwrap_or(false);
        if p.logical("dirichlet point") == Some(true) {
            return Err(SolverError::Unsupported("StreamSolver: Dirichlet Point".into()));
        }
        let anchor = p.real("stream function penalty").unwrap_or(0.0);
        let cs = ctx.simulation_params().string("coordinate system").map(elmer_core::normalize_key).unwrap_or_default();
        let is_axis = cs == "axi symmetric" || cs == "cylindric symmetric";
        let stokes = p.logical("stokes stream function").unwrap_or(is_axis);
        if ctx.mesh.dim() != 2 {
            return Err(SolverError::Unsupported("StreamSolver: 3D meshes".into()));
        }
        let dim = 2;
        let (sperm, sdofs) = {
            let v = ctx.primary().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
            (v.perm.map(|p| p.to_vec()), v.dofs)
        };
        if sdofs != 1 {
            return Err(SolverError::Unsupported("StreamSolver: a variable of more than one dof".into()));
        }
        let in_perm = |k: u32| sperm.as_ref().is_none_or(|p| p[k as usize] != NO_DOF);

        ctx.default_initialize()?;
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let mut be = BasisEval::default();
        for e in elems {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let load = nodal_load(el.nodes, fperm.as_deref(), &fvals, fdofs, 1.0)?;
            let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no FE backend"))?;
            let mut local = LocalSystem::zeros(n);
            for qp in fe.quadrature.gauss_points(el.etype, None) {
                fe.basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let (b, g) = (&be.basis[..n], &be.dbasis_dx[..n]);
                let s = qp.weight * be.det_j;
                let mut lip = [0.0; 2];
                if stokes {
                    let radius = (0..n).fold(0.0, |acc, k| acc + xyz[k][0] * b[k]);
                    for (i, l) in lip.iter_mut().enumerate() {
                        *l = radius * (0..n).fold(0.0, |acc, k| acc + g[k][i] * load[k][i]);
                    }
                    lip[0] += (0..n).fold(0.0, |acc, k| acc + b[k] * load[k][0]);
                } else {
                    for (i, l) in lip.iter_mut().enumerate() {
                        *l = (0..n).fold(0.0, |acc, k| acc + g[k][i] * load[k][i]);
                    }
                }
                for i in 0..n {
                    for j in 0..n {
                        // MATMUL(dBasisdx, TRANSPOSE(dBasisdx)) over the 3 columns.
                        let m = (0..3).fold(0.0, |acc, a| acc + g[i][a] * g[j][a]);
                        local.stiff[i * n + j] += s * m;
                    }
                }
                for i in 0..n {
                    for j in 0..n {
                        local.stiff[i * n + j] += s * anchor * b[i] * b[j];
                    }
                }
                let sl = lip[..dim].iter().sum::<f64>();
                for i in 0..n {
                    local.force[i] += s * sl * b[i];
                }
            }
            ctx.default_update(&local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;

        // Every boundary element whose nodes are all in the system.
        let nb = ctx.mesh.n_bulk();
        for t in nb..nb + ctx.mesh.n_boundary() {
            let e = ElemId(t as u32);
            let el = ctx.mesh.element(e);
            if !el.nodes.iter().all(|&k| in_perm(k)) {
                continue;
            }
            if el.etype.0 != 202 {
                return Err(SolverError::Unsupported(format!("StreamSolver: boundary elements of type {}", el.etype.0)));
            }
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let load = nodal_load(el.nodes, fperm.as_deref(), &fvals, fdofs, -1.0)?;
            let fe = ctx.sim.backends.fe.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no FE backend"))?;
            let rule = fe.quadrature.gauss_points(el.etype, None);
            let mut local = LocalSystem::zeros(n);
            for qp in rule {
                fe.basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let b = be.basis[..n].to_vec();
                let mut s = qp.weight * be.det_j;
                if stokes {
                    s *= (0..n).fold(0.0, |acc, k| acc + xyz[k][0] * b[k]);
                }
                let lip: [f64; 2] = std::array::from_fn(|i| (0..n).fold(0.0, |acc, k| acc + b[k] * load[k][i]));
                // NormalVector(Element, Nodes, u, v, .TRUE.).
                let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [qp.uvw[0], qp.uvw[1]])?;
                let l = lip[0] * normal[0] + lip[1] * normal[1];
                for i in 0..n {
                    local.force[i] += s * l * b[i];
                }
            }
            ctx.default_update(&local, e)?;
        }
        ctx.default_finish_boundary_assembly()?;

        // UpdateDirichletDof(Perm(FirstNode), 0), then DefaultDirichletBCs.
        let (mut fixed, rows) = ctx.dirichlet_values_with_rows()?;
        if let Some(f) = first {
            let node = (f - 1) as u32;
            let k = sperm.as_ref().map_or(node, |p| p[node as usize]);
            if k != NO_DOF && !fixed.iter().any(|&(d, _)| d == k) {
                fixed.insert(0, (k, 0.0));
            }
        }
        ctx.apply_dirichlet_values_with_rows(&fixed, &rows)?;
        ctx.default_solve()?;

        let mut too_small = false;
        if shifting || scaling {
            let v = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
            if shifting {
                let min = v.values.iter().copied().fold(f64::INFINITY, f64::min);
                v.values.iter_mut().for_each(|x| *x -= min);
            }
            if scaling {
                let max = v.values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let min = v.values.iter().copied().fold(f64::INFINITY, f64::min);
                let coeff = max - min;
                if coeff < AEPS {
                    too_small = true;
                } else {
                    v.values.iter_mut().for_each(|x| *x /= coeff);
                }
            }
        }
        if too_small {
            ctx.warn(CALLER, "Maximum absolute value smaller than machine epsilon; cannot scale.");
        }
        Ok(())
    }
}
