//! `SaveProjection` (SaveData/SaveProjection.F90): the solver's mortar
//! projectors applied to fields, each result saved as a new field.
//!
//! For `Variable i` (a field, or `Coordinate 1..3` from the mesh) the target
//! `Target Variable i` (default `Projection <name>`) gets, row by row of each
//! boundary's projector, minus the projector times the field's master
//! values on the row's slave node (`Project To Slave`), the transposed
//! product on the master nodes (`Project To Master`), or the plain product,
//! each optionally divided by the projector's slave weights (`Normalize`,
//! default true). The solver's norm is the combined L2 norm of the targets.
//! Refused: `PMatrix` (periodic) projectors.

use std::path::PathBuf;
use std::sync::Arc;

use elmer_core::{FieldKind, NullLogger, ParamMap, Params, Value, NO_DOF};
use elmer_mortar::{generate_mortar_projectors, Projector, ProjectorEnv};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable};

const CALLER: &str = "SaveProjection";
pub(crate) const PROCEDURE: (&str, &str) = ("SaveData", "SaveProjection");

pub(crate) fn is_section(params: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(params)
        .is_some_and(|(l, e)| l.eq_ignore_ascii_case(PROCEDURE.0) && e.eq_ignore_ascii_case(PROCEDURE.1))
}

/// `SaveProjection_init`: `Variable = -nooutput -global SaveProjection_var`
/// unless given.
pub(crate) fn init_section(p: &mut ParamMap) {
    if p.get("variable").is_none() {
        p.insert("Variable", Value::Str("-nooutput -global SaveProjection_var".into()));
    }
}

pub(crate) fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| {
        Ok(Box::new(SaveProjection { coord_dim: 0, case_dir: PathBuf::new(), discont_perm: None }))
    });
}

pub(crate) struct SaveProjection {
    pub coord_dim: usize,
    pub case_dir: PathBuf,
    pub discont_perm: Option<Arc<Vec<u32>>>,
}

/// One source field: its values, permutation (none for the coordinates)
/// and components.
struct Source {
    values: Vec<f64>,
    perm: Option<Arc<[u32]>>,
    dofs: usize,
}

fn source(ctx: &SolverCtx<'_>, name: &str) -> Option<Source> {
    let key = elmer_core::normalize_key(name);
    if let Some(c) = ["coordinate 1", "coordinate 2", "coordinate 3"].iter().position(|k| *k == key) {
        if ctx.sim.vars.get(&key).is_none() {
            let values = (0..ctx.mesh.n_nodes() as u32).map(|n| ctx.mesh.coord(n)[c]).collect();
            return Some(Source { values, perm: None, dofs: 1 });
        }
    }
    let v = ctx.sim.vars.get(&key)?;
    Some(Source { values: v.values.clone(), perm: v.perm.clone(), dofs: v.dofs })
}

impl Solver for SaveProjection {
    fn init(&mut self, _ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "Creating selected projected values as fields", 8);
        let params = ctx.params().clone();
        let model = &ctx.sim.model;
        if model.boundary_conditions.iter().any(|b| b.get("periodic bc").is_some() || b.get("anti periodic bc").is_some())
        {
            return Err(SolverError::Unsupported("SaveProjection: periodic (PMatrix) projectors".into()));
        }
        // The solver's linear projectors (made when it was set up) and
        // GenerateProjectors(Model, Solver, Nonlinear = .TRUE.)'s.
        let projectors: Vec<Option<Projector>> = {
            let bcs: Vec<&dyn Params> = model.boundary_conditions.iter().map(|b| b as &dyn Params).collect();
            let bc_of = |tag: i32| model.bc_of_tag(tag);
            let env = ProjectorEnv {
                mesh: ctx.mesh,
                bcs: &bcs,
                bc_of: &bc_of,
                simulation: &model.simulation,
                solver: Some(ctx.params()),
                coord_dim: self.coord_dim,
                log: &NullLogger,
                discont_perm: self.discont_perm.as_deref().map(Vec::as_slice),
                edges: None,
            };
            let timestep = ctx.sim.vars.scalar("timestep").map_or(0, |t| t.round() as i32);
            let mut all: Vec<Option<Projector>> = vec![None; bcs.len()];
            for nonlinear in [false, true] {
                match generate_mortar_projectors(&env, &[], nonlinear, &self.case_dir, timestep)? {
                    Some((p, false)) => {
                        for (slot, q) in all.iter_mut().zip(p) {
                            if q.is_some() {
                                *slot = q;
                            }
                        }
                    }
                    Some((_, true)) => return Err(SolverError::Unsupported("Save Projector And Stop".into())),
                    None => {}
                }
            }
            all
        };
        let mut nvar = 0;
        while let Some(name) = params.string(&format!("variable {}", nvar + 1)) {
            if source(ctx, name).is_none() {
                ctx.warn(CALLER, "Requested variable does not exist!");
                break;
            }
            nvar += 1;
        }
        ctx.info(CALLER, &format!("Saving projections from {nvar} fields"), 4);
        if nvar == 0 {
            return Ok(());
        }
        let mut nrm = 0.0f64;
        for i in 1..=nvar {
            let name = params.string(&format!("variable {i}")).unwrap_or_default().to_string();
            let src = source(ctx, &name).expect("checked above");
            ctx.info(CALLER, &format!("Doing variable: {}", elmer_core::normalize_key(&name)), 8);
            let target = params
                .string(&format!("target variable {i}"))
                .map_or_else(|| format!("Projection {name}"), str::to_string);
            let normalize = params.logical(&format!("normalize {i}")).unwrap_or(true);
            let to_slave = params.logical(&format!("project to slave {i}")).unwrap_or(false);
            let to_master = params.logical(&format!("project to master {i}")).unwrap_or(false);
            let active: Option<Vec<i32>> = match params.get(&format!("active projectors {i}")) {
                Some(Value::Integer(v)) => Some(v.clone()),
                _ => None,
            };
            let tid = match ctx.sim.vars.id(&target) {
                Some(id) => id,
                None => {
                    // VariableAddVector: the source's permutation, or the
                    // identity over its values.
                    let perm: Arc<[u32]> = match &src.perm {
                        Some(p) => Arc::clone(p),
                        None => (0..src.values.len() as u32).collect::<Vec<_>>().into(),
                    };
                    let var = Variable::new(&target, FieldKind::Nodal, src.dofs, Some(perm), vec![0.0; src.values.len()]);
                    ctx.sim.vars.add(var)?
                }
            };
            let tv = ctx.sim.vars.by_id(tid).values.len();
            let mut out = vec![0.0f64; tv];
            let mut weight = vec![0.0f64; if normalize { tv } else { 0 }];
            let dofs = src.dofs;
            let map = |node: u32| -> Option<usize> {
                match &src.perm {
                    Some(p) => p.get(node as usize).copied().filter(|&s| s != NO_DOF).map(|s| s as usize),
                    None => Some(node as usize),
                }
            };
            for (b, p) in projectors.iter().enumerate() {
                if active.as_ref().is_some_and(|a| !a.contains(&(b as i32 + 1))) {
                    continue;
                }
                let Some(p) = p else { continue };
                let n = p.n_rows();
                if n == 0 {
                    continue;
                }
                // IsInvInvPerm over the columns.
                let m = p.cols.iter().map(|&c| c as usize + 1).max().unwrap_or(0);
                let mut is_slave = vec![false; m];
                for &k in &p.inv_perm {
                    if k != NO_DOF && (k as usize) < m {
                        is_slave[k as usize] = true;
                    }
                }
                for k in 0..dofs {
                    for r in 0..n {
                        let node = p.inv_perm[r];
                        if node == NO_DOF || node as usize >= m {
                            continue;
                        }
                        let Some(pi) = map(node) else { continue };
                        let pi = dofs * pi + k;
                        for (col, val) in p.row(r) {
                            let Some(pj) = map(col) else { continue };
                            let pj = dofs * pj + k;
                            let actj = is_slave[col as usize];
                            if to_slave {
                                if !actj {
                                    out[pi] += -val * src.values[pj];
                                }
                            } else if !to_master {
                                out[pi] += -val * src.values[pj];
                            }
                            if normalize && actj {
                                weight[pi] += val;
                            }
                            if to_master && !actj {
                                out[pj] += val * src.values[pi];
                                if normalize {
                                    weight[pj] += val;
                                }
                            }
                        }
                    }
                }
            }
            if normalize {
                for (o, w) in out.iter_mut().zip(&weight) {
                    if w.abs() > f64::EPSILON {
                        *o /= w;
                    }
                }
            }
            nrm += out.iter().map(|v| v * v).sum::<f64>();
            let var = ctx.sim.vars.by_id_mut(tid);
            var.values = out;
            var.values_changed = true;
        }
        let nrm = nrm.sqrt();
        ctx.info(CALLER, &format!("Combined L2 norm of all projected fields: {}", elmer_solver::fortran_fmt::es(nrm, 12, 3)), 4);
        if let Some(id) = ctx.sim.slots.iter().find(|s| s.id == ctx.solver_id()).and_then(|s| s.variable) {
            let v = ctx.sim.vars.by_id_mut(id);
            if v.values.len() == 1 {
                v.values[0] = nrm;
                v.conv.norm = nrm.abs();
            }
        }
        Ok(())
    }
}
