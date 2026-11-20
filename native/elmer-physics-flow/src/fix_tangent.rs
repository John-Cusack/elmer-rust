//! `FixTangentVelo` (fem/src/modules/FixTangentVelo.F90): on the boundaries
//! with `Tangent Surface = True`, the correction to the `Active Coordinate`
//! component of the velocity that makes each node's velocity tangent to the
//! surface spanned by its two neighbouring edges, `FixVelo`. With `Fix
//! Tangent Velocity` the correction is subtracted from the velocity field.
//! The solver's norm is that of `FixVelo`.

use std::sync::Arc;

use elmer_core::{ElemId, FieldKind, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{compute_norm, NormSpec, Solver, SolverCtx, SolverError, Variable};

const CALLER: &str = "FixTangentVelo";
const VARIABLE: &str = "fixtangent_var";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(FixTangentVelo::default()))
}

#[derive(Default)]
pub struct FixTangentVelo {
    /// The SAVEd setup of the first visit.
    setup: Option<Setup>,
}

struct Setup {
    /// `ActiveCoord`, 0-based.
    active: usize,
    downstream_only: bool,
    velocity: String,
    /// `FixPerm`: node -> `FixVelo` slot.
    perm: Arc<[u32]>,
}

fn slot_index(ctx: &SolverCtx<'_>) -> usize {
    let id = ctx.solver_id();
    ctx.sim.slots.iter().position(|s| s.id == id).expect("running solver has a slot")
}

fn tangent_surface(ctx: &SolverCtx<'_>, e: ElemId) -> bool {
    ctx.bc_params(e).and_then(|bc| bc.logical("tangent surface")) == Some(true)
}

/// `SolveLinSys2x2( A, x, b )` (LinearAlgebra.F90:652).
fn solve2(a: [[f64; 2]; 2], b: [f64; 2]) -> Result<[f64; 2], SolverError> {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det == 0.0 {
        return Err(SolverError::fatal("SolveLinSys2x2", "Singular matrix, sorry!"));
    }
    let det = 1.0 / det;
    Ok([det * (a[1][1] * b[0] - a[0][1] * b[1]), det * (a[0][0] * b[1] - a[1][0] * b[0])])
}

/// `SolveLinSys3x3( A, x, b )` (LinearAlgebra.F90:678): elimination on the
/// largest entry of the first row, then [`solve2`].
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> Result<[f64; 3], SolverError> {
    // The pivot column, and the other two in order.
    let (p, r, s) = if a[0][0].abs() > a[0][1].abs() && a[0][0].abs() > a[0][2].abs() {
        (0, 1, 2)
    } else if a[0][1].abs() > a[0][2].abs() {
        (1, 0, 2)
    } else {
        (2, 0, 1)
    };
    let q = 1.0 / a[0][p];
    let sv = q * a[1][p];
    let t = q * a[2][p];
    let c = [[a[1][r] - sv * a[0][r], a[1][s] - sv * a[0][s]], [a[2][r] - t * a[0][r], a[2][s] - t * a[0][s]]];
    let g = [b[1] - sv * b[0], b[2] - t * b[0]];
    let y = solve2(c, g)?;
    let mut x = [0.0; 3];
    x[r] = y[0];
    x[s] = y[1];
    x[p] = q * (b[0] - a[0][r] * x[r] - a[0][s] * x[s]);
    Ok(x)
}

impl FixTangentVelo {
    /// The first visit: the parameters, `FixPerm` over the tangent-surface
    /// nodes in node order, and the `FixVelo` variable on it.
    fn first_visit(ctx: &mut SolverCtx<'_>) -> Result<Setup, SolverError> {
        let dim = ctx.mesh.dim();
        if dim != 3 {
            return Err(SolverError::fatal(CALLER, "Implemented only for 3D systems!"));
        }
        let p = ctx.params();
        let active = match p.integer("active coordinate") {
            Some(c) if (1..=3).contains(&c) => c as usize - 1,
            Some(c) => return Err(SolverError::fatal(CALLER, format!("Active Coordinate {c} out of range"))),
            None => dim - 1,
        };
        let downstream_only = p.logical("fix downstream only").unwrap_or(false);
        let velocity = p
            .string("velocity variable name")
            .map_or_else(|| "flow solution".to_string(), elmer_core::normalize_key);
        if ctx.variable(&velocity).is_none() {
            return Err(SolverError::fatal(CALLER, format!("Velocity field variable does not exist: {velocity}")));
        }
        let mut mark = vec![false; ctx.mesh.n_nodes()];
        let nb = ctx.mesh.n_bulk() as u32;
        for e in (nb..nb + ctx.mesh.n_boundary() as u32).map(ElemId) {
            if tangent_surface(ctx, e) {
                for &n in ctx.mesh.element(e).nodes {
                    mark[n as usize] = true;
                }
            }
        }
        let mut n = 0u32;
        let perm: Vec<u32> = mark
            .iter()
            .map(|&m| {
                if m {
                    n += 1;
                    n - 1
                } else {
                    NO_DOF
                }
            })
            .collect();
        // ParallelReduction(n).
        let mut total = [f64::from(n)];
        ctx.comm.sum_f64(&mut total);
        let n = total[0] as u64;
        ctx.info(CALLER, &format!("Total number of nodes to fix: {n}"), 4);
        if n == 0 {
            return Err(SolverError::fatal(CALLER, "No nodes to fix!"));
        }
        let perm: Arc<[u32]> = perm.into();
        let size = perm.iter().filter(|&&s| s != NO_DOF).count();
        ctx.sim
            .vars
            .add(Variable::new("FixVelo", FieldKind::Nodal, 1, Some(perm.clone()), vec![0.0; size]))?;
        Ok(Setup { active, downstream_only, velocity, perm })
    }
}

impl Solver for FixTangentVelo {
    /// `FixTangentVelo_init`: `No Matrix` and the `-nooutput -global
    /// FixTangent_var` solver variable.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let ix = slot_index(ctx);
        if ctx.sim.slots[ix].params.get("no matrix").is_none() {
            ctx.sim.slots[ix].params.insert("No Matrix", Value::Logical(true));
        }
        if ctx.sim.slots[ix].variable.is_none() {
            let mut v = Variable::new(VARIABLE, FieldKind::Global, 1, None, vec![0.0]);
            v.output = false;
            let id = ctx.sim.vars.add(v)?;
            ctx.sim.slots[ix].variable = Some(id);
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "Regulating surface velocity not to go out domain", 4);
        if self.setup.is_none() {
            self.setup = Some(Self::first_visit(ctx)?);
        }
        let st = self.setup.as_ref().expect("set above");
        let size = st.perm.iter().filter(|&&s| s != NO_DOF).count();
        let mut fix = vec![0.0f64; size];
        let (vdofs, vperm, velo) = {
            let f = ctx
                .variable(&st.velocity)
                .ok_or_else(|| SolverError::fatal(CALLER, format!("Velocity field variable does not exist: {}", st.velocity)))?;
            (f.dofs, f.perm.map(|p| p.to_vec()), f.values.to_vec())
        };
        if vdofs < 3 {
            return Err(SolverError::fatal(CALLER, "the velocity variable has fewer than 3 components"));
        }
        let vslot = |n: u32| match &vperm {
            Some(p) => p[n as usize],
            None => n,
        };
        let nb = ctx.mesh.n_bulk() as u32;
        for e in (nb..nb + ctx.mesh.n_boundary() as u32).map(ElemId) {
            if !tangent_surface(ctx, e) {
                continue;
            }
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            if !(n == 3 || n == 4) {
                return Err(SolverError::fatal(CALLER, "Only implemented for 3 and 4 nodes!"));
            }
            // NormalVector(Element, Nodes, Check=.TRUE.) at its default point.
            let uv = if el.etype.0 / 100 == 3 { [1.0 / 3.0, 1.0 / 3.0] } else { [0.0, 0.0] };
            let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, uv)?;
            let x: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            for i in 0..n {
                // GetVectorLocalSolution: zero where the velocity has no dof.
                let s = vslot(el.nodes[i]);
                let v: [f64; 3] = if s == NO_DOF {
                    [0.0; 3]
                } else {
                    std::array::from_fn(|k| velo[vdofs * s as usize + k])
                };
                let (i1, i2) = ((i + 1) % n, (i + n - 1) % n);
                let e1: [f64; 3] = std::array::from_fn(|k| x[i1][k] - x[i][k]);
                let e2: [f64; 3] = std::array::from_fn(|k| x[i2][k] - x[i][k]);
                if st.downstream_only {
                    let dot = |a: &[f64; 3]| v[0] * a[0] + v[1] * a[1] + v[2] * a[2];
                    if dot(&e1) < 0.0 && dot(&e2) < 0.0 {
                        continue;
                    }
                }
                let mut a = [[0.0; 3]; 3];
                for k in 0..3 {
                    a[k][0] = e1[k];
                    a[k][1] = e2[k];
                }
                a[st.active][2] = 1.0;
                let c = solve3(a, v)?;
                let j = st.perm[el.nodes[i] as usize] as usize;
                fix[j] = if normal[st.active] > 0.0 { fix[j].max(c[2]) } else { fix[j].min(c[2]) };
            }
        }
        if let Some(f) = ctx.sim.vars.get_mut("fixvelo") {
            f.values.copy_from_slice(&fix);
        }
        if ctx.params().logical("fix tangent velocity") == Some(true) {
            ctx.info(CALLER, "Applying the fix to current velocity field", 4);
            let f = ctx
                .variable_mut(&st.velocity)
                .ok_or_else(|| SolverError::fatal(CALLER, "velocity variable vanished"))?;
            for (node, &jf) in st.perm.iter().enumerate() {
                let jv = vslot(node as u32);
                if jf == NO_DOF || jv == NO_DOF {
                    continue;
                }
                let k = vdofs * jv as usize + st.active;
                f.values[k] -= fix[jf as usize];
            }
        }
        let spec = NormSpec::from_params(ctx.params(), 1);
        let nrm = compute_norm(&spec, &fix, fix.len(), ctx.comm);
        let ix = slot_index(ctx);
        if let Some(id) = ctx.sim.slots[ix].variable {
            let v = ctx.sim.vars.by_id_mut(id);
            v.values.iter_mut().for_each(|x| *x = nrm);
            v.conv.norm = nrm;
        }
        Ok(())
    }
}
