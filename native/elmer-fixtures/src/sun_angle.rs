//! `SunAngleSolver` (fem/src/modules/SunAngle.F90, built into the
//! `SunAngle` test as `./SunAngle`): the critical sun angle of a height
//! field. With a `Test Height` in a material, the nodes' z coordinates
//! first become it (`SetTestHeight`); then, for the sun direction `Test
//! Angle` in the xy plane, each mesh edge (`FindMeshEdges2D`, in first-use
//! order with that element's orientation) shades every other node the ray
//! from it towards the sun crosses, and the node keeps the largest
//! `tan(phi) = dz / c(1)` of them (`SetShadingAngle`); nodes no edge
//! shades get 0.
//!
//! Serial only: the ridges of other partitions (`Enforce Parallel`, more
//! than one rank) and the multi-dof variant (several angles) are refused.

use elmer_core::{ElemId, Params};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "SunAngleSolver";

/// `SunAngleSolver_init`: `Variable` defaults to `CriticalTan`, and the
/// solver has no matrix.
pub fn init_section(p: &mut elmer_core::ParamMap) {
    if p.get("variable").is_none() {
        p.insert("Variable", elmer_core::Value::Str("CriticalTan".into()));
    }
    p.insert("No Matrix", elmer_core::Value::Logical(true));
}

pub fn factory(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(SunAngle))
}

struct SunAngle;

/// `FindMeshEdges2D` of the bulk elements: consecutive node pairs, numbered
/// in first-use order, each with the orientation of its first element.
fn edges(ctx: &SolverCtx<'_>) -> Vec<[u32; 2]> {
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for t in 0..ctx.mesh.n_bulk() {
        let nodes = ctx.mesh.element(ElemId(t as u32)).nodes;
        let n = nodes.len();
        for k in 0..n {
            let (a, b) = (nodes[k], nodes[(k + 1) % n]);
            if seen.insert((a.min(b), a.max(b))) {
                out.push([a, b]);
            }
        }
    }
    out
}

impl Solver for SunAngle {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        if p.logical("enforce parallel") == Some(true) || ctx.comm.size() > 1 {
            return Err(SolverError::Unsupported("SunAngleSolver in parallel (the ridges of other partitions)".into()));
        }
        if ctx.mesh.dim() != 2 {
            return Err(SolverError::Unsupported("SunAngleSolver on a 3D mesh".into()));
        }
        let name = ctx.variable_name().ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?.to_string();
        let id = ctx.sim.vars.id(&name).ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
        if ctx.sim.vars.by_id(id).dofs != 1 {
            return Err(SolverError::Unsupported("SunAngleSolver with several angles (a multi-dof variable)".into()));
        }
        let nn = ctx.mesh.n_nodes();
        let mut xyz: Vec<[f64; 3]> = (0..nn as u32).map(|k| ctx.mesh.coord(k)).collect();
        let edges = edges(ctx);

        // SetTestHeight: the bulk elements of a 2D mesh.
        if ctx.sim.model.materials.iter().any(|m| m.get("test height").is_some()) {
            for t in 0..ctx.mesh.n_bulk() {
                let e = ElemId(t as u32);
                let n = ctx.mesh.element(e).nodes.len();
                let h = match ctx.body_params(e).and_then(|b| b.material) {
                    Some(m) => ctx.real_at_nodes(m, "test height", e)?.unwrap_or_else(|| vec![0.0; n]),
                    None => vec![0.0; n],
                };
                let v = ctx.sim.vars.by_id(id);
                let nodes = ctx.mesh.element(e).nodes;
                if nodes.iter().any(|&k| v.slot(k).is_none()) {
                    continue;
                }
                for (i, &k) in nodes.iter().enumerate() {
                    xyz[k as usize][2] = h[i];
                }
            }
            if !elmer_physics_meshtools::motion::set_coords(&xyz) {
                return Err(SolverError::NotWired("node coordinates the solver can move".into()));
            }
        }

        let eps = 1.0e-6;
        let alpha = ctx.creal(&p, "test angle")?.unwrap_or(0.0);
        let (a11, a21) = (alpha.cos(), alpha.sin());
        let v = ctx.sim.vars.by_id(id);
        let slots: Vec<Option<usize>> = (0..nn as u32).map(|k| v.slot(k)).collect();
        let mut sun = vec![-100.0f64; v.values.len()];
        for &[i0, i1] in &edges {
            let (i0, i1) = (i0 as usize, i1 as usize);
            if slots[i0].is_none() || slots[i1].is_none() {
                continue;
            }
            let (x0, y0, z0) = ([xyz[i0][0], xyz[i1][0]], [xyz[i0][1], xyz[i1][1]], [xyz[i0][2], xyz[i1][2]]);
            let dz0 = z0[1] - z0[0];
            let a12 = x0[0] - x0[1];
            let a22 = y0[0] - y0[1];
            for j in 0..nn {
                let Some(k) = slots[j] else { continue };
                if j == i0 || j == i1 {
                    continue;
                }
                let b = [x0[0] - xyz[j][0], y0[0] - xyz[j][1]];
                // Solve2x2
                let det = a11 * a22 - a12 * a21;
                if det.abs() < f64::EPSILON {
                    continue;
                }
                let det = 1.0 / det;
                let c1 = det * (a22 * b[0] - a12 * b[1]);
                let c2 = det * (a11 * b[1] - a21 * b[0]);
                if c1 < eps || c2 < -eps || c2 > 1.0 + eps {
                    continue;
                }
                let dz = (z0[0] + c2 * dz0) - xyz[j][2];
                let tanphi = dz / c1;
                sun[k] = sun[k].max(tanphi);
            }
        }
        for s in &mut sun {
            if *s < -99.0 {
                *s = 0.0;
            }
        }
        ctx.sim.vars.by_id_mut(id).values = sun;
        Ok(())
    }
}
