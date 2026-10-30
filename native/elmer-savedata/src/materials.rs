//! `SaveMaterials` (fem/src/modules/SaveData/SaveMaterials.F90): material
//! (or body force) keywords as nodal fields.
//!
//! On every visit, for each `Parameter i`:
//! - the first time, a nodal variable of that name is created, unless one
//!   exists, over the nodes of the elements whose material (body force for
//!   the first `Body Force Parameters` names) has the keyword, numbered in
//!   node order;
//! - its values are set element by element, in mesh order, from the
//!   keyword at the element's nodes, so a node shared by two bodies keeps
//!   the later element's value.

use std::sync::Arc;

use elmer_core::{ElemId, FieldKind, ParamMap, Params, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, Variable};

const CALLER: &str = "SaveMaterials";

pub struct SaveMaterials {
    names: Vec<String>,
    body_force_params: usize,
}

impl SaveMaterials {
    pub fn from_params(params: &ParamMap) -> Result<Self, SolverError> {
        let names: Vec<String> = (1..)
            .map_while(|i| params.string(&format!("parameter {i}")).map(str::to_string))
            .collect();
        let body_force_params = params.integer("body force parameters").unwrap_or(0).max(0) as usize;
        Ok(SaveMaterials { names, body_force_params })
    }
}

/// The section `ParamNo` reads on element `e`: its body's material, or its
/// body force for the first `Body Force Parameters` names.
fn section<'a>(ctx: &'a SolverCtx<'_>, e: ElemId, body_force: bool) -> Option<&'a ParamMap> {
    let b = ctx.body_params(e)?;
    if body_force {
        b.body_force
    } else {
        b.material
    }
}

impl Solver for SaveMaterials {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if self.names.is_empty() {
            ctx.log.warn(CALLER, "No parameters found: No fields will be created");
            return Ok(());
        }
        let n_elems = ctx.mesh.n_bulk() + ctx.mesh.n_boundary();
        let n_nodes = ctx.mesh.n_nodes();
        for (no, name) in self.names.iter().enumerate() {
            if ctx.sim.vars.resolve(name).is_some() {
                continue;
            }
            let body_force = no < self.body_force_params;
            let mut perm = vec![NO_DOF; n_nodes];
            for e in (0..n_elems as u32).map(ElemId) {
                if section(ctx, e, body_force).is_some_and(|p| p.get(name).is_some()) {
                    for &node in ctx.mesh.element(e).nodes {
                        perm[node as usize] = 0;
                    }
                }
            }
            let mut j = 0u32;
            for p in perm.iter_mut().filter(|p| **p != NO_DOF) {
                *p = j;
                j += 1;
            }
            if j == 0 {
                ctx.log.warn(CALLER, &format!("Parameter {name} not present in any material"));
                continue;
            }
            ctx.log.info(CALLER, &format!("Parameter > {name} < defined with {j} dofs"), 5);
            let var = Variable::new(name, FieldKind::Nodal, 1, Some(Arc::from(perm)), vec![0.0; j as usize]);
            ctx.sim.vars.add(var)?;
        }
        for (no, name) in self.names.iter().enumerate() {
            let body_force = no < self.body_force_params;
            let Some(id) = ctx.sim.vars.id(name) else { continue };
            let Some(perm) = ctx.sim.vars.by_id(id).perm.clone() else { continue };
            for e in (0..n_elems as u32).map(ElemId) {
                let nodes = ctx.mesh.element(e).nodes.to_vec();
                if nodes.iter().any(|&n| perm[n as usize] == NO_DOF) {
                    continue;
                }
                let Some(p) = section(ctx, e, body_force) else { continue };
                let Some(values) = ctx.real_at_nodes(p, name, e)? else { continue };
                let var = ctx.sim.vars.by_id_mut(id);
                for (&n, v) in nodes.iter().zip(values) {
                    var.values[perm[n as usize] as usize] = v;
                }
            }
        }
        Ok(())
    }
}
