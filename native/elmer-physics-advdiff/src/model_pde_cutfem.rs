//! `AdvDiffSolver` of ModelPDECutFEM.F90 (`Procedure = "ModelPDECutFEM"
//! "AdvDiffSolver"`): ModelPDE's advection-diffusion-reaction equation on
//! the part of the mesh a zero levelset cuts out, the elements it crosses
//! split on the fly (`elmer_solver::cutfem`): the pieces and the interface
//! elements between them are assembled at the cut dofs, the unsplit
//! elements as usual, each with the material and body force of the body
//! CutFEM gives it.
//!
//! Without `CutFEM` it is ModelPDE's solver with the same local matrices.
//! Refused: coefficients that are not constant on a split element (their
//! cut nodes are not mesh nodes), more than one nonlinear iteration (the
//! module never re-initialises the matrix between them) and bubbles.

use elmer_core::{ElemId, ElementType, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_fe::{ElmerQuadrature, LagrangeBasis};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::fe::Fe;
use crate::model_pde::{model_pde_boundary_local, model_pde_bulk_local, ModelPdeBcNodal, ModelPdeNodal};

/// `Procedure = "ModelPDECutFEM" "AdvDiffSolver"`.
pub const PROCEDURE: (&str, &str) = ("ModelPDECutFEM", "AdvDiffSolver");

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| Ok(Box::new(ModelPdeCutFem { local: LocalSystem::zeros(0) })));
}

pub fn is_section(p: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(p)
        .is_some_and(|(l, e)| l.eq_ignore_ascii_case(PROCEDURE.0) && e.eq_ignore_ascii_case(PROCEDURE.1))
}

/// `AdvDiffSolver_Init`: with `CutFEM`, `Use Global Mass Matrix` and
/// `CutFEM Solver` (`ListAddNewLogical`: a keyword already given wins).
pub fn init_section(p: &mut ParamMap) {
    if p.logical("cutfem") == Some(true) {
        for k in ["Use Global Mass Matrix", "CutFEM Solver"] {
            if p.get(k).is_none() {
                p.insert(k, Value::Logical(true));
            }
        }
    }
}

struct ModelPdeCutFem {
    local: LocalSystem,
}

/// `GetReal(section, key, Found, Element)` on a split element: constants
/// only, zero when absent.
fn constant(p: Option<&ParamMap>, key: &str, n: usize) -> Result<Vec<f64>, SolverError> {
    let Some(p) = p else { return Ok(vec![0.0; n]) };
    match p.get(key) {
        None => Ok(vec![0.0; n]),
        Some(Value::Real { data, .. }) if data.len() == 1 => Ok(vec![data[0]; n]),
        Some(_) => Err(SolverError::Unsupported(format!("CutFEM: \"{key}\" that is not a constant on a split element"))),
    }
}

/// The coefficients of `LocalMatrix` from the body `body` (its material
/// and body force): at the nodes of mesh element `e`, or constants for a
/// piece of a split element.
pub(crate) fn nodal(ctx: &SolverCtx<'_>, body: i32, e: Option<ElemId>, n: usize) -> Result<ModelPdeNodal, SolverError> {
    let bp = usize::try_from(body).ok().and_then(|b| ctx.sim.model.body_params(b));
    let (m, bf) = (bp.as_ref().and_then(|b| b.material), bp.as_ref().and_then(|b| b.body_force));
    let get = |p: Option<&ParamMap>, key: &str| -> Result<Vec<f64>, SolverError> {
        match (e, p) {
            (Some(e), Some(p)) => Ok(ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0; n])),
            (Some(_), None) => Ok(vec![0.0; n]),
            (None, p) => constant(p, key, n),
        }
    };
    let dim = ctx.mesh.dim();
    let mut velocity = vec![[0.0; 3]; n];
    for i in 0..dim {
        for (vel, x) in velocity.iter_mut().zip(get(m, &format!("convection velocity {}", i + 1))?) {
            vel[i] = x;
        }
    }
    Ok(ModelPdeNodal {
        diffusion: get(m, "diffusion coefficient")?,
        convection: get(m, "convection coefficient")?,
        reaction: get(m, "reaction coefficient")?,
        time: get(m, "time derivative coefficient")?,
        velocity,
        source: get(bf, "field source")?,
    })
}

pub(crate) fn bc_nodal(ctx: &SolverCtx<'_>, p: Option<&ParamMap>, e: Option<ElemId>, n: usize) -> Result<ModelPdeBcNodal, SolverError> {
    let get = |key: &str| -> Result<Vec<f64>, SolverError> {
        match (e, p) {
            (Some(e), Some(p)) => Ok(ctx.real_at_nodes(p, key, e)?.unwrap_or_else(|| vec![0.0; n])),
            (Some(_), None) => Ok(vec![0.0; n]),
            (None, p) => constant(p, key, n),
        }
    };
    Ok(ModelPdeBcNodal { flux: get("field flux")?, robin: get("robin coefficient")?, external: get("external field")? })
}

fn etype_of(n: usize) -> ElementType {
    ElementType(match n {
        2 => 202,
        3 => 303,
        _ => 404,
    })
}

impl Solver for ModelPdeCutFem {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if ctx.primary().is_none_or(|f| f.dofs != 1) {
            return Err(SolverError::Unsupported("ModelPDECutFEM solves a scalar variable".into()));
        }
        if ctx.solver_dofs().is_some() {
            return Err(SolverError::Unsupported("ModelPDECutFEM with p elements".into()));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.default_start()?;
        let do_cut = ctx.params().logical("cutfem") == Some(true);
        if !do_cut {
            return Err(SolverError::Unsupported("ModelPDECutFEM without CutFEM (DefaultInitialize is never called)".into()));
        }
        let maxiter = ctx.params().integer("nonlinear system max iterations").map_or(1, |n| n.max(1));
        if maxiter > 1 {
            return Err(SolverError::Unsupported(
                "ModelPDECutFEM with more than one nonlinear iteration (the matrix is not re-initialised)".into(),
            ));
        }
        let dim = ctx.mesh.dim();
        let basis = LagrangeBasis::new(dim);
        let quad = ElmerQuadrature;
        let fe = Fe::new(&basis, &quad);
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let bnd: Vec<ElemId> = ctx.boundary_elements().collect();

        for &e in &elems {
            let (_, active) = ctx.cutfem_check(e)?;
            if !active {
                continue;
            }
            let mut split = false;
            loop {
                let sub = ctx.cutfem_bulk(e)?;
                if !sub.cut {
                    break;
                }
                split = true;
                if ctx.cutfem_in_perm(&sub.nodes)? {
                    let xyz = ctx.cutfem_xyz(&sub.nodes)?;
                    let nod = nodal(ctx, sub.body, None, sub.nodes.len())?;
                    model_pde_bulk_local(&fe, etype_of(sub.nodes.len()), &xyz, &nod, dim, None, &mut self.local)?;
                    let mass = self.local.mass.take();
                    if transient {
                        if let Some(m) = mass {
                            ctx.cutfem_update_mass(&m, &sub.nodes)?;
                        }
                    }
                    ctx.cutfem_update(&self.local, &sub.nodes)?;
                }
                if !sub.more {
                    break;
                }
            }
            if split {
                continue;
            }
            let el = ctx.mesh.element(e);
            let (etype, nodes) = (el.etype, el.nodes.to_vec());
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let body = ctx.cutfem_body(e)?;
            let nod = nodal(ctx, body, Some(e), nodes.len())?;
            model_pde_bulk_local(&fe, etype, &xyz, &nod, dim, None, &mut self.local)?;
            if transient {
                ctx.default_first_order_time(&mut self.local, e, dt)?;
            }
            self.local.mass = None;
            ctx.default_update(&self.local, e)?;
        }
        ctx.default_finish_bulk_assembly()?;

        // The interface elements of the split elements.
        let interface = ctx.cutfem_interface_bc().cloned();
        for &e in &elems {
            while let Some(ie) = ctx.cutfem_interface(e)? {
                if ctx.cutfem_in_perm(&ie.nodes)? {
                    let xyz = ctx.cutfem_xyz(&ie.nodes)?;
                    let bc = bc_nodal(ctx, interface.as_ref(), None, 2)?;
                    model_pde_boundary_local(&fe, ElementType(202), &xyz, &bc, &mut self.local)?;
                    ctx.cutfem_update(&self.local, &ie.nodes)?;
                }
                if !ie.more {
                    break;
                }
            }
        }

        let mut dofs = Vec::new();
        for &e in &bnd {
            ctx.element_dofs(e, &mut dofs);
            if dofs.contains(&NO_DOF) {
                continue;
            }
            let el = ctx.mesh.element(e);
            let (etype, nodes) = (el.etype, el.nodes.to_vec());
            let Some(p) = ctx.bc_params(e).cloned() else { continue };
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let bc = bc_nodal(ctx, Some(&p), Some(e), nodes.len())?;
            model_pde_boundary_local(&fe, etype, &xyz, &bc, &mut self.local)?;
            ctx.default_update(&self.local, e)?;
        }
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        ctx.default_finish()
    }
}
