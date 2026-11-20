//! `Poisson.f90` of the `CurvedBndryPFEM` and `SD_CurvedBndryPFEM` tests
//! (`Procedure = "./Poisson" "PoissonSolver"`): the square mesh's nodes are
//! first moved onto the harmonic extension of `sqrt(2) (x, y) / r`, which
//! maps the square's boundary onto the circle of radius `sqrt(2)`: two
//! p-element Poisson solves without a load, the boundary condition's
//! `Potential` the procedures `circx` and `circy` (`ListAddConstReal(BC,
//! 'Potential', 0, GetProcAddr(...))`), the nodes taking the solutions'
//! values at the vertices. Then the Poisson equation with the body force's
//! `Source` and `Potential = 0` on the moved mesh.
//!
//! The procedure-valued Dirichlet condition of the p elements is that of
//! `DefaultDirichletBCs` for a 2D parent (DefUtils.F90:6520): the nodes take
//! the procedure at their coordinates, and the edge dofs of each boundary
//! element the `LocalBcBDOFs` L2 projection of it with the nodal values
//! moved to the right-hand side (`SolveLinSys` of the edge block). The local
//! matrices take the element's global dofs only (`GetElementNOFDOFs`).
//! Refused: edge blocks beyond 2x2 (p > 3) and 3D.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use elmer_core::{BasisEval, ElemId, LocalSystem, Value};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "PoissonSolver";

pub fn factory(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(CurvedPoisson))
}

struct CurvedPoisson;

/// `CircX` / `CircY`: `sqrt(2) x / r`, `sqrt(2) y / r` (0 at the origin).
fn circ(c: usize, x: f64, y: f64) -> f64 {
    let r = (x * x + y * y).sqrt();
    if r != 0.0 {
        2f64.sqrt() * [x, y][c] / r
    } else {
        0.0
    }
}

/// The element loop with `LocalMatrix(STIFF, FORCE, LOAD, Element, n, nd)`,
/// then `DefaultFinishAssembly`. `geometry` holds the coordinates of every
/// dof slot once the mesh is moved (a p mesh's `Mesh % Nodes` has an entry
/// per dof: `GetElementNodes` fills the element's dofs beyond its nodes from
/// them, and `ElementMetric` maps the element with all its functions).
fn assemble(ctx: &mut SolverCtx<'_>, with_source: bool, geometry: Option<&[Vec<f64>; 2]>) -> Result<(), SolverError> {
    ctx.default_initialize()?;
    let active: Vec<ElemId> = ctx.active_elements().collect();
    let mut slots = Vec::new();
    let mut ev = BasisEval::default();
    for e in active {
        let d = ctx.solver_dofs().ok_or_else(|| SolverError::fatal(CALLER, "no p elements"))?;
        d.element_slots(e, &mut slots);
        let nd = slots.len();
        let nodes = d.element_nodes(e);
        let basis = d.element_basis_scalar(e).map_err(|m| SolverError::fatal("ElementInfo", m))?;
        let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
        let n = nodes.len();
        let etype = ctx.mesh.element(e).etype;
        let xyz: Vec<[f64; 3]> = match geometry {
            Some(g) => slots.iter().map(|&s| [g[0][s as usize], g[1][s as usize], 0.0]).collect(),
            None => nodes.iter().map(|&k| ctx.mesh.coord(k)).collect(),
        };
        let load = match (with_source, ctx.body_params(e).and_then(|b| b.body_force)) {
            (true, Some(bf)) => ctx.real_at_nodes(bf, "source", e)?.unwrap_or_else(|| vec![0.0; n]),
            _ => vec![0.0; n],
        };
        let mut local = LocalSystem::zeros(nd);
        for qp in &rule {
            basis.eval(etype, &xyz, qp.uvw, &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let w = qp.weight * ev.det_j;
            let mut load_at_ip = 0.0;
            for i in 0..n {
                load_at_ip += ev.basis[i] * load[i];
            }
            for i in 0..nd {
                let gi = ev.dbasis_dx[i];
                for j in 0..nd {
                    let gj = ev.dbasis_dx[j];
                    local.stiff[i * nd + j] += w * (gi[0] * gj[0] + gi[1] * gj[1] + gi[2] * gj[2]);
                }
                local.force[i] += w * load_at_ip * ev.basis[i];
            }
        }
        ctx.default_update(&local, e)?;
    }
    ctx.default_finish_boundary_assembly()
}

/// The Dirichlet values the procedure `circ(c, ...)` gives the boundary
/// elements of condition `bc` on p parents: the nodes at their
/// coordinates, then each element's edge dofs (`LocalBcBDOFs`).
fn procedure_dirichlet(ctx: &SolverCtx<'_>, bc: usize, c: usize) -> Result<Vec<(u32, f64)>, SolverError> {
    let d = ctx.solver_dofs().ok_or_else(|| SolverError::fatal(CALLER, "no p elements"))?;
    let id = ctx.variable_name().and_then(|n| ctx.sim.vars.id(n)).ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
    let v = ctx.sim.vars.by_id(id);
    let nb = ctx.mesh.n_bulk();
    let bnd: Vec<ElemId> = (nb..nb + ctx.mesh.n_boundary()).map(|i| ElemId(i as u32)).filter(|&e| ctx.bc_id(e) == Some(bc)).collect();
    let mut nodal = std::collections::BTreeMap::new();
    for &e in &bnd {
        for &k in ctx.mesh.element(e).nodes {
            if let Some(s) = v.slot(k) {
                let p = ctx.mesh.coord(k);
                nodal.insert(s as u32, circ(c, p[0], p[1]));
            }
        }
    }
    let mut fixed: Vec<(u32, f64)> = nodal.iter().map(|(&s, &x)| (s, x)).collect();
    let mut extra = Vec::new();
    let mut ev = BasisEval::default();
    for &e in &bnd {
        d.boundary_extra_slots(e, &mut extra);
        if extra.is_empty() {
            continue;
        }
        // The element's nodes in the order its p basis takes them.
        let etype = ctx.mesh.element(e).etype;
        let nodes = d.element_nodes(e);
        let n = nodes.len();
        let nd = n + extra.len();
        let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
        let basis = d.element_basis_scalar(e).map_err(|m| SolverError::fatal("ElementInfo", m))?;
        let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
        let mut stiff = vec![0.0; nd * nd];
        let mut work = vec![0.0; nd];
        for qp in &rule {
            basis.eval(etype, &xyz, qp.uvw, &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let s = qp.weight * ev.det_j;
            // SUM(Basis(1:nd) * Nodes % x(1:nd)): the node arrays are zero
            // beyond the element's nodes.
            let (mut xip, mut yip) = (0.0, 0.0);
            for i in 0..nd {
                let p = if i < n { xyz[i] } else { [0.0; 3] };
                xip += ev.basis[i] * p[0];
                yip += ev.basis[i] * p[1];
            }
            let load = circ(c, xip, yip);
            for p in 0..nd {
                for q in 0..nd {
                    stiff[p * nd + q] += s * ev.basis[p] * ev.basis[q];
                }
                work[p] += s * load * ev.basis[p];
            }
        }
        for l in 0..n {
            let Some(sl) = v.slot(nodes[l]) else {
                return Err(SolverError::Unsupported("a boundary node outside the p system".into()));
            };
            let sv = nodal[&(sl as u32)];
            for k in n..nd {
                work[k] -= sv * stiff[k * nd + l];
                stiff[k * nd + l] = 0.0;
            }
        }
        if nd == n + 1 {
            work[nd - 1] /= stiff[(nd - 1) * nd + nd - 1];
        } else if nd == n + 2 {
            // SolveLinSys2x2.
            let a = |i: usize, j: usize| stiff[(n + i) * nd + n + j];
            let b = [work[n], work[n + 1]];
            let det = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
            if det == 0.0 {
                return Err(SolverError::fatal("SolveLinSys2x2", "Singular matrix, sorry!"));
            }
            let det = 1.0 / det;
            work[n] = det * (a(1, 1) * b[0] - a(0, 1) * b[1]);
            work[n + 1] = det * (a(0, 0) * b[1] - a(1, 0) * b[0]);
        } else {
            return Err(SolverError::Unsupported("p-element procedure Dirichlet values with more than two edge dofs".into()));
        }
        for (k, &s) in extra.iter().enumerate() {
            fixed.push((s, work[n + k]));
        }
    }
    Ok(fixed)
}

impl Solver for CurvedPoisson {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::Unsupported("the CurvedBndryPFEM Poisson plugin without p elements".into()));
        }
        if ctx.mesh.dim() != 2 {
            return Err(SolverError::Unsupported("the CurvedBndryPFEM Poisson plugin outside 2D".into()));
        }
        // `GetBoundaryElement(1)`'s condition.
        let nb = ctx.mesh.n_bulk();
        let bc = ctx.bc_id(ElemId(nb as u32)).ok_or_else(|| SolverError::fatal(CALLER, "the first boundary element has no condition"))?;
        let name = ctx.variable_name().unwrap_or("").to_string();
        let id = ctx.sim.vars.id(&name).ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
        // The procedure-valued condition replaces the constant one.
        ctx.sim.model.boundary_conditions[bc - 1].0.remove(&name);
        let nn = ctx.mesh.n_nodes();
        let mut new = [vec![0.0; nn], vec![0.0; nn]];
        let mut geometry = [Vec::new(), Vec::new()];
        for c in 0..2 {
            assemble(ctx, false, None)?;
            let fixed = procedure_dirichlet(ctx, bc, c)?;
            ctx.default_dirichlet_after(&fixed)?;
            ctx.default_solve()?;
            let v = ctx.sim.vars.by_id(id);
            for k in 0..nn {
                new[c][k] = v.slot(k as u32).map_or(0.0, |s| v.values[s]);
            }
            // `Mesh % Nodes % x = Values(Perm)` for every dof.
            geometry[c] = v.values.clone();
        }
        let coords: Vec<[f64; 3]> = (0..nn).map(|k| [new[0][k], new[1][k], ctx.mesh.coord(k as u32)[2]]).collect();
        if !elmer_physics_meshtools::motion::set_coords(&coords) {
            return Err(SolverError::NotWired("node coordinates the solver can move".into()));
        }
        ctx.sim.model.boundary_conditions[bc - 1].insert(&name, Value::real(0.0));
        assemble(ctx, true, Some(&geometry))?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}

/// `CircX(Model, x, y, z)`: `sqrt(2) x / r - x`.
pub fn circ_x(p: [f64; 3]) -> f64 {
    let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
    2.0f64.sqrt() * p[0] / r - p[0]
}

/// `CircY(Model, x, y, z)`: `sqrt(2) y / r - y`.
pub fn circ_y(p: [f64; 3]) -> f64 {
    let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
    2.0f64.sqrt() * p[1] / r - p[1]
}

pub fn factory2(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(CurvedPoisson2))
}

/// `CurvedBndryPFEM2/Poisson.f90`: `LocalMatrix` with `GetElementNodes`
/// (the nodes, then the p dofs' coordinates a p MeshSolve has set), the
/// default Dirichlet conditions and solve, then the energy printed with
/// its error against pi/2.
struct CurvedPoisson2;

impl Solver for CurvedPoisson2 {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::Unsupported("the CurvedBndryPFEM2 Poisson plugin without p elements".into()));
        }
        let geometry = ctx.p_geometry();
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut slots = Vec::new();
        let mut ev = BasisEval::default();
        let mut locals = Vec::with_capacity(active.len());
        for &e in &active {
            let d = ctx.solver_dofs().expect("checked");
            d.element_slots(e, &mut slots);
            let nd = slots.len();
            let nodes = d.element_nodes(e);
            let basis = d.element_basis_scalar(e).map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
            let n = nodes.len();
            let etype = ctx.mesh.element(e).etype;
            let mut xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            xyz.extend(slots[n..].iter().map(|s| geometry.as_ref().and_then(|g| g.get(s).copied()).unwrap_or([0.0; 3])));
            let load = match ctx.body_params(e).and_then(|b| b.body_force) {
                Some(bf) => {
                    let el_nodes = ctx.mesh.element(e).nodes.to_vec();
                    let v = ctx.real_at_nodes(bf, "source", e)?.unwrap_or_else(|| vec![0.0; n]);
                    nodes.iter().map(|k| el_nodes.iter().position(|m| m == k).map_or(0.0, |i| v[i])).collect()
                }
                None => vec![0.0; n],
            };
            let mut local = LocalSystem::zeros(nd);
            for qp in &rule {
                basis.eval(etype, &xyz, qp.uvw, &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let w = qp.weight * ev.det_j;
                let mut load_at_ip = 0.0;
                for i in 0..n {
                    load_at_ip += ev.basis[i] * load[i];
                }
                for i in 0..nd {
                    let gi = ev.dbasis_dx[i];
                    for j in 0..nd {
                        let gj = ev.dbasis_dx[j];
                        local.stiff[i * nd + j] += w * (gi[0] * gj[0] + gi[1] * gj[1] + gi[2] * gj[2]);
                    }
                    local.force[i] += w * load_at_ip * ev.basis[i];
                }
            }
            ctx.default_update(&local, e)?;
            locals.push((e, local.stiff.clone(), slots.clone()));
        }
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        // Energy = SUM(POT * MATMUL(STIFF, POT)).
        let name = ctx.variable_name().unwrap_or("").to_string();
        let v = ctx.sim.vars.get(&name).ok_or_else(|| SolverError::fatal(CALLER, "no variable"))?;
        let mut energy = 0.0;
        for (_, stiff, slots) in &locals {
            let nd = slots.len();
            let pot: Vec<f64> = slots.iter().map(|&s| v.values[s as usize]).collect();
            let mut sum = 0.0;
            for i in 0..nd {
                let mut row = 0.0;
                for j in 0..nd {
                    row += stiff[i * nd + j] * pot[j];
                }
                sum += pot[i] * row;
            }
            energy += sum;
        }
        let rows = v.values.len();
        println!("{}", crate::udf::list_real(energy));
        println!(" DOFs, Error in energy: {:>12}{}", rows, crate::udf::list_real((energy - std::f64::consts::PI / 2.0).abs()));
        Ok(())
    }
}
