//! `StatElecForce` (fem/src/modules/ElectricForce.F90): the electric force
//! on the boundaries with `Calculate Electric Force`, by the Maxwell stress
//! tensor `T = (D.E / 2) I - D E^T` of the potential's gradient taken in the
//! boundary element's parent, integrated against the parent-oriented normal
//! (`NormalVector(.., Check = .TRUE.)`). The net force (and, given `Moment
//! About i`, the moment) goes to the Simulation as `res: Electric Force i`
//! (`res: Electric Moment i`); the nodal force density, the element forces
//! lumped over the basis divided by the lumped boundary measure, to the
//! vector field `Electric Force Density` on the potential's numbering.
//!
//! The permittivity is the `Relative Permittivity` (else `Permittivity`) of
//! the left parent's material, else of the right one's, at the boundary
//! element's nodes. Refused: `Calculate Nodal Force` (`Electric Nodal
//! Force`), `All Displaced Boundaries`, tensor permittivities and
//! coordinate systems other than Cartesian.

#![allow(clippy::needless_range_loop)]

use elmer_core::{BasisEval, ElemId, ElementBasis, FieldKind, MeshView, Params, Value};
use elmer_fe::{CoordinateKind, CoordinateSystem, LagrangeBasis};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable};

const CALLER: &str = "StatElecForce";
const PROCEDURE: (&str, &str) = ("ElectricForce", "StatElecForce");
const FIELD: &str = "electric force density";

pub(crate) fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| Ok(Box::new(StatElecForce)));
}

struct StatElecForce;

fn unsupported(m: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {}", m.into()))
}

/// `Element % TYPE % NodeU/V/W` of node `k`.
fn node_uvw(code: u16, k: usize) -> Result<[f64; 3], SolverError> {
    let d = elmer_fe::element_type(code).map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
    let at = |v: &Vec<f64>| v.get(k).copied().unwrap_or(0.0);
    Ok([at(&d.node_u), at(&d.node_v), at(&d.node_w)])
}

fn coords(mesh: &dyn MeshView, nodes: &[u32]) -> Vec<[f64; 3]> {
    nodes.iter().map(|&k| mesh.coord(k)).collect()
}

/// The permittivity of `parent`'s material at `e`'s nodes, if given.
fn permittivity(ctx: &SolverCtx<'_>, parent: ElemId, e: ElemId) -> Result<Option<Vec<f64>>, SolverError> {
    let Some(mat) = ctx.body_params(parent).and_then(|b| b.material) else { return Ok(None) };
    for key in ["relative permittivity", "permittivity"] {
        match mat.get(key) {
            None => continue,
            Some(Value::Real { data, .. }) if data.len() > 1 => return Err(unsupported("a tensor permittivity")),
            Some(_) => return ctx.real_at_nodes(mat, key, e),
        }
    }
    Ok(None)
}

impl Solver for StatElecForce {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let params = ctx.params().clone();
        let csys = ctx.simulation_params().string("coordinate system");
        let cs = CoordinateSystem::from_keyword(csys, ctx.mesh.dim())
            .ok_or_else(|| SolverError::fatal("SetCoordinateSystem", "Unknown Coordinate System"))?;
        if cs.kind != CoordinateKind::Cartesian {
            return Err(unsupported("coordinate systems other than Cartesian"));
        }
        let dim = cs.dim;
        if params.logical("calculate nodal force") == Some(true) || ctx.sim.vars.get("electric nodal force").is_some() {
            return Err(unsupported("Calculate Nodal Force (Electric Nodal Force)"));
        }
        if params.logical("all displaced boundaries") == Some(true) {
            return Err(unsupported("All Displaced Boundaries"));
        }
        let pot_name = params.string("potential field name").unwrap_or("potential").to_string();
        let pot = ctx
            .sim
            .vars
            .get(&pot_name)
            .ok_or_else(|| SolverError::fatal("ElectricForce", format!("Potential field does not exist:{pot_name}")))?;
        let (pot_perm, pot_values) = (pot.perm.clone(), pot.values.clone());
        let pot_slot = |k: u32| -> Option<usize> {
            match &pot_perm {
                Some(p) => p.get(k as usize).copied().filter(|&s| s != elmer_core::NO_DOF).map(|s| s as usize),
                None => Some(k as usize),
            }
        };
        // `VariableAddVector(.., 'Electric Force Density', Dim, Perm =
        // PotentialPerm)` unless it exists.
        if ctx.sim.vars.get(FIELD).is_none() {
            let n = pot_values.len();
            let var = Variable::new(FIELD, FieldKind::Nodal, dim, pot_perm.clone(), vec![0.0; n * dim]);
            let comps: Vec<String> = (1..=dim).map(|i| format!("{FIELD} {i}")).collect();
            ctx.sim.vars.add_with_components(var, &comps)?;
        }
        let (fdofs, flen) = {
            let f = ctx.sim.vars.get(FIELD).expect("added above");
            (f.dofs, f.values.len())
        };
        if fdofs != dim {
            return Err(unsupported(format!("Electric Force Density of {fdofs} components in {dim}D")));
        }

        let eps0 = match ctx.constants().real("permittivity of vacuum") {
            Some(v) => v,
            None => {
                ctx.warn(CALLER, "Permittivity of Vacuum not given, using 1");
                1.0
            }
        };
        let sim = ctx.simulation_params();
        let about: [Option<f64>; 3] = std::array::from_fn(|i| sim.real(&format!("moment about {}", i + 1)));
        let calc_moment = about.iter().any(Option::is_some);
        let about = about.map(|a| a.unwrap_or(0.0));

        let basis = LagrangeBasis::new(dim);
        let (mut ev, mut pev) = (BasisEval::default(), BasisEval::default());
        let mut force = [0.0f64; 3];
        let mut moment = [0.0f64; 3];
        let mut density = vec![0.0f64; flen];
        let mut weight = vec![0.0f64; ctx.mesh.n_nodes().max(flen / dim.max(1))];
        let nb = ctx.mesh.n_bulk();
        for t in nb..nb + ctx.mesh.n_boundary() {
            let e = ElemId(t as u32);
            let el = ctx.mesh.element(e);
            if el.etype.0 == 101 {
                continue;
            }
            if !ctx.bc_params(e).is_some_and(|bc| bc.logical("calculate electric force") == Some(true)) {
                continue;
            }
            let nodes = el.nodes.to_vec();
            let n = nodes.len();
            let [left, right] = el.parents;
            let mut found = None;
            if let Some(p) = left {
                found = permittivity(ctx, p, e)?.map(|eps| (p, eps));
            }
            if found.is_none() {
                let Some(p) = right else {
                    return Err(SolverError::fatal(CALLER, "No permittivity found on specified boundary"));
                };
                let eps = permittivity(ctx, p, e)?
                    .ok_or_else(|| SolverError::fatal(CALLER, "No permittivity found on specified boundary"))?;
                found = Some((p, eps));
            }
            let (parent, eps) = found.expect("set above");
            if nodes.iter().any(|&k| pot_slot(k).is_none()) {
                return Err(SolverError::fatal(CALLER, "No potential available for specified boundary"));
            }
            let pel = ctx.mesh.element(parent);
            let pnodes = pel.nodes.to_vec();
            let pn = pnodes.len();
            let local_pot: Vec<f64> = pnodes.iter().map(|&k| pot_slot(k).map_or(0.0, |s| pot_values[s])).collect();
            let xyz = coords(ctx.mesh, &nodes);
            let pxyz = coords(ctx.mesh, &pnodes);
            // The boundary nodes' reference coordinates in the parent.
            let mut loc = vec![[0.0f64; 3]; n];
            for (i, &k) in nodes.iter().enumerate() {
                if let Some(j) = pnodes.iter().position(|&q| q == k) {
                    loc[i] = node_uvw(pel.etype.0, j)?;
                }
            }
            let def = elmer_fe::element_type(el.etype.0).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
            let rule = elmer_fe::quadrature::gauss_points(el.etype, Some(def.gauss_points))
                .map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
            let mut el_weight = vec![0.0f64; n];
            let mut el_force = vec![[0.0f64; 3]; n];
            for l in 0..rule.len() {
                let uvw = [rule.u[l], rule.v[l], rule.w[l]];
                basis.eval(el.etype, &xyz, uvw, &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let b = ev.basis[..n].to_vec();
                let s = ev.det_j * rule.s[l];
                let normal = elmer_solver::boundary_normal(ctx.mesh, &ctx.sim.model, e, [uvw[0], uvw[1]])?;
                let puvw: [f64; 3] = std::array::from_fn(|c| (0..n).map(|i| b[i] * loc[i][c]).sum());
                basis.eval(pel.etype, &pxyz, puvw, &mut pev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let pdb = &pev.dbasis_dx[..pn];
                let mut efield = [0.0f64; 3];
                let mut dflux = [0.0f64; 3];
                let eps_ip: f64 = (0..n).map(|i| eps[i] * b[i]).sum();
                for i in 0..dim {
                    efield[i] = (0..pn).map(|k| pdb[k][i] * local_pot[k]).sum();
                    for j in 0..dim {
                        // Permittivity(i,j) is diagonal: zero off it.
                        let pij = if i == j { eps_ip } else { 0.0 };
                        dflux[i] += (0..pn).map(|k| pdb[k][j] * local_pot[k]).sum::<f64>() * pij;
                    }
                }
                let dflux = dflux.map(|d| eps0 * d);
                let mut tensor = [[0.0f64; 3]; 3];
                for i in 0..dim {
                    for j in 0..dim {
                        tensor[i][j] = -dflux[i] * efield[j];
                    }
                }
                let de: f64 = (0..3).map(|k| dflux[k] * efield[k]).sum();
                for (i, row) in tensor.iter_mut().enumerate().take(dim) {
                    row[i] += de / 2.0;
                }
                let lforce: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| tensor[i][j] * normal[j]).sum());
                for i in 0..3 {
                    force[i] += s * lforce[i];
                }
                if calc_moment {
                    let r: [f64; 3] = std::array::from_fn(|c| (0..n).map(|i| xyz[i][c] * b[i]).sum::<f64>() - about[c]);
                    let lm = [
                        r[1] * lforce[2] - r[2] * lforce[1],
                        r[2] * lforce[0] - r[0] * lforce[2],
                        r[0] * lforce[1] - r[1] * lforce[0],
                    ];
                    for i in 0..3 {
                        moment[i] += s * lm[i];
                    }
                }
                for i in 0..n {
                    for c in 0..3 {
                        el_force[i][c] += s * lforce[c] * b[i];
                    }
                    el_weight[i] += s * b[i];
                }
            }
            for (i, &node) in nodes.iter().enumerate() {
                let k = pot_slot(node).expect("checked above");
                for j in 0..dim {
                    if let Some(d) = density.get_mut(dim * k + j) {
                        *d += el_force[i][j];
                    }
                }
                weight[k] += el_weight[i];
            }
        }
        // The force density on the potential's numbering.
        for i in 0..ctx.mesh.n_nodes() {
            if weight.get(i).is_some_and(|&w| w > 0.0) {
                for j in 0..dim {
                    if let Some(d) = density.get_mut(dim * i + j) {
                        *d /= weight[i];
                    }
                }
            }
        }
        ctx.sim.vars.get_mut(FIELD).expect("added above").values = density;
        let fmt = |v: f64| elmer_solver::fortran_fmt::es(v, 15, 6);
        let shown: String = force[..if dim == 2 { 2 } else { 3 }].iter().map(|&f| fmt(f)).collect();
        ctx.info(CALLER, &format!("Net electric force  : {shown}"), 4);
        let sf = (force[0] * force[0] + force[1] * force[1] + force[2] * force[2]).sqrt();
        ctx.info(CALLER, &format!("Resultant force  : {}", fmt(sf)), 4);
        let model_sim = &mut ctx.sim.model.simulation;
        if calc_moment {
            for i in [3, 2, 1] {
                elmer_solver::results::add_result(model_sim, &format!("res: Electric Moment {i}"), moment[i - 1]);
            }
        }
        if dim > 2 {
            elmer_solver::results::add_result(model_sim, "res: Electric Force 3", force[2]);
        }
        elmer_solver::results::add_result(model_sim, "res: Electric Force 2", force[1]);
        elmer_solver::results::add_result(model_sim, "res: Electric Force 1", force[0]);
        Ok(())
    }
}
