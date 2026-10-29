//! `RigidMeshMapper` (RigidMeshMapper.F90:46): moves the nodes of the bodies
//! by analytic rotation, scaling and translation read from the `Body Force`
//! (`Mesh Rotate`, `Mesh Scale`, `Mesh Translate`, `Mesh Origin`, `Mesh
//! Matrix`), optionally weighted by `Mesh Relax` or by a relaxation field
//! solved from a Laplace equation on the solver's own variable.
//!
//! Not ported, and refused: `Rotor Mode` (`SetMeshSkew`), `Calculate Mesh
//! Velocity`, `Target Mesh Solver Index`, `mmg remesh` and `Mesh Youngs
//! Modulus`.

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, Value};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::util::{apply_coords, const_array, coords_of, has_prefix, node_real, slot_index};

const CALLER: &str = "RigidMeshMapper";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(RigidMeshMapper::default()))
}

type Mat4 = [[f64; 4]; 4];

const IDENTITY: Mat4 = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
];

/// `MATMUL(a, b)`.
fn matmul(a: &Mat4, b: &Mat4) -> Mat4 {
    let mut c = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            c[i][j] = (0..4).fold(0.0, |s, k| s + a[i][k] * b[k][j]);
        }
    }
    c
}

/// The state `RigidMeshMapper` keeps between nodes (and `SAVE`s between
/// visits): the transformation read for the previous node, reused when
/// `Constant Mapping` says the value list did not change.
struct Map {
    trf: Mat4,
    origin: [f64; 4],
    d_coord: [f64; 3],
    got_matrix: bool,
    got_translate: bool,
}

#[derive(Default)]
pub struct RigidMeshMapper {
    visited: bool,
    /// `Mesh % NodesOrig`.
    orig: Option<Vec<[f64; 3]>>,
    /// The relaxation field (solver variable values and perm), once solved.
    relax: Option<(Vec<f64>, Vec<u32>)>,
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

impl RigidMeshMapper {
    fn refuse(ctx: &SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        if p.integer("target mesh solver index").is_some() {
            return Err(unsupported("Target Mesh Solver Index"));
        }
        if p.logical("rotor mode") == Some(true) {
            return Err(unsupported("Rotor Mode"));
        }
        if p.logical("calculate mesh velocity") == Some(true) {
            return Err(unsupported("Calculate Mesh Velocity"));
        }
        if p.logical("mmg remesh") == Some(true) {
            return Err(unsupported("mmg remesh"));
        }
        if ctx.sim.model.materials.iter().any(|m| m.get("mesh youngs modulus").is_some()) {
            return Err(unsupported("Mesh Youngs Modulus"));
        }
        Ok(())
    }

    /// The Laplace problem for the relaxation field (RigidMeshMapper.F90:259-326).
    fn solve_relaxation(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        ctx.info(CALLER, "Solving mesh relaxation field: ", 5);
        let name = ctx
            .variable_name()
            .ok_or_else(|| SolverError::fatal(CALLER, "no variable for the relaxation field"))?
            .to_string();
        // Moving and fixed BCs / bodies.
        let model = &mut ctx.sim.model;
        for (list, moving, fixed) in [
            (&mut model.boundary_conditions, "moving boundary", "fixed boundary"),
            (&mut model.body_forces, "moving body", "fixed body"),
        ] {
            for sec in list.iter_mut() {
                if sec.logical(moving) == Some(true) {
                    sec.insert(&name, Value::real(1.0));
                } else if sec.logical(fixed) == Some(true) {
                    sec.insert(&name, Value::real(0.0));
                }
            }
        }
        let max_iter = ctx.params().integer("nonlinear system max iterations").unwrap_or(1);
        let coeff = ctx.creal(ctx.params(), "nonlinear conductivity coefficient")?.unwrap_or(0.0);
        let source = ctx.creal(ctx.params(), "mesh relax source")?.unwrap_or(0.0);
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut local = LocalSystem::zeros(0);
        for _ in 0..max_iter.max(1) {
            ctx.default_initialize()?;
            for &e in &active {
                laplace_local(ctx, e, coeff, source, &mut local)?;
                ctx.default_update(&local, e)?;
            }
            ctx.default_finish_bulk_assembly()?;
            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            if ctx.default_converged() {
                break;
            }
        }
        if ctx.params().logical("mesh relax normalize") == Some(true) {
            let f = ctx.primary_mut().expect("relaxation variable");
            let max = f.values.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            ctx.info(CALLER, &format!("Normalizing deformation by:{max:12.3E}"), 6);
            let f = ctx.primary_mut().expect("relaxation variable");
            f.values.iter_mut().for_each(|v| *v /= max);
        }
        if let Some(f) = ctx.primary() {
            if !f.values.is_empty() {
                self.relax = Some((f.values.to_vec(), f.perm.map(<[u32]>::to_vec).unwrap_or_default()));
            }
        }
        Ok(())
    }
}

/// `LocalMatrix` (RigidMeshMapper.F90:620): `-div((1 + c |grad r|) grad r) = s`.
fn laplace_local(
    ctx: &SolverCtx<'_>,
    e: ElemId,
    coeff: f64,
    source: f64,
    out: &mut LocalSystem,
) -> Result<(), SolverError> {
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("element basis and quadrature (elmer-fe)".into()))?;
    let el = ctx.mesh.element(e);
    let n = el.nodes.len();
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&i| ctx.mesh.coord(i)).collect();
    // GetScalarLocalSolution.
    let local_relax: Vec<f64> = {
        let f = ctx.primary().expect("relaxation variable");
        el.nodes
            .iter()
            .map(|&i| match f.perm.map(|p| p[i as usize]) {
                Some(s) if s != elmer_core::NO_DOF => f.values[s as usize],
                _ => 0.0,
            })
            .collect()
    };
    *out = LocalSystem::zeros(n);
    let mut be = BasisEval::default();
    for ip in fe.quadrature.gauss_points(el.etype, None) {
        fe.basis
            .eval(el.etype, &xyz, ip.uvw, &mut be)
            .map_err(|m| SolverError::fatal(CALLER, m))?;
        let mut grad = [0.0; 3];
        for (i, g) in grad.iter_mut().enumerate() {
            *g = (0..n).fold(0.0, |s, p| s + be.dbasis_dx[p][i] * local_relax[p]);
        }
        let cond = 1.0 + coeff * (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
        let s = ip.weight * be.det_j;
        for i in 0..3 {
            for p in 0..n {
                for q in 0..n {
                    out.stiff[p * n + q] += cond * s * be.dbasis_dx[p][i] * be.dbasis_dx[q][i];
                }
            }
        }
        for p in 0..n {
            out.force[p] += source * s * be.basis[p];
        }
    }
    Ok(())
}

impl Solver for RigidMeshMapper {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "---------------------------------------", 4);
        ctx.info(CALLER, "Performing analytic mesh mapping ", 4);
        ctx.info(CALLER, "---------------------------------------", 4);
        Self::refuse(ctx)?;

        let params = ctx.params().clone();
        let cumulative = params.logical("cumulative displacements").unwrap_or(false);
        let use_original = !cumulative;
        let store_original = params.logical("store original coordinates").unwrap_or(false);
        let angle_coeff = if ctx.simulation_params().logical("rotate in radians") == Some(true) {
            1.0
        } else {
            std::f64::consts::PI / 180.0
        };
        let whole_mode = params.logical("whole mesh mode").unwrap_or(false);
        if whole_mode {
            ctx.info(CALLER, "Moving the whole mesh with keywords from Solver section!", 4);
        }
        let rotor_angle = ctx.creal(&ctx.sim.model.simulation, "rotor angle")?;
        let dim = ctx.mesh.dim();

        let current = coords_of(ctx.mesh);
        if !self.visited && (use_original || store_original) {
            ctx.info(CALLER, "Storing original coordinates", 7);
            self.orig = Some(current.clone());
        }
        let orig = if use_original {
            self.orig.clone().unwrap_or_else(|| current.clone())
        } else {
            current.clone()
        };

        let translate_before_rotate = params.logical("translate before rotate").unwrap_or(false);
        // Reversed order of the rotation axes.
        let mut rotate_order = [1usize, 2, 3];
        if let Some(Value::Integer(ia)) = params.get("mesh rotation axis order") {
            if ia.len() != dim {
                return Err(SolverError::fatal(CALLER, "Size of Mesh Rotation Axis Order must match dimension."));
            }
            for i in 0..ia.len() {
                rotate_order[i] = ia[ia.len() - 1 - i] as usize;
            }
        }

        // The relaxation field is solved on the first visit, when the solver
        // has a matrix (a variable).
        let has_matrix = ctx.sim.slots[slot_index(ctx)].variable.is_some_and(|id| {
            ctx.sim.vars.by_id(id).kind == elmer_core::FieldKind::Nodal
        });
        if has_matrix && !self.visited {
            self.solve_relaxation(ctx)?;
        }

        let model = &ctx.sim.model;
        let any = |f: &dyn Fn(&ParamMap) -> bool| -> bool {
            if whole_mode {
                f(&params)
            } else {
                model.body_forces.iter().any(f)
            }
        };
        let any_matrix = any(&|p| p.get("mesh matrix").is_some());
        let any_rotate = any(&|p| has_prefix(p, "mesh rotate"));
        let any_translate = any(&|p| has_prefix(p, "mesh translate"));
        let any_scale = any(&|p| has_prefix(p, "mesh scale"));
        let any_origin = any(&|p| has_prefix(p, "mesh origin"));
        let any_relax = any(&|p| p.get("mesh relax").is_some());
        if model.body_forces.iter().any(|p| has_prefix(p, "mesh displacement")) {
            ctx.info(CALLER, "> Mesh Displacement < is an obsolete keyword", 3);
            ctx.warn(CALLER, "Replace with > Mesh Translate < ");
        }
        let constant_map = params.logical("constant mapping").unwrap_or(false);

        let mut new = current.clone();
        let mut done = vec![false; current.len()];
        let mut map = Map { trf: IDENTITY, origin: [0.0; 4], d_coord: [0.0; 3], got_matrix: false, got_translate: false };
        let mut prev_list: Option<*const ParamMap> = None;
        let mut got_map = false;
        let mut sec_ptr: &ParamMap = &params;

        for t in 0..ctx.mesh.n_bulk() as u32 {
            let e = ElemId(t);
            if !whole_mode {
                let Some(bf) = ctx.body_params(e).and_then(|b| b.body_force) else { continue };
                sec_ptr = bf;
                if constant_map {
                    got_map = prev_list == Some(bf as *const ParamMap);
                    prev_list = Some(bf as *const ParamMap);
                }
            }
            let section = sec_ptr;
            for &node in ctx.mesh.element(e).nodes {
                let ni = node as usize;
                if done[ni] {
                    continue;
                }
                if !got_map {
                    if whole_mode {
                        got_map = true;
                    }
                    self.read_map(
                        ctx, section, node, &mut map,
                        &Flags { any_matrix, any_rotate, any_translate, any_scale, any_origin },
                        angle_coeff, rotor_angle, translate_before_rotate, rotate_order,
                    )?;
                }
                let x0 = [orig[ni][0], orig[ni][1], orig[ni][2], 1.0];
                let mut dx = if map.got_matrix {
                    let v: [f64; 4] = std::array::from_fn(|i| x0[i] - map.origin[i]);
                    let x1: [f64; 4] = std::array::from_fn(|i| {
                        (0..4).fold(0.0, |s, k| s + map.trf[i][k] * v[k]) + map.origin[i]
                    });
                    [x1[0] / x1[3] - x0[0], x1[1] / x1[3] - x0[1], x1[2] / x1[3] - x0[2]]
                } else if map.got_translate {
                    map.d_coord
                } else {
                    continue;
                };
                // Interpolate the displacement between moving and fixed walls.
                if let Some((vals, perm)) = &self.relax {
                    let s = perm[ni];
                    if s != elmer_core::NO_DOF {
                        let r = vals[s as usize];
                        dx.iter_mut().for_each(|d| *d *= r);
                    }
                } else if any_relax {
                    if let Some(r) = node_real(ctx, section, "mesh relax", node)? {
                        dx.iter_mut().for_each(|d| *d *= r);
                    }
                }
                new[ni] = [dx[0] + orig[ni][0], dx[1] + orig[ni][1], dx[2] + orig[ni][2]];
                done[ni] = true;
            }
        }

        if !self.visited {
            ctx.info(CALLER, &format!("Number of nodes mapped: {}", done.iter().filter(|&&d| d).count()), 4);
            ctx.info(CALLER, "All done", 4);
        }
        apply_coords(CALLER, &current, &new)?;
        self.visited = true;
        Ok(())
    }
}

struct Flags {
    any_matrix: bool,
    any_rotate: bool,
    any_translate: bool,
    any_scale: bool,
    any_origin: bool,
}

impl RigidMeshMapper {
    /// The per-node reading of the transformation (RigidMeshMapper.F90:402-513).
    #[allow(clippy::too_many_arguments)]
    fn read_map(
        &self,
        ctx: &SolverCtx<'_>,
        vl: &ParamMap,
        node: u32,
        map: &mut Map,
        f: &Flags,
        angle_coeff: f64,
        rotor_angle: Option<f64>,
        translate_before_rotate: bool,
        rotate_order: [usize; 3],
    ) -> Result<(), SolverError> {
        map.got_matrix = false;
        if f.any_matrix {
            if let Some(Value::Real { rows, cols, data }) = vl.get("mesh matrix") {
                map.got_matrix = true;
                // TrfMatrix(i,j) = Parray(j,i); entries the array lacks keep their old value.
                for i in 0..*rows {
                    for j in 0..*cols {
                        map.trf[i][j] = data[j * cols + i];
                    }
                }
            }
        }
        if !map.got_matrix {
            let mut got_rotate = false;
            let mut angles = [0.0f64; 3];
            if f.any_rotate {
                if let Some(a) = const_array(vl, "mesh rotate") {
                    got_rotate = true;
                    for (i, v) in a.iter().take(3).enumerate() {
                        angles[i] = *v;
                    }
                } else {
                    for (k, key) in ["mesh rotate 1", "mesh rotate 2", "mesh rotate 3"].into_iter().enumerate() {
                        let r = node_real(ctx, vl, key, node)?;
                        angles[k] = r.unwrap_or(0.0);
                        got_rotate |= r.is_some();
                    }
                }
                angles.iter_mut().for_each(|a| *a *= angle_coeff);
            }
            if let Some(ra) = rotor_angle {
                got_rotate = true;
                angles[2] = angle_coeff * ra;
            }

            let mut got_scale = false;
            let mut scaling = [0.0f64; 3];
            if f.any_scale {
                if let Some(a) = const_array(vl, "mesh scale") {
                    got_scale = true;
                    for (i, v) in a.iter().take(3).enumerate() {
                        scaling[i] = *v;
                    }
                } else {
                    for (k, key) in ["mesh scale 1", "mesh scale 2", "mesh scale 3"].into_iter().enumerate() {
                        let r = node_real(ctx, vl, key, node)?;
                        scaling[k] = r.unwrap_or(1.0);
                        got_scale |= r.is_some();
                    }
                }
            }

            map.got_translate = false;
            if f.any_translate {
                if let Some(a) = const_array(vl, "mesh translate") {
                    map.got_translate = true;
                    map.d_coord = [0.0; 3];
                    for (i, v) in a.iter().take(3).enumerate() {
                        map.d_coord[i] = *v;
                    }
                } else {
                    for (k, key) in ["mesh translate 1", "mesh translate 2", "mesh translate 3"].into_iter().enumerate() {
                        let r = node_real(ctx, vl, key, node)?;
                        map.d_coord[k] = r.unwrap_or(0.0);
                        map.got_translate |= r.is_some();
                    }
                }
            }

            map.got_matrix = got_rotate || got_scale;
            if map.got_matrix {
                let mut trs = IDENTITY;
                let mut rot = IDENTITY;
                if got_rotate {
                    for &j in &rotate_order {
                        let alpha = angles[j - 1];
                        if alpha.abs() < f64::MIN_POSITIVE {
                            continue;
                        }
                        let mut t = IDENTITY;
                        let (c, s) = (alpha.cos(), alpha.sin());
                        match j {
                            1 => {
                                t[1][1] = c;
                                t[1][2] = -s;
                                t[2][1] = s;
                                t[2][2] = c;
                            }
                            2 => {
                                t[0][0] = c;
                                t[0][2] = -s;
                                t[2][0] = s;
                                t[2][2] = c;
                            }
                            _ => {
                                t[0][0] = c;
                                t[0][1] = -s;
                                t[1][0] = s;
                                t[1][1] = c;
                            }
                        }
                        rot = matmul(&rot, &t);
                    }
                }
                if map.got_translate {
                    for i in 0..3 {
                        trs[i][3] = map.d_coord[i];
                    }
                }
                map.trf = if translate_before_rotate { matmul(&rot, &trs) } else { matmul(&trs, &rot) };
                if got_scale {
                    let mut scl = IDENTITY;
                    for i in 0..3 {
                        scl[i][i] = scaling[i];
                    }
                    map.trf = matmul(&scl, &map.trf);
                }
            }
        }

        // Mesh origin.
        map.origin = [0.0; 4];
        if map.got_matrix && f.any_origin {
            if let Some(a) = const_array(vl, "mesh origin") {
                for (i, v) in a.iter().take(4).enumerate() {
                    map.origin[i] = *v;
                }
            } else {
                for (k, key) in ["mesh origin 1", "mesh origin 2", "mesh origin 3"].into_iter().enumerate() {
                    map.origin[k] = node_real(ctx, vl, key, node)?.unwrap_or(0.0);
                }
            }
        }
        Ok(())
    }
}
