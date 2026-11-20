//! `Radiation = Diffuse Gray`: the radiating surfaces, their Gebhart
//! factors (`RadiationFactors` of RadiationFactors.F90 with the `ViewFactors`
//! prestep done in process by `elmer-viewfactors`) and the boundary
//! assembly of `LocalMatrixDiffuseGray` (HeatSolve.F90:2774).
//!
//! Supported: axisymmetric two-node boundary elements, explicit Gebhart
//! factors, open boundaries, the Newton and the Picard linearisation.
//! Factors are computed once, on the first solve, from the emissivities
//! and geometry of that moment (the Fortran's `FirstTime` call).

use std::collections::HashMap;

use elmer_core::{BasisEval, CsrMatrix, ElemId, ElementBasis, Params, Quadrature, NO_DOF};
use elmer_solver::{SolverCtx, SolverError};
use elmer_viewfactors::{
    gebhart_factors_with, line_direction, normalize_factors_with, GebhartMatrix, saved_rows, view_factors_axis,
    GebhartInput, ParentInfo,
};

const CALLER: &str = "RadiationFactors";

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(what.to_string())
}

/// `ElementArea` (ElementUtils.F90:3346) of a two-node boundary element in
/// axisymmetric coordinates: `SUM r detJ s` over `GaussPoints(Element)`
/// (its analytic branch is compiled out).
fn element_area_axis(a: [f64; 3], b: [f64; 3]) -> Result<f64, SolverError> {
    let fe = |e: elmer_fe::FeError| SolverError::fatal("ElementArea", e.to_string());
    let def = elmer_fe::element_type(202).map_err(fe)?;
    let ip = elmer_fe::gauss_points(elmer_core::ElementType(202), None).map_err(fe)?;
    let nodes = [[a[0], a[1], 0.0], [b[0], b[1], 0.0]];
    let mut st = elmer_fe::MetricStatus::default();
    let mut ev = BasisEval::default();
    let mut area = 0.0;
    for t in 0..ip.s.len() {
        elmer_fe::element_info(def, &nodes, [ip.u[t], ip.v[t], ip.w[t]], 2, &mut st, &mut ev).map_err(fe)?;
        let x = nodes[0][0] * ev.basis[0] + nodes[1][0] * ev.basis[1];
        area += x * ev.det_j * ip.s[t];
    }
    Ok(area)
}

/// `fm_MatVec`: the dense row-major Jacobian of `NormalizeFactors`.
struct DenseOp<'a> {
    n: usize,
    a: &'a [f64],
}

impl elmer_core::LinearOperator<f64> for DenseOp<'_> {
    fn nrows(&self) -> usize {
        self.n
    }
    fn ncols(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..self.n {
            let mut s = 0.0;
            for j in 0..self.n {
                s += self.a[i * self.n + j] * x[j];
            }
            y[i] = s;
        }
    }
}

/// `fm_DiagPrec`: `u = v * Jdiag`.
struct DiagPrec<'a>(&'a [f64]);

impl elmer_core::Preconditioner<f64> for DiagPrec<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for i in 0..r.len() {
            z[i] = r[i] * self.0[i];
        }
    }
}

/// ViewFactors' `IterSolv`: HUTI CG (the port), `HUTI_TRESID_SCALED_BYB` to
/// 1e-6 in at most 100 iterations from the given `x`, the diagonal as the
/// left preconditioner.
fn hut_itersolv(jac: &[f64], jdiag: &[f64], x: &mut [f64], b: &[f64]) {
    use elmer_core::IterativeSolver;
    use elmer_solver::elmer_krylov::{Huti, HutiConfig, HutiMethod, HutiStop};
    let mut cfg = HutiConfig::new(HutiMethod::Cg);
    cfg.max_iter = 100;
    cfg.tol = 1.0e-6;
    cfg.max_tol = 1.0e20;
    cfg.stop_crit = HutiStop::TrueResidualScaledByB;
    cfg.left_precond = true;
    cfg.abort_not_converged = false;
    let op = DenseOp { n: b.len(), a: jac };
    // A non-converged step is used as it stands, as IterSolv does.
    let _ = Huti::new(cfg).solve(&op, &DiagPrec(jdiag), b, x, &elmer_core::SerialComm);
}

/// `Radiation = Diffuse Gray` on the boundary element's condition.
pub(crate) fn is_diffuse_gray(ctx: &SolverCtx<'_>, e: ElemId) -> bool {
    ctx.bc_params(e)
        .and_then(|bc| bc.string("Radiation"))
        .is_some_and(|s| s.eq_ignore_ascii_case("diffuse gray"))
        && ctx.mesh.element(e).etype.family() > 1
}

/// The diffuse gray radiation of the model, set up on the first solve.
pub(crate) struct DiffuseGray {
    /// Boundary element -> Gebhart factors `(radiating element, factor)`.
    factors: HashMap<u32, Vec<(ElemId, f64)>>,
    /// Every radiating element, by element index: this partition's and, in
    /// a partitioned run, those imported from the others.
    radiating: Vec<ElemId>,
    topology_done: bool,
}

/// `TabulateBoundaryAverages`: element means of the emissivity, the
/// absorptivity and the reflectivity of a radiating boundary element.
fn optical(ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, f64, f64), SolverError> {
    let n = ctx.mesh.element(e).nodes.len() as f64;
    let mean = |v: Option<Vec<f64>>| v.map(|v| v.iter().sum::<f64>() / n);
    let bc = ctx.bc_params(e);
    let from =
        |sec: Option<&elmer_core::ParamMap>, key: &str| -> Result<Option<f64>, SolverError> {
            match sec {
                Some(s) => Ok(mean(ctx.real_at_nodes(s, key, e)?)),
                None => Ok(None),
            }
        };
    if let Some(emis) = from(bc, "Emissivity")? {
        let abso = from(bc, "Absorptivity")?.unwrap_or(emis);
        let refl = from(bc, "Reflectivity")?.unwrap_or(1.0 - abso);
        return Ok((emis, abso, refl));
    }
    // The material of the parent: Left first, then Right.
    for parent in ctx.mesh.element(e).parents.into_iter().flatten() {
        let mat = ctx.body_params(parent).and_then(|b| b.material);
        if let Some(emis) = from(mat, "Emissivity")? {
            let abso = from(mat, "Absorptivity")?.unwrap_or(emis);
            let refl = from(mat, "Reflectivity")?.unwrap_or(1.0 - abso);
            return Ok((emis, abso, refl));
        }
    }
    Ok((0.0, 0.0, 1.0))
}

/// `TabulateEmissivity`'s reflectivity, `1 - absorptivity - transmissivity`.
fn reflectivity_for_factors(ctx: &SolverCtx<'_>, e: ElemId) -> Result<(f64, f64), SolverError> {
    let n = ctx.mesh.element(e).nodes.len() as f64;
    let mean = |v: Option<Vec<f64>>| v.map(|v| v.iter().sum::<f64>() / n);
    let mut list = ctx.bc_params(e).filter(|b| b.get("Emissivity").is_some());
    if list.is_none() {
        for parent in ctx.mesh.element(e).parents.into_iter().flatten() {
            let mat = ctx.body_params(parent).and_then(|b| b.material);
            if mat.is_some_and(|m| m.get("Emissivity").is_some()) {
                list = mat;
                break;
            }
        }
    }
    let Some(list) = list else {
        return Err(SolverError::fatal(
            "TabulateEmissivity",
            "Emissivity list not associated!".to_string(),
        ));
    };
    let get = |key: &str| -> Result<Option<f64>, SolverError> {
        Ok(mean(ctx.real_at_nodes(list, key, e)?))
    };
    let emis = get("Emissivity")?.unwrap_or(0.0);
    let trans = get("Transmissivity")?.unwrap_or(0.0);
    let abso = get("Absorptivity")?.unwrap_or(emis);
    Ok((emis, 1.0 - abso - trans))
}

/// `GetParentInfo`: the body of a parent, the last node not on the boundary
/// element, and whether the body is radiative.
fn parent_info(
    ctx: &SolverCtx<'_>,
    parent: Option<ElemId>,
    nodes: &[u32],
) -> Option<(ParentInfo, bool)> {
    let p = parent?;
    let pe = ctx.mesh.element(p);
    let node = pe
        .nodes
        .iter()
        .rev()
        .find(|n| !nodes.contains(n))
        .copied()?;
    let body = ctx.sim.model.body_of_tag(pe.tag)?;
    let bp = ctx.sim.model.body_params(body)?;
    let rad = bp
        .body
        .logical("Radiative Body")
        .or_else(|| bp.material.and_then(|m| m.logical("Radiative Body")))
        .unwrap_or_else(|| bp.material.is_some_and(|m| m.get("Emissivity").is_some()));
    Some((
        ParentInfo {
            body: body as i32,
            node: ctx.mesh.coord(node),
        },
        rad,
    ))
}

impl DiffuseGray {
    /// Computes the Gebhart factors of every diffuse gray boundary element.
    pub(crate) fn setup(
        ctx: &SolverCtx<'_>,
        bnd: &[ElemId],
        axisym: bool,
    ) -> Result<Option<Self>, SolverError> {
        // GetMeshRadiationSurfaceInfo: in a partitioned run the surfaces
        // imported by RadiationParallelMeshDistribute too, by global element
        // id (the serial order).
        let foreign: Vec<ElemId> = ctx.foreign_boundary_elements().filter(|&e| is_diffuse_gray(ctx, e)).collect();
        let mut radiating: Vec<ElemId> = bnd
            .iter()
            .copied()
            .filter(|&e| is_diffuse_gray(ctx, e))
            .chain(foreign.iter().copied())
            .collect();
        radiating.sort_by_key(|&e| ctx.element_gid(e));
        // Only this partition's surfaces get factors (`PartIndex == myPE`);
        // every rank joins the exchange of the orientations below.
        let parallel = ctx.comm.size() > 1;
        // Every partition with surfaces must see all of them; agreed, so
        // that no rank goes on alone into the exchange below.
        let mut n_max = [radiating.len() as f64, -(radiating.len() as f64)];
        if parallel {
            ctx.comm.max_f64(&mut n_max);
        }
        let (n_hi, n_lo) = (n_max[0] as usize, (-n_max[1]) as usize);
        if n_hi == 0 {
            return Ok(None);
        }
        if parallel && n_lo != n_hi {
            return Err(unsupported(&format!(
                "diffuse gray radiation with surfaces missing on some partitions ({n_lo} to {n_hi}; this rank {}, {} imported)",
                radiating.len(),
                foreign.len()
            )));
        }
        if !axisym {
            return Err(unsupported(
                "diffuse gray radiation outside axisymmetric geometry: the Cartesian view \
                 factors (view3d) are not ported",
            ));
        }
        for key in [
            "Update View Factors",
            "Update Gebhart Factors",
            "Update Gebhardt Factors",
            "Compute View Factors",
            "Use ViewFactors As Gebhart Factors",
            "Implicit Gebhart Factor Fraction",
            "Radiation Factors Solver Full",
            "Gebhart Factors Solver Full",
            "Viewfactor Rigid Mesh Mapping",
            "Save Gebhart Factors",
        ] {
            if ctx.params().get(key).is_some() {
                return Err(unsupported(&format!("solver keyword \"{key}\"")));
            }
        }
        if ctx
            .params()
            .string("radiation: Linear System Solver")
            .is_some_and(|s| s.eq_ignore_ascii_case("iterative"))
        {
            return Err(unsupported("an iterative solver for the Gebhart factors"));
        }

        // Radiation bodies (`Radiation Boundary`), each its own cavity.
        let body_of = |e: ElemId| -> usize {
            ctx.bc_params(e)
                .and_then(|b| b.integer("Radiation Boundary"))
                .map_or(1, |k| k.max(1) as usize)
        };
        let max_body = radiating.iter().map(|&e| body_of(e)).max().unwrap_or(1);

        // RadiationFactors rewrites `mesh.nodes` as `(i7,i3,3f20.12)` for the
        // `ViewFactors` program (`Update Factors Geometry`, default on when
        // it computes the factors), which reloads the mesh from it: the
        // view factors see the coordinates rounded to 12 decimals.
        // A `ViewFactors` run before the solver (`ELMER_VIEWFACTORS_PRESTEP`,
        // which elmer-coverage sets for such a test) leaves the factors file
        // in place: RadiationFactors then neither computes nor rewrites.
        let prestep = std::env::var_os("ELMER_VIEWFACTORS_PRESTEP").is_some();
        let update_geometry = !prestep
            && ctx.params().logical("Update Factors Geometry").unwrap_or(true)
            && ctx.params().logical("Viewfactor Rigid Mesh Mapping") != Some(true);
        if update_geometry && ctx.simulation_params().get("Coordinate Scaling").is_some() {
            return Err(unsupported("Coordinate Scaling with diffuse gray radiation"));
        }
        let vf_coord = |n: u32| -> Result<[f64; 3], SolverError> {
            let x = ctx.mesh.coord(n);
            if !update_geometry {
                return Ok(x);
            }
            let mut r = [0.0; 3];
            for (d, v) in x.iter().enumerate() {
                if v.abs() >= 1.0e7 {
                    return Err(unsupported("coordinates too large for the f20.12 mesh.nodes rewrite"));
                }
                r[d] = format!("{v:.12}").parse().expect("a formatted number");
            }
            Ok(r)
        };
        // ViewFactors runs on the whole mesh: its surfaces and coordinates
        // are by global node number (the local index in a serial run), which
        // fixes the order in which it combines aligned elements.
        let n_nodes = ctx.mesh.n_nodes();
        let vf_index = |n: u32| (ctx.node_gid(n) - 1) as usize;
        let n_global = (0..n_nodes as u32).map(|n| vf_index(n) + 1).max().unwrap_or(0);
        let mut coord: Vec<f64> = vec![0.0; 2 * n_global];
        for n in 0..n_nodes as u32 {
            let x = vf_coord(n)?;
            let g = vf_index(n);
            coord[2 * g] = x[0];
            coord[2 * g + 1] = x[1];
        }
        let divide = ctx
            .params()
            .integer("Viewfactor divide")
            .map_or(1, |d| d.max(1) as usize);
        let combine = ctx
            .params()
            .logical("Viewfactor combine elements")
            .unwrap_or(ctx.mesh.dim() == 2);
        let min_vf = ctx.params().real("Minimum View Factor").unwrap_or(1.0e-20);
        let single = ctx.params().logical("Viewfactor Binary Output") == Some(true)
            && ctx.params().logical("Viewfactor single precision") == Some(true);
        let min_gebhart = ctx
            .params()
            .real("Minimum Gebhart Factor")
            .or_else(|| ctx.params().real("Minimum Gebhardt Factor"))
            .unwrap_or(1.0e-20);

        // ViewFactors' orientation of each surface from its parents. An
        // imported surface has none here: its own partition orients it, and
        // the ranks sum the orientations (each surface has one owner).
        let mut dirs = vec![0.0; radiating.len()];
        for (k, &e) in radiating.iter().enumerate() {
            if foreign.contains(&e) {
                continue;
            }
            let el = ctx.mesh.element(e);
            if el.nodes.len() != 2 {
                return Err(unsupported("diffuse gray radiation on higher order elements"));
            }
            let bc = ctx.bc_params(e).expect("radiating elements have a condition");
            let (a, b) = (vf_coord(el.nodes[0])?, vf_coord(el.nodes[1])?);
            let left = parent_info(ctx, el.parents[0], el.nodes);
            let right = parent_info(ctx, el.parents[1], el.nodes);
            let target = bc.integer("Radiation Target Body").or_else(|| bc.integer("Normal Target Body"));
            let dir = line_direction(
                a,
                b,
                left.map(|l| l.0),
                right.map(|r| r.0),
                target,
                (left.is_some_and(|l| l.1), right.is_some_and(|r| r.1)),
            )
            .map_err(|err| SolverError::fatal("ViewFactors", err.to_string()))?;
            dirs[k] = if dir > 0.0 { 1.0 } else { -1.0 };
        }
        if parallel {
            ctx.comm.sum_f64(&mut dirs);
        }
        let dir_of: HashMap<u32, f64> = radiating.iter().zip(&dirs).map(|(e, &d)| (e.0, d)).collect();

        let mut factors: HashMap<u32, Vec<(ElemId, f64)>> = HashMap::new();
        for rb in 1..=max_body {
            let surfaces: Vec<ElemId> = radiating
                .iter()
                .copied()
                .filter(|&e| body_of(e) == rb)
                .collect();
            let n = surfaces.len();
            if n == 0 {
                continue;
            }
            let mut surf = Vec::with_capacity(2 * n);
            // `ViewFactors`' areas (rewritten coordinates) normalise the
            // factors; RadiationFactors' own (`ElementArea` on the mesh)
            // scale them.
            let mut areas = Vec::with_capacity(n);
            let mut vf_areas = Vec::with_capacity(n);
            let mut open = false;
            for &e in &surfaces {
                let el = ctx.mesh.element(e);
                if el.nodes.len() != 2 {
                    return Err(unsupported(
                        "diffuse gray radiation on higher order elements",
                    ));
                }
                let bc = ctx
                    .bc_params(e)
                    .expect("radiating elements have a condition");
                open |= bc.logical("Radiation Boundary Open") == Some(true);
                let (n1, n2) = (el.nodes[0], el.nodes[1]);
                let (a, b) = (vf_coord(n1)?, vf_coord(n2)?);
                let (ea, eb) = (ctx.mesh.coord(n1), ctx.mesh.coord(n2));
                if dir_of[&e.0] > 0.0 {
                    surf.extend([vf_index(n2), vf_index(n1)]);
                } else {
                    surf.extend([vf_index(n1), vf_index(n2)]);
                }
                vf_areas.push(element_area_axis(a, b)?);
                areas.push(element_area_axis(ea, eb)?);
            }

            let mut vf = view_factors_axis(&surf, &coord, divide, combine);
            normalize_factors_with(&mut vf, &vf_areas, open, &mut hut_itersolv);
            let rows = saved_rows(&vf, n, n, min_vf, single);

            let max_area = areas.iter().copied().fold(f64::MIN, f64::max);
            let rel_areas: Vec<f64> = areas.iter().map(|a| a / max_area).collect();
            let rows: Vec<Vec<(usize, f64)>> = rows
                .into_iter()
                .enumerate()
                .map(|(i, r)| r.into_iter().map(|(j, v)| (j, rel_areas[i] * v)).collect())
                .collect();
            let (mut emis, mut refl) = (Vec::with_capacity(n), Vec::with_capacity(n));
            for &e in &surfaces {
                let (em, re) = reflectivity_for_factors(ctx, e)?;
                emis.push(em);
                refl.push(re);
            }
            // InitRadiationSolver: the `radiation:` keywords, with
            // `Linear System Direct Method = Umfpack` unless given.
            let mut rp = elmer_core::ParamMap::default();
            rp.insert(
                "Linear System Direct Method",
                elmer_core::Value::Str(
                    ctx.params().string("radiation: Linear System Direct Method").unwrap_or("umfpack").to_string(),
                ),
            );
            let mut direct = elmer_solver::linsys::direct_solver(&rp)?;
            let mut factored = false;
            let mut solve = |g: &GebhartMatrix, b: &mut [f64]| -> Result<(), String> {
                if !factored {
                    let n = g.row_ptr.len() - 1;
                    let a = CsrMatrix {
                        nrows: n,
                        ncols: n,
                        row_ptr: g.row_ptr.clone(),
                        col_idx: g.cols.clone(),
                        values: g.values.clone(),
                        diag: None,
                    };
                    direct.factor(&a).map_err(|e| e.to_string())?;
                    factored = true;
                }
                direct.solve_in_place(b).map_err(|e| e.to_string())
            };
            let g = gebhart_factors_with(
                &GebhartInput {
                    rows: &rows,
                    rel_areas: &rel_areas,
                    emissivity: &emis,
                    reflectivity: &refl,
                    min_factor: min_gebhart,
                },
                &mut solve,
            )
            .map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
            ctx.info(
                CALLER,
                &format!(
                    "Minimum Gebhart factors sum {:15.6E}, maximum {:15.6E}",
                    g.min_row_sum, g.max_row_sum
                ),
                5,
            );
            for (t, row) in g.rows.into_iter().enumerate() {
                if foreign.contains(&surfaces[t]) {
                    continue;
                }
                factors.insert(
                    surfaces[t].0,
                    row.into_iter().map(|(i, v)| (surfaces[i], v)).collect(),
                );
            }
        }
        radiating.sort_unstable_by_key(|e| e.0);
        Ok(Some(Self {
            factors,
            radiating,
            topology_done: false,
        }))
    }

    /// The list-matrix rows `MakeListMatrix` adds for the factors
    /// (ElementUtils.F90:744-830): each node of a radiating element coupled
    /// to every node of the elements it sees.
    ///
    /// Every radiating element is listed, in element order, an imported one
    /// (without factors on this partition) with no couplings: its nodes are
    /// what InitialPermutation's radiation branch numbers.
    pub(crate) fn couplings(&self, mesh: &dyn elmer_core::MeshView) -> Vec<(u32, Vec<u32>)> {
        let mut out = Vec::new();
        for &e in &self.radiating {
            let cols: Vec<u32> = self
                .factors
                .get(&e.0)
                .into_iter()
                .flatten()
                .flat_map(|&(j, _)| mesh.element(j).nodes.iter().copied())
                .collect();
            for &n in mesh.element(e).nodes {
                out.push((n, cols.clone()));
            }
        }
        out
    }

    /// `UpdateMatrixTopologyWithFactors`: the couplings between the nodes of
    /// every radiating element and those of the elements it sees.
    pub(crate) fn extend_topology(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if self.topology_done {
            return Ok(());
        }
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        let mut a = Vec::new();
        let mut b = Vec::new();
        for (&e, row) in &self.factors {
            radiation_slots(ctx, ElemId(e), &mut a)?;
            for &(j, _) in row {
                radiation_slots(ctx, j, &mut b)?;
                for &s in a.iter().filter(|&&s| s != NO_DOF) {
                    for &t in b.iter().filter(|&&t| t != NO_DOF) {
                        pairs.push((s as usize, t as usize));
                    }
                }
            }
        }
        let sys = ctx
            .system_mut()
            .ok_or_else(|| SolverError::fatal(CALLER, "no system to extend".to_string()))?;
        extend_pattern(&mut sys.matrix, pairs);
        ctx.system_pattern_changed();
        self.topology_done = true;
        Ok(())
    }

    /// `LocalMatrixDiffuseGray` for boundary element `e`, added straight to
    /// the system. `temps4` holds the element means of `T^4`, by element id.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn assemble(
        &self,
        ctx: &mut SolverCtx<'_>,
        basis: &dyn ElementBasis,
        quad: &dyn Quadrature,
        e: ElemId,
        temps4: &HashMap<u32, f64>,
        newton: bool,
        stefan: f64,
    ) -> Result<(), SolverError> {
        let fact = &self.factors[&e.0];
        let el = ctx.mesh.element(e);
        let n = el.nodes.len();
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&nd| ctx.mesh.coord(nd)).collect();
        let (emis1, _, _) = optical(ctx, e)?;
        let bc = ctx
            .bc_params(e)
            .expect("radiating elements have a condition");
        let open = bc.logical("Radiation Boundary Open") == Some(true);
        let angle_fraction = if open {
            fact.iter().map(|f| f.1).sum::<f64>() / emis1
        } else {
            1.0
        };
        let open_text: Option<Vec<f64>> = if open {
            match ctx.real_at_nodes(bc, "Radiation External Temperature", e)? {
                Some(t) => Some(t),
                None => Some(
                    ctx.real_at_nodes(bc, "External Temperature", e)?
                        .ok_or_else(|| {
                            SolverError::fatal(
                                CALLER,
                                "External Temperature is not defined".to_string(),
                            )
                        })?,
                ),
            }
        } else {
            None
        };

        let ips = quad.gauss_points(el.etype, None);
        let mut slots = Vec::new();
        radiation_slots(ctx, e, &mut slots)?;
        let temps: Vec<f64> = {
            let f = ctx.primary().expect("heat solver has a variable");
            slots.iter().map(|&s| f.values[s as usize]).collect()
        };
        let axisym = {
            let cs = ctx.simulation_params().string("coordinate system");
            cs.is_some_and(|c| {
                let c = elmer_core::normalize_key(c);
                c == "axi symmetric" || c == "cylindric symmetric"
            })
        };

        let mut stiff = vec![0.0; n * n];
        let mut force = vec![0.0; n];
        let mut base = vec![0.0; n];
        let mut be = BasisEval::default();
        let mut weights = Vec::with_capacity(ips.len());
        for qp in &ips {
            basis
                .eval(el.etype, &xyz, qp.uvw, &mut be)
                .map_err(|err| SolverError::fatal(CALLER, format!("{err:?}")))?;
            let mut s = be.det_j * qp.weight;
            if axisym {
                s *= be.basis[..n]
                    .iter()
                    .zip(&xyz)
                    .map(|(b, x)| b * x[0])
                    .sum::<f64>();
            }
            weights.push((s, be.basis[..n].to_vec()));
        }
        let mean_t = temps.iter().sum::<f64>() / n as f64;

        let mut text = 0.0;
        let mut couplings: Vec<(u32, u32, f64)> = Vec::new();
        let mut rhs_add: Vec<(u32, f64)> = Vec::new();

        if newton {
            for (s, phi) in &weights {
                let coeff = 4.0 * emis1 * mean_t.powi(3) * stefan;
                let load = 3.0 * emis1 * mean_t.powi(4) * stefan;
                for p in 0..n {
                    for q in 0..n {
                        stiff[p * n + q] += s * phi[p] * phi[q] * coeff;
                    }
                    force[p] += s * phi[p] * load;
                    base[p] += s * phi[p];
                }
            }
            let mut other = Vec::new();
            for &(j, fj) in fact {
                let t4 = temps4[&j.0];
                let tj = t4.powf(0.25);
                radiation_slots(ctx, j, &mut other)?;
                let k = other.len() as f64;
                let coeff = -4.0 * fj * tj.powi(3) * stefan;
                let load = -3.0 * fj * tj.powi(4) * stefan;
                for p in 0..n {
                    for &q in &other {
                        couplings.push((slots[p], q, coeff * base[p] / k));
                    }
                    rhs_add.push((slots[p], load * base[p]));
                }
            }
        } else {
            for &(j, fj) in fact {
                text += fj * temps4[&j.0] / emis1;
            }
        }

        if let Some(at) = &open_text {
            if angle_fraction < 1.0 {
                let topen = at.iter().map(|t| t.powi(4)).sum::<f64>().powf(0.25);
                if newton {
                    let load = (1.0 - angle_fraction) * emis1 * topen.powi(4) * stefan;
                    for p in 0..n {
                        force[p] += base[p] * load;
                    }
                } else {
                    text += (1.0 - angle_fraction) * topen.powi(4);
                }
            }
        }

        if !newton {
            base.iter_mut().for_each(|b| *b = 0.0);
            let text = text.powf(0.25);
            for (s, phi) in &weights {
                let t0: f64 = phi.iter().zip(&temps).map(|(b, t)| b * t).sum();
                let coeff = emis1
                    * stefan
                    // T0**3 + T0**2*Text + T0*Text**2 + Text**3.
                    * (t0.powi(3) + t0 * t0 * text + t0 * (text * text) + text.powi(3));
                for p in 0..n {
                    for q in 0..n {
                        stiff[p * n + q] += s * phi[p] * phi[q] * coeff;
                    }
                    force[p] += s * phi[p] * coeff * text;
                    base[p] += s * phi[p];
                }
            }
        }

        let sys = ctx
            .system_mut()
            .ok_or_else(|| SolverError::fatal(CALLER, "no assembled system".to_string()))?;
        for (r, c, v) in couplings {
            add_entry(&mut sys.matrix, r as usize, c as usize, v)?;
        }
        for (r, v) in rhs_add {
            sys.rhs[r as usize] += v;
        }
        for p in 0..n {
            for q in 0..n {
                add_entry(
                    &mut sys.matrix,
                    slots[p] as usize,
                    slots[q] as usize,
                    stiff[p * n + q],
                )?;
            }
            sys.rhs[slots[p] as usize] += force[p];
        }
        Ok(())
    }
}

/// The system slots of a radiating boundary element: its nodes', or for a
/// DG solver `DgRadiationIndexes(Element, n, ElemInds, .TRUE.)`, the DG
/// indexes of its nodes in the first parent whose material has an
/// `Emissivity`.
pub(crate) fn radiation_slots(ctx: &SolverCtx<'_>, e: ElemId, out: &mut Vec<u32>) -> Result<(), SolverError> {
    match ctx.solver_dofs().and_then(|d| d.dg_layout()) {
        Some(ix) => {
            *out = crate::solver::dg_diffuse_gray_slots(ctx, ix, e)?;
            Ok(())
        }
        None => {
            ctx.element_slots(e, out);
            Ok(())
        }
    }
}

/// `TabulateBoundaryAverages`' `Temps4`, by boundary element id.
pub(crate) fn temps4(ctx: &SolverCtx<'_>, elems: &[ElemId]) -> Result<HashMap<u32, f64>, SolverError> {
    let mut slots = Vec::new();
    let mut out = HashMap::new();
    for &e in elems.iter().filter(|&&e| is_diffuse_gray(ctx, e)) {
        radiation_slots(ctx, e, &mut slots)?;
        let f = ctx.primary().expect("heat solver has a variable");
        // SUM(NodalTemp(1:n)**4) / n.
        let sum: f64 = slots.iter().map(|&s| f.values[s as usize].powi(4)).sum();
        out.insert(e.0, sum / slots.len() as f64);
    }
    Ok(out)
}

fn add_entry(a: &mut CsrMatrix<f64>, row: usize, col: usize, v: f64) -> Result<(), SolverError> {
    let (lo, hi) = (a.row_ptr[row], a.row_ptr[row + 1]);
    match a.col_idx[lo..hi].binary_search(&(col as u32)) {
        Ok(p) => {
            a.values[lo + p] += v;
            Ok(())
        }
        Err(_) => Err(SolverError::fatal(
            CALLER,
            format!("no matrix entry ({row}, {col}) for a radiation coupling"),
        )),
    }
}

/// Adds the `pairs` to the pattern of `a`, keeping the columns of a row
/// sorted and the values and the diagonal index consistent.
pub(crate) fn extend_pattern(a: &mut CsrMatrix<f64>, mut pairs: Vec<(usize, usize)>) {
    pairs.sort_unstable();
    pairs.dedup();
    let mut row_ptr = Vec::with_capacity(a.nrows + 1);
    let mut col_idx = Vec::with_capacity(a.col_idx.len() + pairs.len());
    let mut values = Vec::with_capacity(a.values.len() + pairs.len());
    let mut at = 0;
    row_ptr.push(0);
    for r in 0..a.nrows {
        let mut row: Vec<(u32, f64)> = (a.row_ptr[r]..a.row_ptr[r + 1])
            .map(|k| (a.col_idx[k], a.values[k]))
            .collect();
        while at < pairs.len() && pairs[at].0 == r {
            let c = pairs[at].1 as u32;
            if !row.iter().any(|&(x, _)| x == c) {
                row.push((c, 0.0));
            }
            at += 1;
        }
        row.sort_by_key(|&(c, _)| c);
        for (c, v) in row {
            col_idx.push(c);
            values.push(v);
        }
        row_ptr.push(col_idx.len());
    }
    let diag = (0..a.nrows)
        .map(|r| {
            col_idx[row_ptr[r]..row_ptr[r + 1]]
                .binary_search(&(r as u32))
                .map_or(usize::MAX, |p| row_ptr[r] + p)
        })
        .collect();
    a.row_ptr = row_ptr;
    a.col_idx = col_idx;
    a.values = values;
    a.diag = Some(diag);
}
