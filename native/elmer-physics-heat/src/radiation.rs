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
    gebhart_factors, line_direction, normalize_factors, saved_rows, view_factors_axis,
    GebhartInput, ParentInfo,
};

const CALLER: &str = "RadiationFactors";

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(what.to_string())
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
        let radiating: Vec<ElemId> = bnd
            .iter()
            .copied()
            .filter(|&e| is_diffuse_gray(ctx, e))
            .collect();
        if radiating.is_empty() {
            return Ok(None);
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

        let n_nodes = ctx.mesh.n_nodes();
        let coord: Vec<f64> = (0..n_nodes as u32)
            .flat_map(|n| {
                let x = ctx.mesh.coord(n);
                [x[0], x[1]]
            })
            .collect();
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
            let mut areas = Vec::with_capacity(n);
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
                let (a, b) = (ctx.mesh.coord(n1), ctx.mesh.coord(n2));
                let left = parent_info(ctx, el.parents[0], el.nodes);
                let right = parent_info(ctx, el.parents[1], el.nodes);
                let target = bc
                    .integer("Radiation Target Body")
                    .or_else(|| bc.integer("Normal Target Body"));
                let dir = line_direction(
                    a,
                    b,
                    left.map(|l| l.0),
                    right.map(|r| r.0),
                    target,
                    (left.is_some_and(|l| l.1), right.is_some_and(|r| r.1)),
                )
                .map_err(|err| SolverError::fatal("ViewFactors", err.to_string()))?;
                if dir > 0.0 {
                    surf.extend([n2 as usize, n1 as usize]);
                } else {
                    surf.extend([n1 as usize, n2 as usize]);
                }
                let len = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
                areas.push(len * (a[0] + b[0]).abs() / 2.0);
            }

            let mut vf = view_factors_axis(&surf, &coord, divide, combine);
            normalize_factors(&mut vf, &areas, open);
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
            let g = gebhart_factors(&GebhartInput {
                rows: &rows,
                rel_areas: &rel_areas,
                emissivity: &emis,
                reflectivity: &refl,
                min_factor: min_gebhart,
            })
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
                factors.insert(
                    surfaces[t].0,
                    row.into_iter().map(|(i, v)| (surfaces[i], v)).collect(),
                );
            }
        }
        Ok(Some(Self {
            factors,
            topology_done: false,
        }))
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
            ctx.element_slots(ElemId(e), &mut a);
            for &(j, _) in row {
                ctx.element_slots(j, &mut b);
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
        ctx.element_slots(e, &mut slots);
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
                ctx.element_slots(j, &mut other);
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
                    * (t0.powi(3) + t0 * t0 * text + t0 * text * text + text.powi(3));
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

/// `TabulateBoundaryAverages`' `Temps4`, by boundary element id.
pub(crate) fn temps4(ctx: &SolverCtx<'_>, elems: &[ElemId]) -> HashMap<u32, f64> {
    let f = ctx.primary().expect("heat solver has a variable");
    elems
        .iter()
        .filter(|&&e| is_diffuse_gray(ctx, e))
        .map(|&e| {
            let nodes = ctx.mesh.element(e).nodes;
            let sum: f64 = nodes
                .iter()
                .map(|&nd| {
                    let s = f.perm.map_or(nd, |p| p[nd as usize]);
                    f.values[s as usize].powi(4)
                })
                .sum();
            (e.0, sum / nodes.len() as f64)
        })
        .collect()
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
fn extend_pattern(a: &mut CsrMatrix<f64>, mut pairs: Vec<(usize, usize)>) {
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
