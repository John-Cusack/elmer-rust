//! Element integrals of SaveScalars: `BulkIntegrals` (SaveScalars.F90:2153),
//! `BoundaryIntegrals` (2438) and `BoundaryStatistics` (2862), plus the
//! nodal statistics `VectorStatistics` (1846) and `VectorMeanDeviation`
//! (2043).
//!
//! Integrals need an element basis and a quadrature rule, which come from
//! `elmer-fe` through [`crate::install_fe`]. Everything is serial; there are
//! no MPI reductions.

use std::f64::consts::PI;

use elmer_core::{normalize_key, ElemId, ElementBasis, BasisEval, ElementType, Params, Quadrature};
use elmer_solver::{SolverCtx, SolverError};

use crate::plan::MaskKind;
use crate::view::VarView;

/// Basis and quadrature of the front end's element library.
pub(crate) struct Fe {
    pub basis: Box<dyn ElementBasis + Send>,
    pub quad: Box<dyn Quadrature + Send>,
}

/// One integration point of an element: `s = DetJ * weight` (with the
/// axisymmetric `2 pi r` already in), the basis values and their physical
/// gradients.
pub(crate) struct Ip {
    pub s: f64,
    pub basis: Vec<f64>,
    pub grad: Vec<[f64; 3]>,
}

/// The `Coordinate System` as far as the integrals care.
pub(crate) fn axisymmetric(ctx: &SolverCtx<'_>) -> Result<bool, SolverError> {
    match ctx.simulation_params().string("coordinate system") {
        None => Ok(false),
        Some(s) => {
            let s = normalize_key(s);
            if s.starts_with("cartesian") {
                Ok(false)
            } else if s == "axi symmetric" || s == "cylindric symmetric" {
                Ok(true)
            } else {
                Err(SolverError::Unsupported(format!("SaveScalars integrals with Coordinate System = {s}")))
            }
        }
    }
}

/// Basis, `DetJ` and quadrature weight at every integration point of `e`
/// (`GaussPoints` + `ElementInfo`). `passive` is Elmer's `Passive
/// Coordinate`, 1-based: that coordinate of every node is taken as zero.
pub(crate) fn element_points(
    ctx: &SolverCtx<'_>,
    fe: &Fe,
    e: ElemId,
    passive: Option<usize>,
    axisym: bool,
) -> Result<Vec<Ip>, SolverError> {
    let el = ctx.mesh.element(e);
    let mut xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
    if let Some(p) = passive.filter(|p| (1..=3).contains(p)) {
        xyz.iter_mut().for_each(|x| x[p - 1] = 0.0);
    }
    let mut out = BasisEval::default();
    let mut pts = Vec::new();
    for q in fe.quad.gauss_points(el.etype, None) {
        fe.basis
            .eval(el.etype, &xyz, q.uvw, &mut out)
            .map_err(|m| SolverError::fatal("ElementInfo", format!("element {}: {m}", e.0 + 1)))?;
        let mut s = out.det_j * q.weight;
        if axisym {
            let x: f64 = xyz.iter().zip(&out.basis).map(|(p, b)| p[0] * b).sum();
            s *= 2.0 * PI * x;
        }
        pts.push(Ip {
            s,
            basis: out.basis.clone(),
            grad: out.dbasis_dx.clone(),
        });
    }
    Ok(pts)
}

fn is_point(t: ElementType) -> bool {
    t.0 == 101
}

/// The section a masked operator reads its logical from.
fn mask_holds(ctx: &SolverCtx<'_>, e: ElemId, kind: MaskKind, name: &str) -> bool {
    let Some(bp) = ctx.body_params(e) else { return false };
    let list = match kind {
        MaskKind::Body => Some(bp.body),
        MaskKind::BodyForce => bp.body_force,
        MaskKind::Material => bp.material,
        MaskKind::None => return true,
    };
    list.is_some_and(|l| l.logical(name) == Some(true))
}

/// `CreateNodeMask`: nodes of the bulk elements whose body, body force or
/// material has `name` true.
pub(crate) fn node_mask(ctx: &SolverCtx<'_>, kind: MaskKind, name: &str, len: usize) -> Vec<bool> {
    let mut mask = vec![false; len];
    for t in 0..ctx.mesh.n_bulk() as u32 {
        let e = ElemId(t);
        if mask_holds(ctx, e, kind, name) {
            for &n in ctx.mesh.element(e).nodes {
                if let Some(m) = mask.get_mut(n as usize) {
                    *m = true;
                }
            }
        }
    }
    mask
}

/// Nodes per `Nodal Variable` flag: the first `n_nodes` permutation entries
/// are nodal, the rest belong to edge and p-element dofs.
fn index_range(view: &VarView<'_>, nodal: Option<bool>, n_nodes: usize) -> std::ops::Range<usize> {
    let n = view.n_perm();
    match nodal {
        None => 0..n,
        Some(true) => 0..n_nodes.min(n),
        Some(false) => n_nodes.min(n)..n,
    }
}

/// `VectorStatistics`: `sum`, `max`, `mean`, `variance`, ... over the nodes
/// of the field, the Euclidean norm per node for a vector field.
pub(crate) fn vector_statistics(
    view: &VarView<'_>,
    oper: &str,
    nodal: Option<bool>,
    mask: Option<&[bool]>,
    n_nodes: usize,
) -> f64 {
    let (mut count, mut sumx, mut sumxx, mut sumabs) = (0usize, 0.0, 0.0, 0.0);
    let (mut max, mut min) = (f64::MAX.copysign(-1.0), f64::MAX);
    let (mut abs_max, mut abs_min) = (0.0f64, f64::MAX);
    let mut range = index_range(view, nodal, n_nodes);
    if let Some(m) = mask {
        range.end = range.end.min(m.len());
    }
    for i in range {
        if mask.is_some_and(|m| !m[i]) {
            continue;
        }
        let Some(slot) = view.slot(i) else { continue };
        let x = view.magnitude(slot);
        count += 1;
        sumx += x;
        sumxx += x * x;
        sumabs += x.abs();
        max = max.max(x);
        min = min.min(x);
        if x.abs() > abs_max.abs() {
            abs_max = x;
        }
        if x.abs() < abs_min.abs() {
            abs_min = x;
        }
    }
    if count == 0 {
        return 0.0;
    }
    let c = count as f64;
    match oper {
        "sum" => sumx,
        "sum square" => sumxx,
        "sum abs" => sumabs,
        "max" => max,
        "min" => min,
        "range" => max - min,
        "max abs" => abs_max,
        "min abs" => abs_min,
        "mean" => sumx / c,
        "mean square" => sumxx / c,
        "rms" => (sumxx / c).sqrt(),
        "mean abs" => sumabs / c,
        "variance" => {
            // Despite the name this is the standard deviation (F90:2025).
            let mean = sumx / c;
            let v = sumxx / c - mean * mean;
            if v > 0.0 {
                v.sqrt()
            } else {
                0.0
            }
        }
        _ => unreachable!("statistic {oper} is checked when the plan is built"),
    }
}

/// `VectorMeanDeviation`: the mean absolute deviation from the mean.
pub(crate) fn vector_mean_deviation(view: &VarView<'_>, mask: Option<&[bool]>) -> f64 {
    let mut n = view.n_perm();
    if let Some(m) = mask {
        n = n.min(m.len());
    }
    let values: Vec<f64> = (0..n)
        .filter(|&i| mask.is_none_or(|m| m[i]))
        .filter_map(|i| view.slot(i))
        .map(|s| view.magnitude(s))
        .collect();
    if values.is_empty() {
        return 0.0;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    values.iter().map(|x| (x - mean).abs()).sum::<f64>() / values.len() as f64
}

/// What `BulkIntegrals` is asked for.
pub(crate) struct BulkRequest<'a> {
    pub oper: &'a str,
    pub mask: MaskKind,
    pub mask_name: &'a str,
    pub passive: Option<usize>,
}

/// `BulkIntegrals`. Without a material coefficient, so the `coeff` of the
/// Fortran is 1 everywhere.
pub(crate) fn bulk_integral(
    ctx: &SolverCtx<'_>,
    fe: &Fe,
    view: &VarView<'_>,
    req: &BulkRequest<'_>,
) -> Result<f64, SolverError> {
    let axisym = axisymmetric(ctx)?;
    let dim = ctx.mesh.dim();
    if req.oper == "diffusive energy" && view.ncomp != 1 {
        return Err(SolverError::Unsupported(
            "SaveScalars: diffusive energy of a vector field".into(),
        ));
    }
    let (mut integral1, mut integral2, mut vol) = (0.0, 0.0, 0.0);
    let mut hits = 0usize;
    let nb = ctx.mesh.n_bulk() as u32;
    let total = nb + ctx.mesh.n_boundary() as u32;
    for t in 0..total {
        // Boundary elements are visited only when no bulk element took part.
        if t == nb && hits > 0 {
            break;
        }
        let e = ElemId(t);
        let el = ctx.mesh.element(e);
        if is_point(el.etype) {
            continue;
        }
        if view.perm.is_some() && el.nodes.iter().any(|&n| view.slot(n as usize).is_none()) {
            continue;
        }
        hits += 1;
        if req.mask != MaskKind::None && !mask_holds(ctx, e, req.mask, req.mask_name) {
            continue;
        }
        for ip in element_points(ctx, fe, e, req.passive, axisym)? {
            vol += ip.s;
            let n = el.nodes.len();
            let comp = |k: usize| -> f64 {
                (0..n)
                    .map(|i| ip.basis[i] * view.slot(el.nodes[i] as usize).map_or(0.0, |s| view.at(s, k)))
                    .sum()
            };
            let func = if view.ncomp > 1 {
                (0..view.ncomp).map(|k| comp(k).powi(2)).sum::<f64>().sqrt()
            } else {
                comp(0)
            };
            match req.oper {
                "volume" => integral1 += ip.s,
                "int" | "int mean" | "potential energy" => integral1 += ip.s * func,
                "int square" | "int square mean" | "int rms" => integral1 += ip.s * func * func,
                "int abs" | "int abs mean" => integral1 += ip.s * func.abs(),
                "int variance" => {
                    integral1 += ip.s * func;
                    integral2 += ip.s * func * func;
                }
                "diffusive energy" => {
                    let mut grad = [0.0; 3];
                    for i in 0..n {
                        let v = view.slot(el.nodes[i] as usize).map_or(0.0, |s| view.at(s, 0));
                        for (g, d) in grad.iter_mut().zip(&ip.grad[i]) {
                            *g += v * d;
                        }
                    }
                    integral1 += ip.s * grad[..dim].iter().map(|g| g * g).sum::<f64>();
                }
                _ => unreachable!("bulk operator {} is checked when the plan is built", req.oper),
            }
        }
    }
    if hits == 0 {
        return Ok(0.0);
    }
    Ok(match req.oper {
        "volume" | "int" | "int abs" | "int square" | "potential energy" => integral1,
        "int mean" | "int square mean" | "int abs mean" => integral1 / vol,
        "int rms" => (integral1 / vol).sqrt(),
        "int variance" => (integral2 / vol - (integral1 / vol).powi(2)).sqrt(),
        "diffusive energy" => 0.5 * integral1,
        _ => unreachable!(),
    })
}

/// Result of `BoundaryIntegrals` for the boundary conditions that took part.
pub(crate) struct BoundaryTotals {
    pub fluxes: Vec<f64>,
    pub areas: Vec<f64>,
    /// `Save Flux Range`: the smallest and largest point flux of a flux
    /// operator.
    pub range: (f64, f64),
}

/// Whether BC `bc` (a section) is one of the `ActiveBC` of SaveScalars.F90:641.
pub(crate) fn bc_is_active(bc: &elmer_core::ParamMap, mask_name: &str) -> bool {
    bc.logical("flux integrate") == Some(true) || bc.logical(mask_name) == Some(true)
}

/// The parent `BoundaryIntegrals` takes the gradient from for `diffusive
/// flux`: the one of `Flux Integrate Body`, else the left parent unless the
/// field misses one of its nodes, then the right one.
fn flux_parent(ctx: &SolverCtx<'_>, e: ElemId, bc: &elmer_core::ParamMap, v: &VarView<'_>) -> Result<ElemId, SolverError> {
    const CALLER: &str = "SaveScalars";
    let [left, right] = ctx.mesh.element(e).parents;
    if let Some(body) = bc.integer("flux integrate body") {
        let of = |p: Option<ElemId>| p.filter(|&p| ctx.body_params(p).is_some_and(|b| b.body_id as i32 == body));
        return of(left).or_else(|| of(right)).ok_or_else(|| {
            SolverError::fatal(CALLER, format!("No such flux integrate body on bc {}", ctx.bc_id(e).unwrap_or(0)))
        });
    }
    let has_field = |p: ElemId| v.perm.is_none() || ctx.mesh.element(p).nodes.iter().all(|&n| v.slot(n as usize).is_some());
    left.filter(|&p| has_field(p))
        .or_else(|| right.filter(|&p| has_field(p)))
        .ok_or_else(|| SolverError::fatal(CALLER, "No solution available for specified boundary"))
}

/// `NormalVector(Element, Nodes, u, v, Check=.TRUE.)` of a line (2D) or a
/// triangle/quadrilateral (3D) boundary element, pointing away from
/// `parent`.
fn outward_normal(xyz: &[[f64; 3]], etype: ElementType, parent_xyz: &[[f64; 3]], point: [f64; 3]) -> Result<[f64; 3], SolverError> {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let mut n = match etype.0 / 100 {
        2 => {
            let t = sub(xyz[1], xyz[0]);
            [t[1], -t[0], 0.0]
        }
        3 | 4 => {
            let (a, b) = if etype.0 / 100 == 3 { (sub(xyz[1], xyz[0]), sub(xyz[2], xyz[0])) } else { (sub(xyz[2], xyz[0]), sub(xyz[3], xyz[1])) };
            [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
        }
        _ => return Err(SolverError::Unsupported(format!("SaveScalars diffusive flux on boundary element {}", etype.0))),
    };
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    n = [n[0] / len, n[1] / len, n[2] / len];
    let k = parent_xyz.len() as f64;
    let c = parent_xyz.iter().fold([0.0; 3], |a, x| [a[0] + x[0] / k, a[1] + x[1] / k, a[2] + x[2] / k]);
    let out = sub(point, c);
    if n[0] * out[0] + n[1] * out[1] + n[2] * out[2] < 0.0 {
        n = [-n[0], -n[1], -n[2]];
    }
    Ok(n)
}

/// `diffusive flux` at one integration point of boundary element `e`
/// (`EvaluateVariableAtGivenPoint(..., DoGrad, Parent)`, SaveScalars.F90
/// :2754): the normal component of the field's gradient taken in `parent`
/// at the point. Without `Coefficient` the energy tensor is the identity.
fn diffusive_flux(ctx: &SolverCtx<'_>, fe: &Fe, e: ElemId, parent: ElemId, v: &VarView<'_>, basis: &[f64]) -> Result<f64, SolverError> {
    let el = ctx.mesh.element(e);
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
    let mut x = [0.0; 3];
    for (b, p) in basis.iter().zip(&xyz) {
        for k in 0..3 {
            x[k] += b * p[k];
        }
    }
    let pel = ctx.mesh.element(parent);
    let pxyz: Vec<[f64; 3]> = pel.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
    let shape = elmer_interp::shape::Shape::get(pel.etype.0)
        .ok_or_else(|| SolverError::Unsupported(format!("SaveScalars diffusive flux on element {}", pel.etype.0)))?;
    let dim = ctx.mesh.dim();
    let uvw = elmer_interp::geom::global_to_local(shape, &pxyz, x, dim);
    let mut out = BasisEval::default();
    fe.basis
        .eval(pel.etype, &pxyz, uvw, &mut out)
        .map_err(|m| SolverError::fatal("ElementInfo", format!("element {}: {m}", parent.0 + 1)))?;
    let mut grad = [0.0; 3];
    for (k, &n) in pel.nodes.iter().enumerate() {
        let u = v.slot(n as usize).map_or(0.0, |s| v.at(s, 0));
        for d in 0..3 {
            grad[d] += out.dbasis_dx[k][d] * u;
        }
    }
    let nrm = outward_normal(&xyz, el.etype, &pxyz, x)?;
    Ok((0..dim).map(|d| nrm[d] * grad[d]).sum())
}

/// `BoundaryIntegrals` for `boundary int`, `boundary int mean`, `area` and
/// `diffusive flux`; per boundary condition, 0-based. `convective flux` is
/// refused earlier.
pub(crate) fn boundary_integrals(
    ctx: &SolverCtx<'_>,
    fe: &Fe,
    view: Option<&VarView<'_>>,
    oper: &str,
    mask_name: &str,
    passive: Option<usize>,
) -> Result<BoundaryTotals, SolverError> {
    let axisym = axisymmetric(ctx)?;
    let nbc = ctx.sim.model.boundary_conditions.len();
    let mut totals = BoundaryTotals {
        fluxes: vec![0.0; nbc],
        areas: vec![0.0; nbc],
        range: (f64::MAX, -f64::MAX),
    };
    let area_only = oper == "area";
    let nb = ctx.mesh.n_bulk() as u32;
    for t in nb..nb + ctx.mesh.n_boundary() as u32 {
        let e = ElemId(t);
        let Some(bc_ix) = ctx.bc_id(e) else { continue };
        let bc = &ctx.sim.model.boundary_conditions[bc_ix - 1];
        if !bc_is_active(bc, mask_name) {
            continue;
        }
        let el = ctx.mesh.element(e);
        if is_point(el.etype) {
            continue;
        }
        let field = if area_only { None } else { view };
        if let Some(v) = field {
            if v.perm.is_some() && el.nodes.iter().any(|&n| v.slot(n as usize).is_none()) {
                continue;
            }
        }
        let parent = match (oper, field) {
            ("diffusive flux", Some(v)) => {
                if v.ncomp != 1 {
                    return Err(SolverError::fatal("SaveScalars", "diffusive flux & NoDofs /= 1?"));
                }
                Some(flux_parent(ctx, e, bc, v)?)
            }
            _ => None,
        };
        for ip in element_points(ctx, fe, e, passive, axisym)? {
            let flux = match field {
                None => ip.s,
                Some(v) if parent.is_some() => {
                    let flux = diffusive_flux(ctx, fe, e, parent.unwrap(), v, &ip.basis)?;
                    totals.range = (totals.range.0.min(flux), totals.range.1.max(flux));
                    ip.s * flux
                }
                Some(v) => {
                    let comp = |k: usize| -> f64 {
                        el.nodes
                            .iter()
                            .zip(&ip.basis)
                            .map(|(&n, b)| b * v.slot(n as usize).map_or(0.0, |s| v.at(s, k)))
                            .sum()
                    };
                    let func = if v.ncomp > 1 {
                        (0..v.ncomp).map(|k| comp(k).powi(2)).sum::<f64>().sqrt()
                    } else {
                        comp(0)
                    };
                    ip.s * func
                }
            };
            totals.fluxes[bc_ix - 1] += flux;
            totals.areas[bc_ix - 1] += ip.s;
        }
    }
    Ok(totals)
}

/// `BoundaryStatistics`: one value over the distinct nodes of the boundary
/// elements whose BC has `mask_name` true. Returns the value and the number
/// of nodes visited.
pub(crate) fn boundary_statistics(
    ctx: &SolverCtx<'_>,
    view: &VarView<'_>,
    oper: &str,
    mask_name: &str,
) -> (f64, usize) {
    let (find_min_max, mut value) = match oper {
        "boundary min" | "boundary min abs" => (true, f64::MAX),
        "boundary max" | "boundary max abs" => (true, -f64::MAX),
        _ => (false, 0.0),
    };
    let mut seen = vec![false; view.values.len() / view.stride.max(1)];
    let mut count = 0usize;
    let nb = ctx.mesh.n_bulk() as u32;
    for t in nb..nb + ctx.mesh.n_boundary() as u32 {
        let e = ElemId(t);
        let Some(bc) = ctx.bc_params(e) else { continue };
        if bc.logical(mask_name) != Some(true) {
            continue;
        }
        for &n in ctx.mesh.element(e).nodes {
            let Some(slot) = view.slot(n as usize) else { continue };
            if std::mem::replace(&mut seen[slot], true) {
                continue;
            }
            let val = if view.ncomp <= 1 {
                view.at(slot, 0)
            } else {
                // BoundaryStatistics sums the components and takes the root
                // without squaring them first (SaveScalars.F90:2965-2968).
                (0..view.ncomp).map(|k| view.at(slot, k)).sum::<f64>().sqrt()
            };
            if find_min_max {
                match oper {
                    "boundary min" => value = value.min(val),
                    "boundary max" => value = value.max(val),
                    "boundary min abs" => {
                        if value.abs() < val.abs() {
                            value = val;
                        }
                    }
                    _ => {
                        if value.abs() > val.abs() {
                            value = val;
                        }
                    }
                }
            } else {
                value += val;
            }
            count += 1;
        }
    }
    if oper == "boundary dofs" {
        value = count as f64;
    }
    (value, count)
}
