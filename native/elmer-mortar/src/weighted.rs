//! The Galerkin projector (`WeightedProjector`, InterpolateMeshToMesh.F90)
//! for nodal degrees of freedom.

use std::collections::BTreeMap;

use elmer_core::{BasisEval, ElemId, Logger, MeshView};
use elmer_fe::{element_info, element_type, gauss_points_rel, LagrangeBasis, MetricStatus};
use elmer_interp::{build_projector, Extras, InterpOptions, Locator};

use crate::interface::InterfaceMesh;
use crate::projector::Projector;
use crate::MortarError;

/// `List_AddToMatrixElement`: a new entry takes the value as is.
fn add(row: &mut BTreeMap<u32, f64>, col: u32, v: f64) {
    row.entry(col).and_modify(|x| *x += v).or_insert(v);
}

fn fe(e: elmer_fe::FeError) -> MortarError {
    MortarError::Fatal(format!("WeightedProjector: {e}"))
}

/// `WeightedProjector(BMesh2, BMesh1, ..., PeriodicScale)` without edge
/// DOFs, repeating (rotational) meshes or axisymmetric weighting: row `i`
/// integrates `N_p (u_this - scale * u_target)` over the elements of `b1`
/// with test function `N_p` of the `i`-th node of `b1` met, in element
/// order. Target values come from interpolating `b2` at the Gauss points
/// of `b1` (`MeshProjector`). `rel_order` is `Projector Relative
/// Integration Order`.
pub fn weighted_projector(
    b2: &InterfaceMesh,
    b1: &InterfaceMesh,
    opts: InterpOptions,
    coord_dim: usize,
    rel_order: i32,
    scale: f64,
    log: &dyn Logger,
) -> Result<Projector, MortarError> {
    log.info(
        "WeightedProjector",
        "Creating Galerkin projector between two boundaries",
        7,
    );
    let n_elem = b1.n_bulk();
    let mut status = vec![MetricStatus::default(); n_elem];
    let mut out = BasisEval::default();
    let xyz_of = |e: usize| -> Vec<[f64; 3]> {
        b1.element(ElemId(e as u32))
            .nodes
            .iter()
            .map(|&n| b1.coords[n as usize])
            .collect()
    };

    // The global coordinates of every Gauss point, in element order.
    let mut gauss: Vec<[f64; 3]> = Vec::new();
    for (e, st) in status.iter_mut().enumerate() {
        let el = b1.element(ElemId(e as u32));
        let def = element_type(el.etype.0).map_err(fe)?;
        let rule = gauss_points_rel(el.etype, rel_order).map_err(fe)?;
        let xyz = xyz_of(e);
        for j in 0..rule.s.len() {
            let uvw = [rule.u[j], rule.v[j], rule.w[j]];
            element_info(def, &xyz, uvw, coord_dim, st, &mut out).map_err(fe)?;
            let mut x = [0.0; 3];
            for (d, xd) in x.iter_mut().enumerate() {
                let mut s = 0.0;
                for (b, p) in out.basis.iter().zip(&xyz) {
                    s += b * p[d];
                }
                *xd = s;
            }
            gauss.push(x);
        }
    }
    log.info(
        "WeightedProjector",
        &format!(
            "Number of nodes and gauss points:{} and {}",
            b1.n_nodes(),
            gauss.len()
        ),
        10,
    );

    // MeshProjector(BMesh2, Gauss points).
    let points = InterfaceMesh::points(gauss);
    let loc = Locator::new(b2, opts)?;
    let gp = build_projector(
        &loc,
        &points,
        &LagrangeBasis::new(coord_dim),
        &Extras::default(),
        log,
    )?;

    // The list matrix: rows numbered as the nodes of b1 are met, columns
    // parent nodes.
    let mut eq: Vec<Option<usize>> = vec![None; b1.n_nodes()];
    let mut inv_perm: Vec<u32> = Vec::new();
    let mut list: Vec<BTreeMap<u32, f64>> = Vec::new();
    let mut g = 0usize;
    for (e, st) in status.iter_mut().enumerate() {
        let el = b1.element(ElemId(e as u32));
        let def = element_type(el.etype.0).map_err(fe)?;
        let rule = gauss_points_rel(el.etype, rel_order).map_err(fe)?;
        let xyz = xyz_of(e);
        let nodes = el.nodes;
        for j in 0..rule.s.len() {
            let uvw = [rule.u[j], rule.v[j], rule.w[j]];
            element_info(def, &xyz, uvw, coord_dim, st, &mut out).map_err(fe)?;
            let weight = out.det_j * rule.s[j];
            for &n in nodes {
                if eq[n as usize].is_none() {
                    eq[n as usize] = Some(list.len());
                    list.push(BTreeMap::new());
                    inv_perm.push(b1.inv_perm[n as usize]);
                }
            }
            for (p, &np) in nodes.iter().enumerate() {
                let row = &mut list[eq[np as usize].unwrap()];
                let val = weight * out.basis[p];
                for (q, &nq) in nodes.iter().enumerate() {
                    add(row, b1.inv_perm[nq as usize], out.basis[q] * val);
                }
                for k in gp.rows[g]..gp.rows[g + 1] {
                    let c = b2.inv_perm[gp.cols[k] as usize];
                    add(row, c, -(scale * gp.values[k] * val));
                }
            }
            g += 1;
        }
    }

    let mut rows = vec![0usize];
    let mut cols = Vec::new();
    let mut values = Vec::new();
    for row in list {
        for (c, v) in row {
            cols.push(c);
            values.push(v);
        }
        rows.push(cols.len());
    }
    Ok(Projector {
        rows,
        cols,
        values,
        inv_perm,
        projector_bc: 0,
        galerkin: true,
        slave_scale: 1.0,
        master_scale: 1.0,
        mortar_perm: false,
        rhs: None,
    })
}

/// `WeightedProjectorDiscont(Mesh, bc)` for nodal degrees of freedom of a
/// serial mesh split by `CreateDiscontMesh`: row `j` (the `j`-th copy of
/// `discont_perm`) integrates `N_p (u - scale * u_copy)` over the boundary
/// elements of `this` (`Mortar BC Scaling`, default -1). With `Mortar BC
/// Diag` (or `Use Biorthogonal Basis`, or the solver's `Eliminate Linear
/// Constraints`) the row is lumped onto the node and its copy. Columns
/// ascending (`CRS_SortMatrix`). `None` when nothing was integrated.
#[allow(clippy::too_many_arguments)]
pub fn weighted_projector_discont(
    mesh: &dyn MeshView,
    bc_of: &dyn Fn(i32) -> Option<usize>,
    this: usize,
    discont_perm: &[u32],
    bc: &dyn elmer_core::Params,
    solver: Option<&dyn elmer_core::Params>,
    simulation: &dyn elmer_core::Params,
    coord_dim: usize,
    log: &dyn Logger,
) -> Result<Option<Projector>, MortarError> {
    const CALLER: &str = "WeightedProjectorDiscont";
    log.info(
        CALLER,
        &format!("Creating projector for discontinuous boundary {this}"),
        7,
    );
    let sflag = |k: &str| solver.is_some_and(|s| s.logical(k) == Some(true));
    let skip = |k: &str| sflag(k) || bc.logical(k) == Some(true);
    if sflag("hcurl basis") && !skip("projector skip edges") {
        return Err(MortarError::Unsupported(
            "WeightedProjectorDiscont edge degrees of freedom".into(),
        ));
    }
    if skip("projector skip nodes") {
        return Ok(None);
    }
    let scale = bc.real("mortar bc scaling").unwrap_or(-1.0);
    // Axisymmetric weighting depends on Projector Metrics and the nodal
    // jump keywords (`ListCheckPrefix`).
    let cs = simulation
        .string("coordinate system")
        .unwrap_or("")
        .to_ascii_lowercase();
    if cs.contains("symmetric") {
        return Err(MortarError::Unsupported(
            "WeightedProjectorDiscont in axisymmetric coordinates".into(),
        ));
    }
    let set_diag = bc
        .logical("mortar bc diag")
        .or_else(|| bc.logical("use biorthogonal basis"))
        .unwrap_or_else(|| sflag("eliminate linear constraints"));

    let n_orig = discont_perm.len();
    let m = discont_perm.iter().filter(|&&j| j > 0).count();
    if m == 0 {
        return Ok(None);
    }
    let mut list: Vec<BTreeMap<u32, f64>> = vec![BTreeMap::new(); m];
    let mut inv_perm = vec![elmer_core::NO_DOF; m];
    let mut out = BasisEval::default();
    let nb = mesh.n_bulk();
    for t in 0..mesh.n_boundary() {
        let el = mesh.element(ElemId((nb + t) as u32));
        if bc_of(el.tag) != Some(this) {
            continue;
        }
        // Serial: every parent is active, so the weight is 1 (0/0 without
        // parents, as in the Fortran).
        let pos = el.parents.iter().flatten().count();
        let coeff = if pos == 0 { f64::NAN } else { 1.0 };
        let nodes = el.nodes;
        if nodes.iter().all(|&i| discont_perm[i as usize] == 0) {
            continue;
        }
        let copy = |i: u32| -> u32 {
            match discont_perm[i as usize] {
                0 => i,
                j => (n_orig + j as usize - 1) as u32,
            }
        };
        let def = element_type(el.etype.0).map_err(fe)?;
        let rule = gauss_points_rel(el.etype, 0).map_err(fe)?;
        let xyz: Vec<[f64; 3]> = nodes.iter().map(|&i| mesh.coord(i)).collect();
        let mut st = MetricStatus::default();
        for j in 0..rule.s.len() {
            let uvw = [rule.u[j], rule.v[j], rule.w[j]];
            element_info(def, &xyz, uvw, coord_dim, &mut st, &mut out).map_err(fe)?;
            let weight = coeff * out.det_j * rule.s[j];
            for (p, &np) in nodes.iter().enumerate() {
                let indp = discont_perm[np as usize];
                if indp == 0 {
                    continue;
                }
                let indp = indp as usize - 1;
                let val = weight * out.basis[p];
                inv_perm[indp] = np;
                let row = &mut list[indp];
                if set_diag {
                    add(row, np, val);
                    add(row, copy(np), scale * val);
                } else {
                    for (q, &nq) in nodes.iter().enumerate() {
                        if discont_perm[nq as usize] == 0 {
                            continue;
                        }
                        add(row, nq, out.basis[q] * val);
                        add(row, copy(nq), scale * out.basis[q] * val);
                    }
                }
            }
        }
    }
    log.info(
        CALLER,
        &format!("Created projector for {m} discontinuous nodes"),
        10,
    );
    // List_ToCRSMatrix keeps the rows up to the last one touched.
    let n_rows = list.iter().rposition(|r| !r.is_empty()).map_or(0, |i| i + 1);
    if n_rows == 0 {
        return Ok(None);
    }
    let mut rows = vec![0usize];
    let mut cols = Vec::new();
    let mut values = Vec::new();
    for row in list.into_iter().take(n_rows) {
        for (c, v) in row {
            cols.push(c);
            values.push(v);
        }
        rows.push(cols.len());
    }
    inv_perm.truncate(n_rows);
    Ok(Some(Projector {
        rows,
        cols,
        values,
        inv_perm,
        projector_bc: this,
        galerkin: true,
        slave_scale: 1.0,
        master_scale: 1.0,
        mortar_perm: false,
        rhs: None,
    }))
}
