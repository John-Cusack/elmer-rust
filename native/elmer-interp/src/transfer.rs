//! Field transfer: direct interpolation, projector construction and
//! projector application (the second half of `InterpolateMeshToMeshQ`).

use elmer_core::{BasisEval, ElemId, ElementBasis, Logger, MeshView, NO_DOF};

use crate::field::{DgIndexes, FieldData, Location};
use crate::locate::{BodyMask, Located, Locator};
use crate::shape::Shape;
use crate::InterpError;

const CALLER: &str = "InterpolateMeshToMesh";
/// `Info` without an explicit level is level 4.
const DEF_LEVEL: u8 = 4;

/// Optional inputs of [`interpolate`] that most callers leave empty.
#[derive(Default)]
pub struct Extras<'a> {
    /// `MaskName`: keep only bodies for which this returns true.
    pub body_mask: Option<BodyMask<'a>>,
    /// `NewMaskPerm`: skip target nodes whose entry is [`NO_DOF`].
    pub new_mask_perm: Option<&'a [u32]>,
    /// DG indexes of the old mesh (needed for `NodesOnElements` sources).
    pub old_dg: Option<&'a DgIndexes>,
    /// DG indexes of the new mesh (needed for `NodesOnElements` targets).
    pub new_dg: Option<&'a DgIndexes>,
}

/// Outcome of one transfer.
#[derive(Clone, Debug)]
pub struct InterpReport {
    /// `FoundNodes`: target nodes located in the old mesh.
    pub found: Vec<bool>,
    pub n_not_found: usize,
    /// Nodes the quadrant tree missed but the full scan found.
    pub qtree_fails: usize,
}

/// Where a located target node lands in the old mesh; `None` for nodes that
/// were skipped by `NewMaskPerm` or not found.
fn locate_all(
    loc: &Locator<'_>,
    new_mesh: &dyn MeshView,
    extras: &Extras<'_>,
    logger: &dyn Logger,
) -> Result<(Vec<Option<Located>>, InterpReport), InterpError> {
    let n = new_mesh.n_nodes();
    let mut hits = vec![None; n];
    let mut found = vec![false; n];
    let mut tot_fails = 0usize;
    let start_fails = loc.qtree_fails();
    for i in 0..n {
        if let Some(m) = extras.new_mask_perm {
            if m[i] == NO_DOF {
                continue;
            }
        }
        let p = loc.search_point(new_mesh.coord(i as u32));
        match loc.locate(p, extras.body_mask) {
            Some(l) => {
                hits[i] = Some(l);
                found[i] = true;
            }
            None => {
                logger.info(
                    CALLER,
                    &format!(
                        "Point {} at {:10.2E}{:10.2E}{:10.2E} not found!",
                        i + 1,
                        p[0],
                        p[1],
                        p[2]
                    ),
                    30,
                );
                tot_fails += 1;
            }
        }
    }
    let qtree_fails = loc.qtree_fails() - start_fails;
    if qtree_fails > 0 {
        logger.info(
            CALLER,
            &format!("Number of points not found in quadtree: {qtree_fails}"),
            DEF_LEVEL,
        );
        if tot_fails == 0 {
            logger.info(
                CALLER,
                "All nodes still found by N^2 dummy search!",
                DEF_LEVEL,
            );
        }
    }
    if tot_fails == 0 {
        logger.info(CALLER, "Found all nodes in the target mesh", 6);
    } else {
        logger.warn(
            CALLER,
            &format!("Points not found: {} (found {})", tot_fails, n - tot_fails),
        );
        if loc.options().abort_not_found {
            return Err(InterpError::NotFound(tot_fails));
        }
    }
    Ok((
        hits,
        InterpReport {
            n_not_found: found.iter().filter(|f| !**f).count(),
            found,
            qtree_fails,
        },
    ))
}

/// `Sum_n N_n(uvw) f_n`, the `InterpolateInElement` evaluation.
fn interpolate_in_element(shape: &Shape, f: &[f64], uvw: [f64; 3]) -> f64 {
    let mut basis = vec![0.0; shape.n_nodes()];
    shape.basis(uvw, &mut basis);
    basis.iter().zip(f).map(|(b, f)| b * f).sum()
}

/// `CreateOneDGIndex`: one representative DG index per target node.
fn one_dg_index(new_mesh: &dyn MeshView, mask: Option<&[u32]>, dg: &DgIndexes) -> Vec<u32> {
    let mut out = vec![NO_DOF; new_mesh.n_nodes()];
    let total = new_mesh.n_bulk() + new_mesh.n_boundary();
    for t in 0..total {
        let e = ElemId(t as u32);
        let el = new_mesh.element(e);
        if let Some(m) = mask {
            if el.nodes.iter().any(|&n| m[n as usize] == NO_DOF) {
                continue;
            }
        }
        let own = dg.get(e);
        let p: Vec<u32> = if !own.is_empty() {
            own.to_vec()
        } else if t >= new_mesh.n_bulk() {
            let mut q = vec![NO_DOF; el.nodes.len()];
            for parent in el.parents.iter().flatten() {
                let pe = new_mesh.element(*parent);
                let pdg = dg.get(*parent);
                for (i, &nd) in el.nodes.iter().enumerate() {
                    if q[i] != NO_DOF {
                        continue;
                    }
                    if let Some(j) = pe.nodes.iter().position(|&x| x == nd) {
                        if let Some(&d) = pdg.get(j) {
                            q[i] = d;
                        }
                    }
                }
            }
            q
        } else {
            continue;
        };
        for (i, &nd) in el.nodes.iter().enumerate() {
            if out[nd as usize] == NO_DOF {
                out[nd as usize] = p[i];
            }
        }
    }
    out
}

/// Interpolate every legitimate old-mesh variable onto the same-named target
/// variable, node by node. Target values at unfound nodes are left alone.
pub fn interpolate(
    loc: &Locator<'_>,
    new_mesh: &dyn MeshView,
    old_fields: &[FieldData],
    new_fields: &mut [FieldData],
    extras: &Extras<'_>,
    logger: &dyn Logger,
) -> Result<InterpReport, InterpError> {
    let old_mesh = loc.mesh();
    if old_mesh.n_nodes() == 0 {
        return Ok(InterpReport {
            found: vec![false; new_mesh.n_nodes()],
            n_not_found: new_mesh.n_nodes(),
            qtree_fails: 0,
        });
    }
    let (hits, report) = locate_all(loc, new_mesh, extras, logger)?;
    let partial = loc.options().partial_hit;

    // Pair each legitimate old variable with its target by name.
    let pairs: Vec<(usize, usize)> = old_fields
        .iter()
        .enumerate()
        .filter(|(_, v)| v.is_legit())
        .filter_map(|(oi, v)| {
            new_fields
                .iter()
                .position(|n| n.name == v.name)
                .map(|ni| (oi, ni))
        })
        .collect();

    let mut one_dg: Option<Vec<u32>> = None;
    let mut nodal_xyz: Vec<f64> = Vec::new();
    for (i, hit) in hits.iter().enumerate() {
        let Some(hit) = hit else { continue };
        let el = old_mesh.element(hit.element);
        let shape = Shape::get(el.etype.0).ok_or(InterpError::UnsupportedElement(el.etype.0))?;
        for &(oi, ni) in &pairs {
            let old = &old_fields[oi];
            let new = &mut new_fields[ni];
            let indexes: &[u32] = if old.location == Location::NodesOnElements {
                extras
                    .old_dg
                    .ok_or(InterpError::MissingDgIndexes("old"))?
                    .get(hit.element)
            } else {
                el.nodes
            };
            let valid = indexes.iter().filter(|&&ix| old.slot(ix).is_some()).count();
            if !(valid == indexes.len() || (partial && valid > 0)) {
                continue;
            }
            let dest = if new.location == Location::NodesOnElements {
                let new_dg = extras.new_dg.ok_or(InterpError::MissingDgIndexes("new"))?;
                let table = one_dg
                    .get_or_insert_with(|| one_dg_index(new_mesh, extras.new_mask_perm, new_dg));
                match table[i] {
                    NO_DOF => None,
                    d => new.slot(d),
                }
            } else {
                new.slot(i as u32)
            };
            let Some(k) = dest else { continue };
            let n = shape.n_nodes();
            nodal_xyz.clear();
            nodal_xyz.extend(
                indexes
                    .iter()
                    .take(n)
                    .map(|&ix| old.slot(ix).map_or(0.0, |s| old.values[s])),
            );
            new.values[k] = interpolate_in_element(shape, &nodal_xyz, hit.uvw);
            for j in 0..old.prev_values.len() {
                nodal_xyz.clear();
                nodal_xyz.extend(
                    indexes
                        .iter()
                        .take(n)
                        .map(|&ix| old.slot(ix).map_or(0.0, |s| old.prev_values[j][s])),
                );
                let v = interpolate_in_element(shape, &nodal_xyz, hit.uvw);
                if let Some(col) = new.prev_values.get_mut(j) {
                    col[k] = v;
                }
            }
        }
    }
    Ok(report)
}

/// Sparse interpolation operator `new = P * old` on node numbering
/// (`Projector_t`): row `i` holds the normalised element basis values at
/// target node `i`, columns are old-mesh node indices (0-based, in element
/// order, not sorted).
#[derive(Clone, Debug)]
pub struct Projector {
    pub rows: Vec<usize>,
    pub cols: Vec<u32>,
    pub values: Vec<f64>,
    /// Transpose (rows per old node, columns ascending); `None` when no
    /// target node was found.
    pub t_rows: Vec<usize>,
    pub t_cols: Vec<u32>,
    pub t_values: Vec<f64>,
    pub found: Vec<bool>,
}

impl Projector {
    pub fn n_rows(&self) -> usize {
        self.rows.len() - 1
    }

    /// `CRS_ApplyProjector` for every variable present in both sets.
    pub fn apply(&self, old_fields: &[FieldData], new_fields: &mut [FieldData]) {
        for old in old_fields.iter().filter(|v| v.is_legit()) {
            let Some(new) = new_fields.iter_mut().find(|n| n.name == old.name) else {
                continue;
            };
            self.apply_one(
                &old.values,
                old.perm.as_deref(),
                &mut new.values,
                new.perm.as_deref(),
            );
            for (j, col) in old.prev_values.iter().enumerate() {
                if let Some(dst) = new.prev_values.get_mut(j) {
                    self.apply_one(col, old.perm.as_deref(), dst, new.perm.as_deref());
                }
            }
        }
    }

    fn apply_one(&self, u: &[f64], uperm: Option<&[u32]>, v: &mut [f64], vperm: Option<&[u32]>) {
        let n = self.n_rows();
        let span = |i: usize| self.rows[i]..self.rows[i + 1];
        match (uperm, vperm) {
            (Some(up), Some(vp)) => {
                for i in 0..n {
                    if vp[i] != NO_DOF && self.values[span(i)].iter().any(|&x| x != 0.0) {
                        v[vp[i] as usize] = 0.0;
                    }
                }
                for i in 0..n {
                    if vp[i] == NO_DOF {
                        continue;
                    }
                    for j in span(i) {
                        let k = up[self.cols[j] as usize];
                        if k != NO_DOF {
                            v[vp[i] as usize] += u[k as usize] * self.values[j];
                        }
                    }
                }
            }
            _ => {
                for (i, vi) in v.iter_mut().enumerate().take(n) {
                    for j in span(i) {
                        *vi += u[self.cols[j] as usize] * self.values[j];
                    }
                }
            }
        }
    }
}

/// Locate every target node and build the projector from the caller's
/// element basis (`ElementInfo`).
pub fn build_projector(
    loc: &Locator<'_>,
    new_mesh: &dyn MeshView,
    basis: &dyn ElementBasis,
    extras: &Extras<'_>,
    logger: &dyn Logger,
) -> Result<Projector, InterpError> {
    let old_mesh = loc.mesh();
    let (hits, report) = locate_all(loc, new_mesh, extras, logger)?;
    let keep = loc.options().keep_unfound_nodes;
    let mut rows = vec![0usize];
    let mut cols: Vec<u32> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    let mut out = BasisEval::default();
    let mut found_cnt = 0usize;
    let entries: Vec<Option<&Located>> = if keep {
        hits.iter().map(|h| h.as_ref()).collect()
    } else {
        hits.iter().filter_map(|h| h.as_ref()).map(Some).collect()
    };
    for hit in &entries {
        if let Some(hit) = hit {
            found_cnt += 1;
            let el = old_mesh.element(hit.element);
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| old_mesh.coord(n)).collect();
            basis
                .eval(el.etype, &xyz, hit.uvw, &mut out)
                .map_err(InterpError::Basis)?;
            let np = el.nodes.len();
            let rowsum: f64 = out.basis[..np].iter().sum();
            for j in 0..np {
                cols.push(el.nodes[j]);
                values.push(out.basis[j] / rowsum);
            }
        }
        rows.push(cols.len());
    }
    let (t_rows, t_cols, t_values) = if found_cnt > 0 {
        transpose(&rows, &cols, &values, old_mesh.n_nodes())
    } else {
        (Vec::new(), Vec::new(), Vec::new())
    };
    Ok(Projector {
        rows,
        cols,
        values,
        t_rows,
        t_cols,
        t_values,
        found: report.found,
    })
}

/// Transpose with the Fortran's row count `max(n_old, max col)`.
fn transpose(
    rows: &[usize],
    cols: &[u32],
    vals: &[f64],
    n_old: usize,
) -> (Vec<usize>, Vec<u32>, Vec<f64>) {
    let n = n_old.max(cols.iter().map(|&c| c as usize + 1).max().unwrap_or(0));
    let mut count = vec![0usize; n];
    for &c in cols {
        count[c as usize] += 1;
    }
    let mut t_rows = vec![0usize; n + 1];
    for i in 0..n {
        t_rows[i + 1] = t_rows[i] + count[i];
    }
    let mut fill = vec![0usize; n];
    let mut t_cols = vec![0u32; cols.len()];
    let mut t_vals = vec![0.0; cols.len()];
    for i in 0..rows.len() - 1 {
        for j in rows[i]..rows[i + 1] {
            let k = cols[j] as usize;
            let l = t_rows[k] + fill[k];
            fill[k] += 1;
            t_cols[l] = i as u32;
            t_vals[l] = vals[j];
        }
    }
    (t_rows, t_cols, t_vals)
}
