//! Geometric multigrid (`Linear System Solver = Multigrid`, `MG
//! Method = geometric`, `MG Equal Split`): `GMGSolve`/`GMGSweep`
//! (Multigrid.F90:136-660) and `MGSmooth` (Smoothers.F90:51) over the mesh
//! hierarchy the caller hands over ([`crate::Simulation::set_multigrid_meshes`]:
//! the `SplitMeshEqual` parents of the simulation mesh, which is the
//! finest level, MainUtils.F90:2730-2770).
//!
//! The coarse permutations follow the solver's (`Optimize Bandwidth`), the
//! transfer operators are `MeshProjector(Mesh2, Mesh1)`'s linear
//! interpolation (triangles only so far) and the coarse matrices the
//! Galerkin product of `ProjectMatrix` ([`galerkin`]; its summation order is
//! not reproduced). Refused: other `MG Method`s, `MG Project Matrix = False`
//! (reassembly through `MSolverActivate`), smoothers other than `jacobi`,
//! `sgs` and `bicgstab`, `MG Lowest Linear Solver` other than a direct
//! solve, vector fields, p elements and parallel runs. Without the coarse
//! meshes the solve is refused.

use elmer_core::{CsrMatrix, DirectSolver, ElemId, ElementType, LinearOperator, Logger, MeshView, ParamMap, Params, Preconditioner, NO_DOF};
use elmer_sparse::Ilu;

use crate::error::SolverError;
use crate::sim::Simulation;

const CALLER: &str = "GMGSolve";

/// One level, finest last. `interp` maps this level's coarse neighbour
/// (`index - 1`) onto it; the coarsest has none.
pub(crate) struct Level {
    pub a: CsrMatrix<f64>,
    pub interp: Option<CsrMatrix<f64>>,
    ilu: Option<Ilu<f64>>,
}

/// The `MG*` keywords `GMGSolve` and `MGSmooth` read.
#[derive(Clone, Debug)]
pub(crate) struct MgOptions {
    /// `MG Max Iterations`, else `Linear System Max Iterations`.
    pub max_iter: usize,
    /// `MG Convergence Tolerance`, else `Linear System Convergence Tolerance`.
    pub tolerance: f64,
    /// `MG Level Max Iterations` (1) and `MG Level Convergence Tolerance`.
    pub level_max_iter: usize,
    pub level_tolerance: f64,
    /// `MG Sweeps` per level.
    pub sweeps: Vec<usize>,
    /// `MG Pre/Post Smoothing Iterations`, indexed by `1 + total - level`.
    pub pre: Vec<usize>,
    pub post: Vec<usize>,
    /// `MG Smoother`, else `Linear System Iterative Method`, else `sgs`.
    pub smoother: String,
}

/// The coarsest-level solve (`MG Lowest Linear Solver = direct`).
pub(crate) type CoarseSolve<'a> = dyn Fn(&CsrMatrix<f64>, &[f64], &mut [f64]) -> Result<(), SolverError> + 'a;

/// Which multigrid's cycle: `GMGSolve` normalises by the right-hand side,
/// tests the top level against the relative residual and prolongs onto the
/// residual (`CRS_ApplyProjector`); `AMGSolve` does none of these
/// (`CRS_ProjectVector` zeroes its output).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Kind {
    Geometric,
    Algebraic,
}

pub(crate) struct Hierarchy {
    pub levels: Vec<Level>,
    pub kind: Kind,
    /// Components per node, for the `bsgs` smoother.
    pub dofs: usize,
}

/// `B = I^T A I`: the coarse operator of `ProjectMatrix` as a Galerkin
/// product (rows and columns in increasing order).
pub(crate) fn galerkin(a: &CsrMatrix<f64>, interp: &CsrMatrix<f64>) -> CsrMatrix<f64> {
    let nc = interp.ncols;
    // AI = A * I, row by row.
    let mut ai: Vec<std::collections::BTreeMap<u32, f64>> = vec![Default::default(); a.nrows];
    for (i, row) in ai.iter_mut().enumerate() {
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            let (j, v) = (a.col_idx[k] as usize, a.values[k]);
            for l in interp.row_ptr[j]..interp.row_ptr[j + 1] {
                *row.entry(interp.col_idx[l]).or_default() += v * interp.values[l];
            }
        }
    }
    // B = I^T * AI.
    let mut b: Vec<std::collections::BTreeMap<u32, f64>> = vec![Default::default(); nc];
    for (i, row) in ai.iter().enumerate() {
        for l in interp.row_ptr[i]..interp.row_ptr[i + 1] {
            let (c, w) = (interp.col_idx[l] as usize, interp.values[l]);
            for (&j, &v) in row {
                *b[c].entry(j).or_default() += w * v;
            }
        }
    }
    let mut m = CsrMatrix { nrows: nc, ncols: nc, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
    for row in b {
        for (j, v) in row {
            m.col_idx.push(j);
            m.values.push(v);
        }
        m.row_ptr.push(m.values.len());
    }
    m
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}

fn residual(a: &CsrMatrix<f64>, x: &[f64], b: &[f64], r: &mut [f64]) {
    a.apply(x, r);
    for (ri, bi) in r.iter_mut().zip(b) {
        *ri = bi - *ri;
    }
}

/// Fortran `Ew.dEe` (mantissa in [0.1, 1)), right-aligned in `w`.
pub(crate) fn e_fmt(x: f64, w: usize, d: usize, e: usize) -> String {
    let body = if x == 0.0 {
        format!("0.{}E+{}", "0".repeat(d), "0".repeat(e))
    } else {
        let s = format!("{:.*e}", d - 1, x.abs());
        let (m, ex) = s.split_once('e').expect("exponent");
        let ex: i32 = ex.parse().expect("exponent");
        let digits: String = m.chars().filter(|c| *c != '.').collect();
        let ex = ex + 1;
        let sign = if ex < 0 { '-' } else { '+' };
        format!("{}0.{}E{}{:0>e$}", if x < 0.0 { "-" } else { "" }, digits, sign, ex.abs(), e = e)
    };
    format!("{body:>w$}")
}

impl Hierarchy {
    /// Builds the ILU(0) factors the `bicgstab` smoother uses on every
    /// level above the coarsest (`CRS_IncompleteLU`, GMGSolve).
    pub(crate) fn new(levels: Vec<(CsrMatrix<f64>, Option<CsrMatrix<f64>>)>, ilu: Option<usize>) -> Result<Self, SolverError> {
        let mut out = Vec::with_capacity(levels.len());
        for (k, (a, interp)) in levels.into_iter().enumerate() {
            let ilu = match ilu {
                Some(level) if k > 0 => Some(Ilu::ilu_n(&a, level)?),
                _ => None,
            };
            out.push(Level { a, interp, ilu });
        }
        Ok(Self { levels: out, kind: Kind::Geometric, dofs: 1 })
    }

    fn total(&self) -> usize {
        self.levels.len()
    }

    /// `GMGSolve` at `level` (1 = coarsest, `total` = finest).
    pub(crate) fn solve(
        &self,
        level: usize,
        x: &mut [f64],
        b: &mut [f64],
        opts: &MgOptions,
        coarse: &CoarseSolve<'_>,
        log: &dyn Logger,
    ) -> Result<(), SolverError> {
        let lv = &self.levels[level - 1];
        let n = lv.a.nrows;
        let geometric = self.kind == Kind::Geometric;
        let mut rhs_norm = norm(&b[..n]);
        if rhs_norm < f64::MIN_POSITIVE.sqrt() {
            if geometric {
                log.info(CALLER, "Right-hand side is almost zero, not normalizing by it", 8);
            }
            rhs_norm = 1.0;
        }
        // AMGSolve leaves the system as it is; GMGSolve normalises it.
        let scale = if geometric { rhs_norm } else { 1.0 };
        x[..n].iter_mut().for_each(|v| *v /= scale);
        b[..n].iter_mut().for_each(|v| *v /= scale);
        if level <= 1 {
            log.info(CALLER, "Starting lowest linear solver using: direct", 10);
            coarse(&lv.a, &b[..n], &mut x[..n])?;
        } else {
            let (max_iter, tol) = if level == self.total() {
                (opts.max_iter, opts.tolerance)
            } else {
                (opts.level_max_iter, opts.level_tolerance)
            };
            let mut r = vec![0.0; n];
            residual(&lv.a, x, b, &mut r);
            for iter in 1..=max_iter {
                let rnorm = self.sweep(level, x, b, &mut r, opts, coarse, log)?;
                let done = if geometric {
                    log.info(
                        CALLER,
                        &format!(
                            "MG Residual at level: {level} iter: {iter} is:{}{}",
                            e_fmt(rnorm / rhs_norm, 20, 12, 3),
                            e_fmt(rnorm, 20, 12, 3)
                        ),
                        5,
                    );
                    if rnorm.is_nan() || rnorm > 1.0e50 {
                        return Err(SolverError::fatal(CALLER, "We seem to have diverged"));
                    }
                    if level == self.total() {
                        rnorm / rhs_norm < tol
                    } else {
                        rnorm < tol
                    }
                } else {
                    // AMGSolve: the relative residual at every level.
                    let rel = rnorm / rhs_norm;
                    log.info(
                        "AMGSolve",
                        &format!("MG Residual at level: {level} iter: {iter} is:{}", e_fmt(rel, 20, 12, 3)),
                        5,
                    );
                    if rel.is_nan() || rel > 1.0e50 {
                        return Err(SolverError::fatal("AMGSolve", "We seem to have diverged"));
                    }
                    rel < tol
                };
                if done {
                    break;
                }
            }
        }
        x[..n].iter_mut().for_each(|v| *v *= scale);
        b[..n].iter_mut().for_each(|v| *v *= scale);
        Ok(())
    }

    /// `GMGSweep`: pre-smooth, restrict the residual, solve the coarse
    /// correction `MG Sweeps` times from zero, prolong it onto the residual
    /// and add both to the solution (as the Fortran does), post-smooth. Returns the residual norm after post-smoothing.
    #[allow(clippy::too_many_arguments)]
    fn sweep(
        &self,
        level: usize,
        x: &mut [f64],
        b: &[f64],
        r: &mut [f64],
        opts: &MgOptions,
        coarse: &CoarseSolve<'_>,
        log: &dyn Logger,
    ) -> Result<f64, SolverError> {
        let lv = &self.levels[level - 1];
        let interp = lv.interp.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "Matrix2 not associated!"))?;
        let sweeps = opts.sweeps.get(level.min(opts.sweeps.len()).saturating_sub(1)).copied().unwrap_or(1);
        self.smooth(level, x, b, r, opts, true)?;
        let nc = interp.ncols;
        let mut r2 = vec![0.0; nc];
        for i in 0..interp.nrows {
            for k in interp.row_ptr[i]..interp.row_ptr[i + 1] {
                r2[interp.col_idx[k] as usize] += interp.values[k] * r[i];
            }
        }
        let mut x2 = vec![0.0; nc];
        for _ in 0..sweeps {
            self.solve(level - 1, &mut x2, &mut r2, opts, coarse, log)?;
        }
        // CRS_ApplyProjector(Trans = .TRUE.) adds the prolonged correction
        // into Residual, which still holds the pre-smoothed residual, and
        // GMGSweep adds all of it: x += r + I x2.
        let mut e = vec![0.0; r.len()];
        interp.apply(&x2, &mut e);
        let adds_residual = self.kind == Kind::Geometric;
        for ((xi, ri), ei) in x.iter_mut().zip(r.iter_mut()).zip(&e) {
            *ri = if adds_residual { *ri + ei } else { *ei };
            *xi += *ri;
        }
        self.smooth(level, x, b, r, opts, false)
    }

    /// `MGSmooth` with the `jacobi`, `sgs` and `bicgstab` smoothers; the
    /// residual `b - A x` is left in `r` and its norm returned.
    fn smooth(&self, level: usize, x: &mut [f64], b: &[f64], r: &mut [f64], opts: &MgOptions, pre: bool) -> Result<f64, SolverError> {
        let lv = &self.levels[level - 1];
        let a = &lv.a;
        let n = a.nrows;
        let inv_level = (1 + self.total()).saturating_sub(level).max(1);
        let pick = |v: &[usize]| v.get(inv_level.min(v.len()).saturating_sub(1)).copied().unwrap_or(1);
        let rounds = if pre { pick(&opts.pre) } else { pick(&opts.post) };
        let inv_diag: Vec<f64> = (0..n)
            .map(|i| {
                let d = (a.row_ptr[i]..a.row_ptr[i + 1]).find(|&k| a.col_idx[k] as usize == i).map_or(0.0, |k| a.values[k]);
                if d != 0.0 {
                    1.0 / d
                } else {
                    0.0
                }
            })
            .collect();
        match opts.smoother.as_str() {
            _ if rounds == 0 => {}
            "jacobi" => {
                for _ in 0..rounds {
                    a.apply(x, r);
                    for j in 0..n {
                        r[j] = b[j] - r[j];
                        x[j] += r[j] * inv_diag[j];
                    }
                }
            }
            "sgs" => {
                let row = |x: &[f64], i: usize| (a.row_ptr[i]..a.row_ptr[i + 1]).map(|k| x[a.col_idx[k] as usize] * a.values[k]).sum::<f64>();
                for _ in 0..rounds {
                    for i in 0..n {
                        r[i] = (b[i] - row(x, i)) * inv_diag[i];
                        x[i] += r[i];
                    }
                    for i in (0..n).rev() {
                        r[i] = (b[i] - row(x, i)) * inv_diag[i];
                        x[i] += r[i];
                    }
                }
            }
            "bsgs" => {
                // BSGS (Smoothers.F90:1054): a node's components together.
                let d = self.dofs.max(1);
                let rowsum = |x: &[f64], i: usize| (a.row_ptr[i]..a.row_ptr[i + 1]).map(|k| x[a.col_idx[k] as usize] * a.values[k]).sum::<f64>();
                for _ in 0..rounds {
                    for node in (0..n / d).chain((0..n / d).rev()) {
                        let s: Vec<f64> = (0..d).map(|c| rowsum(x, node * d + c)).collect();
                        for (c, sc) in s.iter().enumerate() {
                            let id = node * d + c;
                            r[id] = (b[id] - sc) * inv_diag[id];
                            x[id] += r[id];
                        }
                    }
                }
            }
            "bicgstab" => {
                let m = lv.ilu.as_ref().ok_or_else(|| SolverError::fatal("MGSmooth", "no ILU factors for the BiCGStab smoother"))?;
                residual(a, x, b, r);
                let ri = r.to_vec();
                let (mut pr, mut v, mut t) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
                let (mut t1, mut t2, mut s) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
                let (mut omega, mut alpha, mut oldrho) = (1.0, 0.0, 1.0);
                for _ in 0..rounds {
                    let rho = dot(r, &ri);
                    let beta = alpha * rho / (oldrho * omega);
                    for j in 0..n {
                        pr[j] = r[j] + beta * (pr[j] - omega * v[j]);
                    }
                    m.apply(&pr, &mut t1);
                    a.apply(&t1, &mut v);
                    alpha = rho / dot(&ri, &v);
                    for j in 0..n {
                        s[j] = r[j] - alpha * v[j];
                    }
                    m.apply(&s, &mut t2);
                    a.apply(&t2, &mut t);
                    omega = dot(&t, &s) / dot(&t, &t);
                    oldrho = rho;
                    for j in 0..n {
                        r[j] = s[j] - omega * t[j];
                        x[j] += alpha * t1[j] + omega * t2[j];
                    }
                }
            }
            other => return Err(SolverError::Unsupported(format!("MG Smoother = {other}"))),
        }
        residual(a, x, b, r);
        Ok(norm(r))
    }
}

/// Bulk elements of `mesh` on bodies where solver `id` is active.
fn active_on(sim: &Simulation, id: usize, mesh: &dyn MeshView) -> Vec<ElemId> {
    (0..mesh.n_bulk() as u32)
        .map(ElemId)
        .filter(|&e| {
            sim.model
                .body_of_tag(mesh.element(e).tag)
                .is_some_and(|b| sim.model.solver_active_on_body(id, b))
        })
        .collect()
}

/// `MeshProjector(Mesh2, Mesh1)` as linear interpolation: for each node of
/// `fine`, the weights of the `coarse` triangle that contains it.
fn node_projector(coarse: &dyn MeshView, active: &[ElemId], fine: &dyn MeshView) -> Result<Vec<Vec<(u32, f64)>>, SolverError> {
    let mut tris = Vec::with_capacity(active.len());
    for &e in active {
        let el = coarse.element(e);
        if el.etype != ElementType(303) {
            return Err(SolverError::Unsupported(format!("multigrid MeshProjector for {:?} elements", el.etype)));
        }
        let p: Vec<[f64; 3]> = el.nodes.iter().map(|&n| coarse.coord(n)).collect();
        let (lo, hi) = (
            [p[0][0].min(p[1][0]).min(p[2][0]), p[0][1].min(p[1][1]).min(p[2][1])],
            [p[0][0].max(p[1][0]).max(p[2][0]), p[0][1].max(p[1][1]).max(p[2][1])],
        );
        tris.push((el.nodes.to_vec(), p, lo, hi));
    }
    let eps = 1.0e-10;
    let mut out = Vec::with_capacity(fine.n_nodes());
    for n in 0..fine.n_nodes() as u32 {
        let q = fine.coord(n);
        let mut found = None;
        for (nodes, p, lo, hi) in &tris {
            let tol = eps * (hi[0] - lo[0]).max(hi[1] - lo[1]);
            if q[0] < lo[0] - tol || q[0] > hi[0] + tol || q[1] < lo[1] - tol || q[1] > hi[1] + tol {
                continue;
            }
            let det = (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]);
            let l1 = ((q[0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (q[1] - p[0][1])) / det;
            let l2 = ((p[1][0] - p[0][0]) * (q[1] - p[0][1]) - (q[0] - p[0][0]) * (p[1][1] - p[0][1])) / det;
            let w = [1.0 - l1 - l2, l1, l2];
            if w.iter().all(|&x| x >= -eps) {
                found = Some(nodes.iter().copied().zip(w).filter(|(_, x)| x.abs() > eps).collect());
                break;
            }
        }
        out.push(found.unwrap_or_default());
    }
    Ok(out)
}

/// The interpolation operators (level 2 up, fine dofs x coarse dofs) of a
/// `Linear System Solver = Multigrid` solve of slot `ix`, or `None` for
/// another linear solver.
pub(crate) fn transfers(sim: &Simulation, ix: usize, mesh: &dyn MeshView) -> Result<Option<Vec<CsrMatrix<f64>>>, SolverError> {
    let slot = &sim.slots[ix];
    let p = &slot.params;
    if !p.string("linear system solver").is_some_and(|s| elmer_core::normalize_key(s) == "multigrid") {
        return Ok(None);
    }
    if let Some(m) = p.string("mg method").map(elmer_core::normalize_key) {
        if m != "geometric" {
            return Err(SolverError::Unsupported(format!("MG Method = {m}")));
        }
    }
    for k in ["mg algebraic", "mg cluster", "mg pelement"] {
        if p.logical(k) == Some(true) {
            return Err(SolverError::Unsupported(format!("{k} multigrid")));
        }
    }
    if p.logical("mg equal split") != Some(true) {
        return Err(SolverError::Unsupported("geometric multigrid without MG Equal Split".into()));
    }
    if p.logical("mg project matrix") == Some(false) {
        return Err(SolverError::Unsupported("MG Project Matrix = False".into()));
    }
    if p.get("mg lowest linear solver").is_some_and(|v| !matches!(v, elmer_core::Value::Str(s) if elmer_core::normalize_key(s) == "direct"))
        || p.logical("mg lowest linear solver iterative") == Some(true)
    {
        return Err(SolverError::Unsupported("MG Lowest Linear Solver other than direct".into()));
    }
    let levels = p
        .integer("mg levels")
        .or_else(|| p.integer("multigrid levels"))
        .ok_or_else(|| SolverError::fatal("AddEquationBasics", "> MG Levels < must be defined for geometric multigrid!"))?
        .max(1) as usize;
    let id = slot.variable.ok_or_else(|| SolverError::fatal("MultiGridSolve", "solver has no variable"))?;
    let var = sim.vars.by_id(id);
    if var.dofs != 1 || slot.dofs.is_some() {
        return Err(SolverError::Unsupported("geometric multigrid of a vector field or p elements".into()));
    }
    let coarse = slot
        .mg_meshes
        .as_ref()
        .ok_or_else(|| SolverError::Unsupported("Linear System Solver = Multigrid without its mesh hierarchy".into()))?;
    if coarse.len() + 1 != levels {
        return Err(SolverError::fatal(
            "MultiGridSolve",
            format!("{} coarse meshes for MG Levels = {levels}", coarse.len()),
        ));
    }
    let fine_perm: Vec<u32> = var.perm.as_deref().map(<[u32]>::to_vec).unwrap_or_default();
    let mut perms = Vec::with_capacity(levels);
    let mut meshes: Vec<&dyn MeshView> = coarse.iter().map(|m| m.as_ref()).collect();
    meshes.push(mesh);
    for (k, m) in meshes.iter().enumerate() {
        if k + 1 == levels {
            perms.push((fine_perm.clone(), var.values.len()));
            continue;
        }
        let active = active_on(sim, slot.id, *m);
        let perm = match &sim.backends.assembly {
            Some(b) => b.build_perm_for(*m, &active, p)?,
            None => crate::sim::node_order_perm(*m, &active),
        };
        perms.push((perm.perm, perm.n_slots));
    }
    let mut out = Vec::with_capacity(levels - 1);
    for k in 1..levels {
        let active = active_on(sim, slot.id, meshes[k - 1]);
        let proj = node_projector(meshes[k - 1], &active, meshes[k])?;
        let (pf, nf) = (&perms[k].0, perms[k].1);
        let (pc, nc) = (&perms[k - 1].0, perms[k - 1].1);
        let mut rows: Vec<Vec<(u32, f64)>> = vec![Vec::new(); nf];
        for (node, w) in proj.iter().enumerate() {
            let r = pf.get(node).copied().unwrap_or(NO_DOF);
            if r == NO_DOF {
                continue;
            }
            for &(cn, x) in w {
                let c = pc.get(cn as usize).copied().unwrap_or(NO_DOF);
                if c != NO_DOF {
                    rows[r as usize].push((c, x));
                }
            }
        }
        let mut m = CsrMatrix { nrows: nf, ncols: nc, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
        for mut row in rows {
            row.sort_by_key(|&(c, _)| c);
            for (c, x) in row {
                m.col_idx.push(c);
                m.values.push(x);
            }
            m.row_ptr.push(m.values.len());
        }
        out.push(m);
    }
    Ok(Some(out))
}

/// `MultiGridSolve` -> `GMGSolve` at the finest level for the (scaled)
/// system `a x = b`, the coarse matrices projected from `a`.
pub(crate) fn solve(
    p: &ParamMap,
    a: &CsrMatrix<f64>,
    b: &mut [f64],
    x: &mut [f64],
    interps: &[CsrMatrix<f64>],
    log: &dyn Logger,
) -> Result<(), SolverError> {
    let ints = |k: &str| match p.get(k) {
        Some(elmer_core::Value::Integer(v)) => Some(v.iter().map(|&i| i.max(0) as usize).collect::<Vec<_>>()),
        _ => None,
    };
    let multigrid = p.string("linear system solver").is_some_and(|s| elmer_core::normalize_key(s) == "multigrid");
    let max_iter = p
        .integer("mg max iterations")
        .or_else(|| if multigrid { p.integer("linear system max iterations") } else { Some(1) })
        .unwrap_or(1)
        .max(0) as usize;
    let tolerance = p
        .real("mg convergence tolerance")
        .or_else(|| p.real("linear system convergence tolerance"))
        .unwrap_or(0.0);
    let smoother = p
        .string("mg smoother")
        .or_else(|| p.string("linear system iterative method"))
        .map(elmer_core::normalize_key)
        .unwrap_or_else(|| "sgs".into());
    let precond = p
        .string("mg preconditioning")
        .or_else(|| p.string("linear system preconditioning"))
        .map(elmer_core::normalize_key)
        .unwrap_or_else(|| "none".into());
    let ilu = match precond.as_str() {
        s if s.starts_with("ilut") => return Err(SolverError::Unsupported("MG Preconditioning = ILUT".into())),
        s if s.starts_with("ilu") => Some(s[3..].parse::<usize>().ok().filter(|&k| k <= 9).unwrap_or(0)),
        _ => None,
    };
    if smoother == "bicgstab" && ilu.is_none() {
        return Err(SolverError::Unsupported("the BiCGStab multigrid smoother without ILU".into()));
    }
    let opts = MgOptions {
        max_iter,
        tolerance,
        level_max_iter: p.integer("mg level max iterations").map_or(1, |i| i.max(0) as usize),
        level_tolerance: p.real("mg level convergence tolerance").unwrap_or(f64::MAX),
        sweeps: ints("mg sweeps").unwrap_or_else(|| vec![1]),
        pre: ints("mg pre smoothing iterations").unwrap_or_else(|| vec![1]),
        post: ints("mg post smoothing iterations").unwrap_or_else(|| vec![1]),
        smoother,
    };
    // ProjectMatrix down from the finest.
    let mut mats = vec![a.clone()];
    for interp in interps.iter().rev() {
        let c = galerkin(mats.last().expect("finest"), interp);
        mats.push(c);
    }
    mats.reverse();
    let mut levels: Vec<(CsrMatrix<f64>, Option<CsrMatrix<f64>>)> = Vec::with_capacity(mats.len());
    for (k, m) in mats.into_iter().enumerate() {
        levels.push((m, if k == 0 { None } else { Some(interps[k - 1].clone()) }));
    }
    let h = Hierarchy::new(levels, ilu)?;
    log.info("MultiGridSolve", "*********************************", 7);
    log.info("MultiGridSolve", "Performing multigrid solution: geometric", 7);
    let coarse = |a: &CsrMatrix<f64>, b: &[f64], x: &mut [f64]| -> Result<(), SolverError> {
        let mut d = elmer_direct::FaerLu::<f64>::new();
        d.factor(a)?;
        let mut y = b.to_vec();
        d.solve_in_place(&mut y)?;
        x.copy_from_slice(&y);
        Ok(())
    };
    let total = h.levels.len();
    h.solve(total, x, b, &opts, &coarse, log)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `-u'' = 1` on (0, 1) with u(0) = u(1) = 0 on `2^k` elements, the
    /// interior nodes only.
    fn laplace(k: u32) -> CsrMatrix<f64> {
        let n = (1usize << k) - 1;
        let h = 1.0 / (1usize << k) as f64;
        let mut m = CsrMatrix { nrows: n, ncols: n, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
        for i in 0..n {
            if i > 0 {
                m.col_idx.push(i as u32 - 1);
                m.values.push(-1.0 / h);
            }
            m.col_idx.push(i as u32);
            m.values.push(2.0 / h);
            if i + 1 < n {
                m.col_idx.push(i as u32 + 1);
                m.values.push(-1.0 / h);
            }
            m.row_ptr.push(m.values.len());
        }
        m
    }

    /// Linear interpolation from `2^(k-1)` to `2^k` elements.
    fn interp(k: u32) -> CsrMatrix<f64> {
        let (nf, nc) = ((1usize << k) - 1, (1usize << (k - 1)) - 1);
        let mut m = CsrMatrix { nrows: nf, ncols: nc, row_ptr: vec![0], col_idx: Vec::new(), values: Vec::new(), diag: None };
        for i in 0..nf {
            if i % 2 == 1 {
                m.col_idx.push((i / 2) as u32);
                m.values.push(1.0);
            } else {
                if i >= 2 {
                    m.col_idx.push((i / 2 - 1) as u32);
                    m.values.push(0.5);
                }
                if i / 2 < nc {
                    m.col_idx.push((i / 2) as u32);
                    m.values.push(0.5);
                }
            }
            m.row_ptr.push(m.values.len());
        }
        m
    }

    fn dense_solve(a: &CsrMatrix<f64>, b: &[f64], x: &mut [f64]) -> Result<(), SolverError> {
        let n = a.nrows;
        let mut m = vec![0.0; n * n];
        for i in 0..n {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                m[i * n + a.col_idx[k] as usize] = a.values[k];
            }
        }
        x.copy_from_slice(b);
        for c in 0..n {
            for r in c + 1..n {
                let f = m[r * n + c] / m[c * n + c];
                for k in c..n {
                    m[r * n + k] -= f * m[c * n + k];
                }
                x[r] -= f * x[c];
            }
        }
        for c in (0..n).rev() {
            let s: f64 = (c + 1..n).map(|k| m[c * n + k] * x[k]).sum();
            x[c] = (x[c] - s) / m[c * n + c];
        }
        Ok(())
    }

    fn hierarchy(levels: u32, finest: u32, ilu: bool) -> Hierarchy {
        let mut v = Vec::new();
        let coarsest = finest + 1 - levels;
        // Stiffness times h: O(1) entries as in 2D, where the prolongation
        // that also adds the residual stays stable.
        let mut a = laplace(finest);
        let h = 1.0 / (1usize << finest) as f64;
        a.values.iter_mut().for_each(|v| *v *= h);
        let mut stack = vec![(a.clone(), Some(interp(finest)))];
        for k in (coarsest..finest).rev() {
            let p = stack.last().unwrap().1.clone().unwrap();
            a = galerkin(&a, &p);
            stack.push((a.clone(), if k > coarsest { Some(interp(k)) } else { None }));
        }
        while let Some(l) = stack.pop() {
            v.push(l);
        }
        Hierarchy::new(v, ilu.then_some(0)).unwrap()
    }

    fn opts(smoother: &str) -> MgOptions {
        MgOptions {
            max_iter: 50,
            tolerance: 1.0e-10,
            level_max_iter: 1,
            level_tolerance: f64::MAX,
            sweeps: vec![1],
            pre: vec![2],
            post: vec![2],
            smoother: smoother.into(),
        }
    }

    #[test]
    fn galerkin_of_the_1d_laplacian_is_the_coarse_laplacian() {
        let b = galerkin(&laplace(4), &interp(4));
        let c = laplace(3);
        assert_eq!((b.row_ptr.clone(), b.col_idx.clone()), (c.row_ptr.clone(), c.col_idx.clone()));
        for (x, y) in b.values.iter().zip(&c.values) {
            assert!((x - y).abs() < 1e-12);
        }
    }

    #[test]
    fn v_cycles_reach_the_tolerance_and_the_exact_solution() {
        // Plain Jacobi does not damp the highest mode of the 1D Laplacian.
        for smoother in ["sgs", "bicgstab"] {
            let h = hierarchy(4, 6, smoother == "bicgstab");
            let n = h.levels.last().unwrap().a.nrows;
            let mut b = vec![1.0 / (64.0 * 64.0); n];
            let mut x = vec![0.0; n];
            let log = elmer_core::NullLogger;
            h.solve(4, &mut x, &mut b, &opts(smoother), &dense_solve, &log).unwrap();
            // u = x (1 - x) / 2 at the nodes, exact for linear elements.
            for (i, xi) in x.iter().enumerate() {
                let t = (i + 1) as f64 / 64.0;
                assert!((xi - t * (1.0 - t) / 2.0).abs() < 1e-8, "{smoother} node {i}: {xi}");
            }
        }
    }

    #[test]
    fn fortran_e_format() {
        assert_eq!(e_fmt(1.5e-3, 20, 12, 3), " 0.150000000000E-002");
        assert_eq!(e_fmt(-2.0, 20, 12, 3), "-0.200000000000E+001");
        assert_eq!(e_fmt(0.0, 20, 12, 3), " 0.000000000000E+000");
    }
}
