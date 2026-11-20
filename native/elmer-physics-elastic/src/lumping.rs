//! Model lumping (`ModelLumping.F90`): the 6x6 spring matrix an elastic body
//! shows at the boundary flagged `Model Lumping Boundary`, from six load
//! cases (three translations, three rotations), plus the rigid-body mass.
//!
//! `Fix Displacement = True` (the default) prescribes each case as a
//! displacement of the lumping boundary and reads the reactions off the bulk
//! matrix; `False` applies a pure force or moment there and integrates the
//! response. The cases ride the nonlinear iteration counter (`iter` is the
//! case), which the solver forces to six.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use std::io::Write;

use elmer_core::{BasisEval, ElemId, ElementBasis, Params, Quadrature, Value};
use elmer_solver::fortran_fmt::es;
use elmer_solver::{SolverCtx, SolverError};

use crate::nodes::node_uvw;

const CALLER: &str = "ModelLumping";

fn fatal(msg: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, msg.into())
}

/// `ES17.8E3`.
fn es17_8e3(x: f64) -> String {
    let s = es(x, 16, 8);
    match s.rfind('E') {
        Some(i) if s.len() == i + 4 => format!("{}{}0{}", &s[..i + 1], &s[i + 1..i + 2], &s[i + 2..]),
        _ => format!("{s:>17}"),
    }
}

fn write_rows(path: &str, rows: &[[f64; 6]], ncols: usize) -> Result<(), SolverError> {
    let mut f = std::fs::File::create(path).map_err(|e| fatal(format!("{path}: {e}")))?;
    for r in rows {
        let line: String = r[..ncols].iter().map(|&v| es17_8e3(v)).collect();
        writeln!(f, "{line}").map_err(|e| fatal(format!("{path}: {e}")))?;
    }
    Ok(())
}

/// `LUDecomp` + `InvertMatrix` (LinearAlgebra.F90) on a 6x6, in their order.
fn invert_matrix(a: &mut [[f64; 6]; 6], n: usize) -> Result<(), SolverError> {
    let mut pivot = [0usize; 6];
    for i in 0..n {
        let mut j = i;
        for k in i + 1..n {
            if a[i][k].abs() > a[i][j].abs() {
                j = k;
            }
        }
        if a[i][j].abs() == 0.0 {
            return Err(SolverError::fatal("InvertMatrix", "inversion needs successful LU-decomposition"));
        }
        pivot[i] = j;
        if j != i {
            for k in 0..=i {
                a[k].swap(i, j);
            }
        }
        for k in i + 1..n {
            a[i][k] /= a[i][i];
        }
        for k in i + 1..n {
            if j != i {
                a[k].swap(i, j);
            }
            for l in i + 1..n {
                a[k][l] -= a[k][i] * a[i][l];
            }
        }
    }
    pivot[n - 1] = n - 1;
    for i in 0..n {
        if a[i][i].abs() == 0.0 {
            return Err(SolverError::fatal("InvertMatrix", "inversion needs successful LU-decomposition"));
        }
        a[i][i] = 1.0 / a[i][i];
    }
    for i in (0..n - 1).rev() {
        for j in (i + 1..n).rev() {
            let mut s = -a[i][j];
            for k in i + 1..j {
                s -= a[i][k] * a[k][j];
            }
            a[i][j] = s;
        }
    }
    for i in (0..n - 1).rev() {
        for j in (i + 1..n).rev() {
            let mut s = 0.0;
            for k in i + 1..=j {
                s -= a[j][k] * a[k][i];
            }
            a[j][i] = a[i][i] * s;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in i.max(j)..n {
                if k != i {
                    s += a[i][k] * a[k][j];
                } else {
                    s += a[k][j];
                }
            }
            a[i][j] = s;
        }
    }
    for i in (0..n).rev() {
        if pivot[i] != i {
            a.swap(i, pivot[i]);
        }
    }
    Ok(())
}

/// `ModelLumping_t`: kept across the six cases.
pub struct Lump {
    area: f64,
    center: [f64; 3],
    moments: [[f64; 3]; 3],
    kmat: [[f64; 6]; 6],
    kmat_min: [[f64; 6]; 6],
    pub fix_displacement: bool,
    file: String,
    symmetrisize: bool,
}

/// The boundary elements flagged `Model Lumping Boundary`.
fn lumping_boundary(ctx: &SolverCtx<'_>) -> Vec<ElemId> {
    ctx.boundary_elements()
        .filter(|&e| ctx.bc_params(e).is_some_and(|b| b.logical("model lumping boundary") == Some(true)))
        .collect()
}

impl Lump {
    /// `ModelLumpingInit`: the boundary geometry and the rigid-body mass.
    pub fn init(ctx: &mut SolverCtx<'_>, basis: &dyn ElementBasis, quad: &dyn Quadrature) -> Result<Self, SolverError> {
        let p = ctx.params();
        let fix_displacement = p.logical("fix displacement").unwrap_or(true);
        let file = p.string("model lumping filename").unwrap_or("Kmat.dat").trim_matches('"').to_string();
        let symmetrisize = p.logical("symmetrisize") == Some(true);
        ctx.info(
            CALLER,
            if fix_displacement {
                "Using six fixed displacement to compute the spring matrix"
            } else {
                "Using six pure forces and moments to compute the spring matrix"
            },
            5,
        );
        let mut lump = Lump {
            area: 0.0,
            center: [0.0; 3],
            moments: [[0.0; 3]; 3],
            kmat: [[0.0; 6]; 6],
            kmat_min: [[0.0; 6]; 6],
            fix_displacement,
            file,
            symmetrisize,
        };
        lump.coordinate_integrals(ctx, basis, quad)?;
        lump.cartesian_mass(ctx, basis, quad)?;
        Ok(lump)
    }

    /// `CoordinateIntegrals`: area, centre of area, second moments.
    fn coordinate_integrals(&mut self, ctx: &SolverCtx<'_>, basis: &dyn ElementBasis, quad: &dyn Quadrature) -> Result<(), SolverError> {
        let elems = lumping_boundary(ctx);
        if elems.is_empty() {
            return Err(fatal("Model lumping boundary must be defined"));
        }
        let mut be = BasisEval::default();
        for power in 1..=2 {
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                for qp in quad.gauss_points(el.etype, None) {
                    basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(fatal)?;
                    let s = be.det_j * qp.weight;
                    let mut c = [0.0; 3];
                    for (d, cd) in c.iter_mut().enumerate() {
                        let mut acc = 0.0;
                        for k in 0..n {
                            acc += be.basis[k] * xyz[k][d];
                        }
                        *cd = acc;
                    }
                    if power == 1 {
                        self.area += s;
                        for d in 0..3 {
                            self.center[d] += s * c[d];
                        }
                    } else {
                        for d in 0..3 {
                            c[d] -= self.center[d];
                        }
                        for i in 0..3 {
                            for j in 0..3 {
                                self.moments[i][j] += s * c[i] * c[j];
                            }
                        }
                    }
                }
            }
            if power == 1 {
                for d in 0..3 {
                    self.center[d] /= self.area;
                }
            }
        }
        Ok(())
    }

    /// `CartesianMass`: mass and moments of inertia about the lumping
    /// centre, as `res:` values in the Simulation section and `<file>.mass`.
    fn cartesian_mass(&self, ctx: &mut SolverCtx<'_>, basis: &dyn ElementBasis, quad: &dyn Quadrature) -> Result<(), SolverError> {
        let (mut vol, mut m0) = (0.0, 0.0);
        let mut m2 = [[0.0f64; 3]; 3];
        let c = self.center;
        let mut be = BasisEval::default();
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        for e in elems {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let mat = ctx
                .body_params(e)
                .and_then(|b| b.material.cloned())
                .ok_or_else(|| fatal(format!("element {} has no material", e.0)))?;
            let density = ctx.real_at_nodes(&mat, "density", e)?.unwrap_or_else(|| vec![0.0; n]);
            for qp in quad.gauss_points(el.etype, None) {
                basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(fatal)?;
                let s = be.det_j * qp.weight;
                let at = |d: usize| {
                    let mut acc = 0.0;
                    for k in 0..n {
                        acc += xyz[k][d] * be.basis[k];
                    }
                    acc
                };
                let (x, y, z) = (at(0) - c[0], at(1) - c[1], at(2) - c[2]);
                vol += s;
                let mut dens = 0.0;
                for k in 0..n {
                    dens += be.basis[k] * density[k];
                }
                m0 += s * dens;
                m2[0][0] += s * (y * y + z * z) * dens;
                m2[1][1] += s * (x * x + z * z) * dens;
                m2[2][2] += s * (x * x + y * y) * dens;
                m2[0][1] -= s * x * y * dens;
                m2[0][2] -= s * x * z * dens;
                m2[1][2] -= s * y * z * dens;
            }
        }
        if vol < 10.0 * f64::EPSILON {
            return Ok(());
        }
        m2[1][0] = m2[0][1];
        m2[2][0] = m2[0][2];
        m2[2][1] = m2[1][2];
        let sim = &mut ctx.sim.model.simulation;
        for (k, v) in [
            ("res: Mass", m0),
            ("res: Lumped Center X", c[0]),
            ("res: Lumped Center Y", c[1]),
            ("res: Lumped Center Z", c[2]),
            ("res: Moment of inertia XX", m2[0][0]),
            ("res: Moment of inertia YY", m2[1][1]),
            ("res: Moment of inertia ZZ", m2[2][2]),
            ("res: Moment of inertia XY", m2[0][1]),
            ("res: Moment of inertia XZ", m2[0][2]),
            ("res: Moment of inertia YZ", m2[1][2]),
        ] {
            sim.insert(k, Value::real(v));
        }
        let mut mass = [[0.0; 6]; 6];
        for i in 0..3 {
            mass[i][i] = m0;
            for j in 0..3 {
                mass[3 + i][3 + j] = m2[i][j];
            }
        }
        write_rows(&format!("{}.mass", self.file), &mass, 6)
    }

    /// `ModelLumpingLoads`: the boundary load of `case` at the nodes `xyz`
    /// of one lumping boundary element (replaces its surface load).
    pub fn loads(&self, case: usize, xyz: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = xyz.len();
        let mut f = vec![[0.0; 3]; n];
        if case <= 3 {
            for fk in f.iter_mut() {
                fk[case - 1] = 1.0 / self.area;
            }
        } else if case <= 6 {
            let ix = (case - 4) % 3;
            let iy = (case - 3) % 3;
            let iz = (case - 2) % 3;
            let (zc, yc) = match case {
                4 => (2, 1),
                5 => (0, 2),
                _ => (1, 0),
            };
            let m = &self.moments;
            let eps = 1.0e-6;
            let isy = m[iy][ix].abs() < eps * m[iy][iy];
            let isz = m[iz][ix].abs() < eps * m[iz][iz];
            for (k, fk) in f.iter_mut().enumerate() {
                let (y, z) = (xyz[k][yc], xyz[k][zc]);
                if isy {
                    let c = 1.0 / m[iy][iy];
                    fk[iz] = c * (y - self.center[iy]);
                } else if isz {
                    let c = -1.0 / m[iz][iz];
                    fk[iy] = c * (z - self.center[iz]);
                } else {
                    let c = 1.0 / (m[iy][iy] + m[iz][iz]);
                    fk[iy] = -c * (z - self.center[iz]);
                    fk[iz] = c * (y - self.center[iy]);
                }
            }
        }
        f
    }

    /// `ModelLumpingDisplacements`: the prescribed displacement of `case` on
    /// the lumping boundary, as `(dof, value)` for `DefaultDirichletBCs` to
    /// override where the sif's own conditions speak.
    pub fn displacements(&self, ctx: &SolverCtx<'_>, case: usize) -> Result<Vec<(u32, f64)>, SolverError> {
        let f = ctx.primary().ok_or_else(|| fatal("solver has no displacement variable"))?;
        let (perm, dofs) = (f.perm.ok_or_else(|| fatal("displacement has no permutation"))?, f.dofs);
        let dx = 1.0e-2 * self.area.sqrt();
        let dfii = 1.0e-2;
        let mut out = Vec::new();
        for e in lumping_boundary(ctx) {
            for &node in ctx.mesh.element(e).nodes {
                let Some(&k) = perm.get(node as usize) else { continue };
                if k == elmer_core::NO_DOF {
                    continue;
                }
                let mut d = [0.0; 3];
                if case <= 3 {
                    d[case - 1] = dx;
                } else {
                    let x = ctx.mesh.coord(node);
                    let c = [x[0] - self.center[0], x[1] - self.center[1], x[2] - self.center[2]];
                    match case {
                        4 => {
                            d[1] = -dfii * c[2];
                            d[2] = dfii * c[1];
                        }
                        5 => {
                            d[0] = dfii * c[2];
                            d[2] = -dfii * c[0];
                        }
                        _ => {
                            d[0] = -dfii * c[1];
                            d[1] = dfii * c[0];
                        }
                    }
                }
                for (l, &v) in d.iter().enumerate().take(dofs.min(3)) {
                    out.push(((dofs * k as usize + l) as u32, v));
                }
            }
        }
        Ok(out)
    }

    /// `ModelLumpingSprings`: one row or column of the spring matrix from
    /// `case`; after the sixth, write it out.
    pub fn springs(
        &mut self,
        ctx: &SolverCtx<'_>,
        basis: &dyn ElementBasis,
        quad: &dyn Quadrature,
        case: usize,
    ) -> Result<(), SolverError> {
        let f = ctx.primary().ok_or_else(|| fatal("solver has no displacement variable"))?;
        let (perm, disp, dofs) = (f.perm.ok_or_else(|| fatal("displacement has no permutation"))?, f.values, f.dofs);
        let dfii = 1.0e-2;
        let dx = 1.0e-2 * self.area.sqrt();
        let pc = case - 1;
        if case == 1 {
            self.kmat = [[0.0; 6]; 6];
            self.kmat_min = [[f64::MAX; 6]; 6];
        }
        let elems = lumping_boundary(ctx);
        if self.fix_displacement {
            // NodalLoads = BulkValues * Displacement.
            let sys = ctx.system().ok_or_else(|| fatal("no assembled system"))?;
            let (bulk, _) = ctx
                .saved_bulk()
                .ok_or_else(|| fatal("Model lumping needs the bulk matrix (Constant Bulk System)"))?;
            let mut a = sys.matrix.clone();
            a.values.copy_from_slice(bulk);
            let mut loads = vec![0.0; disp.len()];
            elmer_solver::crs_matrix_vector_multiply(&a, dofs, disp, &mut loads);
            let mut visited = vec![false; disp.len() / dofs];
            for e in &elems {
                for &node in ctx.mesh.element(*e).nodes {
                    let j = perm[node as usize] as usize;
                    if visited[j] {
                        continue;
                    }
                    visited[j] = true;
                    let x = ctx.mesh.coord(node);
                    let c = [x[0] - self.center[0], x[1] - self.center[1], x[2] - self.center[2]];
                    let mut fa = [0.0; 3];
                    for k in 0..3 {
                        fa[k] = loads[3 * j + k];
                    }
                    let mo = [
                        -fa[1] * c[2] + fa[2] * c[1],
                        -fa[2] * c[0] + fa[0] * c[2],
                        -fa[0] * c[1] + fa[1] * c[0],
                    ];
                    for k in 0..3 {
                        self.kmat[k][pc] += fa[k];
                        self.kmat[3 + k][pc] += mo[k];
                    }
                }
            }
        } else {
            let in_perm = |p: ElemId| {
                ctx.mesh.element(p).nodes.iter().all(|&k| perm.get(k as usize).is_some_and(|&s| s != elmer_core::NO_DOF))
            };
            let (mut be, mut pbe) = (BasisEval::default(), BasisEval::default());
            for &e in &elems {
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let parent = match (el.parents[0], el.parents[1]) {
                    (Some(l), _) if in_perm(l) => l,
                    (_, Some(r)) if in_perm(r) => r,
                    _ => return Err(fatal("Cannot find proper parent for side element")),
                };
                let pel = ctx.mesh.element(parent);
                let pn = pel.nodes.len();
                let pxyz: Vec<[f64; 3]> = pel.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                let ld: Vec<[f64; 3]> = pel
                    .nodes
                    .iter()
                    .map(|&k| {
                        let s = perm[k as usize] as usize;
                        std::array::from_fn(|i| if i < dofs { disp[dofs * s + i] } else { 0.0 })
                    })
                    .collect();
                let table = node_uvw(pel.etype).ok_or_else(|| fatal(format!("parent element type {}", pel.etype.0)))?;
                let mut pu = vec![[0.0; 3]; n];
                for (i, bn) in el.nodes.iter().enumerate() {
                    if let Some(j) = pel.nodes.iter().position(|pn| pn == bn) {
                        pu[i] = table[j];
                    }
                }
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                for qp in quad.gauss_points(el.etype, None) {
                    basis.eval(el.etype, &xyz, qp.uvw, &mut be).map_err(fatal)?;
                    let s = be.det_j * qp.weight;
                    let mut uvw = [0.0; 3];
                    for d in 0..3 {
                        let mut acc = 0.0;
                        for i in 0..n {
                            acc += pu[i][d] * be.basis[i];
                        }
                        uvw[d] = acc;
                    }
                    basis.eval(pel.etype, &pxyz, uvw, &mut pbe).map_err(fatal)?;
                    let sum = |g: &dyn Fn(usize) -> f64| {
                        let mut acc = 0.0;
                        for k in 0..pn {
                            acc += g(k);
                        }
                        acc
                    };
                    let mut fa = [0.0; 3];
                    for i in 0..3 {
                        fa[i] = sum(&|k| pbe.basis[k] * ld[k][i]);
                    }
                    let db = &pbe.dbasis_dx;
                    let mo = [
                        0.5 * (sum(&|k| db[k][1] * ld[k][2]) - sum(&|k| db[k][2] * ld[k][1])),
                        0.5 * (sum(&|k| db[k][2] * ld[k][0]) - sum(&|k| db[k][0] * ld[k][2])),
                        0.5 * (sum(&|k| db[k][0] * ld[k][1]) - sum(&|k| db[k][1] * ld[k][0])),
                    ];
                    for i in 0..3 {
                        self.kmat[pc][i] += s * fa[i];
                        self.kmat[pc][3 + i] += s * mo[i];
                    }
                    for i in 0..3 {
                        if self.kmat_min[pc][i].abs() > fa[i].abs() {
                            self.kmat_min[pc][i] = fa[i];
                        }
                        if self.kmat_min[pc][3 + i].abs() > mo[i].abs() {
                            self.kmat_min[pc][3 + i] = mo[i];
                        }
                    }
                }
            }
        }

        if case == 6 {
            ctx.info(CALLER, &format!("Saving lumped elastic spring to file {}", self.file), 4);
            let sym = |m: &mut [[f64; 6]; 6]| {
                let t = *m;
                for i in 0..6 {
                    for j in 0..6 {
                        m[i][j] = (t[i][j] + t[j][i]) / 2.0;
                    }
                }
            };
            if self.fix_displacement {
                for row in self.kmat.iter_mut() {
                    for j in 0..3 {
                        row[j] /= dx;
                    }
                    for j in 3..6 {
                        row[j] /= dfii;
                    }
                }
                if self.symmetrisize {
                    sym(&mut self.kmat);
                }
            } else {
                for row in self.kmat.iter_mut() {
                    for v in row.iter_mut() {
                        *v /= self.area;
                    }
                }
                write_rows(&format!("{}.inv", self.file), &self.kmat, 6)?;
                write_rows(&format!("{}.min-inv", self.file), &self.kmat_min, 6)?;
                if self.symmetrisize {
                    sym(&mut self.kmat);
                    sym(&mut self.kmat_min);
                }
                invert_matrix(&mut self.kmat, 6)?;
                invert_matrix(&mut self.kmat_min, 6)?;
                write_rows(&format!("{}.min", self.file), &self.kmat_min, 6)?;
            }
            write_rows(&self.file, &self.kmat, 6)?;
            let center = [[self.center[0], self.center[1], self.center[2], 0.0, 0.0, 0.0]];
            write_rows(&format!("{}.center", self.file), &center, 3)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invert_matrix_inverts() {
        let mut a = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                a[i][j] = if i == j { 4.0 + i as f64 } else { 1.0 / (1.0 + (i + 2 * j) as f64) };
            }
        }
        let a0 = a;
        invert_matrix(&mut a, 6).unwrap();
        for i in 0..6 {
            for j in 0..6 {
                let s: f64 = (0..6).map(|k| a0[i][k] * a[k][j]).sum();
                assert!((s - if i == j { 1.0 } else { 0.0 }).abs() < 1e-13, "{i},{j}: {s}");
            }
        }
    }

    #[test]
    fn es17_8e3_has_a_three_digit_exponent() {
        assert_eq!(es17_8e3(1.0), "  1.00000000E+000");
        assert_eq!(es17_8e3(-2.5e-7), " -2.50000000E-007");
    }
}
