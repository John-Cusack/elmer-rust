//! Element kernels of `WhitneyAVSolver` (steady state, lowest-order
//! Whitney basis without the Piola transform): `LocalMatrix`,
//! `LocalMatrixBC` and `LocalMatrixAirGapBC`.
//!
//! Local DOFs are ordered as Elmer's `GetElementDOFs`: the `np` nodal
//! (scalar potential) DOFs first, then one DOF per element edge in
//! [`crate::topology::edge_map`] order. Matrices are `nd x nd` row-major.
//! Nodal coefficients are interpolated with the basis as the Fortran does
//! (`SUM(Basis(1:n) * c(1:n))`), so a constant keyword still goes through
//! that sum.

use elmer_core::{BasisEval, ElementType};
use elmer_fe::metric::MetricStatus;
use elmer_fe::quadrature::{gauss_points, GaussRule};

use crate::error::MagnetoError;

/// The element integration rule of edge-basis solvers: `GaussPointsAdapt`
/// with `EdgeBasis` for the lowest-order non-Piola basis picks
/// `ElementalNp = [1,2,3,4,4,27,6,8]` by family (SolverBasics.F90:4780),
/// on the standard (non-p) reference element. Boundary elements use
/// `EdgeElementGaussPoints` (Integration.F90:3358), which gives the same
/// counts for families 1-5.
pub fn edge_rule(code: u16) -> Result<GaussRule, MagnetoError> {
    const NP: [usize; 8] = [1, 2, 3, 4, 4, 27, 6, 8];
    let family = (code / 100) as usize;
    if !(1..=8).contains(&family) {
        return Err(MagnetoError::Input(format!("element type {code}")));
    }
    gauss_points(ElementType(code), Some(NP[family - 1])).map_err(|e| MagnetoError::Input(e.to_string()))
}

/// `ElementInfo(..., EdgeBasis=WBasis, RotBasis=RotWBasis)` at one point.
pub struct PointEval {
    pub det_j: f64,
    pub basis: Vec<f64>,
    pub dbasis: Vec<[f64; 3]>,
    pub w: Vec<[f64; 3]>,
    pub curl: Vec<[f64; 3]>,
}

pub fn eval_point(
    code: u16,
    xyz: &[[f64; 3]],
    gindex: &[i64],
    uvw: [f64; 3],
    status: &mut MetricStatus,
) -> Result<PointEval, MagnetoError> {
    let def = elmer_fe::element_type(code).map_err(|e| MagnetoError::Input(e.to_string()))?;
    let mut be = BasisEval::default();
    elmer_fe::info::element_info(def, xyz, uvw, 3, status, &mut be).map_err(|e| MagnetoError::Input(e.to_string()))?;
    let ev = elmer_fe_hp::edge::get_edge_basis(ElementType(code), xyz, gindex, uvw, 3).map_err(MagnetoError::Input)?;
    Ok(PointEval {
        det_j: be.det_j,
        basis: be.basis,
        dbasis: be.dbasis_dx,
        w: ev.edge,
        curl: ev.curl,
    })
}

fn interp(basis: &[f64], c: &[f64]) -> f64 {
    let mut s = 0.0;
    for (b, c) in basis.iter().zip(c) {
        s += c * b;
    }
    s
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let mut s = 0.0;
    for i in 0..3 {
        s += a[i] * b[i];
    }
    s
}

fn matvec(c: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    let mut r = [0.0; 3];
    for (i, ri) in r.iter_mut().enumerate() {
        let mut s = 0.0;
        for j in 0..3 {
            s += c[i][j] * v[j];
        }
        *ri = s;
    }
    r
}

/// An element's local system.
pub struct Local {
    pub nd: usize,
    pub stiff: Vec<f64>,
    pub force: Vec<f64>,
}

impl Local {
    fn new(nd: usize) -> Self {
        Self {
            nd,
            stiff: vec![0.0; nd * nd],
            force: vec![0.0; nd],
        }
    }

    fn add(&mut self, p: usize, q: usize, v: f64) {
        self.stiff[p * self.nd + q] += v;
    }

    pub fn at(&self, p: usize, q: usize) -> f64 {
        self.stiff[p * self.nd + q]
    }
}

/// Nodal coefficients of a bulk element (`n` entries each).
pub struct BulkCoeffs {
    /// Reluctivity (`Acoef`).
    pub nu: Vec<f64>,
    /// Conductivity tensor (`Tcoef(:,:,k)`).
    pub sigma: Vec<[[f64; 3]; 3]>,
    /// `LOAD(1:7, k)`: current density, magnetization, electric potential.
    pub load: Vec<[f64; 7]>,
}

/// `LocalMatrix` (WhitneyAVSolver.F90:2068), steady branch, isotropic
/// reluctivity, no laminate stack, coil or velocity terms. `np` is the
/// number of nodal DOFs (`n` with the scalar potential, else 0).
pub fn bulk_matrix(
    code: u16,
    xyz: &[[f64; 3]],
    gindex: &[i64],
    np: usize,
    co: &BulkCoeffs,
) -> Result<Local, MagnetoError> {
    let n = (code % 100) as usize;
    let ne = crate::topology::edge_map(code / 100).len();
    let nd = np + ne;
    let mut loc = Local::new(nd);
    let rule = edge_rule(code)?;
    let mut status = MetricStatus::default();
    for t in 0..rule.len() {
        let pt = eval_point(code, xyz, gindex, [rule.u[t], rule.v[t], rule.w[t]], &mut status)?;
        let (basis, db, s) = (&pt.basis[..n], &pt.dbasis[..n], rule.s[t]);
        let det_j = pt.det_j;
        let mu = interp(basis, &co.nu);
        let mut c = [[0.0; 3]; 3];
        for (i, ci) in c.iter_mut().enumerate() {
            for (j, cij) in ci.iter_mut().enumerate() {
                let mut v = 0.0;
                for k in 0..n {
                    v += co.sigma[k][i][j] * basis[k];
                }
                *cij = v;
            }
        }
        let mut l = [0.0; 3];
        let mut m = [0.0; 3];
        for i in 0..3 {
            for k in 0..n {
                l[i] += co.load[k][i] * basis[k];
                m[i] += co.load[k][3 + i] * basis[k];
            }
        }
        // L = L - C grad(V^s).
        let mut gv = [0.0; 3];
        for (j, g) in gv.iter_mut().enumerate() {
            for k in 0..n {
                *g += co.load[k][6] * db[k][j];
            }
        }
        let cg = matvec(&c, &gv);
        for i in 0..3 {
            l[i] -= cg[i];
        }

        // Conductor: <C grad V, grad v> and <C grad V, eta>.
        let csum: f64 = c.iter().flatten().map(|x| x.abs()).sum();
        if csum > elmer_sparse::AEPS {
            for p in 0..np {
                for q in 0..np {
                    loc.add(p, q, dot(&matvec(&c, &db[q]), &db[p]) * det_j * s);
                }
                let cdp = matvec(&c, &db[p]);
                for (j, wj) in pt.w.iter().enumerate().take(ne) {
                    loc.add(np + j, p, dot(&cdp, wj) * det_j * s);
                }
            }
        }

        // curl(nu curl A) and the sources.
        for i in 0..ne {
            let p = np + i;
            loc.force[p] += (dot(&l, &pt.w[i]) + dot(&m, &pt.curl[i])) * det_j * s;
            for j in 0..ne {
                loc.add(p, np + j, mu * dot(&pt.curl[i], &pt.curl[j]) * det_j * s);
            }
        }
    }
    Ok(loc)
}

/// The real parts of the harmonic `LocalMatrix` (WhitneyAVHarmonicSolver.F90:1234)
/// for real coefficients and without electrodynamics: the stiffness `K`
/// and damping `D` matrices, combined by the caller into `K + i omega D`
/// (the `MASS` term is zero). Same layout and coefficients as
/// [`bulk_matrix`]; `D` carries `<C A, eta>` and `<C A, grad v>`.
pub fn harmonic_bulk_matrix(
    code: u16,
    xyz: &[[f64; 3]],
    gindex: &[i64],
    np: usize,
    co: &BulkCoeffs,
) -> Result<(Local, Vec<f64>), MagnetoError> {
    let n = (code % 100) as usize;
    let ne = crate::topology::edge_map(code / 100).len();
    let nd = np + ne;
    let mut loc = Local::new(nd);
    let mut damp = vec![0.0; nd * nd];
    let rule = edge_rule(code)?;
    let mut status = MetricStatus::default();
    for t in 0..rule.len() {
        let pt = eval_point(code, xyz, gindex, [rule.u[t], rule.v[t], rule.w[t]], &mut status)?;
        let (basis, db, s) = (&pt.basis[..n], &pt.dbasis[..n], rule.s[t]);
        let det_j = pt.det_j;
        let mut c = [[0.0; 3]; 3];
        for (i, ci) in c.iter_mut().enumerate() {
            for (j, cij) in ci.iter_mut().enumerate() {
                let mut v = 0.0;
                for k in 0..n {
                    v += co.sigma[k][i][j] * basis[k];
                }
                *cij = v;
            }
        }
        let mu = interp(basis, &co.nu);
        let mut l = [0.0; 3];
        let mut m = [0.0; 3];
        for i in 0..3 {
            for k in 0..n {
                l[i] += co.load[k][i] * basis[k];
                m[i] += co.load[k][3 + i] * basis[k];
            }
        }
        let mut gv = [0.0; 3];
        for (j, g) in gv.iter_mut().enumerate() {
            for k in 0..n {
                *g += co.load[k][6] * db[k][j];
            }
        }
        let cg = matvec(&c, &gv);
        for i in 0..3 {
            l[i] -= cg[i];
        }

        let csum: f64 = c.iter().flatten().map(|x| x.abs()).sum();
        if csum > elmer_sparse::AEPS {
            for p in 0..np {
                for q in 0..np {
                    loc.add(p, q, dot(&matvec(&c, &db[q]), &db[p]) * det_j * s);
                }
                let cdp = matvec(&c, &db[p]);
                for (j, wj) in pt.w.iter().enumerate().take(ne) {
                    damp[p * nd + np + j] += dot(&matvec(&c, wj), &db[p]) * det_j * s;
                    loc.add(np + j, p, dot(&cdp, wj) * det_j * s);
                }
            }
        }

        for i in 0..ne {
            let p = np + i;
            loc.force[p] += (dot(&l, &pt.w[i]) + dot(&m, &pt.curl[i])) * det_j * s;
            for j in 0..ne {
                let q = np + j;
                let cj = &pt.curl[j];
                loc.add(p, q, dot(&[mu * cj[0], mu * cj[1], mu * cj[2]], &pt.curl[i]) * det_j * s);
                damp[p * nd + q] += dot(&matvec(&c, &pt.w[j]), &pt.w[i]) * det_j * s;
            }
        }
    }
    Ok((loc, damp))
}

/// Nodal coefficients of `LocalMatrixBC`.
pub struct BcCoeffs {
    /// `Magnetic Transfer Coefficient` (`Bcoef`).
    pub transfer: Vec<f64>,
    /// `LOAD(1:3)` magnetic field strength, `LOAD(4)` electric current
    /// density, `LOAD(5)` electric transfer coefficient.
    pub load: Vec<[f64; 5]>,
}

/// `LocalMatrixBC` (WhitneyAVSolver.F90:2691) without the perfect electric
/// conductor and outflow terms.
pub fn bc_matrix(code: u16, xyz: &[[f64; 3]], gindex: &[i64], np: usize, co: &BcCoeffs) -> Result<Local, MagnetoError> {
    let n = (code % 100) as usize;
    let ne = crate::topology::edge_map(code / 100).len();
    let mut loc = Local::new(np + ne);
    let rule = edge_rule(code)?;
    let mut status = MetricStatus::default();
    for t in 0..rule.len() {
        let pt = eval_point(code, xyz, gindex, [rule.u[t], rule.v[t], rule.w[t]], &mut status)?;
        let (basis, s, det_j) = (&pt.basis[..n], rule.s[t], pt.det_j);
        let b = interp(basis, &co.transfer);
        let mut l = [0.0; 3];
        for (i, li) in l.iter_mut().enumerate() {
            for k in 0..n {
                *li += co.load[k][i] * basis[k];
            }
        }
        let f = interp(basis, &co.load.iter().map(|x| x[3]).collect::<Vec<_>>());
        let tc = interp(basis, &co.load.iter().map(|x| x[4]).collect::<Vec<_>>());
        for p in 0..np {
            loc.force[p] += f * basis[p] * det_j * s;
            for q in 0..np {
                loc.add(p, q, tc * basis[p] * basis[q] * det_j * s);
            }
        }
        for i in 0..ne {
            let p = np + i;
            loc.force[p] -= dot(&l, &pt.w[i]) * det_j * s;
            for j in 0..ne {
                loc.add(p, np + j, b * dot(&pt.w[i], &pt.w[j]) * det_j * s);
            }
        }
    }
    Ok(loc)
}

/// `LocalMatrixAirGapBC` (WhitneyAVSolver.F90:2828):
/// `d / (mu_r mu_0) <curl A, curl eta>` over the gap surface.
pub fn air_gap_matrix(
    code: u16,
    xyz: &[[f64; 3]],
    gindex: &[i64],
    np: usize,
    gap_length: &[f64],
    gap_mu: &[f64],
) -> Result<Local, MagnetoError> {
    let n = (code % 100) as usize;
    let ne = crate::topology::edge_map(code / 100).len();
    let mut loc = Local::new(np + ne);
    let mu_vacuum = 4.0 * crate::material::PI * 1e-7;
    let rule = edge_rule(code)?;
    let mut status = MetricStatus::default();
    for t in 0..rule.len() {
        let pt = eval_point(code, xyz, gindex, [rule.u[t], rule.v[t], rule.w[t]], &mut status)?;
        let (basis, s, det_j) = (&pt.basis[..n], rule.s[t], pt.det_j);
        let len = interp(basis, gap_length);
        let mu_air = interp(basis, gap_mu);
        for i in 0..ne {
            for j in 0..ne {
                loc.add(
                    np + i,
                    np + j,
                    len / (mu_air * mu_vacuum) * dot(&pt.curl[i], &pt.curl[j]) * det_j * s,
                );
            }
        }
    }
    Ok(loc)
}

/// `LocalMatrixThinSheet` (WhitneyAVHarmonicSolver.F90:1828): a thin
/// conducting, permeable sheet of thickness `d` on a boundary face. Returns
/// the stiffness and the damping matrix (`<d C A, grad v>`, `<d C A, eta>`)
/// the caller combines into `K + i omega D`.
pub fn thin_sheet_matrix(
    code: u16,
    xyz: &[[f64; 3]],
    gindex: &[i64],
    np: usize,
    thickness: &[f64],
    permeability: &[f64],
    conductivity: &[f64],
) -> Result<(Local, Vec<f64>), MagnetoError> {
    let n = (code % 100) as usize;
    let ne = crate::topology::edge_map(code / 100).len();
    let nd = np + ne;
    let mut loc = Local::new(nd);
    let mut damp = vec![0.0; nd * nd];
    let mu_vacuum = 4.0 * crate::material::PI * 1e-7;
    let rule = edge_rule(code)?;
    let mut status = MetricStatus::default();
    for t in 0..rule.len() {
        let pt = eval_point(code, xyz, gindex, [rule.u[t], rule.v[t], rule.w[t]], &mut status)?;
        let (basis, db, s, det_j) = (&pt.basis[..n], &pt.dbasis[..n], rule.s[t], pt.det_j);
        let d = interp(basis, thickness);
        let mu = interp(basis, permeability);
        let c = interp(basis, conductivity);
        let conductor = c.abs() > elmer_sparse::AEPS;
        if conductor {
            for p in 0..np {
                for q in 0..np {
                    loc.add(p, q, d * c * dot(&db[q], &db[p]) * det_j * s);
                }
                for (j, wj) in pt.w.iter().enumerate().take(ne) {
                    damp[p * nd + np + j] += d * c * dot(wj, &db[p]) * det_j * s;
                    loc.add(np + j, p, d * c * dot(&db[p], wj) * det_j * s);
                }
            }
        }
        for i in 0..ne {
            let p = np + i;
            for j in 0..ne {
                let q = np + j;
                loc.add(p, q, d / (mu * mu_vacuum) * dot(&pt.curl[i], &pt.curl[j]) * det_j * s);
                if conductor {
                    damp[p * nd + q] += d * c * dot(&pt.w[j], &pt.w[i]) * det_j * s;
                }
            }
        }
    }
    Ok((loc, damp))
}

#[cfg(test)]
mod tests {
    use super::*;

    const TET: [[f64; 3]; 4] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    fn coeffs(nu: f64, sigma: f64) -> BulkCoeffs {
        let mut c = [[0.0; 3]; 3];
        for (i, ci) in c.iter_mut().enumerate() {
            ci[i] = sigma;
        }
        BulkCoeffs {
            nu: vec![nu; 4],
            sigma: vec![c; 4],
            load: vec![[0.0; 7]; 4],
        }
    }

    #[test]
    fn curl_curl_block_is_symmetric_and_kills_gradients() {
        let loc = bulk_matrix(504, &TET, &[1, 2, 3, 4], 0, &coeffs(2.0, 0.0)).unwrap();
        assert_eq!(loc.nd, 6);
        for p in 0..6 {
            for q in 0..6 {
                assert!((loc.at(p, q) - loc.at(q, p)).abs() < 1e-14);
            }
        }
        // A nodal gradient has edge DOFs u_b - u_a along each edge (a -> b
        // towards the larger global number): its curl, and so K x, is zero.
        let phi = [0.3, -1.2, 0.7, 2.0];
        let x: Vec<f64> = crate::topology::edge_map(5)
            .iter()
            .map(|&[a, b]| phi[a.max(b)] - phi[a.min(b)])
            .collect();
        for p in 0..6 {
            let r: f64 = (0..6).map(|q| loc.at(p, q) * x[q]).sum();
            assert!(r.abs() < 1e-13, "row {p}: {r}");
        }
    }

    #[test]
    fn conductor_couples_potential_into_the_edge_rows_only() {
        let loc = bulk_matrix(504, &TET, &[1, 2, 3, 4], 4, &coeffs(1.0, 3.0)).unwrap();
        assert_eq!(loc.nd, 10);
        // <C grad V, grad v> is the P1 Laplacian times sigma: rows sum to 0.
        for p in 0..4 {
            let r: f64 = (0..4).map(|q| loc.at(p, q)).sum();
            assert!(r.abs() < 1e-13);
        }
        // Steady state: the V rows carry no A coupling.
        for p in 0..4 {
            for q in 4..10 {
                assert_eq!(loc.at(p, q), 0.0);
            }
        }
        assert!((4..10).any(|p| (0..4).any(|q| loc.at(p, q) != 0.0)));
    }

    #[test]
    fn current_density_bc_loads_the_potential_rows() {
        let tri = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let co = BcCoeffs {
            transfer: vec![0.0; 3],
            load: vec![[0.0, 0.0, 0.0, 5.0, 0.0]; 3],
        };
        let loc = bc_matrix(303, &tri, &[1, 2, 3], 3, &co).unwrap();
        // J * area = 5 * 2, split evenly over the three nodes.
        let total: f64 = loc.force[..3].iter().sum();
        assert!((total - 10.0).abs() < 1e-13);
        assert!(loc.force[3..].iter().all(|&f| f == 0.0));
    }
}
