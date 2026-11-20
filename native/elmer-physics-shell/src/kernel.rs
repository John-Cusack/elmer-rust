//! `ShellLocalMatrix` (ShellSolver.F90:3454) for a four-node
//! quadrilateral: the Reissner-Mindlin shell of the lines-of-curvature
//! patch with the membrane strains reduced by the curl kernel with edge
//! DOFs (operator 4) and the transverse shear by MITC (operator 2), the
//! large-deflection (von Karman type) terms, and the normal pressure load.
//!
//! The local system is in the patch's surface basis and is rotated to the
//! Cartesian components by `Q` at the end, as in the Fortran.

use elmer_fe::basis::nodal_basis;
use elmer_fe::element_type;

use crate::geom::{
    cross, dot, geometry_basis, quad_p_basis, surface_basis, Basis3, Frame, MAX_PATCH_NODES,
};
use crate::ShellError;

/// The unknowns per node.
pub(crate) const M: usize = 6;
const ND: usize = 4;
const DOFS: usize = M * ND;

/// Reduction operators (`MITC`, `CurlKernelWithEdgeDOFs`).
#[derive(Clone, Copy, PartialEq)]
enum Reduction {
    Mitc,
    CurlKernelEdge,
}

/// Options of the element matrix that come from the solver section.
#[derive(Clone, Copy, Debug)]
pub(crate) struct KernelOpts {
    pub large_deflection: bool,
    pub nonlinear_bending: bool,
    pub shear_alpha: f64,
    pub membrane_alpha: f64,
    pub stretch_alpha: f64,
    pub benchmark: bool,
}

/// Nodal material and load data of the element.
pub(crate) struct KernelInput<'a> {
    pub thickness: &'a [f64],
    pub poisson: &'a [f64],
    pub youngs: &'a [f64],
    pub pressure: &'a [f64],
    /// `LocalSol(k, j)`: the previous solution, node-major.
    pub local_sol: &'a [f64],
    /// Gauss points of the 404 element: `(u, v, s)`.
    pub points: &'a [(f64, f64, f64)],
}

/// The element's local system: `Stiff`, `Force` and the load part
/// `RHSForce` (large deflection only), row-major `DOFS x DOFS`.
pub(crate) struct KernelOutput {
    pub stiff: Vec<f64>,
    pub force: Vec<f64>,
    pub rhs_force: Vec<f64>,
}

type Mat = Vec<Vec<f64>>;

fn zeros(r: usize, c: usize) -> Mat {
    vec![vec![0.0; c]; r]
}

/// `MATMUL( A, x )`.
fn matvec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter().map(|row| dot(row, x)).collect()
}

/// `MATMUL( TRANSPOSE( A ), x )`.
fn matvec_t(a: &[Vec<f64>], x: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .map(|j| {
            let mut s = 0.0;
            for (row, xi) in a.iter().zip(x) {
                s += row[j] * xi;
            }
            s
        })
        .collect()
}

/// `StrainEnergyDensity( A, B, C, m, n, s )`: `A += s C^T (B C)`.
fn strain_energy(a: &mut [f64], b: &[Vec<f64>], c: &[Vec<f64>], m: usize, n: usize, s: f64) {
    let mut t = zeros(m, n);
    for i in 0..m {
        for j in 0..n {
            let mut acc = 0.0;
            for k in 0..m {
                acc += b[i][k] * c[k][j];
            }
            t[i][j] = acc;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0;
            for k in 0..m {
                acc += c[k][i] * t[k][j];
            }
            a[i * n + j] += s * acc;
        }
    }
}

/// `MATMUL( TRANSPOSE( x(1:1,:) ), y(1:1,:) )`: the outer product.
fn outer(x: &[f64], y: &[f64]) -> Vec<f64> {
    let mut r = vec![0.0; DOFS * DOFS];
    for i in 0..DOFS {
        for j in 0..DOFS {
            r[i * DOFS + j] = 0.0 + x[i] * y[j];
        }
    }
    r
}

fn add_rows(a: &[Vec<f64>], b: &[Vec<f64>]) -> Mat {
    a.iter()
        .zip(b)
        .map(|(x, y)| x.iter().zip(y).map(|(p, q)| p + q).collect())
        .collect()
}

/// `ShearCorrectionFactor( Kappa, h, x, y, 4, alpha )`.
fn shear_correction(h: f64, x: &[f64; 4], y: &[f64; 4], stab: f64) -> f64 {
    let alpha = 0.10 * stab;
    let l = |i: usize, j: usize| ((x[i] - x[j]).powi(2) + (y[i] - y[j]).powi(2)).sqrt();
    let hh = l(1, 0).max(l(2, 1)).max(l(3, 2)).max(l(0, 3));
    h.powi(2) / (h.powi(2) + alpha * hh.powi(2))
}

/// `ElasticityMatrix( CMat, GMat, A1, A2, E, nu )` without drilling DOFs.
fn elasticity(a1: f64, a2: f64, e: f64, nu: f64) -> ([[f64; 4]; 4], [[f64; 2]; 2]) {
    let mut c = [[0.0f64; 4]; 4];
    c[0][0] = 1.0;
    c[0][1] = nu;
    c[1][0] = nu;
    c[1][1] = 1.0;
    c[2][2] = (1.0 - nu) / 2.0;
    for row in &mut c {
        for v in row.iter_mut() {
            *v = *v * e / (1.0 - nu.powi(2));
        }
    }
    c[0][0] /= a1.powi(4);
    c[0][1] /= a1.powi(2) * a2.powi(2);
    c[1][0] /= a1.powi(2) * a2.powi(2);
    c[1][1] /= a2.powi(4);
    c[2][2] /= a1.powi(2) * a2.powi(2);
    c[3][3] = (1.0 - nu) * e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mut g = [[0.0f64; 2]; 2];
    g[0][0] = e / (2.0 * (1.0 + nu) * a1.powi(2));
    g[1][1] = e / (2.0 * (1.0 + nu) * a2.powi(2));
    (c, g)
}

/// `ReductionOperatorInfo( Element, Nodes, u, v, StrainBasis, method,
/// ApplyPiolaTransform=.TRUE., ... )`: the reduced strain basis, the
/// Piola-mapped basis gradients, the bilinear basis and `|det F|`.
fn reduction_info(
    nodes: &[[f64; 2]; 4],
    u: f64,
    v: f64,
    method: Reduction,
    edge_dir: usize,
) -> ([[f64; 2]; 4], [[f64; 2]; 4], [f64; 4], f64) {
    let (lbasis, dlbasis) = quad_p_basis(u, v);
    // PiolaTransformationData: F(1:2, i) and det F (z = 0).
    let mut lf = [[0.0f64; 2]; 2];
    for i in 0..2 {
        let (mut sx, mut sy) = (0.0, 0.0);
        for k in 0..4 {
            sx += nodes[k][0] * dlbasis[k][i];
        }
        for k in 0..4 {
            sy += nodes[k][1] * dlbasis[k][i];
        }
        lf[0][i] = sx;
        lf[1][i] = sy;
    }
    let det = lf[0][0] * lf[1][1] - lf[0][1] * lf[1][0];
    let g = [
        [1.0 / det * lf[1][1], -1.0 / det * lf[0][1]],
        [-1.0 / det * lf[1][0], 1.0 / det * lf[0][0]],
    ];
    // LG = TRANSPOSE( LG )
    let lg = [[g[0][0], g[1][0]], [g[0][1], g[1][1]]];
    let mut sb = [[0.0f64; 2]; 4];
    let dofs = match method {
        Reduction::Mitc => {
            sb[0] = [(1.0 - v) / 4.0, 0.0];
            sb[1] = [0.0, (1.0 + u) / 4.0];
            sb[2] = [-(1.0 + v) / 4.0, 0.0];
            sb[3] = [0.0, (-1.0 + u) / 4.0];
            4
        }
        Reduction::CurlKernelEdge => {
            if edge_dir == 1 {
                sb[0] = [(1.0 - v) / 4.0, -1.0 / 4.0 * u];
                sb[1] = [-(1.0 + v) / 4.0, -1.0 / 4.0 * u];
                sb[2] = [0.0, 0.5];
            } else {
                sb[0] = [1.0 / 4.0 * v, -(1.0 - u) / 4.0];
                sb[1] = [1.0 / 4.0 * v, (1.0 + u) / 4.0];
                sb[2] = [0.5, 0.0];
            }
            3
        }
    };
    let mut db = [[0.0f64; 2]; 4];
    for i in 0..4 {
        for j in 0..2 {
            let mut s = 0.0;
            for k in 0..2 {
                s += dlbasis[i][k] * lg[j][k];
            }
            db[i][j] = s;
        }
    }
    for j in 0..dofs {
        let mut b = [0.0f64; 2];
        for k in 0..2 {
            b[k] = dot(&lg[k], &sb[j]);
        }
        sb[j] = b;
    }
    (sb, db, lbasis, det.abs())
}

/// `ReductionOperatorDofs( Element, Nodes, A, nd, 4, method, ModelPars,
/// EdgeDirection )` of a quadrilateral: the DOFs of the reduced strain in
/// terms of the nodal values, `A(i, 2(j-1)+k)`.
fn reduction_dofs(
    nodes: &[[f64; 2]; 4],
    method: Reduction,
    pars: Option<&[[[f64; 2]; 2]; 4]>,
    edge_dir: usize,
) -> Result<[[f64; 8]; 4], ShellError> {
    let def = element_type(404).map_err(|e| ShellError::Fatal(e.to_string()))?;
    let par_at = |u: f64, v: f64| -> [[f64; 2]; 2] {
        match pars {
            None => [[1.0, 0.0], [0.0, 1.0]],
            Some(mp) => {
                let mut b = [0.0f64; 4];
                nodal_basis(def, u, v, 0.0, &mut b);
                let mut p = [[0.0f64; 2]; 2];
                for i in 0..2 {
                    for j in 0..2 {
                        let mut s = 0.0;
                        for k in 0..4 {
                            s += mp[k][i][j] * b[k];
                        }
                        p[i][j] = s;
                    }
                }
                p
            }
        }
    };
    let entry = |p: &[[f64; 2]; 2], tau: &[f64; 3], col: usize, half: f64| {
        half * (p[0][col] * tau[0] + p[1][col] * tau[1]) * tau[2]
    };
    let unit = |t: [f64; 2]| {
        let l = (t[0].powi(2) + t[1].powi(2)).sqrt();
        [t[0] / l, t[1] / l, l]
    };
    let x = |i: usize| nodes[i][0];
    let y = |i: usize| nodes[i][1];
    let mut a = [[0.0f64; 8]; 4];
    match method {
        Reduction::Mitc => {
            let mut tau = [[0.0f64; 3]; 4];
            for i in 1..4 {
                tau[i - 1] = unit([x(i) - x(i - 1), y(i) - y(i - 1)]);
            }
            tau[3] = unit([x(0) - x(3), y(0) - y(3)]);
            let pts = [(0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)];
            let cols = [[0usize, 1, 2, 3], [2, 3, 4, 5], [4, 5, 6, 7], [6, 7, 0, 1]];
            for r in 0..4 {
                let p = par_at(pts[r].0, pts[r].1);
                let c = cols[r];
                a[r][c[0]] = entry(&p, &tau[r], 0, 0.5);
                a[r][c[1]] = entry(&p, &tau[r], 1, 0.5);
                a[r][c[2]] = entry(&p, &tau[r], 0, 0.5);
                a[r][c[3]] = entry(&p, &tau[r], 1, 0.5);
            }
        }
        Reduction::CurlKernelEdge => {
            let raw = if edge_dir == 1 {
                [
                    [x(1) - x(0), y(1) - y(0)],
                    [x(3) - x(2), y(3) - y(2)],
                    [
                        (x(3) + x(2)) / 2.0 - (x(1) + x(0)) / 2.0,
                        (y(3) + y(2)) / 2.0 - (y(1) + y(0)) / 2.0,
                    ],
                ]
            } else {
                [
                    [x(0) - x(3), y(0) - y(3)],
                    [x(2) - x(1), y(2) - y(1)],
                    [
                        (x(2) + x(1)) / 2.0 - (x(3) + x(0)) / 2.0,
                        (y(2) + y(1)) / 2.0 - (y(3) + y(0)) / 2.0,
                    ],
                ]
            };
            let tau = [unit(raw[0]), unit(raw[1]), unit(raw[2])];
            if edge_dir == 1 {
                let p = par_at(0.0, -1.0);
                for (c, k) in [(0, 0), (1, 1), (2, 0), (3, 1)] {
                    a[0][c] = entry(&p, &tau[0], k, 0.5);
                }
                let p = par_at(0.0, 1.0);
                for (c, k) in [(4, 0), (5, 1), (6, 0), (7, 1)] {
                    a[1][c] = entry(&p, &tau[1], k, 0.5);
                }
            } else {
                let p = par_at(-1.0, 0.0);
                for (c, k) in [(6, 0), (7, 1), (0, 0), (1, 1)] {
                    a[0][c] = entry(&p, &tau[0], k, 0.5);
                }
                let p = par_at(1.0, 0.0);
                for (c, k) in [(2, 0), (3, 1), (4, 0), (5, 1)] {
                    a[1][c] = entry(&p, &tau[1], k, 0.5);
                }
            }
            let p = par_at(0.0, 0.0);
            for c in 0..8 {
                a[2][c] = entry(&p, &tau[2], c % 2, 0.25);
            }
        }
    }
    Ok(a)
}

/// `WriteElementNodesVariables`: the element's nodes in the patch
/// coordinates of the 416 geometry element.
fn local_nodes(frame: &Frame) -> Result<[[f64; 2]; 4], ShellError> {
    let def = element_type(404).map_err(|e| ShellError::Fatal(e.to_string()))?;
    let mut out = [[0.0f64; 2]; 4];
    for j in 0..4 {
        let g = geometry_basis(def.node_u[j], def.node_v[j])?;
        let (mut sx, mut sy) = (0.0, 0.0);
        for k in 0..MAX_PATCH_NODES {
            sx += g[k] * frame.patch[k][0];
        }
        for k in 0..MAX_PATCH_NODES {
            sy += g[k] * frame.patch[k][1];
        }
        out[j] = [sx, sy];
    }
    Ok(out)
}

/// `SUM( a(1:n) * Basis(1:n) )`.
fn interp(a: &[f64], b: &[f64; 4]) -> f64 {
    let mut s = 0.0;
    for k in 0..4 {
        s += a[k] * b[k];
    }
    s
}

/// `ShellLocalMatrix( BGElement, 4, 4, 6, LocalSol, ... )`.
pub(crate) fn shell_local_matrix(
    frame: &Frame,
    inp: &KernelInput<'_>,
    o: &KernelOpts,
) -> Result<KernelOutput, ShellError> {
    let nodes = local_nodes(frame)?;
    let plate = frame.planar;
    let membrane = Reduction::CurlKernelEdge;
    let shear = Reduction::Mitc;
    let strain_dim = 4;
    let membrane_dim = 3;

    // The Christoffel symbols and curvatures at the nodes.
    let mut chris1 = [[[0.0f64; 2]; 2]; 4];
    let mut chris2 = [[[0.0f64; 2]; 2]; 4];
    let mut bpar = [[[0.0f64; 2]; 2]; 4];
    let mut bpar1 = [[[0.0f64; 2]; 2]; 4];
    let mut bpar2 = [[[0.0f64; 2]; 2]; 4];
    let mut shear_par = [[[0.0f64; 2]; 2]; 4];
    for j in 0..ND {
        let s = surface_basis(nodes[j][0], nodes[j][1], frame);
        chris1[j] = [[s.c111, s.c112], [s.c211, s.c212]];
        chris2[j] = [[s.c211, s.c212], [s.c221, s.c222]];
        bpar[j][0][0] = s.b11;
        bpar[j][1][1] = s.b22;
        bpar1[j][0][0] = s.b11;
        bpar2[j][1][1] = s.b22;
        shear_par[j][0][0] = s.b11 / s.a11;
        shear_par[j][1][1] = s.b22 / s.a22;
    }
    let mut prev = vec![0.0f64; DOFS];
    prev.copy_from_slice(&inp.local_sol[..DOFS]);
    // Q: the surface basis at the nodes, rows a1, a2, a3.
    let mut q = vec![0.0f64; DOFS * DOFS];
    for j in 0..ND {
        let s = surface_basis(nodes[j][0], nodes[j][1], frame);
        let block = [s.a1, s.a2, s.a3];
        let i0 = j * M;
        for r in 0..3 {
            for c in 0..3 {
                q[(i0 + r) * DOFS + i0 + c] = block[r][c];
                q[(i0 + 3 + r) * DOFS + i0 + 3 + c] = block[r][c];
            }
        }
    }
    {
        let mut out = vec![0.0f64; DOFS];
        for i in 0..DOFS {
            let mut s = 0.0;
            for j in 0..DOFS {
                s += q[i * DOFS + j] * prev[j];
            }
            out[i] = s;
        }
        prev = out;
    }

    let mut stiff = vec![0.0f64; DOFS * DOFS];
    let mut force = vec![0.0f64; DOFS];
    let mut rhs_force = vec![0.0f64; DOFS];
    let xs = [nodes[0][0], nodes[1][0], nodes[2][0], nodes[3][0]];
    let ys = [nodes[0][1], nodes[1][1], nodes[2][1], nodes[3][1]];
    let ld = o.large_deflection;
    let nb = o.nonlinear_bending;

    for &(uq, vq, sq) in inp.points {
        let mut bm = zeros(4, DOFS);
        let mut bb = zeros(3, DOFS);
        let mut bs = zeros(4, DOFS);
        let mut nl_bm = zeros(4, DOFS);
        let mut nl_bs = zeros(4, DOFS);
        let mut nl_bb = zeros(3, DOFS);
        let mut bw = zeros(13, DOFS);
        let mut pf = [0.0f64; 13];
        let (strain_basis, dbasis, basis, det_j) = reduction_info(&nodes, uq, vq, shear, 1);
        let (msb1, _, _, _) = reduction_info(&nodes, uq, vq, membrane, 1);
        let (msb2, _, _, _) = reduction_info(&nodes, uq, vq, membrane, 2);
        let y1 = interp(&xs, &basis);
        let y2 = interp(&ys, &basis);
        let s: Basis3 = surface_basis(y1, y2, frame);
        let a1n = s.a11.sqrt();
        let a2n = s.a22.sqrt();
        let h = interp(inp.thickness, &basis);
        let nu = interp(inp.poisson, &basis);
        let e = interp(inp.youngs, &basis);
        let traction = interp(inp.pressure, &basis);
        let (mut cmat, gmat) = elasticity(a1n, a2n, e, nu);
        let kappa_s = shear_correction(h, &xs, &ys, o.shear_alpha);
        let kappa_m = shear_correction(h, &xs, &ys, o.membrane_alpha);
        let kappa_str = shear_correction(h, &xs, &ys, o.stretch_alpha);
        cmat[3][3] = kappa_str / kappa_m * cmat[3][3];
        let mut weight = h * s.sqrt_det_a * det_j * sq;
        let (c111, c112, c221, c222, c211, c212) = (s.c111, s.c112, s.c221, s.c222, s.c211, s.c212);
        let (b11, b22, a11, a22) = (s.b11, s.b22, s.a11, s.a22);

        // Membrane strains.
        if !plate {
            let r1 = reduction_dofs(&nodes, membrane, Some(&chris1), 1)?;
            for p in 0..ND {
                for j in 0..membrane_dim {
                    bm[0][p * M] -= r1[j][2 * p] * msb1[j][0];
                    bm[2][p * M] -= r1[j][2 * p] * msb1[j][1];
                    bm[0][p * M + 1] -= r1[j][2 * p + 1] * msb1[j][0];
                    bm[2][p * M + 1] -= r1[j][2 * p + 1] * msb1[j][1];
                }
            }
            let r2 = reduction_dofs(&nodes, membrane, Some(&chris2), 2)?;
            for p in 0..ND {
                for j in 0..membrane_dim {
                    bm[2][p * M] -= r2[j][2 * p] * msb2[j][0];
                    bm[1][p * M] -= r2[j][2 * p] * msb2[j][1];
                    bm[2][p * M + 1] -= r2[j][2 * p + 1] * msb2[j][0];
                    bm[1][p * M + 1] -= r2[j][2 * p + 1] * msb2[j][1];
                }
            }
            let r3 = reduction_dofs(&nodes, membrane, Some(&bpar1), 1)?;
            for p in 0..ND {
                for j in 0..membrane_dim {
                    bm[0][p * M + 2] -= r3[j][2 * p] * msb1[j][0];
                    bm[2][p * M + 2] -= r3[j][2 * p] * msb1[j][1];
                }
            }
            let r4 = reduction_dofs(&nodes, membrane, Some(&bpar2), 2)?;
            for p in 0..ND {
                for j in 0..membrane_dim {
                    bm[2][p * M + 2] -= r4[j][2 * p + 1] * msb2[j][0];
                    bm[1][p * M + 2] -= r4[j][2 * p + 1] * msb2[j][1];
                }
            }
            for p in 0..ND {
                bm[0][p * M] += dbasis[p][0];
                bm[1][p * M + 1] += dbasis[p][1];
                bm[2][p * M] += dbasis[p][1];
                bm[2][p * M + 1] += dbasis[p][0];
            }
        } else {
            for p in 0..ND {
                bm[0][p * M] = dbasis[p][0] - c111 * basis[p];
                bm[0][p * M + 1] = -c112 * basis[p];
                bm[0][p * M + 2] = -b11 * basis[p];
                bm[1][p * M] = -c221 * basis[p];
                bm[1][p * M + 1] = dbasis[p][1] - c222 * basis[p];
                bm[1][p * M + 2] = -b22 * basis[p];
                bm[2][p * M] = dbasis[p][1] - 2.0 * c211 * basis[p];
                bm[2][p * M + 1] = dbasis[p][0] - 2.0 * c212 * basis[p];
            }
        }
        for p in 0..ND {
            bm[3][p * M + 5] = -basis[p];
        }
        {
            let f1 = nu / ((1.0 - nu) * a1n.powi(2));
            let f2 = nu / ((1.0 - nu) * a2n.powi(2));
            for c in 0..DOFS {
                bm[3][c] = bm[3][c] + f1 * bm[0][c] + f2 * bm[1][c];
            }
        }
        let mut strain = [0.0f64; 6];
        let mut stress = [0.0f64; 6];
        let bm_prev = |row: &[f64]| dot(row, &prev);
        if ld {
            for p in 0..ND {
                bw[0][p * M] = -c211 * basis[p];
                bw[0][p * M + 1] = dbasis[p][0] - c212 * basis[p];
            }
            pf[0] = dot(&bw[0], &prev);
            for p in 0..ND {
                bw[1][p * M] = b11 / a11 * basis[p];
                bw[1][p * M + 2] = dbasis[p][0];
            }
            pf[1] = dot(&bw[1], &prev);
            let s1 = bm_prev(&bm[0]);
            for c in 0..DOFS {
                nl_bm[0][c] =
                    s1 / a1n.powi(2) * bm[0][c] + pf[0] / a2n.powi(2) * bw[0][c] + pf[1] * bw[1][c];
            }
            strain[0] = 0.5 * dot(&nl_bm[0], &prev);
            for p in 0..ND {
                bw[2][p * M] = dbasis[p][1] - c211 * basis[p];
                bw[2][p * M + 1] = -c212 * basis[p];
            }
            pf[2] = dot(&bw[2], &prev);
            for p in 0..ND {
                bw[3][p * M + 1] = b22 / a22 * basis[p];
                bw[3][p * M + 2] = dbasis[p][1];
            }
            pf[3] = dot(&bw[3], &prev);
            let s2 = bm_prev(&bm[1]);
            for c in 0..DOFS {
                nl_bm[1][c] =
                    s2 / a2n.powi(2) * bm[1][c] + pf[2] / a1n.powi(2) * bw[2][c] + pf[3] * bw[3][c];
            }
            strain[1] = 0.5 * dot(&nl_bm[1], &prev);
            for c in 0..DOFS {
                nl_bm[2][c] = bm[0][c] * pf[2] / a1n.powi(2)
                    + s1 / a1n.powi(2) * bw[2][c]
                    + bm[1][c] * pf[0] / a2n.powi(2)
                    + s2 / a2n.powi(2) * bw[0][c]
                    + bw[1][c] * pf[3]
                    + pf[1] * bw[3][c];
            }
            strain[2] = 0.5 * dot(&nl_bm[2], &prev);
            for p in 0..ND {
                bw[4][p * M + 3] = basis[p];
                bw[5][p * M + 4] = basis[p];
                bw[6][p * M + 5] = basis[p];
            }
            pf[4] = dot(&bw[4], &prev);
            pf[5] = dot(&bw[5], &prev);
            pf[6] = dot(&bw[6], &prev);
            for c in 0..DOFS {
                nl_bm[3][c] = nu / (1.0 - nu) * nl_bm[0][c] / a1n.powi(2)
                    + nu / (1.0 - nu) * nl_bm[1][c] / a2n.powi(2)
                    + bw[4][c] * pf[4] / a1n.powi(2)
                    + bw[5][c] * pf[5] / a2n.powi(2)
                    + bw[6][c] * pf[6];
            }
            strain[3] = 0.5 * dot(&nl_bm[3], &prev);
        }
        let cm: Vec<Vec<f64>> = cmat.iter().map(|r| r.to_vec()).collect();
        let bm_t = add_rows(&bm, &nl_bm);
        strain_energy(&mut stiff, &cm, &bm_t, 4, DOFS, kappa_m * weight);
        let lin = matvec(&bm, &prev);
        for i in 0..4 {
            strain[i] += lin[i];
        }
        let st = matvec(&cm, &strain[..4]);
        stress[..4].copy_from_slice(&st);
        let fm = matvec_t(&bm_t, &stress[..4], DOFS);
        for c in 0..DOFS {
            force[c] -= fm[c] * weight * kappa_m;
        }
        if ld {
            let f1 = stress[0] + nu / (1.0 - nu) * stress[3] / a1n.powi(2);
            let f2 = stress[1] + nu / (1.0 - nu) * stress[3] / a2n.powi(2);
            let (o_bm0, o_bw0, o_bw1) = (
                outer(&bm[0], &bm[0]),
                outer(&bw[0], &bw[0]),
                outer(&bw[1], &bw[1]),
            );
            for k in 0..DOFS * DOFS {
                stiff[k] += o_bm0[k] / a1n.powi(2) * f1 * kappa_m * weight
                    + o_bw0[k] / a2n.powi(2) * f1 * kappa_m * weight
                    + o_bw1[k] * f1 * kappa_m * weight;
            }
            let (o_bm1, o_bw2, o_bw3) = (
                outer(&bm[1], &bm[1]),
                outer(&bw[2], &bw[2]),
                outer(&bw[3], &bw[3]),
            );
            for k in 0..DOFS * DOFS {
                stiff[k] += o_bm1[k] / a2n.powi(2) * f2 * kappa_m * weight
                    + o_bw2[k] / a1n.powi(2) * f2 * kappa_m * weight
                    + o_bw3[k] * f2 * kappa_m * weight;
            }
            let (x1, x2) = (outer(&bw[2], &bm[0]), outer(&bm[0], &bw[2]));
            let (x3, x4) = (outer(&bw[0], &bm[1]), outer(&bm[1], &bw[0]));
            let (x5, x6) = (outer(&bw[1], &bw[3]), outer(&bw[3], &bw[1]));
            for k in 0..DOFS * DOFS {
                stiff[k] += (x1[k] + x2[k]) * stress[2] / a1n.powi(2) * kappa_m * weight
                    + (x3[k] + x4[k]) * stress[2] / a2n.powi(2) * kappa_m * weight
                    + (x5[k] + x6[k]) * stress[2] * kappa_m * weight;
            }
            let (y4, y5, y6) = (
                outer(&bw[4], &bw[4]),
                outer(&bw[5], &bw[5]),
                outer(&bw[6], &bw[6]),
            );
            for k in 0..DOFS * DOFS {
                stiff[k] += y4[k] * stress[3] / a1n.powi(2) * kappa_m * weight
                    + y5[k] * stress[3] / a2n.powi(2) * kappa_m * weight
                    + y6[k] * stress[3] * kappa_m * weight;
            }
        }

        // Transverse shear strains (MITC).
        let rs = reduction_dofs(&nodes, shear, None, 1)?;
        for p in 0..ND {
            for j in 0..strain_dim {
                for r in 0..2 {
                    bs[r][p * M + 3] -= rs[j][2 * p] * strain_basis[j][r];
                    bs[r][p * M + 4] -= rs[j][2 * p + 1] * strain_basis[j][r];
                }
            }
        }
        let rp = reduction_dofs(&nodes, shear, Some(&shear_par), 1)?;
        for p in 0..ND {
            for j in 0..strain_dim {
                for r in 0..2 {
                    bs[r][p * M] += rp[j][2 * p] * strain_basis[j][r];
                    bs[r][p * M + 1] += rp[j][2 * p + 1] * strain_basis[j][r];
                }
            }
            bs[0][p * M + 2] = dbasis[p][0];
            bs[1][p * M + 2] = dbasis[p][1];
        }
        if ld {
            let s1 = bm_prev(&bm[0]);
            let s2 = bm_prev(&bm[1]);
            for c in 0..DOFS {
                nl_bs[0][c] = -bm[0][c] * pf[4] / a1n.powi(2)
                    - s1 / a1n.powi(2) * bw[4][c]
                    - bw[0][c] * pf[5] / a2n.powi(2)
                    - pf[0] / a2n.powi(2) * bw[5][c]
                    - bw[1][c] * pf[6]
                    - pf[1] * bw[6][c];
            }
            strain[4] = 0.5 * dot(&nl_bs[0], &prev);
            for c in 0..DOFS {
                nl_bs[1][c] = -bw[2][c] * pf[4] / a1n.powi(2)
                    - pf[2] / a1n.powi(2) * bw[4][c]
                    - bm[1][c] * pf[5] / a2n.powi(2)
                    - s2 / a2n.powi(2) * bw[5][c]
                    - bw[3][c] * pf[6]
                    - pf[3] * bw[6][c];
            }
            strain[5] = 0.5 * dot(&nl_bs[1], &prev);
        }
        let gm: Vec<Vec<f64>> = gmat.iter().map(|r| r.to_vec()).collect();
        let bs_t = add_rows(&bs[..2], &nl_bs[..2]);
        strain_energy(&mut stiff, &gm, &bs_t, 2, DOFS, kappa_s * weight);
        let lin = matvec(&bs[..2], &prev);
        strain[4] += lin[0];
        strain[5] += lin[1];
        let st = matvec(&gm, &strain[4..6]);
        stress[4] = st[0];
        stress[5] = st[1];
        let fs = matvec_t(&bs_t, &stress[4..6], DOFS);
        for c in 0..DOFS {
            force[c] -= fs[c] * kappa_s * weight;
        }
        if ld {
            let (x1, x2) = (outer(&bw[4], &bm[0]), outer(&bm[0], &bw[4]));
            let (x3, x4) = (outer(&bw[0], &bw[5]), outer(&bw[5], &bw[0]));
            let (x5, x6) = (outer(&bw[1], &bw[6]), outer(&bw[6], &bw[1]));
            for k in 0..DOFS * DOFS {
                stiff[k] = stiff[k]
                    - (x1[k] + x2[k]) * stress[4] / a1n.powi(2) * kappa_s * weight
                    - (x3[k] + x4[k]) * stress[4] / a2n.powi(2) * kappa_s * weight
                    - (x5[k] + x6[k]) * stress[4] * kappa_s * weight;
            }
            let (x1, x2) = (outer(&bw[4], &bw[2]), outer(&bw[2], &bw[4]));
            let (x3, x4) = (outer(&bm[1], &bw[5]), outer(&bw[5], &bm[1]));
            let (x5, x6) = (outer(&bw[3], &bw[6]), outer(&bw[6], &bw[3]));
            for k in 0..DOFS * DOFS {
                stiff[k] = stiff[k]
                    - (x1[k] + x2[k]) * stress[5] / a1n.powi(2) * kappa_s * weight
                    - (x3[k] + x4[k]) * stress[5] / a2n.powi(2) * kappa_s * weight
                    - (x5[k] + x6[k]) * stress[5] * kappa_s * weight;
            }
        }

        // Bending strains.
        weight = h.powi(3) / 12.0 * s.sqrt_det_a * det_j * sq;
        for p in 0..ND {
            let (bp, dx, dy) = (basis[p], dbasis[p][0], dbasis[p][1]);
            bb[0][p * M + 3] = dx - c111 * bp;
            bb[0][p * M + 4] = -c112 * bp;
            bb[1][p * M + 3] = -c221 * bp;
            bb[1][p * M + 4] = dy - c222 * bp;
            bb[2][p * M + 3] = dy - 2.0 * c211 * bp;
            bb[2][p * M + 4] = dx - 2.0 * c212 * bp;
            if o.benchmark {
                bb[2][p * M] = b11 / a11 * (dy - c211 * bp) - b22 / a22 * c211 * bp;
                bb[2][p * M + 1] = b11 / a11 * (-c212) * bp + b22 / a22 * (dx - c212 * bp);
            } else {
                bb[2][p * M] = b11 / a11 * c211 * bp - b22 / a22 * (dy - c211 * bp);
                bb[2][p * M + 1] = -b11 / a11 * (dx - c212 * bp) + b22 / a22 * c212 * bp;
            }
            bb[0][p * M + 5] = -b11 * bp;
            bb[1][p * M + 5] = -b22 * bp;
        }
        if ld && nb {
            bw[7] = bb[0].clone();
            bw[11] = bb[1].clone();
            pf[7] = dot(&bw[7], &prev);
            pf[11] = dot(&bw[11], &prev);
        }
        for c in 0..DOFS {
            bb[0][c] -= b11 / a11 * bm[0][c];
            bb[1][c] -= b22 / a22 * bm[1][c];
        }
        strain[0] = 0.0;
        strain[1] = 0.0;
        strain[2] = 0.0;
        if ld && nb {
            let s1 = bm_prev(&bm[0]);
            let s2 = bm_prev(&bm[1]);
            for p in 0..ND {
                bw[8][p * M + 3] = -c211 * basis[p];
                bw[8][p * M + 4] = dbasis[p][0] - c212 * basis[p];
            }
            pf[8] = dot(&bw[8], &prev);
            for p in 0..ND {
                bw[9][p * M + 3] = b11 / a11 * basis[p];
                bw[9][p * M + 5] = dbasis[p][0];
            }
            pf[9] = dot(&bw[9], &prev);
            for c in 0..DOFS {
                nl_bb[0][c] = bm[0][c] * pf[7] / a1n.powi(2)
                    + s1 / a1n.powi(2) * bw[7][c]
                    + bw[0][c] * pf[8] / a2n.powi(2)
                    + pf[0] / a2n.powi(2) * bw[8][c]
                    + bw[1][c] * pf[9]
                    + pf[1] * bw[9][c];
            }
            for p in 0..ND {
                bw[10][p * M + 3] = dbasis[p][1] - c211 * basis[p];
                bw[10][p * M + 4] = -c212 * basis[p];
            }
            pf[10] = dot(&bw[10], &prev);
            for p in 0..ND {
                bw[12][p * M + 4] = b22 / a22 * basis[p];
                bw[12][p * M + 5] = dbasis[p][1];
            }
            pf[12] = dot(&bw[12], &prev);
            for c in 0..DOFS {
                nl_bb[1][c] = bw[2][c] * pf[10] / a1n.powi(2)
                    + pf[2] / a1n.powi(2) * bw[10][c]
                    + bm[1][c] * pf[11] / a2n.powi(2)
                    + s2 / a2n.powi(2) * bw[11][c]
                    + bw[3][c] * pf[12]
                    + pf[3] * bw[12][c];
            }
            for c in 0..DOFS {
                nl_bb[2][c] = bm[0][c] * pf[10] / a1n.powi(2)
                    + s1 / a1n.powi(2) * bw[10][c]
                    + bw[7][c] * pf[2] / a1n.powi(2)
                    + pf[7] / a1n.powi(2) * bw[2][c]
                    + bw[0][c] * pf[11] / a2n.powi(2)
                    + pf[0] / a2n.powi(2) * bw[11][c]
                    + bm[1][c] * pf[8] / a2n.powi(2)
                    + s2 / a2n.powi(2) * bw[8][c]
                    + bw[1][c] * pf[12]
                    + pf[1] * bw[12][c]
                    + bw[9][c] * pf[3]
                    + pf[9] * bw[3][c];
            }
            for r in 0..3 {
                for c in 0..DOFS {
                    nl_bb[r][c] -= (b11 / a11 + b22 / a22) * nl_bm[r][c];
                }
            }
            let nlv = matvec(&nl_bb, &prev);
            for i in 0..3 {
                strain[i] = 0.5 * nlv[i];
            }
        }
        let lin = matvec(&bb, &prev);
        for i in 0..3 {
            strain[i] += lin[i];
        }
        let bb_t = add_rows(&bb, &nl_bb);
        strain_energy(&mut stiff, &cm, &bb_t, 3, DOFS, weight);
        let c3: Vec<Vec<f64>> = cmat[..3].iter().map(|r| r[..3].to_vec()).collect();
        let st = matvec(&c3, &strain[..3]);
        stress[..3].copy_from_slice(&st);
        let fb = matvec_t(&bb_t, &stress[..3], DOFS);
        for c in 0..DOFS {
            force[c] -= fb[c] * weight;
        }
        if ld && nb {
            let pairs = |a: &[f64], b: &[f64]| {
                let (x, y) = (outer(a, b), outer(b, a));
                x.iter().zip(&y).map(|(p, q)| p + q).collect::<Vec<f64>>()
            };
            let (p1, p2, p3) = (
                pairs(&bw[7], &bm[0]),
                pairs(&bw[8], &bw[0]),
                pairs(&bw[9], &bw[1]),
            );
            for k in 0..DOFS * DOFS {
                stiff[k] += p1[k] * stress[0] / a1n.powi(2) * weight
                    + p2[k] * stress[0] / a2n.powi(2) * weight
                    + p3[k] * stress[0] * weight;
            }
            let (p1, p2, p3) = (
                pairs(&bw[10], &bw[2]),
                pairs(&bw[11], &bm[1]),
                pairs(&bw[12], &bw[3]),
            );
            for k in 0..DOFS * DOFS {
                stiff[k] += p1[k] * stress[1] / a1n.powi(2) * weight
                    + p2[k] * stress[1] / a2n.powi(2) * weight
                    + p3[k] * stress[1] * weight;
            }
            // (BW11'BM1 + BM1'BW11 + BW3'BW8 + BW8'BW3), likewise the others.
            let q1 = {
                let (a, b, c, d) = (
                    outer(&bw[10], &bm[0]),
                    outer(&bm[0], &bw[10]),
                    outer(&bw[2], &bw[7]),
                    outer(&bw[7], &bw[2]),
                );
                (0..DOFS * DOFS)
                    .map(|k| a[k] + b[k] + c[k] + d[k])
                    .collect::<Vec<f64>>()
            };
            let q2 = {
                let (a, b, c, d) = (
                    outer(&bw[11], &bw[0]),
                    outer(&bw[0], &bw[11]),
                    outer(&bm[1], &bw[8]),
                    outer(&bw[8], &bm[1]),
                );
                (0..DOFS * DOFS)
                    .map(|k| a[k] + b[k] + c[k] + d[k])
                    .collect::<Vec<f64>>()
            };
            let q3 = {
                let (a, b, c, d) = (
                    outer(&bw[12], &bw[1]),
                    outer(&bw[1], &bw[12]),
                    outer(&bw[3], &bw[9]),
                    outer(&bw[9], &bw[3]),
                );
                (0..DOFS * DOFS)
                    .map(|k| a[k] + b[k] + c[k] + d[k])
                    .collect::<Vec<f64>>()
            };
            for k in 0..DOFS * DOFS {
                stiff[k] += q1[k] * stress[2] / a1n.powi(2) * weight
                    + q2[k] * stress[2] / a2n.powi(2) * weight
                    + q3[k] * stress[2] * weight;
            }
            let curv = b11 / a11 + b22 / a22;
            let (r1, r2, r3) = (
                outer(&bm[0], &bm[0]),
                outer(&bw[0], &bw[0]),
                outer(&bw[1], &bw[1]),
            );
            for k in 0..DOFS * DOFS {
                stiff[k] = stiff[k]
                    - (r1[k] / a1n.powi(2) + r2[k] / a2n.powi(2) + r3[k])
                        * curv
                        * stress[0]
                        * weight;
            }
            let (r1, r2, r3) = (
                outer(&bm[1], &bm[1]),
                outer(&bw[2], &bw[2]),
                outer(&bw[3], &bw[3]),
            );
            for k in 0..DOFS * DOFS {
                stiff[k] = stiff[k]
                    - (r1[k] / a2n.powi(2) + r2[k] / a1n.powi(2) + r3[k])
                        * curv
                        * stress[1]
                        * weight;
            }
            let (t1, t2) = (pairs(&bw[2], &bm[0]), pairs(&bw[0], &bm[1]));
            let (x, y) = (outer(&bw[1], &bw[3]), outer(&bw[3], &bw[1]));
            let cs = curv * stress[2] * weight;
            for k in 0..DOFS * DOFS {
                stiff[k] = stiff[k]
                    - cs * (t1[k] * 1.0 / a1n.powi(2) + t2[k] * 1.0 / a2n.powi(2) + x[k] + y[k]);
            }
        }

        // The drilling rotation: its gradient against the shear modulus.
        if !o.benchmark {
            for p in 0..ND {
                bs[2][p * M + 5] = dbasis[p][0];
                bs[3][p * M + 5] = dbasis[p][1];
            }
            strain[0] = 0.0;
            strain[1] = 0.0;
            if ld && nb {
                for c in 0..DOFS {
                    nl_bs[2][c] = -bw[4][c] * pf[7] / a1n.powi(2)
                        - bw[7][c] * pf[4] / a1n.powi(2)
                        - bw[5][c] * pf[8] / a2n.powi(2)
                        - bw[8][c] * pf[5] / a2n.powi(2)
                        - bw[6][c] * pf[9]
                        - bw[9][c] * pf[6];
                    nl_bs[3][c] = -bw[4][c] * pf[10] / a1n.powi(2)
                        - bw[10][c] * pf[4] / a1n.powi(2)
                        - bw[5][c] * pf[11] / a2n.powi(2)
                        - bw[11][c] * pf[5] / a2n.powi(2)
                        - bw[6][c] * pf[12]
                        - bw[12][c] * pf[6];
                }
                for c in 0..DOFS {
                    nl_bs[2][c] -= b11 / a11 * nl_bs[0][c];
                    nl_bs[3][c] -= b22 / a22 * nl_bs[1][c];
                }
                let nlv = matvec(&nl_bs[2..4], &prev);
                strain[0] = 0.5 * nlv[0];
                strain[1] = 0.5 * nlv[1];
            }
            let lin = matvec(&bs[2..4], &prev);
            strain[0] += lin[0];
            strain[1] += lin[1];
            let bd = add_rows(&bs[2..4], &nl_bs[2..4]);
            // Stiff + Weight * MATMUL( TRANSPOSE(BD), MATMUL( GMat, BD ) )
            let mut t = zeros(2, DOFS);
            for i in 0..2 {
                for j in 0..DOFS {
                    let mut acc = 0.0;
                    for k in 0..2 {
                        acc += gm[i][k] * bd[k][j];
                    }
                    t[i][j] = acc;
                }
            }
            for i in 0..DOFS {
                for j in 0..DOFS {
                    let mut acc = 0.0;
                    for k in 0..2 {
                        acc += bd[k][i] * t[k][j];
                    }
                    stiff[i * DOFS + j] += weight * acc;
                }
            }
            let st = matvec(&gm, &strain[..2]);
            stress[0] = st[0];
            stress[1] = st[1];
            let fd = matvec_t(&bd, &stress[..2], DOFS);
            for c in 0..DOFS {
                force[c] -= fd[c] * weight;
            }
            if ld && nb {
                let pairs = |a: &[f64], b: &[f64]| {
                    let (x, y) = (outer(a, b), outer(b, a));
                    x.iter().zip(&y).map(|(p, q)| p + q).collect::<Vec<f64>>()
                };
                // Note the Fortran's signs: (P1)/A1^2 - (P2)/A2^2 - X + Y.
                let p1 = pairs(&bw[7], &bw[4]);
                let p2 = pairs(&bw[8], &bw[5]);
                let (x, y) = (outer(&bw[9], &bw[6]), outer(&bw[6], &bw[9]));
                let sw = stress[0] * weight;
                for k in 0..DOFS * DOFS {
                    stiff[k] = stiff[k]
                        - sw * (p1[k] * 1.0 / a1n.powi(2) - p2[k] * 1.0 / a2n.powi(2) - x[k]
                            + y[k]);
                }
                let p1 = pairs(&bw[10], &bw[4]);
                let p2 = pairs(&bw[11], &bw[5]);
                let (x, y) = (outer(&bw[12], &bw[6]), outer(&bw[6], &bw[12]));
                let sw = stress[1] * weight;
                for k in 0..DOFS * DOFS {
                    stiff[k] = stiff[k]
                        - sw * (p1[k] * 1.0 / a1n.powi(2) - p2[k] * 1.0 / a2n.powi(2) - x[k]
                            + y[k]);
                }
                let p1 = pairs(&bw[4], &bm[0]);
                let p2 = pairs(&bw[0], &bw[5]);
                let (x, y) = (outer(&bw[1], &bw[6]), outer(&bw[6], &bw[1]));
                let sw = b11 / a11 * stress[0] * weight;
                for k in 0..DOFS * DOFS {
                    stiff[k] +=
                        sw * (p1[k] * 1.0 / a1n.powi(2) + p2[k] * 1.0 / a2n.powi(2) + x[k] + y[k]);
                }
                let p1 = pairs(&bw[4], &bw[2]);
                let p2 = pairs(&bm[1], &bw[5]);
                let (x, y) = (outer(&bw[3], &bw[6]), outer(&bw[6], &bw[3]));
                let sw = b22 / a22 * stress[1] * weight;
                for k in 0..DOFS * DOFS {
                    stiff[k] +=
                        sw * (p1[k] * 1.0 / a1n.powi(2) + p2[k] * 1.0 / a2n.powi(2) + x[k] + y[k]);
                }
            }
        }

        // The normal pressure.
        if ld {
            let comp = |k: usize| -> [f64; ND] { std::array::from_fn(|p| prev[p * M + k]) };
            let (u1, u2, u3) = (comp(0), comp(1), comp(2));
            let v1 = dot(&basis, &u1);
            let v2 = dot(&basis, &u2);
            let v3 = dot(&basis, &u3);
            let d1: [f64; ND] = std::array::from_fn(|p| dbasis[p][0]);
            let d2: [f64; ND] = std::array::from_fn(|p| dbasis[p][1]);
            let k11 = dot(&d1, &u1) - c111 * v1 - c112 * v2 - b11 * v3;
            let k12 = dot(&d1, &u2) - c211 * v1 - c212 * v2;
            let k13 = dot(&d1, &u3) + b11 / a11 * v1;
            let k21 = dot(&d2, &u1) - c211 * v1 - c212 * v2;
            let k22 = dot(&d2, &u2) - c221 * v1 - c222 * v2 - b22 * v3;
            let k23 = dot(&d2, &u3) + b22 / a22 * v2;
            let mut n1 = [0.0f64; 3];
            let mut n2 = [0.0f64; 3];
            for k in 0..3 {
                n1[k] = s.a1[k] + k11 * s.dual1[k] + k12 * s.dual2[k] + k13 * s.a3[k];
                n2[k] = s.a2[k] + k21 * s.dual1[k] + k22 * s.dual2[k] + k23 * s.a3[k];
            }
            let new_det = dot(&n1, &n1) * dot(&n2, &n2) - dot(&n1, &n2).powi(2);
            let n3 = cross(n1, n2);
            let nr = crate::geom::norm(&n3);
            let n3 = [n3[0] / nr, n3[1] / nr, n3[2] / nr];
            let w = new_det.sqrt() * det_j * sq;
            let t1 = traction * dot(&n3, &s.dual1);
            let t2 = traction * dot(&n3, &s.dual2);
            let t3 = traction * dot(&n3, &s.a3);
            for p in 0..ND {
                rhs_force[p * M] += t1 * basis[p] * w;
                rhs_force[p * M + 1] += t2 * basis[p] * w;
                rhs_force[p * M + 2] += t3 * basis[p] * w;
            }
            for p in 0..ND {
                force[p * M] += t1 * basis[p] * w;
                force[p * M + 1] += t2 * basis[p] * w;
                force[p * M + 2] += t3 * basis[p] * w;
            }
        } else {
            let w = s.sqrt_det_a * det_j * sq;
            for p in 0..ND {
                force[p * M + 2] += traction * basis[p] * w;
            }
        }
    }

    // Stiff = Q^T (Stiff Q), Force = Q^T Force.
    let mut sq_ = vec![0.0f64; DOFS * DOFS];
    for i in 0..DOFS {
        for j in 0..DOFS {
            let mut acc = 0.0;
            for k in 0..DOFS {
                acc += stiff[i * DOFS + k] * q[k * DOFS + j];
            }
            sq_[i * DOFS + j] = acc;
        }
    }
    let mut out = vec![0.0f64; DOFS * DOFS];
    for i in 0..DOFS {
        for j in 0..DOFS {
            let mut acc = 0.0;
            for k in 0..DOFS {
                acc += q[k * DOFS + i] * sq_[k * DOFS + j];
            }
            out[i * DOFS + j] = acc;
        }
    }
    let rot = |f: &[f64]| -> Vec<f64> {
        (0..DOFS)
            .map(|i| {
                let mut acc = 0.0;
                for k in 0..DOFS {
                    acc += q[k * DOFS + i] * f[k];
                }
                acc
            })
            .collect()
    };
    let force = rot(&force);
    let rhs_force = if ld { rot(&rhs_force) } else { rhs_force };
    Ok(KernelOutput {
        stiff: out,
        force,
        rhs_force,
    })
}
