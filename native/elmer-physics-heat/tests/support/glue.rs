//! Serial assembly glue and a CG solve over the kernels' local systems.

use std::collections::BTreeMap;
use std::path::Path;

use elmer_core::{ElementBasis, LocalSystem, Quadrature};
use elmer_physics_heat::{
    heat_boundary_local, heat_bulk_local, HeatBcAtNodes, HeatBoundaryConfig, HeatBulkConfig,
    HeatMaterialAtNodes,
};

use super::meshio::{Elem, Mesh};

/// Global system as sorted triplets. Rows and columns are dof numbers, i.e.
/// already through `perm`.
#[derive(Default, Clone, Debug)]
pub struct System {
    pub a: BTreeMap<(usize, usize), f64>,
    pub b: Vec<f64>,
}

impl System {
    pub fn new(n: usize) -> Self {
        Self {
            a: BTreeMap::new(),
            b: vec![0.0; n],
        }
    }

    pub fn add(&mut self, dofs: &[usize], local: &LocalSystem) {
        let n = local.n;
        for p in 0..n {
            for q in 0..n {
                *self.a.entry((dofs[p], dofs[q])).or_insert(0.0) += local.stiff[p * n + q];
            }
            self.b[dofs[p]] += local.force[p];
        }
    }

    /// Elmer's `SetDirichletBoundaries` on a non-symmetric CRS matrix: zero the
    /// row except the diagonal, set `b = diag * value`; other rows keep their
    /// entries in this column.
    pub fn dirichlet_row_scaled(&mut self, row: usize, value: f64) {
        let keys: Vec<_> = self
            .a
            .range((row, 0)..=(row, usize::MAX))
            .map(|(k, _)| *k)
            .collect();
        let diag = self.a.get(&(row, row)).copied().unwrap_or(0.0);
        for k in keys {
            if k.1 != row {
                self.a.insert(k, 0.0);
            }
        }
        self.b[row] = diag * value;
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.a.get(&(i, j)).copied().unwrap_or(0.0)
    }
}

/// Assembles bulk then boundary elements. `perm[node]` is the dof of a node.
pub fn assemble(
    mesh: &Mesh,
    perm: &[usize],
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    bulk_cfg: &HeatBulkConfig,
    bulk_mat: &dyn Fn(&Elem) -> HeatMaterialAtNodes,
    bnd_cfg: &HeatBoundaryConfig,
    bc: &dyn Fn(&Elem) -> Option<HeatBcAtNodes>,
) -> System {
    let mut sys = System::new(perm.iter().max().map_or(0, |m| m + 1));
    let mut local = LocalSystem::zeros(0);
    for e in &mesh.bulk {
        let xyz = mesh.xyz(e);
        heat_bulk_local(
            basis,
            quad,
            e.etype,
            &xyz,
            &bulk_mat(e),
            bulk_cfg,
            &mut local,
        )
        .unwrap();
        let dofs: Vec<usize> = e.nodes.iter().map(|&i| perm[i]).collect();
        sys.add(&dofs, &local);
    }
    for e in &mesh.boundary {
        let Some(bc) = bc(e) else { continue };
        let xyz = mesh.xyz(e);
        heat_boundary_local(basis, quad, e.etype, &xyz, &bc, bnd_cfg, &mut local).unwrap();
        let dofs: Vec<usize> = e.nodes.iter().map(|&i| perm[i]).collect();
        sys.add(&dofs, &local);
    }
    sys
}

/// Reads a `Linear System Save` matrix dump (`row col value`, 1-based).
pub fn read_matrix(path: &Path) -> BTreeMap<(usize, usize), f64> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter_map(|l| {
            let t: Vec<&str> = l.split_whitespace().collect();
            (t.len() >= 3).then(|| {
                (
                    (
                        t[0].parse::<usize>().unwrap() - 1,
                        t[1].parse::<usize>().unwrap() - 1,
                    ),
                    t[2].parse::<f64>().unwrap(),
                )
            })
        })
        .collect()
}

/// Reads an `i value` dump (rhs, solution) into a 0-based vector.
pub fn read_vector(path: &Path) -> Vec<f64> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter_map(|l| {
            let t: Vec<&str> = l.split_whitespace().collect();
            (t.len() >= 2).then(|| t[1].parse::<f64>().unwrap())
        })
        .collect()
}

/// Reads `Linear System Save Perm`: 0-based node -> 0-based dof, absent
/// (perm 0) nodes panic.
pub fn read_perm(path: &Path) -> Vec<usize> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter_map(|l| {
            let t: Vec<&str> = l.split_whitespace().collect();
            (t.len() >= 2).then(|| {
                let p: usize = t[1].parse().unwrap();
                assert!(p > 0, "node not in system");
                p - 1
            })
        })
        .collect()
}

/// Largest `|x - y|` over both matrices (missing entries are zero), and the
/// largest `|y|`.
pub fn matrix_diff(
    x: &BTreeMap<(usize, usize), f64>,
    y: &BTreeMap<(usize, usize), f64>,
) -> (f64, f64) {
    let mut diff = 0.0f64;
    let mut scale = 0.0f64;
    for (k, v) in x {
        diff = diff.max((v - y.get(k).copied().unwrap_or(0.0)).abs());
        scale = scale.max(v.abs());
    }
    for (k, v) in y {
        diff = diff.max((v - x.get(k).copied().unwrap_or(0.0)).abs());
        scale = scale.max(v.abs());
    }
    (diff, scale)
}

/// Solves `A x = b` for an SPD system with conjugate gradients on a
/// row-compressed copy. Returns `x`.
pub fn cg(n: usize, a: &BTreeMap<(usize, usize), f64>, b: &[f64]) -> Vec<f64> {
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (&(i, j), &v) in a {
        rows[i].push((j, v));
    }
    let mul = |x: &[f64]| -> Vec<f64> {
        rows.iter()
            .map(|r| r.iter().map(|&(j, v)| v * x[j]).sum())
            .collect()
    };
    let dot = |x: &[f64], y: &[f64]| -> f64 { x.iter().zip(y).map(|(p, q)| p * q).sum() };
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rs = dot(&r, &r);
    let tol = 1e-26 * dot(b, b).max(1e-300);
    for _ in 0..20 * n.max(10) {
        if rs < tol {
            break;
        }
        let ap = mul(&p);
        let alpha = rs / dot(&p, &ap);
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rs_new = dot(&r, &r);
        for i in 0..n {
            p[i] = r[i] + rs_new / rs * p[i];
        }
        rs = rs_new;
    }
    x
}
