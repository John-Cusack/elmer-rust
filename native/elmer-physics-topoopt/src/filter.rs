//! The element filters of TopoOpt.F90: `CreateSimpleFilter`,
//! `CreateDistanceFilter`, `NormalizeFilter` and `ApplyFilter`, on a CRS
//! matrix indexed by bulk element.

use std::collections::{BTreeMap, HashMap};

use elmer_core::{ElemId, MeshView};
use elmer_solver::SolverError;

/// A CRS matrix as `List_ToCRSMatrix` leaves it (columns ascending,
/// `ndeg = -1`).
#[derive(Clone, Debug, Default)]
pub(crate) struct Crs {
    pub rows: Vec<usize>,
    pub cols: Vec<usize>,
    pub values: Vec<f64>,
}

impl Crs {
    fn from_list(list: Vec<BTreeMap<usize, f64>>) -> Self {
        let mut c = Crs {
            rows: vec![0],
            ..Default::default()
        };
        for row in list {
            for (k, v) in row {
                c.cols.push(k);
                c.values.push(v);
            }
            c.rows.push(c.cols.len());
        }
        c
    }

    pub fn n_rows(&self) -> usize {
        self.rows.len() - 1
    }

    /// `CRS_MatrixVectorMultiply` with `ndeg = -1`: four partial sums.
    fn mv(&self, u: &[f64], v: &mut [f64]) {
        for i in 0..self.n_rows() {
            let (r0, r1) = (self.rows[i], self.rows[i + 1]);
            let k = r1 - r0;
            let (mut s1, mut s2, mut s3, mut s4) = (0.0, 0.0, 0.0, 0.0);
            let mut j = r0;
            while j < r0 + 4 * (k / 4) {
                s1 += u[self.cols[j]] * self.values[j];
                s2 += u[self.cols[j + 1]] * self.values[j + 1];
                s3 += u[self.cols[j + 2]] * self.values[j + 2];
                s4 += u[self.cols[j + 3]] * self.values[j + 3];
                j += 4;
            }
            for j in r0 + 4 * (k / 4)..r1 {
                s1 += u[self.cols[j]] * self.values[j];
            }
            v[i] = s1 + s2 + s3 + s4;
        }
    }

    /// `CRS_TransposeMatrixVectorMultiply`.
    fn mv_t(&self, u: &[f64], v: &mut [f64]) {
        let m = self.cols.iter().copied().max().map_or(0, |c| c + 1);
        v[..m].iter_mut().for_each(|x| *x = 0.0);
        for i in 0..self.n_rows() {
            for j in self.rows[i]..self.rows[i + 1] {
                let k = self.cols[j];
                v[k] += u[i] * self.values[j];
            }
        }
    }

    /// `NormalizeFilter`: unit row sums.
    pub fn normalize(&mut self) {
        for i in 0..self.n_rows() {
            let r = self.rows[i]..self.rows[i + 1];
            let mut rsum = 0.0;
            for j in r.clone() {
                rsum += self.values[j];
            }
            for j in r {
                self.values[j] /= rsum;
            }
        }
    }
}

/// The corner-node edges (2D) or faces (3D) of a bulk element, sorted.
fn sides(etype: u16, nodes: &[u32]) -> Result<Vec<Vec<u32>>, SolverError> {
    let faces: &[&[usize]] = match etype / 100 {
        3 => &[&[0, 1], &[1, 2], &[2, 0]],
        4 => &[&[0, 1], &[1, 2], &[2, 3], &[3, 0]],
        5 => &[&[0, 1, 2], &[0, 1, 3], &[1, 2, 3], &[0, 2, 3]],
        8 => &[
            &[0, 1, 2, 3],
            &[4, 5, 6, 7],
            &[0, 1, 5, 4],
            &[1, 2, 6, 5],
            &[2, 3, 7, 6],
            &[3, 0, 4, 7],
        ],
        _ => {
            return Err(SolverError::Unsupported(format!(
                "TopoOpt: element filter on element type {etype}"
            )))
        }
    };
    Ok(faces
        .iter()
        .map(|f| {
            let mut s: Vec<u32> = f.iter().map(|&k| nodes[k]).collect();
            s.sort_unstable();
            s
        })
        .collect())
}

/// `CreateSimpleFilter`: for every edge (2D) or face (3D) shared by two
/// bulk elements, ones at the pair and on both diagonals.
pub(crate) fn simple_filter(mesh: &dyn MeshView) -> Result<Crs, SolverError> {
    let nb = mesh.n_bulk();
    let mut owners: HashMap<Vec<u32>, Vec<usize>> = HashMap::new();
    for e in 0..nb {
        let el = mesh.element(ElemId(e as u32));
        for s in sides(el.etype.0, el.nodes)? {
            owners.entry(s).or_default().push(e);
        }
    }
    let mut list: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); nb];
    if nb > 0 {
        list[nb - 1].insert(nb - 1, 0.0);
    }
    let mut add = |i: usize, j: usize| {
        *list[i].entry(j).or_insert(0.0) += 1.0;
    };
    for els in owners.values() {
        if let [a, b] = els[..] {
            add(a, b);
            add(b, a);
            add(a, a);
            add(b, b);
        }
    }
    Ok(Crs::from_list(list))
}

/// `CreateDistanceFilter`: each element's neighbours (through the simple
/// filter) within `rfilter` of its centre, weighted `rfilter - distance`.
pub(crate) fn distance_filter(mesh: &dyn MeshView, rfilter: f64) -> Result<Crs, SolverError> {
    let dim = mesh.dim();
    let nb = mesh.n_bulk();
    let rfilter2 = rfilter * rfilter;
    let emat = simple_filter(mesh)?;
    let centers: Vec<[f64; 3]> = (0..nb)
        .map(|e| {
            let el = mesh.element(ElemId(e as u32));
            let n = el.nodes.len() as f64;
            let mut c = [0.0f64; 3];
            for k in 0..dim {
                let mut s = 0.0;
                for &i in el.nodes {
                    s += mesh.coord(i)[k];
                }
                c[k] = s / n;
            }
            c
        })
        .collect();
    let mut list: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); nb];
    if nb > 0 {
        list[nb - 1].insert(nb - 1, 0.0);
    }
    for i in 0..nb {
        let mut inds = vec![i];
        let mut dist = vec![rfilter];
        let mut k1 = 0;
        loop {
            let k2 = inds.len();
            for k in k1..k2 {
                let row = inds[k];
                for j in emat.rows[row]..emat.rows[row + 1] {
                    let i2 = emat.cols[j];
                    if inds.contains(&i2) {
                        continue;
                    }
                    let mut rik2 = 0.0;
                    for d in 0..dim {
                        let x = centers[i][d] - centers[i2][d];
                        rik2 += x * x;
                    }
                    if rik2 < rfilter2 {
                        if inds.len() >= 1000 {
                            return Err(SolverError::fatal("TopoOpt", "Too many neighbours!"));
                        }
                        inds.push(i2);
                        dist.push(rfilter - rik2.sqrt());
                    }
                }
            }
            if inds.len() == k2 {
                break;
            }
            k1 = k2;
        }
        for (k, d) in inds.into_iter().zip(dist) {
            *list[i].entry(k).or_insert(0.0) += d;
        }
    }
    Ok(Crs::from_list(list))
}

/// `ApplyFilter( Fmat, x, xf, niter, w, wmin, Trans )`: `niter` products
/// with the filter (or its transpose), of `x * w` scaled back by `max(w,
/// wmin)` when weighted. `xf` keeps its entries the products do not reach.
pub(crate) fn apply_filter(
    f: &Crs,
    x: &[f64],
    xf: &mut [f64],
    niter: usize,
    w: Option<(&[f64], f64)>,
    trans: bool,
) {
    let mut xtmp: Vec<f64> = match w {
        Some((w, _)) => x.iter().zip(w).map(|(a, b)| a * b).collect(),
        None => x.to_vec(),
    };
    for i in 1..=niter {
        if trans {
            f.mv_t(&xtmp, xf);
        } else {
            f.mv(&xtmp, xf);
        }
        if i < niter {
            xtmp.copy_from_slice(xf);
        }
    }
    if let Some((w, wmin)) = w {
        for (v, wi) in xf.iter_mut().zip(w) {
            *v /= wi.max(wmin);
        }
    }
}
