//! Gebhart factors of diffuse gray radiation: `CalculateGebhartFactors` of
//! `fem/src/RadiationFactors.F90` for the default (explicit, direct solver)
//! settings. The Fortran solves one sparse system per radiating surface with
//! UMFPACK; here the matrix is factored once, densely, and solved for every
//! right-hand side, which gives the same factors.

use crate::Error;

/// `AEPS = 10 * EPSILON(1._dp)` (Types.F90).
const AEPS: f64 = 10.0 * f64::EPSILON;

/// What `CalculateGebhartFactors` reads.
pub struct GebhartInput<'a> {
    /// The view factors as `ReadViewFactorsFromFile` leaves them: row `i`
    /// holds `(column, RelAreas(i) * factor)`; columns are 0-based surface
    /// numbers.
    pub rows: &'a [Vec<(usize, f64)>],
    /// Areas of the surfaces over the largest one (`RelAreas`).
    pub rel_areas: &'a [f64],
    pub emissivity: &'a [f64],
    /// `1 - absorptivity - transmissivity`.
    pub reflectivity: &'a [f64],
    /// `Minimum Gebhart Factor`, 1e-20 by default.
    pub min_factor: f64,
}

/// The factors and the sums the Fortran reports.
pub struct Gebhart {
    /// Row `t`: `(surface, factor)` of the surfaces whose factor exceeds the
    /// minimum, in surface order.
    pub rows: Vec<Vec<(usize, f64)>>,
    pub min_row_sum: f64,
    pub max_row_sum: f64,
}

/// The matrix `G` of `CalculateGebhartFactors` as `CreateRadiationMatrix`
/// lays it out: row `i` holds the columns of its view factors and the
/// diagonal, sorted (CRS, 0-based).
pub struct GebhartMatrix {
    pub row_ptr: Vec<usize>,
    pub cols: Vec<u32>,
    pub values: Vec<f64>,
}

impl GebhartMatrix {
    fn at(&mut self, i: usize, c: usize) -> &mut f64 {
        let span = self.row_ptr[i]..self.row_ptr[i + 1];
        let k = self.cols[span.clone()].binary_search(&(c as u32)).expect("entry in the pattern");
        &mut self.values[span.start + k]
    }
}

/// A direct solve of `G x = b` in place, `G` factored on the first call
/// (Elmer's `DirectSolver` with `Linear System Refactorize` false after it).
pub type GebhartSolve<'a> = dyn FnMut(&GebhartMatrix, &mut [f64]) -> Result<(), String> + 'a;

/// The Gebhart factors of `input`. See [`GebhartInput`]. The system is
/// solved by a dense LU; [`gebhart_factors_with`] takes the solver.
pub fn gebhart_factors(input: &GebhartInput<'_>) -> Result<Gebhart, Error> {
    let mut lu: Option<Lu> = None;
    let mut dense = |g: &GebhartMatrix, b: &mut [f64]| -> Result<(), String> {
        if lu.is_none() {
            let n = g.row_ptr.len() - 1;
            let mut a = vec![0.0; n * n];
            for i in 0..n {
                for k in g.row_ptr[i]..g.row_ptr[i + 1] {
                    a[i * n + g.cols[k] as usize] = g.values[k];
                }
            }
            lu = Some(Lu::new(n, a).map_err(|e| e.to_string())?);
        }
        lu.as_ref().expect("factored").solve(b);
        Ok(())
    };
    gebhart_factors_with(input, &mut dense)
}

/// [`gebhart_factors`] with the sparse `G` of `CreateRadiationMatrix`,
/// assembled and scaled in the Fortran's order, solved by `solve` (Elmer's
/// `DirectSolver`, UMFPACK by default).
pub fn gebhart_factors_with(input: &GebhartInput<'_>, solve: &mut GebhartSolve<'_>) -> Result<Gebhart, Error> {
    let n = input.rows.len();
    let (rel, emis, refl) = (input.rel_areas, input.emissivity, input.reflectivity);
    if rel.len() != n || emis.len() != n || refl.len() != n {
        return Err(Error::Invalid("per-surface arrays differ in length".into()));
    }

    let triv = refl.iter().all(|r| r.abs() <= AEPS);
    let symm = refl.iter().all(|r| r.abs() > AEPS) || triv;

    // The matrix, scaled to unit diagonal: CRS_AddToMatrixElement of the
    // view factors and then of the diagonal, `G(i,j) = G(i,j)*Diag(i)*Diag(j)`.
    let mut row_ptr = vec![0usize];
    let mut cols: Vec<u32> = Vec::new();
    for (i, row) in input.rows.iter().enumerate() {
        let mut c: Vec<u32> = row.iter().map(|&(c, _)| c as u32).chain(std::iter::once(i as u32)).collect();
        c.sort_unstable();
        c.dedup();
        cols.extend(c);
        row_ptr.push(cols.len());
    }
    let mut g = GebhartMatrix { row_ptr, values: vec![0.0; cols.len()], cols };
    let mut diag = vec![0.0; n];
    if !triv {
        for i in 0..n {
            let r = if symm { refl[i] } else { 1.0 };
            for &(c, v) in &input.rows[i] {
                let s = if symm { refl[i] * refl[c] } else { refl[i] };
                *g.at(i, c) += -s * v;
            }
            diag[i] = r * rel[i];
            *g.at(i, i) += diag[i];
        }
    }
    let scale: Vec<f64> = diag
        .iter()
        .map(|d| (1.0 / d.abs().max(1.0e-12)).sqrt())
        .collect();
    if !triv {
        for i in 0..n {
            for k in g.row_ptr[i]..g.row_ptr[i + 1] {
                g.values[k] = g.values[k] * scale[i] * scale[g.cols[k] as usize];
            }
        }
    }

    let mut rows = Vec::with_capacity(n);
    let mut row_sums = vec![0.0; n];
    let mut fac = vec![0.0; n];
    for t in 0..n {
        if triv {
            fac.iter_mut().for_each(|f| *f = 0.0);
            for &(c, v) in &input.rows[t] {
                fac[c] = v / rel[t];
            }
            for i in 0..n {
                row_sums[i] += fac[i] * rel[t] / (rel[i] * emis[i]);
            }
        } else {
            let mut sol = vec![0.0; n];
            sol[t] = scale[t];
            solve(&g, &mut sol).map_err(Error::Invalid)?;
            for i in 0..n {
                sol[i] *= scale[i];
            }
            for i in 0..n {
                let mut s = 0.0;
                for &(c, v) in &input.rows[i] {
                    s += if symm {
                        refl[c] * v * sol[c]
                    } else {
                        v * sol[c]
                    };
                }
                fac[i] = s * emis[t] * emis[i];
                row_sums[i] += fac[i] * rel[t] / (rel[i] * emis[i]);
            }
        }
        rows.push(
            (0..n)
                .filter(|&i| fac[i] > input.min_factor)
                .map(|i| (i, fac[i]))
                .collect(),
        );
    }

    Ok(Gebhart {
        rows,
        min_row_sum: row_sums.iter().copied().fold(f64::INFINITY, f64::min),
        max_row_sum: row_sums.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    })
}

/// Dense LU with partial pivoting.
struct Lu {
    n: usize,
    lu: Vec<f64>,
    piv: Vec<usize>,
}

impl Lu {
    fn new(n: usize, mut lu: Vec<f64>) -> Result<Self, Error> {
        let mut piv = vec![0; n];
        for c in 0..n {
            let mut p = c;
            for r in c + 1..n {
                if lu[r * n + c].abs() > lu[p * n + c].abs() {
                    p = r;
                }
            }
            if lu[p * n + c] == 0.0 {
                return Err(Error::Singular(c));
            }
            piv[c] = p;
            if p != c {
                for j in 0..n {
                    lu.swap(c * n + j, p * n + j);
                }
            }
            for r in c + 1..n {
                let f = lu[r * n + c] / lu[c * n + c];
                lu[r * n + c] = f;
                for j in c + 1..n {
                    lu[r * n + j] -= f * lu[c * n + j];
                }
            }
        }
        Ok(Self { n, lu, piv })
    }

    fn solve(&self, b: &mut [f64]) {
        let n = self.n;
        for c in 0..n {
            b.swap(c, self.piv[c]);
            for r in c + 1..n {
                b[r] -= self.lu[r * n + c] * b[c];
            }
        }
        for r in (0..n).rev() {
            for j in r + 1..n {
                b[r] -= self.lu[r * n + j] * b[j];
            }
            b[r] /= self.lu[r * n + r];
        }
    }
}
