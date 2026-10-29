//! `MatrixScaling.F90`: `ScaleLinearSystem`, `BackScaleLinearSystem` and
//! `ScaleLinearSystemVectors` for a serial matrix.
//!
//! The complex methods act on the native complex matrix: every real component
//! of an entry is scaled by the same real factor, which is what the interleaved
//! Fortran storage does to the (re, im) pair of a complex row.

use crate::scalar::{ElmerScalar, TINY};
use elmer_core::{CsrMatrix, LinalgError};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScalingMethod {
    /// `(D A D) y = D b`, `D = 1/sqrt(|diag A|)` (the default).
    Diagonal,
    /// Rows scaled to unit absolute row sum (`row equilibration`, `rowsum`).
    RowEquilibration,
    /// `A / mean|diag A|`.
    Constant,
    /// No scaling.
    None,
}

impl ScalingMethod {
    /// Parse `Linear System Scaling Method` (`Linear System Row Equilibration`
    /// maps to [`ScalingMethod::RowEquilibration`] before this is consulted).
    pub fn from_keyword(s: &str) -> Result<Self, LinalgError> {
        match s.trim().to_ascii_lowercase().as_str() {
            "diagonal" => Ok(Self::Diagonal),
            "row equilibration" | "rowsum" => Ok(Self::RowEquilibration),
            "constant" => Ok(Self::Constant),
            "none" => Ok(Self::None),
            other => Err(LinalgError::Unsupported(format!(
                "unknown scaling method: {other}"
            ))),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ScaleOptions {
    pub method: ScalingMethod,
    /// `ApplyScaling`: when false only the factors are computed.
    pub apply: bool,
    /// `RhsScaling`: normalise the right-hand side as well.
    pub rhs_scaling: bool,
}

impl Default for ScaleOptions {
    fn default() -> Self {
        ScaleOptions {
            method: ScalingMethod::Diagonal,
            apply: true,
            rhs_scaling: true,
        }
    }
}

/// The factors a scaling produced (`DiagScaling`, `RhsScaling`, `AveScaling`);
/// keep it to undo the scaling.
#[derive(Clone, Debug)]
pub struct RowScaling {
    method: ScalingMethod,
    diag: Vec<f64>,
    rhs_scaling: f64,
    ave_scaling: f64,
}

fn scale_vals<T: ElmerScalar>(a: &CsrMatrix<T>, vals: &mut [T], d: &[f64]) {
    for i in 0..a.nrows {
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            vals[j] = vals[j].scale_real(d[i] * d[a.col_idx[j] as usize]);
        }
    }
}

fn unscale_vals<T: ElmerScalar>(a: &CsrMatrix<T>, vals: &mut [T], d: &[f64]) {
    for i in 0..a.nrows {
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            vals[j] = vals[j].unscale_real(d[i] * d[a.col_idx[j] as usize]);
        }
    }
}

fn scale_rows<T: ElmerScalar>(a: &CsrMatrix<T>, vals: &mut [T], d: &[f64]) {
    for i in 0..a.nrows {
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            vals[j] = vals[j].scale_real(d[i]);
        }
    }
}

fn unscale_rows<T: ElmerScalar>(a: &CsrMatrix<T>, vals: &mut [T], d: &[f64]) {
    for i in 0..a.nrows {
        for j in a.row_ptr[i]..a.row_ptr[i + 1] {
            vals[j] = vals[j].unscale_real(d[i]);
        }
    }
}

/// `ScaleLinearSystem` for `a`, with optional right-hand side `b` and initial
/// guess `x`.
pub fn scale_system<T: ElmerScalar>(
    a: &mut CsrMatrix<T>,
    mut b: Option<&mut [T]>,
    mut x: Option<&mut [T]>,
    opts: &ScaleOptions,
) -> RowScaling {
    let n = a.nrows;
    let mut rs = RowScaling {
        method: opts.method,
        diag: Vec::new(),
        rhs_scaling: 1.0,
        ave_scaling: 1.0,
    };
    match opts.method {
        ScalingMethod::None => {}
        ScalingMethod::Diagonal => {
            rs.diag = T::diagonal_scaling(a);
            if !opts.apply {
                return rs;
            }
            let mut v = std::mem::take(&mut a.values);
            scale_vals(a, &mut v, &rs.diag);
            a.values = v;
            if let Some(b) = b.as_deref_mut() {
                for i in 0..n {
                    b[i] = b[i].scale_real(rs.diag[i]);
                }
                if opts.rhs_scaling {
                    let mut s = 0.0;
                    for bi in &b[..n] {
                        bi.for_each_component(&mut |c| s += c * c);
                    }
                    let bnorm = s.sqrt();
                    if bnorm >= TINY.sqrt() {
                        rs.rhs_scaling = bnorm;
                        for bi in b[..n].iter_mut() {
                            *bi = bi.unscale_real(bnorm);
                        }
                    }
                }
            }
            if let Some(x) = x.as_deref_mut() {
                for i in 0..n {
                    x[i] = x[i].unscale_real(rs.diag[i] * rs.rhs_scaling);
                }
            }
        }
        ScalingMethod::RowEquilibration => {
            let mut d = vec![0.0; n];
            for (i, di) in d.iter_mut().enumerate() {
                let mut tmp = 0.0;
                for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                    tmp += a.values[k].abs();
                }
                *di = tmp;
            }
            for di in d.iter_mut() {
                *di = if *di > TINY { 1.0 / *di } else { 1.0 };
            }
            rs.diag = d;
            if !opts.apply {
                return rs;
            }
            let mut v = std::mem::take(&mut a.values);
            scale_rows(a, &mut v, &rs.diag);
            a.values = v;
            if let Some(b) = b.as_deref_mut() {
                for i in 0..n {
                    b[i] = b[i].scale_real(rs.diag[i]);
                }
            }
        }
        ScalingMethod::Constant => {
            let mut diag_sum = 0.0;
            for i in 0..n {
                if let Some(j) = crate::scalar::diag_index(a, i) {
                    a.values[j].diag_abs_terms(&mut |t| diag_sum += t);
                }
            }
            let n_sum = (n * T::COMPONENTS) as f64;
            let ascl = diag_sum / n_sum;
            rs.ave_scaling = ascl;
            if !opts.apply || ascl.is_nan() || ascl <= 0.0 {
                return rs;
            }
            for v in a.values.iter_mut() {
                *v = v.unscale_real(ascl);
            }
            if let Some(b) = b {
                if opts.rhs_scaling {
                    let mut bsum = 0.0;
                    for bi in &b[..n] {
                        bi.for_each_component(&mut |c| bsum += c.abs());
                    }
                    let bscl = bsum / n_sum;
                    if bscl > 0.0 && bscl.is_finite() {
                        for bi in b[..n].iter_mut() {
                            *bi = bi.unscale_real(bscl);
                        }
                        rs.rhs_scaling = bscl;
                    }
                }
            }
            if let Some(x) = x {
                let xscl = rs.rhs_scaling / ascl;
                for xi in x[..n].iter_mut() {
                    *xi = xi.unscale_real(xscl);
                }
            }
        }
    }
    rs
}

impl RowScaling {
    pub fn method(&self) -> ScalingMethod {
        self.method
    }
    /// `DiagScaling` (diagonal and row-equilibration methods).
    pub fn diag(&self) -> &[f64] {
        &self.diag
    }
    pub fn rhs_scaling(&self) -> f64 {
        self.rhs_scaling
    }
    /// `AveScaling` (constant method).
    pub fn ave_scaling(&self) -> f64 {
        self.ave_scaling
    }

    /// Apply the scaling to a coefficient array sharing `a`'s structure
    /// (`PrecValues`, `MassValues`, `DampValues`).
    pub fn scale_values<T: ElmerScalar>(&self, a: &CsrMatrix<T>, vals: &mut [T]) {
        match self.method {
            ScalingMethod::Diagonal => scale_vals(a, vals, &self.diag),
            ScalingMethod::RowEquilibration => scale_rows(a, vals, &self.diag),
            ScalingMethod::Constant => {
                if self.ave_scaling > 0.0 {
                    for v in vals.iter_mut() {
                        *v = v.unscale_real(self.ave_scaling);
                    }
                }
            }
            ScalingMethod::None => {}
        }
    }

    /// Undo [`RowScaling::scale_values`].
    pub fn unscale_values<T: ElmerScalar>(&self, a: &CsrMatrix<T>, vals: &mut [T]) {
        match self.method {
            ScalingMethod::Diagonal => unscale_vals(a, vals, &self.diag),
            ScalingMethod::RowEquilibration => unscale_rows(a, vals, &self.diag),
            ScalingMethod::Constant => self.unscale_constant(vals),
            ScalingMethod::None => {}
        }
    }

    fn unscale_constant<T: ElmerScalar>(&self, vals: &mut [T]) {
        let ascl = self.ave_scaling;
        if ascl > 0.0 && (ascl - 1.0).abs() > f64::EPSILON {
            for v in vals.iter_mut() {
                *v = v.scale_real(ascl);
            }
        }
    }

    /// `BackScaleLinearSystem`: restore `a` and transform `b` and `x` back to
    /// the unscaled system.
    pub fn unscale<T: ElmerScalar>(
        &self,
        a: &mut CsrMatrix<T>,
        b: Option<&mut [T]>,
        x: Option<&mut [T]>,
    ) {
        let n = a.nrows;
        match self.method {
            ScalingMethod::None => {}
            ScalingMethod::Diagonal => {
                if let Some(x) = x {
                    for i in 0..n {
                        x[i] = x[i].scale_real(self.diag[i]).scale_real(self.rhs_scaling);
                    }
                }
                if let Some(b) = b {
                    for i in 0..n {
                        b[i] = b[i].unscale_real(self.diag[i]).scale_real(self.rhs_scaling);
                    }
                }
                let mut v = std::mem::take(&mut a.values);
                unscale_vals(a, &mut v, &self.diag);
                a.values = v;
            }
            ScalingMethod::RowEquilibration => {
                if let Some(b) = b {
                    for i in 0..n {
                        b[i] = b[i].unscale_real(self.diag[i]);
                    }
                }
                let mut v = std::mem::take(&mut a.values);
                unscale_rows(a, &mut v, &self.diag);
                a.values = v;
            }
            ScalingMethod::Constant => {
                let xscl = self.rhs_scaling / self.ave_scaling;
                if let Some(x) = x {
                    for xi in x[..n].iter_mut() {
                        *xi = xi.scale_real(xscl);
                    }
                }
                if let Some(b) = b {
                    for bi in b[..n].iter_mut() {
                        *bi = bi.scale_real(self.rhs_scaling);
                    }
                }
                let mut v = std::mem::take(&mut a.values);
                self.unscale_constant(&mut v);
                a.values = v;
            }
        }
    }

    /// `ScaleLinearSystemVectors`: transform `b` (and `x`) like the right-hand
    /// side and solution of the scaled system, or back (`back = true`).
    pub fn scale_vectors<T: ElmerScalar>(&self, b: &mut [T], x: Option<&mut [T]>, back: bool) {
        let n = b.len();
        let d = &self.diag;
        let r = self.rhs_scaling;
        match (self.method, back) {
            (ScalingMethod::None, _) => {}
            (ScalingMethod::Diagonal, false) => {
                for i in 0..n {
                    b[i] = b[i].scale_real(d[i]).unscale_real(r);
                }
                if let Some(x) = x {
                    for i in 0..n {
                        x[i] = x[i].unscale_real(d[i] * r);
                    }
                }
            }
            (ScalingMethod::Diagonal, true) => {
                for i in 0..n {
                    b[i] = b[i].unscale_real(d[i]).scale_real(r);
                }
                if let Some(x) = x {
                    for i in 0..n {
                        x[i] = x[i].scale_real(d[i]).scale_real(r);
                    }
                }
            }
            (ScalingMethod::RowEquilibration, false) => {
                for i in 0..n {
                    b[i] = b[i].scale_real(d[i]);
                }
            }
            (ScalingMethod::RowEquilibration, true) => {
                for i in 0..n {
                    b[i] = b[i].unscale_real(d[i]);
                }
            }
            (ScalingMethod::Constant, false) => {
                for bi in b.iter_mut() {
                    *bi = bi.unscale_real(r);
                }
                if let Some(x) = x {
                    for xi in x.iter_mut() {
                        *xi = xi.scale_real(self.ave_scaling).unscale_real(r);
                    }
                }
            }
            (ScalingMethod::Constant, true) => {
                for bi in b.iter_mut() {
                    *bi = bi.scale_real(r);
                }
                if let Some(x) = x {
                    for xi in x.iter_mut() {
                        *xi = xi.scale_real(r / self.ave_scaling);
                    }
                }
            }
        }
    }
}
