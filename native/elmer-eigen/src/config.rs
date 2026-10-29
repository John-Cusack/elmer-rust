//! `Eigen System *` keyword semantics of `ArpackEigenSolve` (EigenSolve.F90:140),
//! `ArpackSetWhich` (EigenSolve.F90:2296) and the normalisation switches read by
//! `SolveCore.F90:1439`.

use elmer_core::{LinalgError, Params};

/// Part of the spectrum, in terms of the eigenvalues `lambda` of `K x = lambda M x`.
/// The six `Eigen System Select` / `Eigen System Sorting` strings.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum Which {
    #[default]
    SmallestMagnitude,
    LargestMagnitude,
    SmallestReal,
    LargestReal,
    SmallestImag,
    LargestImag,
}

impl Which {
    /// Elmer's `SELECT CASE`: an absent or unrecognised string takes the
    /// `CASE DEFAULT` branch, smallest magnitude.
    pub fn from_keyword(s: Option<&str>) -> Which {
        match s.map(|s| s.trim().to_ascii_lowercase()).as_deref() {
            Some("largest magnitude") => Which::LargestMagnitude,
            Some("smallest real part") => Which::SmallestReal,
            Some("largest real part") => Which::LargestReal,
            Some("smallest imag part") => Which::SmallestImag,
            Some("largest imag part") => Which::LargestImag,
            _ => Which::SmallestMagnitude,
        }
    }
}

/// ARPACK's `WHICH` for the operator the iteration actually sees.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ArpackWhich {
    LM,
    SM,
    LR,
    SR,
    LI,
    SI,
}

impl ArpackWhich {
    /// `ArpackSetWhich`: mode 2 (lumped mass, `OP = M^-1 K`) keeps the natural
    /// order; mode 3 (shift-invert, `OP = (K - sigma M)^-1 M`, eigenvalues
    /// `1/(lambda - sigma)`) inverts it.
    pub fn for_mode(select: Which, shift_invert: bool) -> ArpackWhich {
        use ArpackWhich::*;
        match (select, shift_invert) {
            (Which::SmallestMagnitude, false) | (Which::LargestMagnitude, true) => SM,
            (Which::LargestMagnitude, false) | (Which::SmallestMagnitude, true) => LM,
            (Which::SmallestReal, false) | (Which::LargestReal, true) => SR,
            (Which::LargestReal, false) | (Which::SmallestReal, true) => LR,
            (Which::SmallestImag, false) | (Which::LargestImag, true) => SI,
            (Which::LargestImag, false) | (Which::SmallestImag, true) => LI,
        }
    }
}

/// Eigenvector scaling applied after the solve (`ScaleEigenVectors`,
/// EigenSolve.F90:622).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum Normalize {
    /// Vectors as the iteration returns them: unit `M`-norm, arbitrary sign.
    #[default]
    None,
    /// `Eigen System Mass Normalize`: rescale to `x^H M x = 1`.
    Mass,
    /// `Eigen System Normalize To Unity`: largest component magnitude 1.
    Unity,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EigenConfig {
    /// `Eigen System Values`.
    pub nev: usize,
    /// `Eigen System Lanczos Vectors`, default `3*nev + 1`.
    pub ncv: usize,
    /// `Eigen System Select`.
    pub select: Which,
    /// `Eigen System Sorting`: order of the returned pairs.
    pub sorting: Which,
    /// `Eigen System Shift`; only used with a consistent mass matrix.
    pub shift: f64,
    /// `Eigen System Convergence Tolerance`, default
    /// `100 * Linear System Convergence Tolerance`. `<= 0` means machine epsilon.
    pub tol: f64,
    /// `Eigen System Max Iterations` (restarts), default 300.
    pub max_iter: usize,
    /// `Matrix % Symmetric`: thick-restart Lanczos (`DSAUPD`) when true,
    /// Arnoldi (`DNAUPD`) otherwise.
    pub symmetric: bool,
    pub normalize: Normalize,
    /// `Eigen System Compute Residuals` (`CheckResiduals`).
    pub compute_residuals: bool,
}

impl EigenConfig {
    /// Defaults of `ArpackEigenSolve` for `nev` eigenpairs.
    pub fn new(nev: usize, symmetric: bool) -> Self {
        EigenConfig {
            nev,
            ncv: 3 * nev + 1,
            select: Which::default(),
            sorting: Which::default(),
            shift: 0.0,
            tol: 0.0,
            max_iter: 300,
            symmetric,
            normalize: Normalize::None,
            compute_residuals: false,
        }
    }

    /// Read a solver section. `symmetric` is the matrix flag the caller
    /// assembled with. The damped (`Eigen System Damped`) and stability
    /// (`Stability Analysis`) branches of `ArpackEigenSolve` are refused.
    pub fn from_params(p: &dyn Params, symmetric: bool) -> Result<Self, LinalgError> {
        if p.logical("eigen system damped") == Some(true) {
            return Err(LinalgError::Unsupported(
                "Eigen System Damped (ArpackDampedEigenSolve)".into(),
            ));
        }
        if p.logical("stability analysis") == Some(true) {
            return Err(LinalgError::Unsupported(
                "Stability Analysis (ArpackStabEigenSolve)".into(),
            ));
        }
        let nev = match p.integer("eigen system values") {
            Some(n) if n > 0 => n as usize,
            _ => {
                return Err(LinalgError::Dim(
                    "Eigen Analysis needs a positive 'Eigen System Values'".into(),
                ))
            }
        };
        let mut c = EigenConfig::new(nev, symmetric);
        if let Some(ncv) = p.integer("eigen system lanczos vectors") {
            c.ncv = ncv.max(0) as usize;
        }
        // ListGetConstReal without Found returns 0 for a missing keyword.
        c.tol = match p.real("eigen system convergence tolerance") {
            Some(t) => t,
            None => 100.0 * p.real("linear system convergence tolerance").unwrap_or(0.0),
        };
        if let Some(m) = p.integer("eigen system max iterations") {
            c.max_iter = m.max(0) as usize;
        }
        c.shift = p.real("eigen system shift").unwrap_or(0.0);
        c.select = Which::from_keyword(p.string("eigen system select"));
        c.sorting = Which::from_keyword(p.string("eigen system sorting"));
        c.normalize = if p.logical("eigen system normalize to unity") == Some(true) {
            Normalize::Unity
        } else if p.logical("eigen system mass normalize") == Some(true) {
            Normalize::Mass
        } else {
            Normalize::None
        };
        c.compute_residuals = p.logical("eigen system compute residuals") == Some(true);
        Ok(c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::{ParamMap, Value};

    fn map(entries: &[(&str, Value)]) -> ParamMap {
        let mut m = ParamMap::default();
        for (k, v) in entries {
            m.insert(k, v.clone());
        }
        m
    }

    #[test]
    fn defaults_follow_arpack_eigen_solve() {
        let p = map(&[
            ("Eigen System Values", Value::Integer(vec![10])),
            ("Linear System Convergence Tolerance", Value::real(1e-8)),
        ]);
        let c = EigenConfig::from_params(&p, true).unwrap();
        assert_eq!(c.nev, 10);
        assert_eq!(c.ncv, 31);
        assert_eq!(c.max_iter, 300);
        assert!((c.tol - 1e-6).abs() < 1e-20);
        assert_eq!(c.shift, 0.0);
        assert_eq!(c.select, Which::SmallestMagnitude);
        assert_eq!(c.sorting, Which::SmallestMagnitude);
        assert_eq!(c.normalize, Normalize::None);
    }

    #[test]
    fn explicit_keywords() {
        let p = map(&[
            ("Eigen System Values", Value::Integer(vec![4])),
            ("Eigen System Lanczos Vectors", Value::Integer(vec![20])),
            ("Eigen System Convergence Tolerance", Value::real(1e-9)),
            ("Eigen System Max Iterations", Value::Integer(vec![7])),
            ("Eigen System Shift", Value::real(1e3)),
            (
                "Eigen System Select",
                Value::Str("largest real part".into()),
            ),
            (
                "Eigen System Sorting",
                Value::Str("largest magnitude".into()),
            ),
            ("Eigen System Mass Normalize", Value::Logical(true)),
        ]);
        let c = EigenConfig::from_params(&p, false).unwrap();
        assert_eq!((c.ncv, c.max_iter, c.shift), (20, 7, 1e3));
        assert_eq!(c.tol, 1e-9);
        assert_eq!(c.select, Which::LargestReal);
        assert_eq!(c.sorting, Which::LargestMagnitude);
        assert_eq!(c.normalize, Normalize::Mass);
    }

    #[test]
    fn unity_wins_over_mass_and_unknown_select_is_default() {
        let p = map(&[
            ("Eigen System Values", Value::Integer(vec![2])),
            ("Eigen System Normalize To Unity", Value::Logical(true)),
            ("Eigen System Mass Normalize", Value::Logical(true)),
            (
                "Eigen System Select",
                Value::Str("smallest algebraic".into()),
            ),
        ]);
        let c = EigenConfig::from_params(&p, true).unwrap();
        assert_eq!(c.normalize, Normalize::Unity);
        assert_eq!(c.select, Which::SmallestMagnitude);
        // no tolerance anywhere: 100 * 0, i.e. machine precision
        assert_eq!(c.tol, 0.0);
    }

    #[test]
    fn refused_branches() {
        let damped = map(&[
            ("Eigen System Values", Value::Integer(vec![2])),
            ("Eigen System Damped", Value::Logical(true)),
        ]);
        assert!(matches!(
            EigenConfig::from_params(&damped, true),
            Err(LinalgError::Unsupported(_))
        ));
        let stab = map(&[
            ("Eigen System Values", Value::Integer(vec![2])),
            ("Stability Analysis", Value::Logical(true)),
        ]);
        assert!(matches!(
            EigenConfig::from_params(&stab, true),
            Err(LinalgError::Unsupported(_))
        ));
        assert!(EigenConfig::from_params(&ParamMap::default(), true).is_err());
    }

    #[test]
    fn set_which_inverts_for_shift_invert() {
        use ArpackWhich::*;
        let all = [
            (Which::SmallestMagnitude, SM, LM),
            (Which::LargestMagnitude, LM, SM),
            (Which::SmallestReal, SR, LR),
            (Which::LargestReal, LR, SR),
            (Which::SmallestImag, SI, LI),
            (Which::LargestImag, LI, SI),
        ];
        for (w, lumped, si) in all {
            assert_eq!(ArpackWhich::for_mode(w, false), lumped);
            assert_eq!(ArpackWhich::for_mode(w, true), si);
        }
    }
}
