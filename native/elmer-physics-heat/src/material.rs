//! Nodal inputs of the heat kernels.
//!
//! Every per-node vector has one entry per element node (`xyz.len()`), in
//! Elmer local node order. A property absent from the SIF is `None`, which the
//! kernels treat as Elmer does (`ListGetElementReal` returns not-found, and
//! the term using it drops out) rather than as an implicit default.

use crate::error::HeatError;

/// `Heat Conductivity`, by the shape it was given in the SIF (`Fdim` in
/// Lists.F90: 1x1 is rank 0, an n-by-1 column is rank 1, n-by-n is rank 2).
#[derive(Clone, Debug, PartialEq)]
pub enum Conductivity {
    /// `Heat Conductivity = k`.
    Scalar(Vec<f64>),
    /// `Heat Conductivity(3) = kxx kyy kzz`; only the first `dim` are used.
    Diagonal(Vec<[f64; 3]>),
    /// `Heat Conductivity(3,3) = ...`, row-major `k[i][j]`: the stiffness
    /// couples `dphi_p/dx_i` with `dphi_q/dx_j`.
    Full(Vec<[[f64; 3]; 3]>),
}

impl Conductivity {
    pub(crate) fn len(&self) -> usize {
        match self {
            Conductivity::Scalar(v) => v.len(),
            Conductivity::Diagonal(v) => v.len(),
            Conductivity::Full(v) => v.len(),
        }
    }
}

/// How convection is stabilised. Elmer picks SUPG when the solver has
/// `Stabilize = True` (or `Stabilization Method = stabilized`) and a
/// residual-free bubble otherwise; the bubble path is not implemented here.
#[derive(Copy, Clone, Debug, PartialEq, Default)]
pub enum Stabilization {
    /// Plain Galerkin weighting. Elmer never runs this on convected elements
    /// (it adds a bubble instead), so it is for tests and diffusion-dominated
    /// cases only.
    #[default]
    None,
    /// Equal-order SUPG with the Franca parameter built from `Element % hK`
    /// and `Element % StabilizationMK`, which the caller supplies.
    Supg { h_k: f64, m_k: f64 },
}

/// Element-level switches for [`crate::heat_bulk_local`].
#[derive(Clone, Debug, PartialEq)]
pub struct HeatBulkConfig {
    /// Coordinate-system dimension (`CoordinateSystemDimension()`), 1 to 3.
    pub dim: usize,
    /// Axisymmetric or cylindric coordinates: the measure becomes `r dr dz`
    /// with `r` the first coordinate.
    pub axisymmetric: bool,
    /// Build the mass matrix (`rho*cp`, the time-derivative coefficient).
    pub transient: bool,
    pub stabilization: Stabilization,
    /// Requested Gauss points, or `None` for the element default.
    pub quad_points: Option<usize>,
}

impl HeatBulkConfig {
    /// Steady, Cartesian, unstabilised.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            axisymmetric: false,
            transient: false,
            stabilization: Stabilization::None,
            quad_points: None,
        }
    }
}

/// Material and body-force values at the element nodes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HeatMaterialAtNodes {
    /// `Heat Conductivity`. `None` drops the diffusion term (and the SUPG
    /// term with it), as `HaveCond = .FALSE.` does in `LocalMatrixVec`.
    pub conductivity: Option<Conductivity>,
    /// `Density`. Missing reads as zero, so `rho*cp` and `Heat Source*rho`
    /// vanish.
    pub density: Option<Vec<f64>>,
    /// `Heat Capacity`. Only read when the element is transient or convected.
    pub heat_capacity: Option<Vec<f64>>,
    /// Body force `Heat Source`, per unit mass: the force is `rho * source`.
    pub heat_source: Option<Vec<f64>>,
    /// Body force `Volumetric Heat Source`. When present it replaces
    /// `Heat Source` and is not scaled by density.
    pub volumetric_heat_source: Option<Vec<f64>>,
    /// Convection velocity, already minus any mesh velocity. `Some` turns
    /// convection on (`Convection = Constant` or `Computed`); components
    /// beyond `dim` are ignored.
    pub convection_velocity: Option<Vec<[f64; 3]>>,
    /// `Heat Transfer Multiplier`, scaling the whole local system (for
    /// lower-dimensional elements embedded in a larger mesh).
    pub heat_transfer_multiplier: Option<Vec<f64>>,
}

/// Element-level switches for [`crate::heat_boundary_local`].
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct HeatBoundaryConfig {
    /// Axisymmetric or cylindric coordinates: the measure gains a factor `r`,
    /// the first coordinate.
    pub axisymmetric: bool,
    /// Requested Gauss points, or `None` for the element default.
    pub quad_points: Option<usize>,
}

/// Boundary-condition values at the nodes of a boundary element.
///
/// Elmer assembles `Heat Flux` as a force and `Heat Transfer Coefficient`
/// with `External Temperature` as a Robin term `c*(T - Text)`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HeatBcAtNodes {
    pub heat_flux: Option<Vec<f64>>,
    pub heat_transfer_coefficient: Option<Vec<f64>>,
    /// Read only together with the coefficient; absent reads as zero.
    pub external_temperature: Option<Vec<f64>>,
}

pub(crate) fn check_len(what: &'static str, expected: usize, got: usize) -> Result<(), HeatError> {
    if expected == got {
        Ok(())
    } else {
        Err(HeatError::NodalLength {
            what,
            expected,
            got,
        })
    }
}

pub(crate) fn check_opt(
    what: &'static str,
    n: usize,
    v: &Option<Vec<f64>>,
) -> Result<(), HeatError> {
    match v {
        Some(v) => check_len(what, n, v.len()),
        None => Ok(()),
    }
}

/// `sum_i basis[i] * nodal[i]`, the interpolation `ListGetElementReal` does.
pub(crate) fn interp(basis: &[f64], nodal: &[f64]) -> f64 {
    basis.iter().zip(nodal).map(|(b, v)| b * v).sum()
}
