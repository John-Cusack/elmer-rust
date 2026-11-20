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
/// residual-free bubble otherwise: then the element is a p element with
/// bubble functions ([`HeatBulkConfig::extra_functions`]) and no
/// stabilization term, and the solver condenses the bubbles.
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
    /// Basis functions after the element's nodal ones: a p element's edge,
    /// face and bubble functions (`nd - n` of `LocalMatrixVec`). Only
    /// [`crate::heat_bulk_local_vec`] takes them; the local system then has
    /// `n + extra_functions` rows, and the nodal inputs are interpolated
    /// with the first `n` functions.
    pub extra_functions: usize,
}

impl HeatBulkConfig {
    /// Steady, Cartesian, unstabilised, nodal.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            axisymmetric: false,
            transient: false,
            stabilization: Stabilization::None,
            quad_points: None,
            extra_functions: 0,
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
    /// Convection velocity. `Some` turns convection on (`Convection =
    /// Constant` or `Computed`); components beyond `dim` are ignored.
    pub convection_velocity: Option<Vec<[f64; 3]>>,
    /// `Convection = Computed`: the velocity is a solution field, which
    /// [`crate::heat_bulk_local_vec`] always interpolates, where a
    /// `Convection Velocity i` constant over the element is used as given.
    pub convection_computed: bool,
    /// The `Mesh Velocity` field, subtracted from the convection velocity
    /// (`MeshVeloField_h`).
    pub mesh_velocity: Option<Vec<[f64; 3]>>,
    /// `Heat Transfer Multiplier`, scaling the whole local system (for
    /// lower-dimensional elements embedded in a larger mesh).
    pub heat_transfer_multiplier: Option<Vec<f64>>,
    /// The material's `Compressibility Model`; `None` is incompressible.
    pub compressibility: Option<Compressibility>,
    /// `Joule Heat = True` in the body force, with the field it is computed
    /// from; `None` adds no Joule heat.
    pub joule_heat: Option<JouleHeat>,
}

/// The Joule heat `JouleHeat` (Differentials.F90:284) computes at each
/// integration point, which `LocalMatrixVec` adds to the force after the
/// source. `Electric Conductivity` and the field are nodal values, always
/// interpolated.
#[derive(Clone, Debug, PartialEq)]
pub enum JouleHeat {
    /// The `Joule Field` variable: `sigma * field`.
    Field {
        electric_conductivity: Vec<f64>,
        field: Vec<f64>,
    },
    /// The `Potential` variable: `sigma * |grad V|^2`.
    Potential {
        electric_conductivity: Vec<f64>,
        potential: Vec<f64>,
    },
}

/// A `Compressibility Model` of `LocalMatrixVec` (HeatSolve.F90:1300): the
/// density is replaced at each integration point, and with computed
/// convection the reversible pressure work `Pcoeff*(dp/dt + v.grad p)` is
/// added to the force.
#[derive(Clone, Debug, PartialEq)]
pub struct Compressibility {
    pub model: CompressModel,
    /// The solver's own `Temperature`, the current iterate, at the nodes.
    pub temperature: Vec<f64>,
    /// The flow solution's pressure (its last component) at the nodes. Set
    /// only with `Convection = Computed`; it turns the pressure work on.
    pub pressure: Option<Vec<f64>>,
    /// The flow pressure of the previous time step at the nodes, with the
    /// step size: the `dp/dt` part of the pressure work in a transient run.
    pub prev_pressure: Option<(Vec<f64>, f64)>,
}

/// The density law of a [`Compressibility`]. A keyword absent from the SIF
/// is `None` and takes Elmer's default.
#[derive(Clone, Debug, PartialEq)]
pub enum CompressModel {
    /// `Thermal`: `rho * (1 - beta*(T - Tref))` with `Heat Expansion
    /// Coefficient` beta and `Reference Temperature` Tref (both zero when
    /// absent).
    Thermal {
        heat_expansion: Option<Vec<f64>>,
        reference_temperature: Option<Vec<f64>>,
    },
    /// `Perfect Gas`: `(p + p0) / ((gamma-1)/gamma * cp * T)` with `Specific
    /// Heat Ratio` gamma (default 5/3) and `Reference Pressure` p0 (default
    /// 0), both read at the first integration point. Needs
    /// [`Compressibility::pressure`].
    PerfectGas {
        specific_heat_ratio: Option<Vec<f64>>,
        reference_pressure: Option<Vec<f64>>,
    },
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
    /// `Heat Flux` at the boundary's Gauss points instead of `heat_flux`,
    /// when it depends on an integration-point variable (`HeatFlux_h`
    /// evaluated at each point, HeatSolve.F90:2445).
    pub heat_flux_at_ips: Option<Vec<f64>>,
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
