//! Generic transport solvers of Elmer: `AdvectionDiffusionSolver`
//! (`fem/src/modules/AdvectionDiffusion.F90`, scalar species transport with
//! diffusion, convection, Soret and Nernst-Planck migration) and the ModelPDE
//! `AdvDiffSolver` (`fem/src/modules/ModelPDE.F90`, the prototype
//! advection-diffusion-reaction equation).
//!
//! The kernels ([`advdiff_bulk_local`], [`advdiff_boundary_local`],
//! [`model_pde_bulk_local`], [`model_pde_boundary_local`]) are pure functions
//! of element geometry and nodal property values; time integration is left to
//! the caller (they return the mass matrix next to stiffness and force).
//! [`AdvDiffSolver`] and [`ModelPdeSolver`] wire them into `elmer-solver`'s
//! `Solver` trait over `elmer-fe`'s Lagrange bases, Elmer's Gauss rules and
//! `StabParam` (`hK`, `mK`); [`register`] binds both procedures.
//!
//! Supported: Cartesian coordinates; relative and absolute mass units;
//! constant and computed convection velocity; residual-free bubbles
//! (the AdvectionDiffusion default, also condensed in transient runs exactly
//! as Elmer does) and SUPG on linear simplex elements; Soret diffusion;
//! Nernst-Planck migration; mesh velocity; mass transfer and flux boundary
//! terms. Everything else is refused by [`check_advdiff_supported`] or
//! [`check_model_pde_supported`] instead of being approximated.
#![forbid(unsafe_code)]

mod advdiff_solver;
mod compose;
mod error;
mod fe;
mod model_pde;
mod model_pde_solver;
mod preflight;

pub use advdiff_solver::{advdiff_factory, register_advdiff, AdvDiffSolver};
pub use compose::{
    advdiff_boundary_local, advdiff_bulk_local, AdvDiffConfig, AdvDiffNodal, Stabilization,
};
pub use error::AdvDiffError;
pub use fe::{bubble_gauss_points, condensate, eval_with_bubbles, BasisHessian, Fe};
pub use model_pde::{
    model_pde_boundary_local, model_pde_bulk_local, ModelPdeBcNodal, ModelPdeNodal,
};
pub use model_pde_solver::{model_pde_factory, register_model_pde, ModelPdeSolver};
pub use preflight::{
    check_advdiff_supported, check_model_pde_supported, AdvDiffSections, ModelPdeSections,
    ADVDIFF_PROCEDURE, MODEL_PDE_PROCEDURE,
};

use elmer_solver::SolverRegistry;

/// Registers both `("AdvectionDiffusion", "AdvectionDiffusionSolver")` and
/// `("ModelPDE", "AdvDiffSolver")`.
pub fn register(registry: &mut SolverRegistry) {
    register_advdiff(registry);
    register_model_pde(registry);
}
