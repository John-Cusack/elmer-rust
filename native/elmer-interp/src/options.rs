//! Keyword-controlled behaviour of the interpolation, read from the
//! `Simulation` (and, where the Fortran does, the `Solver`) section.

use elmer_core::Params;

/// Settings of `InterpolateMeshToMeshQ`. Defaults are the Fortran's.
#[derive(Clone, Debug)]
pub struct InterpOptions {
    /// `Interpolation Global Epsilon` (default `2.0e-10`).
    pub global_epsilon: f64,
    /// `Interpolation Local Epsilon` (default `1.0e-10`).
    pub local_epsilon: f64,
    /// `Interpolation Max Iterations`: tolerance-widening attempts (default 12).
    pub max_iterations: usize,
    /// `Interpolation Numeric Epsilon` (default `1.0e-10`).
    pub numeric_epsilon: f64,
    /// `Interpolation Passive Coordinate`: 1..=3 flattens that axis, 0 is off.
    pub passive_coordinate: u8,
    /// `Interpolation Cylindric`: search at `(sqrt(x^2+y^2), z)`.
    pub cylindric: bool,
    /// `Interpolation Partial Hit`: accept elements with some inactive nodes.
    pub partial_hit: bool,
    /// `Try Linear Search If Qtree Fails` (default true).
    pub try_linear_search: bool,
    /// `UseQuadrantTree` argument (default true).
    pub use_quadrant_tree: bool,
    /// `Interpolation Abort Not Found`: make unfound nodes an error.
    pub abort_not_found: bool,
    /// `CoordinateSystemDimension()` of the model. The octree and
    /// `GlobalToLocal` use `max(this, mesh.dim())` and this respectively;
    /// `None` means "take it from the old mesh".
    pub coordinate_system_dim: Option<usize>,
    /// `Use Piola Transform` with `Edge Basis`: remap local coordinates.
    pub piola_transform: bool,
    /// `KeepUnfoundNodes` of the projector build (default true).
    pub keep_unfound_nodes: bool,
}

impl Default for InterpOptions {
    fn default() -> Self {
        InterpOptions {
            global_epsilon: 2.0e-10,
            local_epsilon: 1.0e-10,
            max_iterations: 12,
            numeric_epsilon: 1.0e-10,
            passive_coordinate: 0,
            cylindric: false,
            partial_hit: false,
            try_linear_search: true,
            use_quadrant_tree: true,
            abort_not_found: false,
            coordinate_system_dim: None,
            piola_transform: false,
            keep_unfound_nodes: true,
        }
    }
}

impl InterpOptions {
    /// Read the keywords from a `Simulation` section, with `Interpolation
    /// Passive Coordinate` and `Interpolation Cylindric` falling back to the
    /// `Solver` section as in the Fortran.
    pub fn from_params(simulation: &dyn Params, solver: Option<&dyn Params>) -> Self {
        let mut o = InterpOptions::default();
        if let Some(v) = simulation.real("interpolation global epsilon") {
            o.global_epsilon = v;
        }
        if let Some(v) = simulation.real("interpolation local epsilon") {
            o.local_epsilon = v;
        }
        if let Some(v) = simulation.integer("interpolation max iterations") {
            o.max_iterations = v.max(0) as usize;
        }
        if let Some(v) = simulation.real("interpolation numeric epsilon") {
            o.numeric_epsilon = v;
        }
        let pc = simulation
            .integer("interpolation passive coordinate")
            .or_else(|| solver.and_then(|s| s.integer("interpolation passive coordinate")));
        o.passive_coordinate = pc.unwrap_or(0).clamp(0, 3) as u8;
        o.cylindric = simulation
            .logical("interpolation cylindric")
            .or_else(|| solver.and_then(|s| s.logical("interpolation cylindric")))
            .unwrap_or(false);
        o.partial_hit = simulation
            .logical("interpolation partial hit")
            .unwrap_or(false);
        o.try_linear_search = simulation
            .logical("try linear search if qtree fails")
            .unwrap_or(true);
        o.abort_not_found = simulation
            .logical("interpolation abort not found")
            .unwrap_or(false);
        if let Some(s) = solver {
            if s.logical("edge basis").unwrap_or(false) {
                o.piola_transform = s.logical("use piola transform").unwrap_or(false);
            }
        }
        o
    }
}
