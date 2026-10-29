# elmer-physics-scalar

Three scalar physics solvers for the Rust-native Elmer rewrite, ported onto
`elmer-solver`'s `Solver` trait. They are the electrostatics and acoustics
physics that ElmerGUI offers (decision D2), plus the static current
conduction solver.

| Procedure | Reference | Quantity |
| --- | --- | --- |
| `"StatElecSolve" "StatElecSolver"` (also `"StatElecSolveVec"`) | `fem/src/modules/StatElecSolve.F90` | electric potential |
| `"StatCurrentSolve" "StatCurrentSolver"` | `fem/src/modules/StatCurrentSolve.F90` | electric potential in a conductor |
| `"HelmholtzSolve" "HelmholtzSolver"` | `fem/src/modules/HelmholtzSolve.F90` | complex pressure (two dofs) |

The element kernels are pure functions of geometry and nodal values
(`stat_elec_bulk_local`, `stat_elec_boundary_local`, `stat_current_bulk_local`,
`stat_current_boundary_local`, `potential_post_local`, `helmholtz_bulk_local`,
`helmholtz_boundary_local`). `register(&mut SolverRegistry)` binds the three
solvers; the element basis and quadrature come from `install_fe(provider)`
because registry factories are plain `fn` pointers, and until `elmer-fe` calls
it the factories return `NotWired`.

## Implemented

StatElecSolve (the scalar-element `LocalMatrix` path, which serial runs take):
`Relative Permittivity` (scalar, or field dependent through `elmer-expr`),
`Charge Density`, `Permittivity Of Vacuum` (default 8.854187817e-12),
`Electric Flux` / `Surface Charge Density`, the electrostatic layer condition
(`Layer Relative Permittivity`, `Layer Thickness`, `Electrode Potential`,
`Layer Charge Density`), axisymmetric coordinates, the nonlinear loop. The
post solver runs at the end of the solve: `Electric Energy Density`,
`Electric Flux`, `Electric Field`, `Nodal Energy Density`, `Nodal Electric
Flux`, `Constant Weights`, `res: Electric Energy`, `res: Potential Difference`,
`RES: Effective Capacitance`, `Energy Control` / `Charge Control`, and
`Calculate Average Potential`.

StatCurrentSolve: `Electric Conductivity` as scalar, diagonal or full tensor,
`Current Source`, `Current Density`, the Robin condition (`External
Conductivity`, `External Potential`), the transient displacement-current mass
term (`Relative Permittivity`) through the solver's time scheme, axisymmetric
coordinates, and the post fields `Joule Heating`, `Volume Current`, `Electric
Field`, `Nodal Joule Heating`, `Nodal Current`, with `RES: Total Joule
Heating`, `RES: Effective Resistance` and `Current Control` / `Power Control`.

HelmholtzSolve: `Sound Speed`, `Density`, `Sound Damping`, `Pressure Source`,
`Wave Flux`, `Wave Impedance`, `Plane Wave BC`, `Pressure 1` / `Pressure 2`
Dirichlet values, axisymmetric coordinates, `Frequency` / `Angular Frequency`
looked up as `GetAngularFrequency` does, and `res: Frequency`.

## Refused, not approximated

Each of these returns `SolverError::Unsupported` from `init`:

- StatElecSolve: `Legacy Assembly`, complex permittivity (`Relative
  Permittivity Im`, `Variable DOFs = 2`), `Electric Infinity BC` and
  `Infinity BC`, `Calculate Capacitance Matrix` and the `Capacitance Body`
  constraint modes, `Piezo Material`, elemental (`-dg`) fields, `Average
  Within Materials`, `Adaptive Mesh Refinement`, `Relative Integration Order`.
- StatCurrentSolve: `Legacy Assembly`, `Harmonic Mode`, `Farfield Potential`,
  `Calculate Conductivity Matrix` and `Conductivity Body`, the same DG,
  refinement and integration-order keywords.
- HelmholtzSolve: `Bubbles`, mean-flow `Convection Velocity`, `Flow Interface`
  and `Structure Interface` conditions. `Linear System Preconditioning Damp
  Coefficient` is accepted with a warning: the damped preconditioning matrix is
  not built, which changes the iteration but not the solution.

## What the solvers need from `elmer-solver`

- **Boundary normals.** `Electric Infinity BC`, `Farfield Potential` and the
  Helmholtz convection terms call `NormalVector(Element, Nodes, u, v,
  .TRUE.)`. Needed: `ElementBasis::normal(&self, etype: ElementType,
  nodal_xyz: &[[f64; 3]], uvw: [f64; 3], parent_centre: Option<[f64; 3]>,
  out: &mut [f64; 3]) -> Result<(), String>`, the unit normal at a boundary
  point, flipped to point away from the parent element's centre.
- **`Linear System Complex`.** `SolverCtx::default_solve` refuses it. The
  Helmholtz solver assembles the interleaved real system `[[Re, -Im], [Im,
  Re]]` per entry and does not set the flag, so the same linear system is
  solved by the real direct and Krylov methods. What is missing is Elmer's
  complex scaling and the complex Krylov methods and preconditioners
  (`Linear System Preconditioning Damp Coefficient`), plus a complex
  `Relative Permittivity Im` path for StatElecSolve.
- **`Initialize Dirichlet Conditions`.** `SetInitialConditions`
  (ElmerSolver.F90:1862) writes the Dirichlet values into every solver
  variable before the first time step. `Simulation::setup` does not, so each
  solver here does it in `init` for its own variable (`common::init_dirichlet`);
  without it a transient run starts with zeros on the Dirichlet nodes.
- **Steady state detection.** `Simulation::run` reports the steady state as
  reached after the first of two steady-state iterations for
  `fem/tests/current`, where ElmerSolver takes two (`RunReport::steps == 1`).
  The test states `Steady State Min Iterations = 2` to get the pristine
  result.

## Tests

`tests/data/README.md` records the pristine ElmerSolver runs the end-to-end
tests compare against. They use thin stand-ins under `tests/support` for
`elmer-fe`, `elmer-mesh` and `elmer-sif` until the S1 lane merges, and run the
default `elmer-solver` backends (assembly, time integration, direct and Krylov
solvers).
