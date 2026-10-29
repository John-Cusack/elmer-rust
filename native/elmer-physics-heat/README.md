# elmer-physics-heat

HeatSolve element kernels for the Rust-native Elmer rewrite (lane L13).
Reference: `fem/src/modules/HeatSolve.F90`, `LocalMatrixVec` (the default
vectorised path) and `LocalMatrixBC`.

## API

- `heat_bulk_local(basis, quad, etype, xyz, mat, cfg, out)`: stiffness, force
  and (transient) mass for one bulk element.
- `heat_boundary_local(basis, quad, etype, xyz, bc, cfg, out)`: flux and Robin
  terms for one boundary element.
- `check_supported(&HeatSections)`: preflight that refuses unsupported
  keywords instead of silently ignoring them.
- `PROCEDURE = ("HeatSolve", "HeatSolver")`.
- `HeatSolver` (an `elmer_solver::Solver`) and `register(&mut SolverRegistry)`.
  It loops `Nonlinear System Max Iterations` times over `default_initialize`,
  bulk and boundary kernels, `default_first_order_time` (transient),
  `default_dirichlet` and `default_solve`. Properties must be constants; a
  table, MATC or procedure value is refused. The element basis and quadrature
  come from `install_fe(provider)` because registry factories are plain `fn`
  pointers; until `elmer-fe` calls it, `factory` returns `NotWired`.

`cfg` is an addition to the lane spec's signature: dimension, axisymmetry,
transient and stabilisation are element switches the Fortran reads from
solver state. Properties arrive already evaluated at the element nodes;
the kernels interpolate with the element basis, as `ListGetElementReal` does.
Time integration is not applied: a transient element returns `mass` and the
time-integration lane forms `M/dt`.

## Supported

Scalar, diagonal and full-tensor `Heat Conductivity`; `Density`, `Heat
Capacity`; `Heat Source` (times density) and `Volumetric Heat Source`;
constant or computed convection velocity with Galerkin or SUPG weighting
(`hK` and `StabilizationMK` supplied by the caller); `Heat Transfer
Multiplier`; axisymmetric weighting; `Heat Flux`; `Heat Transfer Coefficient`
with `External Temperature`.

## Refused (see `preflight.rs`)

Radiation, `Farfield Temperature`, `Infinity BC`, phase change and latent-heat
substepping, compressibility, `Draw Velocity`, Joule heating, perfusion,
radiators, DG and `Heat Gap`, Nitsche conditions, bubbles (so convection
without `Stabilize = True`), `Legacy Assembly`, adaptive refinement.

## Not yet

- Convection through `HeatSolver`: the kernels take SUPG `hK` and
  `StabilizationMK`, but `MeshView` has no element sizes, so the solver
  refuses `Convection` at setup.
- Temperature-dependent properties (needs `elmer-expr` wiring and the Picard
  loop's property re-evaluation).
- Bitwise identity with Fortran: matrix parity here is 1e-13 relative through
  a test-local glue loop, basis and quadrature (`tests/support/`); bitwise
  needs `elmer-fe` and the oracle lane.

## Tests

```
cargo test --manifest-path native/Cargo.toml -p elmer-physics-heat
```

- `tests/fortran_reference.rs`: global matrix and rhs for MinimalExample and
  ThermalPipe2D against the pristine ElmerSolver's `Linear System Save`
  (`tests/data/README.md` records how the dumps were made).
- `tests/heat_solver.rs`: MinimalExample, ThermalPipe2D and a transient step
  through `Simulation` with test-local assembly, backward-Euler and dense LU
  backends; solutions against the Fortran run.
- `tests/kernels.rs`: closed-form local matrices.
- `tests/mms.rs`: L2 convergence rates for P1/P2 triangles and quads.
- `tests/preflight.rs`: refusals.
