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
  bulk and boundary kernels, `default_first_order_time` (transient, with the
  solver's `elmer-timeint` scheme), `default_dirichlet` and `default_solve`,
  stopping when `default_converged` reports the nonlinear change below
  tolerance. Scalar properties and boundary values are evaluated at the
  element nodes with `elmer-expr` (`SolverCtx::real_at_nodes`) at the current
  iterate each iteration, so tables, MATC and procedures of `Temperature`
  give a Picard iteration; tensor conductivity must be constant. The element
  basis and quadrature come from `install_fe(provider)` because registry
  factories are plain `fn` pointers; until `elmer-fe` calls it, `factory`
  returns `NotWired`.

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
without `Stabilize = True`), `Legacy Assembly`, adaptive refinement,
`Nonlinear System Linesearch`, and `<property> At IP` evaluation.

## Not yet

- Convection through `HeatSolver` (blocked on `elmer-fe`): the kernels take
  SUPG `hK` and `StabilizationMK`, but `HeatSolver` cannot get them. Elmer
  computes both per element at mesh load (`MeshStabParams` ->
  `StabParam`, Stabilize.F90). Needed from `elmer-fe`, for example
  `fn stab_param(etype: ElementType, xyz: &[[f64; 3]]) -> (f64 /* mK */, f64 /* hK */)`
  passed through `install_fe` with the basis, or `hK`/`mK` on `MeshView`'s
  element. Until then the solver refuses `Convection` at setup.
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
  through `Simulation`; solutions against the Fortran run.
- `tests/heat_breadth.rs`: S2 features against full-precision results of the
  pristine ElmerSolver: `heateq_bdf2` (BDF2, variable step, nodal max diff
  below 1e-13), `heateq_spline_table` (cubic `k(T)` table, Picard converging
  in Elmer's 10 iterations, norm to 1e-10), and Robin and axisymmetric
  variants of ThermalPipe2D (norms to 1e-13).
  `Simulation` runs with its default `elmer-assembly` and `elmer-timeint`
  backends, `elmer-krylov` BiCGStab with `elmer-sparse` ILU(0), and a dense LU
  standing in for UMFPACK (`tests/support/cases.rs`).
- `tests/kernels.rs`: closed-form local matrices.
- `tests/mms.rs`: L2 convergence rates for P1/P2 triangles and quads.
- `tests/preflight.rs`: refusals.
