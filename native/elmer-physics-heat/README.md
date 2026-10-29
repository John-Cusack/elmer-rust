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
convection in the kernels: constant or computed velocity, with Galerkin or
SUPG weighting (`hK` and `StabilizationMK` supplied by the caller).
`HeatSolver` runs `Convection = Constant` (the material's `Convection Velocity
i`, a missing component being zero) with SUPG. It takes `hK` and
`StabilizationMK` per element from `elmer_fe::stab_param`, which is what
`MeshStabParams` computes, honouring `Stabilization Use Longest Element
Edge`. `Heat Transfer
Multiplier`; axisymmetric weighting; `Heat Flux`; `Heat Transfer Coefficient`
with `External Temperature`.

## Refused (see `preflight.rs`)

Radiation, `Farfield Temperature`, `Infinity BC`, phase change and latent-heat
substepping, compressibility, `Draw Velocity`, Joule heating, perfusion,
radiators, DG and `Heat Gap`, Nitsche conditions, bubbles (so convection
without `Stabilize = True`), `Legacy Assembly`, adaptive refinement,
`Nonlinear System Linesearch`, and `<property> At IP` evaluation.

## Not yet

- `Convection = Computed` through `HeatSolver`: the kernels support it, but
  the solver does not yet read the flow solution as the convection field, so
  it refuses the keyword at setup. ALE (`Mesh Velocity`) subtraction is not
  applied either.
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
  variants of ThermalPipe2D (norms to 1e-13), and `angle_convection`
  (constant convection, conductivity 0.01, `Stabilize = True`: SUPG with
  `elmer-fe`'s `hK`/`mK`, nodal max diff below 1e-12, norm to 1e-13).
  `Simulation` runs with its default `elmer-assembly` and `elmer-timeint`
  backends, `elmer-krylov` BiCGStab with `elmer-sparse` ILU(0), a dense LU
  standing in for UMFPACK, and `elmer-fe`'s basis and Gauss rules
  (`tests/support/cases.rs`).
- `tests/kernels.rs`: closed-form local matrices.
- `tests/mms.rs`: L2 convergence rates for P1/P2 triangles and quads.
- `tests/preflight.rs`: refusals.
