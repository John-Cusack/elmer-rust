# elmer-physics-flux

`FluxSolver` for the Rust-native Elmer rewrite (lane `w3flux`, stage S3):
nodal recovery of fluxes and gradients of a scalar field by Galerkin
projection. Reference: `fem/src/modules/FluxSolver.F90`.

## What it computes

For a target field `u` (`Target Variable`, default `Temperature`) and a
material tensor `C` (`Flux Coefficient`, default `Heat Conductivity`), each
requested field is the nodal L2 projection of an integrand onto the solver's
mass-matrix system `M x = b`:

| keyword | exported variable | integrand |
|---|---|---|
| `Calculate Flux` | `<u> Flux[<u> Flux:dim]` | `-C grad u` |
| `Calculate Grad` | `<u> Grad[<u> Grad:dim]` | `grad u` |
| `Calculate Flux Magnitude` | `<u> Flux_mag` | `abs(C grad u)` pointwise |
| `Calculate Grad Magnitude` | `<u> Grad_mag` | `abs(grad u)` pointwise |
| `Calculate Flux Abs` / `Calculate Grad Abs` | `<u> Flux_abs` / `<u> Grad_abs` | norm of the projected components |

Every component reuses one matrix with its own right-hand side. The solver
variable is the hidden `-nooutput <Equation>_temp`, and the defaults of
`FluxSolver_Init` apply: CG with ILU0, `1e-10`, 500 iterations,
`Time derivative order = 0`. The reported norm is that of the last projected
component, and only that one gets a nonlinear-change evaluation (`Skip Compute
Nonlinear Change` is managed as in the Fortran).

## API

- `flux_bulk_local(basis, quad, etype, xyz, potential, coeff, cfg, out)`: the
  element mass matrix and one right-hand side per field. A pure function of
  geometry and nodal values.
- `FluxSolver` (an `elmer_solver::Solver`), `register(&mut SolverRegistry)` for
  `("FluxSolver", "FluxSolver")`. As in `elmer-physics-heat`, the element basis
  and quadrature arrive through `install_fe(provider)`; until `elmer-fe` calls
  it, `factory` returns `NotWired`.

The framework has no hook that runs before solver variables are created (the
role of `FluxSolver_Init`), so `FluxSolver::init` creates the hidden variable,
the exported variables and the time-order and linear-system defaults itself.

## Supported

Scalar, diagonal and full-tensor `Flux Coefficient` (constant arrays, or a
scalar-valued dependent property through `elmer-expr`); `Target Variable`,
`Flux Variable` and `Target Expression` (a material keyword); `Enforce Positive
Magnitude`; `Constant Bulk Matrix` (the matrix is simply reassembled, same
result); axisymmetric weighting; 2D and 3D.

## Refused

- `Discontinuous Galerkin` (the `fluxsolver3` quick test): needs an elemental
  DG solver variable and the face-term assembly. `elmer-solver` refuses
  non-nodal solver variables today.
- Coordinate systems other than Cartesian, `Axi Symmetric` and `Cylindric
  Symmetric`.
- A target variable that is not a scalar nodal field.

## Tests

```
export CARGO_BUILD_JOBS=1 RUST_TEST_THREADS=1 OMP_NUM_THREADS=1 RAYON_NUM_THREADS=1
cargo test --manifest-path native/Cargo.toml -p elmer-physics-flux
```

- `tests/kernels.rs`: the element kernel against closed forms (exact
  projection of linear fields, tensors, interpolated coefficients, the
  axisymmetric weight).
- `tests/flux_solver.rs`: `Simulation` end to end on the `angle.grd` mesh (300
  `404` elements, 341 nodes) of the pristine quick tests, against fields and
  norms captured from the pristine ElmerSolver.

Test-local stand-ins under `tests/support/` (Lagrange basis, Gauss rules, mesh
reader and `MeshView`) are copies of `elmer-physics-heat`'s and are replaced by
`elmer-fe` and `elmer-mesh` when S1 lands. CG and ILU0 come from the real
`elmer-krylov` and `elmer-sparse` (dev-dependencies); assembly and time
integration are the `elmer-solver` defaults.

### Pristine references

Captured with the pristine `build/fem/src/ElmerSolver` (serial, `nice -n 10`,
`OMP_NUM_THREADS=1`), meshes from `ElmerGrid 1 2 angle`. The tests feed the
solver the pristine run's `Temperature` (VTU point data, full precision) so the
comparison isolates the flux recovery from HeatSolve.

| test | pristine case | `Result Norm` (last component) |
|---|---|---|
| `fluxsolver_quick_test_flux` | quick test `fluxsolver` (`flux.sif`, `Calculate Flux`) | `0.19881406507522117` (sif reference `0.198814065076`) |
| `fluxsolver2_quick_test_grad` | quick test `fluxsolver2` (`flux2.sif`, `Calculate Grad`) | `0.066271355025073725` (sif reference `0.66271355E-01`) |
| `all_fields_match_the_pristine_run` | `flux.sif` with all four requests, `Flux Abs`, `Grad Abs`, `Enforce Positive Magnitude` | `0.093238023871151321` |
| `axisymmetric_weighting_matches_the_pristine_run` | `flux.sif` with `Coordinate System = Axi Symmetric`, flux and grad | `0.072722197017180842` |
| `full_tensor_coefficient_matches_the_pristine_run` | `flux.sif` with `Heat Conductivity(2,2) = 3 0.5 0.5 1.5`, flux and flux magnitude | `0.29823103835624498` |

Fields (every exported variable, every node) and norms agree with the pristine
run to about `1e-15`; the tests assert `1e-12`.

`fluxsolver3` is not covered (Discontinuous Galerkin, see above).
