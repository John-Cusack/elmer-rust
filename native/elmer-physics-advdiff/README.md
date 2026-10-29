# elmer-physics-advdiff

The two generic transport solvers of the Rust-native Elmer rewrite:

- `AdvectionDiffusionSolver` (`Procedure = "AdvectionDiffusion"
  "AdvectionDiffusionSolver"`, `fem/src/modules/AdvectionDiffusion.F90`): scalar
  species transport with diffusion, constant or computed convection, Soret
  diffusion, Nernst-Planck migration and mesh velocity, in relative or
  absolute mass units.
- ModelPDE's `AdvDiffSolver` (`Procedure = "ModelPDE" "AdvDiffSolver"`,
  `fem/src/modules/ModelPDE.F90`): the prototype advection-diffusion-reaction
  equation.

## API

- Kernels, pure functions of element geometry and nodal values:
  `advdiff_bulk_local`, `advdiff_boundary_local` (`DiffuseConvectiveCompose`,
  `DiffuseConvectiveBoundary`), `model_pde_bulk_local`,
  `model_pde_boundary_local` (`LocalMatrix`, `LocalMatrixBC`). They return the
  mass matrix next to stiffness and force; time integration stays with the
  caller (`Default1stOrderTime`).
- `AdvDiffSolver`, `ModelPdeSolver` (`elmer_solver::Solver`) and
  `register(&mut SolverRegistry)`. Both build `elmer-fe`'s `LagrangeBasis` and
  `ElmerQuadrature` from the mesh dimension themselves, so there is no
  `install_fe`.
- `check_advdiff_supported`, `check_model_pde_supported`: preflight refusals.

## What follows the Fortran closely

- **Bubbles** (the AdvectionDiffusion default when neither `Stabilize` nor
  `Bubbles` is given): the element gets `n` extra functions `N_i * b` with
  `b` the product of the linear nodal functions (`EvalBubbleBasis`,
  ElemInfo.F90:571), integrated with the element type's second rule
  (`GaussPoints2`, e.g. 16 points on 404), then condensed with `CondensateP`.
  A non-linear element takes its `detJ` from the linear element of its family.
- **Time terms on bubble elements**: Elmer's local time schemes loop over the
  whole allocated local matrix, so the mass matrix enters the bubble rows and
  columns, with the bubbles' previous values zero. Condensation comes after.
  (`bubble_first_order_time` reproduces this through `elmer-timeint`.)
- **SUPG** (`Stabilize = True`): `hK` and `mK` from `elmer_fe::stab_param`.
  The residual needs second derivatives of the basis, which `elmer-fe` does
  not provide, so it works on linear simplices only (202, 303, 504, where they
  are zero). Other element types stop with an error naming the missing
  piece; see the blocker below.
- The zero-order term is always off (`NodalC0 = 0` in the Fortran), `C0`
  carries only the Soret coefficient.
- Boundary elements go through the time scheme with a zero mass matrix, as
  `Add1stOrderTime` is called on them.

## Refused

Non-Cartesian coordinate systems (`DiffuseConvectiveGenCompose`), compressible
species transport, `Concentration Units = Mass To Max Solubility` and the
solubility-jump boundary, FCT (`Use Global Mass Matrix`); for ModelPDE
p-elements (`Element = "p:2"`), CutFEM and the adaptive integration keywords of
`GaussPointsAdapt`.

## Open need from `elmer-fe`

SUPG on elements other than linear simplices needs `ddBasisddx`
(`ElementInfo(..., SecondDerivatives=.TRUE.)`, `GlobalSecondDerivatives`).
`BasisHessian` in `src/fe.rs` is the interim shape of that request:

```rust
fn eval_hessian(&self, etype: ElementType, xyz: &[[f64; 3]], uvw: [f64; 3],
                out: &mut Vec<[[f64; 3]; 3]>) -> Result<(), String>;
// out[p][i][j] = d2 N_p / dx_i dx_j for the xyz.len() nodal functions
```

Once `elmer-fe` implements it (or `BasisEval` grows the field), `AdvDiffSolver`
passes it in `Fe::hessian` in place of `SimplexHessian`.

## Tests

`tests/quick_tests.rs` runs hand-built models through `elmer-solver`'s
`Simulation` and compares with the pristine ElmerSolver's full-precision
solutions (see `tests/data/README.md`): `adv_diff1`, `adv_diff2`, `AdvDiffSoret`
and `AdvDiffSoretAbsMass` (as `fem/tests` ships them, modulo the one change
noted there), a linear-element variant of `ModelPDE`, and Galerkin, bubble and
SUPG variants of `adv_diff2`. `tests/kernels.rs` checks single elements by
hand. The mesh reader and dense direct solver in `tests/support/` stand in for
`elmer-mesh` and UMFPACK until those are wired.
