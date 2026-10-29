# elmer-physics-flow

FlowSolve for the Rust-native Elmer: incompressible Navier-Stokes with
equal-order velocity and pressure and FlowSolve's default residual-based
stabilization (`Stabilization Method = Stabilized`, Franca et al.).
`Procedure = "FlowSolve" "FlowSolver"`.

Ported from `fem/src/modules/FlowSolve.F90` (`FlowSolver`, `FlowSolver_init`)
and `fem/src/NavierStokes.F90` (`NavierStokesCompose`, cartesian
incompressible stabilized path; `NavierStokesBoundary`, cartesian). `hK` and
`mK` come from `elmer-fe`'s `stab_param` (`StabParam`/`ElementDiameter`), and
so do the default basis (`LagrangeBasis`) and Gauss rules (`ElmerQuadrature`).

## Supported

- `Flow Model` = full (default), `No Convection`, `Stokes`; cartesian 2D and
  3D (`Coordinate System = Cartesian*`).
- Picard iteration with FlowSolve's Newton switching: `Nonlinear System
  Newton After Iterations / Tolerance`, `... Newton Max Iterations / Max
  Tolerance` (back to Picard), `... Reset Newton`. The Newton state persists
  across calls, as FlowSolve's implicitly SAVEd locals do. Also `Nonlinear
  System Max Norm Return`.
- `Div`, `Gradp` and `Laplace Discretization`, `H Scale`, `Stabilization Use
  Longest Element Edge`.
- Density and viscosity at nodes through `elmer-expr` (constants, tables,
  MATC, procedures), `Flow BodyForce i`, `Hydrostatic Pressure` (body force or
  equation) with the Constants `Gravity`.
- Boundaries: `Pressure i` traction, `External Pressure` / `Normal Pressure`
  along the outward normal, `Mass Flux`, `Slip Coefficient i`, `Flow Force BC
  = False`, Dirichlet on `Velocity i` / `Pressure` (component names), `Noslip
  Wall BC`.
- Transient runs: the stabilized mass matrix goes through `Default1stOrderTime`
  (`elmer-timeint`, all its first-order schemes).
- FlowSolve's keyword side effects: `Nonlinear System Norm DOFs = dim`,
  `Relative Pressure Relaxation` when every material is incompressible,
  `Time Derivative Order = 1`. Also Elmer's 1e-6 start value for an
  `Equation = Navier-Stokes` variable, and mK = 0 (no stabilization) unless
  some solver names `Stabilize`/`Stabilization Method` (`Mesh % Stabilize`).
  An absent `Nonlinear System Max Iterations` means zero iterations, as in
  FlowSolve.

Refused (`SolverError::Unsupported`), not approximated: bubbles, P2/P1, VMS,
`Element` p-bubbles, compressibility models other than incompressible,
`Viscosity Model` (non-newtonian, anisotropic), porous media, Boussinesq,
Lorentz and potential forces, rotating frames, `Mesh Velocity`, `Solver
Coordinate System`, non-cartesian coordinates, out-of-plane flow,
normal-tangential, free-surface, surface-tension, wall-law, VMS-wall,
periodic, implicit-friction and point (`Target Coordinates`/`Nodes`)
boundaries, `Calculate Loads`, adaptive refinement, line search.

## Tests

- `tests/fortran_reference.rs`: six cases against the pristine ElmerSolver
  (see `tests/data/README.md`). Final norm and every nodal value agree to
  1e-14 or better (the test bound is 1e-9).
- `tests/kernels.rs`: normal orientation, the pressure block equals `Tau
  (grad q, grad p)` (and vanishes with mK = 0), Newton and Picard residuals
  coincide at the linearization point.
- `tests/solver_options.rs`: refusals, and FlowSolve's zero-iteration default.

`tests/support/` still stands in for `elmer-mesh` (mesh files with boundary
parents) and `elmer-sif` (hand-built sections).

## Gaps in the shared crates

These are handled locally; each would be better in its owner's crate.

- `elmer-solver`: `SolverCtx::dirichlet_values` keys components as `"<name>
  k"`, while DefaultDirichletBCs uses `ComponentName` (`velocity 1`,
  `pressure`). FlowSolver builds its own rows from `parse_variable_spec(..)
  .components` and calls `AssemblyBackend::try_apply_dirichlet`. Wanted:
  `dirichlet_values` keyed by the variable's component names.
- `elmer-solver`: `AddEquationBasics`' `InitValue = 1e-6` for `Equation =
  Navier-Stokes` (MainUtils.F90:1241). FlowSolver sets it in `init` for the
  entries no `Initial Condition` sets.
- `elmer-solver`: no `SolverCtx::params_mut(&mut self) -> &mut ParamMap`.
  FlowSolver reaches its slot through `ctx.sim.slots` by solver id to add
  the keywords FlowSolve adds with `ListAddInteger`/`ListAddNewLogical`.
- `elmer-core` `MeshView`: no boundary `OutBody`. Normal orientation follows
  `CheckNormalDirection` for one parent and for two parents of different
  bodies. With two parents of the same body the normal is left unflipped.
