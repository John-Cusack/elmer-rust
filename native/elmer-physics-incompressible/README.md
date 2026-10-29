# elmer-physics-incompressible

`IncompressibleNSSolver` for the Rust-native Elmer: incompressible
Navier-Stokes with equal-order nodal velocity and pressure held stable by
`Pressure Stabilization` (a Brezzi-Pitkaranta pressure block, with SUPG/PSPG
terms for convection and the time derivative). Registered as `Procedure =
"IncompressibleNS" "IncompressibleNSSolver"` and the front `"IncompressibleNSVec"
"IncompressibleNSSolver"`.

Ported from `fem/src/modules/IncompressibleNS.F90` (`LocalBulkMatrix`,
`LocalBoundaryMatrix`, `IncompressibleNSSolver`, `_Init0`, `_init`);
`IncompressibleNSVec.F90` only delegates to it. `hK` comes from `elmer-fe`'s
`stab_param`, and so do the default basis and Gauss rules. The boundary
integrals are `elmer-physics-flow`'s `flow_boundary_local` (the same
`NavierStokesBoundary` element integrals) with no mass flux.

## Supported

- `Pressure Stabilization = True` (or `Stabilize = True` / `Stabilization
  Method = Stabilized`, its synonyms), which makes the element `Element =
  "n:1"`; `Pressure Stabilization Coefficient` and `Pressure Stabilization
  Element Size = shortest edge | longest edge | volume`; the per-family
  constants of the Fortran.
- `Stokes Flow` (also `Flow Model = Stokes`), steady and transient Navier-Stokes
  with Picard or Newton linearization (`ComputeChange` decides, through
  `Nonlinear System Newton After Iterations / Tolerance`), `Linear Equation`,
  `GradP Discretization`, `Div-Curl Discretization`, `Specific Load`,
  `Nonlinear System Min/Max Iterations`.
- 2D and 3D Cartesian, element-constant `Density`, `Viscosity` and `Flow
  Bodyforce i` at the nodes (expressions through `elmer-expr`), `Noslip Wall
  BC`, Dirichlet on `Velocity i` / `Pressure`.
- Boundaries: `Surface Traction i`, `Normal Surface Traction` / `External
  Pressure` along the outward normal, `Slip Coefficient i`.
- Transient runs: the mass matrix (Galerkin plus the PSPG/SUPG time terms)
  goes through `Default1stOrderTime`; `Stokes Flow` has none, as in Fortran.
- Elmer's 1e-6 start value for an `Equation = Navier-Stokes` variable.

Refused (`SolverError::Unsupported`), not approximated: the default MINI
element (`p:1` with condensed velocity bubbles, `Element = "p:1 b:..."`) and any
`Element` but `n:1`, Taylor-Hood, `Bubbles in Global System`, `Compressibility
Model` (thermal, perfect gas), `Viscosity Model` and turbulence viscosity,
`Boussinesq`, axisymmetric and other non-Cartesian coordinates, `Normal Slip
Coefficient`, friction / `Wall Law` / FSSA / normal-tangential boundaries,
`Pressure 1..3` (the Fortran itself asks for `Surface Traction`), `Save
Viscosity`, `Save Slip`, `Implicit Friction Coefficient`, the block
preconditioner and Schur operators, `Nonlinear System Linesearch`, density that
varies within an element.

## Wiring (for the integrator)

- `solver_defaults(&mut ParamMap, dim, transient)` applies the keywords of
  `_Init0` and `_init` (`ListAddNew` semantics): the `Element = "n:1"` string,
  `Variable = "Flow Solution[Velocity:dim Pressure:1]"`, `Nonlinear System Norm
  DOFs = dim`, `Relative Pressure Relaxation`, `Stokes Flow` from `Flow Model`.
  `Variable` and `Element` have to be in the section before the mesh and the
  solver variable are created, so call it on each solver section that names this
  procedure before `Simulation::setup`; `init` repeats it for what still takes
  effect.
- The implicit legacy solver `Equation = "Navier-Stokes"` with no `Procedure`
  maps (`AddEquationBasics`, MainUtils.F90:1231) onto **FlowSolve**, not this
  crate: `Procedure = "FlowSolve" "FlowSolver"`, `Variable = "Flow
  Solution[Velocity:dim Pressure:1]"` (`Velocity:3` for `Cylindric Symmetric`
  as well, with `dofs = dim` + 1 there), the variable start value 1e-6, and
  `Variable Dofs` left unset. `elmer-physics-flow` already applies the 1e-6; the
  sif adapter has to add the `Procedure` and `Variable` keywords when they are
  absent, plus the same for `magnetic induction`, `stress analysis`, `mesh
  update` and `heat equation` (their own procedures).

## Tests

- `tests/fortran_reference.rs`: six cases against the pristine ElmerSolver (see
  `tests/data/README.md`), three of them the `StokesStabilized` quick test
  (tri, quad, tet). Every nodal value agrees to 1e-11 or better (the test
  bound is 1e-9, 1e-8 for the transient case).
- `tests/solver_options.rs`: refusals and the keyword defaults.

`tests/support/` stands in for `elmer-mesh` (mesh files with boundary parents)
and `elmer-sif` (hand-built sections).

## Not done

The default MINI element (`Element = "p:1 -tri b:1 ..."`) with `LCondensate`,
which most of the 58 sif files naming this solver use, needs the p-hierarchic
bubble basis and bubble history; Taylor-Hood needs P2 elements. Non-newtonian
viscosity (`EffectiveViscosityVec`), compressibility, axisymmetry, friction and
wall-law boundaries, and the block preconditioner are ported nowhere yet.
