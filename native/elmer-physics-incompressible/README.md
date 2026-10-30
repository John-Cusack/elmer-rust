# elmer-physics-incompressible

`IncompressibleNSSolver` for the Rust-native Elmer: incompressible (and
thermally or perfect-gas compressible) Navier-Stokes on equal-order nodal
elements held stable by `Pressure Stabilization` (a Brezzi-Pitkaranta
pressure block, with SUPG/PSPG terms for convection and the time
derivative), on the MINI element (velocity bubbles condensed by
`LCondensate`) and on Taylor-Hood p elements. Registered as `Procedure =
"IncompressibleNS" "IncompressibleNSSolver"` and the front `"IncompressibleNSVec"
"IncompressibleNSSolver"`.

Ported from `fem/src/modules/IncompressibleNS.F90` (`LocalBulkMatrix`,
`LCondensate`, `LocalBoundaryMatrix`, `IncompressibleNSSolver`, `_Init0`,
`_init`); `IncompressibleNSVec.F90` only delegates to it. `hK` comes from
`elmer-fe`'s `stab_param`, and so do the default basis and Gauss rules on
nodal elements; p elements take their dofs, basis, rule and node order from
the solver's `SolverDofs` (`elmer-pdofs`). The nodal boundary integrals are
`elmer-physics-flow`'s `flow_boundary_local` (the same `NavierStokesBoundary`
element integrals) with no mass flux.

## Supported

- `Pressure Stabilization = True` (or `Stabilize = True` / `Stabilization
  Method = Stabilized`, its synonyms), which makes the element `Element =
  "n:1"` (`"p:1"` too on `SolverDofs`); `Pressure Stabilization Coefficient`
  and `Pressure Stabilization Element Size = shortest edge | longest edge |
  volume`; the per-family constants of the Fortran.
- Without it, p elements from `SolverDofs`: `_Init0`'s MINI element (`p:1
  -tri b:1 -quad b:3 ...`, its velocity bubbles condensed, with their history
  in a transient run: `LCondensate`'s BDF1) and Taylor-Hood (`p:2`, the edge
  and face pressures pinned to the nodal one; full elements' bubbles
  condensed). `Element Integration Points` and `Number of Integration
  Points` pick the bulk rule.
- `Stokes Flow` (also `Flow Model = Stokes`), steady and transient Navier-Stokes
  with Picard or Newton linearization (`ComputeChange` decides, through
  `Nonlinear System Newton After Iterations / Tolerance`), `Linear Equation`,
  `GradP Discretization`, `Div-Curl Discretization`, `Specific Load`,
  `Nonlinear System Min/Max Iterations`.
- 2D and 3D Cartesian and `Axi Symmetric` (no swirl: the r measure, the hoop
  stress and the `-p v_r` pressure term), element-constant `Density`,
  `Viscosity` and `Flow Bodyforce i` at the nodes (expressions through
  `elmer-expr`), `Noslip Wall BC`, Dirichlet on `Velocity i` / `Pressure`
  (zero on the edge and face dofs of p elements).
- `Boussinesq = True` in the body force: `-rho g beta (T - Tref)` on the
  momentum load, `g` the Constants' `Gravity`.
- `Compressibility Model = Thermal` and `Perfect Gas`, steady and transient,
  from the `Temperature` field: rho(T) or rho(p, T) at the integration
  points, the `u.grad(rho)` continuity term and, transient, the d(rho)/dt
  load (and a perfect gas's pressure mass).
- Boundaries, nodal and p, Cartesian and axisymmetric: `Surface Traction i`,
  `Normal Surface Traction` / `External Pressure` along the outward normal,
  `Slip Coefficient i`.
- Transient runs: the mass matrix (Galerkin plus the PSPG/SUPG time terms)
  goes through `Default1stOrderTime`; `Stokes Flow` has none, as in Fortran.
- Elmer's 1e-6 start value for an `Equation = Navier-Stokes` variable.
- `Viscosity Model = power law` (`Critical Shear Rate`), `carreau`, `cross`,
  `powell eyring`, `glen` (`Set Arrhenius Factor`, or the Arrhenius law of a
  `Relative Temperature` with `Glen Enhancement Factor`), as in
  `EffectiveViscosityVec`: the viscosity at the points from the previous
  iterate's shear rate, and with Newton its derivative in the Jacobian.
- `Save Viscosity`: the exported `Viscosity` (nodal average weighted by the
  hidden `Viscosity Weight`) and `Shearrate` (the nodal integral of the
  squared shear rate, which the Fortran leaves unnormalized) of the last
  assembly.
- `Save Slip`: the exported `Slip Coefficient` (the largest component's,
  averaged with the hidden `Slip Weight`) and `Slip Speed` (zero: friction
  is not ported) on the slip boundaries.

Refused (`SolverError::Unsupported`), not approximated: `Bubbles in Global
System`,
compressibility with `Pressure Stabilization` or in a transient Stokes run,
cylindric coordinates, `Normal Slip Coefficient`, friction / `Wall Law` /
FSSA / normal-tangential boundaries, `Pressure 1..3` (the Fortran itself asks
for `Surface Traction`), `Implicit Friction
Coefficient`, the block preconditioner and Schur operators, `Nonlinear System
Linesearch`, relative and adaptive integration orders, density that varies
within an element, the `power law too` and turbulence viscosity models.

## Wiring (for the integrator)

- `solver_defaults_with(&mut ParamMap, dim, transient, serendipity)` (or
  `solver_defaults`, serendipity elements) applies the keywords of `_Init0`
  and `_init` (`ListAddNew` semantics): the `Element = "n:1"` string with
  stabilization, else the MINI element string and its `Element Integration
  Points`,
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

- `tests/fortran_reference.rs`: six nodal cases against the pristine
  ElmerSolver (see `tests/data/README.md`), three of them the
  `StokesStabilized` quick test (tri, quad, tet). Every nodal value agrees to
  1e-11 or better (the test bound is 1e-9, 1e-8 for the transient case).
- `tests/taylorhood.rs`: the p element, axisymmetric and compressible cases
  (Taylor-Hood on quads, triangles and tetrahedra, full elements, transient,
  with boundary loads; MINI Stokes, Navier-Stokes and transient;
  `IncompressibleNSAxi`; the thermal and perfect gas cavities), every slot
  against the pristine `case.result` to 1e-10.
- `tests/compressible_transient.rs`: the transient compressible cavities, the
  temperature solved by ModelPDE (`elmer-physics-advdiff`).
- `tests/solver_options.rs`: refusals and the keyword defaults.

`tests/support/` stands in for `elmer-mesh` (mesh files with boundary parents)
and `elmer-sif` (hand-built sections).

## Not done

Turbulence models,
friction and wall-law boundaries, normal-tangential conditions, and the block
preconditioner (its Schur matrix needs a child solver `elmer-solver` cannot
create yet).
