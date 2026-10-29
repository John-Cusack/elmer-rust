# Pristine ElmerSolver references

Each directory holds the `case.sif` run by the pristine ElmerSolver
(`build/fem/src/ElmerSolver`, Version 26.2-devel, Rev 9f6af2f85), serially, on
a mesh written by `ElmerGrid 1 2 <grd>` (Rev 9f6af2f85) from the `.grd` next
to it, and the resulting reference:

- `mesh.*`: the ElmerGrid output the run used.
- `solution.dat`: the last time level of the ASCII `.result` that
  `Output File = "case.result"` writes (17 significant digits), as
  `node u v [w] p` lines. SaveResult writes `Values(Perm(i))` in node order.
- The norm in `tests/fortran_reference.rs` is the last `FlowSolve: Result
  Norm` line of the run's stdout. FlowSolve sets `Nonlinear System Norm DOFs =
  dim`, so it is the RMS of the velocity components only.

The cases reuse the geometry, boundary numbering (1 bottom, 2 right, 3 top,
4 left) and wall conditions of the quick test `fem/tests/ChannelFlowCoulomb`
(`rect.grd`, a 4 x 1 channel here split 12 x 4). That test itself needs
implicit friction and load computation, which this crate does not port yet.
The 3D case uses a 2 x 1 x 1 box of 4 x 2 x 2 bricks, with tags 1/3 at y = 0/1,
2/4 at x = 2/0 and 5/6 at z = 0/1.

| case | elements | what it exercises |
|---|---|---|
| `stokes_channel` | 404 | `Flow Model = Stokes`, `Flow BodyForce 1`, `External Pressure` outlet |
| `ns_channel_tri` | 303 | Picard then Newton after 3 iterations, `Slip Coefficient 1` walls, 8 iterations to 1e-10 |
| `ns_quad_newton` | 404 | Newton from iteration 1, `Div Discretization`, `Noslip Wall BC`, `Hydrostatic Pressure` + `Flow BodyForce 2`, `Pressure 1` traction |
| `nc_tri_gradp` | 303 | `Flow Model = No Convection`, `Stabilization Method = Stabilized`, `Gradp` and `Laplace Discretization`, `Normal Pressure`, `Mass Flux`, a `Pressure` Dirichlet |
| `transient_tri` | 303 | BDF2, 4 steps of 0.1 from rest, Newton after 2 iterations each step |
| `ns_hex3d` | 808 / 404 faces | 3D, slip walls, `External Pressure` along the face normal |

All runs use `Linear System Solver = Direct` / `UMFPack`; the Rust side solves
with the default faer LU.

Two Elmer defaults only show up in some of these cases:

- `Equation = Navier-Stokes` creates the flow variable filled with 1e-6
  (`AddEquationBasics`, MainUtils.F90:1241). The steady results are
  independent of the start value. The transient one is not: the 1e-6 enters
  the first history, and leaving it out costs 1e-6 relative in the solution.
- `Element % StabilizationMK` is filled only when some solver says
  `Stabilize = True` or `Stabilization Method = Stabilized|VMS`
  (`MeshStabParams`, MeshBasics.F90:2338). With FlowSolve's implicit
  default, mK = 0 and the run is unstabilized (singular at an inflow corner
  in `nc_tri_gradp`'s geometry), so every case names one of the keywords.
