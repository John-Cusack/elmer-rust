# Pristine ElmerSolver references

Each directory holds the `case.sif` run by the pristine ElmerSolver
(`build/fem/src/ElmerSolver`), serially, on a mesh written by `ElmerGrid` and
the resulting reference:

- `mesh.*`: the ElmerGrid output the run used.
- `solution.dat`: the last time level of the ASCII `.result` that `Output File =
  "case.result"` writes (17 significant digits), as `node u v [w] p` lines in
  node order.
- The norm in `tests/fortran_reference.rs` is the last `ComputeChange` norm
  of the run's stdout (velocity components only: `Nonlinear System Norm DOFs =
  dim`).

| case | source | what it exercises |
|---|---|---|
| `stokes_tri` | `fem/tests/StokesStabilized/tri.sif` (`ElmerGrid 1 2 Step.grd -triangles`) | Stokes, 303, Brezzi-Pitkaranta block |
| `stokes_quad` | `StokesStabilized/quad.sif` (`ElmerGrid 1 2 Step.grd`) | Stokes, 404 |
| `stokes_tet` | `StokesStabilized/tet.sif` (`ElmerGrid 14 2 tet.msh`) | Stokes, 504, 3D |
| `ns_tri` | own | Navier-Stokes, Picard then Newton after 3, `Pressure Stabilization Element Size = longest edge` |
| `ns_quad` | own | Newton from iteration 1, `GradP` and `Div-Curl` forms, `volume` size, coefficient 1.5, `Flow Bodyforce 2`, `Normal Surface Traction`, `Surface Traction 2`, `Slip Coefficient 1` |
| `transient_tri` | own | BDF2, 4 steps of 0.1, Newton after 2 iterations: the PSPG/SUPG mass terms |

The three `StokesStabilized` sifs are the stock ones with the iterative solve
replaced by `Linear System Solver = Direct` / `UMFPack` (the Rust side solves
with the default faer LU) and `Output File = "case.result"` added. Their MATC
inflow profile has no MATC engine in the Rust test, so the test hands the
solver a linear table through the inflow nodes' coordinates, which reproduces
MATC exactly there. In `ns_quad` the pressure level is a null space of the
`GradP` form, so the pressure is compared up to a constant.
