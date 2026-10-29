# Fortran reference dumps

`linsys_a.dat`, `linsys_b.dat`, `linsys_perm.dat` and `linsys_sol.dat` were written by the
pristine ElmerSolver (`build/fem/src/ElmerSolver`, Version 26.2-devel, Rev
9f6af2f85) running `case.sif` with `Linear System Save = True`, `Linear System Save Solution` and `Linear
System Save Perm = True`, serially on copies of `fem/tests/MinimalExample` and
`fem/tests/ThermalPipe2D` (the latter after `ElmerGrid 14 2 OnePipe.msh
-autoclean`). The dump is taken in `DefaultSolve`, so Dirichlet rows are
already applied (diagonal kept, rest of the row zeroed, `b = diag * value`).
Only `case.sif`'s `Linear System Save*` lines differ from the upstream test.

`linsys_sol.dat` comes from a second run of the same case with `Linear System
Save Slot = "finish"`, so it holds the solution (the default slot dumps before
the solve). The dump prints it with 6 significant digits (`ES15.6`), so
comparisons against it use a 1e-5 tolerance. ThermalPipe2D's solution is
BiCGStab to 1e-10 relative.

## Full-precision solutions (`tests/heat_breadth.rs`)

`heateq_bdf2/`, `heateq_spline_table/`, `thermal_pipe2d_robin/` and
`thermal_pipe2d_axi/` come from the same pristine binary, run serially on
copies of the upstream tests with one added Simulation line, `Output File =
"case.result"`, which writes an ASCII `.result` with 17 significant digits.
`solution.dat` is the last time level as `node value` lines (SaveResult
writes `Values(Perm(i))` in node order), and the norms in the tests are
`sqrt(sum(x^2)/n)` over those values. Meshes come from `ElmerGrid 1 2 angle`
as each test's `runtest.cmake` does; `cond.dat` is copied from the upstream
test.

- `heateq_bdf2/`: `fem/tests/heateq_bdf2` unchanged otherwise. Elmer's norm
  0.05214225320136147 (the sif's Reference Norm is 5.214225320136149E-002).
- `heateq_spline_table/`: `fem/tests/heateq_spline_table` unchanged otherwise.
  Elmer converges in 10 nonlinear iterations to norm 1.513099967071337. The
  solution file is not stored (2821 nodes); the test checks the norm and the
  iteration count.
- `thermal_pipe2d_robin/`: ThermalPipe2D with UMFPACK instead of BiCGStab and
  the pipe wall (boundaries 1-4) changed from `Temperature = 0` to `Heat
  Transfer Coefficient = 20.0`, `External Temperature = 2.0`. Norm
  6.6452169458174035. Mesh: `thermal_pipe2d/`.
- `thermal_pipe2d_axi/`: the Robin variant with `Coordinate System = Axi
  Symmetric`, every node moved by x + 6 (so r is in [1, 11]), and `Heat Source
  = 0.01` on body 1. Norm 35.21274966573315. Mesh: `thermal_pipe2d/` moved in
  the test.

Each directory keeps the sif that was run.

## Convection (`angle_convection/`)

`case.sif` is a new steady case on the `heateq_bdf2` `angle` mesh, which the
test reads from `heateq_bdf2/`. It sets `Convection = Constant` with velocity
(2, 1), `Heat Conductivity = 0.01`, `Density = Heat Capacity = 1`, `Heat
Source = 1`, `Temp = 0` on boundary 1, and `Stabilize = True`. The pristine
ElmerSolver ran it serially with `Output File = "case.result"`.
`solution.dat` holds the nodal temperatures (norm 0.3278430808661644; the log
prints 0.32784308). A gdb breakpoint in `LocalMatrixVec` confirmed that the
SUPG block reads `Element % StabilizationMK = 1/3` and `Element % hK` equal to
`ElementDiameter`, as `elmer_fe::stab_param` returns them.
