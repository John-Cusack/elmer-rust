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
