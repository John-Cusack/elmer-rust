# Fortran reference eigenproblems

Stiffness and mass matrices dumped by the pristine Fortran `ElmerSolver`
(v26.2, revision 9f6af2f85), run serially on scratch copies of two quick tests
with `Eigen Analysis = True` and one small non-quick one. All three use
Elmer's default `Linear System Symmetric = False`, so the reference went
through ARPACK's `DNAUPD`/`DNEUPD`.

| Directory | Source test (`fem/tests/`) | Solver | n | nnz | Eigen System Values |
|---|---|---|---:|---:|---:|
| `plates_eigen` | `platesEigen` | Smitc (3 dofs, rotations massless) | 1023 | 16542 | 10 |
| `geomstiff` | `geomstiff` | StressSolve with geometric stiffness (2 dofs) | 1134 | 31775 | 10 |
| `wave_equ_eigen` | `WaveEquEigen` (not labelled quick) | WaveSolver on the unit square | 441 | 3329 | 10 |

All three ask for `Eigen System Select = smallest magnitude`, a consistent
mass and no shift, with `Eigen System Convergence Tolerance = 1.0e-6`.

## Runs

Each case was run twice. Each run added these lines to the end of `Solver 1`:

```
Linear System Save = Logical True
Linear System Save Slot = String "linear solve"
Linear System Save Mass = Logical True
Linear System Save Perm = Logical True
```

It also set `Filename` on the SaveScalars solver and `Output File` in
`Simulation`.

* **asis**: the test's own solver settings. platesEigen solves its
  shift-invert systems with BiCGStab + ILU0 to `1e-8`, WaveEquEigen with
  BiCGStabl + ILU2 to `1e-8` (and `Eigen System Max Iterations = 1`);
  geomstiff uses the default direct solver.
* **tight**: also `Linear System Solver = Direct`,
  `Linear System Direct Method = Umfpack`,
  `Eigen System Convergence Tolerance = 1.0e-12` and
  `Eigen System Max Iterations = 1000`.

The two runs of a case dump identical matrices (`cmp` of the two `_a.dat` files). The
inexact inner solves of the platesEigen asis run move its eigenvalues by up to
about 5e-10 relative to the tight run. The geomstiff runs agree in every
printed digit.

WaveEquEigen has double eigenvalues (the mesh is symmetric). The tight run
returns 130.286 and 173.064 once each and then 180.96 and 202.98, so it
misses one copy of each double. The as-shipped run finds both copies but
only to about 1e-9. The test checks that the native solver returns every copy.
Only its first, simple mode is stored in `modes.txt`.

## Files

* `system.txt`: the `_a.dat` dump (stiffness, mass columns) with entries whose
  stiffness and mass are both zero dropped. The first data line is `n nnz`,
  then one line per row: `row (col stiff mass)...`, 1-based, written with
  Python `repr` so every double round-trips exactly. The matrices are not
  symmetric: Dirichlet rows are reduced to their diagonal (zero mass row) but
  their columns stay, as Elmer leaves them for the eigensolver.
* `eigenvalues.txt`: the SaveScalars `Save EigenValues` columns of each run
  (`ES20.12`, 13 significant digits), in Elmer's output order.
* `modes.txt`: the first three eigenvectors (one for WaveEquEigen) of the
  tight run, from its `.result` file. Each was mapped back to system dof order through
  `Linear System Save Perm`. Each satisfies `||K x - lambda M x|| / ||K x|| < 4e-13`
  and `x^T M x = 1` (checked when the fixture was built).
