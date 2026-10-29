# Fortran reference systems

Assembled linear systems dumped by the pristine Fortran `ElmerSolver`
(`Linear System Save = Logical True`, `Linear System Save Slot = String "after"`,
`Linear System Save Solution = Logical True`), run serially on scratch copies of
three quick tests that use `Linear System Direct Method = UMFPack`:

| Directory | Source test (`fem/tests/`) | n | nnz |
|---|---|---:|---:|
| `bodydir` | `bodydir` (StatElecSolve, conditional Dirichlet) | 441 | 3721 |
| `bodyload` | `bodyload` (StatElecSolve, nodal load) | 441 | 3721 |
| `pfem` | `PoissonPFEM` (Poisson, `Element = "p:2"`) | 441 | 6561 |

* `linsys_a.dat`: `row col value` (1-based, full precision).
* `linsys_b.dat`: `row value` (the right-hand side).
* `linsys_sol.dat`: `row value`, the UMFPACK solution written by Elmer with
  `ES15.6` (six significant digits), so it is a reference to about 1e-5
  relative, not 1e-12.
