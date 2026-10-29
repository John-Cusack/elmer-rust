# HUTI reference fixtures

Captured 2026-10-29 from the pristine ElmerSolver v26.2 build at
`/home/john/runs/elmer-fork/orig_src/build` (fhuti built `-g`, no `-O`, no
OpenMP; reference BLAS 3.12). Each run was serial, `OMP_NUM_THREADS=1`.

## Systems

| Directory | Source | Change to the quick test | Size |
|---|---|---|---|
| `heataniso/` | `fem/tests/HeatAniso` | FluxSolver (Solver 2) removed | 1001 rows, 8401 nnz, symmetric |
| `heataniso_conv/` | `fem/tests/HeatAniso` | as above, plus `Convection = Constant` with `Convection Velocity 1 = 20`, `Convection Velocity 2 = 5`, `Heat Capacity = 1` | 1001 rows, nonsymmetric |

The mesh was generated with `ElmerGrid 1 2 cross`. Solver 1's linear-system
keywords were replaced by:

```
Linear System Solver = String "Iterative"
Linear System Iterative Method = String "<method>"
Linear System Preconditioning = String "<None|Diagonal|ILU0|ILU1>"
Linear System Scaling = Logical False
Linear System Residual Output = Integer 1
Linear System Max Iterations = Integer <2000 or 500>
Linear System Abort Not Converged = Logical False
Linear System Convergence Tolerance = Real 1.0e-8
Linear System Save = Logical True
Linear System Save Slot = String "linear solve"
Linear System Save Solution = Logical True
```

Scaling is off so the dumped matrix is exactly what HUTI receives; the save
slot dumps it inside `SolveLinearSystem` before the iterative call. The saved
initial guess (`linsys_sol.dat`) was all zeros in every run, and the dumped
matrix and right-hand side were byte-identical across all methods of a case.

- `a.txt`: `linsys_a.dat` as `row col value` (1-based, storage order), values
  re-printed as shortest round-trip decimals.
- `b.txt`: `linsys_b.dat`, one value per row.
- `<method>_<preconditioning>.out`: every HUTI iteration line from the run's
  stdout, verbatim, including the line repeated at exit.

Max Iterations was 2000 for the `heataniso` runs except TFQMR and CG/ILU0
(500), and 500 for every `heataniso_conv` run. The ILU runs are checked with
elmer-sparse's `Ilu` as the preconditioner.
