# Pristine ElmerSolver references

Each directory holds the `case.sif` run by the pristine ElmerSolver
(`build/fem/src/ElmerSolver`), serially, with the mesh `ElmerGrid` wrote and
the resulting reference:

- `mesh.*`: the ElmerGrid output the run used (a 10x10 quadrilateral unit
  square, boundaries 1 bottom, 2 right, 3 top, 4 left).
- `solution.dat`: the last time level of the ASCII `.result` that `Output File =
  "case.result"` writes (17 significant digits), as `node temperature
  potential` lines in node order. For the bubble cases only the first 121 lines,
  the nodes, are kept: the `.result` also lists the bubble dofs.

| case | source | what it exercises |
|---|---|---|
| `plain` | `fem/tests/ThermoElectric` | Seebeck material, current source, Newton from the second iteration |
| `aniso_bc` | `ThermoElectricAnisoBC` | anisotropic `Electric Transfer Coefficient` and `Electric Flux Vector` |
| `stabilize` | `ThermoElectricStabilize` | SUPG stabilization of the Peltier drift |
| `bubble` | `ThermoElectricBubble` | `Element = "p:1 b:4"`, steady |
| `bubble_transient` | `ThermoElectricBubbleTransient` | the same, BDF 1, 5 steps of 0.1 |

The Rust side solves with the default faer LU, so the sifs use `Linear System
Solver = Direct` / `UMFPack` (the stock ones iterate). The bubble cases set
`Nonlinear System Convergence Tolerance = 1e-12`: Elmer solves the bubbles
with the nodal dofs in one system and this crate recovers them after each
solve, so the two agree at convergence, not iterate by iterate. Elmer's
reported norm includes the bubble dofs, so the bubble tests compare the
nodal values only.

The two bubble directories also hold `perm.dat` (`slot row`) and `slots.dat`
(`slot temperature potential`, `x(Perm(slot))` slot by slot, the last time
level): the full `.result` of the same pristine runs, nodes and bubbles. The
tests that take the p dofs from `elmer-pdofs` keep the bubbles in the
global system, as Elmer does, and compare every slot and Elmer's norm
(549.0240106292581 steady, 343.1731238726303 after the fifth step).
