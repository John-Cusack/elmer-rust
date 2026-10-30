# Bitwise grader validation (S0, report §4.3 item 2, §5.2 item 4)

Oracle data. `quick.tsv` has one row per registered quick test.

`elmer-oracle refnorm-validate` checks the grader against Elmer itself:
1. For every solver a test checks with `Reference Norm`, the grader
   recomputes the norm from the golden `.result` (last saved step,
   `Values` and `Perm` rebuilt in dof order, serial `ComputeNorm`).
2. It writes that value back with 17 significant digits and
   `Reference Norm Tolerance = 0`. The override goes into the `Solver N`
   section, or right after the solver's last `Solver N :: Reference Norm…`
   one-liner, ahead of any `RUN`.
3. It reruns the test with the same O2 instrumentation as the golden.

Elmer's own `CompareToReferenceSolution` passes only on exact equality. A
pass therefore proves four things: that the captured binary field holds
the solver's exact values, that the grader maps them back through the
Perm, that its `ComputeNorm` sums them in Fortran's order, and that the
O2 injection does not perturb the run.

## Result

| Verdict | Tests |
|---|---:|
| **bitwise:** Elmer accepts the grader's norm at tolerance 0 | **223** (238 solvers) |
| grader gap: edge or DG field (MgDynDiscontBCs, mgdyn_steady_2ndorder, PoissonDG) | 3 |
| unexplained mismatch (PhaseChange3, solver 3 `PhaseSurface`) | 1 |
| out of the validated class | 255 |

**Kill check (§5.2).** Every test in the class that ran without a bitwise
match is counted: 4/482 = 0.8% of the quick tests, under the 2% bar. In
the 3 grader-gap tests the values are all in the file, but the solver's
norm over edge or DG unknowns follows a path the grader does not rebuild
yet (S5 work). PhaseChange3's transient phase-change solver reports a norm
5× the saved field's. It is the one open item for this check.

## Out of class (not validated, not failed)

| Reason | Tests |
|---|---:|
| solver has no explicit `Variable` (legacy equation or default name) | 45 |
| variable not in any result file (not an output variable) | 43 |
| multi-dof variable (`-dofs N`, N > 1) | 29 |
| MPI | 28 |
| several `RUN` blocks | 25 |
| `Nonlinear System Norm *` / `Norm Permutation` options | 23 |
| vector variable (`name[a:1 b:2]`) | 18 |
| custom `runtest.cmake` (replay) | 13 |
| sif uses `Include` | 12 |
| variable in several result files (several meshes) | 9 |
| duplicate `Solver N` sections, or `Show Norm Index` (SaveScalars norms) | 7 |
| other (`-nooutput`, several solver runs, unsaved non-nodal values) | 3 |

So item 4 is **not** done for every captured test. The validated class
covers 223 of the 482 quick tests. Widening it means reconstructing
multi-dof and vector norms (interleaved dofs), legacy default variable
names, the norm options, and per-run checks for multi-`RUN` sifs.

## Reproduce

```sh
S=/home/john/firstmate/data/elmer-s0-oracle
OMP_NUM_THREADS=1 target/release/elmer-oracle refnorm-validate \
  $S/build-release $S/src/fem/tests <golden-root> <out-root>
```

The table is written to `<out-root>/refnorm.tsv`.
