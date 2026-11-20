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

## Which solvers are in the class

A solver is checked when its sif sets `Reference Norm`, the run is serial
and single-`RUN`, and the grader can rebuild the solver's `Values`:
- **Variable.** The explicit `Variable`, scalar or vector. Only one-dof
  variables are saved, so a vector (`-dofs N name`, or `name[a:2 b:1]`)
  is rebuilt from its saved components: component `j` at Perm index `k`
  is `Values(dofs*(k-1)+j)`. The component names and order come from the
  result preamble.
- **Default variable.** A solver without `Variable` takes Elmer's legacy
  default when its equation name selects one (`navier-stokes`,
  `heat equation`, `stress analysis`, `mesh update`,
  `magnetic induction`; `MainUtils.F90`). Otherwise it takes the single
  saved root variable that the preamble assigns to its `Equation`.
  Coupled solvers, equations shared by several solvers, and several
  candidates stay out of class.
- **Includes.** A top-level `include <file>` that only defines MATC
  variables (`$` lines: `params.sif`, the `mesh.names` ElmerGrid writes)
  is expanded as Elmer reads it. The file is looked up in the golden run
  dir, then the source dir. An include that defines sections stays out of
  class, because the override goes into the main sif.
- **Norm options.** Default, plus the velocity-only norm that FlowSolve,
  IncompressibleNS and FilmFlowSolver set internally
  (`Nonlinear System Norm DOFs = DOFs - 1`). `ComputeNorm` then sums each
  component on its own (`Norm = Norm + SUM(x(j::Dofs)**2)`), and the
  grader follows that order.

## Result

| Verdict | Tests |
|---|---:|
| **bitwise:** Elmer accepts the grader's norm at tolerance 0 | **338** (370 solvers) |
| grader gap: the solver sets its own norm (FluxSolver, VectorHelmholtzNodal, WhitneyAVSolver) | 7 |
| grader gap: DG field (PoissonDG) | 1 |
| unexplained mismatch (PhaseChange3, solver 3 `PhaseSurface`) | 1 |
| out of the validated class | 135 |

Before this widening the class held 223 tests (238 solvers). All 223 are
still bitwise, each with at least as many solvers as before.

**Kill check (§5.2).** Every test in the class that ran without a bitwise
match is counted: 9/482 = 1.9% of the quick tests, under the 2% bar.
Eight of them are grader gaps. In each, the values are all in the file, but
the solver's norm is not `ComputeNorm` of them:
- 34 solver entries assign `Solver % Variable % Norm` themselves (grep of
  the pristine `fem/src/modules`, listed in `SELF_NORM` in
  `src/refnorm.rs`). When only such solvers fail, the test is rerun
  without their overrides (`<test>.without-self-norm`), so its other
  solvers still count. CoilSolver1 and fluxsolver-multimesh pass that way.
- PoissonDG's norm runs over DG unknowns that the nodal rebuild does not
  follow.

PhaseChange3's transient phase-change solver reports a norm 5× the saved
field's. It is the one open item for this check.

## Out of class (not validated, not failed)

Each test is listed under the first reason in its row.

| Reason | Tests |
|---|---:|
| MPI | 28 |
| several `RUN` blocks | 25 |
| `Show Norm` / `Show Norm Index` (a SaveScalars value is the norm) | 23 |
| custom `runtest.cmake` (replay) | 13 |
| variable in several result files (several meshes) | 12 |
| variable not in any result file (global, or not output) | 9 |
| duplicate `Solver N` sections (merged by Elmer) | 7 |
| no `Variable`: coupled solver | 5 |
| no `Variable`: several saved variables of its equation | 6 |
| no `Variable`: no saved variable of its equation | 4 |
| other (`-nooutput`, several solver runs) | 2 |
| non-nodal unknowns | 1 |

So item 4 is **not** done for every captured test. The validated class
covers 338 of the 482 quick tests. Widening it further means per-run
checks for multi-`RUN` sifs, the SaveScalars-style norms, and per-mesh
result files.

## Reproduce

```sh
S=/home/john/firstmate/data/elmer-s0-oracle
OMP_NUM_THREADS=1 target/release/elmer-oracle refnorm-validate \
  $S/build-release $S/src/fem/tests <golden-root> <out-root>
```

The table is written to `<out-root>/refnorm.tsv`.
