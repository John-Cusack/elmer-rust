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
- **Several `RUN` blocks.** Elmer compares after every RUN with the
  keywords read so far. So run *j* checks the solvers that have a
  `Reference Norm` in the sif up to RUN *j*. Their golden norm comes from
  the per-RUN golden (`--per-run`: `case.run<j>.result`, with the last RUN
  in `case.result`), and the override is written as a one-liner just
  before RUN *j*. The rerun uses the same per-RUN capture. A multi-RUN sif
  that names its own `Output File`, or sets `Output Intervals = 0`, stays
  out of class: in the second case each saved field is not the solution
  of its RUN (RotatingFlow's holds FlowSolve's 1e-6 initial value).
- **MPI (`_npN`).** Each rank writes `<file>.<rank>`, and the norm is the
  parallel branch of `ComputeNorm`. By default every rank sums all its
  local values, interface nodes included, and the sums and counts are
  `MPI_ALLREDUCE`d. With `Nonlinear System Consistent Norm` a rank sums
  only the rows it owns. A node's owner is the first partition its line
  in ElmerGrid's `part.<rank+1>.shared` lists (Elmer's `Neighbours(1)`,
  `ReadSharedFile`); any other node belongs to the rank itself. Rows map
  to nodes through the Perm. The rank sums combine in Open MPI's
  recursive-doubling grouping (coll/tuned for one double): the first
  `2·(p − 2^k)` ranks pair up, then a doubling tree. Every bitwise MPI
  test passed with that grouping, the first one tried; the left-to-right
  fallback was never needed. Without a `partitioning.<np>` directory
  (Parallel Slices) each rank solves a copy of the whole mesh, and Elmer's
  norm is rank 0's serial one. Out of class: a local mesh whose node count
  differs from the partition files (extruded, refined, or with
  view-factor nodes added at load), and a Perm beyond the mesh nodes
  (edge, p or DG unknowns), because the partition files give the owners
  of mesh nodes only.
- **Several meshes.** When the variable is in several result files, the
  grader picks the solver's own:
  - the directory of its solver-level `Mesh` (multimesh,
    OtherMeshResults);
  - after adaptive refinement, the last refined mesh, `RefinedMesh<N>`
    with the highest N;
  - for a CutFEM solver, its cut mesh `<mesh>-cutfem` (CutFEMPoisson1:
    Elmer's 0.374375920 is the norm there, 0.234 on the background mesh);
  - otherwise the `Mesh DB` mesh.

  The rerun creates each `RefinedMesh<N>/` it needs, as `capture` does.
  A CutFEM solver whose cut mesh is not saved (the library CutFEM of
  `ModelPDE`) stays out of class.
- **Global variables** (`Variable = -global ...`). `SaveResult` skips
  globals unless `Output Global Variables` is set. Their golden is
  captured with `--globals` (O2 plus that keyword), and the rerun matches
  it. A global's `Values` is its NULL-Perm block as written. This applies
  to OdeSolver1stOrder, OdeSolver2ndOrder and SlaveSolverSlots, whose
  `--globals` golden sits in the refnorm golden root beside the per-RUN
  one.
- **Replayed tests (custom `runtest.cmake`).** These run through
  `cmake -P` with the arm's recorded registration, as `capture` does, with
  the overrides in the startinfo sif. A script that runs further sifs
  writes their fields to `<stem>.result`, so the grader reads the startinfo
  sif's `case.result` (all files when there is none).
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
| **bitwise:** Elmer accepts the grader's norm at tolerance 0 | **383** (490 solver checks) |
| grader gap: the solver sets its own norm (FluxSolver, VectorHelmholtzNodal, WhitneyAVSolver) | 7 |
| grader gap: DG field (PoissonDG) | 1 |
| unexplained mismatch (PhaseChange3, solver 3 `PhaseSurface`) | 1 |
| out of the validated class | 90 |

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
- FixTangentVelo writes into the flow solution after the flow solver
  took its norm (`FixTangentVelo.F90:208`), so that solver's saved field
  is not the one Elmer normed. As with the self-norming solvers, the test
  is rerun without that override (`<test>.without-modified`), so its other
  solvers still count (`MODIFIES_OTHER` in `src/refnorm.rs`).
- PoissonDG's norm runs over DG unknowns that the nodal rebuild does not
  follow.

PhaseChange3's transient phase-change solver reports a norm 5× the saved
field's. It is the one open item for this check.

## Out of class (not validated, not failed)

Each test is listed under the first reason in its row.

| Reason | Tests |
|---|---:|
| `Show Norm` / `Show Norm Index` (a SaveScalars value is the norm) | 28 |
| no `Variable`: several saved variables of its equation | 10 |
| no `Variable`: no saved variable of its equation | 10 |
| variable not in any result file (edge fields, a coupled vector Elmer does not output, matrix-only tests) | 6 |
| MPI: non-nodal unknowns (no owners in the partition files) | 8 |
| duplicate `Solver N` sections (merged by Elmer) | 7 |
| MPI: local mesh changed after load (extruded, refined, view factors) | 6 |
| no `Variable`: coupled solver | 5 |
| several `RUN` blocks with `Output Intervals = 0` | 3 |
| CutFEM: the cut-mesh solution is not saved (library CutFEM) | 2 |
| other (`-nooutput`, several solver runs, a multi-RUN sif with no field solver, an MPI variable in two per-rank result sets) | 4 |
| non-nodal unknowns | 1 |

So item 4 is **not** done for every captured test. The validated class
covers 383 of the 482 quick tests. Widening it further means the
SaveScalars-style norms, and parallel ownership of
non-nodal unknowns and of meshes Elmer changes after loading.

## Reproduce

```sh
S=/home/john/firstmate/data/elmer-s0-oracle
OMP_NUM_THREADS=1 target/release/elmer-oracle refnorm-validate \
  $S/build-release $S/src/fem/tests <refnorm golden-root> <out-root>
```

The table is written to `<out-root>/refnorm.tsv`. The refnorm golden root
is the per-RUN golden (`--per-run`), with the three global-variable tests
captured with `--per-run --globals`:
`/home/john/elmer-laptop-lanes/elmer-lane-oracle/refnorm-golden/` on the
oracle host.
