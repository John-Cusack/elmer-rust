# lane-fem-a evidence (batch 1)

Reference: Elmer @ `9f6af2f85`, fork `lane/fork`, `Debug`, defaults
(`-DRESET_ELMER_REVISION=FALSE -DELMER_FEM_REVISION=9f6af2f85`,
branch `devel`; `ELMER_FEM_COMPILATIONDATE` is build-date by nature).
Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)`.
ccache: shared `elmer-lanes-cache`, `CCACHE_BASEDIR=$LANE_ROOT`, sloppiness
per `data/elmer-lanes/ccache.env` (ADR-045).

## Backstop (pristine)

`ctest --output-on-failure -L quick` (the oracle invocation: 482 tests).
Log: `lane/ctest_backstop.log`. Result: **479/482, failures exactly
K=3** (`ConstantBCTemperature`, `ProfileBCTemperature`,
`ProfileBCTemperatureRobin`) — campaign baseline reproduced.
Failure set: `lane/backstop_failures.txt`.
Parallel-validation: `ctest -L quick -j16` on the pristine tree gives the
identical failure set (`lane/ctest_j16.log`) in ~7 min; all per-unit
suites use `-j16` on this evidence.

## Batch-1 units (11)

| unit | crate | true symbols | local differential |
|---|---|---|---|
| c:fem/src/view3d/LUDecomp.c | ludecomp | ludecomp, mtrinv | byte-identical stdout+stderr, 5 cases |
| c:fem/src/view3d/VectorUtil.c | vectorutil | RotateVector, RotateVectorTranspose, GetMatrixToRotateVectorToZAxis, CylinderNormal | byte-identical, incl. exit(1) discriminant path |
| c:fem/src/view3d/second.c | second | second | same-syscall functional |
| c:fem/src/f_stubs.c | f_stubs | rename_c (HAVE_F_ETIME/HAVE_F_FLUSH compile the rest out) | round-trip green |
| c:fem/src/CPUTime.c | cputime | cputime, realtime, cpumemory | same-syscall functional |
| c:fem/src/GFortranLocaleFix.c | gfortranlocalefix | elmer_fix_numeric_locale | TBD |
| c:fem/src/binio/binio.c | binio | 16 × bin*_c + strerrorf_c | byte-identical, all fns + swap/EOF/errno paths |
| fortran:fem/src/minpack_kinds.F90 | minpack_kinds | (none — parameter-only module) | n/a (no symbols) |
| fortran:fem/src/ElementDescription.F90 | elementdescription | (none — facade) | n/a |
| fortran:fem/src/MeshUtils.F90 | meshutils | (none — facade) | n/a |
| fortran:fem/src/SolverUtils.F90 | solverutils | (none — facade) | n/a |

## Batch-2 units (5, local differential green, machine grade pending)

| unit | crate | true symbols | local differential |
|---|---|---|---|
| c:fem/src/umf4_f77wrapper.c | umf4_f77wrapper | 14 × umf4*/umf4_l_* | byte-identical vs vendored UMFPACK (di path) |
| c:fem/src/view3d/2d_4node.c | _2d_4node | 18 shape-fn symbols | byte-identical incl. libc-fprintf %g paths |
| c:fem/src/view3d/16node_quad.c | _16node_quad | 9 shape-fn symbols | byte-identical |
| c:fem/src/view3d/ShaftCull.c | shaftcull | 6 shaft symbols | byte-identical (Walk, planes, counting) |
| c:fem/src/view3d/ContourInteg.c | contourinteg | 5 contour symbols | byte-identical (prepare/eval/poly/counts) |

Machine grades (splice + rebuild + quick suite + matc probe + nm scan),
rows in `lane/store.db` run `fem-a-batch1`:
- `c:fem/src/view3d/LUDecomp.c`: suite failures == K=3
  (`lane/suite_ludecomp.log`, `fails_ludecomp.txt`), matc probe `3`
  (`probe_ludecomp.txt`), wall 198 s. Remaining batch-1 units grading.
- Held-out: no held-out workloads are frozen for the Elmer package
  (`heldout_workloads.json` is empty), so no HeldoutDivergence row is
  recorded — never a vacuous pass.
