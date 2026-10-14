# lane-fem-a evidence (batch 1)

Reference: Elmer @ `9f6af2f85`, fork `lane/fork`, `Debug`, defaults
(`-DRESET_ELMER_REVISION=FALSE -DELMER_FEM_REVISION=9f6af2f85`,
branch `devel`; `ELMER_FEM_COMPILATIONDATE` is build-date by nature).
Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)`.
ccache: shared `elmer-lanes-cache`, `CCACHE_BASEDIR=$LANE_ROOT`, sloppiness
per `data/elmer-lanes/ccache.env` (ADR-045).

## Backstop (pristine)

`ctest --output-on-failure -L quick` (the oracle invocation: 482 tests).
Log: `lane/ctest_backstop.log`. Result: TBD.

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

Machine grades (splice + rebuild + quick suite + matc probe + nm scan):
TBD per `lane/grade_<crate>.log`, rows in `lane/store.db` run `fem-a-batch1`.
