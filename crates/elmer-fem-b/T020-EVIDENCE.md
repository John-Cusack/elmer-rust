# T020 evidence (solver-entry pilot, report §8)

Reference: ElmerCSC/elmerfem @ fork `c0b4bc6` (campaign fork; backstop
tree is a serial `-DWITH_MPI=FALSE -DCMAKE_BUILD_TYPE=Debug` build of
that fork — the campaign MPI config is unbuildable here, see below).
Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)` (asserted in
the proof script) + gfortran 13.3.0. Grade driver:
`ports/poisson/proof-poisson.sh` (lane-owned); proof scratch under the
lane data dir (`proof-poisson/`).

## Graded green (1 unit)

| unit | crate | true symbols | differentials |
|---|---|---|---|
| fortran:fem/src/modules/Poisson.F90 | poisson | `poissonsolver_` (only `T` in `nm -D`, both arms) | A (Source branch, 1 iter) + rerun + B (`w` branch + electrodes, 2 iters): byte-identical stdout/stderr/scalars.dat modulo wall-clock lines |

Corpus (clean-room `caseA.sif`/`caseB.sif` on a generated unit-square
mesh): config A exercises the body-force `Source` path with Dirichlet
BCs; config B exercises the `Solver % Variable % name == 'w'` path —
the `BoundaryCondition` internal (ElementInfo, GaussPoints MATMUL/SUM
order, allocatable `Nodes` hand-off), `Positive/Negative Electrode`
`GetLogical` flags, and repeated invocation through the SAVE state (2
steady iterations). Branch-liveness control: the same config with
`Variable = Potential` gives different scalars, proving the `w` path
executes in both arms.

## Refined symbol gate

- `nm -D --defined-only` on the reference `Poisson.so` and on the port
  `libpoisson.so`: exactly `T poissonsolver_` in both; diff of the
  symbol sets is empty. No other `T` exports; accessor objects are
  hidden (`-fvisibility=hidden` + version script).
- The port's SAVE state (`STIFF`/`LOAD`/`FORCE`, first-mesh sizing,
  `AllocationsDone`) is replicated in Rust; the internal procedures'
  `SAVE Nodes` live as wrapper-side SAVE with identical lifetimes; no
  Rust threads (the solver calls entries serially).

## Backstop

- Serial-config baseline (`ctest -L quick`, pristine tree): 451/454,
  3 failures — all SIGSEGV in temperature-BC tests
  (`ConstantBCTemperature`, `ProfileBCTemperature`,
  `ProfileBCTemperatureRobin`), consistent with the K=3 upstream
- Spliced run (port `Poisson.so` in place of the reference, same serial
  tree): `ctest -L quick` 451/454 with the same 3 SIGSEGV failures —
  identical to the pristine baseline (no suite-wide change from the
  splice; no Poisson quick-test regressed).
