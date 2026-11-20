# Grading a candidate through the external grader (report §4.3 item 1)

`scripts/grade-candidate.sh` runs each listed test with a candidate
`ElmerSolver` (`run --solver <path>`). Everything else comes from the
pristine arm: the registration, ElmerGrid, staging, and the `cmake -P`
replay, where the candidate replaces the recorded `ELMERSOLVER_BIN`. The
script injects the O2 `Output File`, then `elmer-oracle grade-tree` grades
every golden result file against the candidate's, using the committed
Tier-A bars (`tolerances.json`). The candidate's `TEST.PASSED` and PASS
lines are never read.

Per-test verdicts: `tier-a-pass`, `tier-a-FAIL` (a field over its bar,
or a structural mismatch), `no-candidate-output` (no result file to grade),
and `no-golden-fields` (the golden carries no field variables, so the
test is graded by O1 and its logs only).

## Control: the pristine solver as candidate

`control-pristine.tsv`: the Release arm's own `ElmerSolver`, passed as a
candidate, on 7 tests that include two `cmake -P` replays (ChannelFlowCoulomb,
ElasticStabilized) and an MPI variant (AngleMetisLayer_np2). The result is
**7/7 `tier-a-pass`**. Grading the FMA arm's capture as a candidate gives
467/467 `tier-a-pass`, as it must, since σ_F is built from that arm.

## The Rust path: 189 serial quick tests graded tier-a-pass

`rust-elmer-run-60d64f8.tsv` grades the `elmer-run` binary built from main
`60d64f8`. It runs `ElmerSolver`'s process behaviour with `NativeEngine`
and `SifFrontend`, under `ELMER_NATIVE=1` with `ELMER_FORTRAN_SOLVER`
unset, so no case can fall back to Fortran. It ran on all 454 serial
registered quick tests:

| Verdict | Tests |
|---|---:|
| `tier-a-pass` | **189** |
| `tier-a-FAIL` | 6 |
| `no-golden-fields` | 15 |
| `no-candidate-output` (refused, or wrote no result file) | 244 |

Against `elmer-coverage/results.tsv` (main `fc75eec`, in-process harness,
self-graded `TEST.PASSED` plus its field check):
- **Every serial coverage match is confirmed by the external grader:**
  182 `tier-a-pass` and 7 `no-golden-fields` (ElementUnitTest_*,
  L2ndDerivatives).
- 5 tests coverage refused and 2 it called mismatches now pass. They are
  candidate PRs merged after `fc75eec`: NonnewtonianChannelFlow_vec (and
  its SD_ variant), ThermoElectricBubble(Transient), ElasticConstantSystem,
  BlockLinElast2bMany and HeatGap.
- The 5 coverage mismatches still fail (1dtests_MeshLevels,
  ElasticAnisotropic2D, ElasticStressLoad, ElmerGridExtrudeOffset,
  bentonite), on the same fields.
- **ElasticIntegrationOrder:** coverage refused it at `fc75eec`. The binary
  now runs it, and `displacement 1` is 4.6% off (bar 6.1e-14). Its own
  `Reference Norm` check fails as well. It is a wrong answer, not a
  refusal, for the candidate lanes to fix or refuse.

**Test-local fixtures need the ctest directory name.** `elmer-fixtures`
keys its Rust ports of the D12 test-local plugins by the case
directory's name, which under ctest is the test's own directory. The
grader therefore runs the candidate in `<out>/<test>/<test dir>`, with
`run` linking to it. In a directory named `run`, the 25 fixture cases
(1dtests, BlockPoisson*, CoupledPoisson*, ModelPDE*, PoissonPFEM, …)
refused with "procedure not registered".

MPI variants (`_npN`, 28 tests) are not graded here. They need the
candidate's `elmer-mpi-run` under the arm's `mpiexec`.

```sh
printf '#!/bin/sh\nexec env -u ELMER_FORTRAN_SOLVER ELMER_NATIVE=1 <target>/release/elmer-run "$@"\n' > ElmerSolver-rust
awk -F'\t' '$3==1 {print $1}' <golden-root>/registry.tsv > serial.txt
scripts/grade-candidate.sh <oracle-bin> ElmerSolver-rust \
  $S/build-release $S/src/fem/tests <golden-root> tolerances.json <out> 6 < serial.txt
```
