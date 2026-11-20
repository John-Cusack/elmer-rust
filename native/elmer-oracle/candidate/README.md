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

## The partitioned Rust path: MPI quick tests

`rust-elmer-mpi-run-9008121.tsv` grades `elmer-mpi-run` built from main
`9008121` on the 28 registered `_npN` quick tests. The runner starts each
under the arm's `mpiexec -n N` (Open MPI 4.1.6), with `ELMER_NATIVE=1`
and `ELMER_FORTRAN_SOLVER` unset. `grade-tree` grades every per-partition
golden file (`<mesh>/case.result.<rank>`) against the candidate's file of
the same name.

| Verdict | Tests |
|---|---:|
| `tier-a-pass` | **4** (AngleMetisLayer_np2/3/4, ParallelSlices_np2) |
| `tier-a-FAIL` | 2 (passive_np3, passive_np8) |
| `no-candidate-output` (refused) | 20 |
| `candidate-timeout` | 2 (radiation_np2, radiation_np4) |

- **Agreement with `elmer-coverage`:** exact. Its 4 MPI matches pass. Its
  2 mismatches fail on the same fields (passive_np3 `temperature` step 0
  at 1.7e-11 against a 6.7e-12 bar; passive_np8 step 1 at 1.5e-9 against
  1.46e-9). Its refusals refuse.
- **HeatControlExplicit_np3** stops on "procedure AlternatingSource is not
  registered". It runs in the test's own directory name, as the serial
  run does, but the serial binary resolves the D12 fixture and
  `elmer-mpi-run` does not. This is a candidate-side gap.
- **radiation_np2/np4** pass `elmer-coverage`'s static screen, but the
  candidate wrote nothing in 588 s and was killed (pristine: 1.9 s and
  2.6 s). Those two rows were added by hand: `grade-tree` skips a
  test whose run was killed before the `run` link was made. Pristine
  writes no result file for either, so neither could be field-graded.

```sh
printf '#!/bin/sh\nexec env -u ELMER_FORTRAN_SOLVER ELMER_NATIVE=1 HWLOC_COMPONENTS=-gl <target>/release/elmer-mpi-run "$@"\n' > ElmerSolver-rust-mpi
awk -F'\t' '$3>1 {print $1}' <golden-root>/registry.tsv > mpi.txt
scripts/grade-candidate.sh <oracle-bin> ElmerSolver-rust-mpi \
  $S/build-release $S/src/fem/tests <golden-root> tolerances.json <out> 2 < mpi.txt
```

```sh
printf '#!/bin/sh\nexec env -u ELMER_FORTRAN_SOLVER ELMER_NATIVE=1 <target>/release/elmer-run "$@"\n' > ElmerSolver-rust
awk -F'\t' '$3==1 {print $1}' <golden-root>/registry.tsv > serial.txt
scripts/grade-candidate.sh <oracle-bin> ElmerSolver-rust \
  $S/build-release $S/src/fem/tests <golden-root> tolerances.json <out> 6 < serial.txt
```
