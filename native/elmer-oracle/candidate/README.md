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

## The Rust path today: nothing gradable yet

`rust-elmer-run-77e7338.tsv` covers the `elmer-run` binary built from main
`77e7338` (`ElmerSolver` behaviour with `NativeEngine`, `ELMER_NATIVE=1`,
no Fortran fallback). It was run on the 163 serial tests that
`elmer-coverage` reports as matching, plus one replay test. The result is
**157 `no-candidate-output`, 7 `no-golden-fields`, 0 graded.** Two candidate-side gaps block
external grading:

1. **The process-boundary binary has no front end.** `elmer-run`'s `main`
   uses `NativeEngine::default()`, which refuses every case with "no .sif
   and mesh front end (elmer-sif, elmer-mesh) in this build". The 169
   matches in `elmer-coverage/COVERAGE.md` come from the in-process
   harness, which installs `SifFrontend` itself, and from its
   self-graded `TEST.PASSED`.
2. **The native engine does not write `Output File`.** It writes no `.result`
   files, so even with a front end there are no fields for the grader to
   read. Its VTU output (`Post File`) reaches only 57 of the 163 tests.

These belong to the candidate lanes (S1 front end and output), not the
oracle. Once the binary writes `Output File` with `Binary Output`, the same
command grades it with no oracle change:

```sh
scripts/grade-candidate.sh <oracle-bin> <candidate ElmerSolver> \
  $S/build-release $S/src/fem/tests <golden-root> tolerances.json <out> 8 < tests.txt
```
