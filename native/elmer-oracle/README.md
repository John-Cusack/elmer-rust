# elmer-oracle

The external oracle for the Rust-native Elmer rewrite (report §4, S0 in
§5.2). It drives pristine Fortran Elmer builds outside ctest, captures
golden data, and grades candidates against it. It never depends on
candidate crates (`tests/isolation.rs`). Changes need oracle sign-off
through firstmate (D11).

## S0 status by item (§5.2)

| Item | Status | Where |
|---|---|---|
| 1 Pinned arms | Release (golden), FMA, OpenBLAS and LISTUSAGE arms built and published | `oracle-builds.json`, `../oracle-arms/`, `../docker/` |
| 2 Direct runner | Registration-driven: 482 quick, 618 non-quick and 65 ice tests run, with `cmake -P` replay for custom `runtest.cmake` | `src/runner.rs`, `src/registry.rs` |
| 3 Golden capture | O2 for 482/482 quick tests on all three arms, and per RUN (`--per-run`) for the 24 multi-RUN tests, with σ_F. O2m (`--o2m`) and iteration counts (`--max-output`) are supported but not captured corpus-wide. | `golden/` |
| 4 Bitwise grader validation | 350 quick tests bitwise through `Reference Norm Tolerance = 0`, vectors, default variables, MATC includes and every RUN of multi-RUN sifs included; 8 grader gaps (7 self-norming solvers, 1 DG), 1 unexplained, 123 out of class (reasons listed) | `refnorm/` |
| 5 σ_F and tolerance table | 5,665 fields, variable-scale floor for near-zero components (exact zeros stay exact); kill check 24/482 (5.0%) with eigen mode | `tolerances.json`, `golden/README.md` |
| Mutation suite | 14,634/14,634 detectable mutations caught (near-zero floor); below-floor and vacuous counts kept visible | `golden/mutation-summary.txt` |
| 6 Capability matrix | 482/482 LISTUSAGE runs | `capability/` |
| 7 O6 restart round trip | 5/5 restart tests bit-stable | `golden/o6.txt`, `scripts/o6-proof.sh` |
| 8 O7 kernel harness | Rust tests call the pristine `GaussPoints` rule builders and `ElementInfo` through a `BIND(C)` shim; Release reproduces 89 rules and 792 `ElementInfo` calls on 21 types (all S1 types) at 0 ulp, the FMA arm does not | `o7/`, `src/o7.rs`, `tests/o7.rs` |
| 9 Held-out seal and determinism | 190 directories sealed with a hash commitment; pristine-vs-pristine determinism of quick, non-quick and ice tests | `heldout/`, `determinism/` |
| 10 Benchmark baselines | **deferred** (see below) | — |
| 11 ABI harvest | 256 staged campaign cases from 108 solver proofs as oracle tests; 236 golden-captured on Release, FMA and OpenBLAS (bitwise deterministic, σ_F per field), 20 pristine failures explained; 15 known-divergence entries from the documented upstream defects | `harvest/` |
| External grading of a candidate | `--solver` plus `grade-tree`; the Rust path (`elmer-run` at main `60d64f8`, no Fortran fallback) grades 189/454 serial quick tests `tier-a-pass`, confirming every serial `elmer-coverage` match; 6 `tier-a-FAIL`. MPI: `elmer-mpi-run` at `9008121` grades 4/28 `_npN` tests `tier-a-pass`, 2 FAIL, 20 refused, 2 hang | `candidate/` |

## Benchmarks: deferred

§4.5 asks for Release+OpenMP and MPI baselines on 4–6 scaled cases: the
median of at least 5 runs on a quiesced, pinned host, with RSS. The only
host available runs the other lanes at load averages of 12–20, so any
timing taken now would measure contention, not Elmer. Numbers from a
loaded host would not meet §4.5 and would mislead later best-vs-best
comparisons, so firstmate deferred the item (decision of 2026-10-30) and
will schedule a quiet window. Nothing in the oracle depends on the
baselines until the first candidate speed claim (S3).
