# Pristine-vs-pristine determinism (S0, report §4.3 item 9, §5.2 item 9)

Oracle data. `report.tsv` has one row per registered test.

Each set was captured twice on the Release arm with the O2 injection
(`scripts/sweep-capture.sh`, `OMP_NUM_THREADS=1`), in separate run
directories and under different host load. `elmer-oracle sigma-table`
then compared every field of every result file in node space; σ = 0
means the two runs agree bitwise. Verdicts come from
`scripts/determinism-report.py`.

| Set | Tests | Deterministic | Nondeterministic | No fields | Fails in pristine |
|---|---:|---:|---:|---:|---:|
| quick (golden at 8 jobs vs rerun at 4 jobs) | 482 | 466 | 1 | 15 | 0 |
| non-quick `fem/tests` (two runs, 6 jobs each) | 618 | 547 | 0 | 69 | 2 |
| `elmerice/Tests` (two runs) | 65 | 56 | 0 | 2 | 7 |

- **Nondeterministic: `OtherMeshResults_np4` (quick), at 1.95e-16.** It is
  an np=4 MPI test, so the difference is consistent with reduction order
  across ranks, at 1 ulp. Every other test that produces fields is bitwise
  reproducible, including the other 27 MPI quick variants and the 109 parallel
  non-quick ones.
- **Pristine failures, the same in both runs:**
  - **Ice.** Calving2D, Calving3D (timeout), Damage, FSSA_perlin2d,
    Hydro_Coupled, Hydro_SedOnly and InvMeth_AdjRobin. These are exactly
    the ice failures of the pristine ctest baseline, and the held-out seal
    already excludes them.
  - **EM_port_eigen_3D and EM_port_eigen_3D_prisms** segfault (SIGSEGV)
    with `OMP_NUM_THREADS=1` and pass with 16 threads, which is what ctest
    sets. That is a thread-count-sensitive crash in pristine Elmer, so it
    belongs in the known-divergence registry (§4.3 item 8), alongside the
    K=3 class. The candidate may pass these tests where pristine crashes;
    O1 grades them meanwhile.
- **ctest timeouts that are not failures.** The four non-quick tests that
  timed out at 1200 s under ctest `-j` pass in both oracle runs:
  ContactBlunt3DLevelProj (360 s), ContactBlunt3DNormalProj (358 s),
  EM_port_eigen_2ndorder (378 s) and ParStokes_ISMIP_HOM_A010 (616 s).
  They are load-sensitive, not broken.

## Known exceptions

Signed off with https://github.com/John-Cusack/elmer-rust/pull/239.
Grading treats these as exceptions; they are not regressions.

| Test | Exception | Consequence |
|---|---|---|
| `OtherMeshResults_np4` (quick, np=4) | The only test that is not bitwise reproducible: 1.95e-16 between runs (MPI reduction order). | Tier-A bars on it must not rely on bitwise agreement. Its σ_F already exceeds that difference. |
| `EM_port_eigen_3D`, `EM_port_eigen_3D_prisms` (non-quick) | Pristine Elmer segfaults with `OMP_NUM_THREADS=1` and passes with ctest's thread count. | Capture and grade them with ctest's thread count; a candidate that passes at 1 thread does not count as diverging. |

## Runner changes found on the way

- ElmerIce tests use `runTest.cmake` and `RUN_ELMERICE_TEST` (run through
  the `cmake -P` replay), and are looked up in the `elmerice`
  registration.
- Their plugins come from `elmerice/Tests/<dir>`, and their libraries link
  `elmerice/Utils`.
- Directories a test's `CMakeLists.txt` creates with
  `FILE(MAKE_DIRECTORY ...)` (GlaDS `results/`) are created at staging,
  because git keeps no empty directories.

## Reproduce

```sh
S=/home/john/firstmate/data/elmer-s0-oracle; B=target/release/elmer-oracle
for r in a b; do
  scripts/sweep-capture.sh $B $S/build-release $S/src/fem/tests nonquick-$r 6 '!quick'
  SUITE=elmerice scripts/sweep-capture.sh $B $S/build-release $S/src/elmerice/Tests ice-$r 3 elmerice
done
$B sigma-table $S/src/fem/tests nq.json nonquick-a b=nonquick-b > nq.tsv
$B sigma-table $S/src/elmerice/Tests ice.json ice-a b=ice-b > ice.tsv
scripts/determinism-report.py report.tsv nonquick=nonquick-a:nonquick-b:nq.tsv elmerice=ice-a:ice-b:ice.tsv
```
