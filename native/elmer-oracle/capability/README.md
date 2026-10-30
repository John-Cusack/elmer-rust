# Capability matrix of the quick suite (S0, report §5.2 item 6)

Oracle data. Changes only by ADR with captain sign-off.

Built from the LISTUSAGE arm: `-DDEVEL_LISTUSAGE` plus the one-line
ListFind-miss log (`oracle-arms/patches/listfind-miss-log.patch`, D6,
proven numerically inert in https://github.com/John-Cusack/elmer-rust/pull/128).
Every registered quick test (482) is run with `Max Output Level` raised to
12. At that level the logs carry Elmer's unused-keyword report (level 10)
and the miss trace (level 12). **482/482 runs exit cleanly.**

| File | Content |
|---|---|
| `quick.json` | One line per test: procedures loaded, injected libraries, sif keywords used and unused per section kind, queried-absent keyword ids, grading classes (`golden/classes.tsv`). |
| `keywords.tsv` | Per (section kind, keyword): tests whose sifs set it, and tests where Elmer never read it. |
| `queried-absent.tsv` | Keywords the code queried that no sif set, with the number of tests (row `id` = the ids in `quick.json`). |
| `modules.tsv` | Per solver library: tests that load it, and tests where it was injected rather than named in a sif `Procedure`. |

## What a present-keyword check would miss

- **Injected solvers: 192 tests.**
  - `Post File` adds `ResultOutputSolve` in 168 tests.
  - Legacy `Equation = ...` names inject FlowSolve (15), StressSolve (11), MeshSolve (5), HeatSolve (3), HelmholtzSolve (1) and HelmholtzBEM (1). That is 36 tests; the report's count is 39.
- **Defaults: 9,976 distinct keywords queried while unset.** That is 498 per test on average, including variable-named ones (`<var> load`, `<var> condition`). These are the hidden defaults a Rust capability check has to answer: the code reads them, so their absence means "use the default".
- **Unused sif keywords: 1,530 (section, keyword) pairs in 350 tests** were set but never read: dead input that a candidate may ignore. The raw report lists 2,001 pairs in 473 tests; the rest is the `Max Output Level` artifact below.

## Caveats

- **`Max Output Level` is listed as unused in nearly every test.** Elmer reads it before the LISTUSAGE counters are reset at the end of `LoadInputFile`. That is an artifact of the instrumentation, not dead input. When a sif lacked the keyword and only the sweep injected it, the entry is dropped.
- **Present keywords come from the oracle's sif parse** (sections, `::` one-liners, tables and `/* */` blocks skipped, as in `LoadInputFile`). Keywords that reach the lists another way (MATC-generated, `Include` files, solver-added) can show as unused without being present. In aggregate only two keywords do this (`keywords.tsv` rows with `tests_unused` > `tests_present`): `max output level` (above) and `linear system abort not converged` (156 vs 155).
- **The miss trace has no section name** (`ListFindMiss: <name>`), so queried-absent keywords are per test, not per section.

## Reproduce

```sh
B=target/release/elmer-oracle; S=/home/john/firstmate/data/elmer-s0-oracle
scripts/capability-sweep.sh $B $S/build-listusage $S/src/fem/tests cap 8 quick
scripts/capability-matrix.py $S/src/fem/tests cap golden/classes.tsv capability
```
