# femcore lane A evidence (claim A: 177 of the 353 refused/skipped fem core units)

Reference: Elmer fork `c0b4bc6` (`/home/john/runs/elmer-fork`).
Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)` + gfortran
13.3.0 (the campaign pins, ADR-045).

Gate: the **merged** substitution gate — the T014 module profile plus
T017 shape A (the F77-external shape, rustsmith#105 / ADR-049). The
re-grade driver is `rustsmith-adapters`' `elmer-regrade` binary (added by
the lane's tool PR); it runs the campaign's own gate
(`fortran_substitution_verdict`, the reason-carrying form of
`CmakeBridge::check_substitutable`) over the claimed units and records
the verdict plus the refusal reason per offending export.

Claim: `data/elmer-lanes/femcore-claim-a.json` (rustsmith) — the fem
core's refused/skipped units from the lane-fem-a and lane-fem-b ledgers,
sorted by unit id; lane A claims the first half (177 of 353), lane B the
rest (176). The lane ledger for the claimed units is
`crates/elmer-fem-a/FEMCORE-A-LEDGER.md`; the shared `LEDGER.md` is left
untouched on purpose (the sibling lane owns the other half of the same
file — runbook §3 conflict rule).

## Re-grade result (177 claimed units)

| verdict | units |
|---|---|
| admit | 25 |
| refuse | 152 |

| skip reason (current) | offending exports | units |
|---|---|---|
| gate refuses: derived-type dummy (needs whole-program layout) | 1562 | 135 |
| gate refuses: procedure-not-found-in-gate-parse | 385 | 70 |
| gate refuses: explicit-shape dummy (later shape) | 175 | 30 |
| gate refuses: optional dummy (hidden presence arg) | 120 | 23 |
| gate refuses: character dummy (hidden length param) | 114 | 23 |
| gate refuses: non-REAL(8) assumed-shape dtype | 76 | 30 |
| gate refuses: undeclared function result (implicit typing) | 56 | 27 |
| gate refuses: assumed-size dummy (later shape) | 43 | 11 |
| gate refuses: procedure dummy (code address) | 24 | 4 |
| gate refuses: pointer dummy (libgfortran realloc contract) | 22 | 8 |
| gate refuses: undeclared dummy (implicit typing) | 22 | 7 |
| gate refuses: assumed-shape rank > 2 | 15 | 8 |
| gate refuses: allocatable dummy (libgfortran realloc contract) | 9 | 5 |
| gate refuses: external scope owns no procedure | 3 | 3 |
| gate refuses: non-8-byte REAL kind | 1 | 1 |
| not in build (never compiled) | — | 8 |
| compiled as an empty stub (optional solver lib off) | — | 5 |
| header-only helper (never ports) | — | 6 |
| gate admits, but the body needs Elmer whole-program state | — | 2 |
| gate admits vacuously (the scope owns no portable procedure) | — | 1 |

(A unit is refused once per offending export, so the "offending exports"
column sums over exports and the "units" column counts the units that
carry that reason at least once. The per-unit ledger carries the reason
set of each unit.)

## Ports (3 units, byte-identical corpora)

| unit | crate | true symbols | differential |
|---|---|---|---|
| fortran:fem/src/SolveBand.F90 | solveband | `solvebandlapack_`, `solvecomplexbandlapack_` | 5/5 byte-identical stdout+stderr+exit status (32 printed lines) |
| fortran:fem/src/SolveSBand.F90 | solvesband | `solvesbandlapack_`, `solvecomplexsbandlapack_` | 4/4 byte-identical stdout+stderr+exit status (17 printed lines) |
| fortran:fem/src/Walls.F90 | walls | `solve_ufric_`, `wall_law_`, `d_wall_law_`, `kewall_` | 12/12 byte-identical stdout (84 values at `ES25.17E3`) — the shape A port (rustsmith#105), corpus re-verified in this lane |

Proof script: rustsmith `tests/femcore_lapack_band_proof.sh`; machine
rows in the grading store (run `femcore-a-lapack-band-grade`).

### SolveBand / SolveSBand corpora

`SolveBandLapack` / `SolveComplexBandLapack` and `SolveSBandLapack` /
`SolveComplexSBandLapack` are LAPACK band-solver wrappers: the
factorization and the triangular solves stay in LAPACK (`dgbtrf_`/
`dgbtrs_`, `zgbtrf_`/`zgbtrs_`, `dpbtrf_`/`dpbtrs_`, `zpbtrf_`/`zpbtrs_`),
so the numerics are bit-identical by construction. What the port owns is
the wrapper: the argument marshalling, the `N <= 0` early return, the
`info` check and the failure path.

Corpus (clean-room driver, one reference object and one Rust provider
per binary, the Fortran object dropped and the Rust archive linked
`--whole-archive`):

| case | shape | what it pins |
|---|---|---|
| 1 | real general band, `KL=KU=1`, `LDAB=4`, N=5, 1 RHS | the band storage (`AB(KL+KU+1+i-j,j) = A(i,j)`) and the solve |
| 2 | real general band, `KL=KU=2`, `LDAB=7`, N=8, 2 RHS | the multi-right-hand-side path (the Fortran passes `LDB=N` for an `X(M,N)` dummy; the port passes the same values) |
| 3 | `N=0` | the early return (no output, no LAPACK call) |
| 4 | complex general band, N=4, 1 RHS | the `COMPLEX(KIND=dp)` path (`ZGBTRF`/`ZGBTRS`) |
| 5 | exactly singular band (zero pivot) | the failure path: the list-directed `PRINT*` line on stdout, `STOP 1` on stderr, exit status 1 |
| 1 (sband) | symmetric band (lower), `KD=1`, `LDAB=2`, N=6 | the `UPLO='L'` storage (`AB(1+i-j,j) = A(i,j)`) and the Cholesky solve |
| 2 (sband) | `N=0` | the early return |
| 3 (sband) | complex symmetric band, N=5 | the `ZPBTRF`/`ZPBTRS` path |
| 4 (sband) | exactly singular symmetric band | the failure path |

The reference objects are built against a stub `Types` module that carries
the campaign's own `dp` definition (`SELECTED_REAL_KIND(12)`, 8 bytes —
the same value the real `Types` module gives). The two units take
exactly one entity from `Types`, that compile-time constant, so the
reference object is the unit's own code; the real `Types` module needs
the generated `config.h` and the Lua module, which the lane does not
build.

Fidelity notes (the corpus pins both):

- `PRINT*` list-directed output is reproduced byte for byte: a leading
  blank, the literal verbatim, then the integer right-justified in a
  field of 12 columns (gfortran's list-directed integer width).
- `STOP 1` writes `STOP 1` to stderr and exits with status 1. The
  corpus compares stdout, stderr **and** the exit status.
- The Fortran runtime and the Rust provider write to the same file
  descriptor with separate buffers, so the driver flushes unit 6 before
  the failure case; without that flush the two streams interleave
  differently on the two sides (a corpus artifact, not a port defect).

Refined merge gate (runbook §2): each port defines all of its unit
symbols as `T`, and the intersection with the reference object is
exactly those symbols — **zero STRONG port-level dups** beyond the
replaced ones, so each port replaces its unit's object **whole**
(whole-file portable).

## Skips that are not gate refusals

- `fortran:fem/src/BackwardError.F90` and `fortran:fem/src/MGPrec.F90`:
  the gate admits the dummy shape, but the procedure bodies read Elmer's
  whole-program state (`Types.GlobalMatrix`, `CRS_MatrixVectorMultiply`,
  `Multigrid.MultiGridSolve`, `IterSolve.FirstCall`). The file-level
  port must define every symbol of the unit, so a Rust port would
  re-implement the solver core, not these units — recorded as skips,
  never forced.
- `fortran:fem/src/modules/EMWaveSolver.F90#emwavesolverutils`: the gate
  admits it vacuously (the scope owns no portable procedure) — the
  ADR-047 precedent, `skipped-vacuous-admit`.
- The C/C++ units in the claim admit the gate vacuously (no non-BIND(C)
  Fortran export) but are not in the build: eight are never compiled
  (`fem/obsolete` has no `CMakeLists.txt`; `elmer_lua_iface.c` is added
  only under `IF(WITH_LUA)`, which defaults to FALSE; `fft_test.c` and
  the four `view3d/LinearSolve*.c` files are orphan sources; the
  `minpack_kinds.F90` module is added only under `IF(WITH_EXTOPTIM)`),
  five compile as empty stubs (the optional solver libraries HYPRE,
  SuperLU, AMGX, CHOLMOD and Rocalution are all off by default, so the
  guarded bodies are never built), and six are header-only helpers.

## Timing (ADR-045)

lane-fem-a has no `bench_sample_unit` (null at dispatch) and the
lane-fem-b sample unit (`fortran:fem/src/modules/HeatSolve.F90`) is not
in this claim, so **no timing sample runs** for these units: the bench
row stays unmeasured, never a guess. The headroom rows for the Elmer
campaign (`data/headroom/elmer-*.json`) are not in the tree at the lane's
done-gate, so no unit was skipped for "will not pay".
