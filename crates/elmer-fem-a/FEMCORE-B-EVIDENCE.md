# femcore lane B evidence (claim B: 176 of the 353 refused/skipped fem core units)

Reference: Elmer fork `c0b4bc6` (`/home/john/runs/elmer-fork`), verified
byte-identical to the pristine campaign tree @ `9f6af2f85` over all of
`fem/src` (`diff -rq` clean).
Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)` + gfortran
13.3.0 (the campaign pins, ADR-045).

Gate: the **merged** substitution gate — the T014 module profile plus
T017 shape A (the F77-external shape, rustsmith#105 / ADR-049) **and**
T017 shape B (the module-dummy shape, rustsmith#108 / ADR-050). The
re-grade driver is `rustsmith-adapters`' `elmer-regrade` binary; it runs
the campaign's own gate (`fortran_substitution_verdict`, the
reason-carrying form of `CmakeBridge::check_substitutable`) over the
claimed units and records the verdict plus the refusal reason per
offending export.

Claim: `data/elmer-lanes/femcore-claim-a.json` (rustsmith) — the fem
core's refused/skipped units from the lane-fem-a and lane-fem-b ledgers,
sorted by unit id; lane B claims the last half (176 of 353). The lane
ledger for the claimed units is
`crates/elmer-fem-a/FEMCORE-B-LEDGER.md`; the shared ledgers are left
untouched on purpose (this lane owns only its two port crates and its
ledger file — runbook §3 conflict rule).

## Re-grade result (176 claimed units)

| verdict | units |
|---|---|
| admit | 6 |
| refuse | 170 |

| skip reason (current) | offending exports | units |
|---|---|---|
| gate refuses: procedure-not-found-in-gate-parse | 581 | 107 |
| gate refuses: derived-type dummy (needs whole-program layout) | 448 | 166 |
| gate refuses: character dummy (hidden length param) | 29 | 19 |
| gate refuses: non-intrinsic dummy | 7 | 4 |
| gate refuses: undeclared function result (implicit typing) | 2 | 2 |
| gate refuses: assumed-shape rank > 2 | 2 | 2 |
| gate refuses: undeclared dummy (implicit typing) | 2 | 1 |
| gate refuses: optional dummy (hidden presence arg) | 1 | 1 |
| gate admits, but the body needs Elmer whole-program state | — | 1 |
| gate admits vacuously (the scope owns no portable procedure) | — | 2 |
| build-time code generator helper (never linked) | — | 1 |

(A unit is refused once per offending export, so the "offending exports"
column sums over exports and the "units" column counts the units that
carry that reason at least once. The per-unit ledger carries the reason
set of each unit.)

## Ports (2 units, byte-identical corpora)

| unit | crate | true symbols | differential |
|---|---|---|---|
| fortran:fem/src/solve_real.F90 | solvereal | `solvelapack_`, `solvelapacksym_` | 5/5 byte-identical stdout+stderr+exit status |
| fortran:fem/src/mpif_stubs.F90 | mpistubs | 39 `mpi_*_`/`pd*upd_`/`pz*upd_` as `T` + `__mpi_stub_MOD_mpi_stub_init` as `B` | byte-identical stdout+stderr+exit status |

Proof script: rustsmith `tests/femcore_b_proof.sh`.

### solvereal corpus

`SolveLapack` / `SolveLapackSym` are LAPACK dense-solver wrappers: the
LU factorization and the triangular solves stay in LAPACK
(`dgetrf_`/`dgetrs_`, `dpotrf_`/`dpotrs_`), so the numerics are
bit-identical by construction. What the port owns is the wrapper: the
argument marshalling, the `N <= 0` early return and the `info`
diagnostic print.

Corpus (clean-room driver, one reference object and one Rust provider
per binary, the Fortran object dropped and the Rust archive linked
`--whole-archive`):

| case | shape | what it pins |
|---|---|---|
| 1 | 3×3 dense, `b = A*(1,2,3)` | the column-major marshalling and the `DGETRF`/`DGETRS` solve |
| 2 | `N=0` | the early return (no output, no LAPACK call) |
| 3 | singular zero matrix | the diagnostic path: the list-directed `PRINT*` line on stdout, then the `DGETRS` continuation (this unit has no `STOP` — execution continues, and so does the port) |
| 4 | symmetric 2×2 SPD, `b = A*(1,2)` | the `UPLO='L'` Cholesky solve (`DPOTRF`/`DPOTRS`) |
| 5 | singular symmetric (zero matrix) | the `DPOTRF` diagnostic print, then the `DPOTRS` continuation |

Fidelity notes (the corpus pins both):

- `PRINT*` list-directed output is reproduced byte for byte: a leading
  blank, the literal verbatim, then the integer right-justified in a
  field of 12 columns (gfortran's list-directed integer width).
- The Fortran runtime and the Rust provider write to the same file
  descriptor with separate buffers, so the driver flushes unit 6 around
  every wrapper call (ADR-051 §5 precedent); without the flush the two
  streams interleave differently on the two sides (a corpus artifact, not
  a port defect).

### mpistubs corpus

`mpif_stubs.F90` is Elmer's single-rank MPI fallback: seven subroutines
with bodies implement the degenerate communicator (serial success, flag
set/cleared by `init`/`finalize`, size 1, rank 0, split identity, the
requested thread level passed through, `initialized` always `.FALSE.`)
and 32 empty bodies only return. What the port owns is exactly that
behavior — plus the `mpi_stub` module flag as BSS storage.

Corpus (clean-room driver, same substitute-and-compare shape):

| case | what it pins |
|---|---|
| `INIT` | `ierr=0`, flag `.TRUE.` after `mpi_init` |
| `SIZE` / `RANK` | `csize=1`, `rank=0`, `ierr=0` |
| `THREAD` | `prov=req` passthrough, `ierr=0` |
| `SPLIT` | `newcomm=comm` identity, `ierr=0` |
| `INITIALIZED` | `.FALSE.` regardless of the flag, `ierr=0` |
| all 32 empty stubs | each returns (no observable effect, must not trap) |
| `FINAL` | `ierr=0`, flag `.FALSE.` after `mpi_finalize` |

The port performs no I/O, so no unit-6 flush is needed — every printed
byte comes from the driver, identical on both sides by construction.
gfortran `LOGICAL` is 4 bytes (`.FALSE.` = 0, `.TRUE.` = 1); the port's
module flag is a zero-initialized `c_int` in BSS.

Refined merge gate (runbook §2): each port defines all of its unit
symbols (`T` for procedures, `B` for the module flag), and the
intersection with the reference object is exactly those symbols —
**zero STRONG port-level dups** beyond the replaced symbols in both, so
each port replaces its unit's object whole (whole-file portable).

## Backstop (T013-style, lane B)

Fresh tree @ `9f6af2f85`, both lane-B ports spliced over their reference
objects (`ld -r --whole-archive` per object), full relink, `ctest
--output-on-failure -L quick -j16`: **479/482, failures exactly K=3**
(`ConstantBCTemperature`, `ProfileBCTemperature`,
`ProfileBCTemperatureRobin`). `LastTest.log` location + walls in the
rustsmith PR body. Build fact, recorded not hidden: the default
`WITH_MPI=TRUE` build does not link `libmpi_stubs.so` into `elmersolver`,
so the backstop exercises that splice mechanically (link-clean) while
the differential corpus carries the behavioral grade.
