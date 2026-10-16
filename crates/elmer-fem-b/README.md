# elmer-fem-b

Lane-fem-b ports: Elmer `fem/src/modules/` (184 files, 218 recon scopes:
all Fortran).

Round 2 (T014 relaxation landed): 1 scope graded green (`HarmUtils`),
217 scopes recorded skips with the relaxed gate's reason.

## Layout

```text
crates/elmer-fem-b/
  ports/harmutils/ # staticlib port of HarmonicSignals.F90#harmutils
                   # (SinSum/CosSum: scalars + assumed-shape REAL(8))
  LEDGER.md        # all 218 scopes: graded / skipped with reasons
  EVIDENCE.md      # graded-run evidence (corpus, gates, splice, timings)
```

## Rules (campaign standing pins)

- Toolchain: pinned nightly `1.100.0-nightly feaadeeac 2026-09-19`
  (`rust-toolchain.toml` + hash assert in the grade driver).
- The port exports the unit's **true object symbols** (from
  `nm -g --defined-only` on the reference build:
  `__harmutils_MOD_sinsum`, `__harmutils_MOD_cossum`), never the recon
  module name — see EVIDENCE.md.
- Fortran: only the T014-admitted shape ports (COMMON-free,
  real-`MODULE` procedures whose dummies are scalar intrinsics or
  assumed-shape `REAL(8)` rank 1–2). Everything else is a recorded skip
  (ADR-027; `check_substitutable` pinned message), never forced —
  including the 3 vacuous admits (empty / data-only modules own no
  portable procedures).
- Grading is per recon scope (`file#scope`); the splice replaces the
  whole `.o`, so a file whose sibling scopes stay refused cannot splice
  in a full build (ADR-027 honest limitation). The grade splices at
  driver level: drop the Fortran `.o`, link exactly one `--whole-archive`
  provider.
