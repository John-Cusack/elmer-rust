# elmer-fem-a

Lane-fem-a ports: Elmer `fem/` core **excluding** `src/modules/`
(168 units: c 42, cxx 2, fortran 123, python 1).

## Layout

```text
crates/elmer-fem-a/
  ports/<crate>/   # one staticlib crate per ported unit; crate name =
                   # rustsmith scaffold_crate_name(unit stem), lowercased
  LEDGER.md        # all 168 units: ported / refused / skipped with reasons
  EVIDENCE.md      # graded-run evidence (corpus, gates, backstop, timings)
```

## Rules (campaign standing pins)

- Toolchain: pinned nightly `1.100.0-nightly feaadeeac 2026-09-19`
  (`rust-toolchain.toml` + hash assert in the grade driver).
- Every port exports the unit's **true object symbols** (from
  `nm -g --defined-only` on the reference build), not the recon module
  name — see EVIDENCE.md per-unit symbol tables.
- Fortran: only `BIND(C)`-clean units port today; `COMMON`-block files
  are recorded skips (ADR-027); other non-`BIND(C)` units are refused by
  `check_substitutable` (no T014 relaxation landed) and recorded.
- Headers / code generators / orphan sources / ext-dep stubs are
  recorded skips with reasons in LEDGER.md, never forced.
