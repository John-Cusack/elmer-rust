# lane-gui-c: ElmerGUI shell (non-netgen) + GUI-matc + logger + tester

Lane `lane-gui-c` of the Elmer campaign (`data/elmer-lanes/lanes.json` in
rustsmith): 148 units — `ElmerGUI/` 142 + `ElmerGUIlogger/` 3 +
`ElmerGUItester/` 3 (c 107, cxx 41). Upstream checkout `9f6af2f85`.

## Verdicts

| verdict | count | meaning |
|---|---|---|
| `graded-build` | 21 | `ElmerGUI/matc/src/*.c` ported to Rust (this workspace); build + substitution partial-link closure green; honest `no differential probes` halt (below) |
| `skipped-cxx-mangled` | 122 | Qt/C++-only surface (class, `Q_OBJECT`, or C++ TU); Rust cannot provide C++-mangled symbols — recorded skip, never forced |
| `skipped-helper` | 5 | declaration-only headers (`ElmerGUI/matc` `config.h`, `str.h`, `elmer/*.h`); no compilable claim, covered by the shared prelude |

Machine-readable rollup: `units-manifest.json` (unit → verdict/reason/crate).
Grading rows live in the campaign store (`store.db`), rendered via
`combine-stores`; this file is the lane record, not the board.

## Why no differential probes (honest halt)

The oracle baseline builds with `WITH_ELMERGUI*=FALSE` (spec § baseline),
so no lane-gui-c TU participates in any built binary: nothing to substitute
into, no probe binary to diff against. Value delivered is porting-rule
coverage: every C TU builds as a C-ABI staticlib and closes as a drop-in
for its native object (see Gate below). No probe path was invented
(`fhutiter` precedent, T010).

## Ports (`gui_c_*`, pinned nightly `1.100.0 feaadeeac 2026-09-19`)

One crate per C TU. Each `src/lib.rs` = shared T005 prelude (C scalar
aliases, `LIST`/`MATRIX`/`VARIABLE`/graphics types with layout asserts,
libc externs) + `extern` cross-unit imports + transcribed unit body with
`#[no_mangle] pub unsafe extern "C"` exports and per-use `SAFETY` comments.

Provenance: bodies for the 13 code-identical TUs match the T005 reference
ports (rustsmith `fm/rs-elmer-c6` @ `685bdd4`), re-targeted here because the
GUI sources are byte-identical modulo license headers. Eight TUs are the
older GUI variant and were adapted:

- `variable`, `eval`, `oper`: no scalar/wrapper pools, no embedded-data or
  combined-block paths; shared MATRIXes wrapped inline.
- `str`, `files`: `str_p`/`str_pstr` are `str.h` globals (defined by the
  `matc` crate via `MODULE_MATC`), not locals.
- `matc`: defines the `str_p`/`str_pstr` globals; no `setlocale`, no
  `mtc_compile`/`mtc_eval`/`mtc_set_real_array`/`mtc_free_compiled`.
- `parser`: single `doit` (no `doit_compile` split).
- `funcs`: identical modulo brace style (kept as transcribed).

Deliberate fidelity decisions (transcription beats lint):

- `gra_rotate_pip180 = 3.1415926535898 / 180.0` keeps upstream's truncated
  literal (differs from `f64::PI` by ~2 ulp); `approx_constant` allowed
  on that item only.
- Every `unsafe` item carries a `SAFETY` comment (port convention);
  `clippy::missing_safety_doc` allowed crate-wide for that reason.

## Gate (runbook §2 item 2)

`gate-check.py` (campaign side, not shipped) over release archives with
**named `ar x` of unit members only** (never whole archives):

1. Zero STRONG port-level dups across all unit members (`nm -g
   --defined-only`, toolchain classes filtered).
2. Per-unit coverage: every native-defined C global is defined by the
   port; no port-only C globals.
3. `ld -r` closure: all 21 unit members link with only libc/libm
   undefined; each unit also closes as a drop-in against the other 20
   native objects (US4 acceptance adapted to a library lane).

`cargo build --workspace` (debug + release), `cargo test --workspace`,
and `cargo clippy --workspace -- -D warnings` are green on the pinned
nightly. Stable 1.89.0 cannot build the ports (`c_variadic` E0658).

## Merge note

`ElmerGUI/matc` is an older duplicate of top-level `matc/`: 13 TUs are
code-identical, so these archives export the same C globals as the future
`elmer-matc` lane ports. Both providers cannot link into one binary; merge
lane decides (link-instead-of-port re-scope candidate, cf. ADR-045 netgen
note). This lane asserts uniqueness only within its own splice inputs.

## Tight-loop recraft attempt (2026-10-20): `jacobi` — no gain, reverted

The craft pattern (rustsmith ADR-050) was applied to the ported unit
`c:ElmerGUI/matc/src/jacobi.c` (crate `gui_c_jacobi`, symbol
`matc_jacobi`, the generalized Jacobi eigensolver) and measured with the
repo's A/B bench (`bench/`, one corpus, rep-interleaved arms, fresh child
per sample, median + MAD, peak RSS via `wait4`, digest equality).

| workload | before | after | speedup | verdict | peak RSS before -> after |
|---|---|---|---|---|---|
| generalized Jacobi eigen n=128 | 210.291 ms | 208.028 ms | 1.01x | **no-change** | 3096 KiB -> 3080 KiB (-0.5%) |

15 reps; MAD 1.6% (before) / 4.7% (after); digest identical on both arms
(`fb5d76d6feb6dc1c`).

**The recraft was reverted.** The variant kept the value caching (the six
matrix entries an iteration reads are not modified by the update loops, so
they are read once) and hoisted the loop-invariant bounds; it measured
1.01x, inside the 3% noise band, so it is not kept on faith. A second
variant that additionally walked row pointers in the three column-update
loops measured **0.96x** (slower) and was dropped as well — the compiler
unrolls the original shape better.

The port source in this lane is therefore **unchanged**; what the lane
keeps is the bench harness
(`crates/elmer-gui-c/jacobi/benches/jacobi_bench.rs`, machine row under
`bench/jacobi-nogain/`), so the next attempt on this unit is measured
against the same corpus instead of guessed. The no-gain row is the point:
it is the honest record of an attempt that did not pay off.
