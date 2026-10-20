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

## Tight-loop recraft (2026-10-20): `mtr_sum` — 4.8x–19.7x

The craft pattern (rustsmith ADR-050) applied to the ported unit
`c:ElmerGUI/matc/src/matrix.c` (crate `gui_c_matrix`, symbol `mtr_sum`),
measured with the repo's A/B bench (`bench/`).

**Anti-pattern D (per-item lookup a bulk transform replaces).** The C
original walked the matrix **one column at a time**, so every element was
read with a stride of `ncola` doubles — for a 4096x4096 matrix that is
one cache line per element, and the whole 128 MiB is re-read per column
pass. The recraft makes **one sequential pass** over the row-major data
and accumulates the columns as it goes. The accumulation order per column
is unchanged (rows ascending), so the sums are the same values added in
the same order: bit-identical results, same digest on both arms.

| workload | before | after | speedup | verdict | peak RSS before -> after |
|---|---|---|---|---|---|
| column sums n=256 | 2.079 ms | 0.437 ms | **4.76x** | faster | 2268 KiB -> 2264 KiB (-0.2%) |
| column sums n=1024 | 75.509 ms | 3.830 ms | **19.72x** | faster | 9908 KiB -> 9928 KiB (+0.2%) |
| column sums n=4096 | 200.674 ms | 16.126 ms | **12.44x** | faster | 132164 KiB -> 132240 KiB (+0.1%) |

11 reps, arms interleaved, one corpus per size, peak RSS through `wait4`,
digest identical on both arms. Machine rows: `bench/mtr-sum-n256/`,
`bench/mtr-sum-n1024/`, `bench/mtr-sum-n4096/`.

The win grows with the matrix because the original is memory-bound on the
strided read; the honest floor is the smallest size measured (4.76x at
n=256, where the matrix already fits in cache). Peak RSS is unchanged
(the corpus dominates the child's footprint on both sides).

The bench for this unit is `crates/elmer-gui-c/matrix/benches/matrix_bench.rs`
(the matc temp allocator is implemented in the bench, because the measured
result is read back from it).
