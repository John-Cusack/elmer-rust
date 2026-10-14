# elmer-gui-a — lane-gui-a scope record

Rust-side record for campaign lane `lane-gui-a`: `ElmerGUI/netgen/`
sorted-id half 1, 124 units (c 4, cxx 120), upstream `9f6af2f85`.

## Verdict: 0 graded / 124 recorded skips

No Rust substitution is possible for this lane: every unit's defined
surface is C++-mangled, declaration-only, or dead upstream, and Rust
cannot provide C++-mangled (`_Z*`) symbols.

| reason | units | evidence |
|---|---|---|
| `cxx-mangled-only` | 48 | each `.cpp` compiles; `nm -g --defined-only` shows 0 C-linkage symbols, only `_Z` mangled + toolchain refs |
| `cxx-header-only` | 70 | `.hpp` C++ class/template headers; no translation unit, no linkable surface |
| `c-header-decl-only` | 4 | declaration-only headers (`extern "C"` `Ng_*` defs live in `interface/*.cpp`, owned by lane-gui-b) |
| `dead-upstream-uncompilable` | 2 | `geom2d/spline2d.cpp`, `geom2d/splinegeometry2.cpp`: bare-text file head, absent from the `ng` library `SOURCES`, do not compile |

## Method (batched fixed cost)

One compile+`nm` pass over the lane, not one grade per unit
(interactive units carry the 10-min flat estimate; the whole scope pass
ran in ~47 s wall, per-unit seconds in `scope-run.log` held by the lane
workdir):

```sh
g++ -DNO_PARALLEL_THREADS -IElmerGUI/netgen/libsrc/include -c <unit> -o <unit>.o
nm -g --defined-only <unit>.o   # classify: non-`_` symbols = C-linkage
```

Refined-gate inputs: 984 unique strong (`T/D/B/R`) symbols across the 48
lane objects, zero cross-TU collisions, all `_Z`-mangled; with no Rust
replacement objects there are zero port-level strong dups. No binary
exists for this slice in the oracle build (`WITH_ELMERGUI*=FALSE`), so
differential probes keep the honest `no differential probes` halt, as for
`fhutiter`-class slices.

## Timing sample

`bench_sample_unit` (dispatch choice, recorded here rather than in the
landed `lanes.json`, which no lane edits):
`cxx:ElmerGUI/netgen/libsrc/csg/solid.cpp` — scope-probe wall 0.67 s
(from the shared pass log; quiesced-window showcase bench N/A — the lane
ships no ported code, so no speedup claim exists).

## Layout

```text
crates/elmer-gui-a/
  Cargo.toml          # standalone package (no workspace root touched: merge lane owns it)
  README.md           # this file
  src/lib.rs          # lane constants + reason counts
  scope/units.json    # machine-readable per-unit verdicts (124 rows)
  tests/registry.rs   # registry consistency test (zero-dep JSON subset parser)
```
