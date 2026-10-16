# lane-gui-b grade report — 123/123 skipped, 0 ported

Lane `lane-gui-b` (`ElmerGUI/netgen/`, sorted-id half 2): 123 cxx units.
Upstream pin `9f6af2f85`. Toolchain for proof compiles: g++ 11.5.0;
Rust work pinned to `1.100.0-nightly feaadeeac 2026-09-19`.

## Method (batched, fixed-cost)

1. Compiled all 72 lane `.cpp` TUs with the in-tree netgen flags
   (`-DNO_PARALLEL_THREADS -I libsrc/{include,interface,}`): 69 ok, 3 fail.
   Per-unit walls in `evidence/compile-walls.csv` (total 46.0 s).
2. Ran `nm -g --defined-only` over every compiled object: 11,329 defined
   globals, **0 unmangled non-toolchain symbols** — every export is `_Z*`
   C++-mangled. Rust cannot provide C++-mangled symbols (vtable/exception/
   class layout), so all 69 are recorded skips (`cxx_mangled_only`), never
   forced — exactly the brief's rule 2.
3. Archived the 69 objects (`liblane-guib-native.a`, sha256
   `0f0b18b7…088b1d`) and ran the refined-gate assertion
   (`gate-strong-dups.sh` logic, output in `evidence/gate-output.txt`):
   69 members, `test -f` sweep clean, 0 unmangled candidates,
   **0 STRONG port-level dups** (replacement set empty — zero ports).

## Verdicts

| reason | n | meaning |
|---|---|---|
| `cxx_mangled_only` | 69 | compiles; only `_Z*` globals; no C-linkable surface |
| `header_no_tu` | 51 | `.hpp` helper/header; no compilable claim (spec US2 acc.2) |
| `not_built_in_tree` | 3 | uncompilable in-tree, excluded from upstream `ng` target |

The 3 `not_built_in_tree` units, each verified by direct compile:
- `interface/importsolution.cpp`, `interface/nginterface.cpp`: `fatal error:
  ../visualization/*.hpp` — `libsrc/visualization/` is pruned from the Elmer
  snapshot; neither is in upstream `ng` SOURCES. (nginterface.cpp holds the
  only real `extern "C"` Ng_* decls in the lane, but its definitions need the
  missing kernel headers *and* the full C++ mesh kernel — not a separable
  C-linkable unit.)
- `meshing/findip.cpp`: redefines the `FindInnerPoint` template already
  defined `inline` in `findip.hpp:71`; excluded from `ng` SOURCES (dead TU).

Notable negative result: `interface/nglib.cpp` *looks* like a C API
(`Ng_Init`, `Ng_NewMesh`, …) but `nglib.h` has no `extern "C"` — all its
symbols are mangled. Verified in the object inventory, pinned by
`spot_checks_match_measured_proof`.

## Done-gate mapping (runbook §2)

1. **Equivalence corpus**: n/a — HONEST HALT `no_differential_probes`. No
   binary exists: campaign baseline builds with `WITH_ELMERGUI*=FALSE`, the
   netgen kernel builds lib-only, and the GUI shell needs Qt/X (ADR-045
   dispatch observation: GUI lanes need a headless strategy; none exists).
   Same halt class as the `fhutiter` precedent (spec clarification 3).
2. **Refined merge gate**: green — 0 STRONG port-level dups; native
   inventory digest recorded above; selective extraction n/a (no members
   replaced); minimal splice vacuous (no ports).
3. **Machine rows**: `store.db` run `lane-gui-b-r1` (123 `skipped` unit rows +
   `unit_skip` events with walls, `lane_scope`/`lane_gate`/`lane_halt`
   events, run halt `no_differential_probes`) at
   `/home/john/firstmate/data/rs-elmer-lane-gui-b/store.db` + `events.jsonl`,
   written via the `rustsmith-store` API. Merge lane: roll up with
   `rustsmith combine-stores`.
4. **Timing sample**: per-unit compile walls in every `unit_skip` event and
   `evidence/compile-walls.csv`. No showcase bench: `bench_sample_unit` is
   null at dispatch and a mesher-kernel bench on an all-skip lane would be a
   timing claim over ported-nothing — sample intentionally not run.

## Ownership / merge notes

- This branch touches ONLY `crates/elmer-gui-b/` (lane-exclusive paths).
- Shared `README.md` POC-board row intentionally NOT added: the board section
  does not exist yet in `elmer-rust` and eight lanes appending to one file
  would collide; the merge lane should render the board from `store.db` rows
  via `combine-stores` (Ported-gate rule). This crate's `REPORT.md` +
  `evidence/units.json` carry the measured numbers.
- `bench_sample_unit` left null (dispatch-owned, in rustsmith `lanes.json` —
  untouched; this worktree stays clean per lane isolation).
