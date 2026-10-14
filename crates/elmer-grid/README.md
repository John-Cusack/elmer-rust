# elmer-grid: lane-grid-io ports (`elmergrid/` own C units)

Rust transcriptions of Elmer's own `elmergrid/src/*.c` (8 units) at upstream
`9f6af2f85`, workspace layout: `eg-prelude` rlib + one `staticlib` crate per
unit. Toolchain: pinned nightly `1.100.0-nightly feaadeeac` (campaign pin).

## Units (lane `lane-grid-io`, 8 of 166)

| unit | crate | splice |
|---|---|---|
| `egutils.c` | `egutils` | provider, whole-archive, both targets (config-clean) |
| `egexport.c` | `egexport` | exe-only minimal |
| `fempre.c` | `fempre` | exe-only minimal |
| `egextra.c` | `egextra` | exe-only minimal |
| `egparallel.c` | `egparallel` | exe-only minimal (`USE_METIS` blocks are exe-only) |
| `egnative.c` | `egnative` | both targets minimal (config-clean under both) |
| `egconvert.c` | `egconvert` | both targets minimal (zero live config conditionals) |
| `egmesh.c` | `egmesh` | both targets minimal (zero live config conditionals) |

Minimal splice = unit object + alloc-shim members + one whole-archive
`egutils` provider per binary; `test -f` every splice target.

## The other 158 lane units (recorded skips, never forced)

- `metis-5.1.0/` (109): vendored link target (ADR-031), not port scope.
- Headers (13): no compilable claim.
- `fhutiter/` (13 Fortran + `dummymain.c`): honest `no differential probes`
  halt — no binary exists; non-`BIND(C)` shapes unproven under the ADR-027
  shim. No `COMMON` blocks anywhere in the lane (verified by grep).
- `misc/` (19) + `ElmerWorkflows/` (4): not in the campaign build (no
  target, no binary, no tests).

Machine rows: lane `store.db` (`run_id = lane-grid-io`): 166 unit rows
(8 `passed`, 158 `parked` with reasons) + 40 gate rows.

## Proof (this lane)

- Equivalence corpus 27/27 byte-identical (exit + stdout + stderr +
  produced tree): usage, error/exit paths, all oracle conversion inputs
  (grd/unv/comsol/fvcom/gmsh/clone/extrude/`.eg`-script), 4 export formats,
  mirror, cartesian/METIS/optim/halo partitions. Build-stamp date and METIS
  timing block normalized (build identity, not behavior).
- Refined gate: zero STRONG port-level dups across all 8 archives
  (`nm -g --defined-only`, all types, weak + toolchain classes filtered);
  weak libm dups permanent (same-toolchain byte-identity).
- Oracle: `ctest -R ElmerGrid` 10/10 on the spliced tree.
- `cargo build` (dev + release) clean, `cargo clippy -- -D warnings` 0 errors.
- `cargo test`: one pre-existing harness-only failure (`fempre`
  entry-symbol `main`, also present pre-lane); no test fns otherwise.
- No showcase bench: `bench_sample_unit` unset at dispatch (ADR-045 §4).

## Lane errata vs prior slices

Corpus caught two transcription bugs in the migrated slices (fixed here,
both predate this lane; proven by the green corpus above):

- `egmesh` parent-search: the `body*`→`bnry*` rename was merged with the
  allocation guard (`&&`), routing existing names to the exceeds-warning
  (UNV `box.unv`: `body3` vs `bnry3` + warning flood). Restored the C nest.
- `egnative` `SUBCELL SIZES 1`: the X-direction info print was omitted
  behind a false NOTE; C prints it. Restored.
