# lane-fhutiter evidence (1 graded, 16 recorded skips/refusals)

Reference: ElmerCSC/elmerfem @ `9f6af2f85`, `huti_aux.F90` md5
`75fde042f30d996b957313fafe1f4061`.
Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)` (asserted in the
grade driver) + gfortran 13.3.0 + cmake 3.28.3 (Makefiles).
Grade driver: `grade/grade_hutiaux.sh` (lane-owned; regenerates
`grade/probe_{ref,port}.{out,err,rc,time}` on demand — generated logs are
not committed); reference logs in the rustsmith T010 evidence dir
(`proof/gate-verdicts.json`, `proof/selection.json`, `proof/grade-logs/`).

## Graded green (1 unit)

| unit | crate | true symbols | differential |
|---|---|---|---|
| fortran:fhutiter/src/huti_aux.F90 | hutiaux | 12 `__huti_aux_MOD_*` | byte-identical stdout+stderr+exit over the corpus below |

Corpus (driver `grade/probe_hutiaux.F90`, prints `E25.17` explicit):
pcondfun copy s/d/c/z (3–4 elems); LU solve + factor s/d 4x4, c/z 3x3
(17-sig-figure round-trip print); tiny-positive and tiny-negative pivots
s (pins the list-directed field incl. sign), tiny pivot d/c/z (pins the
3-digit exponent expansion + complex 10/18-sig parts + field widths);
`n = 1` and `n = 0` edges; randvec s/d/c/z shape checks (`[0,1)` range,
successive calls differ, complex parts pure-imaginary per the observed
upstream shape — stream identity excluded by construction: the orig
`random_number` intrinsic is process-seeded and non-reproducible run to
run, so no port can match its stream).

## Findings the corpus forced (all fixed + re-proven here)

- gfortran complex division is the inline Smith algorithm, not the naive
  `(ac+bd)/d` formula (decoded from `-O0` asm; naive diverges by ulps,
  e.g. `0x3f1be41c` vs naive). The port reproduces Smith bit-exactly.
- gfortran list-directed complex parts print one more digit than scalars
  (10/18 vs 9/17 sig) with 2/3-digit exponents; scalar fields are
  width 17/26, complex fields 36/54. The port emulates the E-range + zero
  rules actually reachable from the pivot branch (all pivots satisfy
  `|pivot| < 1e-16`); the F-range switch is documented, not emulated.
- Dual-language buffering: libgfortran buffers probe writes while Rust
  `println!` flushes per line, so interleaving is meaningless unless the
  probe runs unbuffered — the driver exports `GFORTRAN_UNBUFFERED_ALL=Y`
  on both sides (harness design, not a port change).

## Refined merge gate (runbook §2)

- Port defines all 12 unit symbols as `T`; claims no sibling symbol
  (`comm` of defined-name sets vs the reference `huti_aux.F90.o` = exactly
  the 12 replaced symbols) → **zero STRONG port-level dups**.
- Port-only `T` symbols are toolchain classes (compiler builtins
  `__*vdi2`-family, `rust_eh_personality`, `__rust_*` alloc shims,
  monomorphized std/gimli internals — the archive bundles 335 members
  like every lane precedent) — resolved by the one-provider /
  byte-identity rule at merge, never asserted away. User namespace is
  otherwise clean (gate checks the full list, not just the unit prefix).
- Substitution is driver-level by recorded necessity: the tool's
  `substitute` halt (`no built object`) misfires on Makefiles Fortran
  targets (DependInfo `CHECK_Fortran` format gap, LEDGER.md note), so the
  driver links `probe.o` + exactly one `--whole-archive` `libhutiaux.a`
  provider (plus system `libgfortran`); the `ld -r` composition the tool
  would emit is identical. No full-tree build was modified, so the
  campaign backstop is unaffected by construction.

## Oracle / backstop position (honest)

- The 17-unit affected-set computation (shipped `affected_ctest_tests`
  over 963 enumerated tests, all with texts) returns **0 affected tests
  for every unit**: no quick-test input references `fhutiter` or any
  `huti_*` stem (verified by independent repo grep). The shipped path
  therefore degrades to fallback-full for the whole slice — never a
  vacuous pass — by construction.
- ctest execution is environmentally blocked in this sandbox (MPI_Init
  hangs even for a trivial singleton; 8 independent reproductions). The
  K=3 anchor is referenced from the frozen records (T002 baseline.md:
  479/482 exactly K=3; c9 re-verified on the pinned stack), not
  re-measured here; per-unit scoped-vs-fallback paths above are pure
  selection computation, no test executed.

## Timing (32-core host)

- Scratch configure (WITH_MPI=OFF, ccache): 6.51 s wall, peak RSS 224 MB.
- Scratch full build: 32.39 s wall, peak RSS 211 MB.
- Port crate rebuild: 0.26 s wall, peak RSS 120 MB.
- Probe runs: ref 0.00 s / 2.6 MB, port 0.00 s / 3.1 MB (differential
  corpus is unit-scale; no showcase bench — `bench_sample_unit` unset at
  dispatch and no fhutiter unit is a bench unit).
