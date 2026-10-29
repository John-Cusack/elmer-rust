# T021-B evidence (solver-entry plugins M–Z, lane rs-elmer-t021b)

Reference: ElmerCSC/elmerfem @ fork `9f6af2f85` (campaign rev; backstop
tree is a fresh `-DWITH_MPI` `Debug` build of that rev — the campaign
configure). Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)`
+ gfortran 13.3.0 (asserted in every proof script).

## Method (T020 Poisson boundary, × 56)

Each unit: whole plugin `.so` → Rust cdylib exporting the reference
`.so`'s exact dynamic symbol set (`nm -D --defined-only` set equality:
`awk '{print $2, $3}' | sort` diff empty). Derived dummies cross as
opaque handles; every component access and Elmer call goes through
generated-style Fortran `BIND(C)` accessors compiled against the built
`.mod`s (no pinned offsets); SAVE semantics replicated; no Rust
threads. Differentials: ≥2 distinct clean-room `.sif` configs exercising
the solver (byte-identical stdout modulo wall-clock / stderr /
`scalars.dat`) + repeated-invocation run + branch-liveness control.

## Graded green (56 units)

See `T021B-LEDGER.md` for the per-unit table (entries, LOC, walls).
Every proof re-ran GREEN by the lane runner against the serial tree
(`ELMER_BUILD=work/build-serial`): `nm` equality + differentials +
rerun + liveness, `cargo +nightly build --release` clean,
`cargo +nightly clippy` warning-free.

## Campaign backstop

Fresh tree @ `9f6af2f85`, campaign configure, all 56 T021-B ports
rebuilt against the campaign `.mod`s and spliced over their reference
`.so` files (entry-presence checked per splice; Poisson.F90 left
pristine — T020 evidence stands), full rebuild, `ctest -L quick`:
**479/482, failures exactly the K=3 upstream segfaults #80/#536/#537**
(`ConstantBCTemperature`, `ProfileBCTemperature`,
`ProfileBCTemperatureRobin`). Log: lane `proof/LastTest-t021b.log`.

## Refused with reasons (17)

14 PluginEntry-gate refusals (module/character/assumed-shape/USF motifs;
full reasons in lane `scratch/regrade-mz.json`) + 3 triage refusals:
ParallelRemesh (no-ParMMG serial tree → unconditional FATAL; USF-shape
co-entry), SimpleCircuits (dead legacy module, zero callers, NaN-halt),
ProjectToPlane (2 of 3 entries need real MPI; zero exercising configs).

## Known limitations (disclosed, out of proof scope)

- ParticleDynamics: config B + MaxIter=2 diverges ONLY on ComputeChange
  NRM summary lines (all physical outputs identical; reference
  self-consistent orig-vs-orig). Config A proves 2-iteration re-entry.
- Probes: no-probe-names input segfaults the reference itself
  (never-allocated handle sized unconditionally); the port yields zero
  probes instead. No valid differential exists for a crashing input.
- Several ports contain upstream quirks bug-compatibly or contain
  undefined reference paths on defined inputs only (see lane UNITS.md
  notes and port code comments).
