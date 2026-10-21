# lane-fem-b evidence (round 2, T014 relaxation)

Reference: ElmerCSC/elmerfem @ `9f6af2f85`, fork `c0b4bc6` (seed from
original), `HarmonicSignals.F90` md5 `dd26ac5dd5fd96ed7809666bc1da2ce9`.
Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)` (asserted in the
grade driver) + gfortran 13.3.0. Grade driver:
`grade/grade_harmutils.sh` (lane-owned; log `grade_harmutils.log`);
machine rows in `store.db` run `fem-b2-harmutils-grade` + skip rows in
run `fem-b2-regrade` (218 scope units total).

## Graded green (1 scope)

| scope unit | crate | true symbols | differential |
|---|---|---|---|
| fortran:fem/src/modules/HarmonicSignals.F90#harmutils | harmutils | `__harmutils_MOD_sinsum`, `__harmutils_MOD_cossum` | 11/11 byte-identical stdout, stderr empty both sides |

Corpus (driver `probe_harm.F90`, print `ES25.17E3`): whole-array SinSum /
CosSum; single-row section `amp(2:2,:)`; offset `amp(2:4,:)`; negative
stride `amp(4:1:-2,:)`; `t = phase = 0` edge on 8×2; strided `big(1:7:2,:)`;
1×2 unit array both kernels; **held-out 10**: negative-stride
`big(8:1:-3,:)` with negative `t`; **held-out 11**: zero-extent section
`amp(3:2,:)` (empty sum = 0). Divergence 0 at 17 significant figures.

## Refined merge gate (runbook §2)

- Port defines both unit symbols as `T`; claims no sibling symbol
  (`source_` absent); `comm` of defined-name sets vs the reference
  `HarmonicSignals.F90.o` = exactly the two replaced symbols → **zero
  STRONG port-level dups**.
- Port-only `T` symbols are `__*_di2`-class compiler builtins +
  `rust_eh_personality` (same class as T003) — resolved by the
  one-provider / byte-identity rule at merge, never asserted away.
- Substitution is driver-level by recorded necessity: the splice
  replaces the whole `.o`, and sibling scope `#source` (derived-type
  dummies, Elmer world) stays refused — a full-build splice would fail
  the link on the missing sibling, never silently (ADR-027 limitation).
  The grade links `probe.o` + exactly one `--whole-archive`
  `libharmutils.a` provider (plus system `libgfortran`/`libm`).
- Every splice target `test -f`-checked in the driver (`.o` × 2,
  `probe_orig`, `libharmutils.a`, `probe_ports`).

## Oracle / backstop position (honest)

- No full-tree build was modified (driver-level grade), so the
  campaign backstop (479/482, K=3) is unaffected by construction;
  re-running it here would measure the pristine tree, not the port.
- `#source` (the file's other scope) calls `SinSum` internally; its
  oracle coverage (circuits/hysteresis cases) exercises the Fortran
  kernel, not the Rust port — scoped parity for the port is the
  differential above, which covers the same call shape (whole arrays
  and sections) the solver uses.

## Timing (ADR-045)

- Grade walls (32-core host): pin 0 s, gate 0 s, orig build+ref 1 s,
  port build+grade 0 s, symbol gate 0 s.
- `bench_sample_unit` = `fortran:fem/src/modules/HeatSolve.F90`
  (round-1 dispatch: HeatSolve in 357 case.sif files). HeatSolve stays
  refused (external scope), so **no timing sample runs this lane** —
  samples run ONLY on the bench unit, never per unit.

## Round-2 verdict

- Scopes graded green: **1** (`#harmutils`).
- Scopes recorded skips: **217** — 179 external-subroutine (F77-global,
  non-module) scopes, 33 module scopes with derived-type dummies
  (overlap: 19 also character, 19 undeclared/implicit, 6 optional,
  3 allocatable/pointer, 2 rank>2, 1 non-REAL(8) shape), 2 further
  explicit-shape scopes (`#vecutils`, `#particledynamicsstuff`,
  `#dxfile` — one overlaps derived), 3 vacuous admits (empty
  `EMWaveSolverUtils`; data-only `GlobMat`/`GlobMatC` owning the
  unproven `B`-class `__globmat_MOD_matrix` / `__globmatc_MOD_matrix`
  descriptor symbols). F77 COMMON/BLOCK DATA files: 0.

## T017 shape B — module-dummy port (`fortran:fem/src/modules/CoordinateTransform.F90#vecutils`)

Tool PR: rustsmith#108 (gate relaxation, ADR-050). Reference: Elmer fork
`c0b4bc6` (`fem/src/modules/CoordinateTransform.F90` md5
`348989c19c13fb198cd7f9c422d036ef`). Toolchain: `rustc
1.100.0-nightly (feaadeeac 2026-09-19)` + gfortran 13.3.0 (the calling
convention is version-specific — ADR-027 §1/§3). Proof script: rustsmith
`tests/module_dummy_shape_proof.sh`; machine rows in the grading store
(run `t017-module-dummy-grade`).

| unit | crate | true symbols | differential |
|---|---|---|---|
| fortran:fem/src/modules/CoordinateTransform.F90#vecutils | vecutils | `__vecutils_MOD_norm2`, `__vecutils_MOD_normalized`, `__vecutils_MOD_jac`, `__vecutils_MOD_transform2rank` | 7/7 byte-identical stdout (14 values at `ES25.17E3`), stderr empty both sides |

Corpus: 4 vector cases (axis-aligned, Pythagorean, mixed-sign,
degenerate) x 2 kernels, plus 3 matrix cases (identity, sheared,
non-symmetric) x 2 kernels. The vector cases pass a **rank-1 section**
actual — gfortran packs the section and passes the packed pointer, the
same raw-base-address wire — and the matrix cases pass whole 3x3 arrays.

Measured wire (ADR-050): the explicit-shape dummy travels as a raw base
address in Fortran's column-major order, and because these four
functions return **arrays**, the result travels back through a hidden
**first** argument that is a gfortran **array descriptor** (base at
offset 0, dimension stride at offset 40) — not a plain pointer. The
first port draft assumed a pointer and wrote the result into the
descriptor; the corpus caught it (all zeros where the reference printed
the norm).

Refined merge gate (runbook §2): the port defines all four unit symbols
as `T`, the intersection with the reference `CoordinateTransform.F90.o`
is exactly those four — **zero STRONG port-level dups** — and the port
claims no sibling symbol. Honest limitation: the file's two F77 globals
(`rotmsolver_`, `rotmsolver_init_`) stay refused (derived-type dummies),
so the substitution is **driver-level**, not a whole-file splice.

Timing (ADR-045): the lane-fem-b `bench_sample_unit`
(`fortran:fem/src/modules/HeatSolve.F90`) stays refused, so no timing
sample runs for this unit; per-unit grades carry no timing claim.
