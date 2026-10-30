# Golden capture of the quick suite (S0, report §5.2 items 2–5)

Oracle data. Changes only by ADR with captain sign-off.

## What is here

| File | Content |
|---|---|
| `quick-release.tsv` | One row per registered quick test (482): verdict, runner mode (`plan` or `replay`), the result files captured, and the sha256 of its golden manifest. |
| `classes.tsv` | Grading class per test: eigen, adaptive, Tier-A-unfit, structure-differs, no-field-output, lumped-matrix, test-local plugin (D12), MPI, replay. |
| `mutation-summary.txt` | Mutation suite outcome over every captured result file. |
| `o6.txt` | O6 restart round trip (`scripts/o6-proof.sh`): 5/5 restart tests bit-stable. |
| `../tolerances.json` | Per (test, result file, field, saved step) σ_F and Tier-A bar. |

The golden files themselves (207 MB per arm) stay outside the repo, under
`/home/john/firstmate/data/elmer-s0-oracle-opus/golden/<arm>/<test>/`, each
with its `golden.manifest`; the index pins every manifest by hash.

## Denominator

The arm's own ctest registration, read from its generated
`CTestTestfile.cmake` files (`elmer-oracle registry <build> --label quick`):
482 tests, 454 serial plus 28 `_npN` MPI variants. Grepping test
directories for `quick` is wrong: it also finds tests behind
`IF(HAVE_LUA)`/MUMPS guards that this build never registers.

## Result

- **482/482 captured and passing** on each of the Release (golden), FMA
  and OpenBLAS arms, with `Output File` + `Binary Output` injected (O2).
  469 run through the declarative plan; 13 whose `runtest.cmake` carries
  custom CMake logic are replayed through `cmake -P` with the arguments
  ctest records, remapped from the container to host paths.
- The report's K=3 pristine segfaults do not occur on these arms: no
  known-divergence entry is needed for the quick suite on this build.
- 467 tests carry field variables in their result files (5,665 fields).
  13 write result files without fields (unit tests, ODE solvers, STL
  output) and 2 write a lumped matrix `.dat` instead (`classes.tsv`).

## σ_F and the kill criterion

**Metric.** The error of a field is its Linf difference in node space
(via the Perm), divided by `max(‖g‖∞, scale of its variable)`.
- The variable scale is the largest `‖·‖∞` over the components of the
  solver variable the field belongs to, taken from the result file's
  preamble (`flow solution[velocity:3 pressure:1]` groups `velocity 1..3`
  and `pressure`).
- A scalar variable is its own scale, so its metric is unchanged.
- **Exact zeros stay exact.** A golden field that is identically zero keeps
  the exact-zero rule. Any candidate value beyond one ulp of its variable's
  scale (`ε · scale`) is a mismatch, so the scaled floor never applies to
  it.
- **Why the floor.** Without it, a component that is round-off around zero
  is judged against its own noise. StokesUnitElements' `velocity 2/3` and
  `pressure` are 5e-17 to 1.6e-15 next to `velocity 1` at 0.1, and failed
  on 1e-16 differences. With it, each component is judged against its
  whole variable, as Elmer's own norm and `Reference Norm` judge it.

**Bars.** σ_F is the largest such distance of any sensitivity arm from
Release, per field and saved step. The bar is `max(10·σ_F, 1e-10)`, with
the cap only where every arm agrees bitwise. 1,295 fields agree bitwise; 6
fields in 6 tests differ in structure between arms and have no Tier-A bar.
Both sensitivity arms, graded as candidates, pass all 467 tests that have
fields.

**Kill/escalate check (§5.2), signed off through firstmate (D11).**
- **With the floor: 25/482 (5.2%)** tests have σ_F > 0.1 × their RefNorm
  tolerance once eigen tests are graded by eigen mode; the raw count is
  30/482.
- **Before the floor** (signed off with
  https://github.com/John-Cusack/elmer-rust/pull/212), the count was 46/482
  (9.5%), or 51/482 (10.6%) raw. The 5 eigen tests are excluded because an
  eigenvector is defined only up to sign (and, within a degenerate cluster,
  up to rotation): a legitimate reordering can flip it, which gives O(1)
  raw distance. §4.1 grades them by eigen mode (sort eigenvalues, normalise
  signs, compare clusters by principal angles).
- The 25 are classed `tier-a-unfit` or `structure-differs` in `classes.tsv`
  and graded by O1+O4+O6, as §5.2 prescribes.

## Mutation suite

Noise (1e-6·‖u‖∞), node permutation (swap the extreme values) and one
changed value are applied per field and step; a dropped step and swapped
variables per file. A mutation is detectable when its size in the grader's
own metric (with the floor) exceeds the field's bar. **14,634/14,634
detectable mutations are caught, 0 escape.**

**What the floor changed.** Every mutation was compared with the
pre-floor run. 9 that were caught before are now below the floor, and none
escapes. All 9 are on StokesUnitElements' round-off components in
`unitelems/case.result`. There all three arms agree bitwise, so σ_F = 0 and
the bar is the 1e-10 cap: these mutations sit under the cap, not inside
measured σ_F noise of the same field. Why they are still undetectable:

| Mutation | Size / variable scale | Absolute change | Evidence |
|---|---:|---:|---|
| noise, one value on `velocity 2`, `velocity 3`, `pressure` (6) | 5e-22 to 1.6e-20 | 5e-23 to 1.6e-21 | below one ulp (2.2e-17) of the variable's scale (0.1): no computation at that scale can register it |
| permutation of `velocity 2` | 7.5e-16 | 7.5e-17 | same size as the ~1e-16 differences an independent implementation shows on these fields (#262); below the reorder noise the arms show on `velocity 2` in this test's Taylor-Hood file (σ_F 1.4e-14) |
| permutation of `velocity 3` | 1.3e-15 | 1.3e-16 | as above (Taylor-Hood σ_F 7.2e-15) |
| permutation of `pressure` | 2.4e-14 | 2.4e-15 | below the Taylor-Hood `pressure` σ_F of 1.4e-13 |

The exact-zero rule keeps all 6 Shell_And_Solid mutations on its
identically zero `u loads 1/2/6` components caught. Two ShoeboxStaticShell
fields that differed in structure before now have a bar, so 4 more
mutations are caught.

**Known coverage limit.** Two counts are kept visible, because the suite does
not test them:

- **2,182 below-floor mutations.** The field's own reorder noise (10·σ_F)
  is at least as large as the perturbation, so at Tier A no grader could
  tell that mutation from a legitimate reordering. These fields are
  covered only by O1/O4/O6.
- **1,445 vacuous mutations.** The mutation cannot change the file: for
  example, a permutation of a constant field, or a swap in a step with a
  single distinct field.

## Reproduce

```sh
B=target/release/elmer-oracle; S=/home/john/firstmate/data/elmer-s0-oracle
for arm in release fma openblas; do
  scripts/sweep-capture.sh $B $S/build-$arm $S/src/fem/tests golden-$arm 8 quick
done
$B sigma-table $S/src/fem/tests tolerances.json golden-release fma=golden-fma openblas=golden-openblas
$B mutation-suite $S/src/fem/tests golden-release tolerances.json
scripts/golden-index.py $S/src/fem/tests golden-release $S/build-release tolerances.json golden
```

Runs use `HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1
OMPI_MCA_hwloc_base_binding_policy=none` (set by the sweep script).
