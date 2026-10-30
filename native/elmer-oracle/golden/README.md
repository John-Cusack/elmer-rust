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

σ_F is the largest relative-Linf distance (node space, via the Perm) of
any sensitivity arm from Release, per field and saved step; the bar is
`max(10·σ_F, 1e-10)`, with the cap only where every arm agrees bitwise.
1,293 fields agree bitwise; 8 fields in 6 tests differ in structure between
arms and have no Tier-A bar.

Kill/escalate check (§5.2), signed off through firstmate (D11):

- **Raw count: 51/482 (10.6%)** tests have σ_F > 0.1 × their RefNorm
  tolerance when eigenvector fields are compared as raw nodal values.
- **Graded count: 46/482 (9.5%)**, under the 10% re-plan bar, so there is no
  S1/S2 re-plan. The 5 tests that make the difference (StressEigen2,
  StressEigen2_np4, StressEigen2_np8, EM_port_eigen_unit_voltage,
  PlatesEigenComplex) are eigen analyses. An eigenvector is defined only
  up to sign (and, within a degenerate cluster, up to rotation), so a
  legitimate reordering can flip or rotate it and give a raw distance of
  O(1). §4.1 assigns the eigen class its own grading mode: sort eigenvalues,
  normalise signs, compare clusters by principal angles. Raw σ_F on these
  fields therefore measures the basis choice, not Fortran's sensitivity,
  and they are excluded from the count.
- The 46 are classed `tier-a-unfit` or `structure-differs` in `classes.tsv`
  and graded by O1+O4+O6, as §5.2 prescribes.

Most of the excess comes from fields whose golden is numerical noise
(‖u‖∞ near 1e-16, so any reordering is O(1) relative), contact and shell
cases, and CutFEM level-set fields.

## Mutation suite

Noise (1e-6·‖u‖∞), node permutation (swap the extreme values) and one
changed value are applied per field and step; a dropped step and swapped
variables per file. A mutation is detectable when its size in the grader's
own metric exceeds the field's bar. **14,639/14,639 detectable mutations
are caught, 0 escape.**

**Known coverage limit.** Two counts are kept visible, because the suite does
not test them:

- **2,173 below-floor mutations.** The field's own reorder noise (10·σ_F)
  is at least as large as the perturbation, so at Tier A no grader could
  tell that mutation from a legitimate reordering. These fields are
  covered only by O1/O4/O6.
- **1,449 vacuous mutations.** The mutation cannot change the file: for
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
