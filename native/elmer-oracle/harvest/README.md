# Harvest of the ABI campaign (S0, report §3.6, §5.2 item 11)

Oracle data. The retired ABI campaign proved each solver port with staged
`.sif` differentials: a pristine-vs-port run of two or more configs per
solver, plus branch-liveness controls. This directory turns those cases
into golden-captured oracle cases. It also turns the upstream defects the
campaign documented into known-divergence entries (§4.3 item 8, ADR 001
decision D4).

| Path | Content |
|---|---|
| `cases/<id>/` | One oracle test per staged case: the campaign `.sif` verbatim, its mesh source (ElmerGrid `.grd`, or a link into `cases/_meshes/`), the data files it names, and a generated `CMakeLists.txt` / `runtest.cmake` / `ELMERSOLVER_STARTINFO`. |
| `cases/_meshes/` | Mesh directories (mostly copied from upstream tests by the proofs), stored once by content hash. The runner's staging follows the links. |
| `cases/cases.tsv` | Provenance per case: lane, proof dir, sif name, sha256 of the sif, mesh source, procedures, and why a case was skipped. |
| `index.tsv` | Golden capture per case on the Release arm: verdict, result files, sha256 of its `golden.manifest`, determinism, worst σ_F, and the pristine error for cases that fail. |
| `tolerances.json` | σ_F and Tier-A bar per (case, result file, field, saved step), in the format of `../tolerances.json`. |
| `defects.tsv` | 15 known-divergence entries with their replacement checks. |

`scripts/harvest-cases.py` rebuilds `cases/` from the campaign lanes'
scratch dirs (`<lane>/scratch/proof-<solver>/case*.sif`; the scripts
that wrote them are `crates/elmer-fem-b/ports/<solver>/proof-<solver>.sh`).
It reads those dirs and never writes them. A later lane that re-staged
the same sif bytes adds nothing, so each case appears once.

## Result

258 distinct staged cases from lanes T020, T021-A, T021-B, T022-A and
T023, covering 108 solver proofs:
- 2 skipped: MagnetoDynamics2D's mesh `1381` is not beside its sifs.
- **236 golden-captured and passing** on the Release, FMA and OpenBLAS
  arms (`scripts/sweep-capture.sh` with `REGISTRY=` and `--exit-verdict`).
  The verdicts are identical on all three arms.
- 104 of the 108 solver proofs have at least one captured case, and 102
  have two or more.
- **Determinism:** a second Release capture agrees bitwise on all 212
  cases that write fields. The other 24 write no field variables (scalars
  or files only).
- **σ_F** (worst of FMA and OpenBLAS against Release, relative Linf in
  node space) over the 212 cases with fields:

  | σ_F | Cases |
  |---|---:|
  | 0 | 12 |
  | ≤ 1e-12 | 147 |
  | ≤ 1e-10 | 7 |
  | ≤ 1e-6 | 19 |
  | > 1e-6 | 26 |
  | structure differs (`inf`) | 1 |

  Most of the 27 cases above 1e-6 are eigen or ill-conditioned staging.
  CraigBamptonSolver's eigenvector signs give 2.0, and
  CoordinateTransform's structure differs across arms. Their golden data
  is reproducible (bitwise on Release), but they are Tier-A-unfit:
  grade them O1 or by class, as `../golden/classes.tsv` does for the
  quick suite.

**20 fail on the pristine arm, the same way on all four runs:**

| Cases | Why |
|---|---|
| elementstats-caseC, meshchecksum-caseC, reynoldssolver-caseC | Negative or liveness controls that fail by design (a Fatal threshold, an altered reference value). |
| eliminateperiodic-caseA/A2/B/C | The arm does not build `EliminatePeriodic.so` (linsolve hook). |
| cs-case_12…15 | ConsistentSplitting cases that stop in `AddEquationBasics` on the pristine arm (the variable exists but no equation owns it). |
| nodetoedgesolver-caseC/D, partitionmesh-caseB | They read the output of an earlier proof run (a restart file, a partitioning). A single-run case cannot express that. |
| gmshoutputreader-caseA/B | The proof put the Gmsh input inside the mesh directory that ElmerGrid generates. |
| elmer2openfoamio-caseA/C | The proof created the OpenFOAM output directory before the run. |
| reloaddata-caseA/B | `LoadInputFile` rejects the reloaded `flow.dat` on this build. |

The ElmerIce USF proofs (T025) are not harvested. Their sifs do not follow
the `case*.sif` convention, and they need the ice tree's libraries.

## Golden data

As for the quick suite, the golden files stay outside the repo, under
`/home/john/elmer-laptop-lanes/elmer-lane-oracle/harvest/golden-<arm>/<id>/`
on the oracle host, each with its `golden.manifest`. `index.tsv` pins
every Release manifest by hash.

## Known divergences (`defects.tsv`)

Each entry names the unit, what triggers the defect, what pristine Elmer
does, and the replacement check. Pristine golden data never grades a
candidate on a trigger path. The defined path is graded against the
harvested golden cases of the same unit. On the trigger the candidate
either refuses (a capability error) or follows the author's evident
intent. In the second case it is graded candidate-vs-candidate (a rerun
must agree bitwise), never against Fortran. Sources are the campaign's
ADR-058/059 and the proof-script headers.

## Reproduce

```sh
D=<campaign data root>       # holds rs-elmer-t02*/scratch/proof-*
B=target/release/elmer-oracle
python3 -I scripts/harvest-cases.py $D harvest/cases
awk -F'\t' 'NR>1 && $8=="" {print $1"\t"$1"\t1"}' harvest/cases/cases.tsv > reg.tsv
REGISTRY=reg.tsv scripts/sweep-capture.sh $B <build-release> harvest/cases out-release 4 harvest --exit-verdict
$B sigma-table harvest/cases sigma.json out-release fma=out-fma openblas=out-openblas
```
