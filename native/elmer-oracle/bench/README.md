# S0 benchmark baselines (strategy report 4.5)

These are the reference timings of the Fortran Elmer that later speed claims
are compared against. Report 4.5 sets the rules:
- Release Fortran with OpenMP, and Release with MPI, on the same pinned BLAS;
- 4–6 scaled cases;
- the median of at least 5 runs on a quiesced, pinned host, with RSS.

## Cases

Every mesh is a structured unit box or square built by the arm's `ElmerGrid`
from `cases/*.grd.in`.

| Case | Solver | Full mesh | DOF (full) | Linear solver |
|---|---|---|---:|---|
| `heat3d` | HeatSolve, steady, unit source, one side at 0 | 100³ hexahedra | 1,030,301 | CG + ILU0 |
| `elast3d` | StressSolve, self-weight, one side clamped | 60³ hexahedra | 680,943 | BiCGStab(4) + ILU0 |
| `poisson2d` | HeatSolve 2-D, boundary at 0 | 2000² quads | 4,004,001 | CG + ILU0 |
| `ns2d` | FlowSolve lid cavity, Re 100, Newton after 3 | 300² quads | 271,803 | BiCGStab(4) + ILU1 |
| `eigen3d` | StressSolve, 10 lowest modes (ARPACK) | 30³ hexahedra | 89,373 | BiCGStab(4) + ILU0 |

`--scale small` runs the same cases on small meshes as a smoke test (each
under a second).

## Running

The arm is the oracle's `release` build in `../oracle-builds.json`:
- Release, OpenMP and Open MPI 4.1.6, reference netlib BLAS/LAPACK;
- used read-only;
- no Fortran is built here.

```
A=/home/john/firstmate/data/elmer-s0-oracle/build-release
bench/run.sh --arm $A --out OUT --scale full --runs 5 --mode omp --threads 8 --cpus 0-7
bench/run.sh --arm $A --out OUT --scale full --runs 5 --mode mpi --np 8 --cpus 0-7
```

**Pinning:**
- **OpenMP:** `taskset -c CPUS`, with `OMP_PROC_BIND=true` and `OMP_PLACES=cores`.
- **MPI:** `mpirun --cpu-set CPUS --bind-to core`, on an `ElmerGrid -partition 2 np/2 1` mesh.
- `HWLOC_COMPONENTS=-gl` is set, because hwloc's GL probe hangs without an X server.

**Output:**

| File | Columns |
|---|---|
| `OUT/runs.tsv` (one row per run) | case, mode, scale, workers, run, DOF, wall time, ElmerSolver's own `SOLVER TOTAL TIME` (real), max RSS |
| `OUT/summary.tsv` | the median of each |

**Wall time and RSS:**
- **OpenMP:** wall time is `/usr/bin/time -v`'s.
- **MPI:** wall time is the launcher's, and RSS is the sum over the ranks.
- **Check:** a run that does not reach `ELMER SOLVER FINISHED` stops the script.

## State

- **Smoke test** (`--scale small`, one run, 2 threads or 2 ranks): passes for all five cases in both modes on the super server.
  - Heat, elasticity and Poisson converge.
  - The cavity's Newton iteration converges in 7 nonlinear iterations.
  - ARPACK returns the 10 modes.
- **Full baselines:** taken 2026-11-21 (01:01Z to 03:55Z), see below.

## Results: S0 baseline, 2026-11-21

Raw data: `results/2026-11-21/{omp,mpi}/{runs,summary}.tsv`.

- **Host:** the super server, 16 cores / 32 threads, 62 GB RAM.
- **Arm:** the oracle's Release+OpenMP build, used read-only, from branch commit `3932b36`.
- **Settings:** full scale, 5 runs per case; the tables show medians.
- **Pinning:** both modes ran on cpus 0-7. OpenMP used 8 threads (`OMP_PROC_BIND=true`, `OMP_PLACES=cores`); MPI used 8 ranks (`--bind-to core`) on an 8-part mesh.
- **Order:** the OpenMP set ran first, then the MPI set.
- **Load:** 1-minute load average 4-6 at the start and 14.06 at the end (5- and 15-minute averages 11.59 and 11.18). Other users' jobs were running, so some timings, mostly the later MPI ones, may be inflated. These are a first baseline, not a quiesced measurement.

| Case | DOF | OpenMP wall (s) | OpenMP max RSS (MB) | MPI wall (s) | MPI RSS, summed over ranks (MB) |
|---|---:|---:|---:|---:|---:|
| heat3d | 1,030,301 | 14.6 | 1,293 | 30.3 | 2,166 |
| elast3d | 680,943 | 64.3 | 1,289 | 99.0 | 2,617 |
| poisson2d | 4,004,001 | 129.1 | 3,126 | 272.4 | 4,711 |
| ns2d | 271,803 | 53.1 | 334 | 252.6 | 746 |
| eigen3d | 89,373 | 168.6 | 292 | 1003.2 | 780 |

Elmer's own reported real time agrees with wall time to within 0.3 s in every row.

At 8 workers, MPI is slower than OpenMP on every case: 1.5x for elast3d, 2.1x for heat3d and poisson2d, 4.8x for ns2d and 6.0x for eigen3d. I have not looked into why. Rising host load accounts for some of it. Other likely causes are more linear iterations with rank-local ILU, and the eigen and Navier-Stokes cases depending more on communication. A per-case iteration count would settle it in a later pass.
