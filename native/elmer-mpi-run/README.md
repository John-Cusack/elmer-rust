# elmer-mpi-run

Partitioned runs of the Rust-native Elmer, the native counterpart of
`ElmerSolver_mpi`. It takes a mesh that `ElmerGrid ... -metis N` or
`-partdual -metiskway N` split into `partitioning.N`. Each rank reads
`part.<rank+1>.*` and assembles its own elements with the normal
`elmer-solver` solvers. Where Elmer's parallel code differs from its serial
code, this crate installs `elmer-solver` backends that do the parallel step:

| Rust | Elmer |
|---|---|
| `ParMesh` (partition reader `elmer-mpi`, ghost-first node order) | `LoadMesh` with `ParallelInfo` |
| `ParAssembly::try_apply_dirichlet` | `EnforceDirichletConditions` parallel branch: `CommunicateDirichletBCs`, then rows scaled by the diagonal summed over ranks (`ScaleLinearSystem(ApplyScaling=.FALSE.)`) |
| `linear_registry` → `elmer_dist::dist_solve` | `ParallelIter`: `ScaleLinearSystem` with `ParallelSumVector`, `SplitMatrix`, block-local preconditioner, HUTI on owned dofs |
| `linear_registry` → `itermethod_solve` | the same, with BiCGStab(l), GCR or IDR(s) from `elmer-itermethods` in place of HUTI |
| `levels::split` (`elmer-mesh-levels` per rank, then a collective numbering) | `SplitMeshEqual` with `SParEdgeNumbering` / `SParFaceNumbering` |
| `prepare_slot` + `Norm Permutation` | `ComputeNorm` with `Nonlinear System Consistent Norm` (owned dofs only) |
| `run_partitioned`, `MpiEngine`, `elmer-mpi-run` binary | `ElmerSolver_mpi` |

```sh
export HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1
mpirun -np 3 elmer-mpi-run case.sif     # or no argument: ELMERSOLVER_STARTINFO
```

At np=1 the binary reads the mesh DB itself, as `ElmerSolver_mpi` does.
Rank 0 writes `TEST.PASSED_<np>` (`TEST.PASSED` at np=1). The exit codes
are `elmer-cli`'s: 0 solved, 1 failed, 2 bad input, 78 refused.

## Single Mesh

With `Single Mesh = True`, every rank reads the whole mesh DB, and Elmer
turns its parallel paths off (`Mesh % SingleMesh`). The case therefore runs
on each rank through `elmer-run`'s serial composition, called with
`elmer_run::run_case_with` and a `RunEnv` naming the job:

- rank 0 writes `TEST.PASSED_<np>`;
- VTU output is written per rank as `<base>_<np>np<rank+1>`, and rank 0
  also writes a `.pvtu` that lists every piece;
- the `Partition` global is set, and with `Parallel Slices` each rank gets
  its own slice: `slice`, `slice ratio` and `slice weight`.

Elmer also takes the maximum of the nonlinear change over the slices
(`ComputeChange` when `Number Of Slices > 1`). That is not reproduced. It
only matters when one slice would stop iterating earlier than another.

## Mesh Levels

Each rank splits its own piece with `elmer-mesh-levels`. A new node is
the midpoint of an old edge, the centre of a quadrilateral face, or the
centre of an element. It can be shared only when every parent node is
shared, and only with ranks that hold all of the parents. Such a rank
shares the node exactly when it has the same edge or face, which both
ranks recognise by the sorted global ids of the parent nodes. The lowest
holding rank owns the new node. New global ids continue after the largest
old one: each owner numbers its own new nodes in local order, rank by rank,
and sends the ids to the other holders. Every level repeats this on the
previous level's mesh.

## Quick tests

The named quick MPI tests with np ≤ 4 were staged as ctest stages them:
the test directory is copied and its `runtest.cmake` ElmerGrid command is
run with the pristine ElmerGrid. Each was then run with `mpirun -np N
elmer-mpi-run`.

| Test | Verdict |
|---|---|
| AngleMetisLayer_np2, _np3, _np4 | match (`TEST.PASSED_N` = 1) |
| ParallelSlices_np2 | match (`Single Mesh`, see above) |
| HeatControlExplicit_np3 | error: its AlternatingSource plugin is not ported (the serial path stops the same way) |
| poisson_transient_conforming_anti_np4 | refused: conforming BCs |
| OtherMeshResults_np4 | refused: ResultOutputSolver in parallel; serially also solver-specific meshes |
| passive_np3 | refused: passive elements |
| radiation_np2, _np4 | refused: radiation |
| CapacitanceMatrix3/4_np3, SD_CapacitanceMatrix3_np3 | refused: p-elements |
| CoilSolver1, CoilSolverLoop (np 2, 4), PoissonDG_np3, SunAngle_np3 | refused: procedure not ported |
| ComponentResistance_np4 | refused: Component sections |
| StressEigen2_np4 | refused: legacy solver without Procedure (eigen analysis) |

AngleMetisLayer against the pristine serial norm on the refined mesh, from
the full-precision `Output File` (0.1657089738977438, 11041 nodes):

| np | relative difference |
|---:|---:|
| 1 | 2.4e-14 |
| 2 | 2.7e-14 |
| 3 | 2.6e-14 |
| 4 | 2.5e-14 |

## Evidence

`tests/partitioned.rs` runs `mpirun -np N elmer-mpi-run` for N = 1..4 on
`tests/data/angle`. That case is AngleMetisLayer (heat, CG without a
preconditioner, `Linear System Symmetric`, consistent norm) without `Mesh
Levels`, on a 300-element `angle.grd`. Its partitions were written by the
quick test's own command, `ElmerGrid 1 2 angle.grd -partdual -metis N 3
-connect 1 2`. Every np must reproduce the pristine serial ElmerSolver norm,
0.16133448050272817 (taken from the full-precision `Output File`), to 1e-12
relative. It must do the same with `Mesh Levels` 2 and 3; level 3 splits
the new nodes of level 2 again.

The same fixture is also solved with the IterativeMethods.F90 methods,
each checked against a pristine serial run with that method:

| Method | Preconditioner | np | Tolerance | Measured |
|---|---|---|---:|---:|
| BiCGStab(l), l = 4 | none | 1-4 | 1e-13 | ≤ 1.6e-15 |
| GCR | diagonal | 2, 3 | 1e-13 | 1.7e-16 |
| IDR(s) | none | 2, 3 | 1e-8 | ≤ 1.3e-9 |
| BiCGStab(l), l = 4 | ILU(0) | 2, 4 | 1e-8 | ≤ 1.7e-9 |

Without a preconditioner, or with the diagonal one, the parallel iterates
equal the serial ones up to rounding. IDR(s) draws its random shadow space
per rank. ILU(0) is factored on each rank's own block, so in parallel it
acts as block Jacobi, as Elmer's parallel ILU does. Those two meet the
serial norm only to about the linear tolerance.

`tests/partitioned.rs` calls `mpirun --bind-to none`. With OpenMPI's
default binding, every rank of an np ≤ 2 run, and every singleton's
`orted`, is pinned to the same core. On a shared host that core can be
saturated, and the run then never progresses. Set
`OMPI_MCA_ess_singleton_isolated=1` when running a pristine serial
ElmerSolver next to other MPI jobs.

Measured without `Mesh Levels` on the quick test's own 2700-element mesh, against the pristine
serial norm 0.16480021253114105:

| np | relative difference |
|---:|---:|
| 1 | 1.2e-14 |
| 2 | 1.6e-14 |
| 3 | 1.5e-14 |
| 4 | 1.5e-14 |

## How it fits elmer-solver without changing it

`elmer-solver` has no parallel hooks yet, so this crate works around three
gaps. Each workaround is local to the run:

- **Consistent norm.** `ParMesh` numbers the nodes other ranks own first.
  `prepare_slot` then sets `Norm Permutation = n_ghost + 1`. With that,
  `elmer-solver`'s parallel `ComputeNorm` sums only the owned nodes, so each
  shared node counts once. This only works for one dof per node; a
  consistent norm on a vector variable is refused.
- **Serial scaling.** `elmer-solver` would scale each rank's system by that
  rank's unsummed diagonal. `prepare_slot` moves `Linear System Scaling` to a
  private key and turns it off, and the distributed solve then scales with
  the summed diagonal.
- **Matrix layout.** The solver factories get no context, so the transport
  and dof layouts live in a per-thread `RankContext`. A layout is matched to
  its matrix by the sparsity pattern `ParAssembly::create_matrix` recorded.

Two additive hooks would remove these workarounds:

1. `elmer-solver`: a `Backends` entry that gives row ownership to
   `ComputeChange`/`compute_norm_consistent`, and skips local scaling while a
   parallel linear backend is installed.
2. `elmer-cli`: a way for the driver to ask the engine before it refuses
   native runs under more than one MPI rank. `MpiEngine` already implements
   `Engine`.

## Refused or not done

- `Mesh Levels` on boundaries whose parent element is in another
  partition: `SplitMeshEqual` drops boundary elements without a parent,
  serially as well.
- Halo partitions (`-halo`), `Post File`/`ResultOutputSolver`, `RUN` blocks,
  direct solvers, the stationary methods (Jacobi, Richardson, SGS),
  complex, coupled and block systems, `Coordinate Mapping` other than
  `1 2 3`, and the `Mesh Levels` options `elmer-mesh-levels` refuses.
- The other quick MPI tests at np ≤ 4 are refused by the serial path for
  their own features: passive elements, conforming BCs, the
  AlternatingSource plugin, ResultOutputSolver, ViewFactors and
  RigidMeshMapper.
- Setup and run errors are not agreed between ranks. An error that only
  some ranks hit leaves the others waiting in a collective.
