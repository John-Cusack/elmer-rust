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
| `prepare_slot` + `Norm Permutation` | `ComputeNorm` with `Nonlinear System Consistent Norm` (owned dofs only) |
| `run_partitioned`, `MpiEngine`, `elmer-mpi-run` binary | `ElmerSolver_mpi` |

```sh
export HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1
mpirun -np 3 elmer-mpi-run case.sif     # or no argument: ELMERSOLVER_STARTINFO
```

At np=1 the binary reads the mesh DB itself, as `ElmerSolver_mpi` does.
Rank 0 writes `TEST.PASSED_<np>` (`TEST.PASSED` at np=1). The exit codes
are `elmer-cli`'s: 0 solved, 1 failed, 2 bad input, 78 refused.

## Evidence

`tests/partitioned.rs` runs `mpirun -np N elmer-mpi-run` for N = 1..4 on
`tests/data/angle`. That case is AngleMetisLayer (heat, CG without a
preconditioner, `Linear System Symmetric`, consistent norm) without `Mesh
Levels`, on a 300-element `angle.grd`. Its partitions were written by the
quick test's own command, `ElmerGrid 1 2 angle.grd -partdual -metis N 3
-connect 1 2`. Every np must reproduce the pristine serial ElmerSolver norm,
0.16133448050272817 (taken from the full-precision `Output File`), to 1e-12
relative.

Measured on the quick test's own 2700-element mesh, against the pristine
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

- `Mesh Levels` in parallel. New edge and face nodes need global ids and
  sharing, matched across ranks by their sorted parent node ids
  (`SParEdgeNumbering`). This blocks AngleMetisLayer_np2/3/4.
- Halo partitions (`-halo`), `Post File`/`ResultOutputSolver`, `RUN` blocks,
  direct solvers, IterativeMethods.F90 methods (BiCGStab(l), GCR, IDR(s)),
  complex, coupled and block systems, and `Coordinate Mapping` other than
  `1 2 3`.
- The other quick MPI tests at np ≤ 4 are refused by the serial path for
  their own features: passive elements, conforming BCs, the
  AlternatingSource plugin, ResultOutputSolver, ViewFactors and
  RigidMeshMapper.
- Setup and run errors are not agreed between ranks. An error that only
  some ranks hit leaves the others waiting in a collective.
