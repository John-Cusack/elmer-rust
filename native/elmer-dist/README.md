# elmer-dist

Distributed linear algebra for the Rust-native Elmer (stage S2 np=2 spike).
This crate does what SParIterSolver.F90, SParIterComm.F90 and ParallelUtils.F90
do around one parallel iterative solve. It builds on `elmer-mpi`
(communicator, partition reader), `elmer-sparse` (spmv, ILU) and
`elmer-krylov` (HUTI).

| Rust | Elmer |
|---|---|
| `DofLayout` (owner per dof, copy lists) | `ParallelInfo % NeighbourList`, `GlobalDOFs` |
| `DofLayout::sum_to_owned` | `ExchangeRHSIf` |
| `DofLayout::owned_to_local` | `ExchangeResult` |
| `DofLayout::sum_shared` | `ParallelSumVector` |
| `DistCsr::from_local` | `SplitMatrix` (inside block + interface blocks held by the column owner) |
| `DistCsr::matvec` | `SParMatrixVectorVals` |
| `dist_solve` | `SolveLinearSystem` zero-rhs check, `ScaleLinearSystem` (diagonal), `SolveHutiter`, `BackScaleLinearSystem` |
| `read_rank_system` | the parallel branch of `SaveLinearSystem` |

The Krylov vectors hold owned dofs only. HUTI's dot products and norms then
need only a global sum, and `elmer-krylov` already does that through
`Communicator::sum_f64`. The preconditioners (`diagonal`, `ILUn`, `ILUT`)
factor each rank's inside block, so across ranks they act as block Jacobi,
which is what parallel IterSolve does.

Transports: `SerialComm` (one rank), `ThreadComm` (in-process ranks, used by
the unit-level tests) and `MpiExchange` (MPI, over `elmer_mpi::MpiComm`).

## Evidence

- `tests/threads.rs` covers the 8x8-quad Poisson fixture on in-process ranks:
  `partitioning.1` against `partitioning.2`, with and without halo elements.
  Iteration counts are equal, solutions agree to 1e-8 relative, and the
  matvec equals the summed local operator.
- `tests/mpi.rs` (`harness = false`) relaunches itself with `mpirun -np 1`
  and `-np 2` (with `HWLOC_COMPONENTS=-gl`). CG/none, CG/diagonal and
  BiCGStab/none give the same iteration counts at np=1 and np=2, with
  `max |x1 - x2|` below 1e-16.
- Captured system: with `ELMER_DIST_CAPTURE=<dir>`, the same test solves the
  per-rank system a pristine `ElmerSolver_mpi` saved at np=2 for the
  AngleMetisLayer quick test (CG, no preconditioner, tolerance 1e-8). It
  checks the iteration count and norm against `<dir>/np2.log`. Measured:
  293 iterations and norm 0.165708973898 here, against Elmer's 293 and
  0.16570897 (reference norm 1.65708974E-01). All 15 residual lines printed
  every 20 iterations are identical to Elmer's.

To make the capture, work in a scratch copy, never in the reference tree:

```sh
B=/home/john/runs/elmer-fork/orig_src/build
cp $B/../fem/tests/AngleMetisLayer/{case.sif,angle.grd,ELMERSOLVER_STARTINFO} <dir>/
# add "Linear System Save = True" to Solver 1 in <dir>/case.sif
cd <dir>
export HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1 ELMER_HOME=$B/fem/src ELMER_LIB=$B/fem/src/modules
$B/elmergrid/src/ElmerGrid 1 2 angle.grd -partdual -metis 2 3 -connect 1 2 -nooverwrite
mpirun -np 2 $B/fem/src/ElmerSolver_mpi case.sif > np2.log
ELMER_DIST_CAPTURE=<dir> cargo test -p elmer-dist --test mpi
```

The saved files carry Elmer's `GlobalDOFs` but not the dof owners, so the
test uses `owners_by_lowest_rank`. Any consistent owner choice gives the same
summed operator. Only the floating-point summation order differs from
Elmer's.

## Known gaps

- **Edge and face dofs.** `DofLayout::new` takes any `(global id, owner)` per
  local dof, but nothing produces those ids for edges or faces yet:
  `DofLayout::nodal` covers nodal dofs only. ElmerGrid partition files
  record sharing for nodes only (`part.K.shared`). Elmer numbers edges and
  faces at run time (`SParEdgeNumbering`, `SParFaceNumbering`,
  SParIterComm.F90:750 and 1414): an edge or face is a sharing candidate only
  when all of its nodes are shared, and the candidates are matched across
  neighbours by their sorted global node ids. To support edge/face dofs,
  either the mesh lane or this crate needs:
  1. per-partition edge and face lists keyed by sorted global node ids;
  2. a global numbering for them. A collective dense renumbering would do,
     since a key of two or more node ids does not fit the `u64` global id
     directly;
  3. an owner per shared edge or face. Any consistent rule works, for
     example the lowest holding rank, which `owners_by_lowest_rank` already
     computes once a global id exists.
- Complex systems (`SParCMatrixVector`), the constraint matrix, `Linear
  System Skip Mask`, row equilibration and the BILU, Vanka and multigrid
  preconditioners are refused or not implemented.
- The interface exchange in the matvec starts after the local product, so
  communication and computation do not overlap.
