# Cross-lane interfaces

Every type or trait that crosses a crate boundary lives in `native/elmer-core/`
and has exactly one owner. Lane codes are from the swarm plan
(`/home/john/firstmate/data/elmer-swarm-plan-opus/report.md`, section 3).
S1 is the front-end lane (`elmer-sif`, `elmer-mesh`, `elmer-fe`).

## Owners

- `elmer-core`: L00 until it merges, then firstmate as gatekeeper.
- `elmer-solver` api (`Solver`, `SolverCtx`, `Simulation`): L12.

## Map

| Item (`elmer_core::`) | Module | Implemented by | Consumed by |
|---|---|---|---|
| `Scalar` (`f64`, `Complex64`) | `scalar` | elmer-core | all linear-algebra lanes (L01-L05) |
| `CsrMatrix` + serial `LinearOperator` impl | `csr` | elmer-core (struct); L01 `elmer-sparse` (operations) | L01, L02, L03, L04, L05, L06, L12 |
| `LinearOperator`, `Preconditioner`, `IdentityPrecond` | `linalg` | L01 (ILU family), L02/L03 (wrappers), elmer-core (`IdentityPrecond`) | L02, L03, L05, L12 |
| `IterativeSolver`, `SolveStats`, `SolveStatus` | `linalg` | L02 `elmer-krylov`, L03 `elmer-itermethods` | L05, L12 |
| `DirectSolver` | `linalg` | L04 `elmer-direct` | L05, L12 |
| `LinalgError` | `linalg` | elmer-core | L01-L05, L12 |
| `Communicator`, `SerialComm` | `comm` | elmer-core (`SerialComm`); a future MPI lane | L02, L03, L05, L12 |
| `Value`, `Dependent`, `DependentKind` | `params` | elmer-core; S1 `elmer-sif` converts into it | L07, L08, L11, L12 |
| `Params`, `ParamMap`, `normalize_key` | `params` | elmer-core (`ParamMap`); S1 `elmer-sif` (`impl Params` or `From` into `ParamMap`) | L07, L08, L11, L12, L13, L14 |
| `ElementType`, `ElemId`, `ElementRef`, `MeshView` | `mesh` | S1 `elmer-mesh` | L06, L09, L10, L12, L13, L14, L15 |
| `QuadPoint`, `Quadrature`, `BasisEval`, `ElementBasis` | `fe` | S1 `elmer-fe` (Lagrange), L15 `elmer-fe-hp` (p and edge) | L06, L12, L13, L14 |
| `FieldKind`, `FieldRef`, `FieldSource`, `NO_DOF` | `field` | L12 `elmer-solver` | L09 `elmer-vtu`, L10 `elmer-saveio`, L12 (norms, restart) |
| `LocalSystem`, `GlobalAssembler` | `local` | L06 `elmer-assembly` (`GlobalAssembler`) | L07 `elmer-timeint`, L12, L13, L14 |
| `Logger`, `NullLogger` | `log` | L12 / L11 (real sinks), elmer-core (`NullLogger`) | every lane that reports |

`CsrMatrix` carries a plain serial `y = A x` and `y = A^T x` in `elmer-core`
because the orphan rule lets only this crate implement `LinearOperator` for
it. Faster or parallel kernels belong in `elmer-sparse`, behind a wrapper type
or free functions.

## Drift rules (swarm plan 4.4)

1. **One owner per interface**, as listed above.
2. **Changing `elmer-core`:** open a PR of at most 50 lines that touches only
   `native/elmer-core/`, labelled `core-change`. Firstmate merges these ahead of
   everything else and broadcasts "rebase onto main" to all lanes.
3. **Drift detector:** after each merge, firstmate rebases every
   `origin/fm/elmer-swarm-*` branch onto `origin/main` in one scratch worktree
   it owns, and runs `cargo check -p <crate>`. A lane that fails gets a "rebase
   now" inbox message with the error attached. Lanes never create worktrees.
4. **API diffing:** for `elmer-core` and the `elmer-solver` api, firstmate diffs
   `pub` items between `main` and each `core-change` PR. A removed or renamed
   item keeps a deprecation shim for one merge cycle.
5. **No local copies of shared traits.** A lane that redefines a shared trait
   gets a review reject. Check with:
   `rg -n 'trait (LinearOperator|Preconditioner|MeshView|ElementBasis|Params|Communicator)' native/ | grep -v elmer-core`
6. **Test-only duplicates** (P1 bases under `tests/support/`, mesh-file
   readers) are allowed but must stay under `tests/`. `pub mod support` outside
   `tests/` is a reject.
7. **S1 adapters:** when S1 merges, firstmate enables
   `native/elmer-core/tests/s1_adapters.rs` (landed commented out) and runs
   `cargo check --tests -p elmer-core`. A missing impl becomes one S1 follow-up
   PR, not a change in every consumer.
