# elmer-fixtures

Rust ports of the Fortran plugins that Elmer's own tests compile next to a case
(`ADD_ELMERTEST_MODULE` in `fem/tests/*/CMakeLists.txt`). The Rust-native Elmer
loads no shared objects, so a case that says `Procedure "Source" "Source"` needs
that function registered under that name. Strategy decision D12: these cases
stay **excluded from the held-out set**; the fixtures are test infrastructure,
not product code. Porting a callback does not mean its parent solver is
supported (InitDist is a callback of an advection solver this repo does not
implement).

References below are paths under `/home/john/runs/elmer-fork/orig_src/fem/tests/`.

## What is ported

| Fixture | Source | Hook kind | Module |
|---|---|---|---|
| `MySource` | `ConstantParamFunc/MySource.F90` | scalar function, `Variable "time, <const>"` | `udf::my_source`, `my_source_fn` |
| `InitDist` | `AdvDiffFCT/InitDist.f90` | scalar function reading `Model % Nodes` | `udf::init_dist`, `init_dist_fn` |
| `Source` | `1dtests/Source.f90`, `1dtests_MeshLevels/Source.f90` | scalar function with per-instance `OldStep` and previous-solution assertions | `udf::OneDSource` |
| `BulkAssembly`, `BoundaryAssembly` | `PoissonAssembly.f90` of `CoupledPoisson1`..`9`, `BlockPoisson1`, `BlockPoisson2` | element assembly kernels (`Bulk Assembly Procedure`) | `poisson_assembly::PoissonAssembly` |

`registry::install(case, model, &UdfRegistry)` registers a case's scalar
fixtures under the names its `.sif` uses, as **new instances**, so two runs in
one process never share `SAVE`d state. `registry::lookup` builds a single
instance. Keys are `(case, library, entry)` because plugin names repeat across
cases: `Source` alone is five different files (see the inventory).

```rust
let model: Arc<dyn ModelView> = /* live view of the running simulation */;
let udfs = elmer_expr::UdfRegistry::new();
elmer_fixtures::registry::install("1dtests", model, &udfs);
let prop = elmer_expr::compile_with(&value, &udfs, &Default::default())?;
```

### Fidelity notes (each is checked by a test)

- **InitDist threshold.** The Fortran writes `0.3**2` without `_dp`, so
  gfortran folds it in single precision. The pristine `InitDist.so` stores
  `0x3FB70A3D80000000` for it (`objdump -s -j .rodata`), not the double `0.09`.
  `INIT_DIST_RADIUS_SQ` reproduces that value, so a node at distance `0.3`
  is inside. The inequality is strict: a point whose squared distance equals
  the threshold gets `0`.
- **1dtests tolerances.** `1dtests_MeshLevels/Source.f90` differs from
  `1dtests/Source.f90` only in the polynomial tolerance (`1e-8` against
  `1e-10`). `OneDTolerances` carries both and the registry picks by case.
- **Source state.** `OldStep` (implicitly `SAVE`d through `= -1`) lives in the
  `OneDSource` instance. It advances only after a clean check of all nodes, so
  a failed check fails again on the next call, as in the Fortran. A node with
  no `Potential` value, a step outside `1..=5` and a call at an integration
  point are errors, where the Fortran reads out of bounds or returns an unset
  result.
- **`Fatal`.** Cores return `FixtureError::Fatal` with the Fortran text. The
  `NodeFn` adapter panics with it, because `NodeFn` has no error channel and
  Elmer halts there too.
- **PoissonAssembly.** One file serves all nine `CoupledPoisson` cases (they
  differ only in declaration order and unused `USE` lines); cases 2 to 9 differ
  in the solver around it (coupled variables, integral/equality/floating
  constraints, line search, a MATC source `1.0-20*tx` in case 2). The kernels
  keep the Fortran's association order, so results match bit for bit given the
  same basis and quadrature. Absent keywords are zeros (`GetReal`), the
  boundary routine returns without contribution when neither `Flux` nor
  `Coefficient` is defined, and each routine reports `1st time` once per
  instance.

## Tests

```
export CARGO_BUILD_JOBS=1 RUST_TEST_THREADS=1 OMP_NUM_THREADS=1 RAYON_NUM_THREADS=1 \
       RUSTC_WRAPPER=sccache HWLOC_COMPONENTS=-gl
cargo test --manifest-path <worktree>/native/Cargo.toml -p elmer-fixtures
```

`tests/support/` holds test-only stand-ins (ADR 003, rule 6): a P1 triangle and
two-node line basis with quadrature, and a mutable model. They go away when
`elmer-fe` and the solver's live model view land. Reference behaviour was
checked on the pristine `ElmerSolver`: `1dtests` runs steps 1 to 5 through the
compiled `Source.so` and passes (`Norm = 1.40257363E+00`, relative error
`1.9e-12`).

## What an integration needs (not blockers for this crate)

1. A live model view with this exact shape, implemented by the solver side:
   `pub trait ModelView: elmer_expr::FieldLookup + Send + Sync { fn n_nodes(&self) -> usize; }`
   (defined in this crate; `FieldLookup::value("potential", node)` must be
   `None` for nodes outside the variable's permutation, and `time()` must be the
   simulation `Time`).
2. A solver-side hook that calls `PoissonAssembly::{bulk, boundary}` per
   element for `Bulk Assembly Procedure` / `Boundary Assembly Procedure`.
   `elmer-solver` has none yet. The coupled cases (`Coupled Solver`,
   `Constraint n Type`, `Variable n`) also need the block/constraint machinery.
3. The `1dtests` Poisson **solver** is a plugin, not a scalar function; see the
   inventory.

## Findings about the quick list

Quick status here is read from the pristine build's `CTestTestfile.cmake`
labels. The CMake `LABELS` argument of `ADD_ELMER_TEST` alone misses tests that
get `quick` from a separate `ADD_ELMER_LABEL(<test> quick)` line (the
`CoupledPoisson` cases do). The `fem/tests` directories register 455 tests with
the label (the strategy report's 482 counts the whole build).
Of the cases named for this lane:

- `1dtests`, `1dtests_MeshLevels` and `CoupledPoisson1`, `3`..`9` are `quick`.
- `ConstantParamFunc`, `AdvDiffFCT` and `CoupledPoisson2` are **not** `quick`
  (labels `serial` or `serial;transient`). `MySource` and `InitDist` were
  ported as asked; `CoupledPoisson2` shares the ported `PoissonAssembly`.
- `Lua` and `KeywordUnitTestLua` carry the `quick` argument but register no
  test in the pristine build (no Lua support), so they are not counted.
- 66 quick test directories compile a local plugin, 72 plugin files, 48
  distinct file contents (strategy report counts 67).

## Inventory: remaining quick-test local plugins

Hook kinds: *scalar function* is `FUNCTION F(Model, n, t)` called through a
`Procedure` keyword (an `elmer_expr::NodeFn`); *solver* is
`SUBROUTINE S(Model, Solver, dt, Transient)`; *element assembly* is the
11-argument `Bulk/BoundaryAssembly`; *UMAT* is a 37-argument material routine.
"State" lists `SAVE`d and initialised-on-declaration variables, which are
per-process in Fortran and per-instance here. "default_* solve path" means the
plugin calls `DefaultInitialize` / `DefaultSolve` / `DefaultUpdateEquations`
style services, which map to `SolverCtx::default_*`; "elmer-fe" means it needs
`ElementBasis` and `Quadrature` (`ElementInfo`, `GaussPoints`). Generated from the
sources and `.sif` files by a script, then spot-checked by hand; the blocker
column is a guide to what to read first, not a proof.

### Scalar functions (cheapest next batch)

Needs only `ModelView`. `SourceFun` has three bodies: `32(y-y^2+x-x^2)`
(BDM2D, RaviartThomas2D), `0` (MixedPoisson2D/3D/3D_Brick) and
`128((y^2-1/4)(z^2-1/4)+(x^2-1/4)(z^2-1/4)+(x^2-1/4)(y^2-1/4))`
(RaviartThomas3D). `AlternatingSource` takes `x(2) = (temperature, time)`, flips
a sign when the temperature crosses `+-0.07`, and prints on each flip. Its
`ParEnv % MyPe == 0` guard is rank 0 in serial. `TangentForce` needs
`Model % CurrentElement`, the element normal and the `Displacement Loads`
variable.

| Path (under `fem/tests/`) | Hook kind and entries | State | Blockers |
|---|---|---|---|
| `BDM2D/Source.F90` (+1: RaviartThomas2D) | SourceFun:scalar function | none | none (ModelView only) |
| `HeatControlExplicit/AlternatingSource.F90` | AlternatingSource:scalar function | SAVE sgn; init tlim=0.07_dp, sgn=1 | ParEnv/MPI use |
| `MixedPoisson2D/Source.F90` (+2: MixedPoisson3D, MixedPoisson3D_Brick) | SourceFun:scalar function | none | none (ModelView only) |
| `NormalTangentialDisplacement/TangentForce.f90` | TangentForce:scalar function | SAVE Visited, ElementNodes, LoadPerm, LoadDofs, LoadValues; init Visited=.FALSE. | Model%Variables access; current-element context |
| `RaviartThomas3D/Source.F90` | SourceFun:scalar function | none | none (ModelView only) |

### Element assembly

`BlockLinElast` is the vector-elasticity analogue of `PoissonAssembly` and is
used with `Bulk Assembly Procedure 11..22` block numbering.

| Path (under `fem/tests/`) | Hook kind and entries | State | Blockers |
|---|---|---|---|
| `BlockLinElast1/BlockLinElast.f90` (+1: SD_BlockLinElast1) | BulkAssembly:element assembly; BoundaryAssembly:element assembly | SAVE Visited, Source, Nodes, Basis, dBasisdx, Params, Flux, Text, Coeff; init Visited=.FALSE., Visited=.FALSE.; work arrays | elmer-fe; p/edge/face Element in .sif; coupled/constraint/block .sif; current-element context |
| `BlockPoisson1/PoissonAssembly.f90` (+1: BlockPoisson2) | BulkAssembly:element assembly; BoundaryAssembly:element assembly | SAVE Visited, Source, Cond, React, Nodes, Basis, dBasisdx, Flux, Text, Coeff; init Visited=.FALSE., Visited=.FALSE.; work arrays | elmer-fe; coupled/constraint/block .sif; current-element context |

### Whole solvers

These are not scalar fixtures: each is a `Solver` implementation (assembly,
`DefaultSolve`, often `Update` of norms) and belongs with the physics crates or a
second fixtures batch once `SolverCtx` and `elmer-fe`/`elmer-fe-hp` are on
`main`. The Poisson family (`PoissonSolver`, 7 sources, ~12 tests) is the
largest group and the best next target after the scalar batch; `CountSolver` and
`TrescaSolver` are small.

| Path (under `fem/tests/`) | Hook kind and entries | State | Blockers |
|---|---|---|---|
| `1dtests/Poisson.f90` (+1: 1dtests_MeshLevels) | PoissonSolver:solver | SAVE STIFF, LOAD, FORCE, AllocationsDone, Nodes; init ; work arrays | default_* solve path; elmer-fe |
| `AdvReactDG_P/AdvReactDG_P.F90` | AdvectionReactionSolver:solver | SAVE MASS, STIFF, LOAD, FORCE, Velo, MeshVelo, Gamma, &; init ; work arrays | default_* solve path; elmer-fe; boundary parent lookup; edge/face basis (L15); p/edge/face Element in .sif; writes Solver%Variable |
| `BDM2D/ComputeError.F90` | IntegrateError:solver | SAVE SolDOFs, Indices, AllocationsDone; init ; work arrays | elmer-fe; edge/face basis (L15); writes Solver%Variable |
| `CurvedBndryPFEM/Poisson.f90` (+1: SD_CurvedBndryPFEM) | PoissonSolver:solver | SAVE STIFF, LOAD, FORCE, POT, AllocationsDone, Nodes; init ; work arrays | default_* solve path; elmer-fe; boundary parent lookup; p/edge/face Element in .sif; writes Solver%Variable |
| `CurvedBndryPFEM2/Poisson.f90` (+1: SD_CurvedBndryPFEM2) | PoissonSolver:solver | SAVE STIFF, LOAD, FORCE, POT, AllocationsDone, Nodes; init ; work arrays | default_* solve path; elmer-fe; p/edge/face Element in .sif |
| `ElementUnitTest_edge/SingleEdgeElements.F90` (+2: ElementUnitTest_edge_piola, ElementUnitTest_edge_quadratic) | SingleEdgeElements_init:solver _init; SingleEdgeElements:solver | SAVE :: NoTest, cumrt, cumct, Nodes | elmer-fe; edge/face basis (L15); writes Solver%Variable |
| `ElementUnitTest_h1/SingleElements.F90` (+2: ElementUnitTest_p1, SD_ElementUnitTest_p1) | SingleElements_init:solver _init; SingleElements:solver | SAVE :: NoTest, cumrt, cumct, Nodes | elmer-fe; p/edge/face Element in .sif; writes Solver%Variable |
| `ExtrudedMesh/Poisson.f90` (+1: ExtrudedMeshBaseline1) | PoissonSolver:solver | SAVE STIFF, LOAD, FORCE, AllocationsDone, Nodes; init ; work arrays | default_* solve path; elmer-fe; Extruded Mesh Levels |
| `FaceElement2D_BCs/FaceFEInterpolationTest.F90` | BestApproximationSolver:solver | SAVE STIFF, LOAD, FORCE, Acoef, AllocationsDone, Nodes, Indices; init ; work arrays | default_* solve path; elmer-fe; edge/face basis (L15); p/edge/face Element in .sif; writes Solver%Variable |
| `FaceElement3D_BCs/FaceFEInterpolationTest.F90` | BestApproximationSolver:solver | SAVE STIFF, LOAD, FORCE, Acoef, AllocationsDone, Nodes, Indices; init ; work arrays | default_* solve path; elmer-fe; edge/face basis (L15); p/edge/face Element in .sif; writes Solver%Variable |
| `HelmholtzEdge/Hrot.f90` | HrotSolver:solver | SAVE STIFF, LOAD, MASS, FORCE, Acoef, Bcoef, &, Nodes; init ; work arrays | default_* solve path; elmer-fe; boundary parent lookup; edge/face basis (L15); p/edge/face Element in .sif; direct CSR row surgery; current-element context; writes Solver%Variable |
| `HelmholtzFace/Helmholtz.f90` | HelmholtzHDIVSolver:solver | SAVE STIFF, LOAD, MASS, FORCE, Acoef, Bcoef, &, Nodes; init ; work arrays | default_* solve path; elmer-fe; boundary parent lookup; edge/face basis (L15); p/edge/face Element in .sif; direct CSR row surgery; current-element context; writes Solver%Variable |
| `L2ndDerivatives/L2ndDerivativesTest.f90` | TestSolver:solver | SAVE Nodes | elmer-fe; writes Solver%Variable |
| `MeshColour/MeshColour.F90` | MeshColour_init:solver _init; MeshColour:solver | work arrays | METIS C binding; writes Solver%Variable |
| `ModelPDEhandle/ModelPDEhandle.F90` | AdvDiffSolver:solver | SAVE Nodes | default_* solve path; elmer-fe; boundary parent lookup; writes Solver%Variable |
| `ModelPDEsplitted/ModelPDEsplitted.f90` (+1: SD_ModelPDEsplitted) | AdvDiffSolver_bulk:solver; AdvDiffSolver_boundary:solver | SAVE Nodes | default_* solve path; elmer-fe; boundary parent lookup; p/edge/face Element in .sif |
| `PoissonML/PoissonML.f90` | PoissonSolver:solver | init diag=0._dp, eCount=0; work arrays | elmer-fe; boundary parent lookup; p/edge/face Element in .sif; writes Solver%Variable |
| `PoissonPFEM/Poisson.f90` | PoissonSolver:solver | SAVE STIFF, LOAD, FORCE, AllocationsDone, Nodes; init ; work arrays | default_* solve path; elmer-fe; p/edge/face Element in .sif |
| `Q1Q0/Q1Q0.f90` | StokesSolver:solver | SAVE Nodes; init Norm=0, mu=1.0d0; work arrays | default_* solve path; elmer-fe; p/edge/face Element in .sif; writes Solver%Variable |
| `RaviartThomas2D/ComputeError.F90` | IntegrateError:solver | SAVE SolDOFs, Indices, AllocationsDone; init ; work arrays | elmer-fe; edge/face basis (L15); writes Solver%Variable |
| `RaviartThomas3D/ComputeError.F90` | IntegrateError:solver | SAVE SolDOFs, Indices, AllocationsDone; init ; work arrays | elmer-fe; edge/face basis (L15); writes Solver%Variable |
| `SD_Classical2DShell/BenchmarkShell.F90` (+4: SD_p-FEM_two_solvers, SD_p-FEM_with_varying_p, p-FEM_two_solvers, p-FEM_with_varying_p) | ShellSolver:solver | SAVE Nodes | default_* solve path; elmer-fe; p/edge/face Element in .sif; writes Solver%Variable |
| `SD_P2ndDerivatives/P2ndDerivativesTest.f90` | TestSolver:solver | SAVE Nodes; work arrays | elmer-fe; p/edge/face Element in .sif; writes Solver%Variable |
| `SD_PoissonPFEM/Poisson.f90` (+2: SD_SolverSpecificMeshes, SolverSpecificMeshes) | PoissonSolver:solver | SAVE STIFF, LOAD, FORCE, AllocationsDone, Nodes; init ; work arrays | default_* solve path; elmer-fe; p/edge/face Element in .sif |
| `SecondOrderEdgeElement2D/EdgeElementTest.F90` | EdgeElementSolver:solver | SAVE STIFF, LOAD, FORCE, Acoef, AllocationsDone, Nodes, Indices, PrismBasis, PrismCurlBasis, Br; init BricksVisited=.FALSE., PyramidsVisited=.FALSE.; work arrays | default_* solve path; elmer-fe; edge/face basis (L15); p/edge/face Element in .sif; writes Solver%Variable |
| `SecondOrderEdgeElement3D/EdgeElementTest.F90` | EdgeElementSolver:solver | SAVE STIFF, LOAD, FORCE, Acoef, AllocationsDone, Nodes, Indices, PrismBasis, PrismCurlBasis, Br; init BricksVisited=.FALSE., PyramidsVisited=.FALSE.; work arrays | default_* solve path; elmer-fe; edge/face basis (L15); p/edge/face Element in .sif; writes Solver%Variable |
| `SlaveSolverSlots/CountSolver.F90` | CountSolver:solver | init VisitedTimes=0 | writes Solver%Variable |
| `SlaveSolverSlotsScan/CountSolver2.F90` | CountSolver:solver | init VisitedTimes=0 | writes Solver%Variable |
| `StokesPFEM/Stokes.f90` | StokesSolver:solver | SAVE STIFF, MASS, LOAD, FORCE, rho, mu, Velocity, AllocationsDone, Nodes; init Norm=0, mu=1.0d0; work arrays | default_* solve path; elmer-fe; p/edge/face Element in .sif |
| `StokesProj/Velocity.f90` | VelocitySolver:solver | SAVE STIFF, LOAD, FORCE, rho, mu, Velocity, Pressure, AllocationsDone, Nodes; init Norm=0, mu=1.0d0; work arrays | default_* solve path; elmer-fe |
| `StokesProj/Pressure.f90` | PressureSolver:solver | SAVE STIFF, LOAD, FORCE, Velocity, AllocationsDone, Nodes; init Norm=0; work arrays | default_* solve path; elmer-fe; writes Solver%Variable |
| `bentonite/BentoniteSolver.F90` | BentoniteSolver:solver | SAVE STIFF, MASS, FORCE, ElementNodes, &; init ; work arrays | default_* solve path; elmer-fe; boundary parent lookup; writes Solver%Variable |
| `tresca/Tresca.f90` | TrescaSolver:solver | SAVE TrescaValues; work arrays | Model%Variables access |

### UMAT

| Path (under `fem/tests/`) | Hook kind and entries | State | Blockers |
|---|---|---|---|
| `UMAT_StVenant_2D/MyUMAT.F90` | stvenant_kirchhoff:UMAT | none | none (ModelView only) |
| `UMAT_StVenant_axials/LocalUMAT.F90` | UMAT_template:UMAT; linear_isotropic:UMAT; stvenant_kirchhoff:UMAT; stvenant_kirchhoff_axials:UMAT; hencky_stvenant_kirchhoff:UMAT | init PriLWork=102 | none (ModelView only) |
