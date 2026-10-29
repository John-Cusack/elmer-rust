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
| `BulkAssembly`, `BoundaryAssembly` | `PoissonAssembly.f90` of `CoupledPoisson1`..`9`, `BlockPoisson1`, `BlockPoisson2` | element assembly kernels (`Bulk Assembly Procedure`) | `poisson_assembly::PoissonAssembly`; as `elmer-solver` `ElementAssembly` on `elmer-fe`: `hooks::{PoissonBulk, PoissonBoundary}` |
| `PoissonSolver` | `1dtests/Poisson.f90` (+ `1dtests_MeshLevels`), `ExtrudedMesh/Poisson.f90` (+ `ExtrudedMeshBaseline1`) | solver procedure (`Procedure "Poisson" "PoissonSolver"`) | `poisson_solver` |
| `TrescaSolver` | `tresca/Tresca.f90` | solver procedure; creates `Tresca` on `Stress`'s permutation (`SAVE TrescaValues`) | `tresca` |
| `SingleElements`, `SingleElements_init` | `ElementUnitTest_h1/SingleElements.F90` (also `ElementUnitTest_p1`, `SD_ElementUnitTest_p1`: the same source on `Element = "p:1"`) | solver procedure plus a `_init` section edit (`registry::init_sections`); `SAVE NoTest, NoMode`; p basis and rule from `elmer-pdofs` | `single_elements` |
| `AdvDiffSolver` | `ModelPDEhandle/ModelPDEhandle.F90` | solver procedure, keyword handles at integration points | `model_pde_handle` |
| `PressureSolver`, `VelocitySolver` | `StokesProj/Pressure.f90`, `StokesProj/Velocity.f90` | two solver procedures of a projection scheme; `SAVE Newton, Norm` | `stokes_proj` |
| `BentoniteSolver` | `bentonite/BentoniteSolver.F90` | 3-dof coupled solver with `CompEta` initialisation (`SAVE AllocationsDone`), axisymmetric metric, Robin boundaries | `bentonite` |
| `SingleEdgeElements`, `SingleEdgeElements_init` | `ElementUnitTest_edge*/SingleEdgeElements.F90` | solver procedure on `elmer-fe-hp`'s Whitney and Piola edge bases; the `_init` edit sets the variable, the edge element choice stays in the solver | `single_edge_elements` |
| `AlternatingSource` | `HeatControlExplicit/AlternatingSource.F90` | scalar function, `Variable "cpar, time"`, `SAVE sgn`, prints each flip | `udf::AlternatingSource` (not registered, see below) |

### Cases that match the pristine reference end to end

`tests/end_to_end_*.rs` run the cases' own `.sif` files on ElmerGrid meshes
through `elmer-run`'s `SifFrontend`, `elmer-mesh`, the mesh-load transforms
and `elmer-solver` with `elmer-run`'s solver registry plus these fixtures
(`tests/case/`, a copy of `elmer_run::run_case` that adds them). Each passes
its reference check and prints the `CompareToReferenceSolution` norm the
pristine ElmerSolver (26.2-devel, Rev 9f6af2f85, with the compiled `.so`
plugins) printed:

| Case | Fixtures | Pristine norm |
|---|---|---|
| `CoupledPoisson1` (`Assembly Solver`) | `PoissonBulk` | 7.78240309E-02 |
| `CoupledPoisson3`..`9` (`Coupled Solver`, constraints) | `PoissonBulk` | 6.96052509E-02, 5.10670175E+01, 5.85178077E+01, 7.56321899E-01, 7.84091700E-02, 5.25188761E-01, 2.40813064E+00 |
| `1dtests`, `1dtests_MeshLevels` (5 transient steps, `Mesh Levels = 8`) | `PoissonSolver`, `Source` | 1.40257363E+00, 1.41393619E+00 |
| `ExtrudedMesh` (`Extruded Mesh Layers = 3`) | `PoissonSolver` | 3.78170795E-02 |
| `tresca` (StressSolve, then the plugin) | `TrescaSolver` | 2.19639953E+01 |
| `ModelPDEhandle` | `AdvDiffSolver` | 1.79722562E+00 |
| `ElementUnitTest_h1` (12 `RUN` blocks) | `SingleElements` | all 12 reference values pass |
| `ElementUnitTest_edge`, `_edge_piola`, `_edge_quadratic` (12 `RUN` blocks each) | `SingleEdgeElements` | all 36 reference values pass |
| `StokesProj` (45 coupled iterations, as pristine) | `PressureSolver`, `VelocitySolver` | 3.79706752E-01 |
| `bentonite` (29 nonlinear iterations, as pristine) | `BentoniteSolver` | 1.95292467E+02 |

`StokesProj` names `Before Linsolve = "EliminateDirichlet"`, which nothing
implements and `elmer-solver` ignores. Eliminating Dirichlet rows does not
change the solution of the system, and the run matches the pristine
iteration count and final norm. The intermediate pressure-correction norms
differ from about the sixth digit (`0.59287902E-03` against `0.59287805E-03`
at iteration 45), within the iterative solver's tolerance.

`ElementUnitTest_h1` prints its per-element values. The six `SUM(MASS)`
values and four of the six `SUM(ABS(STIFF))` values are the pristine run's
digits exactly. The stiffness sums of 808 and 404 differ in the last place
(`5.3333333333333321` against `...330`, `...348` against `...339`). The mass
sums agree, so the difference is in `elmer-fe`'s `dBasisdx` rounding, not in
the kernel.

The edge cases print 36 values. All 24 Piola values, and 10 of the 12
Whitney values, are the pristine run's 17 digits (the values in
`elmer-fe-hp`'s `edge_unit_norms` test). The Whitney mass sums of 808 and
706 differ in the last place (`2.3333333333333330` against `...326`,
`0.99999999999999956` against `...967`). `SingleEdgeElements_init` writes an
`Element = "n:0 e:1 ..."` keyword that only the plugin's own `ElementInfo`
reads. The port keeps that choice in the solver instead, so `elmer-run`'s
element check and `elmer-pdofs` never see a non-nodal section for a
solver whose variable is a global.

The `tresca` test supplies one step the composition lacks: `StressSolver_Init0`
(StressSolve.F90:292) gives a section without `Variable` the default
`-dofs <dim> Displacement`. Without it, `elmer-physics-elastic` stops with
`'Variable DOFs' must equal the coordinate dimension`. Through `elmer-run`
today the case is therefore an error, not a mismatch.

Refused with a named reason:

- `ExtrudedMeshBaseline1`: its `PoissonSolver` runs on the baseline, boundary
  elements that a BC's `Body ID = 2` turns into a body. `elmer-solver`'s
  active elements are bulk elements only, so the solver refuses
  (`Boundary Condition 5 Body ID`) instead of solving an empty system.
  Solver 1 (HeatSolve, StructuredMeshMapper) matches the pristine norm.
- `HeatControlExplicit`: `AlternatingSource` is ported and tested against
  the pristine run's `PRINT *` lines, but it is **not registered**. Its
  argument `cpar` exists only through `Apply Explicit Control`, which
  `elmer-solver` does not implement, and `elmer-expr` evaluates a missing
  variable as zero where Elmer stops with `Can't find independent
  variable`. Registered, the case runs to norm 0 (a mismatch) instead of
  erroring. Register it once explicit control exists or `elmer-run`
  refuses `Apply Explicit Control`.
- `ElementUnitTest_p1`, `SD_ElementUnitTest_p1`: registered. The solver's
  variable is a global, so `elmer-solver` builds no dof layout for it, and
  `SingleElements` builds `elmer-pdofs`' layout for its own section in
  `init`. That layout has no tetrahedral p elements yet (the unit meshes
  have a 504), so both are refused by name before the run. On the brick of
  the unit mesh alone (`tests/data/ElementUnitTest_p1_brick`), the p:1 sums
  match the pristine values to 1e-15 (the stiffness sum differs in the last
  place, as for `_h1`). Wedges have no p basis in `elmer-pdofs` yet either. Its parentless 404 boundary face is "not an active p
  element" for `elmer-pdofs`, while the pristine run gives it a p basis (its
  stiffness sum differs from `_h1`'s).
- `BlockPoisson1`, `BlockPoisson2`: registered, but `Block Solver = True` is
  refused by `elmer-solver` (`coupled.rs`).

`registry::install(case, model, &UdfRegistry)` registers a case's scalar
fixtures under the names its `.sif` uses, as **new instances**, so two runs in
one process never share `SAVE`d state. `registry::lookup` builds a single
instance. `registry::register(case, &mut SolverRegistry)` adds the case's
solver procedures and element assembly routines; each routine the solver
builds from them starts with fresh `SAVE`d state. Keys are `(case, library,
entry)` because plugin names repeat across cases: `Source` alone is five
different files, `Poisson` `PoissonSolver` seven (see the inventory).

A composition running case `case` (the `fem/tests` directory name) does,
before `Simulation::setup`:

```rust
let mut registry = elmer_run::solver_registry();
elmer_fixtures::registry::register(case, &mut registry);
let mut sim = Simulation::new(model);
elmer_fixtures::registry::install_live(case, &mut sim, mesh); // Simulation::live_model + Backends::udfs
sim.setup(mesh, &registry, log, comm)?;
```

`elmer_solver::LiveModel` implements `ModelView`.

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

`tests/support/` holds the stand-ins the kernel unit tests use (ADR 003,
rule 6): a P1 triangle and two-node line basis with quadrature, and a mutable
model. The fixtures themselves integrate with `elmer-fe` and read
`elmer_solver::LiveModel`. `tests/case/` is the end-to-end composition;
`install_fe` keeps one basis per process, so there is one test binary per
coordinate dimension (`end_to_end_1d`, `_2d`, `_3d`). Case files and
ElmerGrid meshes are in `tests/data/` (`meshes/` holds the two shared by
several cases).

## What `elmer-run` needs to run these cases (not this crate's directory)

1. In `run_case`: `elmer_fixtures::registry::register(case, &mut registry)`
   after `solver_registry()`, and
   `elmer_fixtures::registry::install_live(case, &mut sim, mesh)` before
   `sim.setup`, with `case` the `fem/tests` directory name. `elmer-coverage`
   stages each test under its test name (`HeatControlExplicit_np3` for the
   `HeatControlExplicit` directory), so it should hand over the source
   directory's name rather than `case_dir.file_name()`.
2. `elmer_fixtures::registry::init_sections(case, &mut model)` before
   `Simulation::new`: the `_init` routines that only edit their section
   (`SingleElements_init` sets `Variable = -nooutput -global dummy_var`).
   `elmer-solver` has no pre-variable init hook, so the edit is applied to
   the model.
3. `preflight` refuses `Solver n: no Procedure` for sections with `Assembly
   Solver = True` or `Coupled Solver = True`, which `elmer-solver` runs itself
   (`coupled::solver_for_mode`). It should let those through (and `Block
   Solver`, which `elmer-solver` refuses by name). `tests/case/mod.rs` has the
   preflight the end-to-end tests use.
5. `StressSolver_Init0`'s default `Variable = -dofs <dim> Displacement`
   (for `tresca`; in `elmer-physics-elastic` or `elmer-run`).
4. `Apply Explicit Control` belongs in `unported.rs` (see `HeatControlExplicit`
   above). `elmer-expr`'s zero for a missing independent variable, where
   Elmer's `ListParseStrToValues` is `Fatal`, is the underlying gap.

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
(RaviartThomas3D). `TangentForce` needs
`Model % CurrentElement`, the element normal and the `Displacement Loads`
variable.

| Path (under `fem/tests/`) | Hook kind and entries | State | Blockers |
|---|---|---|---|
| `BDM2D/Source.F90` (+1: RaviartThomas2D) | SourceFun:scalar function | none | none (ModelView only) |
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
| `AdvReactDG_P/AdvReactDG_P.F90` | AdvectionReactionSolver:solver | SAVE MASS, STIFF, LOAD, FORCE, Velo, MeshVelo, Gamma, &; init ; work arrays | default_* solve path; elmer-fe; boundary parent lookup; edge/face basis (L15); p/edge/face Element in .sif; writes Solver%Variable |
| `BDM2D/ComputeError.F90` | IntegrateError:solver | SAVE SolDOFs, Indices, AllocationsDone; init ; work arrays | elmer-fe; edge/face basis (L15); writes Solver%Variable |
| `CurvedBndryPFEM/Poisson.f90` (+1: SD_CurvedBndryPFEM) | PoissonSolver:solver | SAVE STIFF, LOAD, FORCE, POT, AllocationsDone, Nodes; init ; work arrays | default_* solve path; elmer-fe; boundary parent lookup; p/edge/face Element in .sif; writes Solver%Variable |
| `CurvedBndryPFEM2/Poisson.f90` (+1: SD_CurvedBndryPFEM2) | PoissonSolver:solver | SAVE STIFF, LOAD, FORCE, POT, AllocationsDone, Nodes; init ; work arrays | default_* solve path; elmer-fe; p/edge/face Element in .sif |
| `FaceElement2D_BCs/FaceFEInterpolationTest.F90` | BestApproximationSolver:solver | SAVE STIFF, LOAD, FORCE, Acoef, AllocationsDone, Nodes, Indices; init ; work arrays | default_* solve path; elmer-fe; edge/face basis (L15); p/edge/face Element in .sif; writes Solver%Variable |
| `FaceElement3D_BCs/FaceFEInterpolationTest.F90` | BestApproximationSolver:solver | SAVE STIFF, LOAD, FORCE, Acoef, AllocationsDone, Nodes, Indices; init ; work arrays | default_* solve path; elmer-fe; edge/face basis (L15); p/edge/face Element in .sif; writes Solver%Variable |
| `HelmholtzEdge/Hrot.f90` | HrotSolver:solver | SAVE STIFF, LOAD, MASS, FORCE, Acoef, Bcoef, &, Nodes; init ; work arrays | default_* solve path; elmer-fe; boundary parent lookup; edge/face basis (L15); p/edge/face Element in .sif; direct CSR row surgery; current-element context; writes Solver%Variable |
| `HelmholtzFace/Helmholtz.f90` | HelmholtzHDIVSolver:solver | SAVE STIFF, LOAD, MASS, FORCE, Acoef, Bcoef, &, Nodes; init ; work arrays | default_* solve path; elmer-fe; boundary parent lookup; edge/face basis (L15); p/edge/face Element in .sif; direct CSR row surgery; current-element context; writes Solver%Variable |
| `L2ndDerivatives/L2ndDerivativesTest.f90` | TestSolver:solver | SAVE Nodes | elmer-fe; writes Solver%Variable |
| `MeshColour/MeshColour.F90` | MeshColour_init:solver _init; MeshColour:solver | work arrays | METIS C binding; writes Solver%Variable |
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

### UMAT

| Path (under `fem/tests/`) | Hook kind and entries | State | Blockers |
|---|---|---|---|
| `UMAT_StVenant_2D/MyUMAT.F90` | stvenant_kirchhoff:UMAT | none | none (ModelView only) |
| `UMAT_StVenant_axials/LocalUMAT.F90` | UMAT_template:UMAT; linear_isotropic:UMAT; stvenant_kirchhoff:UMAT; stvenant_kirchhoff_axials:UMAT; hencky_stvenant_kirchhoff:UMAT | init PriLWork=102 | none (ModelView only) |
