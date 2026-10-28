# T021-A lane ledger (solver-entry ports, report §8 lane)

Scope: `fortran:fem/src/modules/**` units whose basename starts with A–L,
graded at the plugin boundary (rustsmith ADR-057 entry shape,
`fortran_plugin_entry_verdict`) + the `nm` symbol surface. Lane T021-B
takes the rest of the alphabet; this ledger owns A–L only.
Statuses: graded = Rust cdylib in `ports/<name>` with byte-identical
differentials green (see T021A-EVIDENCE.md); refused-* = recorded reason,
never forced; deferred-lane-budget = admissible and portable in
principle, not reached in smallest-first order (this lane graded 6/61
admissible A–L units).

Entry-gate regrade: rustsmith `elmer-regrade` driver pointed at
`fortran_plugin_entry_verdict` (temporary bin, removed after the run):
184 `fem/src/modules` file units → 139 admit / 45 refuse; A–L slice →
61 admit / 21 refuse. The per-file `.so` mapping below is `nm -D` of the
serial reference build (not source scan): five subdirectory groups
bundle several files into one `.so`, which the per-unit port boundary
cannot split (see refused-shared-lib rows).

| unit | status | reason |
|---|---|---|
| fortran:fem/src/modules/AllocateSolver.F90 | graded | T021-A entry shape: F77-global `allocatesolver_` + `allocatesolver_init_`, derived dummies opaque, no CONTAINS internals. Port crate `ports/allocatesolver`: Rust cdylib exporting exactly the reference `.so` dynamic symbols (`nm -D` equality); `_init` list mutation + `Info` print through generated-style Fortran `BIND(C)` wrappers. Differentials: 3 distinct `.sif` configs (default / `always` + transient x2 / `never`) byte-identical + rerun + dispatch-count liveness. See T021A-EVIDENCE.md |
| fortran:fem/src/modules/DivProjection.F90 | graded | T021-A entry shape: F77-global `divprojection_` entry, empty body (stub). Port crate `ports/divprojection`: entry-only Rust cdylib, exact `nm -D` surface, no accessors (no Elmer calls in the original). Differentials: steady + transient-x2 byte-identical + rerun; `never`-toggle dispatch-line control. See T021A-EVIDENCE.md |
| fortran:fem/src/modules/ElementSizeSolver.F90 | graded | T021-A entry shape: F77-global `elementsizesolver_` entry + 1 internal. Port crate `ports/elementsizesolver`: Galerkin assembly transcribed into Fortran wrappers (same-compiler floats, `ES12.4` reports in Fortran); per-call allocate/deallocate mirrored (vestigial `AllocationsDone` documented as unreplicated). Differentials: defaults vs exponent+weight transient-x2 byte-identical + rerun; size-report liveness. See T021A-EVIDENCE.md |
| fortran:fem/src/modules/DivergenceSolver.F90 | graded | T021-A entry shape: F77-global `divergencesolver_` + `divergencesolver_init_` + 1 internal. Port crate `ports/divergencesolver`: full assembly transcribed into wrappers; guard checks map to identical `Warn`+return / `Fatal`; `Level=5` timing reports structurally suppressed via standard `Max Output Level = 4` (deterministic `Level=4` norm kept). Differentials: cartesian / relative+constbulk transient-x2 / axi-symmetric + SaveScalars byte-identical + rerun; norm liveness. See T021A-EVIDENCE.md |
| fortran:fem/src/modules/ElementStats.F90 | graded | T021-A entry shape: F77-global `elementstats_` entry + 3 internals. Port crate `ports/elementstats`: mesh-statistics sweeps transcribed into wrappers incl. cohort `GOTO 100` re-pass; parallel branch transcribed as written (dead in serial by the same guard). Differentials: defaults / histogram+skip+warn / rerun byte-identical; histogram liveness; `Fatal` negative control byte-identical in both arms. See T021A-EVIDENCE.md |
| fortran:fem/src/modules/EnergyRelease.F90 | graded | T021-A entry shape: F77-global `releaseratesolver_` entry + 1 internal. Port crate `ports/energyrelease`: fracture-integral assembly transcribed into wrappers; SAVE work arrays first-mesh sized and kept across calls. Differentials: linear vs quadratic MATC source fields (+transient x2) byte-identical + rerun; `2G` liveness. Staged via owning AllocateSolver solvers + MATC ICs. See T021A-EVIDENCE.md |
| fortran:fem/src/modules/ADIOS2OutputSolver.F90 | refused-no-reference | gate-refused (module derived-type dummies) AND no reference `.so` in the build (optional ADIOS2 dependency): no port boundary to prove against |
| fortran:fem/src/modules/AmsPrecSolver.F90 | refused-gate | `__amsvectorsolver_init_MOD_parallelprojtransposemultiply`: pointer dummy (libgfortran realloc contract) |
| fortran:fem/src/modules/ArteryOutlet.F90 | refused-gate | `__outletcompute_init_MOD_outletinit`: untyped function result (no type prefix) |
| fortran:fem/src/modules/DCRComplexSolve.F90 | refused-gate | `__dcrcomplexsolver_MOD_dcrcomplexsolver_boundary_residual`: non-REAL(8) assumed-shape dtype |
| fortran:fem/src/modules/ResultOutputSolve/DXOutputSolver.F90 | refused-shared-lib | ships inside bundled `ResultOutputSolve.so` (with `dxoutputsolver_`, `gidoutputsolver_`, `elmerpostoutputsolver_`, `resultoutputsolver_`, `vtkoutputsolver_`, `vtuoutputsolver_` + VTK module procs): no per-file `.so` boundary; a byte-identical port must replicate the whole library |
| fortran:fem/src/modules/ResultOutputSolve/ElmerPostOutputSolver.F90 | refused-shared-lib | same bundled `ResultOutputSolve.so` as DXOutputSolver (also gate-refused: optional dummy + module derived dummies) |
| fortran:fem/src/modules/ElasticSolve.F90 | refused-gate | `__elasticsolver_init0_MOD_elasticsolver_boundary_residual`: non-REAL(8) assumed-shape dtype |
| fortran:fem/src/modules/FindOptimum.F90 | refused-gate | `__findoptimum_MOD_name`, `__findoptimum_MOD_or`: non-intrinsic dummy |
| fortran:fem/src/modules/FlowSolve.F90 | refused-gate | `__flowsolver_init_MOD_flowsolver_boundary_residual`: non-REAL(8) assumed-shape dtype |
| fortran:fem/src/modules/HydrostaticNSVec.F90 | refused-gate | `__hydrostaticnsvecfront_MOD_delegateto`: character dummy (hidden length param) |
| fortran:fem/src/modules/HeatSolveVec.F90 | refused-gate | `__heatsolvevecfront_MOD_delegate`, `__heatsolvevecfront_MOD_delegateto`: character dummy |
| fortran:fem/src/modules/HydrostaticNS.F90 | refused-gate | `__hydrostaticnsutils_MOD_api`: undeclared function result (implicit typing) |
| fortran:fem/src/modules/HeatSolveLegacy.F90 | refused-gate | `__heatsolverlegacy_init_MOD_heatsolverlegacy_boundary_residual`: non-REAL(8) assumed-shape dtype |
| fortran:fem/src/modules/HeatSolve.F90 | refused-gate | `__heatsolvefront_MOD_delegatetoheatlegacy`: character dummy |
| fortran:fem/src/modules/IncompressibleNSVec.F90 | refused-gate | `__incompressiblensvecfront_MOD_delegateto`: character dummy |
| fortran:fem/src/modules/IncompressibleNS.F90 | refused-gate | `__incompressiblelocalforms_MOD_api`: undeclared function result (implicit typing) |
| fortran:fem/src/modules/MagnetoDynamics/CalcFields.F90 | refused-gate | `__magnetodynamicscalcfields_init0_MOD_as`, `__magnetodynamicscalcfields_init0_MOD_not`: non-intrinsic dummy (also inside bundled `MagnetoDynamics.so`) |
| fortran:fem/src/modules/LevelSet/LevelSetTimestep.F90 | refused-gate | `__levelsettimestep_MOD_levelsettimestep`: untyped function result; `levelsettimestep`: untyped function result (also inside bundled `LevelSet.so`) |
| fortran:fem/src/modules/TurbulenceModels/Komega.F90 | refused-gate | `__komegafront_MOD_delegatetokomegalegacy`: character dummy |
| fortran:fem/src/modules/TurbulenceModels/KESolver.F90 | refused-gate | `__kesolverfront_MOD_delegatetokesolverlegacy`: character dummy |
| fortran:fem/src/modules/BatterySolver/BatteryUtils.F90 | refused-gate | `__batterymodule_MOD_butlervolmer`, `__batterymodule_MOD_butlervolmerupdate`: derived-type dummy (also inside bundled `BatterySolver.so`) |
| fortran:fem/src/modules/ResultOutputSolve/GidOutputSolver.F90 | refused-shared-lib | admits under the entry gate but ships inside bundled `ResultOutputSolve.so` (see DXOutputSolver): cannot splice a per-file `.so` without dropping sibling entries |
| fortran:fem/src/modules/BatterySolver/BatteryPost.F90 | refused-shared-lib | admits under the entry gate but ships inside bundled `BatterySolver.so` (module-procedure surface `__batterymodule_MOD_*`): no per-file `.so` boundary |
| fortran:fem/src/modules/BatterySolver/ElectrolyteCons.F90 | refused-shared-lib | same bundled `BatterySolver.so` as BatteryPost |
| fortran:fem/src/modules/BatterySolver/ElectrolytePot.F90 | refused-shared-lib | same bundled `BatterySolver.so` as BatteryPost |
| fortran:fem/src/modules/ConsistentSplitting/CompressibilitySolver.F90 | refused-shared-lib | admits under the entry gate but ships inside bundled `ConsistentSplitting.so` (4 F77 entries across 4 files, 3 of which belong to lane T021-B): the whole-library port crosses the lane split |
| fortran:fem/src/modules/LevelSet/LevelSetIntegrate.F90 | refused-shared-lib | admits under the entry gate but ships inside bundled `LevelSet.so` (6 entries incl. gate-refused `levelsettimestep_`): the whole-library port is out of unit scope (port drafted during triage, removed unproven) |
| fortran:fem/src/modules/LevelSet/LevelSetSolver.F90 | refused-shared-lib | same bundled `LevelSet.so` as LevelSetIntegrate |
| fortran:fem/src/modules/LevelSet/LevelSetCurvature.F90 | refused-shared-lib | same bundled `LevelSet.so` as LevelSetIntegrate |
| fortran:fem/src/modules/LevelSet/LevelSetDistance.F90 | refused-shared-lib | same bundled `LevelSet.so` as LevelSetIntegrate |
| fortran:fem/src/modules/MagnetoDynamics/JfixPotentialSolver.F90 | refused-shared-lib | admits under the entry gate but ships inside bundled `MagnetoDynamics.so` (module-procedure surface `__magnetodynamicsutils_MOD_*` + gate-refused siblings): no per-file `.so` boundary |
| fortran:fem/src/modules/HarmonicSignals.F90 | refused-new-abi | admits under the entry gate (`__harmutils_MOD_sinsum/cossum` fit the substitution profile + F77 `source_`), but the `.so` ALSO exports the module procedures whose assumed-shape `REAL(8)` dummies arrive as gfortran descriptors: replicating the callee-side descriptor ABI from Rust is an unmeasured caller shape (T020 fixtures cover Rust-as-caller only) — refused, no new gate shapes |
| fortran:fem/src/modules/BeamSolver3D.F90 | deferred-lane-budget | admissible, own `.so` (`timoshenkosolver_` + `timoshenkosolver_init0_`); staging needs a 1-D beam mesh (`Element = p:1 b:1`) + SolidMechanicsUtils optional-arg calls; not reached in smallest-first order |
| fortran:fem/src/modules/CahnHilliard.F90 | deferred-lane-budget | admissible, own `.so`; coupled 2-DOF system with `!$omp` parallel assembly + nonlinear iteration loop; needs coupled concentration/potential staging; not reached |
| fortran:fem/src/modules/CyclicConvergence.F90 | deferred-lane-budget | admissible, own `.so` (2 entries); needs cyclic-simulation staging (`Periodic Cycle` variable, multi-cycle run) for the SAVE valTable state; not reached |
| fortran:fem/src/modules/CraigBamptonSolver.F90 | deferred-lane-budget | admissible, own `.so`; Craig-Bampton eigenmode reduction staging (eigensolves); not reached |
| fortran:fem/src/modules/CoordinateTransform.F90 | deferred-lane-budget | admissible, own `.so`; coordinate-transform + polar-decomposition staging with SAVE rotation state; not reached |
| fortran:fem/src/modules/CompressibleNS.F90 | deferred-lane-budget | admissible, own `.so`; compressible flow staging (density/pressure/temperature coupling); not reached |
| fortran:fem/src/modules/CoilSolver.F90 | deferred-lane-budget | admissible, own `.so` (1620 LOC); coil electromagnetics staging; not reached |
| fortran:fem/src/modules/CircuitsAndDynamics.F90 | deferred-lane-budget | admissible, own `.so` (1940 LOC); circuit-dynamics coupled staging; not reached |
| fortran:fem/src/modules/DirectionSolver.F90 | deferred-lane-budget | admissible, own `.so` (3 entries incl. model-surgery `_Init0` that appends a DG solver to `Model % Solvers`); coil-direction staging; not reached |
| fortran:fem/src/modules/DistributeSource.F90 | deferred-lane-budget | admissible, own `.so`; needs two-mesh nonconforming interpolation staging; not reached |
| fortran:fem/src/modules/DistanceSolve.F90 | deferred-lane-budget | admissible, own `.so`; distance-function staging; not reached |
| fortran:fem/src/modules/DataToFieldSolver.F90 | deferred-lane-budget | admissible, own `.so`; data-to-field staging; not reached |
| fortran:fem/src/modules/DFTSolver.F90 | deferred-lane-budget | admissible, own `.so`; discrete-Fourier staging; not reached |
| fortran:fem/src/modules/EdgeFEInterpolationTest.F90 | deferred-lane-budget | admissible, own `.so`; edge-element interpolation staging; not reached |
| fortran:fem/src/modules/EliminateDirichlet.F90 | deferred-lane-budget | admissible, own `.so`; Dirichlet-elimination staging; not reached |
| fortran:fem/src/modules/ElectricForce.F90 | deferred-lane-budget | admissible, own `.so`; electrostatic force staging; not reached |
| fortran:fem/src/modules/Electrokinetics.F90 | deferred-lane-budget | admissible, own `.so`; electrokinetics staging; not reached |
| fortran:fem/src/modules/Elmer2OpenFoamIO.F90 | deferred-lane-budget | admissible, own `.so`; OpenFOAM-mesh IO staging; not reached |
| fortran:fem/src/modules/EliminatePeriodic.F90 | deferred-lane-budget | admissible, own `.so` (boundary-leak pair with EliminateDirichlet per report §7); periodic-elimination staging; not reached |
| fortran:fem/src/modules/EMWaveSolver.F90 | deferred-lane-budget | admissible, own `.so`; electromagnetic-wave staging; not reached |
| fortran:fem/src/modules/EMPort.F90 | deferred-lane-budget | admissible, own `.so`; port-mode staging; not reached |
| fortran:fem/src/modules/FixTangentVelo.F90 | deferred-lane-budget | admissible, own `.so` (2 entries); staging needs a 3-D mesh + velocity field + `Tangent Surface` BCs (+ `ParallelReduction`, serially dead); not reached |
| fortran:fem/src/modules/FreeSurfaceReduced.F90 | deferred-lane-budget | admissible, own `.so`; reduced free-surface staging; not reached |
| fortran:fem/src/modules/FilterTimeSeries.F90 | deferred-lane-budget | admissible, own `.so`; time-series filter staging; not reached |
| fortran:fem/src/modules/FluidicForce.F90 | deferred-lane-budget | admissible, own `.so`; fluidic-force staging; not reached |
| fortran:fem/src/modules/FilmCompressibility.F90 | deferred-lane-budget | admissible, own `.so`; film-compressibility staging; not reached |
| fortran:fem/src/modules/FDiffusion.F90 | deferred-lane-budget | admissible, own `.so`; diffusion staging; not reached |
| fortran:fem/src/modules/FluxSolver.F90 | deferred-lane-budget | admissible, own `.so`; flux staging; not reached |
| fortran:fem/src/modules/FDiffusion3D.F90 | deferred-lane-budget | admissible, own `.so`; 3-D diffusion staging; not reached |
| fortran:fem/src/modules/FilmFlowSolver.F90 | deferred-lane-budget | admissible, own `.so`; film-flow staging; not reached |
| fortran:fem/src/modules/FourierLoss.F90 | deferred-lane-budget | admissible, own `.so` (boundary-leak pair with MagnetoDynamics per report §7); loss staging; not reached |
| fortran:fem/src/modules/FreeSurfaceSolver.F90 | deferred-lane-budget | admissible, own `.so`; free-surface staging; not reached |
| fortran:fem/src/modules/FacetShellSolve.F90 | deferred-lane-budget | admissible, own `.so`; shell-facet staging; not reached |
| fortran:fem/src/modules/GmshOutputReader.F90 | deferred-lane-budget | admissible, own `.so`; Gmsh-output staging; not reached |
| fortran:fem/src/modules/HelmholtzProjection.F90 | deferred-lane-budget | admissible, own `.so`; Helmholtz-projection staging; not reached |
| fortran:fem/src/modules/HelmholtzBEM.F90 | deferred-lane-budget | admissible, own `.so`; BEM staging; not reached |
| fortran:fem/src/modules/HelmholtzSolve.F90 | deferred-lane-budget | admissible, own `.so`; Helmholtz staging; not reached |
| fortran:fem/src/modules/Isosurface.F90 | deferred-lane-budget | admissible, own `.so`; isosurface staging; not reached |
| fortran:fem/src/modules/ArtificialCompressibility.F90 | deferred-lane-budget | admissible, own `.so`; artificial-compressibility staging; not reached |
| fortran:fem/src/modules/AdvectionReaction.F90 | deferred-lane-budget | admissible, own `.so`; advection-reaction staging; not reached |
| fortran:fem/src/modules/AdvectionDiffusion.F90 | deferred-lane-budget | admissible, own `.so`; advection-diffusion staging; not reached |
| fortran:fem/src/modules/Acoustics.F90 | deferred-lane-budget | admissible, own `.so` (4735 LOC); acoustics staging; not reached |
| fortran:fem/src/modules/TurbulenceModels/KomegaLegacy.F90 | deferred-lane-budget | admissible, own `.so` (`komegalegacy_` + `komegalegacy_init_`); k-omega turbulence staging (wall distance, coupled k/omega fields); not reached |
| fortran:fem/src/modules/TurbulenceModels/KESolverLegacy.F90 | deferred-lane-budget | admissible, own `.so` (`kesolverlegacy_` + `kesolverlegacy_init_`); k-epsilon turbulence staging; not reached |
