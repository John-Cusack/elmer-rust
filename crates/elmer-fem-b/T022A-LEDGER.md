# T022-A lane ledger (deferred A–L solver-plugin ports)

Scope: the 44 `deferred-lane-budget` rows of T021A-LEDGER.md (basenames
A–L, each with its own reference `.so` in the serial build). Count note:
the T021-A report text says "41 deferred"; the ledger itself lists 44
`deferred-lane-budget` rows, and this lane dispositioned all 44.
Statuses: ported = Rust cdylib in `ports/<name>` with `nm -D` equality
plus ≥2 byte-identical `.sif` differentials, a repeated-invocation run,
and a branch-liveness control (see T022A-EVIDENCE.md). No refusals: every
deferred row admitted proof, including the two triage-first hazards
(CoordinateTransform module procedures, HelmholtzBEM matvec/precond),
both resolved as plain-ABI on inspection.

Reference build for all proofs: serial Elmer tree
`/home/john/firstmate/data/rs-elmer-t020/scratch/elmer-build` (campaign
configure minus MPI; both arms share it, so every differential is
self-consistent). Toolchain pins asserted by every proof script:
gfortran 13.3.0, rustc 1.100.0-nightly (feaadeeac 2026-09-19).

| unit | status | reason |
|---|---|---|
| fortran:fem/src/modules/Acoustics.F90 | ported | `ports/acoustics`: `acousticssolver_` + `_init_`; driver + LocalMatrix family in Rust, iterative/BEM kernels verbatim in wrappers. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/AdvectionDiffusion.F90 | ported | `ports/advectiondiffusion`: 2 entries; entry/init + loop nest in Rust, 2488-line body verbatim in wrappers. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/AdvectionReaction.F90 | ported | `ports/advectionreaction`: 3 entries incl `_post`; all Elmer traffic via wrappers. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/ArtificialCompressibility.F90 | ported | `ports/artificialcompressibility`: 3 entries; evident-intent divergence documented (pressure-values path unexercised). 2 diffs + reruns + liveness GREEN |
| fortran:fem/src/modules/BeamSolver3D.F90 | ported | `ports/beamsolver3d`: `timoshenkosolver_` + `_init0_`; pure driver port over 1-D beam mesh. 3 diffs + rerun + liveness GREEN |
| fortran:fem/src/modules/CahnHilliard.F90 | ported | `ports/cahnhilliard`: coupled 2-DOF nonlinear loop; FP/I-O in wrappers. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/CircuitsAndDynamics.F90 | ported | `ports/circuitsanddynamics`: 6 entries; harmonic + transient proofs; dropped-`CM =>` transcription bug caught by differential. 2 configs + reruns + liveness GREEN |
| fortran:fem/src/modules/CoilSolver.F90 | ported | `ports/coilsolver` (accessor `coil_acc.F90`): 4 entries; 4 diffs incl UDF entries per node + rerun + liveness GREEN; reference in-process re-entry aborts (double ALLOCATE), proven port-vs-port per bar |
| fortran:fem/src/modules/CompressibleNS.F90 | ported | `ports/compressiblens`: compressible flow coupling; dropped whole-array zeroing caught by differential. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/CoordinateTransform.F90 | ported | `ports/coordinatetransform` (+ `vecutils_rep.F90`): triage-first hazard cleared — 4 `__vecutils_MOD_*` procs are explicit-shape REAL(8), plain F77 ABI, transcribed verbatim. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/CraigBamptonSolver.F90 | ported | `ports/craigbamptonsolver`: eigensolve + reduction staging on cube mesh; reduction files compared. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/CyclicConvergence.F90 | ported | `ports/cyclicconvergence`: SAVE valTable state wrapper-side; cyclic-simulation staging. 2 diffs + rerun + liveness GREEN |
| fortran:fem/src/modules/DataToFieldSolver.F90 | ported | `ports/datatoffieldsolver`: mask + point-data paths; allocate-once repair documented. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/DFTSolver.F90 | ported | `ports/dftm` (dir `dftm`): 3 entries incl `poissonsolver_` (own-.so splice, no conflict). 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/DirectionSolver.F90 | ported | `ports/directionsolver`: 3 entries incl model-surgery `_Init0`; write-only/unobservable skips documented. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/DistanceSolve.F90 | ported | `ports/distancesolve`: 3 entries; PDE + dummy + real branches; `mpi_allreduce_` via libmpi_stubs as reference. 3 diffs + rerun + liveness GREEN |
| fortran:fem/src/modules/DistributeSource.F90 | ported | `ports/distributesource`: two-mesh interpolation staging; SAVE allocate-if-needed (author's intent). 2 diffs + rerun + liveness GREEN |
| fortran:fem/src/modules/EdgeFEInterpolationTest.F90 | ported | `ports/edgefeinterpolationtest`: library-vs-entry name (`BestApproximationSolver`), campaign form. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/ElectricForce.F90 | ported | `ports/electricforce`: campaign form `StatElecForce`; upstream dangling-pointer segfault avoided (defined paths only). 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/Electrokinetics.F90 | ported | `ports/electrokinetics`: 6 UDF entries driven as Real Procedures; ABI verified via `execrealfunction`. 3 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/EliminateDirichlet.F90 | ported | `ports/eliminatedirichlet`: linsolve-hook boundary (INTEGER FUNCTION via `DoLinSolveProcs` C trampoline, 8 plain pointers); pointer-deref + off-by-one caught by differential. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/EliminatePeriodic.F90 | ported | `ports/eliminateperiodic`: hook boundary; delegates to reference `EliminateDirichlet.so` (NEEDED+RUNPATH replicated); non-SAVE ColN kept as wrapper SAVE (author's intent). 3 diffs + reruns + liveness GREEN |
| fortran:fem/src/modules/Elmer2OpenFoamIO.F90 | ported | `ports/elmer2openfoamio`: campaign form `Elmer2OpenFoamWrite`; type-prefixed keywords. 3 diffs + rerun + liveness GREEN |
| fortran:fem/src/modules/EMPort.F90 | ported | `ports/emport`: 4 entries; eigen + potentials formulations. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/EMWaveSolver.F90 | ported | `ports/emwavesolver`: 5 entries; `CONTIG` empty-macro transcribed verbatim; defined-path carve-out documented. 3 diffs + rerun + liveness GREEN |
| fortran:fem/src/modules/FacetShellSolve.F90 | ported | `ports/facetshellsolve`: campaign form `ShellSolver`; 2 transcription slips caught by differential. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FDiffusion3D.F90 | ported | `ports/fdiffusion3d`: reference reads uninitialized Symb/dSymb on Cartesian (self-disagreement proven) — zeroed per library's Cartesian branch, defined Cylindrical path proven. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FDiffusion.F90 | ported | `ports/fdiffusion`: campaign form `FourierDiffusionSolver`; stale-acc_elem bug caught by differential. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FilmCompressibility.F90 | ported | `ports/filmcompressibility`: FullMode + Differential paths; author's stale-Found read replicated via wrapper cell. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FilmFlowSolver.F90 | ported | `ports/filmflowsolver`: 3 entries; darcy/Newton/GradP path pinned (unpinned reference blows to 1e17 deterministically). 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FilterTimeSeries.F90 | ported | `ports/filtertimeseries`: transcendental/order-sensitive math in Fortran; reference-UB paths documented, avoided. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FixTangentVelo.F90 | ported | `ports/fixtangentvelo`: real 3-D mesh staging; reference re-entry is UB (stale Mesh/Params; segfault at Level 10) — port implements author's intent, defined path + port-vs-port rerun per bar. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FluidicForce.F90 | ported | `ports/fluidicforce`: T021-A EnergyRelease shape; shear files compared. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FluxSolver.F90 | ported | `ports/fluxsolver`: scalar + tensor conductivity paths; BulkValues-reuse path. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FourierLoss.F90 | ported | `ports/fourierloss`: curl/inexact-Simpsons multi-Rhs path; end-cycle firing verified. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FreeSurfaceReduced.F90 | ported | `ports/freesurfacereduced`: axi + cartesian mapping paths; FlowSolve `(s)` value-mask. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/FreeSurfaceSolver.F90 | ported | `ports/freesurfacesolver`: 6 entries; USF `_`/`_3` nm-proven but not .sif-reachable; RotInit negative control. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/GmshOutputReader.F90 | ported | `ports/gmshoutputreader`: Python-generated .msh staging; GOTO-10 sweep in wrappers. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/HelmholtzBEM.F90 | ported | `ports/helmholtzbem`: triage-first hazard cleared — `matvec_`/`precond_` are assumed-size COMPLEX/INTEGER (plain ABI), HUTI via ipar/dpar; plugin-private `B __globmatc_MOD_matrix` re-declared + globalized. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/HelmholtzProjection.F90 | ported | `ports/helmholtzprojection`: 4 entries; single-face-singular mesh footgun (reference diverges) fixed with full 6-face Dirichlet. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/HelmholtzSolve.F90 | ported | `ports/helmholtsolve`: complex Galerkin assembly in wrappers; PlaneWave/flux/bubbles paths. 2 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/Isosurface.F90 | ported | `ports/isosurface`: 12 CONTAINS internals stay Fortran-side (by-value/assumed-shape dummies = refusal class avoided by keeping them un-crossed); mesh.geo/stl compared. 3 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/TurbulenceModels/KESolverLegacy.F90 | ported | `ports/kesolverlegacy`: implicit-interface sequence association replicated; 3-body-minimum .sif footgun. 3 configs + rerun + liveness GREEN |
| fortran:fem/src/modules/TurbulenceModels/KomegaLegacy.F90 | ported | `ports/komegalegacy`: legacy-bubble + transient two-solve paths; Symb reads verified Cartesian-clean. 3 configs + rerun + liveness GREEN |

Campaign backstop: fresh tree @ `9f6af2f85`, campaign configure, all 6
T021-A + all 44 T022-A ports spliced - `ctest -L quick` 479/482, failures
exactly K=3 upstream segfaults #80/#536/#537 (log: lane `proof/`).
