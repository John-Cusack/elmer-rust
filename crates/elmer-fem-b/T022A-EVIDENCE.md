# T022-A evidence (deferred A–L ports, lane rs-elmer-t022a)

Reference build for all proofs below: serial Elmer tree
`/home/john/firstmate/data/rs-elmer-t020/scratch/elmer-build`
(campaign configure minus MPI; both arms share it, so every differential
is self-consistent). Toolchain pins asserted by every proof script:
gfortran 13.3.0, rustc 1.100.0-nightly (feaadeeac 2026-09-19).
Elmer checkout (mesh inputs only): `/home/john/runs/elmer-fork`.

Norm class (lane-merged, wall-clock only): drop lines matching
`Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started
at`; mask timing *values* on per-iteration lines via `: +\(s\)`
(structure/iterations compared). Never widened beyond wall-clock.

Proof logs: lane `scratch/proof-<name>/` (per-config `run*/stdout.txt`,
`stderr.txt`, `scalars.dat` where applicable, `ref.syms`/`port.syms`).
Each proof rebuilds its cdylib against the reference `fmodules`
(`ELMER_MODDIR`), checks `nm -D` symbol equality, runs the shadow-dir
differentials (`ELMER_LIB` shadow with the port `.so` swapped in under
its own name, everything else symlinked), and halts on any divergence.
All runs under the lane 3-slot `measure-slot.sh` guard with isolated
`CARGO_TARGET_DIR`s and target-dir-aware `PORT_SO`.

## beamsolver3d (`ports/beamsolver3d`, proof-beamsolver3d.sh)

- Unit `fortran:fem/src/modules/BeamSolver3D.F90` (194 code lines).
- Symbol surface: `T timoshenkosolver_`, `T timoshenkosolver_init0_`. Diff empty.
- Config A (steady baseline), B (`Displace Mesh` + higher load), C (transient Bossak `MassAssembly`): stdout/stderr/scalars byte-identical modulo wall-clock lines.
- Rerun of A on port byte-identical. Liveness: tip values differ across A/B/C in both arms.
- `src/lib.rs` 139, accessor 177. Pure driver port (`BeamStiffnessMatrix` assembles internally); `DisplaceMesh` assumed-shape traffic stays in Fortran. PROOF GREEN.

## fixtangentvelo (`ports/fixtangentvelo`, proof-fixtangentvelo.sh)

- Unit `fortran:fem/src/modules/FixTangentVelo.F90` (215 lines).
- Symbol surface: `T fixtangentvelo_`, `T fixtangentvelo_init_`. Diff empty.
- Config A (defaults, 3-D cube mesh, 242 nodes fixed) + rerun, config B (explicit `Active Coordinate = 3`, `Fix Downstream Only`, `Fix Tangent Velocity`): byte-identical.
- Liveness: A scalars `0 0.2 0.2 0.2` vs B `0 0.2 ~0 0.2` in both arms.
- Staging must be a REAL 3-D mesh (2-D extrusion adds no top/bottom faces; vertical walls make the 3x3 system singular — upstream tags top/bottom only).
- Reference re-entry is UB: non-`SAVE` Mesh/Params set only on first visit → all-zero FixVals at Level 4, segfault at Level 10 (`FixTangentVelo.F90:125`, backtrace captured). Port implements author's intent (current Mesh/Params every call); proof covers the defined single-invocation path + port-vs-port rerun per bar item 3.
- `src/lib.rs` 117, accessor 254. PROOF GREEN.

## cyclicconvergence (`ports/cyclicconvergence`, proof-cyclicconvergence.sh)

- Unit `fortran:fem/src/modules/CyclicConvergence.F90` (265 lines).
- Symbol surface: `T cyclicconvergence_`, `T cyclicconvergence_init_`. Diff empty.
- Config A (constant `Convergence Value`, converges cycle 2 → finish cycle 3) + rerun; config B (`Convergence Value Name = time`, all tolerances, Min Iterations, custom file, 2 production cycles, converges 5 → 7): byte-identical.
- SAVE state (valTable/tsize/iCycle/prev*/ConvCycle/Converged/ConvVisited) wrapper-side, identical lifetimes; convergence-file writes in Fortran.
- Liveness: cycle counts differ across A/B in both arms. `src/lib.rs` 123, accessor 311. PROOF GREEN.

## cahnhilliard (`ports/cahnhilliard`, proof-cahnhilliard.sh)

- Unit `fortran:fem/src/modules/CahnHilliard.F90` (284 lines).
- Symbol surface: exactly `T cahnhilliardsolver_`. Diff empty.
- Config A (steady, 2 nonlinear iters, Mobility + Mixing Energy Density + Interface Thickness, UMFPack): byte-identical. Config B (transient 2x3 iters, Surface-Tension + h²-mobility + ElementDiameter + Convection-Velocity + Default1stOrderTime): byte-identical + rerun byte-identical.
- Divergence-style whole-element body in Fortran (automatic 2·nd arrays inner, bounds bind on entry); Rust drives nonlinear loop/assembly order (serial index order = serial reference build, no OpenMP symbols). Set-but-never-read `NewtonLinearization` replicated as wrapper SAVE; `BoundaryLocalMatrix` empty stub driver-calls-only. Recipe deviation: added `-fPIC` for the `Message` global reference.
- Liveness: `2G`-class fracture-equivalent reports differ across A/B in both arms. `src/lib.rs` 116, accessor 268. PROOF GREEN.

## craigbamptonsolver (`ports/craigbamptonsolver`, proof-craigbamptonsolver.sh)

- Unit `fortran:fem/src/modules/CraigBamptonSolver.F90` (317 lines).
- Symbol surface: `T craigbamptonsolver_`, `T craigbamptonsolver_init_`. Diff empty.
- Config A (stress eigensolve 12 modes + constraint modes + CB with `Calculate Matrix Norm`/`Save Reduction Basis`/`Save Reduction Matrix`) + rerun; config B (norm+basis branches off): byte-identical stdout/stderr/scalars.dat/`ReductionHeader/Indeces/Modes/Stiff/Mass.dat`.
- Liveness: B writes no basis/norms, A does, both arms. Proof wall ≈ 430 s (real 8³ cube mesh from `StressConstraintModes3`). `src/lib.rs` 175, accessor 335. PROOF GREEN.

## distributesource (`ports/distributesource`, proof-distributesource.sh)

- Unit `fortran:fem/src/modules/DistributeSource.F90` (389 lines).
- Symbol surface: exactly `T distributesource_`. Diff empty.
- Config A (two-mesh donut→sphere from `DistributeSourceMgdyn`, `Before all`, SaveScalars probe) + rerun; config B: byte-identical incl `scalars.dat`.
- Wrapper-side SAVE Nodes allocate-if-needed (author's first-mesh-sizing intent). The no-`g fix` Warn branch agrees on both arms but can't serve as exit-0 differential (downstream Solver 2 fails symmetrically).
- Liveness holds. `src/lib.rs` 444, accessor 441. PROOF GREEN wall 4 s.

## directionsolver (`ports/directionsolver`, proof-directionsolver.sh)

- Unit `fortran:fem/src/modules/DirectionSolver.F90` (410 lines).
- Symbol surface: `T directionsolver_`, `T directionsolver_dummy_` (empty), `T directionsolver_init0_` (DG guard + full model surgery: equation widening, DG dummy-solver append with exported `<Variable> Direction`). Diff empty.
- Config A (Source/plain-Dirichlet) + port-vs-port rerun; config B (`w`-gated BoundaryCondition, ±electrode loading, absent body force, namespaced `body 2: w` BC, 2 steady iterations in-process re-entry): byte-identical.
- SAVE (`STIFF/LOAD/FORCE`, `Nodes`, `First/varname/directionvar` incl. undefined-at-first-compare `varname`) replicated wrapper-side; write-only `visited` + dead `getSecondDerivatives` skipped as unobservable (documented).
- Liveness holds. `src/lib.rs` 427, accessor 397. PROOF GREEN.

## eliminatedirichlet (`ports/eliminatedirichlet`, proof-eliminatedirichlet.sh)

- Unit `fortran:fem/src/modules/EliminateDirichlet.F90` (452 lines).
- Symbol surface: exactly `T eliminatedirichlet_`. Diff empty.
- Boundary class: linsolve hook — `INTEGER FUNCTION(Model, Solver, A, b, x, n, DOFs, Norm)` reached through `Load.c`'s `DoLinSolveProcs` C trampoline (8 plain pointers, int result; backtrace-confirmed). Exercised via campaign `Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"` (cf. `PredictorCorrector` test).
- Config A (full-Dirichlet Poisson + hook) + rerun; config B (partial Dirichlet + Direct solver, 2 steady iters): byte-identical.
- Transcription footguns caught by the differential: gfortran passes scalar `TYPE(Matrix_t), POINTER` dummy as the address of the pointer variable (dedicated ABI probe) — entry dereferences `a` once; real-fill off-by-one (`Rows(i+1)-1` bound, caught via `ScaleLinearSystem` backtrace). Complex fill/eigen map-back/PEs>1 transcribed verbatim, unexercised (documented).
- Liveness: Dirichlet 1.0→2.0 differs, both arms agree. `src/lib.rs` 225, accessor 532. PROOF GREEN.

## filtertimeseries (`ports/filtertimeseries`, proof-filtertimeseries.sh)

- Unit `fortran:fem/src/modules/FilterTimeSeries.F90` (455 lines).
- Symbol surface: exactly `T filtertimeseries_`. Diff empty.
- Config A (pass-through mean, 2 steady iterations) + rerun; config B (`square` operator + `Time Filter 1 = Real 2.0` phase branch, 2 iters): byte-identical.
- Pure-`square` B is liveness-dead on Potential ∈ [0,1] (min/max endpoints preserved), so B also takes the phase-scaling branch. All transcendental/order-sensitive math (`EXP/SIN/COS/SQRT/**`, `2·π·freq`, `INT` interval checks) + whole-array updates in Fortran; plain `+-*/` weight arithmetic in Rust. Reference-UB paths (uninitialized `time`/`t0`, unset `MeanDofs`, unassociated `MeanField` on `Reset Interval`) documented, avoided.
- Liveness live. `src/lib.rs` 509, accessor 334. 3 consecutive GREEN runs. PROOF GREEN.

## edgefeinterpolationtest (`ports/edgefeinterpolationtest`, proof-edgefeinterpolationtest.sh)

- Unit `fortran:fem/src/modules/EdgeFEInterpolationTest.F90` (464 lines).
- Symbol surface: `T bestapproximationsolver_`, `T bestapproximationsolver_init0_`. Diff empty. Campaign form `Procedure = "EdgeFEInterpolationTest" "BestApproximationSolver"`.
- Config A (quadratic, Test Mode 5, L2 error) + rerun; config B (linear Piola, Test Mode 2, `Material Param = 1.0`, curl error): byte-identical.
- Dead `np > 0` branch transcribed verbatim in wrapper; unreferenced main-level `Nodes` SAVE documented as skipped.
- Liveness: A/B diverge, both arms agree. `src/lib.rs` 261, accessor 387. PROOF GREEN.

## freesurfacereduced (`ports/freesurfacereduced`, proof-freesurfacereduced.sh)

- Unit `fortran:fem/src/modules/FreeSurfaceReduced.F90` (467 lines).
- Symbol surface: exactly `T freesurfacereduced_`. Diff empty.
- Config A (Axi Symmetric, no mapping) + rerun; config B (Cartesian + mapping + softer relaxation): byte-identical; norm() masks FlowSolve `(s)` timing values only (merged class + `Assembly: \(s\)|Solve: \(s\)` value mask, structure compared).
- Liveness: flux 1.36e-06 vs 4.50e-03 in both arms. `src/lib.rs` 494, accessor 439. PROOF GREEN.

## electricforce (`ports/electricforce`, proof-electricforce.sh)

- Unit `fortran:fem/src/modules/ElectricForce.F90` (504 lines).
- Symbol surface: exactly `T statelecforce_`. Diff empty. Campaign form `Procedure = "ElectricForce" "StatElecForce"`.
- Config A (single force boundary) + rerun; config B (two force boundaries + nodal/field variables + Moment About, 2 steady iters): byte-identical.
- Found-and-avoided upstream defect: `Calculate Nodal Force` without `Calculate Force Density` segfaults the reference (`ElectricForce.F90:172`, dangling pointer); config B sets both. Port mirrors the code statement-for-statement (would fault identically); proofs cover defined paths only.
- Liveness: net-force x `-1.386576E-10` vs `-1.391979E-10`, both arms. `src/lib.rs` 120, accessor 466. PROOF GREEN.

## coordinatetransform (`ports/coordinatetransform` + `vecutils_rep.F90`, proof-coordinatetransform.sh)

- Unit `fortran:fem/src/modules/CoordinateTransform.F90` (505 lines).
- TRIAGE-FIRST hazard cleared: the four `__vecutils_MOD_*` dummies are all explicit-shape `REAL(8)` (`a(3)` for `norm2`/`normalized`, `(3,3)` pairs for `jac`, `(3,3)` trio for `transform2rank`) — plain F77 ABI, zero descriptor-carrying dummies. Full port incl. the 4 module procs transcribed statement-for-statement into `accessors/vecutils_rep.F90` (same compiler/flags → exact ABI); accessor hidden-visibility, vecutils object default-visibility, version script globalizes exactly the 6 symbols.
- Symbol surface: `T rotmsolver_`, `T rotmsolver_init_`, `T __vecutils_MOD_jac`, `T __vecutils_MOD_norm2`, `T __vecutils_MOD_normalized`, `T __vecutils_MOD_transform2rank`. Diff empty.
- Config A (alpha/beta-gradient + polar decomposition, 1 iter) + rerun; config B (gamma path via body-section `Local Coordinate System Beta Reference and Gamma` + `coilcurrent e`, rotated refs, `Use Polar Decomposition = False` with both Warn defaults, body-4 `Beta Reference` CYCLE, 2 iters SAVE re-entry): byte-identical.
- Real bugs fixed en route: refs/flag must live in Body sections (`GetBodyParams` reads `Bodies(i)%Values`, not the force list); Rust `cta_node(j)` vs Fortran `cta_node(j, nn)` (garbage `nn` → 17 GB RSS kill).
- Liveness: A sums `2182.38…/1260…` vs B `1600.41…/924…` + B Warn line, both arms. `src/lib.rs` 142, accessors 330+73. PROOF GREEN.

## fluidicforce (`ports/fluidicforce`, proof-fluidicforce.sh)

- Unit `fortran:fem/src/modules/FluidicForce.F90` (517 lines).
- Symbol surface: exactly `T forcecompute_`. Diff empty. Campaign form `Procedure = "FluidicForce" "ForceCompute"`.
- Config A (lid-driven Stokes cavity, viscous) + rerun; config B (inviscid + Moment About 1/2/3 + Shear Stress Output + Velocity Field Name + Incompressible, 2 steady iters): byte-identical stdout/stderr/scalars.dat/scalars.dat.names/shearstress.dat(+names).
- Thin Rust entry (`ffc_begin` + `ffc_run`); whole `ForceCompute` body + `ForceIntegrate` internal statement-for-statement in the wrapper. Norm class extends within-kind: `File started at` (SaveScalars `.names` header) + `FINISHED AT` (actual line is `ELMER SOLVER FINISHED AT`).
- Liveness holds. `src/lib.rs` 66, accessor 471. PROOF GREEN.

## gmshoutputreader (`ports/gmshoutputreader`, proof-gmshoutputreader.sh)

- Unit `fortran:fem/src/modules/GmshOutputReader.F90` (543 lines).
- Symbol surface: `T gmshoutputreader_`, `T gmshoutputreader_init_`. Diff empty.
- Config A (clean-room square.grd + Python-generated `inputA.msh`: 240 nodes/420 303-triangles `$Nodes/$Elements/$NodeData`, `Temperature` onto solver-declared variable via SaveScalars — `InterpolateMeshToMeshQ` only fills pre-existing targets, same as upstream `Outlet2InletBC` test) + rerun; config B (`Use Quadrant Tree = Logical True`, different field data, 2 steady iters): byte-identical.
- Rust holds no state; label-10 `GOTO` sweep via `gmr_pass`/`gmr_need_alloc`/`gmr_allocate`; 2 CONTAINS internals as accessor module procs over SAVE state; all file I/O/reductions/formatted output in Fortran.
- Liveness holds. `src/lib.rs` 118, accessor 556. PROOF GREEN.

## elmer2openfoamio (`ports/elmer2openfoamio`, proof-elmer2openfoamio.sh)

- Unit `fortran:fem/src/modules/Elmer2OpenFoamIO.F90` (546 lines).
- Symbol surface: exactly `T elmer2openfoamwrite_`. Diff empty.
- Config A + rerun; config B (early-return branch); config C (field differs): all byte-identical incl `scalars.dat.names` (normed) + field files; liveness (A reaches All done, B returns early, C field differs) holds.
- Staging footguns: unlisted solver keywords (`Elmer Center Filename`, `Passive OpenFOAM Coordinate`, `Number of cells`) dropped with desync unless type-prefixed (`= File …`, `= Integer …`); parser lowercases unquoted values (field stages as `fieldc.dat`).
- `src/lib.rs` 201, accessor 621. PROOF GREEN.

## distancesolve (`ports/distancesolve`, proof-distancesolve.sh)

- Unit `fortran:fem/src/modules/DistanceSolve.F90` (547 lines).
- Symbol surface: `T distancesolver_`, `T distancesolver1_`, `T distancesolver_init_`. Diff empty. (`U mpi_allreduce_` resolves via `libmpi_stubs.so` exactly as the reference does.)
- Config A (PDE `DistanceSolver` with RelaxDT terms, `H Scale`, ~50-iter nonlinear loop × 2 steady iters through SAVE) + rerun; config B (geometric `DistanceSolver1` dummy/`distcomp0` branch with BodyForce-condition + BC + Noslip fixed-node loops); config C (real `distcomp` branch with `RANDOM_NUMBER` + `SortR`): all byte-identical.
- Liveness: B-vs-C differs, both arms agree. `src/lib.rs` 178, accessor 521. PROOF GREEN.

## electrokinetics (`ports/electrokinetics`, proof-electrokinetics.sh)

- Unit `fortran:fem/src/modules/Electrokinetics.F90` (561 lines).
- Symbol surface: all 6 UDF entries `T getjouleheat_`, `T helmholtz_smoluchowski_`, `T _1`, `T _2`, `T _3`, `T _comp_`. Diff empty.
- ABI verified as `(Model, Node, dummy-real-slot) -> real` via `execrealfunction` in `LoadMod.F90`; `helmholtz_smoluchowski_comp_` overlays its INTEGER direction on the caller slot like the Fortran callee.
- Config A (computed E-field + EO Mobility + HS1/HS2 + JouleHeat) + rerun; config B (constant E-field + Zeta Potential + tangential HS + HS3 + JouleHeat, `Steady State Min Iterations = 2`, 2 in-process iters); config C (B with doubled Zeta): all byte-identical. All FP (`NormalVector`, `SUM`, mobility quotient, tangential projection) + formatted I/O in accessor; Rust holds control flow only.
- Liveness: B-vs-C differs, both arms agree. `src/lib.rs` 329, accessor 460. PROOF GREEN.

## compressiblens (`ports/compressiblens`, proof-compressiblens.sh)

- Unit `fortran:fem/src/modules/CompressibleNS.F90` (564 lines).
- Symbol surface: exactly `T compressiblens_`. Diff empty. Campaign form `Procedure = "CompressibleNS" "CompressibleNS"`.
- Config A (body-force + bulk-viscosity, transient 2×BDF-1) + rerun; config B (no-force + default bulk viscosity + `p:1 b:1` bubbles + 2 nonlinear iters): byte-identical.
- Real bug caught by the differential: dropped per-element whole-array `LOAD = 0.0d0` zeroing showed as 0.3% `ComputeChange` NRM drift on first solve — fixed by zeroing in Rust orchestration before the associated-body-force fill.
- Liveness holds. `src/lib.rs` 368, accessor 419. PROOF GREEN.

## helmholtzprojection (`ports/helmholtzprojection`, proof-helmholtzprojection.sh)

- Unit `fortran:fem/src/modules/HelmholtzProjection.F90` (584 lines).
- Symbol surface: `T helmholtzprojector_`, `T helmholtzprojector_init_`, `T removekernelcomponent_`, `T removekernelcomponent_init0_`. Diff empty.
- Config A (Constant Bulk Matrix + 2 steady iters, incl ReadySystemMatrix force-only re-entry) + rerun; config B (full reassembly, stronger source): byte-identical.
- Staging footgun: proof mesh needs a real 3-D tet mesh (`cube3`) with phi Dirichlet on all 6 boundaries — single-face leaves the hp Laplacian singular and the reference itself diverges.
- Transposed `PotSol(1,n)` wrapper buffer corrupted the heap — caught, fixed. Liveness live. `src/lib.rs` 554, accessor 573. PROOF GREEN.

## filmcompressibility (`ports/filmcompressibility`, proof-filmcompressibility.sh)

- Unit `fortran:fem/src/modules/FilmCompressibility.F90` (606 lines).
- Symbol surface: `T filmcompressibility_`, `T filmcompressibility_init_`. Diff empty.
- Config A (FullMode + `AC Diffusion`, PDE assembly) + rerun; config B (Differential Mode + `AC Relaxation Factor 0.5`, diagonal assembly + WHERE-relaxation, 2 BDF steps): byte-identical.
- `.sif` needs Iterative linear-system keywords (`_init` sets no defaults, unlike DivergenceSolver); `Params => GetSolverParams()` must not be dropped (null-list segfault in first `ListGetCReal`). Author's stale-`Found` read after Hydrostatic-Pressure lookup replicated via single wrapper-side `Found` cell in exact call order. `Levelset Initialize` path transcribed, unexercised by design.
- Liveness: `AC Average` differs, both arms agree. `src/lib.rs` 186, accessor 589. PROOF GREEN.

## fluxsolver (`ports/fluxsolver`, proof-fluxsolver.sh)

- Unit `fortran:fem/src/modules/FluxSolver.F90` (704 lines).
- Symbol surface: `T fluxsolver_`, `T fluxsolver_init_`. Diff empty.
- Config A (scalar conductivity, flux-only, Rank-1) + rerun; config B (2×2-tensor conductivity Rank-3, flux+grad+magnitudes, Enforce Positive Magnitude, Constant Bulk Matrix with 2 steady iters taking the BulkValues-reuse path, Flux/Grad Abs, DG reads): byte-identical stdout/stderr/scalars.dat.
- Thin Rust orchestration (guards, sweep, element loop, finish, solve); whole BulkAssembly element body, face terms (`AddLocalFaceTerms`/`LocalJumps` incl SAVE nodes), per-component `DefaultSolve` loop, `FluxSolver_Init` fully in wrappers.
- Liveness: `Result Norm` A `3.4144869924308132` vs B `2.5683642073499798`, identical in both arms. `src/lib.rs` 169, accessor 590. PROOF GREEN.

## artificialcompressibility (`ports/artificialcompressibility`, proof-artificialcompressibility.sh)

- Unit `fortran:fem/src/modules/ArtificialCompressibility.F90` (748 lines).
- Symbol surface: `T compressibilityscale_`, `T compressibilityscale_init_`, `T compressibilitysolver_`. Diff empty.
- Config A (CompressibilitySolver, Reference Pressure branch, 2 steady iters) + rerun; config B (CompressibilityScale + `_init`); liveness C (Displaced Shape=True, arms agree, differs from A); liveness D (Weight By Displacement=True, arms agree, differs from B): all byte-identical modulo merged norm class.
- Rust owns control flow/branch defaults only; Gauss integration, SUM reductions, WRITE/PRINT, SAVE state in accessor. Documented evident-intent divergence: original assigns `PressureExists = .FALSE.` immediately after associating pressure pointers, so its element loop could only read uninitialized `ReferencePressure` when a pressure variable exists; port uses pressure values when present. No proof config creates a pressure variable — the defined Reference Pressure path is proven byte-identical.
- `src/lib.rs` 291, accessor 777. PROOF GREEN.

## eliminateperiodic (`ports/eliminateperiodic`, proof-eliminateperiodic.sh)

- Unit `fortran:fem/src/modules/EliminatePeriodic.F90` (748 lines).
- Symbol surface: exactly `T eliminateperiodic_`. Diff empty. Port carries reference-matching `NEEDED EliminateDirichlet.so` + RUNPATH for its two delegation calls (resolved to the reference `.so` exactly like the original).
- Config A (reduction) + rerun; config B (delegation/Direct/2-iter); config A2 (2-iter projection-reuse) + rerun: all byte-identical.
- Rust transcribes control flow + Perm-remap scalar loop (`(ColN-1)/DOFs+1`); matrix surgery/CRS helpers/Elmer calls verbatim in wrappers. Original's non-SAVE `ColN` kept as wrapper SAVE with identical contents (author's evident intent) — A2 orig-vs-port 2-iteration differential proves re-entry byte-identical.
- Liveness variant differs, both arms agree. `src/lib.rs` 165, accessor 681. PROOF GREEN.

## fdiffusion (`ports/fdiffusion`, proof-fdiffusion.sh)

- Unit `fortran:fem/src/modules/FDiffusion.F90` (758 lines).
- Symbol surface: exactly `T fourierdiffusionsolver_`. Diff empty. Campaign form `Procedure = File "FDiffusion" "FourierDiffusionSolver"`.
- Config A (Angular Frequency branch) + rerun; config B (Frequency + Warn + 2-iteration branches): byte-identical.
- Rust drives entry/iteration/element loops; all FP, `SUM`, `WRITE(Message)` timing, SAVE state in wrappers. Stale `acc_elem` in the boundary path (degenerate-element errors) caught by the differential, fixed.
- Liveness holds. `src/lib.rs` 101, accessor 486. PROOF GREEN.

## advectionreaction (`ports/advectionreaction`, proof-advectionreaction.sh)

- Unit `fortran:fem/src/modules/AdvectionReaction.F90` (808 lines).
- Symbol surface: `T advectionreactionsolver_`, `T _init_`, `T _post_`. Diff empty.
- Config A (all 3 entries incl nodal-average `_post`) + rerun; config B (Limit Solution + 2 nonlinear + 2 steady iters): byte-identical.
- Fidelity points: `CPUTIme()` == `CPUTime()` (case-insensitive same external); `FlowDOFs` zero-initialised (undefined only on unobservable mesh-velocity path); `Indexes(128)`/`1:N` broadcast bounds/pre-loop counter increments replicated wrapper-side; unlisted `.sif` keywords (`Limit Solution`, `Tracer Upper/Lower Limit`) type-prefixed.
- Liveness: A-vs-B branch delta confirmed. `src/lib.rs` 184, accessor 706. PROOF GREEN.

## datatoffieldsolver (`ports/datatoffieldsolver`, proof-datatoffieldsolver.sh)

- Unit `fortran:fem/src/modules/DataToFieldSolver.F90` (809 lines).
- Symbol surface: `T datatofieldsolver_`, `T datatofieldsolver_init_`. Diff empty.
- Config A (mask + `Mask Diffusion` + `Continue BC With Boundary Assembly` + Dirichlet, global diffusion, 2 steady iters = SAVE re-entry both arms) + rerun; config B (point-data file + `Normalize Data by Weight` alias + `Set Constant Weight Sum` + material `Data Diffusivity` local path): byte-identical stdout/stderr/scalars.dat.
- Documented repair: `AsciiPointsToMesh`'s never-guarded SAVE allocation made allocate-once; defined first-call path statement-identical; rerun is a fresh process in both arms.
- Liveness: `Mask is active` only in A, scalars differ. `src/lib.rs` 136, accessor 797. PROOF GREEN.

## helmholtzbem (`ports/helmholtzbem`, proof-helmholtzbem.sh)

- Unit `fortran:fem/src/modules/HelmholtzBEM.F90` (824 lines).
- TRIAGE-FIRST hazard cleared: `matvec`/`precond` carry only plain-ABI assumed-size dummies (`COMPLEX(KIND=dp) :: u(*),v(*)`, `INTEGER :: ipar(*)`) — no CHARACTER, no descriptors. `HUTI_*` names are CPP macros over `ipar`/`dpar` slots (verified in `huti_fdefs.h`); whole HUTI sequence transcribes with plain indexing; `IterCallFTNC` (descriptor dummies) reached via a Fortran wrapper.
- Symbol surface: `B __globmatc_MOD_matrix` (plugin-private `GlobMatC` matrix, re-declared same-statement in accessor, globalized in version script) + `T helmholtzbemsolver_`, `T matvec_`, `T precond_`. Diff empty.
- Config A (direct + mixed Dirichlet/flux) + rerun; config B (iterative HUTI + swapped BC mix, 2 steady iters): byte-identical. Reference iterative path verified self-deterministic (fixed HUTI seed) before porting.
- Staging: boundary elements need `Body Id = 1` to join the active set (else 0 boundary nodes); `P1`/`P2` sizing follows author's evident intent (documented).
- Liveness control diverges on both arms. `src/lib.rs` 627, accessor 687. PROOF GREEN.

## helmholtzsolve (`ports/helmholtsolve`, proof-helmholtsolve.sh)

- Unit `fortran:fem/src/modules/HelmholtzSolve.F90` (845 lines).
- Symbol surface: `T helmholtzsolver_`, `T helmholtzsolver_init_`. Diff empty.
- Config A + rerun; config B (PlaneWave-BC, Wave-Flux, Bubbles/condensation, 2 nonlinear iters in-process SAVE re-entry): byte-identical. Flow/structure-interface branches transcribed but correctly dormant (`need_iface=0`, same as reference).
- Complex Galerkin assembly, `MATMUL`/`SUM`, `CPUTime` reports in wrappers; Rust holds control flow. Norm: merged class + `iter: N Assembly/Solve (s)` value mask.
- Liveness: B differs, both arms agree. Proof run twice GREEN. `src/lib.rs` 137, accessor 683. PROOF GREEN.

## fdiffusion3d (`ports/fdiffusion3d`, proof-fdiffusion3d.sh)

- Unit `fortran:fem/src/modules/FDiffusion3D.F90` (937 lines).
- Symbol surface: exactly `T fourierdiffusionsolver_`. Diff empty.
- Config A (Angular Frequency + Electric Conductivity, Cylindrical on x-shifted unit cube) + rerun; config B (Frequency + Electrical→Warn + Electric + 2 nonlinear iters): byte-identical.
- Reference UB found and fixed per defined-path rule: `LocalMatrix` inits Metric but skips `CoordinateSystemInfo` on the Cartesian path, reading uninitialized `Symb`/`dSymb` — provably nondeterministic (same binary disagreed with itself run-to-run, ~1% norm swings). Port zeroes both (matches `CoordinateSystemInfo`'s Cartesian branch: `Symbols=0`); differentials run the defined `Cylindrical` path (documented in proof header + lib.rs). The 3-D `Solve:    (s)` tag (4 spaces) motivated the generalized `: +\(s\)` norm mask.
- Liveness live. `src/lib.rs` 111, accessor 645. PROOF GREEN.

## filmflowsolver (`ports/filmflowsolver`, proof-filmflowsolver.sh)

- Unit `fortran:fem/src/modules/FilmFlowSolver.F90` (1076 lines).
- Symbol surface: `T filmflowsolver_`, `T filmflowsolver_init_`, `T filmflowsolver_init0_`. Diff empty.
- Config A (laminar 2-D channel, 1 nonlinear iter) + rerun; config B (Convect + Gradp Discretization + darcy friction + Newton-after-1 + friction heating with heat source, 2 nonlinear iters, pressure pinned — unpinned reference blows to uniform 1e17 deterministically in both arms, pinned is the meaningful differential): byte-identical.
- All entry SAVE + CONTAINS-internal SAVE (Nodes/Visited/MinPres/MinSpeed) wrapper-side, identical lifetimes; `FrictionLawPraks`'s set-but-never-read `Visited` the one documented skip.
- Liveness: B differs, both arms agree. `src/lib.rs` 206, accessor 974. PROOF GREEN.

## fourierloss (`ports/fourierloss`, proof-fourierloss.sh)

- Unit `fortran:fem/src/modules/FourierLoss.F90` (1133 lines).
- Symbol surface: `T fourierlosssolver_`, `T fourierlosssolver_init_`. Diff empty.
- Config A (Direct/exact/summed, bodyloss.dat/seriesloss.dat compared) + fresh-process rerun; config B (curl/inexact-Simpsons/separate-components multi-Rhs path): byte-identical. End-cycles verified firing at steps 20/30 (loss trajectory 0 → 11.9 in scalars.dat).
- Stateless Rust orchestrator (all SAVE state — entry, DFT, BulkAssembly Nodes — wrapper-side, identical lifetimes); all float work, formatted/file I/O, `Rhs` pointer swaps verbatim in Fortran. Recipe delta: `-fPIC -fallow-argument-mismatch` (accessor references Elmer COMMON `CurrentModel`/`Message`). No FourierLoss-specific UB (no Symb/Metric traffic; rerun identical).
- Liveness: A-vs-B diverges, both arms agree. `src/lib.rs` 147, accessor 1099. PROOF GREEN.

## komegalegacy (`ports/komegalegacy`, proof-komegalegacy.sh)

- Unit `fortran:fem/src/modules/TurbulenceModels/KomegaLegacy.F90` (737 lines).
- Symbol surface: `T komegalegacy_`, `T komegalegacy_init_`. Diff empty.
- Config A (default legacy bubbles) + rerun; config B (p-bubbles + wall BCs, 2 iters); config C (transient, Init-injected 2 solves) + rerun (bx/BubbleTimestep re-entry): all byte-identical.
- Rust pure orchestration; all SAVE state (first-mesh arrays, bx/bxprev, ElementNodes) wrapper-side. `_Init` p-bubble + transient two-solve branches covered. Metric/Symb/dSymb reads verified strictly inside non-Cartesian ELSE branches (no FDiffusion3D-class UB on this path).
- Liveness: Result Norm differs across A/B/C, both arms. `src/lib.rs` 136, accessor 514, proof 315 lines. PROOF GREEN.

## kesolverlegacy (`ports/kesolverlegacy`, proof-kesolverlegacy.sh)

- Unit `fortran:fem/src/modules/TurbulenceModels/KESolverLegacy.F90` (1170 lines).
- Symbol surface: `T kesolverlegacy_`, `T kesolverlegacy_init_`. Diff empty. (`kewall_` stays implicit-external `U` in both arms, resolved from libelmersolver as in the original.)
- Config A (steady/standard-model/legacy-bubbles + Epsilon Wall BC, direct, 2 outer iters) + rerun; config B (transient BDF(1)×2/RNG/`Bubbles=False`, Default1stOrderTime, `_Init` legacy-bubble minimums); config C (transient/standard/default-bubbles, `CondensatePTransient` + bx history): all byte-identical.
- Transcription replicates even the original's implicit-interface array→scalar sequence association (`GetParentUVW`/`NormalVector`/`ElementInfo` scalar dummies called with `(1)` elements — bit-identical); dead stores (`at`/`at0`, `NewtonLinearization`, `ForceVector`) documented skips.
- Staging: all four Body sections declared (mesh carries bulk bodies 1–2 + boundary bodies 3–4; 3-body `.sif` segfaults in core `initcond` in *both* arms). Symb advisory checked: `Symb`/`Metric` read only inside non-Cartesian branch after `CoordinateSystemInfo` fills them — clean.
- Liveness: Result Norms 3.4e11 / 2.1e-3 / C-distinct, all non-vacuous. `src/lib.rs` 134, accessor 850. PROOF GREEN.

## freesurfacesolver (`ports/freesurfacesolver`, proof-freesurfacesolver.sh)

- Unit `fortran:fem/src/modules/FreeSurfaceSolver.F90` (1223 lines).
- Symbol surface: `T freesurfacesolver_`, `T freesurfacesolver_rotinit_`, `T freesurfacetomeshupdate_`, `T _1`, `T _2`, `T _3`. Diff empty.
- Config A (stabilized + global Maximum Displacement) + rerun; config B (bubbles + Apply Dirichlet limiters + local max displacement + Relaxation 0.7 + RotInit entry + `FreeSurfaceToMeshUpdate1/2` user functions, 2 nonlinear iters): byte-identical.
- Exercise evidence: `Result Norm` lines present in both configs; B shows 45 upper-constrained limiter points (tuned `Max FreeSurface = 0.002` so `ZeroRow/SetMatrixElement` fires); loader lines confirm both USFs loaded from the port `.so` in both arms; negative control (RotInit keyword removed) yields `ERROR:: FreeSurfaceSolver_RotInit: Can't find 'Free Surface Variable Name'` + STOP 1, proving entry liveness. `freesurfacetomeshupdate_`/`_3` nm-proven but not `.sif`-reachable (USF protocol supplies no axis argument). USF core keeps the original's uninitialized-`SAVE` `FirstTime`/`t0` verbatim (serial zero-fill → branch never taken). Dead locals documented as skipped.
- Liveness: A/B scalars differ, both arms agree. `src/lib.rs` 299, accessor 1006. PROOF GREEN.

## dftm (`ports/dftm`, proof-dftm.sh)

- Unit `fortran:fem/src/modules/DFTSolver.F90` (1244 lines).
- Symbol surface: `T chargedensitysolver_`, `T poissonsolver_`, `T wavefunctionsolver_`. Diff empty. (`poissonsolver_` also exists in `Poisson.so` — no conflict: each `.so` is spliced independently under its own name.)
- Config A (exponential mixing + Perdew-Zunger + explicit weights) + rerun; config B (constant mixing + XC none + shifting constant + default weights — exercises constant-mix, XC-none, shift-PRINT, default-weights branches): byte-identical. Liveness asserts each branch message present in both arms, A≠B.
- All three LocalMatrix integrations + GRPULAY3 as transcribed Fortran module procedures over wrapper-side SAVE state (complex eigenvector arithmetic, DOT_PRODUCT/LUSolve, formatted I/O bit-identical by construction); Rust transcribes the three entries' control flow. Fixed en route: mass matrix pointer in `DefaultUpdateMass`. Documented deviations: `RF` SAVE-zero-init vs original stack garbage on never-assigned paths (unexercised — both configs always mix); skipped dead `EigsVar` lookup; CDS plain resize (original ALLOCATE unchecked).
- `src/lib.rs` 549, accessor 805. PROOF GREEN.

## emwavesolver (`ports/emwavesolver`, proof-emwavesolver.sh)

- Unit `fortran:fem/src/modules/EMWaveSolver.F90` (1249 lines).
- Symbol surface: `T emwavecalcfields_`, `T emwavecalcfields_init_`, `T emwavecalcfields_init0_`, `T emwavesolver_`, `T emwavesolver_init0_` (5/5, no `__*_MOD_*` symbols). Diff empty.
- Config A (2nd-order time + constant bulk + nodal Elfield) + rerun; config B (1st-order time, no bulk, 2 nonlinear iters, TEM outlet); config C (DG/elemental-field Init0 path + time derivatives): all byte-identical stdout/stderr/scalars.dat.
- `CONTIG` is a `-DCONTIG=` empty macro in the reference build (`elmerf90-helper.h`); build.rs passes it, accessor transcribes both `CONTIG` attributes verbatim. Defined-path divergence documented: reference fails its own iterative solve when the elemental-field path combines with first-order time and no constant bulk matrix — configs B/C isolate the two sides stably.
- Liveness: B, C differ from A on both arms. `src/lib.rs` 366, accessor 988. PROOF GREEN.

## emport (`ports/emport`, proof-emport.sh)

- Unit `fortran:fem/src/modules/EMPort.F90` (1591 lines).
- Symbol surface: `T emportpotential_`, `T emportpotential_init_`, `T emportsolver_`, `T emportsolver_init0_`. Diff empty.
- Config A (base eigen formulation, 30 ARPACK values, Umfpack, SaveScalars eigenvalues) + rerun; config B (potentials formulation + nodal post `EMPortPost` + Eigenfunction Index 2 + Port Ground presets): byte-identical.
- No defined-path divergence observed; unexercised paths documented (MaxPort≥1 filter / multi-matrix surgery need port-tagged 3-D meshes; PEs>1 guarded exactly as original). Proof elapsed 190 s (ARPACK solves).
- Liveness: A/B differ on both arms. `src/lib.rs` 298, accessor 1510. PROOF GREEN.

## isosurface (`ports/isosurface`, proof-isosurface.sh)

- Unit `fortran:fem/src/modules/Isosurface.F90` (2119 lines).
- Symbol surface: `T isosurfacesolver_`, `T isosurfacesolver_init_`. Diff empty.
- Config A (2-D single level) + rerun; config B (2-D fixed-surface remap + gmsh geo, 2 iters — asserts non-empty `mesh.geo` + `Remapping the fields` line, fixed-surface SAVE path on iter 2); config C (3-D isosurface + STL — asserts non-empty `mesh.stl`): all byte-identical stdout/stderr/scalars(+geo/stl).
- Entry control flow (phases, all four GOTO-100 exits, collective reduction, reports, savers, norm) in Rust; the twelve CONTAINS internals stay Fortran-side as whole units (their `TYPE(Element_t)`-by-value / assumed-shape dummies cannot cross the C boundary — the refusal class avoided by keeping them un-crossed). All `MPI_*` behind `ParEnv % PEs > 1` guards, same `libmpi_stubs`. SAVE mapping (`FixedSurface, SurfaceExist, Interpolant, InvPerm, IsoMesh, SaveParTmp, calls`) mirrored wrapper-side; other host state reset per call. No defined-path UB (only cross-call state explicitly SAVEd; `eperm`/temporaries fresh automatics reset in `iso_begin`).
- Liveness: A, B, C pairwise differ, both arms agree. `src/lib.rs` 145, accessor 1685, proof 409 lines. PROOF GREEN.

## facetshellsolve (`ports/facetshellsolve`, proof-facetshellsolve.sh)

- Unit `fortran:fem/src/modules/FacetShellSolve.F90` (2216 lines).
- Symbol surface: `T shellsolver_`, `T shellsolver_init0_`. Diff empty. Campaign form `Procedure = "FacetShellSolve" "ShellSolver"`.
- Config A (linear, triangles, stress) + rerun; config B (large deflection, quads, point load, 2 iters): byte-identical.
- Two real transcription slips caught by mechanical diff vs original (`Sdofs(5)` line in SMITC4 edge-41, `EpszzElement` source in averaging) — fixed, re-verified; all other procedures diff-identical modulo intentional `EXIT`→converged-flag split. `RelChange`/`LoadScale`/`NumberOfElementNodes` defined-path assumptions documented. Norm merged class + `Max Output Level = 4` both configs. Reference self-consistent on all proof paths.
- Liveness holds. `src/lib.rs` 135, accessor 2152. PROOF GREEN.

## coilsolver (`ports/coilsolver`, proof-coilsolver.sh; accessor `coil_acc.F90`)

- Unit `fortran:fem/src/modules/CoilSolver.F90` (2397 lines).
- Symbol surface: `T coilpotential_`, `T coilpotentialnormalized_`, `T coilsolver_`, `T coilsolver_init_`. Diff empty.
- Four clean-room differentials on the CoilSolver1 Curve mesh + separate-process rerun + liveness: A basic open-coil; B `Normalize Coil Current` + `Coil Conductivity Fix` 2nd-iteration + `GotElCond` path; C/D Poisson solve with `Source = Variable CoilPot` + `Real Procedure "CoilSolver" "CoilPotential[_Normalized]"` (both user-function entries per node incl Normalized's `n`-overwrite + per-element gradient cache). Poisson UDF configs use Direct/Umfpack (iterative BiCGStab diverges on this matrix in both arms); UDF Dirichlet merged into the Coil-Start BC (shared-target BCs map to one section). 1 steady iteration per config because the reference itself aborts on in-process re-entry (double `ALLOCATE` of entry work arrays) — re-entry proven port-vs-port per the bar.
- Entry-level control flow (coil loop CYCLE/EXIT, `TestCut` GOTO-100 bypass, `MaxNonlinIter`/`Part`/`dimi` loops, `CalcCurr`/`Fix Input Current Density` branches) in Rust; all data-touching work incl. fifteen CONTAINS internals verbatim in accessor (assembly, `SUM`/`MINVAL`, `WRITE(Message)`, `DSYEV` fit, `LUdecomp`/`LUSolve`); host-shared scratch `j`/`Found` made procedure-locals (clobbering provably unobservable). MPI/`ParallelReduction` transcribed verbatim under same `PEs > 1` guard.
- Liveness: B/C/D scalars differ from A in both arms (`coilpotential_` values C max 198.65 vs D max 208.55). `src/lib.rs` 248, accessor 1818. PROOF GREEN.

## advectiondiffusion (`ports/advectiondiffusion`, proof-advectiondiffusion.sh)

- Unit `fortran:fem/src/modules/AdvectionDiffusion.F90` (2488 lines).
- Symbol surface: `T advectiondiffusionsolver_`, `T advectiondiffusionsolver_init_`. Diff empty.
- Config A + rerun; config B (convection + stabilized + absolute mass + 2 iters): byte-identical.
- Rust owns entry/init + steady/nonlinear/element loop nest only; all 2488 lines of body logic (preamble, bulk, mixed solubility, boundary, solve tail, mass check, six CONTAINS internals verbatim incl. host-associated `Diffusivity`) statement-for-statement in wrappers. `Relax`/`EquationName`/first-`Norm` reads replicated and discarded; `ElementList`/`BackupForceVector` omitted as provably unused.
- Liveness: A/B differ, both arms agree. `src/lib.rs` 138, accessor 1782. PROOF GREEN.

## circuitsanddynamics (`ports/circuitsanddynamics`, proof-circuitsanddynamics.sh)

- Unit `fortran:fem/src/modules/CircuitsAndDynamics.F90` (2932 lines).
- Symbol surface: `T circuitsanddynamics_`, `T _init_`, `T circuitsanddynamicsharmonic_`, `T _init_`, `T circuitsoutput_`, `T _init_`. Diff empty.
- Config A (harmonic stranded coil + series resistor: harmonic driver/resistor branch/reporter) + rerun; config B (transient stranded + resistor over 5 BDF-1 steps: transient driver, timestep-refresh block, transient reporter) + cross-process rerun; liveness (Resistance 2.5→5.0 changes both arms identically): all byte-identical.
- Defect found and fixed by the differential: transcription dropped the `CM =>` association in lifted `hm_AddComponentEquationsAndCouplings` (while adding the `ASolver` host-ref fix) — unassociated pointer writes vanished silently on CRS (valgrind confirmed uninit read at `AddToMatrixElement:152`); one-line restoration → full identity. Two earlier bugs fixed the same way (stale `n_Circuits` read before `AllocateCircuitsList`; missing `Circuits(p) % n_comp` assignment). Massive/foil kernels verbatim Fortran under same compiler/flags (documented uncovered paths); MPI arm compiled, serial-guarded as reference.
- 1381 mesh, proof ≈ 15 s (8 Elmer runs). `src/lib.rs` 444, accessor 2859. PROOF GREEN (run twice).

## acoustics (`ports/acoustics`, proof-acoustics.sh)

- Unit `fortran:fem/src/modules/Acoustics.F90` (4735 lines, largest in lane).
- Symbol surface: `T acousticssolver_`, `T acousticssolver_init_`. Diff empty.
- Config A (driven cavity: BodyForce 1 + zero velocity Dirichlet, `Angular Frequency`, SaveScalars min/max on Flow, 1 iter) + rerun; config B (A + `Bubbles = Logical True` exercising Rust `LCondensate` + different frequency): byte-identical.
- Transcription split: driver + `LocalMatrix` family in Rust; BlockPreconditioning iterative suite, slip/interface/BEM/fluidic-force kernels verbatim-hidden in Fortran with identical statements/compiler. `IF(.FALSE.)` arms omitted (documented); `Bulk Viscosity` leading-space keyword kept verbatim.
- Two real bugs found and fixed via element-level differentials against an instrumented reference copy (scaffolding removed afterwards): Fortran col-major `dBasisdx` read with row-major indexing; complex-chain association order (`s*C` vs `s*bq` first, folded `s` in viscous terms).
- Liveness: A/B differ, both arms agree. `src/lib.rs` 1876, accessor 6121. PROOF GREEN.

## Campaign backstop

Fresh tree @ `9f6af2f85`, campaign configure (Debug,
`HWLOC_COMPONENTS=-gl`, `--allow-multiple-definition`), all 6 T021-A
ports + all 44 T022-A ports spliced (`<Name>.so` ← port cdylib, rebuilt
against the campaign `fmodules`; whole-`.so` replacement, the T020
recipe), full rebuild, `ctest -L quick`: **479/482, failures exactly
K=3 upstream segfaults #80/#536/#537** (log: lane `proof/`).
