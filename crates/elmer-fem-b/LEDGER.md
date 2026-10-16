# lane-fem-b scope ledger (218 scopes in 184 files)
# generated mechanically from the T014-relaxed gate run 2026-10-15 (gate1.json) + store.db; reasons reviewed against gate internals.
# status: graded = Rust port in crates/elmer-fem-b/ports + differential green; skipped-* = recorded skip, never forced.

| scope unit | status | reason |
|---|---|---|
| fortran:fem/src/modules/ADIOS2OutputSolver.F90#adios2outputsolver_init | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ADIOS2OutputSolver.F90#adiosoutputsolverutils | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/Acoustics.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/AdvectionDiffusion.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/AdvectionReaction.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/AllocateSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/AmsPrecSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ArteryOutlet.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ArtificialCompressibility.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/BatterySolver/BatteryPost.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/BatterySolver/BatteryUtils.F90 | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/BatterySolver/ElectrolyteCons.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/BatterySolver/ElectrolytePot.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/BatterySolver/SolidPhaseCons.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/BatterySolver/SolidPhasePot.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/BeamSolver3D.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/CahnHilliard.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/CircuitsAndDynamics.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/CoilSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/CompressibleNS.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ConsistentSplitting/CompressibilitySolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ConsistentSplitting/OptimalSolutionUpdate.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ConsistentSplitting/PressureSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ConsistentSplitting/VelocitySolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/CoordinateTransform.F90#rotmsolver_init | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/CoordinateTransform.F90#vecutils | skipped-explicit-shape | explicit-shape dummy (later shape; T017) |
| fortran:fem/src/modules/CraigBamptonSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/CyclicConvergence.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/DCRComplexSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/DFTSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/DataToFieldSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/DirectionSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/DistanceSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/DistributeSource.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/DivProjection.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/DivergenceSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/EMPort.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/EMWaveSolver.F90#emwavesolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/EMWaveSolver.F90#emwavesolverutils | skipped | gate admits (no non-BIND(C) procedures) but scope owns no portable procedure: empty module or data-only module (module ALLOCATABLE descriptor symbol unproven shape, T017) |
| fortran:fem/src/modules/EdgeFEInterpolationTest.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ElasticSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ElectricForce.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/Electrokinetics.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ElementSizeSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ElementStats.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/EliminateDirichlet.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/EliminatePeriodic.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/Elmer2OpenFoamIO.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/EnergyRelease.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FDiffusion.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FDiffusion3D.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FacetShellSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FilmCompressibility.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FilmFlowSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FilterTimeSeries.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FindOptimum.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FixTangentVelo.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FlowSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FluidicForce.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FluxSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FourierLoss.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FreeSurfaceReduced.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/FreeSurfaceSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/GmshOutputReader.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HarmonicSignals.F90#harmutils | graded | 11-case differential byte-identical (2 held out), strong-symbol gate, one-provider splice; EVIDENCE.md |
| fortran:fem/src/modules/HarmonicSignals.F90#source | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HeatSolve.F90#heatsolvefront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/HeatSolve.F90#heatsolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HeatSolveLegacy.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HeatSolveVec.F90#heatsolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HeatSolveVec.F90#heatsolvevecfront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/HelmholtzBEM.F90#globmatc | skipped | gate admits (no non-BIND(C) procedures) but scope owns no portable procedure: empty module or data-only module (module ALLOCATABLE descriptor symbol unproven shape, T017) |
| fortran:fem/src/modules/HelmholtzBEM.F90#helmholtzbemsolver | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HelmholtzProjection.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HelmholtzSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HydrostaticNS.F90#hydrostaticnssolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HydrostaticNS.F90#hydrostaticnsutils | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/HydrostaticNSVec.F90#hydrostaticnssolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/HydrostaticNSVec.F90#hydrostaticnsvecfront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/IncompressibleNS.F90#incompressiblelocalforms | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/IncompressibleNS.F90#incompressiblenssolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/IncompressibleNSVec.F90#incompressiblenssolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/IncompressibleNSVec.F90#incompressiblensvecfront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/Isosurface.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/LevelSet/LevelSetCurvature.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/LevelSet/LevelSetDistance.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/LevelSet/LevelSetIntegrate.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/LevelSet/LevelSetSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/LevelSet/LevelSetTimestep.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MagneticSolve/MagneticSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MagneticSolve/Maxwell.F90 | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/MagneticSolve/MaxwellAxiS.F90 | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/MagneticSolve/MaxwellGeneral.F90 | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/MagnetoDynamics/CalcFields.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MagnetoDynamics/JfixPotentialSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MagnetoDynamics/PlatiseFerroModel.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MagnetoDynamics/Utils.F90 | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/MagnetoDynamics/WhitneyAVHarmonicSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MagnetoDynamics/WhitneyAVSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MagnetoDynamics2D.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MarchingODESolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/Mesh2MeshSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MeshChecksum.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/MeshSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ModelMixedPoisson.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ModelPDE.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ModelPDECutFEM.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ModelPDEhandle.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/NodeToEdgeSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/NonphysicalMeshSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/NormalSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/OdeSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/OpenFoam2ElmerIO.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ParStokes.F90#anisotropicmaterialmodels | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/ParStokes.F90#stokessolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ParallelRemesh.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ParticleAdvector.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ParticleDynamics.F90#particledynamics | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ParticleDynamics.F90#particledynamicsstuff | skipped-explicit-shape | explicit-shape dummy (later shape; T017) |
| fortran:fem/src/modules/PartitionMesh.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/PhaseChangeSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/Poisson.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/PoissonBEM.F90#globmat | skipped | gate admits (no non-BIND(C) procedures) but scope owns no portable procedure: empty module or data-only module (module ALLOCATABLE descriptor symbol unproven shape, T017) |
| fortran:fem/src/modules/PoissonBEM.F90#poissonbemsolver | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/PoissonBoltzmannSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/PoissonDG.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/PressurePrecond.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/Probes.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ProjectToPlane.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/RateOfChange.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ReloadData.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ReloadInput.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ResultOutputSolve/DXOutputSolver.F90#dxfile | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/ResultOutputSolve/DXOutputSolver.F90#dxoutputsolver | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ResultOutputSolve/ElmerPostOutputSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ResultOutputSolve/GidOutputSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ResultOutputSolve/ResultOutputSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ResultOutputSolve/VtkOutputSolver.F90#vtklegacyfile | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/ResultOutputSolve/VtkOutputSolver.F90#vtkoutputsolver | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ResultOutputSolve/VtuOutputSolver.F90#vtuoutputsolver | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ResultOutputSolve/VtuOutputSolver.F90#vtuxmlfile | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/ReynoldsSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/RichardsSolver.F90#porousmaterials | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/RichardsSolver.F90#richardssolver_init | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/RigidBodyReduction.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/RigidMeshMapper.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SaveData/SaveBoundaryValues.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SaveData/SaveDependence.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SaveData/SaveLine.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SaveData/SaveMaterials.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SaveData/SaveProjection.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SaveData/SaveScalars.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SaveGridData.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SaveMesh.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ScalarPotentialSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ScannedFieldSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ShallowWaterNS.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ShearrateSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ShellSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SimpleCircuits.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/Smitc.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SpringAssembly.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StatCurrentSolve.F90#statcurrentsolvefront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/StatCurrentSolve.F90#statcurrentsolver_init | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StatCurrentSolveLegacy.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StatCurrentSolveVec.F90#statcurrentsolver_init | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StatCurrentSolveVec.F90#statcurrentsolvevecfront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/StatElecBoundary.F90#statelecboundaryenergy | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StatElecBoundary.F90#statelecboundaryutils | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/StatElecSolve.F90#statelecsolvefront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/StatElecSolve.F90#statelecsolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StatElecSolveLegacy.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StatElecSolveVec.F90#statelecsolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StatElecSolveVec.F90#statelecsolvevecfront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/StatMagSolve.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SteadyPhaseChange.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/Stokes.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StreamSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StressSolve.F90#stresssolvefront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/StressSolve.F90#stresssolver_init0 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StructuredFlowLine.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StructuredMeshMapper.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/StructuredProjectToPlane.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/SunAngle.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/ThermoElectricSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TopoOpt.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TransientCost.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TransientPhaseChange.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TransportEquation.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/KESolver.F90#kesolver | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/KESolver.F90#kesolverfront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/KESolver.F90#kesolverlocalforms | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/KESolverLegacy.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/Komega.F90#komega | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/Komega.F90#komegafront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/Komega.F90#komegalocalforms | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/KomegaLegacy.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/SSTKomega.F90#sstkomega | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/SSTKomega.F90#sstkomegafront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/SSTKomega.F90#sstkomegalocalforms | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/SSTKomegaLegacy.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/Spalart-Allmaras.F90#spalartallmaras | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/Spalart-Allmaras.F90#spalartallmarasfront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/Spalart-Allmaras.F90#spalartallmaraslocalforms | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/Spalart-AllmarasLegacy.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/V2FSolver.F90#v2f_ldm | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/TurbulenceModels/V2FSolver.F90#v2fsolverfront | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/V2FSolver.F90#v2fsolverlocalforms | skipped-derived-dummy | module procedure with derived-type dummy (needs whole-program layout; T017) |
| fortran:fem/src/modules/TurbulenceModels/V2FSolverLegacy.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/UMATLib.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/VectorHelmholtz.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/VectorHelmholtzLumping.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/VectorHelmholtzNodal.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/VelocityPrecond.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/VorticitySolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/WPotentialSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/WVectorFix.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
| fortran:fem/src/modules/WaveSolver.F90 | skipped-external-scope | external (F77-global) subroutine scope: no MODULE, synthetic MOD linkage unprovable by gate parse (T017) |
