# femcore lane B unit ledger (claim B: 176 of the 353 refused/skipped fem core units)
# generated mechanically from the relaxed-gate re-grade (`data/elmer-lanes/femcore-regrade-b.json` in rustsmith);
# the reasons are the gate's CURRENT reasons, not the pre-relaxation ledger's.
# status: ported = Rust port in crates/elmer-fem-a/ports + byte-identical corpus;
#   skipped-* = recorded skip, never forced.
# The shared lane-fem-a/b ledgers are left untouched on purpose: this lane owns
# only its two port crates and this ledger file (runbook section 3 conflict rule).

| unit | status | reason |
|---|---|---|
| fortran:fem/src/modules/EliminateDirichlet.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/EliminatePeriodic.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/Elmer2OpenFoamIO.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/EnergyRelease.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/FDiffusion.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/FDiffusion3D.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/FacetShellSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (16); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/FilmCompressibility.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/FilmFlowSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4); procedure-not-found-in-gate-parse (3) |
| fortran:fem/src/modules/FilterTimeSeries.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/FindOptimum.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (9); derived-type dummy (needs whole-program layout) (2); non-intrinsic dummy (2) |
| fortran:fem/src/modules/FixTangentVelo.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/FlowSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/FluidicForce.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/FluxSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/FourierLoss.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/FreeSurfaceReduced.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/FreeSurfaceSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (7); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/GmshOutputReader.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/HarmonicSignals.F90#source | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/HeatSolve.F90#heatsolvefront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/HeatSolve.F90#heatsolver_init0 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (20); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/HeatSolveLegacy.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (14); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/HeatSolveVec.F90#heatsolver_init0 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/HeatSolveVec.F90#heatsolvevecfront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (2); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/HelmholtzBEM.F90#globmatc | skipped-vacuous-admit | gate admits with 0 non-BIND(C) Fortran exports: the scope owns no portable procedure (ADR-047 precedent) |
| fortran:fem/src/modules/HelmholtzBEM.F90#helmholtzbemsolver | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (9); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/HelmholtzProjection.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (3) |
| fortran:fem/src/modules/HelmholtzSolve.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (3) |
| fortran:fem/src/modules/HydrostaticNS.F90#hydrostaticnssolver_init0 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/HydrostaticNS.F90#hydrostaticnsutils | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); undeclared function result (implicit typing) (1); assumed-shape rank > 2 (1) |
| fortran:fem/src/modules/HydrostaticNSVec.F90#hydrostaticnssolver_init0 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/HydrostaticNSVec.F90#hydrostaticnsvecfront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1) |
| fortran:fem/src/modules/IncompressibleNS.F90#incompressiblelocalforms | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); undeclared function result (implicit typing) (1); assumed-shape rank > 2 (1); optional dummy (hidden presence arg) (1); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/IncompressibleNS.F90#incompressiblenssolver_init0 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/IncompressibleNSVec.F90#incompressiblenssolver_init0 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/IncompressibleNSVec.F90#incompressiblensvecfront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1) |
| fortran:fem/src/modules/Isosurface.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (15); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/LevelSet/LevelSetCurvature.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/LevelSet/LevelSetDistance.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/LevelSet/LevelSetIntegrate.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/LevelSet/LevelSetSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/LevelSet/LevelSetTimestep.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/MagneticSolve/MagneticSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/MagneticSolve/Maxwell.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/MagneticSolve/MaxwellAxiS.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/MagneticSolve/MaxwellGeneral.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/MagnetoDynamics/CalcFields.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (9); derived-type dummy (needs whole-program layout) (4); non-intrinsic dummy (2) |
| fortran:fem/src/modules/MagnetoDynamics/JfixPotentialSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/MagnetoDynamics/PlatiseFerroModel.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/MagnetoDynamics/Utils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (16) |
| fortran:fem/src/modules/MagnetoDynamics/WhitneyAVHarmonicSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (22); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/MagnetoDynamics/WhitneyAVSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (23); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/MagnetoDynamics2D.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (21); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/MarchingODESolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/Mesh2MeshSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/MeshChecksum.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/MeshSolve.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/ModelMixedPoisson.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/ModelPDE.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/ModelPDECutFEM.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/ModelPDEhandle.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/NodeToEdgeSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (5); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/NonphysicalMeshSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/NormalSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/OdeSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/OpenFoam2ElmerIO.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ParStokes.F90#anisotropicmaterialmodels | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/ParStokes.F90#stokessolver_init0 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (9); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/ParallelRemesh.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ParticleAdvector.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (6); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ParticleDynamics.F90#particledynamics | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ParticleDynamics.F90#particledynamicsstuff | skipped-gate-admitted-not-portable | gate admits the dummy shape (T017 shape B), but the scope body reads Elmer whole-program state (ValueList_t Params via ListGet*/ListFind, Solver_t/Model_t beside it); a file-level port must define all 4 symbols, so a Rust port would re-implement the particle core instead of this unit (ADR-051 section 4 precedent; T017b judgment) |
| fortran:fem/src/modules/PartitionMesh.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/PhaseChangeSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/Poisson.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/PoissonBEM.F90#globmat | skipped-vacuous-admit | gate admits with 0 non-BIND(C) Fortran exports: the scope owns no portable procedure (ADR-047 precedent) |
| fortran:fem/src/modules/PoissonBEM.F90#poissonbemsolver | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (8); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/PoissonBoltzmannSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/PoissonDG.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/PressurePrecond.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/Probes.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ProjectToPlane.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (10); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/RateOfChange.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/ReloadData.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ReloadInput.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ResultOutputSolve/DXOutputSolver.F90#dxfile | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); character dummy (hidden length param) (3) |
| fortran:fem/src/modules/ResultOutputSolve/DXOutputSolver.F90#dxoutputsolver | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/ResultOutputSolve/ElmerPostOutputSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ResultOutputSolve/GidOutputSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ResultOutputSolve/ResultOutputSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/ResultOutputSolve/VtkOutputSolver.F90#vtklegacyfile | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4); character dummy (hidden length param) (2) |
| fortran:fem/src/modules/ResultOutputSolve/VtkOutputSolver.F90#vtkoutputsolver | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/ResultOutputSolve/VtuOutputSolver.F90#vtuoutputsolver | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ResultOutputSolve/VtuOutputSolver.F90#vtuxmlfile | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); character dummy (hidden length param) (1) |
| fortran:fem/src/modules/ReynoldsSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (8); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/RichardsSolver.F90#porousmaterials | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/RichardsSolver.F90#richardssolver_init | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (6); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/RigidBodyReduction.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (6); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/RigidMeshMapper.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/SaveData/SaveBoundaryValues.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/SaveData/SaveDependence.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/SaveData/SaveLine.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (12); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/SaveData/SaveMaterials.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/SaveData/SaveProjection.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/SaveData/SaveScalars.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (11); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/SaveGridData.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/SaveMesh.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ScalarPotentialSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/ScannedFieldSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/ShallowWaterNS.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ShearrateSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/ShellSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (30); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/SimpleCircuits.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (8); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/Smitc.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/SpringAssembly.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/StatCurrentSolve.F90#statcurrentsolvefront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/StatCurrentSolve.F90#statcurrentsolver_init | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (10); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/StatCurrentSolveLegacy.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (6); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/StatCurrentSolveVec.F90#statcurrentsolver_init | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/StatCurrentSolveVec.F90#statcurrentsolvevecfront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (2); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/StatElecBoundary.F90#statelecboundaryenergy | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (5) |
| fortran:fem/src/modules/StatElecBoundary.F90#statelecboundaryutils | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/StatElecSolve.F90#statelecsolvefront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/StatElecSolve.F90#statelecsolver_init0 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (13); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/StatElecSolveLegacy.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (11); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/StatElecSolveVec.F90#statelecsolver_init0 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (5) |
| fortran:fem/src/modules/StatElecSolveVec.F90#statelecsolvevecfront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/StatMagSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/SteadyPhaseChange.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/Stokes.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (15); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/StreamSolver.F90 | skipped-gate-refused | gate refuses: non-intrinsic dummy (2); procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/StressSolve.F90#stresssolvefront | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); character dummy (hidden length param) (1) |
| fortran:fem/src/modules/StressSolve.F90#stresssolver_init0 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (9); derived-type dummy (needs whole-program layout) (7) |
| fortran:fem/src/modules/StructuredFlowLine.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/StructuredMeshMapper.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/StructuredProjectToPlane.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/SunAngle.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (6); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/ThermoElectricSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/TopoOpt.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (13); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/TransientCost.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); non-intrinsic dummy (1) |
| fortran:fem/src/modules/TransientPhaseChange.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/TransportEquation.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/TurbulenceModels/KESolver.F90#kesolver | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/TurbulenceModels/KESolver.F90#kesolverfront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/TurbulenceModels/KESolver.F90#kesolverlocalforms | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/TurbulenceModels/KESolverLegacy.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/TurbulenceModels/Komega.F90#komega | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/TurbulenceModels/Komega.F90#komegafront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/TurbulenceModels/Komega.F90#komegalocalforms | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/TurbulenceModels/KomegaLegacy.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/TurbulenceModels/SSTKomega.F90#sstkomega | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/TurbulenceModels/SSTKomega.F90#sstkomegafront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/TurbulenceModels/SSTKomega.F90#sstkomegalocalforms | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (5) |
| fortran:fem/src/modules/TurbulenceModels/SSTKomegaLegacy.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/TurbulenceModels/Spalart-Allmaras.F90#spalartallmaras | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/TurbulenceModels/Spalart-Allmaras.F90#spalartallmarasfront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/TurbulenceModels/Spalart-Allmaras.F90#spalartallmaraslocalforms | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/TurbulenceModels/Spalart-AllmarasLegacy.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/TurbulenceModels/V2FSolver.F90#v2f_ldm | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/TurbulenceModels/V2FSolver.F90#v2fsolverfront | skipped-gate-refused | gate refuses: character dummy (hidden length param) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/modules/TurbulenceModels/V2FSolver.F90#v2fsolverlocalforms | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/TurbulenceModels/V2FSolverLegacy.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/UMATLib.F90 | skipped-gate-refused | gate refuses: character dummy (hidden length param) (5) |
| fortran:fem/src/modules/VectorHelmholtz.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (15); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/VectorHelmholtzLumping.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/VectorHelmholtzNodal.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/VelocityPrecond.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/VorticitySolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/WPotentialSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/WVectorFix.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/WaveSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/mpif_stubs.F90 | ported | gate admits (T017 shape A: F77 global external subprograms + module variable); port crate `mpistubs`, driver corpus byte-identical (stdout+stderr+exit), zero STRONG port-level dups beyond the 40 replaced symbols; compiled into libmpi_stubs.so, which the default build (WITH_MPI=TRUE) does not link into elmersolver (recorded, not hidden) |
| fortran:fem/src/solve_cmplx.F90 | skipped-gate-refused | gate refuses: undeclared dummy (implicit typing) (2) |
| fortran:fem/src/solve_real.F90 | ported | gate admits (T017 shape A: F77 global external subprogram); port crate `solvereal`, 5-case corpus byte-identical (stdout+stderr+exit), zero STRONG port-level dups beyond the 2 replaced symbols |
| python:fem/src/split_ElementDescription.py | skipped-helper | gate admits vacuously (no Fortran export); build-time code generator: links, never ports (matc T005 precedent) |
