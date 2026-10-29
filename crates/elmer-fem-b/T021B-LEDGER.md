# T021-B ledger (solver-entry plugins M–Z, lane rs-elmer-t021b)

| unit | status | entries | LOC (rs/acc) | wall |
|---|---|---|---|---|
| fortran:fem/src/modules/MeshChecksum.F90 | graded | meshchecksum_, meshchecksum_init_ | 115/128 | 135s |
| fortran:fem/src/modules/VelocityPrecond.F90 | graded | velocityprecond_, velocityprecond_init0_ | 101/76 | 79s |
| fortran:fem/src/modules/PressurePrecond.F90 | graded | pressureprecond_, pressureprecond_init0_ | 104/83 | 420s |
| fortran:fem/src/modules/SaveMesh.F90 | graded | savemesh_ | 87/85 | 480s |
| fortran:fem/src/modules/ReloadData.F90 | graded | reloadsolution_ | 164/107 | 230s |
| fortran:fem/src/modules/ReloadInput.F90 | graded | reloadinput_ | 79/74 | 340s |
| fortran:fem/src/modules/Probes.F90 | graded | volumeprobe_ | 155/236 | 370s |
| fortran:fem/src/modules/SpringAssembly.F90 | graded | springassembler_ | 270/171 | 450s |
| fortran:fem/src/modules/OdeSolver.F90 | graded | odesolver_, odesolver_init_ | 124/200 | 300s |
| fortran:fem/src/modules/PartitionMesh.F90 | graded | partitionmeshsolver_ | 111/191 | 690s |
| fortran:fem/src/modules/ScannedFieldSolver.F90 | graded | scannedfieldsolver_, scannedfieldsolver_init_ | 150/272 | 236s |
| fortran:fem/src/modules/ModelPDEhandle.F90 | graded | advdiffsolver_ | 367/341 | 2100s |
| fortran:fem/src/modules/ModelPDE.F90 | graded | advdiffsolver_ | 462/351 | 500s |
| fortran:fem/src/modules/ScalarPotentialSolver.F90 | graded | scalarpotentialsolver_, scalarpotentialsolver_init_ | 347/334 | 450s |
| fortran:fem/src/modules/WVectorFix.F90 | graded | wvectorfix_, wvectorfix_init_ | 137/224 | 360s |
| fortran:fem/src/modules/Mesh2MeshSolver.F90 | graded | mesh2meshsolver_, mesh2meshsolver_init0_ | 115/285 | 150s |
| fortran:fem/src/modules/ShearrateSolver.F90 | graded | shearratesolver_, shearratesolver_init_ | 215/365 | 210s |
| fortran:fem/src/modules/VorticitySolver.F90 | graded | vorticitysolver_, vorticitysolver_init_ | 394/445 | 390s |
| fortran:fem/src/modules/WaveSolver.F90 | graded | wavesolver_, wavesolver_init_ | 418/375 | 1010s |
| fortran:fem/src/modules/ModelPDECutFEM.F90 | graded | advdiffsolver_, advdiffsolver_init_ | 646/484 | 570s |
| fortran:fem/src/modules/ShallowWaterNS.F90 | graded | shallowwaternssolver_ | 208/444 | 1240s |
| fortran:fem/src/modules/NormalSolver.F90 | graded | normalsolver_, normalsolver_init_ | 228/490 | 307s |
| fortran:fem/src/modules/SunAngle.F90 | graded | sunanglesolver_, sunanglesolver_init_ | 595/308 | 540s |
| fortran:fem/src/modules/PoissonDG.F90 | graded | poissonsolver_ | 512/443 | 470s |
| fortran:fem/src/modules/RateOfChange.F90 | graded | rateofchangesolver_ | 226/328 | 855s |
| fortran:fem/src/modules/MeshSolve.F90 | graded | meshsolver_, meshsolver_init_ | 563/610 | 770s |
| fortran:fem/src/modules/OpenFoam2ElmerIO.F90 | graded | openfoam2elmerfit_ | 427/682 | 670s |
| fortran:fem/src/modules/NodeToEdgeSolver.F90 | graded | nodetoedgefield_, nodetoedgefield_init0_, extrudedrestart_, extrudedrestart_init0_ | 276/580 | 1230s |
| fortran:fem/src/modules/VectorHelmholtzLumping.F90 | graded | vectorhelmholtzlumping_ | 445/662 | 1600s |
| fortran:fem/src/modules/WPotentialSolver.F90 | graded | wsolve_, wsolve_dummy_, wsolve_init0_ | 401/649 | 1270s |
| fortran:fem/src/modules/StructuredFlowLine.F90 | graded | structuredflowline_, structuredflowline_init_ | 497/697 | 730s |
| fortran:fem/src/modules/ThermoElectricSolver.F90 | graded | thermoelectricsolver_, thermoelectricsolver_init_ | 300/604 | 2030s |
| fortran:fem/src/modules/RigidMeshMapper.F90 | graded | rigidmeshmapper_ | 647/1019 | 960s |
| fortran:fem/src/modules/NonphysicalMeshSolve.F90 | graded | meshsolver_ | 599/768 | 780s |
| fortran:fem/src/modules/ModelMixedPoisson.F90 | graded | mixedpoisson_, mixedpoisson_init0_, mixedpoisson_post_ | 286/694 | 620s |
| fortran:fem/src/modules/MarchingODESolver.F90 | graded | marchingodesolver_, marchingodesolver_init_ | 161/859 | 600s |
| fortran:fem/src/modules/VectorHelmholtzNodal.F90 | graded | vectorhelmholtznodal_, vectorhelmholtznodal_init_ | 257/926 | 1000s |
| fortran:fem/src/modules/StructuredMeshMapper.F90 | graded | structuredmeshmapper_, structuredmeshmapper_init_ | 815/973 | 1320s |
| fortran:fem/src/modules/TransientPhaseChange.F90 | graded | transientphasechange_, transientphasechange_init_ | 438/935 | 930s |
| fortran:fem/src/modules/StructuredProjectToPlane.F90 | graded | structuredprojecttoplane_, structuredprojecttoplane_init_ | 101/1103 | 540s |
| fortran:fem/src/modules/PoissonBEM.F90 | graded | poissonbemsolver_, matvec_, precond_, B __globmat_MOD_matrix | 167/597 | 3130s |
| fortran:fem/src/modules/Smitc.F90 | graded | smitcsolver_, smitcsolver_init0_ | 465/801 | 660s |
| fortran:fem/src/modules/PoissonBoltzmannSolve.F90 | graded | poissonboltzmannsolve_ | 152/713 | 520s |
| fortran:fem/src/modules/SteadyPhaseChange.F90 | graded | steadyphasechange_, meltingheat_ | 231/901 | 970s |
| fortran:fem/src/modules/TransportEquation.F90 | graded | transportequationsolver_ | 507/620 | 1930s |
| fortran:fem/src/modules/StatMagSolve.F90 | graded | statmagsolver_, statmagsolver_init_ | 468/730 | 820s |
| fortran:fem/src/modules/ParticleAdvector.F90 | graded | particleadvector_, particleadvector_init_ | 278/1137 | 1100s |
| fortran:fem/src/modules/TopoOpt.F90 | graded | topoopt_, topoopt_init_, topoopt_init0_ | 189/1486 | 830s |
| fortran:fem/src/modules/PhaseChangeSolve.F90 | graded | phasechangesolve_, phasechangesolve_init_, meltingheat_ | 279/1277 | 1380s |
| fortran:fem/src/modules/RigidBodyReduction.F90 | graded | rigidbody_ (Before-Linsolve hook, INTEGER FUNCTION) | 166/1787 | 1050s |
| fortran:fem/src/modules/ReynoldsSolver.F90 | graded | reynoldssolver_, reynoldssolver_init_, reynoldspostprocess_, reynoldspostprocess_init_ | 254/1404 | 1160s |
| fortran:fem/src/modules/ParticleDynamics.F90 | graded | particledynamics_, particledynamics_init_, 3 module procs | 1228/1821 | 2290s |
| fortran:fem/src/modules/VectorHelmholtz.F90 | graded | 7 entries (solver + dummy + calcfields families) | 296/2458 | 1610s |
| fortran:fem/src/modules/Stokes.F90 | graded | stokessolver_ | 1525/1009 | 2270s |
| fortran:fem/src/modules/MagnetoDynamics2D.F90 | graded | 7 entries (magnetodynamics2d[_init] + harmonic[_init,_init0] + bsolver[_init]) | 333/4200 | 1010s |
| fortran:fem/src/modules/ShellSolver.F90 | graded | shellsolver_, shellsolver_init0_ | 525/6758 | 1750s |

Proofs: `ports/<lower>/proof-<lower>.sh` (nm equality + >=2 .sif differentials + rerun + liveness).
Independently re-verified by lane runner (proof/verify-*.log in lane dir).
Deviations documented in UNITS.md notes (lane dir) and port code comments.
Poisson.F90 (in range) already graded in T020; not redone.
Refused: 14 PluginEntry-gate (see regrade evidence) + 3 triage (ParallelRemesh, SimpleCircuits, ProjectToPlane) with reasons in lane UNITS.md.
