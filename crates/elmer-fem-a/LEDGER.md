# lane-fem-a unit ledger (168 units)
# generated mechanically from lanes.json + frozen recon + source scan; reasons reviewed by hand.
# status: ported = Rust port in crates/elmer-fem-a/ports + local differential green (grade pending backstop);
# refused-* = recorded skip, never forced.

| unit | status | reason |
|---|---|---|
| c:fem/obsolete/SC2Elmer.c | skipped-not-in-build | fem/obsolete has no CMakeLists: never compiled |
| c:fem/src/CPUTime.c | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| c:fem/src/GFortranLocaleFix.c | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| c:fem/src/Load.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/SolveHypre.c | skipped-ext-dep | needs absent optional solver lib (HYPRE/SuperLU/AMGX/CHOLMOD/ROCALUTION all FALSE); stub bodies unlinked |
| c:fem/src/SolveSuperLU.c | skipped-ext-dep | needs absent optional solver lib (HYPRE/SuperLU/AMGX/CHOLMOD/ROCALUTION all FALSE); stub bodies unlinked |
| c:fem/src/amgx.c | skipped-ext-dep | needs absent optional solver lib (HYPRE/SuperLU/AMGX/CHOLMOD/ROCALUTION all FALSE); stub bodies unlinked |
| c:fem/src/binio/binio.c | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| c:fem/src/cholmod.c | skipped-ext-dep | needs absent optional solver lib (HYPRE/SuperLU/AMGX/CHOLMOD/ROCALUTION all FALSE); stub bodies unlinked |
| c:fem/src/elmer_lua_iface.c | skipped-not-in-build | WITH_LUA=FALSE: excluded from solverlib_SOURCES, never compiled |
| c:fem/src/elmer_lua_iface.h | skipped-helper | header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/elmerf90-helper.in.h | skipped-helper | header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/elmerf90.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/elmerld.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/f_stubs.c | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| c:fem/src/fft.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/fft_test.c | skipped-not-in-build | no CMake target (test main never built) |
| c:fem/src/mpif_stub.h | skipped-helper | header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/umf4_f77wrapper.c | ported | batch2: local differential byte-identical (di path vs vendored UMFPACK); machine grade pending |
| c:fem/src/view3d/16node_quad.c | ported | batch2: local differential byte-identical; machine grade pending |
| c:fem/src/view3d/2d_4node.c | ported | batch2: local differential byte-identical (libc fprintf %g paths); machine grade pending |
| c:fem/src/view3d/BiCubicUtil.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/view3d/BiLinearUtil.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/view3d/BiQuadraticUtil.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/view3d/ClipShadow.c | ported | batch2: local differential byte-identical; machine grade pending |
| c:fem/src/view3d/ContourInteg.c | ported | batch2: local differential byte-identical (back-face/0.5*ABS paths); machine grade pending |
| c:fem/src/view3d/Ipoints.h | skipped-helper | header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/view3d/LUDecomp.c | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| c:fem/src/view3d/LineUtil.c | ported | batch2: local differential byte-identical (dead-division preserved); machine grade pending |
| c:fem/src/view3d/LinearSolve.c | skipped-not-in-build | orphan source: no CMake target compiles it (verified: unreferenced in any CMakeLists) |
| c:fem/src/view3d/LinearSolveConjGrad.c | skipped-not-in-build | orphan source: no CMake target compiles it (verified: unreferenced in any CMakeLists) |
| c:fem/src/view3d/LinearSolveGaussSeidel.c | skipped-not-in-build | orphan source: no CMake target compiles it (verified: unreferenced in any CMakeLists) |
| c:fem/src/view3d/LinearSolveJacob.c | skipped-not-in-build | orphan source: no CMake target compiles it (verified: unreferenced in any CMakeLists) |
| c:fem/src/view3d/RayTrace.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/view3d/ShaftCull.c | ported | batch2: local differential byte-identical (driver-UB fixed, Walk/planes verified); machine grade pending |
| c:fem/src/view3d/TestModel.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/view3d/TriangleUtil.c | ported | batch2: local differential byte-identical (pair+subdivide paths); machine grade pending |
| c:fem/src/view3d/VectorUtil.c | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| c:fem/src/view3d/ViewFactors.c | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| c:fem/src/view3d/ViewFactors.h | skipped-helper | header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/view3d/second.c | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| c:fem/src/viewaxis/viewfact.h | skipped-helper | header-only unit: links, never ports (matc T005 precedent) |
| cxx:fem/src/rocalution.cpp | skipped-ext-dep | needs absent optional solver lib (HYPRE/SuperLU/AMGX/CHOLMOD/ROCALUTION all FALSE); stub bodies unlinked |
| cxx:fem/src/viewaxis/viewfact.cpp | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| fortran:fem/examples/turing/Turing.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/obsolete/MovingElstatSolver.F90 | skipped-not-in-build | fem/obsolete has no CMakeLists: never compiled (also non-BIND(C)) |
| fortran:fem/src/ADIOS2Utils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (17); no T014 relaxation landed |
| fortran:fem/src/Adaptive.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (26); no T014 relaxation landed |
| fortran:fem/src/AddrFunc.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (2); no T014 relaxation landed |
| fortran:fem/src/BackwardError.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (6); no T014 relaxation landed |
| fortran:fem/src/BandMatrix.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (12); no T014 relaxation landed |
| fortran:fem/src/BandwidthOptimize.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (6); no T014 relaxation landed |
| fortran:fem/src/BlockSolve.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (38); no T014 relaxation landed |
| fortran:fem/src/BoundaryConditionUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (26); no T014 relaxation landed |
| fortran:fem/src/CRSMatrix.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (77); no T014 relaxation landed |
| fortran:fem/src/CircuitUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (83); no T014 relaxation landed |
| fortran:fem/src/ClusteringMethods.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (11); no T014 relaxation landed |
| fortran:fem/src/Constitutive.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (14); no T014 relaxation landed |
| fortran:fem/src/ContactUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (20); no T014 relaxation landed |
| fortran:fem/src/CoordinateSystems.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (13); no T014 relaxation landed |
| fortran:fem/src/CutFEMUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (31); no T014 relaxation landed |
| fortran:fem/src/DefUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (144); no T014 relaxation landed |
| fortran:fem/src/Differentials.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (9); no T014 relaxation landed |
| fortran:fem/src/DiffuseConvectiveAnisotropic.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/src/DiffuseConvectiveGeneralAnisotropic.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/src/DirectSolve.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (42); no T014 relaxation landed |
| fortran:fem/src/EigenSolve.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (15); no T014 relaxation landed |
| fortran:fem/src/ElemInfo.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (38); no T014 relaxation landed |
| fortran:fem/src/ElementBasis.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (30); no T014 relaxation landed |
| fortran:fem/src/ElementDescription.F90 | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| fortran:fem/src/ElementGeometry.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (21); no T014 relaxation landed |
| fortran:fem/src/ElementUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (37); no T014 relaxation landed |
| fortran:fem/src/ElmerSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (24); no T014 relaxation landed |
| fortran:fem/src/ExchangeCorrelations.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (7); no T014 relaxation landed |
| fortran:fem/src/Feti.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (31); no T014 relaxation landed |
| fortran:fem/src/FreeSurface.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (4); no T014 relaxation landed |
| fortran:fem/src/GeneralUtils.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/GeometryFitting.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (7); no T014 relaxation landed |
| fortran:fem/src/H1Basis.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (104); no T014 relaxation landed |
| fortran:fem/src/HashTable.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (13); no T014 relaxation landed |
| fortran:fem/src/Integration.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (36); no T014 relaxation landed |
| fortran:fem/src/InterpVarToVar.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (5); no T014 relaxation landed |
| fortran:fem/src/InterpolateMeshToMesh.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (11); no T014 relaxation landed |
| fortran:fem/src/Interpolation.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (14); no T014 relaxation landed |
| fortran:fem/src/IpFieldInterface.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (2); no T014 relaxation landed |
| fortran:fem/src/IterSolve.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (21); no T014 relaxation landed |
| fortran:fem/src/IterativeMethods.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (27); no T014 relaxation landed |
| fortran:fem/src/LinearAlgebra.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (25); no T014 relaxation landed |
| fortran:fem/src/LinearForms.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (8); no T014 relaxation landed |
| fortran:fem/src/ListMatrix.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (29); no T014 relaxation landed |
| fortran:fem/src/ListMatrixArray.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (18); no T014 relaxation landed |
| fortran:fem/src/Lists.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (176); no T014 relaxation landed |
| fortran:fem/src/LoadMod.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (48); no T014 relaxation landed |
| fortran:fem/src/Lua.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (19); no T014 relaxation landed |
| fortran:fem/src/LumpingUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (15); no T014 relaxation landed |
| fortran:fem/src/MGPrec.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (2); no T014 relaxation landed |
| fortran:fem/src/MagnetoDynamicsUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (14); no T014 relaxation landed |
| fortran:fem/src/MainUtils.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/MaterialModels.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (8); no T014 relaxation landed |
| fortran:fem/src/MatrixAssembly.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (30); no T014 relaxation landed |
| fortran:fem/src/MatrixScaling.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (10); no T014 relaxation landed |
| fortran:fem/src/MeshAllocations.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (8); no T014 relaxation landed |
| fortran:fem/src/MeshBasics.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (79); no T014 relaxation landed |
| fortran:fem/src/MeshExtrusion.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (8); no T014 relaxation landed |
| fortran:fem/src/MeshGenerate.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (5); no T014 relaxation landed |
| fortran:fem/src/MeshGraph.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (15); no T014 relaxation landed |
| fortran:fem/src/MeshIO.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (35); no T014 relaxation landed |
| fortran:fem/src/MeshLoad.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/MeshPartition.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (40); no T014 relaxation landed |
| fortran:fem/src/MeshRemeshing.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (20); no T014 relaxation landed |
| fortran:fem/src/MeshSplit.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (7); no T014 relaxation landed |
| fortran:fem/src/MeshTagging.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/MeshTransform.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (8); no T014 relaxation landed |
| fortran:fem/src/MeshUtils.F90 | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| fortran:fem/src/Messages.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (7); no T014 relaxation landed |
| fortran:fem/src/ModelDescription.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (43); no T014 relaxation landed |
| fortran:fem/src/ModelLumping.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (8); no T014 relaxation landed |
| fortran:fem/src/MortarUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (45); no T014 relaxation landed |
| fortran:fem/src/Multigrid.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/NavierStokes.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (7); no T014 relaxation landed |
| fortran:fem/src/NavierStokesCylindrical.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/src/NavierStokesGeneral.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/src/OptimizationUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (27); no T014 relaxation landed |
| fortran:fem/src/PElementBase.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (150); no T014 relaxation landed |
| fortran:fem/src/PElementMaps.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (35); no T014 relaxation landed |
| fortran:fem/src/PElementVisual.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/src/ParallelEigenSolve.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (12); no T014 relaxation landed |
| fortran:fem/src/ParallelUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (44); no T014 relaxation landed |
| fortran:fem/src/ParticleUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (84); no T014 relaxation landed |
| fortran:fem/src/PowellOptBobyqa.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (11); no T014 relaxation landed |
| fortran:fem/src/PowellOptNewuoa.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (10); no T014 relaxation landed |
| fortran:fem/src/ProjectorUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (14); no T014 relaxation landed |
| fortran:fem/src/Radiation.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/src/RadiationFactors.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (30); no T014 relaxation landed |
| fortran:fem/src/Radiators.F90 | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| fortran:fem/src/ResultToPost.F90 | portable-pending | gate-clean (all exports BIND(C)), in build; later batch |
| fortran:fem/src/ResultToResult.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/src/SParIterComm.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/SParIterGlobals.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/src/SParIterPrecond.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (7); no T014 relaxation landed |
| fortran:fem/src/SParIterSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (31); no T014 relaxation landed |
| fortran:fem/src/SaveUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (16); no T014 relaxation landed |
| fortran:fem/src/Smoothers.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (35); no T014 relaxation landed |
| fortran:fem/src/SolidMechanicsUtils.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/SolveBand.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (4); no T014 relaxation landed |
| fortran:fem/src/SolveCore.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (58); no T014 relaxation landed |
| fortran:fem/src/SolveSBand.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (4); no T014 relaxation landed |
| fortran:fem/src/Solver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (2); no T014 relaxation landed |
| fortran:fem/src/SolverActivate_x.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (2); no T014 relaxation landed |
| fortran:fem/src/SolverBasics.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/SolverUtils.F90 | ported | batch1: local differential byte-identical; machine grade vs backstop pending |
| fortran:fem/src/Stress.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (45); no T014 relaxation landed |
| fortran:fem/src/StressGeneral.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (3); no T014 relaxation landed |
| fortran:fem/src/TimeIntegrate.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (17); no T014 relaxation landed |
| fortran:fem/src/Types.F90 | portable-pending | gate-clean but LARGE: 108 __types_MOD___copy_* symbols (derived-type copies); later batch |
| fortran:fem/src/VankaCreate.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/VectorHelmholtzUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (4); no T014 relaxation landed |
| fortran:fem/src/ViewFactors.F90 | refused-COMMON | ADR-027: F77 COMMON-block file, refused until follow-up design (brief s2) |
| fortran:fem/src/ViewUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (9); no T014 relaxation landed |
| fortran:fem/src/Walls.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (8); no T014 relaxation landed |
| fortran:fem/src/ZirkaHysteresis.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (31); no T014 relaxation landed |
| fortran:fem/src/binio/biniomod.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (14); no T014 relaxation landed |
| fortran:fem/src/minpack.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (28); no T014 relaxation landed |
| fortran:fem/src/minpack_kinds.F90 | skipped-not-in-build | WITH_EXTOPTIM off: minpack*.F90 excluded from solverlib_SOURCES, never compiled (port kept in tree) |
| fortran:fem/src/mpif_stubs.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (40); no T014 relaxation landed |
| fortran:fem/src/solve_cmplx.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (2); no T014 relaxation landed |
| fortran:fem/src/solve_real.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports (4); no T014 relaxation landed |
| python:fem/src/split_ElementDescription.py | skipped-helper | build-time code generator, not linked |

counts: {'portable-pending': 14, 'ported': 18, 'skipped-ext-dep': 5, 'skipped-not-in-build': 9, 'skipped-helper': 7, 'refused-gate': 105, 'refused-COMMON': 10} total: 168
