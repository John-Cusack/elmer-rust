# femcore lane A unit ledger (claim A: 177 of the 353 refused/skipped fem core units)
# generated mechanically from the relaxed-gate re-grade (`data/elmer-lanes/femcore-regrade-a.json` in rustsmith);
# the reasons are the gate's CURRENT reasons, not the pre-relaxation ledger's.
# status: ported = Rust port in crates/elmer-fem-a/ports + byte-identical corpus;
#   skipped-* = recorded skip, never forced.
# The shared lane-fem-a ledger is left untouched on purpose: the sibling
# femcore lane owns the other half of the same fem core and the campaign
# forbids same-file edits across lanes (runbook §3 conflict rule).

| unit | status | reason |
|---|---|---|
| c:fem/obsolete/SC2Elmer.c | skipped-not-in-build | gate admits vacuously (no non-BIND(C) Fortran export); fem/obsolete has no CMakeLists.txt: never compiled |
| c:fem/src/SolveHypre.c | skipped-stub-only | gate admits vacuously (no non-BIND(C) Fortran export); compiled as an empty stub (HAVE_HYPRE off by default); the HYPRE solver body is not in the build |
| c:fem/src/SolveSuperLU.c | skipped-stub-only | gate admits vacuously (no non-BIND(C) Fortran export); compiled as an empty stub (HAVE_SUPERLU off by default); the SuperLU solver body is not in the build |
| c:fem/src/amgx.c | skipped-stub-only | gate admits vacuously (no non-BIND(C) Fortran export); compiled as an empty stub (HAVE_AMGX off by default); the AMGX solver body is not in the build |
| c:fem/src/cholmod.c | skipped-stub-only | gate admits vacuously (no non-BIND(C) Fortran export); compiled as an empty stub (WITH_CHOLMOD=FALSE by default); the CHOLMOD solver body is not in the build |
| c:fem/src/elmer_lua_iface.c | skipped-not-in-build | gate admits vacuously (no non-BIND(C) Fortran export); added to solverlib_SOURCES only under IF(WITH_LUA), and WITH_LUA defaults to FALSE: never compiled |
| c:fem/src/elmer_lua_iface.h | skipped-helper | gate admits vacuously (no non-BIND(C) Fortran export); header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/elmerf90-helper.in.h | skipped-helper | gate admits vacuously (no non-BIND(C) Fortran export); header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/fft_test.c | skipped-not-in-build | gate admits vacuously (no non-BIND(C) Fortran export); no CMake target compiles it (the test main is never built) |
| c:fem/src/mpif_stub.h | skipped-helper | gate admits vacuously (no non-BIND(C) Fortran export); header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/view3d/Ipoints.h | skipped-helper | gate admits vacuously (no non-BIND(C) Fortran export); header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/view3d/LinearSolve.c | skipped-not-in-build | gate admits vacuously (no non-BIND(C) Fortran export); orphan source: no CMakeLists references it (re-verified 2026-10-20) |
| c:fem/src/view3d/LinearSolveConjGrad.c | skipped-not-in-build | gate admits vacuously (no non-BIND(C) Fortran export); orphan source: no CMakeLists references it (re-verified 2026-10-20) |
| c:fem/src/view3d/LinearSolveGaussSeidel.c | skipped-not-in-build | gate admits vacuously (no non-BIND(C) Fortran export); orphan source: no CMakeLists references it (re-verified 2026-10-20) |
| c:fem/src/view3d/LinearSolveJacob.c | skipped-not-in-build | gate admits vacuously (no non-BIND(C) Fortran export); orphan source: no CMakeLists references it (re-verified 2026-10-20) |
| c:fem/src/view3d/ViewFactors.h | skipped-helper | gate admits vacuously (no non-BIND(C) Fortran export); header-only unit: links, never ports (matc T005 precedent) |
| c:fem/src/viewaxis/viewfact.h | skipped-helper | gate admits vacuously (no non-BIND(C) Fortran export); header-only unit: links, never ports (matc T005 precedent) |
| cxx:fem/src/rocalution.cpp | skipped-stub-only | gate admits vacuously (no non-BIND(C) Fortran export); compiled as an empty stub (WITH_ROCALUTION=FALSE by default); the Rocalution body is not in the build |
| fortran:fem/examples/turing/Turing.f90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/obsolete/MovingElstatSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/ADIOS2Utils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (16) |
| fortran:fem/src/Adaptive.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (16); procedure-not-found-in-gate-parse (6); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/AddrFunc.F90 | skipped-gate-refused | gate refuses: procedure dummy (code address) (2) |
| fortran:fem/src/BackwardError.F90 | skipped-gate-admitted-not-portable | gate admits the dummy shape, but 3 of the 5 procedures read the global matrix (Types.GlobalMatrix) and call the Elmer CRS matrix-vector products; the file-level port must define all five symbols, so a Rust port would re-implement the CRS core, not this unit |
| fortran:fem/src/BandMatrix.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (9); assumed-size dummy (later shape) (2) |
| fortran:fem/src/BandwidthOptimize.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/BlockSolve.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (27); assumed-size dummy (later shape) (3); optional dummy (hidden presence arg) (2) |
| fortran:fem/src/BoundaryConditionUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (17); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/CRSMatrix.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (58); explicit-shape dummy (later shape) (5); non-REAL(8) assumed-shape dtype (3); assumed-size dummy (later shape) (3); procedure-not-found-in-gate-parse (2); optional dummy (hidden presence arg) (1) |
| fortran:fem/src/CircuitUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (26); character dummy (hidden length param) (7); pointer dummy (libgfortran realloc contract) (6); non-REAL(8) assumed-shape dtype (5); optional dummy (hidden presence arg) (2); undeclared function result (implicit typing) (2) |
| fortran:fem/src/ClusteringMethods.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (6); optional dummy (hidden presence arg) (2) |
| fortran:fem/src/Constitutive.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (7); undeclared function result (implicit typing) (3); procedure-not-found-in-gate-parse (3) |
| fortran:fem/src/ContactUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); undeclared function result (implicit typing) (1) |
| fortran:fem/src/CoordinateSystems.F90 | skipped-gate-refused | gate refuses: undeclared function result (implicit typing) (2); explicit-shape dummy (later shape) (1) |
| fortran:fem/src/CutFEMUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (15); explicit-shape dummy (later shape) (3); non-REAL(8) assumed-shape dtype (2); undeclared function result (implicit typing) (1) |
| fortran:fem/src/DefUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (101); character dummy (hidden length param) (13); non-REAL(8) assumed-shape dtype (10); optional dummy (hidden presence arg) (5); assumed-shape rank > 2 (1); undeclared function result (implicit typing) (1); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/Differentials.F90 | skipped-gate-refused | gate refuses: non-REAL(8) assumed-shape dtype (2); assumed-shape rank > 2 (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/DiffuseConvectiveAnisotropic.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (1); assumed-shape rank > 2 (1) |
| fortran:fem/src/DiffuseConvectiveGeneralAnisotropic.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (1); assumed-shape rank > 2 (1) |
| fortran:fem/src/DirectSolve.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (28); procedure-not-found-in-gate-parse (4); explicit-shape dummy (later shape) (4) |
| fortran:fem/src/EigenSolve.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (13) |
| fortran:fem/src/ElemInfo.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (26); explicit-shape dummy (later shape) (8); optional dummy (hidden presence arg) (1); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/ElementBasis.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (24); explicit-shape dummy (later shape) (1); pointer dummy (libgfortran realloc contract) (1) |
| fortran:fem/src/ElementGeometry.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (14); explicit-shape dummy (later shape) (2) |
| fortran:fem/src/ElementUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (30); explicit-shape dummy (later shape) (2); non-REAL(8) assumed-shape dtype (2); undeclared dummy (implicit typing) (1) |
| fortran:fem/src/ElmerSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (22); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/ExchangeCorrelations.F90 | skipped-gate-refused | gate refuses: undeclared dummy (implicit typing) (6) |
| fortran:fem/src/Feti.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (18); explicit-shape dummy (later shape) (2); optional dummy (hidden presence arg) (2); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/FreeSurface.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/GeneralUtils.F90 | skipped-gate-refused | gate refuses: character dummy (hidden length param) (16); optional dummy (hidden presence arg) (10); pointer dummy (libgfortran realloc contract) (9); non-REAL(8) assumed-shape dtype (8); derived-type dummy (needs whole-program layout) (7); undeclared function result (implicit typing) (6); assumed-size dummy (later shape) (1) |
| fortran:fem/src/GeometryFitting.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (5) |
| fortran:fem/src/H1Basis.F90 | skipped-gate-refused | gate refuses: explicit-shape dummy (later shape) (75); non-REAL(8) assumed-shape dtype (2); undeclared function result (implicit typing) (1) |
| fortran:fem/src/HashTable.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (10); character dummy (hidden length param) (2) |
| fortran:fem/src/Integration.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (16); optional dummy (hidden presence arg) (7); explicit-shape dummy (later shape) (5); undeclared function result (implicit typing) (2) |
| fortran:fem/src/InterpVarToVar.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/InterpolateMeshToMesh.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (7); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/Interpolation.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (7); explicit-shape dummy (later shape) (6) |
| fortran:fem/src/IpFieldInterface.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/IterSolve.F90 | skipped-gate-refused | gate refuses: assumed-size dummy (later shape) (7); procedure-not-found-in-gate-parse (6); derived-type dummy (needs whole-program layout) (2); explicit-shape dummy (later shape) (2); character dummy (hidden length param) (1); undeclared function result (implicit typing) (1) |
| fortran:fem/src/IterativeMethods.F90 | skipped-gate-refused | gate refuses: explicit-shape dummy (later shape) (12); derived-type dummy (needs whole-program layout) (9); assumed-size dummy (later shape) (4) |
| fortran:fem/src/LinearAlgebra.F90 | skipped-gate-refused | gate refuses: explicit-shape dummy (later shape) (7); non-REAL(8) assumed-shape dtype (7); procedure-not-found-in-gate-parse (2); optional dummy (hidden presence arg) (1) |
| fortran:fem/src/LinearForms.F90 | skipped-gate-refused | gate refuses: optional dummy (hidden presence arg) (3); assumed-shape rank > 2 (2) |
| fortran:fem/src/ListMatrix.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (24); non-REAL(8) assumed-shape dtype (2) |
| fortran:fem/src/ListMatrixArray.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (17) |
| fortran:fem/src/Lists.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (147); character dummy (hidden length param) (19); procedure-not-found-in-gate-parse (2); non-REAL(8) assumed-shape dtype (1); undeclared function result (implicit typing) (1) |
| fortran:fem/src/LoadMod.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (15); procedure-not-found-in-gate-parse (14); procedure dummy (code address) (7); character dummy (hidden length param) (6); non-REAL(8) assumed-shape dtype (2); explicit-shape dummy (later shape) (2) |
| fortran:fem/src/Lua.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (16); character dummy (hidden length param) (1); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/LumpingUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (8); explicit-shape dummy (later shape) (1); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/MGPrec.F90 | skipped-gate-admitted-not-portable | gate admits the dummy shape, but the procedure body is Elmer's multigrid preconditioner over the global matrix (Types.CurrentModel/GlobalMatrix + Multigrid.MultiGridSolve + IterSolve.FirstCall): a Rust port would re-implement the multigrid core, not this unit |
| fortran:fem/src/MagnetoDynamicsUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (9); explicit-shape dummy (later shape) (4) |
| fortran:fem/src/MainUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (29); procedure-not-found-in-gate-parse (5); explicit-shape dummy (later shape) (1); optional dummy (hidden presence arg) (1) |
| fortran:fem/src/MaterialModels.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); character dummy (hidden length param) (2); explicit-shape dummy (later shape) (1) |
| fortran:fem/src/MatrixAssembly.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (23); optional dummy (hidden presence arg) (3); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/MatrixScaling.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (5) |
| fortran:fem/src/MeshAllocations.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (6); optional dummy (hidden presence arg) (1) |
| fortran:fem/src/MeshBasics.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (53); allocatable dummy (libgfortran realloc contract) (3); explicit-shape dummy (later shape) (3); undeclared function result (implicit typing) (2); character dummy (hidden length param) (1) |
| fortran:fem/src/MeshExtrusion.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (7) |
| fortran:fem/src/MeshGenerate.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/MeshGraph.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (8); non-REAL(8) assumed-shape dtype (5); allocatable dummy (libgfortran realloc contract) (1) |
| fortran:fem/src/MeshIO.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (11); character dummy (hidden length param) (2); undeclared function result (implicit typing) (1); optional dummy (hidden presence arg) (1) |
| fortran:fem/src/MeshLoad.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/MeshPartition.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (19); assumed-size dummy (later shape) (3); allocatable dummy (libgfortran realloc contract) (2); non-REAL(8) assumed-shape dtype (2); undeclared function result (implicit typing) (1) |
| fortran:fem/src/MeshRemeshing.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (18); character dummy (hidden length param) (1) |
| fortran:fem/src/MeshSplit.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (5) |
| fortran:fem/src/MeshTagging.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); undeclared function result (implicit typing) (1) |
| fortran:fem/src/MeshTransform.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (6); character dummy (hidden length param) (1) |
| fortran:fem/src/Messages.F90 | skipped-gate-refused | gate refuses: character dummy (hidden length param) (5); optional dummy (hidden presence arg) (1) |
| fortran:fem/src/ModelDescription.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (15); character dummy (hidden length param) (14); allocatable dummy (libgfortran realloc contract) (2); optional dummy (hidden presence arg) (1); non-REAL(8) assumed-shape dtype (1); pointer dummy (libgfortran realloc contract) (1) |
| fortran:fem/src/ModelLumping.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (7) |
| fortran:fem/src/MortarUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (34); undeclared function result (implicit typing) (1); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/Multigrid.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (28); procedure-not-found-in-gate-parse (1); pointer dummy (libgfortran realloc contract) (1); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/NavierStokes.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); character dummy (hidden length param) (1); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/NavierStokesCylindrical.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (1); character dummy (hidden length param) (1) |
| fortran:fem/src/NavierStokesGeneral.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/OptimizationUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (8); undeclared function result (implicit typing) (5); procedure dummy (code address) (3); explicit-shape dummy (later shape) (2) |
| fortran:fem/src/PElementBase.F90 | skipped-gate-refused | gate refuses: optional dummy (hidden presence arg) (60); undeclared function result (implicit typing) (4); assumed-shape rank > 2 (1) |
| fortran:fem/src/PElementMaps.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (19); optional dummy (hidden presence arg) (2) |
| fortran:fem/src/PElementVisual.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (1); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/ParallelEigenSolve.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (8); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/ParallelUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (24); non-REAL(8) assumed-shape dtype (3); undeclared function result (implicit typing) (3); optional dummy (hidden presence arg) (3); allocatable dummy (libgfortran realloc contract) (1) |
| fortran:fem/src/ParticleUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (69); character dummy (hidden length param) (7); procedure-not-found-in-gate-parse (2); explicit-shape dummy (later shape) (2) |
| fortran:fem/src/PowellOptBobyqa.F90 | skipped-gate-refused | gate refuses: undeclared dummy (implicit typing) (6); non-8-byte REAL kind (1); procedure-not-found-in-gate-parse (1); undeclared function result (implicit typing) (1) |
| fortran:fem/src/PowellOptNewuoa.F90 | skipped-gate-refused | gate refuses: undeclared dummy (implicit typing) (5); procedure-not-found-in-gate-parse (1); assumed-size dummy (later shape) (1); undeclared function result (implicit typing) (1) |
| fortran:fem/src/ProjectorUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (9) |
| fortran:fem/src/Radiation.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/RadiationFactors.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (28); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/ResultToResult.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); external scope owns no procedure (1) |
| fortran:fem/src/SParIterComm.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (25); assumed-size dummy (later shape) (8); undeclared function result (implicit typing) (2); procedure-not-found-in-gate-parse (2); pointer dummy (libgfortran realloc contract) (1); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/SParIterGlobals.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/SParIterPrecond.F90 | skipped-gate-refused | gate refuses: assumed-size dummy (later shape) (4); explicit-shape dummy (later shape) (1); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/SParIterSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (20); assumed-size dummy (later shape) (7) |
| fortran:fem/src/SaveUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (7); undeclared function result (implicit typing) (1); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/Smoothers.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (29); non-REAL(8) assumed-shape dtype (2); undeclared dummy (implicit typing) (1) |
| fortran:fem/src/SolidMechanicsUtils.F90 | skipped-gate-refused | gate refuses: optional dummy (hidden presence arg) (2); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/SolveBand.F90 | ported | gate admits (T017 shape A: F77 global external subprogram); port crate solveband |
| fortran:fem/src/SolveCore.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (39); procedure-not-found-in-gate-parse (4); pointer dummy (libgfortran realloc contract) (2); undeclared function result (implicit typing) (1); optional dummy (hidden presence arg) (1) |
| fortran:fem/src/SolveSBand.F90 | ported | gate admits (T017 shape A: F77 global external subprogram); port crate solvesband |
| fortran:fem/src/Solver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (1); external scope owns no procedure (1) |
| fortran:fem/src/SolverActivate_x.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/SolverBasics.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (58); character dummy (hidden length param) (6); non-REAL(8) assumed-shape dtype (5); undeclared function result (implicit typing) (4); pointer dummy (libgfortran realloc contract) (1); explicit-shape dummy (later shape) (1) |
| fortran:fem/src/Stress.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (10); assumed-shape rank > 2 (6); undeclared function result (implicit typing) (5); explicit-shape dummy (later shape) (4); non-REAL(8) assumed-shape dtype (1); procedure-not-found-in-gate-parse (1); character dummy (hidden length param) (1) |
| fortran:fem/src/StressGeneral.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (1); assumed-shape rank > 2 (1) |
| fortran:fem/src/TimeIntegrate.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (7); undeclared dummy (implicit typing) (1) |
| fortran:fem/src/VankaCreate.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (10); undeclared dummy (implicit typing) (2); derived-type dummy (needs whole-program layout) (1) |
| fortran:fem/src/VectorHelmholtzUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); explicit-shape dummy (later shape) (1) |
| fortran:fem/src/ViewFactors.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (20); external scope owns no procedure (1) |
| fortran:fem/src/ViewUtils.F90 | skipped-gate-refused | gate refuses: explicit-shape dummy (later shape) (3); derived-type dummy (needs whole-program layout) (2); character dummy (hidden length param) (1); non-REAL(8) assumed-shape dtype (1) |
| fortran:fem/src/Walls.F90 | ported | gate admits (T017 shape A: F77 global external subprogram); port crate walls |
| fortran:fem/src/ZirkaHysteresis.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (29) |
| fortran:fem/src/binio/biniomod.f90 | skipped-gate-refused | gate refuses: optional dummy (hidden presence arg) (8); character dummy (hidden length param) (5) |
| fortran:fem/src/minpack.F90 | skipped-gate-refused | gate refuses: procedure dummy (code address) (12); explicit-shape dummy (later shape) (10); procedure-not-found-in-gate-parse (5) |
| fortran:fem/src/minpack_kinds.F90 | skipped-not-in-build | gate admits vacuously (no non-BIND(C) Fortran export); added to solverlib_SOURCES only under IF(WITH_EXTOPTIM), which is off by default: never compiled (the lane-fem-a port is kept in the tree) |
| fortran:fem/src/modules/ADIOS2OutputSolver.F90#adios2outputsolver_init | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (3) |
| fortran:fem/src/modules/ADIOS2OutputSolver.F90#adiosoutputsolverutils | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); character dummy (hidden length param) (1) |
| fortran:fem/src/modules/Acoustics.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (47); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/AdvectionDiffusion.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (6); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/AdvectionReaction.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/AllocateSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/AmsPrecSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (7); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/ArteryOutlet.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (7); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/ArtificialCompressibility.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/BatterySolver/BatteryPost.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/BatterySolver/BatteryUtils.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (11); undeclared function result (implicit typing) (2) |
| fortran:fem/src/modules/BatterySolver/ElectrolyteCons.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/BatterySolver/ElectrolytePot.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/BatterySolver/SolidPhaseCons.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/BatterySolver/SolidPhasePot.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/BeamSolver3D.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/CahnHilliard.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/CircuitsAndDynamics.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (13); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/CoilSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (17); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/CompressibleNS.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/ConsistentSplitting/CompressibilitySolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/ConsistentSplitting/OptimalSolutionUpdate.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (9); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ConsistentSplitting/PressureSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/ConsistentSplitting/VelocitySolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/CoordinateTransform.F90#rotmsolver_init | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/CoordinateTransform.F90#vecutils | skipped-gate-refused | gate refuses: explicit-shape dummy (later shape) (4); the sibling T017 lane's port (elmer-rust#11) covers this scope under the shape B relaxation (rustsmith#108, open at this lane's done-gate): the gate still refuses it under the merged gate |
| fortran:fem/src/modules/CraigBamptonSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/CyclicConvergence.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/DCRComplexSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (7); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/DFTSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/DataToFieldSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/DirectionSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (4); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/DistanceSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (5); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/DistributeSource.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/DivProjection.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/DivergenceSolver.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (2); derived-type dummy (needs whole-program layout) (2) |
| fortran:fem/src/modules/EMPort.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (8); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/EMWaveSolver.F90#emwavesolver_init0 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (12); derived-type dummy (needs whole-program layout) (3) |
| fortran:fem/src/modules/EMWaveSolver.F90#emwavesolverutils | skipped-vacuous-admit | the scope owns no portable procedure (an empty module): the gate admits it vacuously (ADR-047 precedent, skipped-vacuous-admit) |
| fortran:fem/src/modules/EdgeFEInterpolationTest.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (3); procedure-not-found-in-gate-parse (2) |
| fortran:fem/src/modules/ElasticSolve.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (21); derived-type dummy (needs whole-program layout) (4) |
| fortran:fem/src/modules/ElectricForce.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/Electrokinetics.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (7) |
| fortran:fem/src/modules/ElementSizeSolver.F90 | skipped-gate-refused | gate refuses: derived-type dummy (needs whole-program layout) (2); procedure-not-found-in-gate-parse (1) |
| fortran:fem/src/modules/ElementStats.F90 | skipped-gate-refused | gate refuses: procedure-not-found-in-gate-parse (3); derived-type dummy (needs whole-program layout) (2) |

counts: ported=3, skipped-gate-admitted-not-portable=2, skipped-gate-refused=152, skipped-helper=6, skipped-not-in-build=8, skipped-stub-only=5, skipped-vacuous-admit=1 total=177
