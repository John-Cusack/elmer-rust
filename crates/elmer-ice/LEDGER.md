# lane-ice unit ledger (184 units)
# generated mechanically from lane_ice_records.json (T008 triage re-land); reasons reviewed by hand.
# status: graded = Rust port in crates/elmer-ice + differential corpus green;
# refused-gate = recorded skip (non-BIND(C) Fortran or non-linkable python), never forced;
# skipped-* = recorded skip (doc stub / interface fragment), never forced.

| unit | status | reason |
|---|---|---|
| c:elmerice/Meshers/ExtrudeMesh.c | graded | crates/elmer-ice/extrude: differential corpus byte-identical (see README) |
| c:elmerice/Solvers/ScatteredDataInterpolator/csa_interpolate_points.c | graded | crates/elmer-ice/csainterp: differential corpus byte-identical (see README) |
| fortran:elmerice/IceSheet/Greenland/MAKE_MESH/USF_PCond.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/IceSheet/Greenland/SSA/Scalar_OUTPUT.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/IceSheet/Greenland/SSA_FRICTION_INIT/USFs.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/IceSheet/Tools/ConservativeInterpolation/CELL_TO_NODE.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/IceSheet/Tools/ConservativeInterpolation/MakeCDOGridDesc.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/IceSheet/Tools/ConservativeInterpolation/READ_CELL_NETCDF.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Meshers/MshGlacier.F90 | graded | crates/elmer-ice/mshglacier: differential corpus byte-identical (see README) |
| fortran:elmerice/Meshers/MshGlacierDEM.F90 | graded | crates/elmer-ice/mshglacierdem: differential corpus byte-identical (see README) |
| fortran:elmerice/Meshers/MshGlacierSynthetic.F90 | graded | crates/elmer-ice/mshglaciersynth: differential corpus byte-identical (see README) |
| fortran:elmerice/Solvers/AIFlowSolve_nlD2.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AIFlowSolve_nlS2.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Adjoint/Adjoint_CostContSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Adjoint/Adjoint_CostDiscSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Adjoint/Adjoint_CostRegSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Adjoint/Adjoint_GradientValidation.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Adjoint/Adjoint_LinearSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointSSA/AdjointSSA_AdjointSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointSSA/AdjointSSA_CostContSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointSSA/AdjointSSA_CostDiscSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointSSA/AdjointSSA_CostFluxDivSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointSSA/AdjointSSA_CostRegSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointSSA/AdjointSSA_CostTaubSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointSSA/AdjointSSA_GradientSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointSSA/AdjointSSA_SSASolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointStokes/AdjointStokes_GradientBetaSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointStokes/AdjointStokes_GradientMu.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointThickness/AdjointThickness_GradientSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/AdjointThickness/AdjointThickness_ThicknessSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/BasalMelt3D.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CaffeSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Calving.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Calving3D.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Calving3D_lset.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CalvingFrontAdvance3D.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CalvingGeometry.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CalvingGlacierAdvance3D.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CalvingHydroInterp.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CalvingRemesh.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CalvingRemeshMMG.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CalvingRemeshparMMG.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ComputeCalvingNormal.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ComputeDevStress.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ComputeEigenValues.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ComputeGroundingLineFlux.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ComputeNormal.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ComputeStrainRate.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CostSolver_Adjoint.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/CostSolver_Robin.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Covarianceutils/BackgroundErrorCostSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Covarianceutils/CovarianceUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Covarianceutils/CovarianceVectorMultiplySolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Covarianceutils/GaussianSimulationSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/DJDBeta_Adjoint.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/DJDBeta_Robin.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/DJDmu_Adjoint.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/DJDmu_Robin.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/DeformationalHeat.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Documentation/CalvingHydroInterp.F90 | skipped-doc | doc stub: no exports, nothing to port |
| fortran:elmerice/Solvers/EPLSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ElmerIceUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Emergence.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/EnthalpySolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ExportVertically.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Flotation.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Flowdepth.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ForceToStress.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/FrontDisplacement.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/GMValid.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/GetHydrostaticLoads.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/GlaDSCoupledSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/GlaDSchannelSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/GolfLaw.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Grid2DInterpolator.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/GridDataReader/GridDataReader.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/GroundedSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/HydroRestart.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/IDSSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/IcyMaskSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/IntegrateVertically.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/IntegratedVelocity.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/MMG3DSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/MeshAdaptation_2D/Compute2DNodalGradient.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/MeshAdaptation_2D/GetActiveMesh.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/MeshAdaptation_2D/GetConnectedAreas.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/MeshAdaptation_2D/MMG2DSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/MeshAdaptation_2D/MMG2D_MetricAniso.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/MeshAdaptation_2D/MMG2D_MetricIntersect.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Optimize_m1qn3Parallel.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/OutPutSolvers/XIOSOutputSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/OutputStrainHeating.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Permafrost/PermafrostMaterials.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Permafrost/Permafrost_Darcy.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Permafrost/Permafrost_HTEQ.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Permafrost/Permafrost_Utils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Permafrost/Permafrost_solid.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Permafrost/Permafrost_solute.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/PlumeSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/PorousSolve.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ProjectCalving.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/SIASolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/SSASolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/SSAmask.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Scalar_OUTPUT_Glacier.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ScatteredDataInterpolator/Scattered2DDataInterpolator.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ScatteredDataInterpolator/Scattered2D_FInterface.F90#nearestneighbour | skipped-fragment | Scattered2D_FInterface fragment in COMMON file: no independent port surface |
| fortran:elmerice/Solvers/ScatteredDataInterpolator/Scattered2D_FInterface.F90#point_type | skipped-fragment | Scattered2D_FInterface fragment in COMMON file: no independent port surface |
| fortran:elmerice/Solvers/SubShelfMelt.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/SurfaceBoundaryEnthalpy.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/TemperateIce.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/ThicknessSolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/TwoMeshes.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/UGridDataReader.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/UpdateExport.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/Weertman2Coulomb.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/m1qn3.F | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Solvers/pointwise.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Calving2D/PROG/bedrockfunction_flat.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Calving3D/PROG/DummySolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Calving3D/PROG/bedrockfunction_3D.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Calving3D_lset/PROG/DummySolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Calving3D_lset/PROG/bedrockfunction_3D.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Calving3D_lset_parMMG/PROG/DummySolver.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Calving3D_lset_parMMG/PROG/bedrockfunction_3D.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/CovarianceVector/Init.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/CovarianceVector2/Init.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/DGsolver/PROG/DGtoNodalVariable.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/DGsolver/PROG/InitializeDGVariable.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Emergence/buelerprofile.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/FrictionHeat/DummySolver.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/FrictionHeatMasked/DummySolver.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Friction_Weertman2/MyModule.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Grounded/PROG/bedrock.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Grounded/PROG/fbed.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/GroundingLineFlux/gmInit.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/InvMeth_AdjRobin/PROG/GradientValidation.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/InvMeth_AdjRobin/PROG/USF_Init.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/MMG2D_Transient/PROG/USFs.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Proj6/PROJ_VALIDATION.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Teterousse3a/PROG/USF_TR.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Teterousse3a_b/PROG/USF_TR.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Teterousse_DeformHeat/PROG/USF_TR.f90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/UGridDataReader/src/Check.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Xios/Check.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Tests/Xios2/Check_xios2.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/Buoyancy.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/CaffeFlow.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_CoV.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_Contact.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_CouplingGlaDS_SSA.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_Damage.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_GetFrictionHeating.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_GlacierMeshMetric.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_Haf.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_IceProperties.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_LateralFriction.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_ShapeFactor.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_Sliding.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_SourceCalcCalving.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_WaterTransfer.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_Zs.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/UserFunctions/USF_proj.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Utils/ComputeFluxUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Utils/PorousMaterialModels.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Utils/ProjUtils.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Utils/SSAMaterialModels.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/Utils/proj6_interface.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/examples/InverseMethods_OLD/src/GradientValidation.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/examples/InverseMethods_OLD/src/USF_Init.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/examples/Inverse_Methods/MassConservation/src/RAMP.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/examples/Inverse_Methods/MassConservation/src/Random.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/examples/Inverse_Methods/src/MacAyeal_USFs.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/examples/Inverse_Methods/src/USFs_RonneFilchner.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/examples/Test_MshGlacierSynthetic/fbed.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| fortran:elmerice/examples/Test_MshGlacierSynthetic/fsurf.F90 | refused-gate | check_substitutable: non-BIND(C) Fortran exports; T014 relaxation covers assumed-shape REAL(8) file-granularity only |
| python:elmerice/Meshers/GIS/Contour2geo.py | refused-gate | non-BIND(C) gate: python unit, no C-ABI link surface |
| python:elmerice/Meshers/GIS/MeshToShp.py | refused-gate | non-BIND(C) gate: python unit, no C-ABI link surface |
| python:elmerice/Meshers/makemoulin.py | refused-gate | non-BIND(C) gate: python unit, no C-ABI link surface |
| python:elmerice/examples/Inverse_Methods/DATA/AddNoise.py | refused-gate | non-BIND(C) gate: python unit, no C-ABI link surface |
| python:elmerice/examples/Inverse_Methods/MassConservation/src/MakeObs.py | refused-gate | non-BIND(C) gate: python unit, no C-ABI link surface |
| python:elmerice/examples/Inverse_Methods/SCRIPTS/MakeReport.py | refused-gate | non-BIND(C) gate: python unit, no C-ABI link surface |
| python:elmerice/examples/Inverse_Methods/SCRIPTS/PlotHist.py | refused-gate | non-BIND(C) gate: python unit, no C-ABI link surface |

counts: {'graded': 5, 'refused-gate': 176, 'skipped-doc': 1, 'skipped-fragment': 2} total: 184
