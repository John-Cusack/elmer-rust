# Reference cases

`mgdyn_airgap/` is a verbatim copy of `fem/tests/mgdyn_airgap` (`airgap.sif`
and the `airgap/` mesh database) from the Elmer source tree.

The pristine ElmerSolver (`build/fem/src/ElmerSolver`, run serially on a
copy of the test) prints for Solver 1 (`WhitneyAVSolver`):

    CompareToReferenceSolution: Solver 1 PASSED:  Norm = 7.08538973E-05  RefNorm = 7.08538965E-05

The norm is `ComputeNorm` with `Norm Permutation = nNodes+1`: the RMS of the
edge DOFs only. The steady curl-curl system is not gauged, so the edge part
of the solution depends on the Krylov path. This case uses BiCGStab(l)
without a preconditioner, whose path does not depend on the DOF numbering.

`mgdyn_airgap_harmonic/` is a verbatim copy of `fem/tests/mgdyn_airgap_harmonic`
(`WhitneyAVHarmonicSolver`, complex BiCGStab(l) to 1e-6). The pristine run
prints `Norm = 2.42423001E-03` (reference 2.42423000E-03). At that
tolerance the norm hangs on rounding: BiCGStab(l)'s residual history on
this system changes from iteration 15 on under a 1e-15 relative
perturbation of the matrix, so a port that does not reproduce Elmer's
assembly bit for bit lands within about 1e-5 of the reference (this one at
1.4e-5). The test therefore compares the converged solution: with
`Linear System Convergence Tolerance = 1.0e-10` (the only change to the sif)
the pristine run takes 45 iterations and prints
`Norm = 2.42424125E-03`, which this crate matches.

`mgdyn_airgap` Solver 2 (`MagnetoDynamicsCalcFields`, `Calculate Magnetic
Field Strength`): the pristine run prints `Norm = 3.06334699E+03`
(reference 3.06334400E+03), the norm of the last nodal component solved,
`Magnetic Field Strength 3`. The projections are SPD mass-matrix solves
(CG to 1e-8), so this value does not hang on the solver path.

`mgdyn_airgap_harmonic/thinsheet.sif` is `fem/tests/mgdyn_thinsheet_harmonic/thinsheet.sif`
verbatim; that test uses the same `airgap` mesh, so it sits here instead of
copying the mesh again. The pristine run prints `Norm = 1.22780296E-03`
(reference 1.22780252E-03); with `Linear System Convergence Tolerance =
1.0e-10` it prints `Norm = 1.22780291E-03`.

`InductionHeating3/` holds `fem/tests/InductionHeating3/crucible.sif`
verbatim and the mesh `ElmerGrid 1 2 crucible` makes from its
`crucible.grd` (as the test's `runtest.cmake` does). The pristine run
prints `Norm = 8.88220240E-06` for Solver 1 (`MagnetoDynamics2DHarmonic`,
BiCGStab(l) with ILU to 1e-10 on a nonsingular system).

`InductionHeating4/` holds `fem/tests/InductionHeating4/crucible.sif`
verbatim and the mesh `ElmerGrid 1 2 crucible-bc` makes; the pristine run
prints `Norm = 8.85766857E-06` for Solver 1 (surface-current coil).
