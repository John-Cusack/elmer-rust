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
