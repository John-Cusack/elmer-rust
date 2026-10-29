Captured from the pristine ElmerSolver (orig_src build), serial, `nice -n 10`.

Case: Elmer sif: heat equation on a 6x6-division triangulated unit square
(49 nodes), conductivity `1+2x`, source `1+x`, Dirichlet T=0 on x=0 and x=1,
`Optimize Bandwidth = False`, `Linear System Save = True`, solver
`CG` with `Linear System Max Iterations = 1` and
`Linear System Convergence Tolerance = 1e-30`, started from T=0.

* `sq_{a,b,perm,sizes}.dat`: the saved system (identical for every run below).
* `x1_<precond>.txt`: the Temperature after exactly one CG iteration from the
  VTU (`ES15.7`, 8 digits) for `Linear System Preconditioning` =
  ILU0, ILU1, ILU2, ILUT (`Linear System ILUT Tolerance = 0.05`), Diagonal.
  With x0 = 0 that is `alpha * M^-1 b`, which pins the preconditioner.
