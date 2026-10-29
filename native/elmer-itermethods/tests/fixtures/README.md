# Fortran reference data

`linsys_a.dat` / `linsys_b.dat` are the matrix (`i j value`, 1-based, CRS order,
explicit zeros kept) and right-hand side written by `Linear System Save = True`
from the pristine ElmerSolver (`orig_src/build/fem/src/ElmerSolver`, serial) on
the `fem/tests/linearsolvers` mesh with a `Heat Equation`, `Temperature = 1` on all
boundaries, `Linear System Preconditioning = None`, `Linear System Scaling =
False`, `Linear System Convergence Tolerance = 1e-8` (408 unknowns).

`ref_<case>.txt` hold the per-iteration lines Elmer printed for that matrix, as
`iteration relative-residual` (4 significant digits):

| file | method | printed by |
|---|---|---|
| `ref_gcr.txt`, `ref_gcr10.txt` | GCR, restart 100 / 10 | `gcr:` lines, column 2 |
| `ref_bicgstabl.txt`, `ref_bicgstabl4.txt` | BiCGStab(l), l = 2 / 4 | per-round lines (the converging round exits early and prints only the summary `BiCGStabl:  <rounds>`: 20 and 9) |
| `ref_idrs.txt`, `ref_idrs1.txt` | IDR(s), s = 4 / 1 | per-iteration lines (random shadow space, so only the count is compared) |
| `ref_jacobi.txt`, `ref_sgs.txt` | Jacobi, SGS (omega 1.8) | per-iteration lines, column 3 |

Final iteration counts: GCR(100) 53, GCR(10) 88, BiCGStab(2) 20, BiCGStab(4) 9,
IDR(4) 59, IDR(1) 66, Jacobi 687, SGS 130.

`linsys_x0.dat` is the initial temperature (the Dirichlet values already applied),
saved with `Linear System Save Solution`; the solvers start from it.

Complex references come from `fem/tests/linearsolvers_cmplx/HelmholtzFEM.sif` with
the same settings (iterative method, no preconditioner, no scaling, tolerance
1e-8, zero initial guess): `helm_a.dat` (`i j re im`, 1-based complex CSR rebuilt
from the dumped real `[[a, -b], [b, a]]` matrix, 2110 unknowns), `helm_b.dat`
(`i re im`), and `ref_zgcr.txt` (restart 100, 393 iterations), `ref_zgcr10.txt`
(restart 10, 466), `ref_zbic.txt` (BiCGStab(2), 119 rounds), `ref_zidrs.txt`
(IDR(4), 500 iterations).
