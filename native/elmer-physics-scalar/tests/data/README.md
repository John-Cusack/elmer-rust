# Reference runs

Everything here was produced by the pristine ElmerSolver (`build/fem/src/ElmerSolver`,
Version 26.2-devel), run serially (`OMP_NUM_THREADS=1`) on copies of quick tests
in `fem/tests`, with `Output File = "case.result"` added to `Simulation`. The
ASCII `.result` holds 17 significant digits in node order; `solution.dat` /
`fields.dat` are `node v1 v2 ...` lines of its last time level. Meshes come
from `ElmerGrid 1 2 <grd>` as each test's `runtest.cmake` does, except where
noted. `case.sif` is the sif that was run.

Norms are Elmer's `sqrt(sum(x^2)/n)` over all values of the variable (both
components for the complex pressure).

## StatElecSolve

- `elstat_coarse/` is `fem/tests/elstat` on a coarser mesh: `Reference
  Density` 0.3 in `elstat_coarse.grd` raised to 1.0 (1323 nodes, 905
  hexahedra), so the mesh files stay small. The full-size run (32296 nodes)
  has norm 584.1590981586087 (the sif's Reference Norm is 584.159),
  `res: electric energy` 1.067158481342e-10. The sif run here also turns on
  `Calculate Electric Field`, `Calculate Electric Flux`, `Calculate Nodal
  Energy = Logical True` and `Calculate Nodal Flux = Logical True`. Coarse mesh:
  potential norm 631.9878246439413, `res: electric energy` 1.073050209726e-10,
  `res: potential difference` 1000, `res: effective capacitance`
  2.146100419453e-16. `fields.dat` columns: potential, electric energy density,
  nodal energy density, electric field 1-3, electric flux 1-3, nodal electric
  flux 1-3. The reference solve is BiCGStab/ILU1 to 1e-10, the tests solve
  directly, hence the 1e-8 tolerance. `Coordinate Scaling = 1.0e-6` is applied
  by the test when it reads the mesh.
- `elstat_source/` is `fem/tests/elstat_source` unchanged (2D quad8 mesh from
  `Step.grd`). Norm 795.5054410960914.
- `statelec_tet/` is a lane-made `case.sif` for p tetrahedra: `Element =
  p:2` on the StokesStabilized tet mesh (143 nodes, 386 tetrahedra; 803 dofs,
  no bubbles), `Relative Permittivity = 1 + x` (MATC in the sif, a linear
  table in the test), `Charge Density = 0.25`, potential 0 and 1 on
  boundaries 5 and 6. `tet_perm.dat` / `tet_solution.dat` are the pristine
  `case.result` slot by slot; norm 0.26499018387576206.

## StatCurrentSolve

- `current/` is `fem/tests/current` unchanged (`Step` mesh as shipped in the
  test). Conductivity `Variable Time` table, two steady-state iterations, norm
  7.6811928421903835.
- `current_control/` is `fem/tests/current_heat_control` (`Step.grd` meshed
  with ElmerGrid), including its heat solver; only the potential is compared.
  `Max Output Level = 7` shows `Total Heating Power` 1.3140158060369138e-4,
  `Effective Resistance` 7610.2585327037350 and `Control Scaling`
  761.02585327037355 for `Current Control = 0.1`. Potential norm after the
  scaling 378.7865721979089.
- `current_tensor/` is a derived case on the `current` mesh: a full symmetric
  `Electric Conductivity(2,2)`, a Robin condition (`External Conductivity` 3e-3,
  `External Potential` 2) on boundary 1, and all `Calculate ... Joule Heating /
  Volume Current / Electric Field / Nodal Heating / Nodal Current` requests.
  Total heating power 3.3061969337776302e-4. `solution.dat` columns: potential,
  joule heating, nodal joule heating, volume current 1-2, nodal current 1-2,
  electric field 1-2.
- `current_transient/` is a derived case on the same mesh: `Simulation Type =
  Transient`, BDF order 1, six steps of 5e-9 s, constant conductivity 1e-3,
  `Relative Permittivity = 4`, 1 V on boundary 1 and 0 V on boundary 2. The
  potential at step 6 (norm 0.2844367) is compared.

## HelmholtzSolve

- `helmholtz_fem/` is the first pass of `fem/tests/HelmholtzFEM` (`MeshFiles`
  as shipped), with a direct UMFPACK solve instead of BiCGStab/ILU2 so the
  solution is exact: pressure norm 0.012831819 (the sif's Reference Norm is
  0.128318191837E-01 for the iterative solve to 1e-8).
- `helmholtz_plane_waves/` is `fem/tests/HelmholtzPlaneWaves` unchanged apart
  from `Output File`. Its Reference Norm 0.99967001 is the `min: pressure` that
  SaveScalars reports, 0.9996700116143.
