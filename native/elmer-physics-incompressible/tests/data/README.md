# Pristine ElmerSolver references

Each directory holds the `case.sif` run by the pristine ElmerSolver
(`build/fem/src/ElmerSolver`), serially, on a mesh written by `ElmerGrid` and
the resulting reference:

- `mesh.*`: the ElmerGrid output the run used.
- `solution.dat`: the last time level of the ASCII `.result` that `Output File =
  "case.result"` writes (17 significant digits), as `node u v [w] p` lines in
  node order.
- The norm in `tests/fortran_reference.rs` is the last `ComputeChange` norm
  of the run's stdout (velocity components only: `Nonlinear System Norm DOFs =
  dim`).

| case | source | what it exercises |
|---|---|---|
| `stokes_tri` | `fem/tests/StokesStabilized/tri.sif` (`ElmerGrid 1 2 Step.grd -triangles`) | Stokes, 303, Brezzi-Pitkaranta block |
| `stokes_quad` | `StokesStabilized/quad.sif` (`ElmerGrid 1 2 Step.grd`) | Stokes, 404 |
| `stokes_tet` | `StokesStabilized/tet.sif` (`ElmerGrid 14 2 tet.msh`) | Stokes, 504, 3D |
| `ns_tri` | own | Navier-Stokes, Picard then Newton after 3, `Pressure Stabilization Element Size = longest edge` |
| `ns_quad` | own | Newton from iteration 1, `GradP` and `Div-Curl` forms, `volume` size, coefficient 1.5, `Flow Bodyforce 2`, `Normal Surface Traction`, `Surface Traction 2`, `Slip Coefficient 1` |
| `transient_tri` | own | BDF2, 4 steps of 0.1, Newton after 2 iterations: the PSPG/SUPG mass terms |
| `taylorhood_quad` | `StokesStabilized/quad_taylorhood.sif` (`ElmerGrid 1 2 Step.grd`) | Stokes on `Element = "p:2"` serendipity quads, no `Pressure Stabilization`: edge pressures pinned |
| `taylorhood_tri` | `StokesStabilized/tri_taylorhood.sif` (`ElmerGrid 1 2 Step.grd -triangles`) | the same on p:2 triangles |
| `taylorhood_tet` | `StokesStabilized/tet_taylorhood.sif` (the `stokes_tet` mesh) | the same on p:2 tetrahedra (reordered nodes), `Element Integration Points = "-tetra 5"` |
| `taylorhood_quad_full` | `quad_taylorhood.sif` with `Serendipity P Elements = False` | full p:2 quads, one bubble each condensed (`LCondensate`, no pressure pinning in a Stokes run); velocity only (the unpinned pressure has spurious modes) |
| `mini_quad` | `StokesStabilized/quad.sif` without `Pressure Stabilization` | `_Init0`'s MINI element `p:1 -quad b:3`, velocity bubbles condensed; norm 7.5549150E-02 |
| `mini_tri` | `StokesStabilized/tri.sif` without `Pressure Stabilization` | the same, `-tri b:1` |
| `mini_ns_tri` | `ns_tri` without `Pressure Stabilization` | Navier-Stokes on the MINI element, Picard then Newton after 3, bubbles condensed every iteration; norm 7.5920118E-02 |
| `mini_transient_tri` | `transient_tri` without `Pressure Stabilization` | the MINI element over 4 steps of 0.1, Newton after 2; `LCondensate`'s BDF1 time derivative with each element's bubble history (the run asks for BDF2); norm 7.4404965E-02 |
| `taylorhood_transient_tri` | `transient_tri` with `Element = p:2`, no `Pressure Stabilization` | Taylor-Hood Navier-Stokes over 4 BDF2 steps (`Default1stOrderTime`, the pinned pressures' mass rows zero); norm 3.7828652E-02 |
| `thermal_quad` | own (`ThermalCompress_vec` without its HeatSolve, `ElmerGrid 1 2 ThermalCompress_vec.grd`) | `Compressibility Model = Thermal` on the MINI element: the temperature an exported variable at its initial condition `293.5 - x`, gravity, pressure pinned on boundary 1; norm 4.4556902E-03 |
| `perfect_gas_quad` | own (`PerfectGasCompress_vec` without its HeatSolve, the `thermal_quad` mesh) | `Compressibility Model = Perfect Gas` (`Heat Capacity` 1000, `Specific Heat Ratio` 5/3, `Reference Pressure` 1e5), the same fixed temperature; norm 5.1831983E-05 |
| `axi_pipe` | `fem/tests/IncompressibleNSAxi` (`ElmerGrid 1 2 pipe.grd`) | axisymmetric Stokes, Pressure Stabilization on `p:1`, `Specific Load` `-1/r`; norm 1.3600543 |

The `taylorhood_*` and `mini_*` cases have p dofs, so their references are by dof
slot instead of node: `perm.dat` holds the `Perm` of the `.result` as `slot
row` lines and `solution.dat` the values as `slot u1 u2 [u3] p` lines, slot by
slot (`x(Perm(slot))`); `tests/taylorhood.rs` compares every slot.

The three `StokesStabilized` sifs are the stock ones with the iterative solve
replaced by `Linear System Solver = Direct` / `UMFPack` (the Rust side solves
with the default faer LU) and `Output File = "case.result"` added. Their MATC
inflow profile has no MATC engine in the Rust test, so the test hands the
solver a linear table through the inflow nodes' coordinates, which reproduces
MATC exactly there. In `ns_quad` the pressure level is a null space of the
`GradP` form, so the pressure is compared up to a constant.
