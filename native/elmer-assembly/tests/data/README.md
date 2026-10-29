# Pristine Elmer captures for `elmer-assembly`

Each directory holds one case: the mesh as ElmerSolver read it
(`mesh.*`), the case's `case.sif`, and `Linear System Save` dumps from the
pristine solver (Rev 9f6af2f85,
`/home/john/runs/elmer-fork/orig_src/build/fem/src/ElmerSolver`, serial,
`OMP_NUM_THREADS=1`). Every `.dat` file starts with a `#` provenance header;
the rest is Elmer's own output, unchanged.

| Case | Source | Solver | Notes |
|---|---|---|---|
| `minimal_example` | `fem/tests/MinimalExample` | HeatSolve | one 303 triangle |
| `thermal_pipe_2d` | `fem/tests/ThermalPipe2D`, mesh from `ElmerGrid 14 2 OnePipe.msh -autoclean` | HeatSolve | 3 bodies, bandwidth 889 -> 60 |
| `hex_elasticity` | lane-made `box.grd` (2x1x1 subcells, 5x2x2 808 elements), `ElmerGrid 1 2 box.grd` | StressSolve, 3 DOFs | coordinate-dependent BC |
| `hex_heat_partial` | same mesh | HeatSolve | only body 1 active, `Linear System Symmetric = True` |

Each case was run twice with these lines added to `Solver 1`:

```
Linear System Save = True
Linear System Save Perm = True
Linear System Save Prefix = String "pre"      ! or "post"
Linear System Save Slot = String "assembly"   ! "pre" run only
```

The `assembly` slot (DefaultFinishAssembly) dumps the system before
Dirichlet conditions; the default `solve` slot (start of DefaultSolve)
dumps it after `DefaultDirichletBCs`. The two runs have the same Perm and
pattern. Files:

- `perm.dat`: node, 1-based slot (0 = not in the system);
- `pre_a.dat`, `pre_b.dat`: matrix (`row col value`, 1-based, CSR order)
  and rhs before Dirichlet;
- `post_a_changed.dat`: the entries of the post-Dirichlet matrix that differ
  from `pre_a.dat`;
- `post_b.dat`: rhs after Dirichlet.

Values are gfortran list-directed (17 significant digits), so they
round-trip to the exact `f64`. `PrintRHS` writes `0.0` for any
`|v| <= tiny(v)`.
