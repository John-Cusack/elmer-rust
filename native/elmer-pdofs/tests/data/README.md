# Pristine Elmer captures for `elmer-pdofs`

Each directory holds one case: the mesh as ElmerSolver read it
(`mesh.*`), the case's `case.sif`, and dumps from the pristine solver
(Rev 9f6af2f85, `/home/john/runs/elmer-fork/orig_src/build/fem/src/ElmerSolver`,
serial, `OMP_NUM_THREADS=1`). Every `.dat` file starts with a `#`
provenance header; the rest is Elmer's own output, unchanged.

| Case | Source | Solver | Notes |
|---|---|---|---|
| `poisson_pfem` | `fem/tests/PoissonPFEM` (`ElmerGrid 1 2 square.grd`) | the test's `Poisson.f90`, `Element = "p:2"` | 10x10 quads, 2 bodies; 441 DOFs = 121 nodes + 220 edges + 100 global bubbles; bandwidth 342 -> 85 |
| `varying_p` | `fem/tests/p-FEM_with_varying_p` (`ElmerGrid 1 2 eighth.grd`) | `BenchmarkShell`, 5 DOFs, `Element = "p:%bodywisep"` | 2x2 quads, p = 4 (body 1) and 8 (body 2); 189 of 289 slots; bandwidth 187 -> 187 (identity kept) |
| `sd_poisson_pfem` | `fem/tests/SD_PoissonPFEM` (same mesh as `poisson_pfem`, which it reads) | the same `Poisson.f90`, serendipity `p:2` | 341 DOFs (no bubbles), 8-point economic quad rule; bandwidth 233 -> 65; norm 2.23542169E-02 |
| `tri_p3` | lane-made: `poisson_pfem`'s `square.grd` meshed with `ElmerGrid 1 2 square.grd -triangles`, `Element = "p:3"`, no Reference Norm | the same `Poisson.f90` | 200 triangles; 961 DOFs = 121 nodes + 320 edges x 2 + 200 bubbles; 12-point p-reference simplex rule; bandwidth 852 -> 104 |
| `brick_p3` | lane-made: `fem/tests/QuadraturesBrick/cube.grd` with `Element Divisions` 3 3 3, `Element = "p:3"`, `Serendipity P Elements = False` | the same `Poisson.f90` | 27 bricks; 1000 DOFs = 64 nodes + 144 edges x 2 + 108 faces x 4 + 27 bubbles x 8; bandwidth 958 -> 692 |
| `brick_sd_p4` | lane-made: the `brick_p3` mesh (which it reads), `Element = "p:4"`, serendipity | the same `Poisson.f90` | 604 DOFs = 64 + 144 x 3 + 108 x 1; bandwidth 562 -> 399 |

`poisson_pfem` was run twice with these lines added to `Solver 1`:

```
Linear System Save = True
Linear System Save Perm = True
Linear System Save Prefix = String "pre"      ! or "post"
Linear System Save Slot = String "assembly"   ! "pre" run only
```

The `assembly` slot (DefaultFinishAssembly) dumps the system before
Dirichlet conditions; the default `solve` slot (start of DefaultSolve)
dumps it after `DefaultDirichletBCs`. Files:

* `pre_a.dat`, `pre_b.dat`: matrix (`row col value`, 1-based, CSR order)
  and right-hand side before Dirichlet conditions;
* `post_a.dat`, `post_b.dat`: the same after them;
* `perm.dat`: `Perm` (slot -> row, 1-based; identical in both runs);
* `solution.dat`: the solution from `Output File = "case.result"` (17
  digits), one value per slot, `x(Perm(i))`.

`varying_p` was run once (the `pre` run). Its matrix dump is 12 MB, so
only `perm.dat` and `pattern.dat`, the list-matrix pattern derived from
it (rows and columns of the 5x5 blocks, as ranges), are kept; the
shell kernel itself is not ported here.

`sd_poisson_pfem`, `tri_p3`, `brick_p3` and `brick_sd_p4` were run once, with the default `solve`
slot dump, and keep only `case.sif` (the run's), `perm.dat` and
`solution.dat` (from `case.result`); `tri_p3` and `brick_p3` also keep
their meshes.
