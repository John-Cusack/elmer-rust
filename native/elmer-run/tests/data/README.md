# elmer-run test data

Each directory is one quick test from `elmer/fem/tests/<Test>`: the pristine
`case.sif` and the mesh DB that test's `runtest.cmake` generates with
ElmerGrid (the pristine C preprocessor, run once to make these files):

| Directory | Test | ElmerGrid command |
|---|---|---|
| `minimal_example/triangle` | MinimalExample | none (mesh is in the test) |
| `thermal_pipe_2d/OnePipe` | ThermalPipe2D | `ElmerGrid 14 2 OnePipe.msh -autoclean` |
| `coordinate_scaling/square` | CoordinateScaling | `ElmerGrid 1 2 square` |
| `wildcardsize/square` | wildcardsize | `ElmerGrid 1 2 square` |

The tests build each `case.sif` by hand (`tests/support/case.rs`) until
`elmer-sif` lands; the `.sif` files are kept as the reference text.
