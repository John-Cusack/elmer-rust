# Pristine-binary fixtures for `tests/pristine.rs`

`mesh2d/` and `mesh3d/` are captured runs of the unmodified `ElmerSolver`
(`26.2-devel`, rev `9f6af2f85`, serial, `OMP_NUM_THREADS=1`) on
`case.sif`. `Mesh2MeshSolver` runs `InterpolateMeshToMesh` from the solver
mesh `src` to the main mesh `dst`.

| File | What it is |
|---|---|
| `gen2d.py`, `gen3d.py` | write the `src/` and `dst/` meshes (seeded jitter, deterministic) |
| `case.sif` | the case that was run |
| `pristine_run.log` | its stdout (`Points not found: N (found M)`) |
| `expected_tempa.txt` | `TempA` on `dst` after the run, one `repr` per node; `-1000` marks nodes the Fortran left untouched |
| `src_tempa.txt` | `TempA` on `src`, dumped by the second `ResultOutputSolver` of the same run |
| `extract.py`, `extract2.py` | turn the binary `.vtu` files into the two text files above |

To regenerate: run `gen*.py` in a scratch directory, copy `case.sif` next to
the `src/` and `dst/` folders, run `ElmerSolver case.sif` there (the VTU files
land in `dst/case_t0001.vtu` and `src/src_case_t0002.vtu`), then extract.
