# Per-function reference data from the pristine ElmerSolver

`pristine/*.txt` hold basis values captured from the pristine, unmodified
`ElmerSolver` (`/home/john/runs/elmer-fork/orig_src/build`, serial) while it ran
copies of the `ElementUnitTest_h1/p1/edge/edge_piola/edge_quadratic` quick
tests. `capture/dump.py` is a gdb Python script: it breaks on
`__eleminfo_MOD_elementinfo` (full DWARF info is present), and for every call
made from the tests' `LocalMatrix` it records the element, its node numbers and
coordinates, `(u,v,w)`, and after the call returns `detJ`, `Basis`,
`dBasisdx`, `EdgeBasis` and `RotBasis` read through the gfortran array
descriptors, at full double precision. Nothing in the Elmer tree was built or
modified.

* `h_*`/`h2d_*`: p-element H1 bases. The tests' `SingleEdgeElements` solver was
  run with `Element = "p:N"` in the sif (it only supplies its own element
  string when the keyword is absent) and, for these files, the `EdgeBasis` /
  `RotBasis` arguments were nulled in the call frame (`NOEDGE=1`) so
  `ElementInfo` returns the H1 `Basis(200)`/`dBasisdx` instead of the edge
  functions. `h_*` use the 3D unit mesh (808/706/605/504), `h2d_*` a 2D mesh
  (404/303); `_id`/`_sh` are identity/shuffled global node numbering (this also
  makes Elmer pick tetra type 2), `_sd` the serendipity scheme.
* `e_*`: the three edge configurations (`edge` = lowest order, non-Piola;
  `edge_piola`; `edge_quadratic`), identity and shuffled numbering.

Line format: see `tests/h1_vs_pristine.rs` / `tests/edge_vs_pristine.rs`
(whitespace separated numbers; floats printed with Python `repr`, i.e. exact).
`capture/gen_cfg.py` builds the case directories, `mkfix*.py` turns the gdb
JSON lines into these files (H1 files keep every 5th Gauss point to limit size).
The scripts contain the absolute paths of the machine they ran on.
