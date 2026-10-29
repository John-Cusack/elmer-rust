# Reference data from the pristine ElmerSolver

Every directory holds an Elmer mesh (`mesh.*`), the `.sif` that was run
(`case.sif`; the ModelPDE ones are `lin.sif` and `evol.sif`) and node-ordered
solutions `solution.dat` (`node value` per line, 17 significant digits) from
`build/fem/src/ElmerSolver`, run serially on copies of the upstream tests with
one added line, `Output File = "case.result"` in `Simulation`, and with the
`Post File` line dropped. Meshes come from `ElmerGrid 1 2 <name>` as each
test's `runtest.cmake` does. The `.result` values are in node order
(`Values(Perm(i))`), not perm order. Norms are `sqrt(sum(x^2)/n)`.

- `adv_diff1/`: `fem/tests/adv_diff1` unchanged otherwise. 408 elements, BDF2,
  no convection. Norm 0.42696309972986873.
- `adv_diff2/`: `fem/tests/adv_diff2` unchanged (default residual-free
  bubbles; 404). Norm 0.8766881677902325. Variants of the same case:
  `nb_solution.dat` with `Bubbles = False` (0.8691287700487503) and
  `st_solution.dat` with `Simulation Type = Steady State` (0.9980513902691243).
- `adv_diff2_tri/`: `adv_diff2` on a triangle mesh (`Triangles = True` in
  `Step.grd`) with `Stabilize = True` instead of the default bubbles. Norm
  0.8844887712095778.
- `adv_diff_soret/`: `fem/tests/AdvDiffSoret` (`case.sif`, solutions
  `temp.dat`, `conc.dat`) and `AdvDiffSoretAbsMass` (`case_abs.sif`,
  `conc_abs.dat`). One change: the equation's `Convection = "Constant"` is
  `Conc Convection = String "Constant"`, so the heat solver sees no convection
  (its velocity is zero in the case anyway; the pristine norms are unchanged,
  351.2003227031049 and 0.7638989898318287 for relative mass). Mesh `strip`.
- `model_pde/`: `fem/tests/ModelPDE` with `Element = "p:2"` removed (linear
  bilinear 404 elements) and the post-file line replaced by `Output File`.
  `lin_solution.dat` is the steady case (norm 1.7972256169743914; the quick
  test's own p:2 reference norm is 0.923486039). `evol.sif` /
  `evol_solution.dat` add `Simulation Type = Transient` (BDF2, 6 steps of
  0.05), `time derivative coefficient = 2.0` and `reaction coefficient = 0.5`
  (norm 1.7550196309497508).
  The quick tests themselves, with their `Element = "p:2"`: `p2.sif` is
  `fem/tests/ModelPDE` and `sd_p2.sif` `fem/tests/SD_ModelPDE` (serendipity),
  each with the post-file line replaced by `Output File`. `p2_solution.dat` /
  `sd_p2_solution.dat` hold `slot value` per line, the value of every Perm
  slot (nodes, edges, bubbles) in slot order, and `p2_perm.dat` /
  `sd_p2_perm.dat` the `slot row` Perm (1681 dofs, norm 0.9234860388553119;
  1281 dofs, norm 1.057859263392748).
