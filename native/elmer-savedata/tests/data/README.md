# Reference data

`two_cell/` is a 3 x 1 rectangle of two bilinear quadrilaterals (widths 1 and 2, so the
nodal mean and the volume-weighted mean differ), two bodies, and the sections of
`case.sif`/`ops2/case.sif`. `golden/` holds the files the pristine ElmerSolver (Version
26.2-devel, Rev 9f6af2f85, `build/fem/src/ElmerSolver`) wrote when run serially on that
case with `ELMER_HOME` and `ELMER_LIB` pointing at its build tree:

- `two_cell/case.sif`: transient, two BDF(1) steps of 0.5 s with a heat source in body 2.
  `scalars.dat` has 38 columns (statistics, integrals, counts, point values),
  `reader.dat` reads the exported `Temperature Integral` back and sets its norm,
  `line.dat` is the SaveLine of the two flagged boundaries.
- `ops2/case.sif`: the steady solution, the remaining operators, `Line Marker`,
  `File Append`, `Live Graph` (`.csv`, `.marker`), `Partition Numbering`, `Show Norm Name`,
  `Reference Values`, `Cost Function Index`, `Optimize Node Ordering = False`,
  `Skip Boundary Info` and a second `Save Mask`.

The tests drive the hooks from the pristine run's own temperatures (a stand-in for
`HeatSolver`, `Replay` in the tests), so the SaveData columns are compared to
rounding. The transient heat solution itself is not compared here. Elmer's first step of
`two_cell` (node 2: 18/11) is reproduced exactly by a backward-Euler step that starts from
the initial field with the Dirichlet values (0 and 10) already in place; the heat solver
under `elmer-solver` starts from the plain initial condition and gets -0.788. That is a
question for the solver framework, not for SaveData. The steady case is compared: heat and
Elmer agree (`heat_run.rs`).

The tests parse the real `case.sif` with `elmer-sif` (`harness::sif_model`). Two things it
hands over differently from Elmer, worked around in the tests and worth a look in the
front end: `Save Coordinates(1,2) = Real 2.9 0.9` comes out as a 2x2 array holding two values
(`take_reals` returns the element count where the row count belongs), and quoted strings keep
their case where Elmer's `ListAddString` lowercases them (`Comment = String "... SaveData ..."`).

The `.names` files carry the run's date and the full path of the data file; the tests
compare them without those two lines.
