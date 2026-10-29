# elmer-guicompat

The ElmerGUI process boundary as executable checks. ElmerGUI never links the solver; it

1. probes `ElmerSolver -v` and `ElmerSolver_mpi -v` (the GUI's `args` list accumulates, so the `_mpi` probe is `-v -v -v`),
2. writes a case directory with `ELMERSOLVER_STARTINFO`, starts `ElmerSolver` with **no arguments**,
3. scrapes `Time:` and `ComputeChange ... (NRM,RELC)` lines (`mainwindow.cpp:6926-6978`, ported in `logparse`),
4. loads what the case wrote (`.ep` header, VTU through the same fields `vtkpost.cpp` reads).

This crate does the same against any solver executable given by absolute path: the pristine Fortran build, or the `elmer-cli` binaries.

```text
elmer-guicompat fixtures --edf-dir <ElmerGUI/Application/edf> --out DIR
elmer-guicompat check --solver /abs/ElmerSolver [--solver-mpi /abs/ElmerSolver_mpi] \
    --edf-dir <edf dir> --work DIR [--only ID]... [--native-only] [--strict] [--json]
elmer-guicompat parse-log FILE
```

Pristine Fortran needs `ELMER_HOME=<build>/fem/src ELMER_LIB=<build>/fem/src/modules` in the environment. For `elmer-cli` set `ELMER_FORTRAN_SOLVER` to delegate; `--native-only` strips it, so the executable must answer alone.

## Fixtures

Seven cases, one per shipped EDF (`heatequation`, `linearelasticity`, `navier-stokes`, `electrostatics`, `helmholtz`, `meshdeform`, `resultoutput`; the last is heat plus `ResultOutputSolve` writing VTU), on a 2x2-cell unit square. `Procedure` and `Variable` come from the EDF defaults with the GUI's own rewriting (`Name(n)` becomes `-dofs n Name`, `[Velocity:cdim ...]` evaluates `cdim`); check boxes left at their EDF default are not written, as in `SifGenerator::handleCheckBox`. Values are untyped, like the GUI writes them, so the `typing:` step resolves every keyword through `elmer-keywords` and fails on an unlisted keyword or a value that does not fit its type. An EDF `Type` tag that disagrees with the catalog is a *note*: it is a widget hint, not the SIF value type.

Classes: `heat`, `elasticity` and `output` are native targets; `unsupported` (Navier-Stokes, electrostatics, Helmholtz, mesh deformation) must be refused with exit 78 by a native-only solver, never solved. A solver exit of 78 is reported as `unsupported`; it fails the run only under `--strict` for native-target fixtures.

## What is not claimed

- Linux "Compile" is unimplemented in this ElmerGUI fork and is not modelled.
- The GUI itself was not launched: visual validation is outstanding for S3.
- Passing proves the process boundary and file formats, not numerical parity (the oracle owns that).

Tests need no solver (scripts stand in for refusals and crashes). Set `ELMER_GUICOMPAT_SOLVER` (and `_MPI`) to run the whole contract against a real one; `ELMER_GUICOMPAT_EDF_DIR` points at the EDF files.
