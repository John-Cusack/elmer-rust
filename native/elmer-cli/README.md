# elmer-cli

`ElmerSolver` and `ElmerSolver_mpi` of the Rust-native Elmer. Until the native
engine lands, they answer `-v` themselves and hand every case to the pristine
Fortran solver named by `ELMER_FORTRAN_SOLVER` (never linked, only `exec`ed).

| Variable | Meaning |
|---|---|
| `ELMER_FORTRAN_SOLVER` | path of the pristine serial `ElmerSolver` |
| `ELMER_FORTRAN_SOLVER_MPI` | path of `ElmerSolver_mpi` (default: the serial path) |
| `ELMER_NATIVE=1` | try the native engine first; refusal falls back to Fortran |
| `ELMER_SHADOW=1` | run Fortran as a child; native comparison is a stub |
| `ELMER_CLI_DELEGATED` | set on the delegated process; blocks re-entry |

Default is Fortran-first: with `ELMER_FORTRAN_SOLVER` set the case is `exec`ed
before anything is printed, so its stdout, cwd, environment and exit status are
the Fortran ones. Without it the native engine is asked, and a refusal exits 78.

| Exit | Meaning |
|---|---|
| 0 | solved |
| 1 | runtime or numerical failure |
| 2 | bad arguments, missing input file, unusable Fortran path |
| 78 | not supported natively and no Fortran fallback |

`-v` prints the pristine format plus one `Backend:` line. The ElmerGUI parsing
of `Time:` and `ComputeChange … (NRM,RELC)` lines is covered by
`tests/gui_compat.rs`; the capability manifest is described in
`docs/manifest.md`.

Not in scope: Linux "Compile" (unimplemented in this fork's ElmerGUI).
