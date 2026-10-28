# T020 lane ledger (solver-entry pilot, report §8)

Scope: `fortran:fem/src/modules/Poisson.F90` at the plugin boundary
(T020 solver-entry gate shape, rustsmith ADR-057) + the `nm` symbol
surface. Statuses: graded = Rust cdylib in `ports/poisson` with
byte-identical differentials; skipped-* = recorded skip, never forced.

| unit | status | reason |
|---|---|---|
| fortran:fem/src/modules/Poisson.F90 | graded | T020 entry shape: F77-global `poissonsolver_` entry, derived dummies opaque under the accessor contract, `CONTAINS` internals proven non-global (`nm`, rustsmith `tests/t020_phantom_nm_proof.sh`). Port crate `ports/poisson`: Rust cdylib exporting exactly the reference `.so` dynamic symbols (`nm -D` equality); component access + all Elmer calls through generated-style Fortran `BIND(C)` wrappers compiled against the built `.mod` files; SAVE semantics replicated. Differentials: 2 distinct `.sif` configs byte-identical (stdout/stderr/scalars.dat modulo wall-clock lines) + one repeated-invocation run + one rerun; branch-liveness control included. See T020-EVIDENCE.md |
