# T020 lane ledger (solver-entry pilot, report §8)

Scope: `fortran:fem/src/solve_cmplx.F90`, admitted by the T020 `KIND=`
legacy-parse fix (rustsmith ADR-057 §1.2). Statuses: graded = Rust port
in `ports/solvecmplx` with byte-identical corpus; skipped-* = recorded
skip, never forced.

| unit | status | reason |
|---|---|---|
| fortran:fem/src/solve_cmplx.F90 | graded | T020 `KIND=` fix: legacy `COMPLEX(KIND=dp)` dummies parse (was: mistaken for initialization). Port crate `ports/solvecmplx` (clones the merged `solvereal` template: ZGETRF/ZGETRS, `N<=0` early return, list-directed `info` diagnostics). 5-case corpus byte-identical (stdout+stderr+exit); strong-symbol gate exactly `solvelapack_cmplx_` as T, whole-file portable. Proof: rustsmith `tests/femcore_cmplx_proof.sh`. See T020-EVIDENCE.md |
