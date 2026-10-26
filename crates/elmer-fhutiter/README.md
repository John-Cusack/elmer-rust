# lane-fhutiter: 1 fhutiter port graded + 16 recorded skips/refusals

Lane `lane-fhutiter` (`fhutiter/` — the HUTI iterative-solver library:
4 C/header + 13 Fortran units). Upstream `9f6af2f85`; pinned nightly
`1.100.0 feaadeeac`; gfortran 13.3.0; ccache per ADR-045.

## Units: 1 graded, 16 recorded (machine verdicts in LEDGER.md)

- **graded 1** (`hutiaux`): `huti_aux.F90` — 12 auxiliary routines in 4
  precisions (dummy preconditioner copy, PRNG fill, LU solve). Driver-level
  differential byte-identical (see EVIDENCE.md).
- **refused-gate 9**: 7 solver files (procedure-pointer callbacks,
  unproven libgfortran realloc shape) + `huti_interfaces` (abstract
  interfaces only, no code) + `huti_sfe` (28 wrappers with callbacks +
  external BLAS interfaces). Never forced.
- **skipped 7**: 3 headers (no compilable claim) + 4 example units
  (`examples/` not in the campaign build — no objects).

## Layout

- `Cargo.toml` — nested lane workspace (`ports/hutiaux`).
- `ports/hutiaux/` — the port (staticlib, 12 `__huti_aux_MOD_*` exports).
- `grade/` — `grade_hutiaux.sh` + `probe_hutiaux.F90` (driver-level
  differential + strong-symbol gate; logs beside the script).
- `LEDGER.md` — per-unit verdicts + reasons.
- `EVIDENCE.md` — corpus, findings, merge-gate position, timings.
- `build-check.sh` — pinned-toolchain guard (run by the grade driver).
