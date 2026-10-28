# T020 evidence (solve_cmplx via the KIND= fix, report §8)

Reference: ElmerCSC/elmerfem @ fork `c0b4bc6`.
Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)` + gfortran
13.3.0. Grade driver: rustsmith `tests/femcore_cmplx_proof.sh`
(T020-CMPLX PROOF GREEN).

## Graded green (1 unit)

| unit | crate | true symbols | corpus |
|---|---|---|---|
| fortran:fem/src/solve_cmplx.F90 | solvecmplx | `solvelapack_cmplx_` (only `T`) | 5 cases byte-identical stdout+stderr+exit status |

Corpus (clean-room driver `probe_cmplx.F90`): 2×2 dense complex solve;
`N=0` early return; singular zero matrix (`ZGETRF` info print, then
`ZGETRS`); 1×1 complex; 3×3 dense with non-trivial phases. The
reference stages a clean-room `Types` stub with the identical
`dp = SELECTED_REAL_KIND(12)` definition (the real `Types.F90`
unconditionally `USE`s Lua; both arms share the stub, so the
differential is self-consistent by build).

## Refined strong-symbol gate

- Port defines exactly `solvelapack_cmplx_` as `T`; the defined-name
  intersection with `solve_cmplx.F90.o` is exactly the replaced symbol
  → zero STRONG port-level dups; whole-file portable.
