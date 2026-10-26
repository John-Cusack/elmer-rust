# lane-fhutiter unit ledger (17 units)

Lane `lane-fhutiter` (`fhutiter/` — the HUTI iterative-solver library):
4 C/header + 13 Fortran units. Upstream `9f6af2f85`; pinned nightly
`1.100.0 feaadeeac`; gfortran 13.3.0; ccache per ADR-045.

Machine verdicts from the shipped gate (`fortran_substitution_verdict` +
`has_built_object`, `proof/gate-verdicts.json` in the rustsmith T010
evidence dir) plus manual object verification where noted.
Status: `graded` = Rust port in `ports/hutiaux` + driver-level
differential green (see EVIDENCE.md); `skip-*` = recorded skip, never
forced; `refused-*` = gate refusal with reason, never forced.

| unit | status | reason |
|---|---|---|
| fortran:fhutiter/src/huti_aux.F90 | graded | gate ADMIT (scalars by ref + assumed-size/explicit-shape raw arrays, COMMON-free); object exists (`fhuti.dir/huti_aux.F90.o`); crate `hutiaux`, 12/12 symbols, differential byte-identical (EVIDENCE.md) |
| c:fhutiter/examples/ex1/dummymain.c | skip-no-built-object | gate admits (plain C) but `examples/` is not in the campaign build (no `add_subdirectory(examples)`): no object, scheduler-skip class (tool #107) |
| c:fhutiter/include/huti_fdefs.h | skip-helper | header: no compilable claim (tool #107 scheduler skip) |
| c:fhutiter/src/huti_fdefs.h | skip-helper | header: no compilable claim (tool #107 scheduler skip) |
| c:fhutiter/src/huti_intdefs.h | skip-helper | header: no compilable claim (tool #107 scheduler skip) |
| fortran:fhutiter/examples/ex1/globals_module.f90 | skip-no-built-object | gate admits (data-only module) but `examples/` unbuilt: no object |
| fortran:fhutiter/examples/ex1/huti-ex.F | skip-no-built-object | gate admits (example program; `own_matvec` has the proven F77-external profile) but `examples/` unbuilt: no object |
| fortran:fhutiter/examples/ex1/huti-ex.F90 | skip-no-built-object | same as `.F` (dual fixed/free sources of one example): no object |
| fortran:fhutiter/src/huti_bicgstab.F90 | refused-procedure-pointer | solver callbacks are `procedure(...), pointer` dummies (libgfortran realloc contract) in all 4 precision variants; unproven shape, never forced |
| fortran:fhutiter/src/huti_bicgstab_2.F90 | refused-procedure-pointer | same: 4/4 solver procedures carry procedure-pointer callbacks |
| fortran:fhutiter/src/huti_cg.F90 | refused-procedure-pointer | same: 4/4 solver procedures carry procedure-pointer callbacks |
| fortran:fhutiter/src/huti_cgs.F90 | refused-procedure-pointer | same: 4/4 solver procedures carry procedure-pointer callbacks |
| fortran:fhutiter/src/huti_gmres.F90 | refused-procedure-pointer | same: 4/4 solver procedures carry procedure-pointer callbacks |
| fortran:fhutiter/src/huti_qmr.F90 | refused-procedure-pointer | same: 4/4 solver procedures carry procedure-pointer callbacks |
| fortran:fhutiter/src/huti_tfqmr.F90 | refused-procedure-pointer | same: 4/4 solver procedures carry procedure-pointer callbacks |
| fortran:fhutiter/src/huti_interfaces.F90 | refused-interface-only | abstract-interface bodies only (`mv_iface_*`, `pc_iface_*`, ...); the built object defines no text symbols — nothing to substitute |
| fortran:fhutiter/src/huti_sfe.F90 | refused-procedure-pointer | 28 `huti_{s,d,c,z}_*` convenience wrappers take procedure-pointer callbacks + BLAS `INTERFACE` bodies external to the file (unprovable) |

Note on `has_built_object`: the shipped resolver reports `false` for
`huti_aux.F90` although `fhuti.dir/huti_aux.F90.o` exists — its DependInfo
fallback parses only the newer `CMAKE_DEPENDS_DEPENDENCY_FILES` section,
while Makefiles Fortran targets record the mapping in
`CMAKE_DEPENDS_CHECK_Fortran`. Recorded as a follow-up tool gap (no tool
change in this lane); substitution mechanics validated manually with the
identical `ld -r` composition in the grade driver.
