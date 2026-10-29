# ADR 001: S0 oracle design and ABI campaign retirement

- Status: accepted.
- Date: 2026-10-29.
- Decides: report §5.2 items 1-4 and 6 (runner, golden capture, grader, tolerance inputs), D4 (ABI retirement), D7 as applied to new crates, workspace layout, and container pinning.

## Context

- The strategy review concludes the per-file ABI port campaign cannot delete Fortran by construction (pointer-graph core, 81% of solver logic still in Fortran accessors, untouched core).
- The captain accepted all review recommendations, including a Rust-native Elmer graded by an external oracle built first (stage S0).
- This task ships S0 in reviewable increments; this ADR records the oracle design and the retirement so later increments stay consistent.

## Decision D4: retire the ABI campaign as product code

- No new ABI-port slices are opened after this point.
- The staged `.sif` cases and the documented defect registry (UB re-entry in FixTangentVelo and CoilSolver, uninitialised reads in FDiffusion3D, dangling pointer at `ElectricForce.F90:172`) are harvested as oracle inputs: the `.sif` cases join the golden pool and the defects seed the known-divergence registry.
- The `crates/` and `fem-modules/` trees are frozen evidence; the native workspace never modifies or depends on them.

## Workspace layout

- The Rust-native rewrite lives in `native/`; the oracle crate `native/elmer-oracle/` is its own Cargo workspace with its own lockfile, so no shared manifest can smuggle a candidate crate into the oracle build.
- A future `native/` swarm workspace (candidate lanes) explicitly excludes `elmer-oracle`; the root workspace (frozen ABI crates) references neither.
- Rationale over putting the oracle in a shared `native/` workspace: a 16-lane swarm resolves `elmer-*` members together, which would void the crate-local isolation proof below the moment a lane is added.
- New crates are licensed LGPL-2.1-or-later (D7); existing crates keep their licences.

## Oracle crate (`elmer-oracle`)

- Responsibility split: `runner` (direct execution), `norm` (ComputeNorm replication plus the `CompareToReferenceSolution` verdict), `result` (binary `.result` reader), `grade` (field comparison plus mutation suite), `manifest` (hashed golden manifests), `main` (thin CLI).
- The crate has no path or git dependencies and no `elmer-*` dependency other than itself; `tests/isolation.rs` enforces this over the manifest and lockfile on every `cargo test` run.
- External grading only: candidate `TEST.PASSED` lines are ignored by the grader (they are still recorded as diagnostics).

## Direct runner

- The runner replays each test's `runtest.cmake` outside ctest: copy the test source tree to a run dir, run the preparation tools with the same binary variables the harness passes, then run `ElmerSolver` with the `RUN_ELMER_TEST` environment (`ELMER_HOME`, `ELMER_LIB` with extra module dirs, `PATH` prefix, `OMP_*` defaults, inherited `HWLOC_*`).
- Only the declarative subset is covered: `execute_process(COMMAND <known-bin> args...)`, `RUN_ELMER_TEST`, `EXECUTE_ELMER_SOLVER`, `EXECUTE_ELMER_SOLVER_MPI`.
- Anything else (custom macros, conditionals, file comparisons, result-variable verdicts) makes the plan `custom` with the offending construct named; the runner refuses to run it rather than silently skipping steps.
- Measured on the pristine tree: 1042 of 1073 test dirs plan declaratively (97%).
- The 31 refusals are genuinely exotic (adaptive-remesh rerun loops, FETI template generation, multi-corner viewfactor matrices, restart staging into subdirectories), matching the special classes that §4.1 grades differently anyway.
- Known gap: `file(COPY ...)` staging inside `runtest.cmake` (e.g. the `reload` test) is refused today; supporting it as a `Stage` step is follow-up work, not silent coverage.

## Norm replication

- `compute_norm` replicates the serial `ComputeNorm` branches exactly: degree 0/1/2/general, ascending index order with plain `+=`, strided per-component accumulation scaled by `NormDofs*n/Dofs`, nodal truncation, and `Norm Permutation` compaction.
- Two semantics were re-derived from the Fortran rather than assumed: an explicit `Components` list sets `NormDofs` to its length (it is not just a mask), and strided sections sum to the end of the array while only the scale uses the truncated length.
- Residual risk: the bitwise instrumentation-invariance check (§4.3 item 2) needs the grader's 17-digit recomputation to equal the solver's internal norm exactly, which holds only if gfortran's `SUM` compiles to ascending sequential accumulation.
- The S0 item-4 pilot (bitwise validation on every captured test) is the empirical gate for this risk, with the report's kill criterion (>2% mismatch escalates) applying unchanged.

## Binary `.result` reader

- The reader implements the `SaveResult` layout: ASCII `BINARY 3.<endian>` header, NUL separator, per-step `Time:` blocks, NUL-terminated variable names, Perm blocks (NULL as size 0, shared as size -1 with an absolute file offset), and native-endian doubles in Perm order.
- Cross-endian files are refused rather than misread; the host marker is checked against the header.
- Values stay in Perm order in the parsed struct; node-space mapping through the Perm pairs lives in the grader, where edge-DOF sign handling will also live.

## Containers and Fortran arms

- Docker is available on this host, so the spec's preference stands: one pinned image (`native/docker/Dockerfile.elmer-oracle`, Ubuntu 24.04 plus the Elmer toolchain) serves all arms.
- `native/oracle-arms/build-arm.sh` configures and builds one arm (`release`, `fma`, `openblas`, `listusage`) from pristine sources at the upstream pin, refusing dirty trees, wrong pins, existing build dirs, and a third concurrent host build.
- Every arm writes `versions.txt` (compilers, CMake, MPI, BLAS/LAPACK libraries, CPU, image digest) and the golden data commits reference it.
- The `listusage` arm carries the oracle-only `-DDEVEL_LISTUSAGE` flag (D6, no source change).
- The one-line `ListFind`-miss log patch (D6, source change) is not applied yet; it lands as its own change with the numerical-inertness proof the decision requires.
- The existing campaign Debug build remains the reference arm; the oracle never modifies the campaign fork.

## Consequences

- S0 follow-ups (golden capture, σ_F tables, capability matrix, O6 restart proof, O7 harness, held-out seal, benchmarks, harvest) build on these mechanics without changing them.
- Grader changes re-run the mutation suite; tolerance-table changes need an ADR with captain sign-off, never a port PR.
