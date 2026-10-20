# lane-fem-a evidence (batch 1)

Reference: Elmer @ `9f6af2f85`, fork `lane/fork`, `Debug`, defaults
(`-DRESET_ELMER_REVISION=FALSE -DELMER_FEM_REVISION=9f6af2f85`,
branch `devel`; `ELMER_FEM_COMPILATIONDATE` is build-date by nature).
Toolchain: `rustc 1.100.0-nightly (feaadeeac 2026-09-19)`.
ccache: shared `elmer-lanes-cache`, `CCACHE_BASEDIR=$LANE_ROOT`, sloppiness
per `data/elmer-lanes/ccache.env` (ADR-045).

## Backstop (pristine)

`ctest --output-on-failure -L quick` (the oracle invocation: 482 tests).
Log: `lane/ctest_backstop.log`. Result: **479/482, failures exactly
K=3** (`ConstantBCTemperature`, `ProfileBCTemperature`,
`ProfileBCTemperatureRobin`) — campaign baseline reproduced.
Failure set: `lane/backstop_failures.txt`.
Parallel-validation: `ctest -L quick -j16` on the pristine tree gives the
identical failure set (`lane/ctest_j16.log`) in ~7 min; all per-unit
suites use `-j16` on this evidence.

## Batch-1 units (11)

| unit | crate | true symbols | local differential |
|---|---|---|---|
| c:fem/src/view3d/LUDecomp.c | ludecomp | ludecomp, mtrinv | byte-identical stdout+stderr, 5 cases |
| c:fem/src/view3d/VectorUtil.c | vectorutil | RotateVector, RotateVectorTranspose, GetMatrixToRotateVectorToZAxis, CylinderNormal | byte-identical, incl. exit(1) discriminant path |
| c:fem/src/view3d/second.c | second | second | same-syscall functional |
| c:fem/src/f_stubs.c | f_stubs | rename_c (HAVE_F_ETIME/HAVE_F_FLUSH compile the rest out) | round-trip green |
| c:fem/src/CPUTime.c | cputime | cputime, realtime, cpumemory | same-syscall functional |
| c:fem/src/GFortranLocaleFix.c | gfortranlocalefix | elmer_fix_numeric_locale | identical (setlocale path + idempotence) |
| c:fem/src/binio/binio.c | binio | 16 × bin*_c + strerrorf_c | byte-identical, all fns + swap/EOF/errno paths |
| fortran:fem/src/minpack_kinds.F90 | minpack_kinds | (none — parameter-only module) | n/a (no symbols) |
| fortran:fem/src/ElementDescription.F90 | elementdescription | (none — facade) | n/a |
| fortran:fem/src/MeshUtils.F90 | meshutils | (none — facade) | n/a |
| fortran:fem/src/SolverUtils.F90 | solverutils | (none — facade) | n/a |

## Batch-2 units (5, local differential green, machine grade pending)

| unit | crate | true symbols | local differential |
|---|---|---|---|
| c:fem/src/umf4_f77wrapper.c | umf4_f77wrapper | 14 × umf4*/umf4_l_* | byte-identical vs vendored UMFPACK (di path) |
| c:fem/src/view3d/2d_4node.c | _2d_4node | 18 shape-fn symbols | byte-identical incl. libc-fprintf %g paths |
| c:fem/src/view3d/16node_quad.c | _16node_quad | 9 shape-fn symbols | byte-identical |
| c:fem/src/view3d/ShaftCull.c | shaftcull | 6 shaft symbols | byte-identical (Walk, planes, counting) |
| c:fem/src/view3d/ContourInteg.c | contourinteg | 5 contour symbols | byte-identical (prepare/eval/poly/counts) |

Machine grades (splice + rebuild + quick suite + matc probe + nm scan),
rows in `lane/store.db` run `fem-a-batch1`. Batch-1: **10/10 green** —
every suite fails exactly K=3, every matc probe prints `3`:

| unit | suite log | wall |
|---|---|---|
| c:fem/src/view3d/LUDecomp.c | suite_ludecomp.log | 198 s |
| c:fem/src/view3d/second.c | suite_second.log | 193 s |
| c:fem/src/f_stubs.c | suite_f_stubs.log | 181 s |
| c:fem/src/CPUTime.c | suite_cputime.log | 189 s |
| c:fem/src/view3d/VectorUtil.c | suite_vectorutil.log | 185 s |
| c:fem/src/GFortranLocaleFix.c | suite_gfortranlocalefix.log | 191 s |
| c:fem/src/binio/binio.c | suite_binio.log | 185 s |
| fortran:fem/src/ElementDescription.F90 | suite_elementdescription.log | 180 s |
| fortran:fem/src/MeshUtils.F90 | suite_meshutils.log | 180 s |
| fortran:fem/src/SolverUtils.F90 | suite_solverutils.log | 177 s |

Batch-2 (8 units) machine grades, rows in `lane/store.db` run `fem-a-batch1`:
- `c:fem/src/umf4_f77wrapper.c`: K=3, probe `3`, wall 170 s.
- `c:fem/src/view3d/2d_4node.c`: K=3, probe `3`, wall 175 s.
- `c:fem/src/view3d/16node_quad.c`: K=3, probe `3`, wall 174 s.
- `c:fem/src/view3d/ShaftCull.c`: K=3, probe `3`, wall 248 s.
- `c:fem/src/view3d/ContourInteg.c`: K=3, probe `3`, wall 259 s.
- `c:fem/src/view3d/ClipShadow.c`: K=3, probe `3`, wall 268 s.
- `c:fem/src/view3d/LineUtil.c`: K=3, probe `3`, wall 203 s
  (after inlining header-static `FunctionValue`; the extern form failed
  the link — recorded in the port notes).
- `c:fem/src/view3d/TriangleUtil.c`: K=3, probe `3`, wall 183 s (same fix).

Batch-2: **8/8 green**. Lane total machine-green: **18/18 portable
non-facade units** (10 batch-1 + 8 batch-2). Combined splice gate green (see below).

## Refined merge gate (panel corrections 1–3)

Selective member extraction (`lane/combine2.sh`) over all 14 non-facade
port archives: 99 unit symbols unique (zero STRONG port-level dups);
closure keeps 27 members into `combined_unit.o` (1908 strong symbols,
no duplicate); runtime members come from exactly one provider
(`ludecomp`); the 5 co-owned `__rustc_*` alloc shims are byte-different
across crates (audit in `stage_dry/byteid.txt`) and therefore included
from the provider only — never asserted identical. Log:
`lane/stage_dry/combine.log`.
- Held-out: no held-out workloads are frozen for the Elmer package
  (`heldout_workloads.json` is empty), so no HeldoutDivergence row is
  recorded — never a vacuous pass.

## Combined splice grade (SPLICE1)

`lane/grade_combined.sh`: the 14 non-facade units\u2019 selected members
spliced per object (provider runtime members only in the `ludecomp`
object \u2014 exactly one runtime provider per final binary); facade units
keep pristine objects (vacuous-port rule: their ports define no symbols).
Rebuild + `ctest -L quick -j16` + matc probe. Log:
`lane/combined_grade.log`, suite `lane/suite_combined.log`.
Result: **failures exactly K=3** (`lane/fails_combined.txt`,
`COMBINED-SAME-AS-BACKSTOP`), matc probe `3`
(`lane/probe_combined.txt`). The tree was restored to pristine after
(`*.pristine18` round-trip, smoke-tested: probe `3`, sample ctest pass).

## Store

`lane/store.db`: run `fem-a-batch1` (18 `gated` units with
OracleIntegrity/OracleParity/Differential rows), run `fem-a-ledger`
(137 `skipped` rows with ledger reasons + 14 `parked` rows for the
T008 scope-cut units, each with the round-2 deferral reason in
`spec_json`). Ledger and store agree on all 168 units; the cut itself is
recorded as a `decisions` row (`fixforward-triage`).

## Rebase onto merged main (2026-10-16)

Rebased `elmer/lane-fem-a` onto main at `d19e5b1` (merges #3 grid-io,
#6 fem-b-r2, #4 gui-c). Conflicts: `rust-toolchain.toml` add/add +
content — kept the lane's float-fix wording, extended to cite both
`elmer-fem-a` and `elmer-fem-b` `build-check.sh`. Port sources
(`crates/elmer-fem-a/ports/`, LEDGER.md) byte-identical after rebase.
One structural fix required by the rebase: main's root workspace
(grid-io) seized standalone port builds (`believes it's in a workspace`),
so added nested `crates/elmer-fem-a/Cargo.toml` workspace (19 port
members, lane-gui-c precedent); no root-manifest touch, no port edits.
Re-proof on the rebased tree, pinned `rustc 1.100.0-nightly feaadeeac`:
`cargo build --release` 19/19 archives; refined strong-symbol gate over
fresh archives — zero non-toolchain STRONG port-level dups (4188 dup
names, all compiler_builtins/std-monomorphized/weak-libm/alloc-shim
classes), unit true symbols (`ludecomp`, `mtrinv`, `RotateVector`,
`second`, `cputime` family, `rename_c`, `elmer_fix_numeric_locale`,
14 × `umf4*`) each defined-strong in exactly one archive; root
workspace `cargo test --release` green (7 passed, 0 failed).
Full splice+corpus re-grade (ctest backstop K=3, SPLICE, matc probe)
reruns on the grade host — `lane/` state is not in git.

## Tight-loop recraft (2026-10-20): `ludecomp` / `mtrinv`

Craft pattern: rustsmith **ADR-050** (the crc-style tight loop — hoisted
branching, primitive-int state, no hot-path allocation, narrow fast path
with delegation guards). Applied to the ported unit
`c:fem/src/view3d/LUDecomp.c` (crate `ludecomp`, symbols `ludecomp` +
`mtrinv`) and measured with this repo's own A/B bench (`bench/`, ADR-039
conventions: one corpus, rep-interleaved arms, fresh child per sample,
median + MAD, peak RSS via `wait4`, digest equality).

**This is an explicit, per-unit measurement — not a per-unit grade claim.**
ADR-045 keeps per-unit grades free of timing claims; the rows below are a
separate, reproducible bench of one unit, with the commands in
`bench/README.md`.

### Measured (this host, pinned nightly `1.100.0-nightly feaadeeac`)

| workload | before | after | speedup | verdict | peak RSS before -> after | digest |
|---|---|---|---|---|---|---|
| LU decompose n=256 | 255.947 ms | 227.187 ms | **1.13x** | faster | 2760 KiB -> 2772 KiB (+0.4%) | identical |
| matrix inverse n=256 | 187.181 ms | 179.213 ms | **1.04x** | faster | 2760 KiB -> 2748 KiB (-0.4%) | identical |
| LU decompose n=512 | 214.458 ms | 203.146 ms | **1.06x** | faster | 5840 KiB -> 5832 KiB (-0.1%) | identical |
| matrix inverse n=512 | 467.345 ms | 421.481 ms | **1.11x** | faster | 5816 KiB -> 5824 KiB (+0.1%) | identical |

15 reps (11 at n=512), arms interleaved, same corpus both sides. Machine
rows: `bench/ludecomp/results.json`, `bench/mtrinv/results.json`,
`bench/ludecomp-n512/results.json`, `bench/mtrinv-n512/results.json`
(plus `rows.md` per directory). The digest is FNV-1a over the result
matrix bytes: identical on both arms, so the speed numbers are for
bit-identical results.

Harness overhead (the corpus copy that the measured region needs, because
the port mutates its argument in place) is reported per row and is 0.8% of
the LU workload and 0.1% of the inverse workload; it is never subtracted
from the workload number.

### What changed (and why it is bit-exact)

* **LU decomposition** (`port_ludecomp`): the pivot candidate is held in a
  register instead of being re-read from the matrix on every candidate;
  the `j != i` test is invariant over the elimination k-loop, so the swap
  leaves the loop body (two loop shapes, one body each); the elimination
  inner loop walks row pointers, so the C index math `n*i + j` is done
  once per row instead of on every access.
* **Matrix inverse** (`port_mtrinv`): the inverse product wrote its
  partial sum back to `a(i,j)` on every k iteration (an n^3 store chain).
  The only value ever read back from `a(i,j)` is read on the *first*
  iteration, before any store, and no other iteration of the product reads
  `a(i,j)` afterwards — so the single final store is equivalent (n^3 stores
  become n^2). The arithmetic is the same operations in the same order, so
  the results stay bit-identical.

### Rejected after measurement (kept out on purpose)

Phase medians, n=256, 9 reps, same corpus (scratch study; every variant
produced an identical digest):

| phase | variant | median | verdict |
|---|---|---|---|
| LU | original | 1.28 ms | — |
| LU | branch hoist only | 1.21 ms | — |
| LU | **row bases + branch hoist (kept)** | **1.11 ms** | 1.15x |
| INV(U) | **original (kept)** | **3.58 ms** | — |
| INV(U) | row base for row i | 3.80 ms | 0.94x, rejected |
| INV(U) | row bases for both rows | 3.86 ms | 0.93x, rejected |
| INV(L) | original (kept) | 3.29 ms | — |
| INV(L) | row bases + cached diagonal | 3.27 ms | inside noise, not kept |
| INV(AP) | original | 7.06 ms | — |
| INV(AP) | **hoisted store (kept)** | **6.49 ms** | 1.09x |
| INV(AP) | split k-range into two branch-free bodies | 8.72 ms | 0.81x, rejected |
| INV(AP) | hoisted store + split k-range | 7.31 ms | 0.90x, rejected |
| INV(AP) | hoisted store + split k-range + row bases | 7.85 ms | 0.90x, rejected |

The lesson recorded for the next port: **the craft is not a formula, it is
a measurement.** "Hoist the branch" and "walk a row pointer" are not
guaranteed wins in compiled Rust — the compiler already unswitches
loop-invariant branches and strength-reduces the index math, and a
hand-rolled form can block the unroller that the original shape enabled.

### Lane-fem-a bench harness

`crates/elmer-fem-a/ports/ludecomp/benches/ludecomp_bench.rs` + the
shared driver `bench/harness` (see `bench/README.md`). The port crate
gained `rlib` in its `crate-type` so the bench target can link the
package's own lib; the `staticlib` archive — the artifact the campaign
splices — is unchanged.
