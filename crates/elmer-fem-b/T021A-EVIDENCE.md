# T021-A evidence (solver-entry ports, lane rs-elmer-t021a)

Reference build for all proofs below: serial Elmer tree
`/home/john/firstmate/data/rs-elmer-t020/scratch/elmer-build`
(campaign configure minus MPI; both arms share it, so every differential
is self-consistent). Toolchain pins asserted by every proof script:
gfortran 13.3.0, rustc 1.100.0-nightly (feaadeeac 2026-09-19).
Elmer checkout (mesh inputs only): `/home/john/runs/elmer-fork`.

Proof logs: lane `scratch/proof-<name>/` (per-config `run*/stdout.txt`,
`stderr.txt`, `scalars.dat` where applicable, `ref.syms`/`port.syms`).
Each proof rebuilds its cdylib against the reference `fmodules`
(`ELMER_MODDIR`), checks `nm -D` symbol equality, runs the shadow-dir
differentials (`ELMER_LIB` shadow with the port `.so` swapped in under
its own name, everything else symlinked), and halts on any divergence.

## allocatesolver (`ports/allocatesolver`, proof-allocatesolver.sh)

- Unit `fortran:fem/src/modules/AllocateSolver.F90` (18 code lines).
- Symbol surface (`nm -D --defined-only`, both arms):
  `T allocatesolver_`, `T allocatesolver_init_`. Diff empty.
- Config A (default; `_init` appends `Exec Solver = never` at the list
  tail, after setup-time `AddExecWhenFlag` cached `ALWAYS`, so the main
  entry still runs once): stdout/stderr byte-identical modulo
  wall-clock lines (`FINISHED AT`/`STARTED AT`/`Elapsed REAL time`
  filtered — the same filter class as T020).
- Config B (transient BDF, 2 timesteps, `Exec Solver = always`): main
  entry runs twice (in-process repetition). Byte-identical + full-process
  rerun byte-identical.
- Config C (`Exec Solver = never` in the `.sif`: setup-time `NEVER`,
  main entry never dispatched). Byte-identical.
- Liveness: the `Info` line `This solver does nothing...` prints 1x in
  A, 2x in B, 0x in C, in both arms.
- Proof wall ≈ 1 s. PROOF GREEN.

## divprojection (`ports/divprojection`, proof-divprojection.sh)

- Unit `fortran:fem/src/modules/DivProjection.F90` (10 code lines, empty
  body apart from `USE`s).
- Symbol surface: exactly `T divprojection_`. Diff empty. The port
  compiles no Fortran (no Elmer calls in the original); the version
  script alone enforces the surface.
- Config A (steady) + config B (transient, 2 timesteps): stdout/stderr
  byte-identical modulo wall-clock lines + full-process rerun of B
  byte-identical.
- Control C (`Exec Solver = never`): byte-identical; the
  `SingleSolver: Attempting to call solver` dispatch lines are present
  in A and absent in C, in both arms (a missing/unresolvable entry
  aborts the run, so exit 0 in both arms proves dispatch).
- Proof wall ≈ 1 s. PROOF GREEN.

## elementsizesolver (`ports/elementsizesolver`, proof-elementsizesolver.sh)

- Unit `fortran:fem/src/modules/ElementSizeSolver.F90` (85 code lines).
- Symbol surface: exactly `T elementsizesolver_`. Diff empty.
- Config A (defaults: no `Element Size Exponent` → `Power = 1.0`, no `No
  Integration Weight` → `DetJ`-weighted): stdout/stderr byte-identical
  modulo wall-clock lines.
- Config B (transient BDF, 2 timesteps, `Element Size Exponent = 2.0`,
  `No Integration Weight = True`): byte-identical + full-process rerun
  byte-identical.
- Liveness: `Minimum/Maximum Element Size` `ES12.4` reports differ
  across A/B in both arms (exponent + weighting branches change the
  numerics).
- Notes: per-call allocate/deallocate mirrored (the `SAVE`d
  `AllocationsDone` the original sets-but-never-reads is not
  replicated); allocation failure calls the same
  `Fatal('PoissonSolve', ...)` (bug-compatible program name).
- Proof wall ≈ 1 s. PROOF GREEN.

## divergencesolver (`ports/divergencesolver`, proof-divergencesolver.sh)

- Unit `fortran:fem/src/modules/DivergenceSolver.F90` (162 code lines).
- Symbol surface: `T divergencesolver_`, `T divergencesolver_init_`.
  Diff empty.
- Config A (cartesian, MATC linear velocity IC owned by an
  `AllocateSolver` solver, defaults): stdout/stderr/scalars.dat
  byte-identical modulo wall-clock lines.
- Config B (transient BDF 2 steps, `Relative Divergence = True`,
  `Constant Bulk Matrix = True` — exercises both the
  `DefaultUpdateEquations` and the `DefaultUpdateForce` paths across the
  two timesteps): byte-identical + full-process rerun byte-identical.
- Config C (axi-symmetric coordinates: the `CSymmetry` x-weighting
  branch): byte-identical.
- Liveness: `Result Norm` (`Level=4`, deterministic) and the saved
  Divergence min/max differ across A/B/C in both arms. The
  nondeterministic `Level=5` `RealTime` reports are structurally
  suppressed in both arms via the standard `Max Output Level = 4`
  keyword (the deterministic `Level=4` line still prints and is
  compared).
- Incidentally verified while staging: the square test mesh carries bulk
  bodies 1–4, so an `Initial Condition` run must declare all four `.sif`
  bodies (single-`Body 1` + IC reads `Bodies(3)` out of bounds — an
  upstream footgun, deterministic SIGSEGV; all IC proofs in this lane
  declare 4 bodies, matching the T020 Poisson configs).
- Proof wall ≈ 1 s. PROOF GREEN.

## elementstats (`ports/elementstats`, proof-elementstats.sh)

- Unit `fortran:fem/src/modules/ElementStats.F90` (287 code lines).
- Symbol surface: exactly `T elementstats_`. Diff empty.
- Config A (defaults: bulk + boundary sweeps): stdout/stderr
  byte-identical modulo wall-clock lines.
- Config B (`Create Histogram` with 5 intervals, `Skip Boundaries`,
  huge `Element Warn Size` tripping the `Warn` path but not `Fatal`):
  byte-identical + full-process rerun byte-identical.
- Negative control C (`Stop Below Warn Size` with the same threshold):
  both arms abort via `Fatal` with byte-identical stdout/stderr
  (exit 1 both arms).
- Liveness: `Histogram:` sections present in B, absent in A, in both
  arms. Reference accessor recipe matched to the Elmer module build
  (`-fPIC -fallow-argument-mismatch`, from the modules'
  `flags.make`); the legacy-MPI-idiom mismatch errors without it.
- Proof wall ≈ 1 s. PROOF GREEN.

## energyrelease (`ports/energyrelease`, proof-energyrelease.sh)

- Unit `fortran:fem/src/modules/EnergyRelease.F90` (96 code lines).
- Symbol surface: exactly `T releaseratesolver_`. Diff empty. (Campaign
  invocation is `Procedure = "EnergyRelease" "ReleaseRateSolver"` —
  library vs entry names differ; the proof uses the campaign form.)
- Config A (linear MATC source fields): stdout/stderr byte-identical
  modulo wall-clock lines.
- Config B (transient BDF, 2 timesteps, quadratic MATC fields):
  byte-identical + full-process rerun byte-identical (SAVE work arrays
  persist with identical lifetimes in both arms).
- Liveness: the `2G =` fracture report differs across A/B in both
  arms. Fields staged via three owning `AllocateSolver` solvers (one
  per field, distinct `Equation` names — the loader rejects duplicate
  solver equation names) + MATC initial conditions on all four bodies.
- Notes: allocation failure calls the same `Fatal('PoissonSolve',
  ...)` (bug-compatible program name); `Norm` is declared-but-unused in
  the original and absent here.
- Proof wall ≈ 1 s. PROOF GREEN.

## Campaign backstop

Fresh tree @ `9f6af2f85`, campaign configure (Debug, `HWLOC_COMPONENTS=-gl`,
`--allow-multiple-definition`), all six T021-A ports spliced
(`<Name>.so` ← port cdylib, rebuilt against the campaign `fmodules`;
`ld -r` composition not needed — whole-`.so` replacement, the T020
recipe), full rebuild, `ctest -L quick`: **479/482, failures exactly
K=3 upstream segfaults #80/#536/#537** (log: lane `proof/`).
