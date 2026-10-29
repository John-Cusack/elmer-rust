#!/bin/bash
# T022-A KESolverLegacy solver-entry proof: the Rust cdylib replaces
# KESolverLegacy.so (TurbulenceModels unit).
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `kesolverlegacy_` + `kesolverlegacy_init_`);
#  2. three distinct .sif differentials are byte-identical to the
#     reference except wall-clock lines — config A (steady, standard
#     k-epsilon, default legacy-bubble condensation, direct solve), config B
#     (transient BDF(1) x2 steps, RNG model, `Bubbles = False` stark
#     path with Default1stOrderTime, exercising the transient bx
#     history, the RNG secant loop, and the _Init legacy-bubble
#     timestep minimums), config C (transient like B but standard model
#     with default legacy bubbles — the transient condensed-bubble
#     CondensatePTransient path); all three carry a `SaveScalars` solver so
#     the solved Kinetic Energy / Kinetic Dissipation min/max land in
#     scalars.dat;
#  3. a second full process run of config B is byte-identical
#     (SAVE re-entry: AllocationsDone/bx/BubbleTimestep persist);
#  4. branch-liveness control: the `Result Norm` report differs across
#     configs A/B/C in both arms.
# The nondeterministic `Level=5` `RealTime` reports are structurally
# suppressed in both arms via `Max Output Level = 4` (the deterministic
# `Level=4` lines still print and are compared).
#
# Needs a configured+built serial Elmer tree (ELMER_BUILD — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-kesolverlegacy"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A KESolverLegacy proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: proof pinned to gfortran 13.3.0, re-measure on mismatch"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: unexpected rustc $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libkesolverlegacy.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/KESolverLegacy.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the KESolverLegacy module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T kesolverlegacy_" "$WORK/port.syms" \
  || { echo "FAIL: kesolverlegacy_ not the exported entry"; exit 1; }
grep -qx "T kesolverlegacy_init_" "$WORK/port.syms" \
  || { echo "FAIL: kesolverlegacy_init_ not exported"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 2
  Output Intervals = 1
  Max Output Level = Integer 4
End
Body 1
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 2
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 3
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 4
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Initial Condition 1
  Velocity 1 = Variable Coordinate 2
    Real MATC "0.1+0.1*tx"
  Velocity 2 = Real 0.0
  Kinetic Energy = Real 0.003
  Kinetic Dissipation = Real 0.0001
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = Real 1.0
  Viscosity = Real 1.0e-3
  KE SigmaK = Real 1.0
  KE SigmaE = Real 1.3
  KE C1 = Real 1.44
  KE C2 = Real 1.92
  KE Cmu = Real 0.09
  KE Clip = Real 1.0e-6
End
Boundary Condition 1
  Target Boundaries = 1
  Epsilon Wall BC = Logical True
End
Solver 1
  Equation = AllocateStuff
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = Velocity
  Variable DOFs = 2
End
Solver 2
  Equation = KEpsEq
  Procedure = "KESolverLegacy" "KESolverLegacy"
  Nonlinear System Max Iterations = Integer 1
  Nonlinear System Convergence Tolerance = Real 1.0e-5
  Linear System Solver = Direct
  Linear System Direct Method = Umfpack
  Steady State Convergence Tolerance = Real 1.0e-5
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Kinetic Energy
  Operator 1 = min
  Variable 2 = Kinetic Energy
  Operator 2 = max
  Variable 3 = Kinetic Dissipation
  Operator 3 = min
  Variable 4 = Kinetic Dissipation
  Operator 4 = max
End
EOF
sed -e 's/Simulation Type = Steady State/Simulation Type = Transient/' \
    -e 's/Steady State Max Iterations = 2/Timestepping Method = BDF\n  BDF Order = Integer 1\n  Timestep Sizes(1) = Real 0.1\n  Timestep Intervals(1) = Integer 2/' \
    -e 's/Procedure = "KESolverLegacy" "KESolverLegacy"/Procedure = "KESolverLegacy" "KESolverLegacy"\n  KE Model = String rng\n  Bubbles = Logical False/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -c "KE Model" "$WORK/caseB.sif"
grep -c "Bubbles = Logical False" "$WORK/caseB.sif"
grep -c "Timestep Intervals" "$WORK/caseB.sif"
# Config C: transient like B, but standard model with default (legacy)
# bubbles — exercises the transient condensed-bubble path
# (CondensatePTransient + bx history) that B skips.
grep -v -e "KE Model = String rng" -e "Bubbles = Logical False" \
    "$WORK/caseB.sif" > "$WORK/caseC.sif"
grep -c "Timestep Intervals" "$WORK/caseC.sif"
! grep -q "KE Model" "$WORK/caseC.sif" || echo "   config C drops the RNG branch"
! grep -q "Bubbles" "$WORK/caseC.sif" || echo "   config C uses default legacy bubbles"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/KESolverLegacy.so"
cp "$PORT_SO" "$WORK/shadow-port/KESolverLegacy.so"
nm -D --defined-only "$WORK/shadow-port/KESolverLegacy.so" | grep -q " T kesolverlegacy_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" || true
}

echo "-- differentials A/B/C + rerun of B"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseB port runB-port2
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
  norm "$WORK/runB-port2/$f" > "$WORK/b3.txt"
  diff "$WORK/b2.txt" "$WORK/b3.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   configs A+B+C byte-identical (stdout/stderr/scalars) + rerun byte-identical"

echo "-- branch-liveness control (Result Norm differs across A/B/C, both arms)"
for arm in orig port; do
  na=$(grep "Result Norm" "$WORK/runA-$arm/stdout.txt" | head -1)
  nb=$(grep "Result Norm" "$WORK/runB-$arm/stdout.txt" | head -1)
  nc=$(grep "Result Norm" "$WORK/runC-$arm/stdout.txt" | head -1)
  [ -n "$na" ] || { echo "FAIL: no Result Norm in A($arm)"; exit 1; }
  [ -n "$nb" ] || { echo "FAIL: no Result Norm in B($arm)"; exit 1; }
  [ -n "$nc" ] || { echo "FAIL: no Result Norm in C($arm)"; exit 1; }
  [ "$na" != "$nb" ] || { echo "FAIL: A($arm) and B($arm) norms identical"; exit 1; }
  [ "$na" != "$nc" ] || { echo "FAIL: A($arm) and C($arm) norms identical"; exit 1; }
  [ "$nb" != "$nc" ] || { echo "FAIL: B($arm) and C($arm) norms identical"; exit 1; }
done
echo "   control holds: branches change the numerics in both arms"

echo "KESOLVERLEGACY PROOF GREEN wall=$((SECONDS - T0))s"
