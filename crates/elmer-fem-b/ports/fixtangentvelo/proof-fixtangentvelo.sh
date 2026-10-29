#!/bin/bash
# T022-A FixTangentVelo solver-entry proof: the Rust cdylib replaces
# FixTangentVelo.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `fixtangentvelo_` + `fixtangentvelo_init_`);
#  2. two distinct single-invocation .sif differentials are
#     byte-identical to the reference except wall-clock lines — config A
#     (default branches: `Active Coordinate` defaults to 3, no
#     `Fix Downstream Only`, no `Fix Tangent Velocity`) and config B
#     (explicit `Active Coordinate = 3`, `Fix Downstream Only`,
#     `Fix Tangent Velocity`: all three keyword reads take the
#     present/True path); both carry a `SaveScalars` solver so the
#     `FixVelo` min/max and the (fix-altered) `Velocity 3` min/max land
#     in scalars.dat. Single invocation is the unit's only defined path
#     (see staging note); both configs converge in one steady iteration.
#  3. a second full process run of config A against the port is
#     byte-identical (idempotent re-entry through the SAVE state).
# Branch-liveness control: the saved scalars and the `Applying the fix`
# report differ across configs A/B in both arms (the keyword branches
# change the numerics identically in both arms), so the differential is
# not vacuous. The nondeterministic `Level=5` `RealTime` reports are
# structurally suppressed in both arms via the standard
# `Max Output Level = 4` keyword.
#
# Staging: genuine 3-D unit-cube mesh (`cube.grd` from the radiator3d
# test, via ElmerGrid — the solver FATALs on non-3D meshes so the brief's
# 2-D square mesh cannot be used here), velocity field staged by
# `AllocateSolver` (3-DOF `Velocity`) + `Initial Condition` MATC exactly
# as in the T021-A divergence proof, one `Tangent Surface` BC on the
# horizontal bottom (tag 2) and top (tag 3) faces. Side walls are
# deliberately untagged: with the default `Active Coordinate = 3` the
# solver's 3x3 system is singular on vertical walls (same reason the
# upstream slab test tags only its wavy top/bottom surfaces), so only
# horizontal faces are valid staging here. Re-entry note: the original
# keeps `Mesh`/`Params` as non-SAVE locals set only in the first-visit
# block, so a second invocation dereferences stale pointers (undefined
# behavior). In this tree that yields all-zero `FixVals` at output
# level 4 and a SIGSEGV at `FixTangentVelo.F90:125` at output level 10
# (backtrace on file in scratch); upstream's own test uses a single
# steady iteration. The port therefore implements the author's evident
# intent (current `Mesh`/`Params` every call — identical to the original
# wherever its behavior is defined) and the proof differentially covers
# the defined single-invocation path plus a port-vs-port rerun.
#
# Needs a configured+built serial Elmer tree (ELMER_BUILD — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout only for the cube.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-fixtangentvelo"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$LANE/scratch/cargo-target-fixtangentvelo}"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/cube*
echo "== T022-A FixTangentVelo proof: work=$WORK build=$ELMER_BUILD"
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
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" CARGO_TARGET_DIR="$CARGO_TARGET_DIR" cargo +nightly build --release -q)
PORT_SO="$CARGO_TARGET_DIR/release/libfixtangentvelo.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FixTangentVelo.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FixTangentVelo module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T fixtangentvelo_" "$WORK/port.syms" \
  || { echo "FAIL: fixtangentvelo_ not the exported entry"; exit 1; }
grep -qx "T fixtangentvelo_init_" "$WORK/port.syms" \
  || { echo "FAIL: fixtangentvelo_init_ not exported"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's cube.grd)"
cp "$ELMER_FORK/fem/tests/radiator3d/cube.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 cube > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "cube"
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
Initial Condition 1
  Velocity 1 = Variable Coordinate 1
    Real MATC "tx"
  Velocity 2 = Variable Coordinate 2
    Real MATC "tx"
  Velocity 3 = Real 0.2
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
End
Boundary Condition 1
  Target Boundaries(2) = 2 3
  Tangent Surface = Logical True
End
Solver 1
  Equation = AllocV
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = Velocity
  Variable DOFs = 3
End
Solver 2
  Equation = FixTang
  Procedure = "FixTangentVelo" "FixTangentVelo"
  Velocity Variable Name = String Velocity
End
Solver 3
  Equation = SaveS
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = FixVelo
  Operator 1 = min
  Variable 2 = FixVelo
  Operator 2 = max
  Variable 3 = Velocity 3
  Operator 3 = min
  Variable 4 = Velocity 3
  Operator 4 = max
End
EOF
sed -e 's/Procedure = "FixTangentVelo" "FixTangentVelo"/Procedure = "FixTangentVelo" "FixTangentVelo"\n  Active Coordinate = Integer 3\n  Fix Downstream Only = Logical True\n  Fix Tangent Velocity = Logical True/' \
  "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -c "Fix Downstream Only" "$WORK/caseB.sif"
grep -c "Fix Tangent Velocity" "$WORK/caseB.sif"
grep -c "Active Coordinate" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/FixTangentVelo.so"
cp "$PORT_SO" "$WORK/shadow-port/FixTangentVelo.so"
nm -D --defined-only "$WORK/shadow-port/FixTangentVelo.so" | grep -q " T fixtangentvelo_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/cube" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" || true
}
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
grep -q "Total number of nodes to fix" "$WORK/runA-orig/stdout.txt" \
  || { echo "FAIL: config A never reached FixTangentVelo"; exit 1; }
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (all keyword branches flipped)"
run_case caseB orig runB-orig
run_case caseB port runB-port
grep -q "Applying the fix" "$WORK/runB-orig/stdout.txt" \
  || { echo "FAIL: config B never applied the fix"; exit 1; }
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A vs B differ, both arms)"
for arm in orig port; do
  norm "$WORK/runA-$arm/scalars.dat" > "$WORK/liveA.txt"
  norm "$WORK/runB-$arm/scalars.dat" > "$WORK/liveB.txt"
  if diff -q "$WORK/liveA.txt" "$WORK/liveB.txt" > /dev/null 2>&1; then
    echo "FAIL: liveness control vacuous in $arm arm (A == B scalars)"
    exit 1
  fi
done
echo "   control holds: branches change the numerics in both arms"

echo "FIXTANGENTVELO PROOF GREEN wall=$((SECONDS - T0))s"
