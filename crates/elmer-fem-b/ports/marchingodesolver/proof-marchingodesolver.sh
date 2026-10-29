#!/bin/bash
# T021-B MarchingODESolver solver-entry proof: the Rust cdylib replaces
# MarchingODESolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `marchingodesolver_` and
#     `marchingodesolver_init_`);
#  2. two distinct .sif differentials (config A: steady Newmark-Beta
#     march with 2 steady iterations, i.e. the repeated-invocation run;
#     config B: transient march, exercising the DoTransient/dtn passage
#     branch) are byte-identical to the reference except wall-clock
#     lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state).
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout for the MarchingODE1 rect.grd mesh input
# (structured extruded grid via the same ElmerGrid flow; the PoissonDB
# triangular mesh is not extruded and cannot march).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="/home/john/firstmate/data/rs-elmer-t021b/scratch/proof-marchingodesolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/rect*
echo "== MarchingODESolver proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libmarchingodesolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/MarchingODESolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the MarchingODESolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T marchingodesolver_" "$WORK/port.syms" \
  || { echo "FAIL: marchingodesolver_ not the exported entry"; exit 1; }
grep -qx "T marchingodesolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: marchingodesolver_init_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's MarchingODE1 rect.grd)"
cp "$ELMER_FORK/fem/tests/MarchingODE1/rect.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 rect > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Check Keywords "Warn"

Header
  Mesh DB "." "rect"
End

Simulation
  Max Output Level = 5
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Output Intervals(1) = 1
  Steady State Max Iterations = 2
End

Constants
  Stefan Boltzmann = 5.67e-8
End

Body 1
  Equation = 1
  Material = 1
  Initial Condition = 1
End

Equation 1
  Active Solvers(1) = 1
End

Initial Condition 1
  Name = "Guess"
  Integrator = Real 0.0
End

Material 1
  Name = "Material"
  Integrator: Source = Equals "Coordinate 1"
End

Solver 1
  Exec Solver = always
  Equation = "MarchingOde"
  Procedure = "MarchingODESolver" "MarchingODESolver"

! Detects structured mesh in x-direction
  Active Coordinate = Integer 1

! This gives 2nd order accuracy
  Newmark Beta = Real 0.5

  Draw Velocity = Real 1.0

! The name of the target field
  Variable = Integrator
End

Boundary Condition 1
  Name = "Sides"
  Target Boundaries(4) = 1 2 3 4
End

Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Integrator
  Operator 1 = min
  Variable 2 = Integrator
  Operator 2 = max
End
EOF
# Config B: transient march (one mesh layer per timestep). Exercises the
# DoTransient/dtn passage branch that config A never takes.
sed -e 's/Simulation Type = Steady State/Simulation Type = Transient/' \
    -e 's/Steady State Max Iterations = 2/Timestep Intervals = 7\n  Timestep Sizes = 0.1\n  Steady State Max Iterations = 1/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -q "Timestep Intervals = 7" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/MarchingODESolver.so"
cp "$PORT_SO" "$WORK/shadow-port/MarchingODESolver.so"
nm -D --defined-only "$WORK/shadow-port/MarchingODESolver.so" | grep -q "T marchingodesolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/rect" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (steady Newmark-Beta march, 2 steady iterations) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
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

echo "-- differential B (transient march, DoTransient branch)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A and B take different branches)"
if diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1; then
  echo "FAIL: branch control produced identical output"; exit 1;
fi
echo "   A and B outputs differ, both arms agree on each"

echo "MARCHINGODESOLVER PROOF GREEN"
echo "wall seconds: $((SECONDS - T0))"
