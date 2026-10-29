#!/bin/bash
# T022-A ElectricForce solver-entry proof: the Rust cdylib replaces
# ElectricForce.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `statelecforce_`, invoked via the campaign
#     `Procedure = "ElectricForce" "StatElecForce"` form);
#  2. two distinct .sif differentials (force-on-one-boundary config A;
#     force-on-two-boundaries + nodal + moment config B with 2 steady
#     iterations, i.e. the repeated-invocation run) are byte-identical to
#     the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state).
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="${RS_ELMER_T022A_WORK:-/home/john/firstmate/data/rs-elmer-t022a/scratch/proof-electricforce}"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A ElectricForce proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libelectricforce.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ElectricForce.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ElectricForce module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T statelecforce_" "$WORK/port.syms" \
  || { echo "FAIL: statelecforce_ not the exported entry"; exit 1; }
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
  Steady State Max Iterations = 1
  Output Intervals = 1
  Max Output Level = 4
End
Constants
  Permittivity Of Vacuum = 8.8542e-12
End
Body 1
  Equation = 1
  Material = 1
End
Body 2
  Equation = 1
  Material = 1
End
Body 3
  Equation = 1
  Material = 1
End
Body 4
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Relative Permittivity = 2.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Potential = Real 0.0
  Calculate Electric Force = Logical True
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Potential = Real 1.0
End
Solver 1
  Equation = Stat Elec Solver
  Procedure = "StatElecSolve" "StatElecSolver"
  Variable = Potential
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-10
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = Electric Force
  Procedure = "ElectricForce" "StatElecForce"
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
End
EOF
# Config B: force on two boundaries + nodal/field variables + moment
# reporting, 2 steady iterations (in-process re-entry through the port).
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    -e 's/Target Boundaries(3) = 2 3 4/Target Boundaries(1) = 2/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
t = t.replace("""Boundary Condition 2
  Target Boundaries(1) = 2
  Potential = Real 1.0
End""", """Boundary Condition 2
  Target Boundaries(1) = 2
  Potential = Real 1.0
  Calculate Electric Force = Logical True
End
Boundary Condition 3
  Target Boundaries(2) = 3 4
  Potential = Real 1.0
End""", 1)
t = t.replace("""  Procedure = "ElectricForce" "StatElecForce"
End""", """  Procedure = "ElectricForce" "StatElecForce"
  Calculate Nodal Force = Logical True
  Calculate Force Density = Logical True
End""", 1)
t = t.replace("""  Max Output Level = 4
End""", """  Max Output Level = 4
  Moment About 1 = Real 0.5
  Moment About 2 = Real 0.5
End""", 1)
open(p, "w").write(t)
PYEOF
grep -c "Calculate Electric Force" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ElectricForce.so"
cp "$PORT_SO" "$WORK/shadow-port/ElectricForce.so"
nm -D --defined-only "$WORK/shadow-port/ElectricForce.so" | grep -q "T statelecforce_" \
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
norm() { # drop wall-clock-only lines for the byte-compare (merged class)
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" || true
}

echo "-- differential A (single force boundary) + rerun (re-entry)"
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
norm "$WORK/runA-orig/scalars.dat.names" > "$WORK/a0.txt"
norm "$WORK/runA-port/scalars.dat.names" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A scalars.dat.names diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat scalars.dat.names; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (two force boundaries + nodal + moment, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat scalars.dat.names; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness control (B takes a different branch than A)"
if diff -q "$WORK/runA-port/stdout.txt" "$WORK/runB-port/stdout.txt" > /dev/null 2>&1; then
  echo "FAIL: liveness control produced identical output"
  exit 1
fi
grep -h "Net electric force" "$WORK/runA-port/stdout.txt" | head -2
grep -h "Net electric force" "$WORK/runB-port/stdout.txt" | head -4
echo "   liveness control differs as expected"

echo "ELECTRICFORCE PROOF GREEN"
echo "wall seconds: $((SECONDS - T0))"
