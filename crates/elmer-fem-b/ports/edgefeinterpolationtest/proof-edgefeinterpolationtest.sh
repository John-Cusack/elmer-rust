#!/bin/bash
# T022-A EdgeFEInterpolationTest solver-entry proof: the Rust cdylib
# replaces EdgeFEInterpolationTest.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `bestapproximationsolver_` and
#     `bestapproximationsolver_init0_`);
#  2. two distinct .sif differentials (config A: quadratic edge elements,
#     Test Mode 5, L2 error; config B: linear Piola elements, Test Mode 2
#     with Material Param, curl error) are byte-identical to the reference
#     except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state).
# The campaign invocation form `Procedure = "EdgeFEInterpolationTest"
# "BestApproximationSolver"` (library-vs-entry names differ) is used in
# both arms. The mesh is the interpolation-test square.grd from the Elmer
# checkout: edge-element interpolation staging needs that mesh.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="/home/john/firstmate/data/rs-elmer-t022a/scratch/proof-edgefeinterpolationtest"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A EdgeFEInterpolationTest proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libedgefeinterpolationtest.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/EdgeFEInterpolationTest.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T bestapproximationsolver_" "$WORK/port.syms" \
  || { echo "FAIL: bestapproximationsolver_ not an exported entry"; exit 1; }
grep -qx "T bestapproximationsolver_init0_" "$WORK/port.syms" \
  || { echo "FAIL: bestapproximationsolver_init0_ not an exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's interpolation-test square.grd)"
cp "$ELMER_FORK/fem/tests/SecondOrderEdgeElement2D_BCs/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Check Keywords "Warn"

Header
  Mesh DB "." "square"
End

Simulation
  Max Output Level = 4
  Coordinate System = "Cartesian 2D"
  Simulation Type = Steady
  Steady State Max Iterations = 1
  Output Intervals(1) = 1
End

Body 1
  Equation = 1
  Material = 1
End

Material 1
  Material Param = Real 0.0
End

Equation 1
  Active Solvers(1) = 1
End

Solver 1
  Equation = "Interpolation test"
  Variable = "P"
  Variable DOFs = 1
  Procedure = "EdgeFEInterpolationTest" "BestApproximationSolver"
  Use Piola Transform = Logical True
  Quadratic Approximation = Logical True
  Second Kind Basis = Logical False
  Compute Curl Error = Logical False
  Error Computation = Logical True
  Test Mode = Integer 5
  Optimize Bandwidth = False
  Bubbles in Global System = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = bicgstab2
  Linear System Preconditioning = none
  Linear System Residual Output = 0
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1e-11
End

Boundary Condition 1
  Target Boundaries(1) = 3
  P {e} 1 = Variable Coordinate
    Real MATC "-tx(0)*tx(1)*tx(1)"
  P {e} 2 = Real 0.0
End

Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = P
  Operator 1 = min
  Variable 2 = P
  Operator 2 = max
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Check Keywords "Warn"

Header
  Mesh DB "." "square"
End

Simulation
  Max Output Level = 4
  Coordinate System = "Cartesian 2D"
  Simulation Type = Steady
  Steady State Max Iterations = 1
  Output Intervals(1) = 1
End

Body 1
  Equation = 1
  Material = 1
End

Material 1
  Material Param = Real 1.0
End

Equation 1
  Active Solvers(1) = 1
End

Solver 1
  Equation = "Interpolation test"
  Variable = "P"
  Variable DOFs = 1
  Procedure = "EdgeFEInterpolationTest" "BestApproximationSolver"
  Use Piola Transform = Logical True
  Second Kind Basis = Logical False
  Compute Curl Error = Logical True
  Error Computation = Logical True
  Test Mode = Integer 2
  Optimize Bandwidth = False
  Bubbles in Global System = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = bicgstab2
  Linear System Preconditioning = none
  Linear System Residual Output = 0
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1e-11
End

Boundary Condition 1
  Target Boundaries(1) = 3
  P {e} 1 = Variable Coordinate
    Real MATC "1-tx(1)"
  P {e} 2 = Variable Coordinate
    Real MATC "1+tx(0)"
End

Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = P
  Operator 1 = min
  Variable 2 = P
  Operator 2 = max
End
EOF

echo "-- shadow module dirs (port .so under its own campaign name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/EdgeFEInterpolationTest.so"
cp "$PORT_SO" "$WORK/shadow-port/EdgeFEInterpolationTest.so"
nm -D --defined-only "$WORK/shadow-port/EdgeFEInterpolationTest.so" | grep -q "T bestapproximationsolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (quadratic, Test Mode 5, L2 error) + rerun (re-entry)"
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

echo "-- differential B (linear Piola, Test Mode 2, curl error, Material Param)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A and B take different code paths)"
norm "$WORK/runA-orig/stdout.txt" > "$WORK/l1.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/l2.txt"
if diff -q "$WORK/l1.txt" "$WORK/l2.txt" > /dev/null 2>&1; then
  echo "FAIL: liveness configs produce identical output (no branch taken)"
  exit 1
fi
echo "   liveness OK: configs differ, both arms agree on each"

echo "EDGEFEINTERPOLATIONTEST PROOF GREEN (wall $((SECONDS - T0))s)"
