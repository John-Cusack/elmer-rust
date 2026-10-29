#!/bin/bash
# T021-B ShearrateSolver solver-entry proof: the Rust cdylib replaces
# ShearrateSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `shearratesolver_` + `shearratesolver_init_`);
#  2. two distinct .sif differentials (lid-driven cavity via FlowSolve;
#     config A: plain shearrate; config B: `Calculate Viscosity = True`
#     branch with extra viscosity solve + exported `Viscosity velocity`
#     variable) are byte-identical to the reference except wall-clock
#     lines (both configs run 2 steady iterations, i.e. the
#     repeated-invocation SAVE re-entry run is inside every arm);
#  3. a second full process run of config A against the port is
#     byte-identical (idempotent re-entry through the SAVE state);
#  4. branch liveness: config B output differs from config A while both
#     arms agree on each config.
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
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-shearratesolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B ShearrateSolver proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libshearratesolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ShearrateSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ShearrateSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T shearratesolver_" "$WORK/port.syms" \
  || { echo "FAIL: shearratesolver_ not an exported entry"; exit 1; }
grep -qx "T shearratesolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: shearratesolver_init_ not an exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's PoissonDB square.grd)"
echo "   (triangular unit-square mesh: 4 bodies, so all 4 are declared in"
echo "   the .sif files below; Elmer core segfaults in initcond otherwise)"
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
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
  Viscosity = 1.0e-3
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Velocity 1 = 0.0
  Velocity 2 = 0.0
End
Boundary Condition 2
  Target Boundaries(2) = 2 4
  Velocity 1 = 0.0
  Velocity 2 = 0.0
End
Boundary Condition 3
  Target Boundaries(1) = 3
  Velocity 1 = 1.0
  Velocity 2 = 0.0
End
Solver 1
  Equation = Navier-Stokes
  Procedure = "FlowSolve" "FlowSolver"
  Variable = Flow Solution[Velocity:2 Pressure:1]
  Exec Solver = Always
  Stabilize = True
  Nonlinear System Convergence Tolerance = 1.0e-8
  Nonlinear System Max Iterations = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStabL
  BiCgStabl Polynomial Degree = 4
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
  Linear System Abort Not Converged = False
End
Solver 2
  Equation = Shearrate
  Procedure = "ShearrateSolver" "ShearrateSolver"
  Target Variable = String velocity
  Exec Solver = Always
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Shearrate velocity
  Operator 1 = min
  Variable 2 = Shearrate velocity
  Operator 2 = max
End
EOF
python3 - "$WORK/caseA.sif" "$WORK/caseB.sif" <<'PYEOF'
import sys
a = open(sys.argv[1]).read()
b = a.replace("""  Target Variable = String velocity
  Exec Solver = Always""", """  Target Variable = String velocity
  Calculate Viscosity = Logical True
  Exec Solver = Always""")
b = b.replace("""  Variable 2 = Shearrate velocity
  Operator 2 = max""", """  Variable 2 = Shearrate velocity
  Operator 2 = max
  Variable 3 = Viscosity velocity
  Operator 3 = max
  Variable 4 = Viscosity velocity
  Operator 4 = min""")
open(sys.argv[2], "w").write(b)
PYEOF
grep -c "Calculate Viscosity" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ShearrateSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/ShearrateSolver.so"
nm -D --defined-only "$WORK/shadow-port/ShearrateSolver.so" | grep -q "T shearratesolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/ShearrateSolver.so" | grep -q "T shearratesolver_init_" \
  || { echo "FAIL: staged port .so lost its init entry"; exit 1; }

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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|SOLVER FINISHED AT|ELMER FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (plain shearrate) + rerun (re-entry)"
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

echo "-- differential B (Calculate Viscosity branch)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B must differ from A, both arms agree)"
norm "$WORK/runA-orig/scalars.dat" > "$WORK/la.txt"
norm "$WORK/runB-orig/scalars.dat" > "$WORK/lb.txt"
diff "$WORK/la.txt" "$WORK/lb.txt" > /dev/null 2>&1 && \
  { echo "FAIL: branch-liveness control: A and B scalars identical"; exit 1; }
echo "   A vs B scalars differ, both arms agree"

echo "SHEARRATESOLVER PROOF GREEN"
echo "elapsed: $((SECONDS - T0))s"
