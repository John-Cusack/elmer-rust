#!/bin/bash
# T021-B ThermoElectric solver-entry proof: the Rust cdylib replaces
# ThermoElectricSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `thermoelectricsolver_` and
#     `thermoelectricsolver_init_`);
#  2. three distinct .sif differentials (Picard steady config A;
#     SUPG-stabilized + anisotropic-BC + Newton config B, 2 steady
#     iterations i.e. the repeated-invocation run; transient
#     condensed-bubble config C exercising the bx/bxprev SAVE roll) are
#     byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch-liveness control: config B takes different code branches
#     (Stabilize, anisotropic BCs, Newton linearization) and produces
#     different output, with both arms agreeing on it.
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
WORK="$LANE/scratch/proof-thermoelectricsolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B ThermoElectric proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0, got: $GFORTRAN_VERSION"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: need nightly feaadeeac 2026-09-19, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libthermoelectricsolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release 2>&1 | tail -3)
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release -- -D warnings) \
  || { echo "FAIL: clippy warnings"; exit 1; }
echo "   build + clippy clean"

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ThermoElectricSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ThermoElectricSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T thermoelectricsolver_" "$WORK/port.syms" \
  || { echo "FAIL: thermoelectricsolver_ not the exported entry"; exit 1; }
grep -qx "T thermoelectricsolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: thermoelectricsolver_init_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the ThermoElectric test's square.grd)"
cp "$ELMER_FORK/fem/tests/ThermoElectric/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Check Keywords "Warn"
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady
  Steady State Max Iterations = 2
  Output Intervals(1) = 1
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(1) = 1
End
Body Force 1
  Current Source = Real 100
End
Material 1
  Electric Conductivity = 1
  Density = 1
  Heat Capacity = 1
  Heat Conductivity = 1
  Seebeck Coefficient = Real -0.005
End
Solver 1
  Equation = "ThermoElectric"
  Variable = POT[Temperature:1 Potential:1]
  Procedure = "ThermoElectricSolver" "ThermoElectricSolver"
  Element = "p:1"
  Nonlinear System Convergence Tolerance = 1e-6
  Nonlinear System Max Iterations = 5
  Nonlinear System Newton After Iterations = 50
  Nonlinear System Newton After Tolerance = 1e-9
  Linear System Solver = "Iterative"
  Linear System Iterative Method = BicgstabL
  Bicgstabl Polynomial Degree = 2
  Linear System Max Iterations = 200
  Linear System Residual Output = 40
  Linear System Preconditioning = Ilu
  Linear System Convergence Tolerance = 1e-8
  Steady State Convergence Tolerance = 1e-6
End
Boundary Condition 1
  Target Boundaries = 1
  Potential = 0
  Temperature = 300
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = min
  Variable 2 = Temperature
  Operator 2 = max
  Variable 3 = Potential
  Operator 3 = min
  Variable 4 = Potential
  Operator 4 = max
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Check Keywords "Warn"
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady
  Steady State Max Iterations = 2
  Output Intervals(1) = 1
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(1) = 1
End
Body Force 1
  Current Source = Real 100
End
Material 1
  Electric Conductivity = 1
  Density = 1
  Heat Capacity = 1
  Heat Conductivity = 1
  Seebeck Coefficient = Real -0.005
End
Solver 1
  Equation = "ThermoElectric"
  Variable = POT[Temperature:1 Potential:1]
  Procedure = "ThermoElectricSolver" "ThermoElectricSolver"
  Stabilize = True
  Nonlinear System Convergence Tolerance = 1e-6
  Nonlinear System Max Iterations = 5
  Nonlinear System Newton After Iterations = 1
  Nonlinear System Newton After Tolerance = 1e-9
  Linear System Solver = "Iterative"
  Linear System Iterative Method = BicgstabL
  Bicgstabl Polynomial Degree = 2
  Linear System Max Iterations = 200
  Linear System Residual Output = 40
  Linear System Preconditioning = Ilu
  Linear System Convergence Tolerance = 1e-8
  Steady State Convergence Tolerance = 1e-6
End
Boundary Condition 1
  Target Boundaries = 1
  Potential = 0
  Temperature = 300
End
Boundary Condition 2
  Target Boundaries(2) = 2 4
  Electric Transfer Coefficient(2,2) = Real 2.0 0.5 0.5 1.0
  External Potential = Real 5.0
End
Boundary Condition 3
  Target Boundaries = 3
  Electric Flux Vector 1 = Real 30.0
  Electric Flux Vector 2 = Real 50.0
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = min
  Variable 2 = Temperature
  Operator 2 = max
  Variable 3 = Potential
  Operator 3 = min
  Variable 4 = Potential
  Operator 4 = max
End
EOF
cat > "$WORK/caseC.sif" <<'EOF'
Check Keywords "Warn"
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Transient
  Timestepping Method = BDF
  BDF Order = 1
  Timestep Sizes(1) = 0.1
  Timestep Intervals(1) = 2
  Output Intervals(1) = 1
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(1) = 1
End
Body Force 1
  Current Source = Real 100
End
Material 1
  Electric Conductivity = 1
  Density = 1
  Heat Capacity = 1
  Heat Conductivity = 1
  Seebeck Coefficient = Real -0.005
End
Solver 1
  Equation = "ThermoElectric"
  Variable = POT[Temperature:1 Potential:1]
  Procedure = "ThermoElectricSolver" "ThermoElectricSolver"
  Element = "p:1 b:4"
  Nonlinear System Convergence Tolerance = 1e-6
  Nonlinear System Max Iterations = 5
  Nonlinear System Newton After Iterations = 1
  Nonlinear System Newton After Tolerance = 1e-9
  Linear System Solver = "Iterative"
  Linear System Iterative Method = BicgstabL
  Bicgstabl Polynomial Degree = 2
  Linear System Max Iterations = 200
  Linear System Preconditioning = Ilu
  Linear System Convergence Tolerance = 1e-8
End
Boundary Condition 1
  Target Boundaries = 1
  Potential = 0
  Temperature = 300
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = min
  Variable 2 = Temperature
  Operator 2 = max
  Variable 3 = Potential
  Operator 3 = min
  Variable 4 = Potential
  Operator 4 = max
End
EOF

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ThermoElectricSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/ThermoElectricSolver.so"
nm -D --defined-only "$WORK/shadow-port/ThermoElectricSolver.so" | grep -q "T thermoelectricsolver_" \
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
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; tail -20 "$WORK/$run/stderr.txt"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER SOLVER.*STARTED AT|ELMER.*FINISHED AT|STARTED AT|Factorization time" "$1" || true
}

echo "-- differential A (Picard steady) + rerun (SAVE re-entry)"
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

echo "-- differential B (stabilized + anisotropic BC + Newton, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- differential C (transient condensed bubbles, SAVE roll)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B takes different branches than A)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
if diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1; then
  echo "FAIL: liveness control configs A and B agree (expected branch divergence)"
  exit 1
fi
echo "   A and B differ as expected (both arms agree per differentials above)"

echo "THERMOELECTRICSOLVER PROOF GREEN ($((SECONDS - T0))s wall)"
