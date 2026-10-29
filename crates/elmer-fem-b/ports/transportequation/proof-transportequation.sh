#!/bin/bash
# T021-B TransportEquation solver-entry proof: the Rust cdylib replaces
# TransportEquation.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `transportequationsolver_`);
#  2. two distinct .sif differentials (transient constant-advection config
#     A; computed-advection + FlowSolver config B, exercising the other
#     `Advection` branch) are byte-identical to the reference except
#     wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state; the 2-timestep transient
#     run itself already re-enters the solver, hitting the even-call
#     no-op gate on the second step).
#
# Defined-path notes: both configs set `Time Derivative Order = 2` on the
# transport solver (Elmer core then sizes `PrevValues` to 7 columns, so the
# `PrevValues(:,3/4/5)` X/V/A reads are in bounds — without it the
# reference itself reads out of bounds). The ` Assembly: (s)` /
# ` Solve:    (s)` Info lines print process CPU times and are normed out
# as wall-clock lines.
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
WORK="$LANE/scratch/proof-transportequation"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B TransportEquation proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libtransportequation.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/TransportEquation.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the TransportEquation module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T transportequationsolver_" "$WORK/port.syms" \
  || { echo "FAIL: transportequationsolver_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's PoissonDB square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Transient
  Timestep Sizes = 0.05
  Timestep Intervals = 2
  Output Intervals = 1
  Steady State Max Iterations = 1
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
  Active Solvers(4) = 1 2 3 4
  Heat Equation = Logical True
  RateOfChange = Logical True
  Transport = Logical True
End
Material 1
  Density = 1.0
  Heat Conductivity = 1.0
  Advection Velocity 1 = Real 0.5
  Advection Velocity 2 = Real 0.25
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Temperature = Real 0.0
  TempRate = Real 0.0
  Transp = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Temperature = Real 1.0
  Transp = Real 1.0
End
Solver 1
  Equation = Heat Equation
  Procedure = "HeatSolve" "HeatSolver"
  Variable = Temperature
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = RateOfChange
  Procedure = "RateOfChange" "RateOfChangeSolver"
  Variable = TempRate
  Variable DOFs = 1
  Advection Variable = String Temperature
  Advection = String constant
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 3
  Equation = Transport
  Procedure = "TransportEquation" "TransportEquationSolver"
  Variable = Transp
  Variable DOFs = 1
  Time Derivative Order = 2
  Advection Variable = String Temperature
  Rate Of Change Equation Variable = String TempRate
  Advection = String constant
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 4
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Transp
  Operator 1 = min
  Variable 2 = Transp
  Operator 2 = max
End
EOF
# Config B: computed-advection branch — FlowSolver provides the
# Flow Solution velocity field, both RateOfChange and TransportEquation
# take the 'computed' path instead of the constant Material velocity.
sed -e 's/Active Solvers(4) = 1 2 3 4/Active Solvers(5) = 1 2 3 4 5/' \
    -e 's/RateOfChange = Logical True/RateOfChange = Logical True\n  Navier-Stokes = Logical True/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
flow = """Solver 2
  Equation = Navier-Stokes
  Procedure = "FlowSolve" "FlowSolver"
  Variable = Flow Solution[velocity:2 pressure:1]
  Stabilization Method = Stabilized
  Steady State Convergence Tolerance = 1.0e-5
  Nonlinear System Convergence Tolerance = 1.0e-8
  Nonlinear System Max Iterations = 1
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Relaxation Factor = 1
  Nonlinear System Newton After Tolerance = 1.0e-3
  Linear System Solver = Iterative
  Linear System Symmetric = True
  Linear System Iterative Method = BicgstabL
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-10
  Linear System Preconditioning = ILU0
End
"""
t = t.replace("Solver 2\n  Equation = RateOfChange", flow + "Solver 3\n  Equation = RateOfChange", 1)
t = t.replace("Solver 3\n  Equation = Transport", "Solver 4\n  Equation = Transport", 1)
t = t.replace("Solver 4\n  Equation = \"SaveScalars\"", "Solver 5\n  Equation = \"SaveScalars\"", 1)
t = t.replace("Advection = String constant", "Advection = String computed")
t = t.replace("""Material 1
  Density = 1.0
  Heat Conductivity = 1.0
  Advection Velocity 1 = Real 0.5
  Advection Velocity 2 = Real 0.25
End""", """Material 1
  Density = 1.0
  Heat Conductivity = 1.0
  Viscosity = 1.0
  Advection Velocity 1 = Real 0.5
  Advection Velocity 2 = Real 0.25
End""")
t = t.replace("""Boundary Condition 1
  Target Boundaries(1) = 1
  Temperature = Real 0.0
  TempRate = Real 0.0
  Transp = Real 0.0
End""", """Boundary Condition 1
  Target Boundaries(1) = 1
  Temperature = Real 0.0
  TempRate = Real 0.0
  Transp = Real 0.0
  Velocity 1 = Real 0.0
  Velocity 2 = Real 0.0
  Pressure = Real 0.0
End""")
t = t.replace("""Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Temperature = Real 1.0
  Transp = Real 1.0
End""", """Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Temperature = Real 2.0
  Transp = Real 1.0
  Velocity 1 = Real 0.1
  Velocity 2 = Real 0.0
End""")
open(p, "w").write(t)
PYEOF
grep -q "FlowSolve" "$WORK/caseB.sif"
grep -q "Advection = String computed" "$WORK/caseB.sif"
! grep -q "Advection = String constant" "$WORK/caseB.sif" || { echo "FAIL: caseB still has constant advection"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/TransportEquation.so"
cp "$PORT_SO" "$WORK/shadow-port/TransportEquation.so"
nm -D --defined-only "$WORK/shadow-port/TransportEquation.so" | grep -q "T transportequationsolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER SOLVER.*STARTED AT|ELMER.*FINISHED AT|STARTED AT|Assembly: \(s\)|Solve: *\(s\)" "$1" || true
}

echo "-- differential A (constant advection) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
done
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (computed advection + flow, 2 steps)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A and B outputs must differ)"
if diff -q "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: A/B scalars identical — branch control dead"; exit 1
fi
echo "   A/B outputs differ and both arms agree"

echo "TRANSPORTEQUATION PROOF GREEN (wall: $((SECONDS - T0))s)"
