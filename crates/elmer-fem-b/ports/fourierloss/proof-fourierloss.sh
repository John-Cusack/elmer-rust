#!/bin/bash
# T022-A FourierLoss solver-entry proof: the Rust cdylib replaces FourierLoss.so.
#
# Reference unit: fem/src/modules/FourierLoss.F90 (two F77-global entries
# fourierlosssolver_ + fourierlosssolver_init_, three CONTAINS internals).
# Differential strategy: transient AdvDiff driver (ModelPDE) produces a
# time-varying nodal target field with PrevValues; the FourierLoss solver
# accumulates its on-the-fly DFT over [t0,t0+tcycle] and computes losses at
# the two end-cycle steps.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: fourierlosssolver_ + fourierlosssolver_init_);
#  2. two distinct .sif differentials (config A: Direct/exact/summed
#     components + loss-file outputs; config B: curl/inexact-Simpsons/
#     separate components i.e. the multi-Rhs OtherRhs path) are
#     byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch liveness: A and B outputs differ on both arms.
#
# SAVE-state note: the reference keeps entry/DFT state in SAVE variables
# that persist per process; the rerun bar item is satisfied with the
# port-vs-port fresh-process rerun the bar explicitly allows.
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
WORK="/home/john/firstmate/data/rs-elmer-t022a/scratch/proof-fourierloss"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A FourierLoss proof: work=$WORK build=$ELMER_BUILD"
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
  *) echo "FAIL: need rustc nightly 1.100.0"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-/home/john/firstmate/data/rs-elmer-t022a/scratch/cargo-target-fourierloss}"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libfourierloss.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FourierLoss.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FourierLoss module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T fourierlosssolver_" "$WORK/port.syms" \
  || { echo "FAIL: fourierlosssolver_ not the exported entry"; exit 1; }
grep -qx "T fourierlosssolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: fourierlosssolver_init_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- configs (clean-room .sif files below)"
python3 - "$WORK/caseA.sif" "$WORK/caseB.sif" <<'PYEOF'
import sys, math
pathA, pathB = sys.argv[1], sys.argv[2]

def pot_table():
    lines = []
    t = 0.0
    while t <= 3.0001:
        v = (2 + math.cos(2*math.pi*t) + 0.5*math.cos(4*math.pi*t)
             + 0.25*math.cos(6*math.pi*t) + 0.4*math.sin(2*math.pi*t)
             + 0.2*math.sin(4*math.pi*t) + 0.1*math.sin(6*math.pi*t))
        lines.append("      %.2f %.12e" % (t, v))
        t += 0.05
    return "\n".join(lines)

table = pot_table()

head = """Header
  Mesh DB "." "square"
End
Simulation
  Max Output Level = 4
  Coordinate System = Cartesian
  Simulation Type = Transient
  Timestepping Method = BDF
  BDF Order = 1
  Timestep Intervals(1) = 30
  Timestep Sizes(1) = 0.1
  Steady State Max Iterations = 1
  Frequency = Real 1.0
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 2
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 3
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 4
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Time Derivative Coefficient = Real 0.0
  Diffusion Coefficient = Real 1.0
  Harmonic Loss Coefficient 1 = Variable frequency
    Real
      0.0 1.0
      3.0 1.0
    End
  Harmonic Loss Coefficient 2 = Variable frequency
    Real
      0.0 2.0
      3.0 2.0
    End
End
Body Force 1
  Field Source = Real 0.0
End
Boundary Condition 1
  Target Boundaries(4) = 1 2 3 4
  Potential = Variable Time
    Real
__POT_TABLE__
    End
End
Solver 1
  Equation = AdvDiff
  Procedure = "ModelPDE" "AdvDiffSolver"
  Variable = String Potential
  Variable DOFs = Integer 1
  Linear System Solver = Iterative
  Linear System Preconditioning = ILU0
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Iterative Method = BiCGStab
End
"""

solverA = """Solver 2
  Equation = Fourier
  Procedure = "FourierLoss" "FourierLossSolver"
  Target Variable = String Potential
  Target Variable Direct = Logical True
  Fourier Series Components = Integer 2
  Fourier Series Output = Logical True
  Calculate Nodal Fields = Logical True
  Calculate Elemental Fields = Logical True
  Separate Loss Components = Logical False
  Calculate Nodal Losses = Logical True
  Harmonic Loss Frequency Exponent(2) = Real 0.0 0.0
  Harmonic Loss Field Exponent(2) = Real 1.0 1.0
  Fourier Start Cycles = Integer 1
  Fourier Integrate Cycles = Integer 1
  Fourier Loss Filename = File bodyloss.dat
  Series Loss Filename = File seriesloss.dat
  Frequency = Real 1.0
  Linear System Solver = Iterative
  Linear System Preconditioning = None
  Linear System Max Iterations = 5000
  Linear System Iterative Method = CG
  Linear System Convergence Tolerance = 1.0e-6
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Exec Solver = after timestep
  Filename = scalars.dat
  Variable 1 = String "Fourier Loss"
  Operator 1 = max
  Variable 2 = String "Fourier Loss"
  Operator 2 = min
  Variable 3 = String Potential
  Operator 3 = mean
End
"""

solverB = """Solver 2
  Equation = Fourier
  Procedure = "FourierLoss" "FourierLossSolver"
  Target Variable = String Potential
  Target Variable Direct = Logical False
  Inexact Integration = Logical True
  Simpsons Rule = Logical True
  Fourier Series Components = Integer 1
  Fourier Series Output = Logical True
  Calculate Nodal Fields = Logical True
  Calculate Elemental Fields = Logical True
  Separate Loss Components = Logical True
  Calculate Nodal Losses = Logical True
  Harmonic Loss Frequency Exponent(2) = Real 1.0 2.0
  Harmonic Loss Field Exponent(2) = Real 1.0 2.0
  Fourier Start Cycles = Integer 1
  Fourier Integrate Cycles = Integer 1
  Frequency = Real 1.0
  Linear System Solver = Iterative
  Linear System Preconditioning = None
  Linear System Max Iterations = 5000
  Linear System Iterative Method = CG
  Linear System Convergence Tolerance = 1.0e-6
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Exec Solver = after timestep
  Filename = scalars.dat
  Variable 1 = String "Fourier Loss 1"
  Operator 1 = max
  Variable 2 = String "Fourier Loss 2"
  Operator 2 = max
  Variable 3 = String Potential
  Operator 3 = mean
End
"""

open(pathA, "w").write((head + solverA).replace("__POT_TABLE__", table))
open(pathB, "w").write((head + solverB).replace("__POT_TABLE__", table))
PYEOF
grep -c "Harmonic Loss Coefficient" "$WORK/caseA.sif" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/FourierLoss.so"
cp "$PORT_SO" "$WORK/shadow-port/FourierLoss.so"
nm -D --defined-only "$WORK/shadow-port/FourierLoss.so" | grep -q "T fourierlosssolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/FourierLoss.so" | grep -q "T fourierlosssolver_init_" \
  || { echo "FAIL: staged port .so lost its init entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/hutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" || true
}

echo "-- differential A (Direct/exact/summed + loss files) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat bodyloss.dat seriesloss.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars/loss files)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (curl/inexact-Simpsons/separate components)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stderr.txt scalars.dat; do
  diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" | head; exit 1; }
done
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config B stdout diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (A vs B differ; both arms agree per config)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: A/B stdout identical, branches not exercised"; exit 1; }
diff "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: A/B scalars identical, branches not exercised"; exit 1; }
echo "   liveness control GREEN (configs diverge, arms agree)"

echo "FOURIERLOSS PROOF GREEN ($((SECONDS - T0))s wall)"
