#!/bin/bash
# T022-A DFT solver-entry proof: the Rust cdylib replaces DFTSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `chargedensitysolver_`, `wavefunctionsolver_`,
#     `poissonsolver_`);
#  2. two distinct .sif differentials (config A: exponential mixing +
#     Perdew-Zunger XC, weights given; config B: constant mixing +
#     XC none + shifting constant set + default weights) are
#     byte-identical to the reference except wall-clock lines/values;
#  3. a second full process run of config A against the port is
#     byte-identical (SAVE re-entry through fresh first-time state);
#  4. branch-liveness control: config B takes different code branches
#     (constant mixing, XC-none, shift PRINT, default weights) and both
#     arms agree on it, while A and B differ from each other.
#
# Defined-path notes (documented, not chased):
# - `PoissonSolver`'s `RF` local is assigned only in the constant/
#   exponential mixing arms and read by the GRPULAY3 gate on every pass.
#   Both proof configs set a mixing method, so the gate only reads a
#   freshly assigned value. `RelaxType` 0/3-on-entry paths (stale RF) are
#   not exercised; the GRPulay activation thresholds are set far out of
#   reach (`Start GRPulay after iterations = 1000000`), so GRPULAY3 only
#   ever takes its store-only path (identical by construction, in the
#   Fortran wrapper).
# - The commented-out MPI alternatives in the original are skipped
#   (serial proof tree); the serial branches are transcribed.
# - `PoissonSolver`'s dead `EigsVar => VariableGet(...)` pointer query
#   (never read afterwards) is skipped: a pure lookup with no side
#   effects.
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
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-dftm"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A DFTM proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libdftm.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/DFTSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the DFTSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in chargedensitysolver_ wavefunctionsolver_ poissonsolver_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
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
End
Body 1
  Equation = 1
End
Body 2
  Equation = 1
End
Body 3
  Equation = 1
End
Body 4
  Equation = 1
End
Equation 1
  Active Solvers(4) = 1 2 3 4
End
Solver 1
  Equation = Poisson
  Procedure = "DFTSolver" "PoissonSolver"
  Variable = Potential
  Variable DOFs = 1
  NOFnuclei = Integer 3
  NucleiTable(3,4) = Real 8.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 1.0 0.0
  Linear System Solver = Iterative
  Linear System Iterative Method = CG
  Linear System Max Iterations = 2000
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = 1.0e-10
  Linear System Symmetric = Logical True
  Linear System Residual Output = 0
  Relaxation Method = String "exponential mixing"
  Relaxation Parameter A = Real 1.0
  Relaxation Parameter B = Real 0.05
  Relaxation Parameter C = Real 0.005
  Start GRPulay after iterations = Integer 1000000
  Start GRPulay if relaxation factor is more than = Real 1.0e8
  Steady State Convergence Tolerance = 1e-6
End
Solver 2
  Equation = EigenProblem
  Variable = WaveFunctions
  Procedure = "DFTSolver" "WaveFunctionSolver"
  Eigen Analysis = Logical True
  Eigen System Values = Integer 4
  Eigen System Select = String "Smallest magnitude"
  Eigen System Convergence Tolerance = Real 1.0e-7
  Eigen System Compute Residuals = Logical True
  Eigen System Lanczos Vectors = Integer 20
  Eigen System Max Iterations = Integer 500
  Linear System Iterative Method = CG
  Linear System Solver = Iterative
  Linear System Max Iterations = 10000
  Linear System Preconditioning = ILU1
  Linear System Convergence Tolerance = 1e-9
  Linear System Residual Output = 0
  Linear System Symmetric = Logical True
  Linear System Abort Not Converged = Logical False
  XC Potential type = String "perdew-zunger"
End
Solver 3
  Equation = Charge density calculation
  Variable = Charge density
  Procedure = "DFTSolver" "ChargeDensitySolver"
  Number of Eigenmodes Included = Integer 3
  Weights of Eigen States(3,1) = Real 2.0 2.0 2.0
  Linear System Solver = Iterative
  Linear System Iterative Method = CG
  Linear System Max Iterations = 10000
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = 1e-9
  Linear System Residual Output = 0
  Linear System Symmetric = Logical True
  Linear System Abort Not Converged = Logical False
  Steady State Convergence Tolerance = 1e-4
End
Solver 4
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
  Variable 3 = Charge density
  Operator 3 = min
  Variable 4 = Charge density
  Operator 4 = max
End
Boundary Condition 1
  Target Boundaries(1) = 1
  WaveFunctions = Real 0.0
  Potential = Real 100.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  WaveFunctions = Real 0.0
  Potential = Real 100.0
End
EOF
# Config B: constant mixing + XC none + shifting constant set + default
# eigenstate weights (exercises the other branch of each choice).
python3 - "$WORK/caseA.sif" "$WORK/caseB.sif" <<'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
t = open(src).read()
t = t.replace('Relaxation Method = String "exponential mixing"',
              'Relaxation Method = String "constant mixing"')
t = t.replace('Relaxation Parameter A = Real 1.0',
              'Relaxation Parameter A = Real 0.1')
t = t.replace('XC Potential type = String "perdew-zunger"',
              'XC Potential type = String "none"\n  Eigen System Shifting Constant = Real 5.0')
t = t.replace('  Weights of Eigen States(3,1) = Real 2.0 2.0 2.0\n', '')
open(dst, 'w').write(t)
PYEOF
grep -q 'Relaxation Method = String "constant mixing"' "$WORK/caseB.sif"
grep -q 'Shifting Constant' "$WORK/caseB.sif"
! grep -q 'Weights of Eigen States' "$WORK/caseB.sif" \
  || { echo "FAIL: caseB weights line not removed"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/modules/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/DFTSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/DFTSolver.so"
for sym in chargedensitysolver_ wavefunctionsolver_ poissonsolver_; do
  nm -D --defined-only "$WORK/shadow-port/DFTSolver.so" | grep -q "T $sym" \
    || { echo "FAIL: staged port .so lost entry $sym"; exit 1; }
done

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
norm() { # drop wall-clock lines, mask wall-clock values for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT" "$1" \
    | sed -E 's/\(s\):[[:space:]]+[0-9.]+/(s): X/' || true
}

echo "-- differential A (exponential mixing + Perdew-Zunger) + rerun (re-entry)"
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
grep -q "Applied exponential mixing" "$WORK/runA-orig/stdout.txt" \
  || { echo "FAIL: config A missed exponential-mixing branch"; exit 1; }
grep -q "Using Perdew-Zunger type XC-potential" "$WORK/runA-orig/stdout.txt" \
  || { echo "FAIL: config A missed XC branch"; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (constant mixing + XC none + shift + default weights)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
grep -q "Applied constant mixing" "$WORK/runB-port/stdout.txt" \
  || { echo "FAIL: config B missed constant-mixing branch"; exit 1; }
grep -q "Not using any type of XC-potential" "$WORK/runB-port/stdout.txt" \
  || { echo "FAIL: config B missed XC-none branch"; exit 1; }
grep -q "Eigen System Shifting constant was set to" "$WORK/runB-port/stdout.txt" \
  || { echo "FAIL: config B missed shift PRINT branch"; exit 1; }
grep -q "All weights were set to 1" "$WORK/runB-port/stdout.txt" \
  || { echo "FAIL: config B missed default-weights branch"; exit 1; }
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness: A and B differ (distinct branches, distinct outputs)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
if diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1; then
  echo "FAIL: configs A and B agree (no branch discrimination)"; exit 1;
fi
echo "   liveness control GREEN"

echo "DFTM PROOF GREEN (wall: $((SECONDS - T0))s)"
