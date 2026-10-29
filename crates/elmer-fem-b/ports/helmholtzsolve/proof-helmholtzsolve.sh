#!/bin/bash
# T022-A HelmholtzSolve solver-entry proof: the Rust cdylib replaces HelmholtzSolve.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `helmholtzsolver_`, `helmholtzsolver_init_`);
#  2. two distinct .sif differentials (base config A; PlaneWave-flux +
#     Bubbles config B with 2 nonlinear iterations, i.e. the
#     repeated-invocation run) are byte-identical to the reference
#     except wall-clock lines/values;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. config B takes different code branches with observably different
#     output agreed on by both arms (liveness control).
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
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-helmholtzsolve"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== HelmholtzSolve proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran 13.3.0 required, got: $GFORTRAN_VERSION"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: rustc nightly-1.100.0 required, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libhelmholtzsolve.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/HelmholtzSolve.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the HelmholtzSolve module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T helmholtzsolver_" "$WORK/port.syms" \
  || { echo "FAIL: helmholtzsolver_ not an exported entry"; exit 1; }
grep -qx "T helmholtzsolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: helmholtzsolver_init_ not an exported entry"; exit 1; }
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
  Max Output Level = 4
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
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
  Active Solvers(1) = Integer 1
End
Material 1
  Density = Real 1.2
  Sound Speed = Real 340.0
  Sound Damping = Real 0.0
End
Body Force 1
  Pressure Source 1 = Real 10.0
  Pressure Source 2 = Real 0.0
End
Boundary Condition 1
  Target Boundaries(1) = Integer 1
  Pressure 1 = Real 1.0
  Pressure 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = Integer 2
  Wave Impedance 1 = Real 340.0
  Wave Impedance 2 = Real 0.0
End
Boundary Condition 3
  Target Boundaries(1) = Integer 3
End
Boundary Condition 4
  Target Boundaries(1) = Integer 4
  Wave Flux 1 = Real 1.0
  Wave Flux 2 = Real 0.5
End
Solver 1
  Equation = Helmholtz
  Procedure = "HelmholtzSolve" "HelmholtzSolver"
  Variable = Pressure
  Variable DOFs = Integer 2
  Frequency = Real 1000.0
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = Integer 500
  Linear System Convergence Tolerance = Real 1.0e-8
  Linear System Preconditioning = ILU0
  Nonlinear System Max Iterations = Integer 1
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Pressure
  Operator 1 = min
  Variable 2 = Pressure
  Operator 2 = max
End
EOF
# Config B: PlaneWave branch (impedance from parent Sound Speed),
# Wave Flux on the idle boundary, Bubbles assembly path, and 2 nonlinear
# iterations (in-process SAVE re-entry across invocations).
sed -e 's/  Wave Impedance 1 = Real 340.0/  Plane Wave BC = Logical True/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
t = t.replace("""Boundary Condition 3
  Target Boundaries(1) = Integer 3
End""", """Boundary Condition 3
  Target Boundaries(1) = Integer 3
  Wave Flux 1 = Real 2.0
  Wave Flux 2 = Real -1.0
End""", 1)
t = t.replace("  Nonlinear System Max Iterations = Integer 1",
              "  Nonlinear System Max Iterations = Integer 2\n  Bubbles = Logical True", 1)
open(p, "w").write(t)
PYEOF
grep -c "Plane Wave BC" "$WORK/caseB.sif"
grep -c "Bubbles = Logical True" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/HelmholtzSolve.so"
cp "$PORT_SO" "$WORK/shadow-port/HelmholtzSolve.so"
nm -D --defined-only "$WORK/shadow-port/HelmholtzSolve.so" | grep -q "T helmholtzsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/HelmholtzSolve.so" | grep -q "T helmholtzsolver_init_" \
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
norm() { # drop wall-clock lines/values for the byte-compare; the
         # Helmholtz `iter: N Assembly/Solve (s)` timing values are
         # normalized with the lane `: +\(s\)` mask (values only,
         # structure/iterations still compared).
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" | \
    sed -E 's/(: +\(s\)).*/\1/' || true
}

echo "-- differential A (base) + rerun (re-entry)"
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

echo "-- differential B (PlaneWave + flux + Bubbles, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness control (B takes different branches with different output)"
if diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: config B output identical to config A (branches not exercised)"
  exit 1
fi
if diff -q "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: reference config B output identical to config A"
  exit 1
fi
echo "   branches live: A and B differ, both arms agree"

echo "HELMHOLTZSOLVE PROOF GREEN (wall $((SECONDS - T0))s)"
