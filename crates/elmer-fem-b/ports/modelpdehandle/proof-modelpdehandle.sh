#!/bin/bash
# T021-B ModelPDEhandle solver-entry proof: the Rust cdylib replaces
# ModelPDEhandle.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `advdiffsolver_`);
#  2. two distinct .sif differentials (config A: diffusion+source+reaction
#     steady branch; config B: convection-velocity + Robin-BC + 2
#     nonlinear iterations branch) are byte-identical to the reference
#     except wall-clock lines;
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
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-modelpdehandle"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== ModelPDEhandle proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran pin 13.3.0 not met"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: rustc pin 1.100.0-nightly not met"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libmodelpdehandle.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ModelPDEhandle.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ModelPDEhandle module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T advdiffsolver_" "$WORK/port.syms" \
  || { echo "FAIL: advdiffsolver_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"
echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
# NOTE: the lane default (PoissonDB square.grd) is triangular, and the
# reference ModelPDEhandle.so itself aborts on it (mGetElementDOFS bounds
# error on the point/line boundary elements), so no differential is
# possible there. This proof uses the solver's own upstream quad square
# mesh (fem/tests/ModelPDEhandle/square.grd), same `ElmerGrid 1 2` flow.
cp "$ELMER_FORK/fem/tests/ModelPDEhandle/square.grd" "$WORK/"
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
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Diffusion Coefficient = Real 1.0
  Reaction Coefficient = Real 0.5
  Time Derivative Coefficient = Real 1.0
End
Body Force 1
  Field Source = Real 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Field = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Field = Real 1.0
End
Solver 1
  Equation = ModelPDE
  Procedure = "ModelPDEhandle" "AdvDiffSolver"
  Variable = Field
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Field
  Operator 1 = min
  Variable 2 = Field
  Operator 2 = max
End
EOF
# Config B: convection + Robin branch, 2 nonlinear iterations, direct
# linear solve (the coarse mesh lets strong convection diverge the
# iterative solver in *both* arms, so B pins the linear path to Umfpack).
sed -e 's/Reaction Coefficient = Real 0.5/Reaction Coefficient = Real 0.25\n  Convection Coefficient = Real 0.5\n  Convection Velocity 1 = Real 0.5\n  Convection Velocity 2 = Real 0.25/' \
    -e 's/  Equation = ModelPDE/  Equation = ModelPDE\n  Nonlinear System Max Iterations = 2/' \
    -e 's/Linear System Solver = Iterative/Linear System Solver = Direct/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
t = t.replace("""Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Field = Real 1.0
End""", """Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Field Flux = Real 0.5
  Robin Coefficient = Real 10.0
  External Field = Real 1.0
End""", 1)
open(p, "w").write(t)
PYEOF
grep -c "Robin Coefficient" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ModelPDEhandle.so"
cp "$PORT_SO" "$WORK/shadow-port/ModelPDEhandle.so"
nm -D --defined-only "$WORK/shadow-port/ModelPDEhandle.so" | grep -q "T advdiffsolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER SOLVER FINISHED AT|ELMER FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (diffusion+source+reaction branch) + rerun (re-entry)"
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
echo "-- differential B (convection + Robin + 2 nonlinear iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A and B must differ)"
if diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: branch-liveness control: A and B scalars identical"
  exit 1
fi
echo "   A/B outputs differ as expected, both arms agree"

echo "MODELPDHANDLE PROOF GREEN in $((SECONDS - T0))s"
