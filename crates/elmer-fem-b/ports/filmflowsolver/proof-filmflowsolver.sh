#!/bin/bash
# T022-A FilmFlowSolver solver-entry proof: the Rust cdylib replaces
# FilmFlowSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: filmflowsolver_, filmflowsolver_init_,
#     filmflowsolver_init0_);
#  2. two distinct .sif differentials (config A: laminar 2-D channel,
#     single nonlinear iteration; config B: Convect + GradP + darcy
#     friction + newton + friction/pressure heating, 2 steady iterations
#     i.e. the repeated-invocation run) are byte-identical to the
#     reference except wall-clock lines;
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
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-filmflowsolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== FilmFlowSolver proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0, have: $GFORTRAN_VERSION"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: need rustc nightly 1.100.0, have: $RUSTC_VERSION"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libfilmflowsolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FilmFlowSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FilmFlowSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T filmflowsolver_" "$WORK/port.syms" \
  || { echo "FAIL: filmflowsolver_ not the exported entry"; exit 1; }
grep -qx "T filmflowsolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: filmflowsolver_init_ not exported"; exit 1; }
grep -qx "T filmflowsolver_init0_" "$WORK/port.syms" \
  || { echo "FAIL: filmflowsolver_init0_ not exported"; exit 1; }
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
  Coordinate System = Cartesian 2D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Max Output Level = 4
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
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0e-3
  Viscosity = 1.5
  Gap Height = 0.5
End
Boundary Condition 1
  Target Boundaries(1) = 1
  FilmVelocity 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 2
  FilmVelocity 2 = Real 0.0
End
Boundary Condition 3
  Target Boundaries(1) = 3
  FilmVelocity 2 = Real 0.0
End
Boundary Condition 4
  Target Boundaries(1) = 4
  FilmVelocity 1 = Real 1.3
End
Solver 1
  Equation = FilmFlow
  Procedure = "FilmFlowSolver" "FilmFlowSolver"
  Model Dimension = Integer 2
  Nonlinear System Max Iterations = 1
  Nonlinear System Convergence Tolerance = 1.0e-6
  Linear System Solver = Direct
  Linear System Direct Method = umfpack
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = FilmPressure
  Operator 1 = min
  Variable 2 = FilmPressure
  Operator 2 = max
  Variable 3 = FilmVelocity 1
  Operator 3 = min
  Variable 4 = FilmVelocity 1
  Operator 4 = max
End
EOF
# Config B: different physics branches — convection + GradP
# discretization + darcy friction model + newton linearization +
# friction heating with heat source, two steady iterations.
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    -e 's/Nonlinear System Max Iterations = 1/Nonlinear System Max Iterations = 2/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"

python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
t = t.replace("""  Model Dimension = Integer 2
""", """  Model Dimension = Integer 2
  Convect = True
  Gradp Discretization = True
  Friction model = String darcy
  Nonlinear System Newton After Iterations = 1
  Calculate Friction Heating = Logical True
  Use Heating Source = Logical True
""", 1)
t = t.replace("""  Gap Height = 0.5
""", """  Gap Height = 0.5
  Darcy Roughness = Real 1.0e-4
  Latent Heat = 1.0e6
""", 1)
t = t.replace("""Boundary Condition 1
  Target Boundaries(1) = 1
  FilmVelocity 2 = Real 0.0
End""", """Boundary Condition 1
  Target Boundaries(1) = 1
  FilmVelocity 2 = Real 0.0
  FilmPressure = Real 0.0
End""", 1)
open(p, "w").write(t)
PYEOF
grep -c "darcy\|Convect\|Heating" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/FilmFlowSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/FilmFlowSolver.so"
nm -D --defined-only "$WORK/shadow-port/FilmFlowSolver.so" | grep -q "T filmflowsolver_" \
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
norm() { # drop wall-clock lines, mask per-iteration (s) timing values
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" \
    | sed -E '/: +\(s\)/ s/[0-9]+\.[0-9]+/T/g' || true
}

echo "-- differential A (laminar branch) + rerun (re-entry)"
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

echo "-- differential B (convect+gradp+darcy+newton+heating, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness control (B takes different branches, output differs)"
if cmp -s "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat"; then
  echo "FAIL: liveness control: A and B scalars identical"
  exit 1
fi
echo "   A vs B outputs differ as expected"

echo "FILMFLOWSOLVER PROOF GREEN (wall $((SECONDS - T0))s)"
