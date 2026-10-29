#!/bin/bash
# T021-B VorticitySolver solver-entry proof: the Rust cdylib replaces
# VorticitySolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `vorticitysolver_` + `vorticitysolver_init_`);
#  2. two distinct .sif differentials (config A: default `Velocity`
#     target variable, lid 1.0, 1 steady iteration; config B: explicit
#     `Vorticity Variable`, `Constant Bulk Matrix`, lid 2.0, 2 steady
#     iterations) are byte-identical to the reference except
#     wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch-liveness control: config B output differs from config A,
#     and both arms agree on it.
#
# Wall-clock note: besides the framework lines the lane strips
# (`Elapsed REAL time`, `SOLVER TOTAL TIME`, `ELMER FINISHED AT`,
# `STARTED AT`), this unit and its driver solver report `RealTime()`- and
# CPU-clock-derived values on stdout (`VorticitySolver: Assembly Time:`,
# `VorticitySolver: Solution Time:`, `FlowSolve: iter: N Assembly: (s)`).
# These call the clock and cannot be bit-identical by construction, in
# either arm; the norm filter drops exactly those clock-derived lines.
# Every computational line (norms, residuals, iteration counts) and
# `scalars.dat` are compared byte-identically.
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
WORK="$LANE/scratch/proof-vorticitysolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B VorticitySolver proof: work=$WORK build=$ELMER_BUILD"
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
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libvorticitysolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/VorticitySolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the VorticitySolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T vorticitysolver_" "$WORK/port.syms" \
  || { echo "FAIL: vorticitysolver_ not the exported entry"; exit 1; }
grep -qx "T vorticitysolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: vorticitysolver_init_ not exported"; exit 1; }
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
  Viscosity = 0.1
End
Boundary Condition 1
  Target Boundaries(3) = 1 2 4
  Velocity 1 = Real 0.0
  Velocity 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  Velocity 1 = Real 1.0
  Velocity 2 = Real 0.0
End
Solver 1
  Equation = Navier-Stokes
  Stabilization Method = Stabilized
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU1
  Nonlinear System Convergence Tolerance = 1.0e-6
  Nonlinear System Max Iterations = 20
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Newton After Tolerance = 1.0e-3
  Nonlinear System Relaxation Factor = 1.0
  Steady State Convergence Tolerance = 1.0e-5
End
Solver 2
  Equation = Vorticity
  Procedure = "VorticitySolver" "VorticitySolver"
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Vorticity
  Operator 1 = max
  Variable 2 = Vorticity
  Operator 2 = min
End
EOF
# Config B: different ValueList branches — explicit target variable name,
# constant bulk matrix path, predeclared solver variable (Init
# ListCheckPresent branch), faster lid, 2 steady iterations (re-entry).
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    -e 's/Velocity 1 = Real 1.0/Velocity 1 = Real 2.0/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
t = t.replace("""Solver 2
  Equation = Vorticity
  Procedure = "VorticitySolver" "VorticitySolver"
End""", """Solver 2
  Equation = Vorticity
  Procedure = "VorticitySolver" "VorticitySolver"
  Variable = Vorticity
  Variable DOFs = 1
  Vorticity Variable = String "Velocity"
  Constant Bulk Matrix = Logical True
End""", 1)
open(p, "w").write(t)
PYEOF
grep -c "Constant Bulk Matrix" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/VorticitySolver.so"
cp "$PORT_SO" "$WORK/shadow-port/VorticitySolver.so"
nm -D --defined-only "$WORK/shadow-port/VorticitySolver.so" | grep -q "T vorticitysolver_" \
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
  # The brief's set plus this tree's spelling (`ELMER SOLVER FINISHED AT`);
  # every remaining clock-derived line in these runs carries `(s)`
  # (FlowSolve per-iteration Assembly/Solve timings, ILU factorization
  # time) or is this unit's own RealTime report (Assembly/Solution Time).
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT|SOLVER FINISHED AT|Assembly Time:|Solution Time:|\(s\)" "$1" || true
}

echo "-- differential A (default branches) + rerun (re-entry)"
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

echo "-- differential B (explicit branches + 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B output must differ from A)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
if diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1; then
  echo "FAIL: config B output identical to config A (control vacuous)"; exit 1;
fi
if diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: config B scalars identical to config A (control vacuous)"; exit 1;
fi
echo "   control live: B differs from A, both arms agree"

echo "-- release build + clippy clean"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release 2>&1 | tee "$WORK/clippy.txt" | grep -E "^(warning|error)" && { echo "FAIL: clippy warnings"; exit 1; } || true)
echo "   build + clippy clean"

echo "VORTICITYSOLVER PROOF GREEN (wall $((SECONDS - T0))s)"
