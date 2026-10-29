#!/bin/bash
# T021-B PressurePrecond solver-plugin proof: the Rust cdylib replaces
# PressurePrecond.so.
#
# Proves, on pinned gfortran 13.3.0 + lane nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `pressureprecond_` + `pressureprecond_init0_`);
#  2. two distinct .sif differentials (config A: Init0 defaults install
#     `Variable = P`; config B: the .sif pre-sets `Variable = Q` plus the
#     linear-system keys so Init0's `ListAddNew*` take the already-present
#     path, with a nonzero Initial Condition) are byte-identical to the
#     reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the no-op solver body);
#  4. branch liveness: config B's output differs from config A's
#     (variable q vs p, nonzero scalars), and both arms agree on it.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial configure —
# both arms share it, so the differential is self-consistent) and the
# Elmer checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-pressureprecond"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B PressurePrecond proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libpressureprecond.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/PressurePrecond.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the PressurePrecond module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T pressureprecond_" "$WORK/port.syms" \
  || { echo "FAIL: pressureprecond_ not an exported entry"; exit 1; }
grep -qx "T pressureprecond_init0_" "$WORK/port.syms" \
  || { echo "FAIL: pressureprecond_init0_ not an exported entry"; exit 1; }
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
  Steady State Max Iterations = 1
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Density = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  P = Real 0.0
End
Solver 1
  Equation = Pressure Preconditioning
  Procedure = "PressurePrecond" "PressurePrecond"
  Variable DOFs = 1
  Linear System Row Equilibration = True
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
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Density = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Q = Real 0.0
End
Solver 1
  Equation = Pressure Preconditioning
  Procedure = "PressurePrecond" "PressurePrecond"
  Variable = "Q"
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 7
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = 1.0e-5
  Linear System Abort Not Converged = Logical True
  Linear System Row Equilibration = True
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Q
  Operator 1 = min
  Variable 2 = Q
  Operator 2 = max
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
rm -f "$WORK/shadow-port/PressurePrecond.so"
cp "$PORT_SO" "$WORK/shadow-port/PressurePrecond.so"
nm -D --defined-only "$WORK/shadow-port/PressurePrecond.so" | grep -q "T pressureprecond_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/PressurePrecond.so" | grep -q "T pressureprecond_init0_" \
  || { echo "FAIL: staged port .so lost its Init0 entry"; exit 1; }

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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (Init0 defaults: Variable P installed by _Init0)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat scalars.dat.names; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
grep -q "Creating standard variable: p" "$WORK/runA-orig/stdout.txt" \
  || { echo "FAIL: config A never installed Variable P (Init0 not observed)"; exit 1; }
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat scalars.dat.names; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (re-entry through the no-op body)"

echo "-- differential B (sif pre-sets Variable Q + linear-system keys: ListAddNew keep-path)"
run_case caseB orig runB-orig
run_case caseB port runB-port
grep -q "Creating standard variable: q" "$WORK/runB-orig/stdout.txt" \
  || { echo "FAIL: config B lost the sif-set Variable Q"; exit 1; }
for f in stdout.txt stderr.txt scalars.dat scalars.dat.names; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (B must differ from A; arms already agree)"
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b1.txt"
if diff -q "$WORK/a0.txt" "$WORK/b1.txt" > /dev/null 2>&1; then
  echo "FAIL: config B output identical to config A (variant not live)"; exit 1;
fi
if diff -q "$WORK/runA-orig/scalars.dat.names" "$WORK/runB-orig/scalars.dat.names" > /dev/null 2>&1; then
  echo "FAIL: config B scalars names identical to config A (Init0 delta not observed)"; exit 1;
fi
echo "   variant is live and both arms agree on it"

echo "PRESSUREPRECOND PROOF GREEN (wall: $((SECONDS - T0))s)"
