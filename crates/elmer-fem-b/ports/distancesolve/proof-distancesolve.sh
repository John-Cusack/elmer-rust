#!/bin/bash
# T022-A DistanceSolve solver-entry proof: the Rust cdylib replaces DistanceSolve.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `distancesolver_`, `distancesolver1_`,
#     `distancesolver_init_`);
#  2. three distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines: config A (DistanceSolver PDE branch, with
#     RelaxDT terms and 2 steady iterations, i.e. the in-process
#     repeated-invocation run through the SAVE state); config B
#     (DistanceSolver1 geometric branch with Dummy Distance Computation,
#     i.e. distcomp0); config C (DistanceSolver1 with the real distcomp
#     branch: RANDOM_NUMBER + SortR search);
#  3. a second full process run of config B against the port is
#     byte-identical (idempotent re-entry);
#  4. branch-liveness control: configs B and C take different code
#     branches and produce different output, and both arms agree on each.
#
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:-/home/john/firstmate/data/rs-elmer-t020/scratch/elmer-build}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-distancesolve"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK/square"*
echo "== DistanceSolve proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libdistancesolve.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- clippy clean"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release -- -D warnings) || \
  { echo "FAIL: clippy warnings"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/DistanceSolve.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the DistanceSolve module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in distancesolver_ distancesolver1_ distancesolver_init_; do
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
  Steady State Min Iterations = 2
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
  Density = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Distance = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Distance = Real 0.0
End
Solver 1
  Equation = DistancePDE
  Procedure = "DistanceSolve" "DistanceSolver"
  Variable = Distance
  Variable DOFs = 1
  Nonlinear System Max Iterations = 50
  Distance Pseudo dt = Real 0.05
  H Scale = Real 1.0
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
  Variable 1 = Distance
  Operator 1 = min
  Variable 2 = Distance
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
  Max Output Level = 4
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
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Body Force 1
  Distance = Real 0.0
  Distance Condition = Real 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Distance = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Noslip Wall BC = Logical True
End
Solver 1
  Equation = DistanceGeom
  Procedure = "DistanceSolve" "DistanceSolver1"
  Variable = Distance
  Variable DOFs = 1
  Dummy Distance Computation = Logical True
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Distance
  Operator 1 = min
  Variable 2 = Distance
  Operator 2 = max
End
EOF
sed -e 's/Dummy Distance Computation = Logical True/Dummy Distance Computation = Logical False/' \
  "$WORK/caseB.sif" > "$WORK/caseC.sif"
grep -q "Logical False" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/DistanceSolve.so"
cp "$PORT_SO" "$WORK/shadow-port/DistanceSolve.so"
nm -D --defined-only "$WORK/shadow-port/DistanceSolve.so" | grep -q "T distancesolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/DistanceSolve.so" | grep -q "T distancesolver1_" \
  || { echo "FAIL: staged port .so lost its geometric entry"; exit 1; }

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
norm() { # merged wall-clock norm class: drop wall-clock lines, mask FlowSolve per-iteration (s) values
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" | \
    sed -E -e 's/Assembly: \(s\)[^0-9]*[0-9.Ee+-]+/Assembly: (s) MASKED/' \
           -e 's/Solve: \(s\)[^0-9]*[0-9.Ee+-]+/Solve: (s) MASKED/' || true
}

echo "-- differential A (PDE branch, RelaxDT terms, 2 steady iterations)"
run_case caseA orig runA-orig
run_case caseA port runA-port
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"

echo "-- differential B (geometric dummy branch: distcomp0 + bulk/bc fixed nodes)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- rerun B against the port (re-entry)"
run_case caseB port runB-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port2/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (re-entry)"

echo "-- differential C (geometric real branch: RANDOM_NUMBER + SortR distcomp)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"

echo "-- liveness control (B dummy vs C real-dist branches differ, both arms agree)"
if diff -q "$WORK/runB-port/scalars.dat" "$WORK/runC-port/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: liveness control produced identical output (branches not exercised)"
  exit 1
fi
if diff -q "$WORK/runB-orig/scalars.dat" "$WORK/runC-orig/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: reference arms agree too much (branches not exercised)"
  exit 1
fi
echo "   liveness control passes (branches differ, arms agree)"

echo "DISTANCESOLVE PROOF GREEN (wall $((SECONDS - T0))s)"
