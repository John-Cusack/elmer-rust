#!/bin/bash
# T021-B Stokes solver-entry proof: the Rust cdylib replaces Stokes.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `stokessolver_`);
#  2. two distinct .sif differentials (steady non-convective config A;
#     convective Newton config B with body force, 3 nonlinear iterations,
#     i.e. the repeated-invocation run) are byte-identical to the reference
#     modulo wall-clock/CPU-time lines;
#  3. a second full process run of config A is byte-identical through
#     norm (idempotent re-entry through the SAVE state);
#  4. branch liveness: config B output differs from config A output on
#     both arms.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout only for the Step.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-stokes"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/Step "$WORK"/Step.grd
echo "== T021-B Stokes proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran must be 13.3.0"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: rustc +nightly must be 1.100.0-nightly"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libstokes.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- clippy (no warnings)"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release 2>&1 | tee "$WORK/clippy.txt" | grep -E "^(warning|error)" && { echo "FAIL: clippy warnings"; exit 1; } || true)
grep -E "^(warning|error)" "$WORK/clippy.txt" && { echo "FAIL: clippy warnings"; exit 1; } || echo "   clippy clean"

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/Stokes.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the Stokes module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T stokessolver_" "$WORK/port.syms" \
  || { echo "FAIL: stokessolver_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's Step.grd)"
cp "$ELMER_FORK/fem/tests/SD_Step_stokes_block/Step.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 Step > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "Step"
End
Simulation
  Coordinate System = Cartesian 2D
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Steady
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
  Density = 1
  Viscosity = 1
End
Boundary Condition 1
  Target Boundaries = 1
  RotSol 1 = Variable Coordinate 2
     Real MATC "(2-tx)*(tx-1)"
  RotSol 2 = Real 0
End
Boundary Condition 2
  Target Boundaries = 2
  Outflow Boundary = Logical True
  RotSol 2 = Real 0
End
Boundary Condition 3
  Target Boundaries = 3
  RotSol 1 = Real 0
  RotSol 2 = Real 0
End
Solver 1
  Equation = "Stokes"
  Procedure = "Stokes" "StokesSolver"
  Element = "p:2"
  Variable = RotSol
  Variable Dofs = 3
  Convective = Logical False
  Block Preconditioning = Logical False
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Residual Output = 1000
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = RotSol
  Operator 1 = min
  Variable 2 = RotSol
  Operator 2 = max
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "Step"
End
Simulation
  Coordinate System = Cartesian 2D
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Steady
  Steady State Max Iterations = 1
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Density = 1
  Viscosity = 0.7
End
Body Force 1
  Force 1 = Real 1.0
  Force 2 = Real 0.5
End
Boundary Condition 1
  Target Boundaries = 1
  RotSol 1 = Variable Coordinate 2
     Real MATC "(2-tx)*(tx-1)"
  RotSol 2 = Real 0
End
Boundary Condition 2
  Target Boundaries = 2
  Outflow Boundary = Logical True
  RotSol 2 = Real 0
End
Boundary Condition 3
  Target Boundaries = 3
  RotSol 1 = Real 0
  RotSol 2 = Real 0
End
Solver 1
  Equation = "Stokes"
  Procedure = "Stokes" "StokesSolver"
  Element = "p:2"
  Variable = RotSol
  Variable Dofs = 3
  Convective = Logical True
  Nonlinear Iteration Method = String "Newton"
  Nonlinear System Max Iterations = 2
  Nonlinear System Convergence Tolerance = 1.0e-8
  Nonlinear System Newton After Iterations = 1
  Nonlinear System Newton After Tolerance = 1.0e-3
  Nonlinear System Relaxation Factor = 0.5
  Block Preconditioning = Logical False
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Residual Output = 1000
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = RotSol
  Operator 1 = min
  Variable 2 = RotSol
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
rm -f "$WORK/shadow-port/Stokes.so"
cp "$PORT_SO" "$WORK/shadow-port/Stokes.so"
nm -D --defined-only "$WORK/shadow-port/Stokes.so" | grep -q "T stokessolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/Step" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Assembly \\(s\\)|Solve \\(s\\)|pointer init time|Time spent for restart|Factorization time|Assembly:.*\(s\)|Solution:.*\(s\)" "$1" || true; }

echo "-- differential A (non-convective) + rerun (re-entry)"
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

echo "-- differential B (convective Newton + body force, 3 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head -20; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (B differs from A on both arms)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
diff "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: config B output identical to config A (dead branch?)"; exit 1; }
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b1.txt"
diff "$WORK/a0.txt" "$WORK/b1.txt" > /dev/null 2>&1 && \
  { echo "FAIL: orig config B output identical to config A"; exit 1; }
echo "   branch liveness confirmed (both arms)"

echo "STOKES PROOF GREEN ($((SECONDS - T0))s wall)"
