#!/bin/bash
# T022-A FluidicForce solver-entry proof: the Rust cdylib replaces FluidicForce.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `forcecompute_`; invoked via the campaign
#     `Procedure = "FluidicForce" "ForceCompute"` library-vs-entry form);
#  2. two distinct .sif differentials (viscous config A; inviscid +
#     moment + shear-output config B with 2 steady iterations, i.e. the
#     repeated-invocation run) are byte-identical to the reference
#     except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry; the original keeps no SAVE state);
#  4. branch liveness: config B output differs from config A output,
#     and both arms agree on each config (incl. the shear-stress files).
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
WORK="$LANE/scratch/proof-fluidicforce"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square* "$WORK"/*.syms
echo "== T022-A FluidicForce proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libfluidicforce.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FluidicForce.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FluidicForce module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T forcecompute_" "$WORK/port.syms" \
  || { echo "FAIL: forcecompute_ not the exported entry"; exit 1; }
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
  NS Convect = False
End
Material 1
  Density = 1.0
  Viscosity = 0.1
End
Solver 1
  Equation = N-S
  Variable = Flow Solution[Velocity:2 Pressure:1]
  Procedure = "FlowSolve" "FlowSolver"
  Stabilize = True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Max Iterations = 200
  Linear System Convergence Tolerance = 1.0e-8
  Nonlinear System Max Iterations = 1
  Nonlinear System Convergence Tolerance = 1.0e-5
  Steady State Convergence Tolerance = 1.0e-5
End
Solver 2
  Equation = Fluidic Force
  Procedure = "FluidicForce" "ForceCompute"
  Calculate Viscous Force = True
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = String "res: fluid force bc 1"
  Operator 1 = min
  Variable 2 = String "res: fluid force area bc 1"
  Operator 2 = min
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Velocity 1 = 0.0
  Velocity 2 = 0.0
  Calculate Fluidic Force = True
End
Boundary Condition 2
  Target Boundaries(2) = 2 4
  Velocity 1 = 0.0
  Velocity 2 = 0.0
End
Boundary Condition 3
  Target Boundaries(1) = 3
  Velocity 1 = 1.0
  Velocity 2 = 0.0
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "square"
End
Simulation
  Max Output Level = 4
  Coordinate System = Cartesian 2D
  Simulation Type = Steady State
  Steady State Max Iterations = 2
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
  NS Convect = False
End
Material 1
  Density = 1.0
  Viscosity = 0.1
  Compressibility Model = String "Incompressible"
End
Solver 1
  Equation = N-S
  Variable = Flow Solution[Velocity:2 Pressure:1]
  Procedure = "FlowSolve" "FlowSolver"
  Stabilize = True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Max Iterations = 200
  Linear System Convergence Tolerance = 1.0e-8
  Nonlinear System Max Iterations = 1
  Nonlinear System Convergence Tolerance = 1.0e-5
  Steady State Convergence Tolerance = 1.0e-5
End
Solver 2
  Equation = Fluidic Force
  Procedure = "FluidicForce" "ForceCompute"
  Calculate Viscous Force = False
  Velocity Field Name = String "Flow Solution"
  Shear Stress Output = True
  Shear Stress Output File = shearstress.dat
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = String "res: fluid force bc 1"
  Operator 1 = min
  Variable 2 = String "res: fluid force area bc 1"
  Operator 2 = min
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Velocity 1 = 0.0
  Velocity 2 = 0.0
  Calculate Fluidic Force = True
  Moment About 1 = Real 0.5
  Moment About 2 = Real 0.25
  Moment About 3 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(2) = 2 4
  Velocity 1 = 0.0
  Velocity 2 = 0.0
End
Boundary Condition 3
  Target Boundaries(1) = 3
  Velocity 1 = 2.0
  Velocity 2 = 0.0
End
EOF
grep -q "Calculate Viscous Force = False" "$WORK/caseB.sif"
grep -q "Shear Stress Output = True" "$WORK/caseB.sif"
grep -q "Moment About 1" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/FluidicForce.so"
cp "$PORT_SO" "$WORK/shadow-port/FluidicForce.so"
nm -D --defined-only "$WORK/shadow-port/FluidicForce.so" | grep -q "T forcecompute_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" || true
}

echo "-- differential A (viscous branch) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat scalars.dat.names; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
grep -q "Resultant Force" "$WORK/a1.txt" \
  || { echo "FAIL: config A shows no force output (solver may not have run)"; exit 1; }
! grep -q "Moment about" "$WORK/a1.txt" \
  || { echo "FAIL: config A unexpectedly shows moment output"; exit 1; }
for f in stdout.txt stderr.txt scalars.dat scalars.dat.names; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (inviscid + moment + shear, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat scalars.dat.names shearstress.dat shearstress.dat.names; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars/shear files)"
[ "$(grep -c "Moment about" "$WORK/runB-port/stdout.txt")" -eq 2 ] \
  || { echo "FAIL: config B lacks the 2 expected moment reports"; exit 1; }

echo "-- branch liveness (B must differ from A; both arms must agree)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
! diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1 \
  || { echo "FAIL: config B output identical to config A (branch dead?)"; exit 1; }
echo "   B differs from A on both arms (liveness control holds)"

echo "FLUIDICFORCE PROOF GREEN (${SECONDS}s wall)"
