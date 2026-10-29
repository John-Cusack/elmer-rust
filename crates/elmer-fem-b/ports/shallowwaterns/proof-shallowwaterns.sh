#!/bin/bash
# T021-B ShallowWaterNS solver-entry proof: the Rust cdylib replaces ShallowWaterNS.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `shallowwaternssolver_`);
#  2. two distinct .sif differentials (config A: Mass Flux Discretization
#     False, upstream-like; config B: Mass Flux Discretization True,
#     i.e. the boundary-flux branch) are byte-identical to the reference
#     except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch-liveness control: config B output differs from config A
#     output while both arms agree on each config.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config —
# both arms share it, so the differential is self-consistent) and the
# Elmer checkout only for the 1d.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-shallowwaterns"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/1d "$WORK"/1d.grd
echo "== T021-B ShallowWaterNS proof: work=$WORK build=$ELMER_BUILD"
echo "File started at: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
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
PORT_SO="$PORT_DIR/target/release/libshallowwaterns.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ShallowWaterNS.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ShallowWaterNS module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T shallowwaternssolver_" "$WORK/port.syms" \
  || { echo "FAIL: shallowwaternssolver_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's ShallowWaterNS 1d.grd)"
cp "$ELMER_FORK/fem/tests/ShallowWaterNS/1d.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 1d > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "1d"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Transient
  Timestep Intervals = 2
  Timestep Sizes = 0.01
  Timestepping Method = BDF
  BDF Order = 1
  Output Intervals(1) = 1
  Steady State Max Iterations = 1
  Max Output Level = 4
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
  Initial Condition = 1
End
Initial Condition 1
  F = Variable Coordinate 1;Real MATC "0.005*exp(-20*((tx-1)^2))"
  Ground = Real -0.03
End
Equation 1
  Active Solvers(2) = 1 2
End
Body Force 1 :: F Lower Limit = Equals Ground
Material 1
  Wind-Stress Coefficient = Real 1
  Wind 1 = Real 1
  Coriolis Coefficient = Real 0
  Graviational Acceleration = Real 0
  Viscosity Coefficient = Real 1
  Inverse Chezy Coefficient = Real 0
End
Solver 1
  Equation = "Shallow"
  Variable = POT[U:1 F:1]
  Procedure = "ShallowWaterNS" "ShallowWaterNSSolver"
  Element = "p:2"
  Exported Variable 1 = Ground
  Apply Limiter = True
  Mass Flux Discretization = Logical False
  Nonlinear System Convergence Tolerance = 1e-6
  Nonlinear System Max Iterations = 25
  Nonlinear System Newton After Iterations = 1
  Nonlinear System Newton After Tolerance = 1e-9
  Nonlinear System Relaxation Factor = 1.0
  Linear System Solver = Iterative
  Linear System Iterative Method = BicgstabL
  Bicgstabl Polynomial Degree = 2
  Linear System Max Iterations = 200
  Linear System Residual Output = 40
  Linear System Preconditioning = Ilu
  Linear System Convergence Tolerance = 1e-8
  Steady State Convergence Tolerance = 1e-6
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = POT
  Operator 1 = max
End
Boundary Condition 1
  Target Boundaries = 1
End
Boundary Condition 2
  Target Boundaries = 2
  Periodic BC = 1
  Periodic BC F = Logical True
  Periodic BC U = Logical True
End
EOF
sed -e 's/Mass Flux Discretization = Logical False/Mass Flux Discretization = Logical True/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -q "Logical True" "$WORK/caseB.sif" \
  || { echo "FAIL: caseB did not take the Mass Flux branch"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ShallowWaterNS.so"
cp "$PORT_SO" "$WORK/shadow-port/ShallowWaterNS.so"
nm -D --defined-only "$WORK/shadow-port/ShallowWaterNS.so" | grep -q "T shallowwaternssolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/1d" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (MassFlux False branch) + rerun (re-entry)"
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

echo "-- differential B (MassFlux True branch, boundary flux live)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A and B must differ; both arms agree)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
if diff "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1; then
  echo "FAIL: branch control produced identical output (MassFlux branch dead?)"
  exit 1
fi
echo "   A-vs-B outputs differ as expected (MassFlux branch live)"

echo "SHALLOWWATERNS PROOF GREEN"
echo "File finished at: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "Wall seconds: $((SECONDS - T0))"
