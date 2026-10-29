#!/bin/bash
# T021-B TransientPhaseChange solver-entry proof: the Rust cdylib replaces
# TransientPhaseChange.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `transientphasechange_` + `transientphasechange_init_`);
#  2. two distinct clean-room .sif differentials (internal-flux + simple
#     update config A; nodal-loads + PDE-smoothing config B, each with 2
#     timesteps for in-run SAVE re-entry) are byte-identical to the
#     reference except wall-clock lines;
#  3. a second full process run of config A against the port is
#     byte-identical (idempotent re-entry through the SAVE state);
#  4. branch-liveness control: config B takes different code branches and
#     produces different output, and both arms agree on it.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-transientphasechange"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/PhaseChange*
echo "== T021B TransientPhaseChange proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libtransientphasechange.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/TransientPhaseChange.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T transientphasechange_" "$WORK/port.syms" \
  || { echo "FAIL: transientphasechange_ missing"; exit 1; }
grep -qx "T transientphasechange_init_" "$WORK/port.syms" \
  || { echo "FAIL: transientphasechange_init_ missing"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout PhaseChange3 test)"
cp "$ELMER_FORK/fem/tests/PhaseChange3/PhaseChange.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 PhaseChange > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "PhaseChange"
End
Simulation
  Max Output Level = 5
  Coordinate System = "Axi Symmetric"
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Transient
  Steady State Min Iterations = 1
  Steady State Max Iterations = 1
  Timestepping Method = Implicit Euler
  Timestep Sizes = 3.0e-2
  Timestep Intervals = 2
  Output Intervals = 1
End
Constants
  Gravity(4) = 0 -1 0 9.82
  Stefan Boltzmann = 5.67e-08
End
Body 1
  Name = "solid"
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 2
  Name = "melt"
  Equation = 1
  Material = 2
  Initial Condition = 1
End
Body 3
  Name = "heater"
  Equation = 2
  Material = 3
  Initial Condition = 1
  Body Force = 1
End
Body 4
  Name = "interface"
  Equation = 3
  Material = 1
  Initial Condition = 1
End
Equation 1
  Active Solvers(2) = 1 2
  Convection = Constant
End
Equation 2
  Active Solvers(1) = 2
  Convection = Constant
End
Equation 3
  Active Solvers(1) = 3
End
Solver 1
  Equation = "Mesh Update"
  Linear System Solver = "Direct"
  Linear System Direct Method = "umfpack"
  Steady State Convergence Tolerance = 1.0e-4
End
Solver 2
  Equation = "Heat Equation"
  Procedure = "HeatSolve" "HeatSolver"
  Linear System Solver = "Direct"
  Linear System Direct Method = "umfpack"
  Nonlinear System Convergence Tolerance = 1.0e-5
  Nonlinear System Max Iterations = 1
  Nonlinear System Relaxation Factor = 1.0
  Calculate Loads = True
  Steady State Convergence Tolerance = 1.0e-4
  Stabilize = Logical True
End
Solver 3
  Variable = PhaseSurface
  Equation = "Phase Surface"
  Procedure = "TransientPhaseChange" "TransientPhaseChange"
  Linear System Solver = "Direct"
  Linear System Direct Method = "umfpack"
  Nonlinear System Relaxation Factor = 1.0
  Steady State Convergence Tolerance = 1.0e-3
  Surface Smoothing Factor = Real 0.0
  Transient Speedup = Real 1.0
  Velocity Smoothing Factor = Real 0.05
  Velocity Relaxation Factor = Real 0.5
End
Solver 4
  Equation = String SaveScalars
  Procedure = File "SaveData" "SaveScalars"
  Filename = File "scalars.dat"
  Variable 1 = String PhaseSurface
  Operator 1 = String max
  Variable 2 = String PhaseSurface
  Operator 2 = String min
  Variable 3 = String Temperature
  Operator 3 = String mean
End
Body Force 1
End
Material 1
  Solid = Logical True
  Melting Point = Real 0.5
  Density = 1
  Heat Capacity = 1.0
  Heat Conductivity = 5.0
  Latent Heat = Real 1.0
  Convection Velocity 1 = 0.0
  Convection Velocity 2 = Real 1.0
End
Material 2
  Liquid = Logical True
  Melting Point = Real 0.5
  Density = 1
  Heat Capacity = 1.0
  Heat Conductivity = 5.0
End
Material 3
  Density = 1
  Heat Capacity = 1.0
  Heat Conductivity = 1.0
End
Initial Condition 1
  Temperature = Variable Coordinate 2
    Real
      0.0 1.0
      2.0 0.0
    End
  PhaseSurface = Real 0.0
  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.0
End
Boundary Condition 1
  Name = "melt_crystal"
  Target Boundaries = 1
  Temperature = Real 0.5
  Mesh Update 1 = 0
  Mesh Update 2 = Equals PhaseSurface
  Normal Target Body = Integer 1
  Body Id = Integer 4
End
Boundary Condition 2
  Name = "melt_heater"
  Target Boundaries = 2
  Mesh Update 1 = 0
  Mesh Update 2 = 0
  Temperature = 1.0
End
Boundary Condition 3
  Name = "solid_top"
  Target Boundaries = 3
  Mesh Update 1 = 0
  Mesh Update 2 = 0
  Temperature = 0.0
End
Boundary Condition 4
  Name = "right_walls"
  Target Boundaries = 4
  Mesh Update 1 = 0
  Heat Flux BC = Logical True
  Heat Transfer Coefficient = Real 1.0
  External Temperature = Real 1.0
  Phase Change Side = Logical True
End
Boundary Condition 5
  Name = "axis"
  Target Boundaries = 5
  Mesh Update 1 = 0
End
Boundary Condition 6
  Name = "heater_bottom"
  Target Boundaries = 6
End
EOF
# Config B: nodal-loads flux path + PDE-smoothed surface update path.
sed -e 's/Surface Smoothing Factor = Real 0.0/Surface Smoothing Factor = Real 0.05/' \
    -e 's/Velocity Smoothing Factor = Real 0.05/Velocity Smoothing Factor = Real 0.10\n  Use Nodal Loads = Logical True/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -q "Use Nodal Loads" "$WORK/caseB.sif" || { echo "FAIL: case B branch flags missing"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/TransientPhaseChange.so"
cp "$PORT_SO" "$WORK/shadow-port/TransientPhaseChange.so"
nm -D --defined-only "$WORK/shadow-port/TransientPhaseChange.so" | grep -q "T transientphasechange_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/PhaseChange" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER SOLVER.*STARTED AT|ELMER.*FINISHED AT" "$1" || true
}

echo "-- differential A (internal flux + simple update) + rerun (re-entry)"
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

echo "-- differential B (nodal loads + PDE smoothing) + branch liveness"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"
# Liveness: B must exercise different branches (different output), agreed by both arms.
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
if diff -q "$WORK/a0.txt" "$WORK/b1.txt" > /dev/null 2>&1; then
  echo "FAIL: branch-liveness control produced identical output"; exit 1;
fi
norm "$WORK/runA-orig/scalars.dat" > "$WORK/sa.txt"
norm "$WORK/runB-orig/scalars.dat" > "$WORK/sb.txt"
if diff -q "$WORK/sa.txt" "$WORK/sb.txt" > /dev/null 2>&1; then
  echo "FAIL: branch-liveness scalars identical"; exit 1;
fi
echo "   branch-liveness control differs, both arms agree"

echo "TRANSIENTPHASECHANGE PROOF GREEN ($((SECONDS - T0))s wall)"
