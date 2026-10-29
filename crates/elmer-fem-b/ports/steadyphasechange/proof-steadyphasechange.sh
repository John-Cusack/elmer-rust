#!/bin/bash
# T021-B SteadyPhaseChange solver-plugin proof: the Rust cdylib replaces
# SteadyPhaseChange.so (both entries: steadyphasechange_ + meltingheat_).
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: steadyphasechange_ + meltingheat_);
#  2. two distinct clean-room .sif differentials that exercise the solver
#     (config A: isotherm-mapping branch; config B: Newton + Triple Point
#     Fixed + Internal Mesh Movement branches, 3 steady iterations each)
#     are byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch liveness: A and B produce different output, and both arms
#     agree on each config.
#
# Mesh note: the PoissonDB triangular square.grd mesh cannot exercise
# this solver (no coupled heat field, no Phase Change boundary, wrong
# body layout), so the proof uses the upstream PhaseChange2 structured
# mesh (fem/tests/PhaseChange2/PhaseChange.grd) through the same
# ElmerGrid 1 2 flow — quadrilateral bulk elements (code 404) and line
# boundary elements (code 202), exactly what CreateIsotherm needs.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout only for the PhaseChange.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-steadyphasechange"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/PhaseChange*
echo "== T021-B SteadyPhaseChange proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libsteadyphasechange.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/SteadyPhaseChange.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the SteadyPhaseChange module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T steadyphasechange_" "$WORK/port.syms" \
  || { echo "FAIL: steadyphasechange_ not the exported entry"; exit 1; }
grep -qx "T meltingheat_" "$WORK/port.syms" \
  || { echo "FAIL: meltingheat_ not the exported co-entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid structured PhaseChange mesh from the Elmer checkout)"
cp "$ELMER_FORK/fem/tests/PhaseChange2/PhaseChange.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 PhaseChange > /dev/null)

echo "-- configs (clean-room .sif files below, adapted from"
echo "   fem/tests/PhaseChange2/PhaseChangeSS.sif: axisymmetric coupled"
echo "   heat + mesh-update + phase-surface solve)"
cat > "$WORK/caseA.sif" <<'EOF'
$ velo = 1.0
$ melt = 0.5
$ trans = 1.0
$ cap = 1.0
$ lat = 1.0
$ cond = 5.0
Header
  Mesh DB "." "PhaseChange"
End
Simulation
  Max Output Level = 5
  Coordinate System = "Axi Symmetric"
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Steady State
  Steady State Min Iterations = 3
  Steady State Max Iterations = 3
  Output Intervals = 0
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
  Nonlinear System Newton After Iterations = 0
  Nonlinear System Newton After Tolerance = 1.0e-2
  Steady State Convergence Tolerance = 1.0e-4
  Stabilize = Logical True
End
Solver 3
  Variable = PhaseSurface
  Equation = "Phase Surface"
  Procedure = "SteadyPhaseChange" "SteadyPhaseChange"
  Nonlinear System Relaxation Factor = 1.0
  Nonlinear System Newton After Iterations = 10
  Steady State Convergence Tolerance = 1.0e-3
End
Solver 4
  Equation = String SaveScalars
  Procedure = File "SaveData" "SaveScalars"
  Filename = File scalars.dat
  Variable 1 = String PhaseSurface
  Operator 1 = String max
  Variable 2 = String PhaseSurface
  Operator 2 = String min
  Variable 3 = String Temperature
  Operator 3 = String mean
  Variable 4 = String time
End
Body Force 1
End
Material 1
  Solid = Logical True
  Melting Point = Real $ melt
  Density = 1
  Heat Capacity = $ cap
  Heat Conductivity = $ cond
  Youngs Modulus = 1.0
  Poisson Ratio = 0.3
  Latent Heat = Real $ lat
  Convection Velocity 1 = 0.0
  Convection Velocity 2 = Real $ velo
End
Material 2
  Liquid = Logical True
  Melting Point = Real $ melt
  Density = 1
  Heat Capacity = $ cap
  Heat Conductivity = $ cond
  Youngs Modulus = 1.0
  Poisson Ratio = 0.3
  Convection Velocity 1 = 0.0
  Convection Velocity 2 = Real $ velo
End
Material 3
  Density = 1
  Heat Capacity = $ cap
  Heat Conductivity = 1.0
End
Initial Condition 1
  Temperature = Real 1.0
  PhaseSurface = Real 0.0
  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.0
End
Boundary Condition 1
  Name = "melt_crystal"
  Target Boundaries = 1
  Heat Flux BC = Logical True
  Phase Change = Logical True
  Heat Flux = Variable Coordinate 1
    Real Procedure "SteadyPhaseChange" "MeltingHeat"
  Mesh Update 1 = 0
  Mesh Update 2 = Equals PhaseSurface
  Latent Heat = Real $ lat
  Melting Point = Real $ melt
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
  Heat Transfer Coefficient = Real $ trans
  External Temperature = Real 1.0
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
# Config B: same physics, different ValueList branches — Newton
# formulation from iteration 2, Triple Point Fixed melting point, and
# internal mesh movement. Must produce different output than A.
cp "$WORK/caseA.sif" "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
old = """  Nonlinear System Newton After Iterations = 10
  Steady State Convergence Tolerance = 1.0e-3"""
new = """  Nonlinear System Newton After Iterations = 1
  Triple Point Fixed = Logical True
  Internal Mesh Movement = Logical True
  Steady State Convergence Tolerance = 1.0e-3"""
assert old in t, "Solver 3 block not found"
open(p, "w").write(t.replace(old, new))
PYEOF
grep -c "Triple Point Fixed" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/SteadyPhaseChange.so"
cp "$PORT_SO" "$WORK/shadow-port/SteadyPhaseChange.so"
nm -D --defined-only "$WORK/shadow-port/SteadyPhaseChange.so" | grep -q "T steadyphasechange_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/SteadyPhaseChange.so" | grep -q "T meltingheat_" \
  || { echo "FAIL: staged port .so lost its co-entry"; exit 1; }

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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER.*STARTED AT|ELMER.*FINISHED AT|STARTED AT|FINISHED AT" "$1" || true
}

echo "-- differential A (isotherm branch, 3 iterations) + rerun (re-entry)"
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

echo "-- differential B (Newton + Triple Point Fixed + mesh movement)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (B takes different branches, different output)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
if cmp -s "$WORK/a1.txt" "$WORK/b2.txt"; then
  echo "FAIL: config B output identical to A (branches not exercised)"
  exit 1
fi
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b1.txt"
if cmp -s "$WORK/a0.txt" "$WORK/b1.txt"; then
  echo "FAIL: reference B output identical to A (bad liveness control)"
  exit 1
fi
grep -q "Steady state newton formulation" "$WORK/runB-port/stdout.txt" \
  || { echo "FAIL: Newton branch not taken in B"; exit 1; }
grep -q "Melting point set" "$WORK/runB-port/stdout.txt" \
  || { echo "FAIL: Triple Point Fixed branch not taken in B"; exit 1; }
grep -q "Steady state isoterm formulation" "$WORK/runA-port/stdout.txt" \
  || { echo "FAIL: isotherm branch not taken in A"; exit 1; }
echo "   liveness control holds (A != B, both arms agree)"

echo "PROOF STEADYPHASECHANGE GREEN ($(($SECONDS - $T0))s wall)"
