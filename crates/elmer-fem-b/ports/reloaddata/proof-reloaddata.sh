#!/bin/bash
# T021-B ReloadData port proof: the Rust cdylib replaces ReloadData.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `reloadsolution_`);
#  2. two distinct .sif differentials exercise the solver (config A:
#     Starting Position 5 / Range 1..7 / Intervals 1 / Continuous True,
#     i.e. the continuous-disabled non-continuous load arm with wrap;
#     config B: Starting Position 1 / Range 1..4 / Intervals 2, i.e. the
#     continuous load arm plus the MOD-skip branch) are byte-identical to
#     the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (SAVE re-entry through Round/Visit);
#  4. the two configs produce different output (branch liveness), and
#     both arms agree on each.
#
# Restart input is the fem/tests/reload fixture (Step.grd mesh + flow.dat);
# the .sif files below are clean-room (own solver list with a SaveScalars
# dump of the reloaded fields so the load is observable in scalars.dat).
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout only for the Step mesh input. Honest halts: any divergence is
# a FAIL. NEVER edit this script to make it pass — fix the port or the
# wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-reloaddata"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/Step*
echo "== T021-B ReloadData proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libreloaddata.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ReloadData.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ReloadData module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T reloadsolution_" "$WORK/port.syms" \
  || { echo "FAIL: reloadsolution_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's Step.grd + flow.dat fixture)"
cp "$ELMER_FORK/fem/tests/reload/Step.grd" "$ELMER_FORK/fem/tests/reload/flow.dat" \
   "$ELMER_FORK/fem/tests/reload/flow.dat.pos" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 Step > /dev/null)
cp "$WORK/flow.dat" "$WORK/flow.dat.pos" "$WORK/Step/"

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "Step"
End
Simulation
  Max Output Level = 5
  Coordinate System = Cartesian 2D
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Transient
  Steady State Max Iterations = 20
  Timestep Intervals = 6
  Timestep Sizes = 2.5
  Output Intervals = 1
End
Body 1
  Name = "Body1"
  Equation = 1
  Material = 1
End
Equation 1
  Name = "Equation1"
  Active Solvers(4) = 1 2 3 4
  Convection = Computed
End
Solver 1
  Exec Solver = Before Timestep
  Equation = restarting
  Procedure = "ReloadData" "ReloadSolution"
  Reload Solution File = "flow.dat"
  Reload Starting Position = Integer 5
  Reload Range Minimum = Integer 1
  Reload Range Maximum = Integer 7
  Reload Reading Intervals = Integer 1
  Continuous Reading = Logical True
End
Solver 2
  Exec Solver = Never
  Equation = "Navier-Stokes"
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Stabilization Method = String Stabilized
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU1
  Linear System Residual Output = 1
  Steady State Convergence Tolerance = 1.0e-5
  Nonlinear System Convergence Tolerance = 1.0e-5
  Nonlinear System Max Iterations = 1
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Newton After Tolerance = 1.0e-2
  Nonlinear System Relaxation Factor = 1.0
End
Solver 3
  Equation = "Mass transfer"
  Procedure = "AdvectionDiffusion" "AdvectionDiffusionSolver"
  Variable = "Species"
  Linear System Solver = Direct
  Linear System Direct Method = UMFPack
  Stabilization Method = Bubbles
  Steady State Convergence Tolerance = 1.0e-5
End
Solver 4
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Species
  Operator 1 = min
  Variable 2 = Species
  Operator 2 = max
  Variable 3 = Velocity 1
  Operator 3 = min
  Variable 4 = Velocity 1
  Operator 4 = max
End
Material 1
  Name = "Material1"
  Density = 1
  Viscosity = 0.01
  Species Diffusivity = 1.0e-4
End
Initial Condition 1
  Velocity 1 = 0.0
End
Boundary Condition 1
  Name = "Constraint1"
  Target Boundaries = 1
  Velocity 1 = 1.0
  Velocity 2 = 0
  Species = Real 0.1
Boundary Condition 2
  Name = "Constraint2"
  Target Boundaries = 2
  Pressure = 0
  Velocity 2 = 0
End
Boundary Condition 3
  Name = "Constraint3"
  Target Boundaries = 3
  Velocity 1 = 0
  Velocity 2 = 0
  Species = Real 0
End
EOF
# Config B: continuous arm (Starting Position 1 keeps Continuous Reading
# enabled), narrower range, and Intervals 2 (MOD-skip branch).
sed -e 's/Reload Starting Position = Integer 5/Reload Starting Position = Integer 1/' \
    -e 's/Reload Range Maximum = Integer 7/Reload Range Maximum = Integer 4/' \
    -e 's/Reload Reading Intervals = Integer 1/Reload Reading Intervals = Integer 2/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -q "Reload Range Maximum = Integer 4" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ReloadData.so"
cp "$PORT_SO" "$WORK/shadow-port/ReloadData.so"
nm -D --defined-only "$WORK/shadow-port/ReloadData.so" | grep -q "T reloadsolution_" \
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
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Time spent for restart" "$1" || true
}

echo "-- differential A (non-continuous arm + wrap) + rerun (re-entry)"
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

echo "-- differential B (continuous arm + MOD-skip branch)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (configs must differ; arms agree on each)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: configs A and B produce identical stdout (no branch delta)"; exit 1; }
diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: configs A and B produce identical scalars (no branch delta)"; exit 1; }
echo "   configs differ on both arms; differentials agree"

echo "RELOADDATA PROOF GREEN (wall: $((SECONDS - T0))s)"
