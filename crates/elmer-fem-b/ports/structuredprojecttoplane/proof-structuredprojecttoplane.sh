#!/bin/bash
# T021-B StructuredProjectToPlane solver-entry proof: the Rust cdylib replaces
# StructuredProjectToPlane.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `structuredprojecttoplane_`,
#     `structuredprojecttoplane_init_`);
#  2. two distinct .sif differentials (config A: geometry-only operators
#     depth/height/distance/index coupled back into the heat source via MATC;
#     config B: Temperature-dependent operators bottom/top/sum/max/isosurface
#     over 2 steady iterations, i.e. the repeated-invocation run) are
#     byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state).
#  4. branch-liveness control (config C: different operator branches) gives
#     different output, and both arms agree on it.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both arms
# share it, so the differential is self-consistent) and the Elmer checkout
# only for the cube.grd mesh input (upstream structured-mesh test mesh).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="${RS_ELMER_T021B_SCRATCH:-/home/john/firstmate/data/rs-elmer-t021b/scratch/proof-structuredprojecttoplane}"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/cube*
echo "== T021-B StructuredProjectToPlane proof: work=$WORK build=$ELMER_BUILD"
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
  *) echo "FAIL: need rustc nightly 1.100.0, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libstructuredprojecttoplane.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/StructuredProjectToPlane.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T structuredprojecttoplane_" "$WORK/port.syms" \
  || { echo "FAIL: structuredprojecttoplane_ not the exported entry"; exit 1; }
grep -qx "T structuredprojecttoplane_init_" "$WORK/port.syms" \
  || { echo "FAIL: structuredprojecttoplane_init_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's StructuredHeightDepth cube.grd)"
cp "$ELMER_FORK/fem/tests/StructuredHeightDepth/cube.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 cube > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "cube"
End
Simulation
  Coordinate System = Cartesian
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Max Output Level = 5
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
  Heat Conductivity = 1.0
End
Body Force 1
  Heat Source = Variable depth_f, height_f
    Real MATC "tx(0)+tx(1)"
End
Boundary Condition 1
  Target Boundaries(4) = 1 2 3 4
End
Boundary Condition 2
  Target Boundaries(1) = 5
  Temperature = Real 0.0
End
Boundary Condition 3
  Target Boundaries(1) = 6
  Temperature = Real 0.0
End
Solver 1
  Equation = HeightDepth
  Procedure = "StructuredProjectToPlane" "StructuredProjectToPlane"
  Active Coordinate = Integer 3
  Operator 1 = depth
  Target Variable 1 = String depth_f
  Operator 2 = height
  Target Variable 2 = String height_f
  Operator 3 = distance
  Target Variable 3 = String dist_f
  Operator 4 = index
  Target Variable 4 = String index_f
End
Solver 2
  Equation = Heat Equation
  Variable = String Temperature
  Variable DOFs = 1
  Procedure = "HeatSolve" "HeatSolver"
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = none
  Linear System Residual Output = 1
  Steady State Convergence Tolerance = 1.0e-5
  Nonlinear System Convergence Tolerance = 1.0e-5
  Nonlinear System Max Iterations = 1
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Newton After Tolerance = 1.0e-2
  Nonlinear System Relaxation Factor = 1.0
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = max
  Variable 2 = depth_f
  Operator 2 = max
  Variable 3 = height_f
  Operator 3 = max
  Variable 4 = dist_f
  Operator 4 = max
  Variable 5 = index_f
  Operator 5 = max
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "cube"
End
Simulation
  Coordinate System = Cartesian
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Steady State
  Steady State Max Iterations = 2
  Output Intervals = 1
  Max Output Level = 5
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
  Heat Conductivity = 1.0
End
Body Force 1
  Heat Source = Real 10.0
End
Boundary Condition 1
  Target Boundaries(4) = 1 2 3 4
End
Boundary Condition 2
  Target Boundaries(1) = 5
  Temperature = Real 1.0
End
Boundary Condition 3
  Target Boundaries(1) = 6
  Temperature = Real 0.0
End
Solver 1
  Equation = Heat Equation
  Variable = String Temperature
  Variable DOFs = 1
  Procedure = "HeatSolve" "HeatSolver"
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = none
  Linear System Residual Output = 1
  Steady State Convergence Tolerance = 1.0e-5
  Nonlinear System Convergence Tolerance = 1.0e-5
  Nonlinear System Max Iterations = 1
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Newton After Tolerance = 1.0e-2
  Nonlinear System Relaxation Factor = 1.0
End
Solver 2
  Equation = HeightDepth
  Procedure = "StructuredProjectToPlane" "StructuredProjectToPlane"
  Active Coordinate = Integer 3
  Variable 1 = Temperature
  Operator 1 = bottom
  Target Variable 1 = String bot_temp
  Variable 2 = Temperature
  Operator 2 = top
  Target Variable 2 = String top_temp
  Target Variable 2 At Bottom = Logical True
  Variable 3 = Temperature
  Operator 3 = sum
  Target Variable 3 = String sum_temp
  Variable 4 = Temperature
  Operator 4 = max
  Target Variable 4 = String max_temp
  Variable 5 = Temperature
  Operator 5 = isosurface
  Isosurface Value 5 = Real 0.5
  Target Variable 5 = String iso_temp
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = max
  Variable 2 = bot_temp
  Operator 2 = max
  Variable 3 = top_temp
  Operator 3 = min
  Variable 4 = sum_temp
  Operator 4 = max
  Variable 5 = max_temp
  Operator 5 = min
  Variable 6 = iso_temp
  Operator 6 = max
End
EOF
# Config C: branch-liveness control — different operator branches, must differ.
sed -e 's/Operator 1 = bottom/Operator 1 = top/' \
    -e 's/Operator 3 = sum/Operator 3 = min/' \
    -e 's/Target Variable 3 = String sum_temp/Target Variable 3 = String min_temp/' \
    -e 's/Variable 4 = sum_temp/Variable 4 = min_temp/' \
    "$WORK/caseB.sif" > "$WORK/caseC.sif"
grep -q "min_temp" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/StructuredProjectToPlane.so"
cp "$PORT_SO" "$WORK/shadow-port/StructuredProjectToPlane.so"
nm -D --defined-only "$WORK/shadow-port/StructuredProjectToPlane.so" | grep -q "T structuredprojecttoplane_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/cube" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; tail -30 "$WORK/$run/stdout.txt"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER SOLVER.*STARTED AT|ELMER.*FINISHED AT|STARTED AT|FINISHED AT|Projection time" "$1" || true
}

echo "-- differential A (geometry operators + MATC coupling) + rerun (re-entry)"
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

echo "-- differential B (variable operators, 2 steady iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head -40; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control C (different operator branches)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head -40; exit 1; }
done
echo "   config C byte-identical (both arms agree)"
if diff "$WORK/runB-port/scalars.dat" "$WORK/runC-port/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: branch-liveness control did not change the output"
  exit 1
fi
echo "   control output differs across branches, as required"

T1=$SECONDS
echo "STRUCTUREDPROJECTTOPLANE PROOF GREEN ($((T1-T0))s wall)"
