#!/bin/bash
# T021-B StructuredMeshMapper solver-entry proof: the Rust cdylib replaces
# StructuredMeshMapper.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `structuredmeshmapper_` and
#     `structuredmeshmapper_init_`);
#  2. two distinct .sif differentials that exercise the solver on different
#     ValueList branches (config A: binary layer mapping on a 2-D
#     structured box mesh with 2 steady iterations, i.e. the
#     repeated-invocation run; config B: multilayer mapping with fixed
#     layers on a 3-D extruded cube mesh) are byte-identical to the
#     reference except wall-clock lines;
#  3. a second full process run of config A against the port is
#     byte-identical (SAVE re-entry);
#  4. branch-liveness control: config B takes the MultiLayerMapper branch
#     (vs config A's BinaryLayerMapper) and produces different output,
#     and both arms agree on it.
#
# Wall-clock lines (`Elapsed REAL time`, `SOLVER TOTAL TIME`,
# `ELMER SOLVER ... STARTED AT`, `ELMER ... FINISHED AT`, the solver's own
# `Active coordinate mapping time` CPU-time line, and
# DetectExtrudedStructure's `Top and bottom pointer init time` CPU-time
# line) are wall-clock by construction and are normalized away; everything
# else must match byte for byte.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout for the box/cube mesh inputs.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-structuredmeshmapper"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/box* "$WORK"/cube*
echo "== StructuredMeshMapper proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*"feaadeeac 2026-09-19"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly (feaadeeac 2026-09-19)"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libstructuredmeshmapper.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/StructuredMeshMapper.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the StructuredMeshMapper module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T structuredmeshmapper_" "$WORK/port.syms" \
  || { echo "FAIL: structuredmeshmapper_ not exported"; exit 1; }
grep -qx "T structuredmeshmapper_init_" "$WORK/port.syms" \
  || { echo "FAIL: structuredmeshmapper_init_ not exported"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- meshes (ElmerGrid structured box + cube from the Elmer checkout)"
cp "$ELMER_FORK/fem/tests/structmap/box.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 box > /dev/null)
cp "$ELMER_FORK/fem/tests/structmap_multilayer/cube.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 cube > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "box"
End
Simulation
  Max Output Level = 5
  Coordinate System = "Cartesian 2D"
  Simulation Type = "Steady State"
  Steady State Max Iterations = 2
  Output Intervals = 1
End
Body 1
  Name = "Body1"
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Name = "Material1"
  Density = 1.0
  Heat Conductivity = 1.0
End
Body Force 1
  Name = "Heating"
  Heat Source = Real 1.0
End
Boundary Condition 1
  Name = "bottom"
  Target Boundaries = 1
  Bottom Surface = Variable "Coordinate 1"
    Real MATC "0.1*cos(5*tx)"
  Temperature = Real 0.0
End
Boundary Condition 2
  Name = "right"
  Target Boundaries = 2
  Temperature = Real 0.0
End
Boundary Condition 3
  Name = "top"
  Target Boundaries = 3
  Top Surface = Variable "Coordinate 1"
    Real MATC "1+0.1*sin(5*tx)"
  Temperature = Real 0.0
End
Boundary Condition 4
  Name = "left"
  Target Boundaries = 4
  Temperature = Real 0.0
End
Solver 1
  Equation = MapCoordinate
  Procedure = "StructuredMeshMapper" "StructuredMeshMapper"
  Active Coordinate = Integer 2
End
Solver 2
  Equation = "Heat Equation"
  Procedure = "HeatSolve" "HeatSolver"
  Variable = "Temperature"
  Variable Dofs = 1
  Linear System Solver = "Direct"
  Linear System Direct Method = "Umfpack"
  Steady State Convergence Tolerance = 1.0e-05
  Nonlinear System Convergence Tolerance = 1.0e-05
  Nonlinear System Max Iterations = 1
  Stabilize = True
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = min
  Variable 2 = Temperature
  Operator 2 = max
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "cube"
End
Simulation
  Max Output Level = 6
  Coordinate System = "Cartesian"
  Simulation Type = "Steady State"
  Steady State Max Iterations = 1
  Output Intervals = 0
End
Body 1
  Name = "Body1"
  Equation = 1
  Material = 1
End
Body 2
  Name = "Surface1"
  Body Force = 1
  Equation = 2
  Material = 1
  Initial Condition = 1
End
Equation 1
  Name = "Equation1"
  Active Solvers(1) = 1
End
Equation 2
  Name = "SurfaceHeat"
  Active Solvers(2) = 2 3
End
Initial Condition 1
  Name = "LayerInit"
  LayerHeight 1 = Variable "Coordinate"
    Real MATC "-0.1"
  LayerHeight 2 = Variable "Coordinate"
    Real MATC "0.2-0.1*tx(0)"
  LayerHeight 3 = Variable "Coordinate"
    Real MATC "0.4+0.1*tx(0)"
  LayerHeight 4 = Variable "Coordinate"
    Real MATC "0.6+0.1*tx(1)"
End
Solver 1
  Equation = "MapCoordinate"
  Procedure = "StructuredMeshMapper" "StructuredMeshMapper"
  Active Coordinate = Integer 3
  Fixed Layer Indexes(4) = Integer 1 6 15 21
  Fixed Layer Variable = String "LayerHeight"
End
Solver 2
  Equation = "Heat Equation"
  Procedure = "HeatSolve" "HeatSolver"
  Variable = "Temperature"
  Linear System Solver = "Direct"
  Linear System Direct Method = "Umfpack"
  Nonlinear System Max Iterations = 1
  Exported Variable 1 = -dofs 4 LayerHeight
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = min
  Variable 2 = Temperature
  Operator 2 = max
End
Material 1
  Name = "Material1"
  Density = 1
  Heat Conductivity = 1
End
Body Force 1
  Name = "Heating"
  Heat Source = 1
End
Boundary Condition 1
  Name = "Sides"
  Target Boundaries(4) = 1 2 3 4
  Temperature = 0.0
End
Boundary Condition 2
  Name = "Bottom"
  Target Boundaries = 5
End
Boundary Condition 3
  Name = "Top"
  Target Boundaries = 6
  Body Id = 2
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
rm -f "$WORK/shadow-port/StructuredMeshMapper.so"
cp "$PORT_SO" "$WORK/shadow-port/StructuredMeshMapper.so"
nm -D --defined-only "$WORK/shadow-port/StructuredMeshMapper.so" | grep -q "T structuredmeshmapper_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir> <mesh>
  local case=$1 arm=$2 run=$3 mesh=$4
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/$mesh" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|Active coordinate mapping time|Top and bottom pointer init time" "$1" || true
}

echo "-- differential A (binary mapping, 2 steady iterations) + rerun (re-entry)"
run_case caseA orig runA-orig box
run_case caseA port runA-port box
run_case caseA port runA-port2 box
for f in stderr.txt scalars.dat; do
  cmp "$WORK/runA-orig/$f" "$WORK/runA-port/$f" \
    || { echo "FAIL: config A $f diverges"; exit 1; }
  cmp "$WORK/runA-port/$f" "$WORK/runA-port2/$f" \
    || { echo "FAIL: config A rerun $f diverges"; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  if [ "$f" = stdout.txt ]; then
    norm "$WORK/runA-port/stdout.txt" > "$WORK/ar1.txt"
    norm "$WORK/runA-port2/stdout.txt" > "$WORK/ar2.txt"
    diff "$WORK/ar1.txt" "$WORK/ar2.txt" > /dev/null 2>&1 || \
      { echo "FAIL: config A rerun stdout diverges"; diff "$WORK/ar1.txt" "$WORK/ar2.txt" | head; exit 1; }
  else
    cmp "$WORK/runA-port/$f" "$WORK/runA-port2/$f" \
      || { echo "FAIL: config A rerun $f diverges"; exit 1; }
  fi
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (multilayer mapping, fixed layers)"
run_case caseB orig runB-orig cube
run_case caseB port runB-port cube
for f in stdout.txt stderr.txt scalars.dat; do
  if [ "$f" = stdout.txt ]; then
    norm "$WORK/runB-orig/stdout.txt" > "$WORK/b0.txt"
    norm "$WORK/runB-port/stdout.txt" > "$WORK/b1.txt"
    diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
      { echo "FAIL: config B stdout diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head -40; exit 1; }
  else
    cmp "$WORK/runB-orig/$f" "$WORK/runB-port/$f" \
      || { echo "FAIL: config B $f diverges"; exit 1; }
  fi
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B takes the multilayer branch, output differs; arms agree)"
grep -q "Mapping using" "$WORK/runB-port/stdout.txt" \
  || { echo "FAIL: multilayer branch marker missing in B"; exit 1; }
grep -q "Mapping using" "$WORK/runA-port/stdout.txt" \
  && { echo "FAIL: multilayer branch marker unexpectedly in A"; exit 1; }
if cmp -s "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat"; then
  echo "FAIL: branch-liveness control produced identical scalars"
  exit 1
fi
echo "   B differs from A (multilayer branch live), arms agree on both"

echo "-- wall seconds: $((SECONDS - T0))"
echo "STRUCTUREDMESHMAPPER PROOF GREEN"
