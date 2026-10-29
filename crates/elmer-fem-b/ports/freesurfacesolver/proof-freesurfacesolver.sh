#!/bin/bash
# T022-A FreeSurfaceSolver solver-entry proof: the Rust cdylib replaces
# FreeSurfaceSolver.so.
#
# Reference symbol surface (nm -D on FreeSurfaceSolver.so): exactly
#   T freesurfacesolver_
#   T freesurfacesolver_rotinit_
#   T freesurfacetomeshupdate_
#   T freesurfacetomeshupdate1_
#   T freesurfacetomeshupdate2_
#   T freesurfacetomeshupdate3_
# (the campaign ticket abbreviated the mesh-update names; the port
# matches the measured `nm` ground truth, proven below by empty diff).
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (nm -D equality over the six entry symbols above);
#  2. two distinct clean-room .sif differentials (config A: stabilized +
#     global Maximum Displacement + flux accumulation, 1 nonlinear
#     iteration; config B: bubbles + Apply Dirichlet limiters + local
#     maximum displacement + relaxation 0.7 + scalar accumulation +
#     RotInit entry + FreeSurfaceToMeshUpdate1/2 user functions, 2
#     nonlinear iterations; both transient, 2 timesteps) are
#     byte-identical between reference and port except wall-clock lines;
#  3. a second full process run of config A against the port is
#     byte-identical (idempotent re-entry through the wrapper-side SAVE
#     state);
#  4. a branch-liveness check: configs A and B produce different
#     FreeSurface max/min output, and both arms agree on each config.
#
# Mesh: the freesurf_maxd slab (10 x 0.5, Cartesian 2D), rebuilt from the
# Elmer checkout's Mesh.grd via ElmerGrid. Needs a configured+built Elmer
# tree (ELMER_BUILD, serial or campaign config — both arms share it, so
# the differential is self-consistent).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="/home/john/firstmate/data/rs-elmer-t022a/scratch/proof-freesurfacesolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/Mesh* "$WORK"/mesh
echo "== T022-A FreeSurfaceSolver proof: work=$WORK build=$ELMER_BUILD"
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
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_MODDIR"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libfreesurfacesolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FreeSurfaceSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FreeSurfaceSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in freesurfacesolver_ freesurfacesolver_rotinit_ freesurfacetomeshupdate_ \
    freesurfacetomeshupdate1_ freesurfacetomeshupdate2_ freesurfacetomeshupdate3_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not exported"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's freesurf_maxd Mesh.grd)"
cp "$ELMER_FORK/fem/tests/freesurf_maxd/Mesh.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 Mesh > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "Mesh"
  Check Keywords Warn
End
Constants
  Gravity(4) = 0 1 0 9.81
  Stefan Boltzmann = 5.67e-8
End
Simulation
  Max Output Level = 4
  Coordinate System = Cartesian 2D
  Simulation Type = Transient
  Steady State Max Iterations = 1
  Timestepping Method = BDF
  BDF Order = 1
  Timestep Sizes = Real 606825.5
  Timestep Intervals = 2
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 2
  Equation = 2
  Material = 1
  Body Force = 1
  Initial Condition = 1
End
Equation 1
  Active Solvers(1) = 2
End
Equation 2
  Active Solvers(1) = 1
  Convection = "constant"
End
Material 1
  Viscosity = 1.0E12
  Density = 918.0
  Convection Velocity 1 = Variable Coordinate
    Real MATC "(tx(0) - 5.0)*tx(1)/(31556926.0)"
  Convection Velocity 2 = Variable Coordinate
    Real MATC "((5.0 - abs(tx(0) - 5.0))*tx(1)/10.0)/(31556926.0)"
  Youngs Modulus = 1
  Poisson Ratio = 0.3
End
Body Force 1
  FreeSurface Accumulation Flux 1 = Real 0.0
  FreeSurface Accumulation Flux 2 = Variable Coordinate 1
    Real MATC "1.0*(4.0 - sqrt((tx-5.0)*(tx-5.0)))/(31556926.0)"
End
Initial Condition 1
  Mesh Update 1 = 0.0
  Mesh Update 2 = 0.0
  FreeSurface = Equals Coordinate 2
End
Solver 1
  Equation = Free
  Procedure = "FreeSurfaceSolver" "FreeSurfaceSolver"
  Variable = FreeSurface
  Variable DOFs = 1
  Stabilization Method = Stabilized
  Linear System Solver = Direct
  Linear System Convergence Tolerance = 1.0e-7
  Nonlinear System Max Iterations = 1
  Nonlinear System Convergence Tolerance = 1.0e-6
  Steady State Convergence Tolerance = 1.0e-6
  Exported Variable 1 = Freesurface Residual
  Exported Variable 1 DOFs = 1
  Maximum Displacement = Real 0.01
End
Solver 2
  Equation = Mesh Update
  Linear System Solver = Direct
  Steady State Convergence Tolerance = 1.0e-6
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Exec Solver = After Timestep
  Filename = scalars.dat
  Variable 1 = FreeSurface
  Operator 1 = max
  Variable 2 = FreeSurface
  Operator 2 = min
End
Boundary Condition 1
  Target Boundaries = 1
  Name = "Bed"
  Mesh Update 1 = 0
  Mesh Update 2 = 0
End
Boundary Condition 2
  Target Boundaries = 3
  Name = "Free"
  Body Id = 2
  Mesh Update 2 = Variable FreeSurface
    Real MATC "tx - 0.4"
  Mesh Update 1 = 0
End
Boundary Condition 3
  Target Boundaries(2) = 2 4
  Name = "sides"
  Mesh Update 1 = 0
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "Mesh"
  Check Keywords Warn
End
Constants
  Gravity(4) = 0 1 0 9.81
  Stefan Boltzmann = 5.67e-8
End
Simulation
  Max Output Level = 4
  Coordinate System = Cartesian 2D
  Simulation Type = Transient
  Steady State Max Iterations = 1
  Timestepping Method = BDF
  BDF Order = 1
  Timestep Sizes = Real 606825.5
  Timestep Intervals = 2
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 2
  Equation = 2
  Material = 1
  Body Force = 1
  Initial Condition = 1
End
Equation 1
  Active Solvers(1) = 2
End
Equation 2
  Active Solvers(2) = 1 4
  Convection = "constant"
End
Material 1
  Viscosity = 1.0E12
  Density = 918.0
  Convection Velocity 1 = Variable Coordinate
    Real MATC "(tx(0) - 5.0)*tx(1)/(31556926.0)"
  Convection Velocity 2 = Variable Coordinate
    Real MATC "((5.0 - abs(tx(0) - 5.0))*tx(1)/10.0)/(31556926.0)"
  Youngs Modulus = 1
  Poisson Ratio = 0.3
  Min FreeSurface = Real 0.0
  Max FreeSurface = Real 0.002
End
Body Force 1
  FreeSurface Accumulation = Real 1.0e-7
  Max Displacement FreeSurface = Real 0.005
End
Initial Condition 1
  Mesh Update 1 = 0.0
  Mesh Update 2 = 0.0
  FreeSurface = Equals Coordinate 2
End
Solver 1
  Equation = Free
  Procedure = "FreeSurfaceSolver" "FreeSurfaceSolver"
  Variable = FreeSurface
  Variable DOFs = 1
  Linear System Solver = Direct
  Linear System Convergence Tolerance = 1.0e-7
  Nonlinear System Max Iterations = 2
  Nonlinear System Convergence Tolerance = 1.0e-6
  Steady State Convergence Tolerance = 1.0e-6
  Exported Variable 1 = Freesurface Residual
  Exported Variable 1 DOFs = 1
  Apply Dirichlet = Logical True
  Compute Local Maximum Displacement = Logical True
  Relaxation Factor = Real 0.7
End
Solver 2
  Equation = Mesh Update
  Linear System Solver = Direct
  Steady State Convergence Tolerance = 1.0e-6
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Exec Solver = After Timestep
  Filename = scalars.dat
  Variable 1 = FreeSurface
  Operator 1 = max
  Variable 2 = FreeSurface
  Operator 2 = min
End
Solver 4
  Equation = Free RotInit
  Procedure = "FreeSurfaceSolver" "FreeSurfaceSolver_RotInit"
  Exec Solver = Before Simulation
  Variable = FSRef
  Variable DOFs = 1
  Free Surface Variable Name = String "FreeSurface"
  FreeSurface Orientation(3) = Real 0.0 0.0 1.0
End
Boundary Condition 1
  Target Boundaries = 1
  Name = "Bed"
  Mesh Update 1 = 0
  Mesh Update 2 = 0
End
Boundary Condition 2
  Target Boundaries = 3
  Name = "Free"
  Body Id = 2
  Mesh Update 2 = Variable FreeSurface, FSRef
    Real Procedure "FreeSurfaceSolver" "FreeSurfaceToMeshUpdate2"
  Mesh Update 1 = Variable FreeSurface, FSRef
    Real Procedure "FreeSurfaceSolver" "FreeSurfaceToMeshUpdate1"
End
Boundary Condition 3
  Target Boundaries(2) = 2 4
  Name = "sides"
  Mesh Update 1 = 0
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
rm -f "$WORK/shadow-port/FreeSurfaceSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/FreeSurfaceSolver.so"
nm -D --defined-only "$WORK/shadow-port/FreeSurfaceSolver.so" | grep -q "T freesurfacesolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/Mesh" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines; mask per-iteration (s) timing values
  # (structure and iteration counts stay compared)
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" | \
    sed -E -e 's/(Assembly: \(s\))( +[0-9.]+)( +[0-9.]+)/\1 T T/' \
           -e 's/(Solve: +\(s\))( +[0-9.]+)( +[0-9.]+)/\1 T T/' || true
}

echo "-- differential A (stabilized + global max displacement) + rerun (re-entry)"
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

echo "-- differential B (bubbles + Dirichlet limiters + RotInit + USFs)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness: A and B exercise different branches"
if diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: liveness control produced identical scalars; branches not distinguished"
  exit 1
fi
echo "   configs A/B differ as expected (both arms agree)"

echo "FREESURFACESOLVER PROOF GREEN (wall $((SECONDS - T0))s)"
