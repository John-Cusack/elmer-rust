#!/bin/bash
# T021-B RigidMeshMapper solver-entry proof: the Rust cdylib replaces
# RigidMeshMapper.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `rigidmeshmapper_`);
#  2. four distinct .sif differentials are byte-identical to the
#     reference except wall-clock lines —
#     A: body-force-gated translate + per-node Mesh Relax, Constant
#        Mapping, no matrix (GotTranslate/AnyRelax path);
#     B: relaxation via Laplace (matrix assembly, Moving/Fixed BCs,
#        normalize, cumulative aliasing, rotate+scale+translate+origin,
#        rotation-axis-order, mesh velocity);
#     C: Whole Mesh Mode from the solver section (GotMap latch);
#     D: Rotor Mode with Rotor Radius (rotor marking, partial mapping);
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch-liveness: A/B/C/D outputs mutually differ while both arms
#     agree on each.
#
# Defined-input notes (upstream quirks, contained, documented in the
# port): `Mesh Scale` is only ever combined with rotation (a lone scale
# reads the original's uninitialized RotMatrix local); `Mesh Matrix`
# arrays are not used (the original recomputes GotMatrix afterwards,
# making a lone array a defined no-op); per-component keywords are
# always fully specified, as in upstream tests.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-rigidmeshmapper"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== RigidMeshMapper proof: work=$WORK build=$ELMER_BUILD"
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
  *"1.100.0-nightly"*"feaadeeac"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly (feaadeeac 2026-09-19)"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/librigidmeshmapper.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/RigidMeshMapper.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the RigidMeshMapper module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T rigidmeshmapper_" "$WORK/port.syms" \
  || { echo "FAIL: rigidmeshmapper_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- clippy clean"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release 2>&1 | tail -5)
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release -- -D warnings 2>&1 | tail -3) \
  || { echo "FAIL: clippy warnings"; exit 1; }

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

HEAT_SOLVER='Solver 2
  Equation = HeatSolver
  Variable = Temperature
  Procedure = "HeatSolve" "HeatSolver"
  Steady State Convergence Tolerance = 1.0e-5
  Nonlinear System Convergence Tolerance = 1.0e-4
  Nonlinear System Max Iterations = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = min
  Variable 2 = Temperature
  Operator 2 = max
End'

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<EOF
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 2
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 2
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 3
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 4
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Density = 1.0
  Heat Conductivity = 1.0
End
Body Force 1
  Mesh Translate 1 = Real 0.1
  Mesh Relax = Real 0.5
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Temperature = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Temperature = Real 1.0
End
Solver 1
  Equation = MeshDeform
  Procedure = "RigidMeshMapper" "RigidMeshMapper"
  Constant Mapping = Logical True
End
$HEAT_SOLVER
EOF

cat > "$WORK/caseB.sif" <<EOF
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 2
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 2
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 3
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 4
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Density = 1.0
  Heat Conductivity = 1.0
End
Body Force 1
  Mesh Rotate 3 = Real 30.0
  Mesh Origin(3) = Real 0.5 0.5 0.0
  Mesh Translate 1 = Real 0.05
  Mesh Scale 1 = Real 1.1
  Heat Source = Real 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Fixed Boundary = Logical True
  Temperature = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Moving Boundary = Logical True
End
Solver 1
  Equation = MeshDeform
  Procedure = "RigidMeshMapper" "RigidMeshMapper"
  Cumulative Displacements = Logical True
  Calculate Mesh Velocity = Logical True
  Mesh Relax Normalize = Logical True
  Mesh Rotation Axis Order(2) = Integer 2 1
  Variable = MeshRelax
  Nonlinear System Relaxation Factor = Real 0.5
  Nonlinear System Relaxation After = Integer 1
  Nonlinear System Max Iterations = Integer 30
  Nonlinear System Convergence Tolerance = Real 1.4e-4
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-12
  Linear System Preconditioning = ILU1
  Linear System ILUT Tolerance = 1.0e-3
  Linear System Abort Not Converged = False
  Linear System Residual Output = 20
  Linear System Precondition Recompute = 1
End
$HEAT_SOLVER
EOF

cat > "$WORK/caseC.sif" <<EOF
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
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
  Active Solvers(2) = 1 2
End
Material 1
  Density = 1.0
  Heat Conductivity = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Temperature = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Temperature = Real 1.0
End
Solver 1
  Equation = MeshDeform
  Procedure = "RigidMeshMapper" "RigidMeshMapper"
  Whole Mesh Mode = Logical True
  Mesh Translate 2 = Real -0.05
End
$HEAT_SOLVER
EOF

cat > "$WORK/caseD.sif" <<EOF
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 2
  Output Intervals = 1
  Rotor Radius = Real 0.4
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
  Active Solvers(2) = 1 2
End
Material 1
  Density = 1.0
  Heat Conductivity = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Temperature = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Temperature = Real 1.0
End
Solver 1
  Equation = MeshDeform
  Procedure = "RigidMeshMapper" "RigidMeshMapper"
  Rotor Mode = Logical True
  Mesh Rotate 3 = Real 15.0
End
$HEAT_SOLVER
EOF

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/RigidMeshMapper.so"
cp "$PORT_SO" "$WORK/shadow-port/RigidMeshMapper.so"
nm -D --defined-only "$WORK/shadow-port/RigidMeshMapper.so" | grep -q "T rigidmeshmapper_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER.*FINISHED AT|STARTED AT|Coordinate mapping time" "$1" || true
}

for cfg in caseA caseB caseC caseD; do
  echo "-- differential $cfg"
  run_case "$cfg" orig "run-${cfg}-orig"
  run_case "$cfg" port "run-${cfg}-port"
  for f in stdout.txt stderr.txt scalars.dat; do
    if ! cmp -s "$WORK/run-${cfg}-orig/$f" "$WORK/run-${cfg}-port/$f"; then
      if [ "$f" = "stdout.txt" ]; then
        norm "$WORK/run-${cfg}-orig/$f" > "$WORK/n0.txt"
        norm "$WORK/run-${cfg}-port/$f" > "$WORK/n1.txt"
        diff "$WORK/n0.txt" "$WORK/n1.txt" > /dev/null 2>&1 \
          || { echo "FAIL: $cfg stdout diverges"; diff "$WORK/n0.txt" "$WORK/n1.txt" | head; exit 1; }
      else
        echo "FAIL: $cfg $f diverges"; exit 1;
      fi
    fi
  done
  echo "   $cfg byte-identical (stdout/stderr/scalars)"
done

echo "-- repeated-invocation run (same config twice against the port)"
run_case caseA port runA-port2
for f in stdout.txt stderr.txt scalars.dat; do
  if ! cmp -s "$WORK/run-caseA-port/$f" "$WORK/runA-port2/$f"; then
    if [ "$f" = "stdout.txt" ]; then
      norm "$WORK/run-caseA-port/$f" > "$WORK/r0.txt"
      norm "$WORK/runA-port2/$f" > "$WORK/r1.txt"
      diff "$WORK/r0.txt" "$WORK/r1.txt" > /dev/null 2>&1 \
        || { echo "FAIL: rerun stdout diverges"; diff "$WORK/r0.txt" "$WORK/r1.txt" | head; exit 1; }
    else
      echo "FAIL: rerun $f diverges"; exit 1;
    fi
  fi
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- branch-liveness control (configs must mutually differ; arms agree)"
for cfg in caseA caseB caseC caseD; do
  norm "$WORK/run-${cfg}-port/stdout.txt" > "$WORK/live-${cfg}.out"
  cp "$WORK/run-${cfg}-port/scalars.dat" "$WORK/live-${cfg}.dat"
done
if cmp -s "$WORK/live-caseA.dat" "$WORK/live-caseB.dat" \
  && cmp -s "$WORK/live-caseA.out" "$WORK/live-caseB.out"; then
  echo "FAIL: A vs B outputs identical (no branch liveness)"; exit 1;
fi
if cmp -s "$WORK/live-caseC.dat" "$WORK/live-caseD.dat" \
  && cmp -s "$WORK/live-caseC.out" "$WORK/live-caseD.out"; then
  echo "FAIL: C vs D outputs identical (no branch liveness)"; exit 1;
fi
echo "   liveness confirmed (A/B and C/D differ; arms agreed above)"

echo "RIGIDMESHMAPPER PROOF GREEN"
echo "wall seconds: $((SECONDS - T0))"
