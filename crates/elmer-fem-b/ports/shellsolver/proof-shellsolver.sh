#!/bin/bash
# T021-B ShellSolver port proof: the Rust cdylib replaces ShellSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `shellsolver_` + `shellsolver_init0_`);
#  2. three distinct .sif differentials (A: blending+nonlinear;
#     B: cartesian+nonlinear; C: blending+linear) are byte-identical to
#     the reference except wall-clock lines (stdout normed, stderr and
#     scalars.dat raw);
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch-liveness: configs A/B/C produce pairwise-different
#     solutions, and both arms agree on each.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config —
# both arms share it, so the differential is self-consistent).
# The strip mesh (with mesh.director) is the upstream Shell_Cantilever
# mesh: a shell solver needs a 3-D shell mesh, so the Poisson square.grd
# note does not apply (per the lane brief's mesh note).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-shellsolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/strip
echo "== T021-B ShellSolver proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran 13.3.0 required"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: nightly 1.100.0 required"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libshellsolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ShellSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ShellSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T shellsolver_" "$WORK/port.syms" \
  || { echo "FAIL: shellsolver_ not the exported entry"; exit 1; }
grep -qx "T shellsolver_init0_" "$WORK/port.syms" \
  || { echo "FAIL: shellsolver_init0_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (upstream Shell_Cantilever strip mesh with mesh.director)"
cp -r "$ELMER_FORK/fem/tests/Shell_Cantilever/strip" "$WORK/"

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "strip"
End

Simulation
  Max Output Level = 5
  Coordinate System = Cartesian 3D
  Simulation Type = Steady
  Output Intervals = 1
  Steady State Max Iterations = 1
End

Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End

Material 1
  Poisson Ratio = Real 0.0
  Youngs Modulus = Real 1.2e+6
  Shell Thickness = Real 1.0e-1
End

Body Force 1
  Normal Pressure = Real 0.0
End

Equation 1 :: Active Solvers(1) = 1

Solver 1
  Equation = "Shell equations"
  Procedure = "ShellSolver" "ShellSolver"

  Large Deflection = True
  Displace Mesh = True
  Nonlinear System Max Iterations = 40
  Nonlinear System Convergence Tolerance = 1e-7

  Linear System Solver = "Iterative"
  Linear System Preconditioning = ILU2
  Linear System Max Iterations = 1000
  Linear System Normwise Backward Error = Logical True
  Linear System Convergence Tolerance = 1e-9
  Linear System Iterative Method = GCR
  Linear System GCR Restart = 100
  Linear System Abort Not Converged = False
  Steady State Convergence Tolerance = 1e-09
End

Solver 2
  Equation = "Save Scalars"
  Exec Solver = After Timestep
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = U
  Save Points(2) = 17 34
End

Boundary Condition 1
  Target Boundaries(1) = 4
  U 1 = Real 0
  U 2 = Real 0
  U 3 = Real 0
  DNU 1 = Real 0
  DNU 2 = Real 0
  DNU 3 = Real 0
End

Boundary Condition 2
  Target Boundaries(1) = 2
  Resultant Force 1 = Real 0.0
  Resultant Force 2 = Real 0.0
  Resultant Force 3 = Variable Time
    Real MATC "0.2*tx"
End

Solver 1 :: Reference Norm = Real 1.39293663E-01
Solver 1 :: Reference Norm Tolerance = Real 1.0e-5
EOF
sed -e 's/Large Deflection = True/Large Deflection = True\n  Cartesian Formulation = True/' \
    -e '/Reference Norm/d' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
sed -e 's/Large Deflection = True/Large Deflection = False/' \
    -e '/Nonlinear System Max Iterations/d; /Nonlinear System Convergence Tolerance/d' \
    -e '/Reference Norm/d' \
    "$WORK/caseA.sif" > "$WORK/caseC.sif"
grep -c "Cartesian Formulation" "$WORK/caseB.sif"
grep -c "Large Deflection = False" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ShellSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/ShellSolver.so"
nm -D --defined-only "$WORK/shadow-port/ShellSolver.so" | grep -q "T shellsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/strip" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Assembly \(s\)|Solve \(s\)|pointer init time|Time spent for restart|Factorization time" "$1" || true; }

check_equal() { # <case> <what>
  local case=$1
  for f in stderr.txt scalars.dat; do
    diff "$WORK/run${case}-orig/$f" "$WORK/run${case}-port/$f" > /dev/null 2>&1 || \
      { echo "FAIL: config $case $f diverges"; diff "$WORK/run${case}-orig/$f" "$WORK/run${case}-port/$f" | head; exit 1; }
  done
  norm "$WORK/run${case}-orig/stdout.txt" > "$WORK/${case}0.txt"
  norm "$WORK/run${case}-port/stdout.txt" > "$WORK/${case}1.txt"
  diff "$WORK/${case}0.txt" "$WORK/${case}1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config $case stdout diverges"; diff "$WORK/${case}0.txt" "$WORK/${case}1.txt" | head; exit 1; }
  echo "   config $case byte-identical (stdout/stderr/scalars)"
}

echo "-- differential A (blending + nonlinear) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
check_equal A
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (cartesian + nonlinear)"
run_case caseB orig runB-orig
run_case caseB port runB-port
check_equal B

echo "-- differential C (blending + linear)"
run_case caseC orig runC-orig
run_case caseC port runC-port
check_equal C

echo "-- branch-liveness control (configs must differ; arms agree)"
norm "$WORK/runA-orig/stdout.txt" > "$WORK/la.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/lb.txt"
norm "$WORK/runC-orig/stdout.txt" > "$WORK/lc.txt"
if cmp -s "$WORK/la.txt" "$WORK/lb.txt"; then echo "FAIL: A/B liveness outputs identical"; exit 1; fi
if cmp -s "$WORK/la.txt" "$WORK/lc.txt"; then echo "FAIL: A/C liveness outputs identical"; exit 1; fi
echo "   A/B/C solutions pairwise differ and both arms agree on each"

echo "-- clippy (no warnings)"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release 2>&1 | tee "$WORK/clippy.txt" | grep -E "^(warning|error)" && { echo "FAIL: clippy warnings"; exit 1; } || true)
echo "   clippy clean"

T1=$SECONDS
echo "SHELLSOLVER PROOF GREEN (wall $((T1-T0))s)"
