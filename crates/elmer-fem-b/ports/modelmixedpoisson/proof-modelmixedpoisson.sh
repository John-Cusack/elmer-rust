#!/bin/bash
# T021-B ModelMixedPoisson solver-entry proof: the Rust cdylib replaces
# ModelMixedPoisson.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `mixedpoisson_`, `mixedpoisson_init0_`,
#     `mixedpoisson_post_`);
#  2. two distinct clean-room .sif differentials (config A: first-kind
#     basis + Material Parameter + Source Field, no BCs, 1 steady
#     iteration; config B: second-kind basis + Material Tensor +
#     Scalar Field boundary condition, 2 steady iterations) are
#     byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch liveness: A and B produce different outputs, and both arms
#     agree on each config (B also re-enters SAVE state in-process via
#     its 2nd steady iteration).
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-modelmixedpoisson"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B ModelMixedPoisson proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libmodelmixedpoisson.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- clippy (no warnings)"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release 2>&1 | tee "$WORK/clippy.log" | grep -E "^(warning|error)" && { echo "FAIL: clippy warnings/errors"; exit 1; } || true)
echo "   clippy clean"

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ModelMixedPoisson.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ModelMixedPoisson module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in mixedpoisson_ mixedpoisson_init0_ mixedpoisson_post_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
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
  Coordinate System = "Cartesian 2D"
  Simulation Type = Steady State
  Steady State Max Iterations = 1
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
End
Body 3
  Equation = 1
  Material = 1
End
Body 4
  Equation = 1
  Material = 1
  Body Force = 2
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Material Parameter = Real 1.0
End
Body Force 1
  Source Field = Real 1.0
End
Body Force 2
  Source Field = Real -1.0
End
Solver 1
  Equation = "Mixed Poisson"
  Variable = "MixedSol"
  Variable DOFs = 1
  Procedure = "ModelMixedPoisson" "MixedPoisson"
  Second Kind Basis = Logical False
  Optimize Bandwidth = False
  Bubbles in Global System = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = GCR
  Linear System Preconditioning = ILU0
  Linear System Row Equilibration = Logical True
  Linear System Max Iterations = 100
  Linear System Convergence Tolerance = 1.0e-6
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = mixedpot
  Operator 1 = min
  Variable 2 = mixedpot
  Operator 2 = max
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = "Cartesian 2D"
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
End
Body 3
  Equation = 1
  Material = 1
End
Body 4
  Equation = 1
  Material = 1
  Body Force = 2
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Material Tensor = Real 0.5
End
Body Force 1
  Source Field = Real 2.0
End
Body Force 2
  Source Field = Real -2.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Scalar Field = Real 2.0
  MixedSol Condition = Real 1.0
End
Solver 1
  Equation = "Mixed Poisson"
  Variable = "MixedSol"
  Variable DOFs = 1
  Procedure = "ModelMixedPoisson" "MixedPoisson"
  Second Kind Basis = Logical True
  Optimize Bandwidth = False
  Bubbles in Global System = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = GCR
  Linear System Preconditioning = ILU0
  Linear System Row Equilibration = Logical True
  Linear System Max Iterations = 100
  Linear System Convergence Tolerance = 1.0e-6
 End
 Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = mixedpot
  Operator 1 = min
  Variable 2 = mixedpot
  Operator 2 = max
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
rm -f "$WORK/shadow-port/ModelMixedPoisson.so"
cp "$PORT_SO" "$WORK/shadow-port/ModelMixedPoisson.so"
nm -D --defined-only "$WORK/shadow-port/ModelMixedPoisson.so" | grep -q "T mixedpoisson_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (first-kind + parameter + source) + rerun (re-entry)"
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
   norm "$WORK/runA-port/$f" > "$WORK/r0.txt"
   norm "$WORK/runA-port2/$f" > "$WORK/r1.txt"
   diff "$WORK/r0.txt" "$WORK/r1.txt" > /dev/null 2>&1 || \
     { echo "FAIL: rerun $f diverges (SAVE re-entry)"; diff "$WORK/r0.txt" "$WORK/r1.txt" | head; exit 1; }
 done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (second-kind + tensor + Scalar Field BC, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
  diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (A and B must differ, on both arms)"
diff "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: configs A and B agree on port (no branch divergence)"; exit 1; }
diff "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: configs A and B agree on orig (no branch divergence)"; exit 1; }
echo "   A/B diverge on both arms (liveness control holds)"
 echo "MODELMIXEDPOISSON PROOF GREEN in $((SECONDS - T0))s"
