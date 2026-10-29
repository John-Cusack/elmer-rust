#!/bin/bash
# T021-B StructuredFlowLine port proof: the Rust cdylib replaces
# StructuredFlowLine.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `structuredflowline_`, `structuredflowline_init_`);
#  2. two distinct clean-room .sif differentials (masked + TrueLines
#     config A; unmasked + displacement-mode + AveMode 1 config B with 2
#     steady iterations, i.e. the repeated-invocation run) are
#     byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. the branch-liveness control (A vs B) produces different output on
#     both arms, and both arms agree on each config.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial configure —
# both arms share it, so the differential is self-consistent).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-structuredflowline"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/slab*
echo "== T021-B StructuredFlowLine proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libstructuredflowline.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/StructuredFlowLine.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the StructuredFlowLine module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T structuredflowline_" "$WORK/port.syms" \
  || { echo "FAIL: structuredflowline_ not an exported entry"; exit 1; }
grep -qx "T structuredflowline_init_" "$WORK/port.syms" \
  || { echo "FAIL: structuredflowline_init_ not an exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- clippy (no warnings)"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release -q -- -D warnings) \
  || { echo "FAIL: clippy warnings"; exit 1; }
echo "   clippy clean"

echo "-- mesh (ElmerGrid from the upstream ExtrusionStructured2D slab.grd)"
cp "$ELMER_FORK/fem/tests/ExtrusionStructured2D/slab.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 slab > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "slab"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
End
Body 1
  Name = "Fluid"
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Initial Condition 1
  Velocity 1 = 0
  Velocity 2 = -0.5
End
Equation 1
  Name = "Fluid"
  Active Solvers(3) = 1 2 3
End
Material 1
  Name = "Liquid"
  Density = 10.0
  Viscosity = 1.0
End
Solver 1
  Equation = "Navier-Stokes"
  Procedure = "FlowSolve" "FlowSolver"
  Variable = Flow Solution[Velocity:2 Pressure:1]
  Stabilization Method = String Stabilized
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
  Linear System Residual Output = 10
  Linear System Abort Not Converged = False
  Nonlinear System Convergence Tolerance = 1.0e-4
  Nonlinear System Max Iterations = 1
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Newton After Tolerance = 1.0e-2
  Nonlinear System Relaxation Factor = 1.0
  Steady State Convergence Tolerance = 1.0e-5
End
Solver 2
  Equation = "FlowLine"
  Procedure = "StructuredFlowLine" "StructuredFlowLine"
  Nonlinear System Relaxation Factor = 0.7
  Active Coordinate = Integer 2
  Velocity Variable Name = String "Flow Solution"
  Displacement Mode = Logical False
  True Flow Line Iterations = Integer 2
  Averaging Order = Integer 1
  Averaging Method = Integer 0
  Dot Product Tolerance = Real 0.001
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Flow Solution
  Operator 1 = norm
  Variable 2 = FlowLineDisp
  Operator 2 = norm
End
Boundary Condition 1
  Name = "Free"
  Target Boundaries = 1
  Flow Line = Logical True
End
Boundary Condition 2
  Name = "Walls"
  Target Boundaries = 2
  Velocity 1 = 0
  Velocity 2 = 0
End
Boundary Condition 3
  Name = "Outlet"
  Target Boundaries = 3
  Velocity 1 = Real 0.0
End
Boundary Condition 4
  Name = "Inlet"
  Target Boundaries = 4
  Velocity 1 = 0
  Velocity 2 = -0.5
End
EOF
# Config B: no Flow Line mask (MaskExist=.FALSE. path), displacement mode
# on (defined: the mask guard is skipped), TrueLines only on the first of
# 2 steady iterations (VisitedTimes re-entry), AveMode 1 retry path.
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    -e 's/Displacement Mode = Logical False/Displacement Mode = Logical True/' \
    -e 's/True Flow Line Iterations = Integer 2/True Flow Line Iterations = Integer 1/' \
    -e 's/Averaging Method = Integer 0/Averaging Method = Integer 1/' \
    -e '/^  Flow Line = Logical True$/d' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -q "Displacement Mode = Logical True" "$WORK/caseB.sif"
grep -q "Steady State Max Iterations = 2" "$WORK/caseB.sif"
! grep -q "Flow Line = Logical True" "$WORK/caseB.sif" \
  || { echo "FAIL: mask line not removed from case B"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/StructuredFlowLine.so"
cp "$PORT_SO" "$WORK/shadow-port/StructuredFlowLine.so"
nm -D --defined-only "$WORK/shadow-port/StructuredFlowLine.so" | grep -q "T structuredflowline_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/StructuredFlowLine.so" | grep -q "T structuredflowline_init_" \
  || { echo "FAIL: staged port .so lost its init entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/slab" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER SOLVER.*STARTED AT|ELMER.*FINISHED AT|STARTED AT|pointer init time|Streamline mapping time|Solve:.*\(s\)" "$1" || true
}

echo "-- differential A (masked + TrueLines) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 \
    || { echo "FAIL: config A $f diverges"; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  if [ "$f" = stdout.txt ]; then
    norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
    norm "$WORK/runA-port2/stdout.txt" > "$WORK/a2.txt"
    diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 \
      || { echo "FAIL: rerun stdout diverges"; exit 1; }
  else
    diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" > /dev/null 2>&1 \
      || { echo "FAIL: rerun $f diverges"; exit 1; }
  fi
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (unmasked + displacement + 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stderr.txt scalars.dat; do
  diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" > /dev/null 2>&1 \
    || { echo "FAIL: config B $f diverges"; exit 1; }
done
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b0.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b1.txt"
diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config B stdout diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A vs B differ on both arms)"
diff "$WORK/a0.txt" "$WORK/b0.txt" > /dev/null 2>&1 \
  && { echo "FAIL: configs A and B agree on orig (control has no branch delta)"; exit 1; }
diff "$WORK/a1.txt" "$WORK/b1.txt" > /dev/null 2>&1 \
  && { echo "FAIL: configs A and B agree on port (control has no branch delta)"; exit 1; }
echo "   A-vs-B outputs differ on both arms, arms agree per config"

echo "ELMER STRUCTUREDFLOWLINE PROOF GREEN (wall $((SECONDS - T0))s)"
