#!/bin/bash
# T022-A CraigBamptonSolver solver-entry proof: the Rust cdylib replaces
# CraigBamptonSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `craigbamptonsolver_` + `craigbamptonsolver_init_`);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines — config A (Craig-Bampton eigenmode
#     reduction staging: stress eigensolve + constraint modes, then the
#     CB solver with `Calculate Matrix Norm`, `Save Reduction Basis` and
#     `Save Reduction Matrix` all true) and config B (same staging with
#     the norm and basis branches off: no `Calculate Matrix Norm`, no
#     `Save Reduction Basis`); both carry a SaveScalars solver so the
#     eigensolve lands in scalars.dat;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry; the port holds no cross-call state — the
#     reduction work arrays are allocated and dropped every call).
# Branch-liveness control: config A writes ReductionModes.dat and the
# `Matrix norm ...` reports while config B writes neither, in both arms
# (the keyword branches change the output identically in both arms), so
# the differential is not vacuous. The nondeterministic wall-clock lines
# (`Elapsed REAL time`, `SOLVER TOTAL TIME`, `Solver time (CPU,REAL)`,
# `STARTED/FINISHED AT`) are structurally filtered in both arms.
#
# Needs a configured+built serial Elmer tree (ELMER_BUILD — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout only for the cube.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-craigbamptonsolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/cube*
echo "== T022-A CraigBamptonSolver proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libcraigbamptonsolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/CraigBamptonSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the CraigBamptonSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T craigbamptonsolver_" "$WORK/port.syms" \
  || { echo "FAIL: craigbamptonsolver_ not the exported entry"; exit 1; }
grep -qx "T craigbamptonsolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: craigbamptonsolver_init_ not exported"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's cube.grd)"
cp "$ELMER_FORK/fem/tests/StressConstraintModes3/cube.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 cube > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  CHECK KEYWORDS Warn
  Mesh DB "." "cube"
End
Simulation
  Max Output Level = 5
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1000.0
  Youngs Modulus = 1.0e9
  Poisson Ratio = 0.3
End
Solver 1
  Equation = "Stress Analysis"
  Variable = "Displacement"
  Variable Dofs = 3
  Constraint Modes Analysis = Logical True
  Constraint Modes Mass Lumping = Logical True
  Eigen Analysis = True
  Eigen System Values = 12
  Eigen System Convergence Tolerance = Real 1.0e-6
  Eigen System Normalize To Unity = Logical False
  Eigen System Max Iterations = 50
  Linear System Solver = string "Iterative"
  Linear System Convergence Tolerance = real 1e-8
  Linear System Iterative Method = BiCGStabl
  BiCGStabl Polynomial Degree = 4
  Linear System Max Iterations = 2000
  Linear System Preconditioning = BILU1
  Linear System Residual Output = 0
  Solver Timing = Logical False
  Displace Mesh = Logical False
End
Solver 2
  Equation = "CB"
  Procedure = "CraigBamptonSolver" "CraigBamptonSolver"
  Elasticity Solver Name = String "Stress Analysis"
  Calculate Matrix Norm = Logical True
  Save Reduction Basis = Logical True
  Save Reduction Matrix = Logical True
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Save Eigenvalues = Logical True
End
Boundary Condition 1
  Name = "Bottom"
  Target Boundaries(1) = 1
  Constraint Modes Displacement = Logical True
End
Boundary Condition 2
  Name = "Top"
  Target Boundaries = 3
  Constraint Modes Displacement = Logical True
End
Boundary Condition 3
  Name = "Else"
  Target Boundaries(4) = 2 4 5 6
End
EOF
# Config B: same eigenmode-reduction staging, but the norm and basis
# branches off (exercises the other side of both ValueList branches).
sed -e 's/  Calculate Matrix Norm = Logical True//' \
    -e 's/  Save Reduction Basis = Logical True//' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -c "CraigBamptonSolver" "$WORK/caseB.sif"
grep -c "Save Reduction Matrix" "$WORK/caseB.sif"
[ "$(grep -c "Calculate Matrix Norm" "$WORK/caseB.sif")" -eq 0 ] \
  || { echo "FAIL: case B still sets Calculate Matrix Norm"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/CraigBamptonSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/CraigBamptonSolver.so"
nm -D --defined-only "$WORK/shadow-port/CraigBamptonSolver.so" | grep -q " T craigbamptonsolver_" \
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
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|Solver time \(CPU,REAL\)" "$1" || true
}

echo "-- differentials A/B + rerun of A"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat ReductionHeader.dat ReductionIndeces.dat ReducedStiff.dat ReducedMass.dat ReductionModes.dat; do
  # Config B never writes ReductionModes.dat (basis branch off); compare
  # only where the reference arm produced the file, and require the port
  # arm to agree on its presence.
  if [ "$f" = "ReductionModes.dat" ]; then
    test -f "$WORK/runA-orig/$f" || { echo "FAIL: config A orig missing $f"; exit 1; }
    test '!' -f "$WORK/runB-orig/$f" || { echo "FAIL: config B orig unexpectedly wrote $f"; exit 1; }
    test -f "$WORK/runA-port/$f" || { echo "FAIL: config A port missing $f"; exit 1; }
    test '!' -f "$WORK/runB-port/$f" || { echo "FAIL: config B port unexpectedly wrote $f"; exit 1; }
  fi
  if [ -f "$WORK/runA-orig/$f" ]; then
    norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
    norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
    diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
      { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
    norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
    diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
      { echo "FAIL: rerun $f diverges"; exit 1; }
  fi
  if [ -f "$WORK/runB-orig/$f" ]; then
    norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
    norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
    diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
      { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
  fi
done
echo "   configs A+B byte-identical (stdout/stderr/scalars.dat/reduction files) + rerun byte-identical"

echo "-- branch-liveness control (A writes basis+norms, B writes neither; both arms agree)"
for arm in orig port; do
  grep -q "Matrix norm for reduced stiffness matrix" "$WORK/runA-$arm/stdout.txt" \
    || { echo "FAIL: no stiffness norm report in A($arm)"; exit 1; }
  grep -q "Matrix norm for reduced stiffness matrix" "$WORK/runB-$arm/stdout.txt" \
    && { echo "FAIL: stiffness norm report leaks into B($arm)"; exit 1; }
  test -f "$WORK/runA-$arm/ReductionModes.dat" \
    || { echo "FAIL: no ReductionModes.dat in A($arm)"; exit 1; }
  test '!' -f "$WORK/runB-$arm/ReductionModes.dat" \
    || { echo "FAIL: ReductionModes.dat leaks into B($arm)"; exit 1; }
done
echo "   control holds: branches change the output in both arms"

echo "CRAIGBAMPTONSOLVER PROOF GREEN wall=$((SECONDS - T0))s"
