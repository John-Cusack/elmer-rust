#!/bin/bash
# T022-A CoordinateTransform solver-entry proof: the Rust cdylib replaces
# CoordinateTransform.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `rotmsolver_` + `rotmsolver_init_` +
#     `__vecutils_MOD_jac/norm2/normalized/transform2rank`);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines — config A (alpha/beta-gradient path with
#     polar decomposition, 1 steady iteration); config B (gamma path via
#     `Local Coordinate System Beta Reference and Gamma` + coilcurrent,
#     rotated references, `Use Polar Decomposition = False` with the two
#     Tolerance/Iteration Warn lines, CYCLE on the body without a Beta
#     Reference, 2 steady iterations exercising SAVE re-entry);
#     both carry a `SaveScalars` solver so the computed fields land in
#     scalars.dat;
#  3. a second full process run of config B is byte-identical
#     (idempotent re-entry through the SAVE state).
# Branch-liveness control: the saved field statistics and the Warn lines
# differ across configs A/B in both arms, so the differential is not
# vacuous. The nondeterministic `Level=5` `RealTime` reports are
# structurally suppressed in both arms via the standard
# `Max Output Level = 4` keyword (the deterministic reports still print
# and are compared).
#
# Needs a configured+built serial Elmer tree (ELMER_BUILD — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-coordinatetransform"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A CoordinateTransform proof: work=$WORK build=$ELMER_BUILD"
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
CTD="${CARGO_TARGET_DIR:-$LANE/scratch/cargo-target-coordinatetransform}"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" CARGO_TARGET_DIR="$CTD" cargo +nightly build --release -q)
PORT_SO="$CTD/release/libcoordinatetransform.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/CoordinateTransform.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the CoordinateTransform module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T rotmsolver_" "$WORK/port.syms" \
  || { echo "FAIL: rotmsolver_ not the exported entry"; exit 1; }
grep -qx "T rotmsolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: rotmsolver_init_ not exported"; exit 1; }
grep -qx "T __vecutils_MOD_jac" "$WORK/port.syms" \
  || { echo "FAIL: __vecutils_MOD_jac not exported"; exit 1; }
grep -qx "T __vecutils_MOD_norm2" "$WORK/port.syms" \
  || { echo "FAIL: __vecutils_MOD_norm2 not exported"; exit 1; }
grep -qx "T __vecutils_MOD_normalized" "$WORK/port.syms" \
  || { echo "FAIL: __vecutils_MOD_normalized not exported"; exit 1; }
grep -qx "T __vecutils_MOD_transform2rank" "$WORK/port.syms" \
  || { echo "FAIL: __vecutils_MOD_transform2rank not exported"; exit 1; }
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
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Max Output Level = Integer 4
End
Body 1
  Equation = 1
  Material = 1
  Initial Condition = 1
  Alpha Reference(3) = Real 1.0 0.0 0.0
  Beta Reference(3) = Real 0.0 1.0 0.0
End
Body 2
  Equation = 1
  Material = 1
  Initial Condition = 1
  Alpha Reference(3) = Real 1.0 0.0 0.0
  Beta Reference(3) = Real 0.0 1.0 0.0
End
Body 3
  Equation = 1
  Material = 1
  Initial Condition = 1
  Alpha Reference(3) = Real 1.0 0.0 0.0
  Beta Reference(3) = Real 0.0 1.0 0.0
End
Body 4
  Equation = 1
  Material = 1
  Initial Condition = 1
  Alpha Reference(3) = Real 1.0 0.0 0.0
  Beta Reference(3) = Real 0.0 1.0 0.0
End
Initial Condition 1
  alpha = Variable Coordinate 1
    Real MATC "tx"
  beta = Variable Coordinate 2
    Real MATC "tx"
  coilcurrent e 1 = Variable Coordinate 1
    Real MATC "tx"
  coilcurrent e 2 = Variable Coordinate 2
    Real MATC "tx"
  coilcurrent e 3 = Real 1.0
End
Equation 1
  Active Solvers(5) = 1 2 3 4 5
End
Material 1
  Density = 1.0
End
Solver 1
  Equation = AllocAlpha
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = alpha
  Variable DOFs = 1
End
Solver 2
  Equation = AllocBeta
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = beta
  Variable DOFs = 1
End
Solver 3
  Equation = AllocCoil
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = String "coilcurrent e"
  Variable DOFs = 3
End
Solver 4
  Equation = Transformation matrix
  Procedure = "CoordinateTransform" "RotMSolver"
  Variable = -nooutput Dummy
  Optimize Bandwidth = False
  Discontinuous Galerkin = True
  Polar Decomposition Determinant Tolerance = Real 1.0e-9
  Polar Decomposition Max Iterations = Integer 100
  Exported Variable 1 = RotM E[RotM E:9]
  Exported Variable 2 = Alpha Vector E[Alpha Vector E:3]
  Exported Variable 3 = Beta Vector E[Beta Vector E:3]
  Exported Variable 4 = Gamma Vector E[Gamma Vector E:3]
End
Solver 5
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = RotM E
  Operator 1 = min
  Variable 2 = RotM E
  Operator 2 = max
  Variable 3 = RotM E
  Operator 3 = sum
  Variable 4 = Alpha Vector E
  Operator 4 = sum
  Variable 5 = Gamma Vector E
  Operator 5 = sum
End
EOF
python3 - "$WORK/caseA.sif" "$WORK/caseB.sif" <<'PYEOF'
import sys
s = open(sys.argv[1]).read()
# Config B: gamma path, rotated body references, no polar decomposition (the
# two Warn defaults), 2 steady iterations (SAVE re-entry), Body 4 without a
# Beta Reference (CYCLE path).
s = s.replace("Steady State Max Iterations = 1",
              "Steady State Max Iterations = 2")
s = s.replace("  Alpha Reference(3) = Real 1.0 0.0 0.0",
              "  Alpha Reference(3) = Real 0.0 1.0 0.0")
s = s.replace("  Beta Reference(3) = Real 0.0 1.0 0.0",
              "  Beta Reference(3) = Real -1.0 0.0 0.0\n  Local Coordinate System Beta Reference and Gamma = Logical True")
s = s.replace("  Polar Decomposition Determinant Tolerance = Real 1.0e-9",
              "  Use Polar Decomposition = Logical False")
s = s.replace("  Polar Decomposition Max Iterations = Integer 100\n", "")
# Body 4 loses its Beta Reference so its elements take the CYCLE skip.
head, sep, tail = s.partition("Body 4\n")
assert sep
tail = tail.replace("  Beta Reference(3) = Real -1.0 0.0 0.0\n", "", 1)
assert "Beta Reference(3)" not in tail.split("End")[0]
s = head + sep + tail
open(sys.argv[2], "w").write(s)
PYEOF
grep -c "Local Coordinate System Beta Reference and Gamma" "$WORK/caseB.sif"
grep -c "Use Polar Decomposition" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/CoordinateTransform.so"
cp "$PORT_SO" "$WORK/shadow-port/CoordinateTransform.so"
nm -D --defined-only "$WORK/shadow-port/CoordinateTransform.so" | grep -q " T rotmsolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" || true
}

echo "-- differential A (alpha/beta path + polar decomposition)"
run_case caseA orig runA-orig
run_case caseA port runA-port
echo "-- differential B (gamma path, no polar decomposition, CYCLE, 2 iterations) + rerun"
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseB port runB-port2
for f in stdout.txt stderr.txt scalars.dat scalars.dat.names; do
  [ -f "$WORK/runA-orig/$f" ] || continue
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
  [ -f "$WORK/runB-orig/$f" ] || continue
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
  norm "$WORK/runB-port2/$f" > "$WORK/b3.txt"
  diff "$WORK/b2.txt" "$WORK/b3.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   configs A+B byte-identical (stdout/stderr/scalars) + rerun byte-identical"

echo "-- branch-liveness control (saved statistics differ across A/B, both arms)"
for arm in orig port; do
  grep -q "[1-9]" "$WORK/runA-$arm/scalars.dat" \
    || { echo "FAIL: A($arm) statistics all zero (solver wrote nothing?)"; exit 1; }
  cmp -s "$WORK/runA-$arm/scalars.dat" "$WORK/runB-$arm/scalars.dat" \
    && { echo "FAIL: A($arm) and B($arm) statistics identical"; exit 1; }
  grep -q "Polar Decomposition Determinant Tolerance not set" "$WORK/runB-$arm/stdout.txt" \
    || grep -q "Polar Decomposition Determinant Tolerance not set" "$WORK/runB-$arm/stderr.txt" \
    || { echo "FAIL: B($arm) missing tolerance Warn line"; exit 1; }
done
echo "   control holds: branches change the numerics in both arms"

echo "COORDINATETRANSFORM PROOF GREEN wall=$((SECONDS - T0))s"
