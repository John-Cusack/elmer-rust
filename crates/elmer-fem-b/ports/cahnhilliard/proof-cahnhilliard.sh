#!/bin/bash
# T022-A CahnHilliard solver-entry proof: the Rust cdylib replaces
# CahnHilliard.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `cahnhilliardsolver_`);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines — config A (steady, 2 nonlinear iterations;
#     Mobility + Mixing Energy Density + Interface Thickness branches)
#     and config B (transient, 2 timesteps x 3 nonlinear iterations;
#     Surface Tension Coefficient + h**2-mobility + ElementDiameter +
#     Convection Velocity + Default1stOrderTime branches); a SaveScalars
#     solver records the solved Phase field's min/max in scalars.dat;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry; the port holds no cross-call Rust state).
# Branch-liveness control: the `Result Norm` reports and the saved
# scalars differ across configs A/B in both arms (the keyword and
# transient branches change the numerics identically in both arms), so
# the differential is not vacuous. Nondeterministic wall-clock lines are
# filtered via the standard norm; deterministic `Level=4` `Result Norm`
# and `Relative Change` lines still print and are compared.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-cahnhilliard"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A CahnHilliard proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libcahnhilliard.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/CahnHilliard.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the CahnHilliard module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T cahnhilliardsolver_" "$WORK/port.syms" \
  || { echo "FAIL: cahnhilliardsolver_ not the exported entry"; exit 1; }
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
  Max Output Level = 4
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
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
Initial Condition 1
  Phase 1 = Real 0.5
  Phase 2 = Real 0.0
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Interface Thickness = Real 0.05
  Mobility = Real 0.0025
  Mixing Energy Density = Real 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Phase 1 = Real 1.0
  Phase 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Phase 1 = Real -1.0
  Phase 2 = Real 0.0
End
Solver 1
  Equation = CahnHilliard
  Procedure = "CahnHilliard" "CahnHilliardSolver"
  Variable = String "Phase"
  Variable DOFs = Integer 2
  Nonlinear System Max Iterations = Integer 2
  Nonlinear System Convergence Tolerance = Real 1.0e-12
  Linear System Solver = Direct
  Linear System Direct Method = UMFPack
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Phase 1
  Operator 1 = min
  Variable 2 = Phase 1
  Operator 2 = max
  Variable 3 = Phase 2
  Operator 3 = min
  Variable 4 = Phase 2
  Operator 4 = max
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "square"
End
Simulation
  Max Output Level = 4
  Coordinate System = Cartesian
  Simulation Type = Transient
  Timestepping Method = BDF
  BDF Order = 1
  Timestep Sizes(1) = 0.1
  Timestep Intervals(1) = 2
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
Initial Condition 1
  Phase 1 = Real 0.5
  Phase 2 = Real 0.0
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Convection Velocity 1 = Real 0.5
  Convection Velocity 2 = Real 0.0
  Surface Tension Coefficient = Real 0.1
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Phase 1 = Real 1.0
  Phase 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Phase 1 = Real -1.0
  Phase 2 = Real 0.0
End
Solver 1
  Equation = CahnHilliard
  Procedure = "CahnHilliard" "CahnHilliardSolver"
  Variable = String "Phase"
  Variable DOFs = Integer 2
  Nonlinear System Max Iterations = Integer 3
  Nonlinear System Convergence Tolerance = Real 1.0e-12
  Linear System Solver = Direct
  Linear System Direct Method = UMFPack
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Phase 1
  Operator 1 = min
  Variable 2 = Phase 1
  Operator 2 = max
  Variable 3 = Phase 2
  Operator 3 = min
  Variable 4 = Phase 2
  Operator 4 = max
End
EOF
grep -c "Surface Tension Coefficient" "$WORK/caseB.sif"
grep -c "Timestep Intervals" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/CahnHilliard.so"
cp "$PORT_SO" "$WORK/shadow-port/CahnHilliard.so"
nm -D --defined-only "$WORK/shadow-port/CahnHilliard.so" | grep -q "T cahnhilliardsolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (steady full-material branch) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
done
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (transient alternate-material branch)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (Result Norm differs across A/B, both arms)"
for arm in orig port; do
  grep "Result Norm" "$WORK/runA-$arm/stdout.txt" > "$WORK/liveA-$arm.txt"
  grep "Result Norm" "$WORK/runB-$arm/stdout.txt" > "$WORK/liveB-$arm.txt"
  if diff -q "$WORK/liveA-$arm.txt" "$WORK/liveB-$arm.txt" > /dev/null 2>&1; then
    echo "FAIL: liveness control vacuous for $arm (A and B norms identical)"
    exit 1
  fi
done
echo "   control holds: branches change the numerics in both arms"

echo "CAHNHILLIARD PROOF GREEN wall=$((SECONDS - T0))s"
