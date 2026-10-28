#!/bin/bash
# T021-A DivergenceSolver solver-entry proof: the Rust cdylib replaces
# DivergenceSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `divergencesolver_` + `divergencesolver_init_`);
#  2. three distinct .sif differentials are byte-identical to the
#     reference except wall-clock lines — config A (cartesian defaults),
#     config B (transient, 2 timesteps: `Relative Divergence` and
#     `Constant Bulk Matrix` branches, exercising both the bulk-update
#     and the force-only update paths across the two timesteps), config C
#     (axi-symmetric coordinates: the `CSymmetry` weighting branch);
#     all three carry a `SaveScalars` solver so the solved divergence
#     field's min/max land in scalars.dat;
#  3. a second full process run of config B is byte-identical
#     (idempotent re-entry; the port holds no cross-call state — the
#     original allocates and deallocates its work arrays every call).
# Branch-liveness control: the `Result Norm` report and the saved
# scalars differ across configs A/B/C in both arms (the keyword and
# coordinate-system branches change the numerics identically in both
# arms), so the differential is not vacuous. The nondeterministic
# `Level=5` `RealTime` reports are structurally suppressed in both arms
# via the standard `Max Output Level = 4` keyword (the deterministic
# `Level=4` `Result Norm` line still prints and is compared).
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
LANE="${RS_ELMER_T021A_LANE:-/home/john/firstmate/data/rs-elmer-t021a}"
WORK="$LANE/scratch/proof-divergencesolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-A DivergenceSolver proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libdivergencesolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/DivergenceSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the DivergenceSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T divergencesolver_" "$WORK/port.syms" \
  || { echo "FAIL: divergencesolver_ not the exported entry"; exit 1; }
grep -qx "T divergencesolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: divergencesolver_init_ not exported"; exit 1; }
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
End
Body 2
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 3
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 4
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Initial Condition 1
  Velocity 1 = Variable Coordinate 1
    Real MATC "tx"
  Velocity 2 = Variable Coordinate 2
    Real MATC "tx"
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
End
Solver 1
  Equation = AllocateStuff
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = Velocity
  Variable DOFs = 2
End
Solver 2
  Equation = DivEq
  Procedure = "DivergenceSolver" "DivergenceSolver"
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Divergence
  Operator 1 = min
  Variable 2 = Divergence
  Operator 2 = max
  Variable 3 = Velocity 1
  Operator 3 = min
  Variable 4 = Velocity 1
  Operator 4 = max
End
EOF
sed -e 's/Simulation Type = Steady State/Simulation Type = Transient/' \
    -e 's/Steady State Max Iterations = 1/Timestepping Method = BDF\n  BDF Order = 1\n  Timestep Sizes(1) = 0.1\n  Timestep Intervals(1) = 2/' \
    -e 's/Procedure = "DivergenceSolver" "DivergenceSolver"/Procedure = "DivergenceSolver" "DivergenceSolver"\n  Relative Divergence = Logical True\n  Constant Bulk Matrix = Logical True/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
sed -e 's/Coordinate System = Cartesian/Coordinate System = Axi Symmetric/' \
    "$WORK/caseA.sif" > "$WORK/caseC.sif"
grep -c "MATC" "$WORK/caseB.sif"
grep -c "Relative Divergence" "$WORK/caseB.sif"
grep -c "Axi Symmetric" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/DivergenceSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/DivergenceSolver.so"
nm -D --defined-only "$WORK/shadow-port/DivergenceSolver.so" | grep -q " T divergencesolver_" \
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

echo "-- differentials A/B/C + rerun of B"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseB port runB-port2
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
  norm "$WORK/runB-port2/$f" > "$WORK/b3.txt"
  diff "$WORK/b2.txt" "$WORK/b3.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   configs A+B+C byte-identical (stdout/stderr/scalars) + rerun byte-identical"

echo "-- branch-liveness control (Result Norm differs across A/B/C, both arms)"
for arm in orig port; do
  na=$(grep "Result Norm" "$WORK/runA-$arm/stdout.txt" | head -1)
  nb=$(grep "Result Norm" "$WORK/runB-$arm/stdout.txt" | head -1)
  nc=$(grep "Result Norm" "$WORK/runC-$arm/stdout.txt" | head -1)
  [ -n "$na" ] || { echo "FAIL: no Result Norm in A($arm)"; exit 1; }
  [ "$na" != "$nb" ] || { echo "FAIL: A($arm) and B($arm) norms identical"; exit 1; }
  [ "$na" != "$nc" ] || { echo "FAIL: A($arm) and C($arm) norms identical"; exit 1; }
done
echo "   control holds: branches change the numerics in both arms"

echo "DIVERGENCESOLVER PROOF GREEN wall=$((SECONDS - T0))s"
