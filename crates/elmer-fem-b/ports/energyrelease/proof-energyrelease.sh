#!/bin/bash
# T021-A EnergyRelease solver-entry proof: the Rust cdylib replaces
# EnergyRelease.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `releaseratesolver_`);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines — config A (linear MATC source fields) and
#     config B (transient, 2 timesteps: quadratic MATC source fields,
#     exercising the per-timestep repeated-invocation path through the
#     SAVE work arrays);
#  3. a second full process run of config B is byte-identical
#     (idempotent re-entry; the SAVE work arrays persist across calls
#     with identical lifetimes in both arms).
# Branch-liveness control: the `2G =` fracture report differs between
# configs A and B in both arms (the source fields change the numerics
# identically in both arms), so the differential is not vacuous. The
# three source fields are staged through owning `AllocateSolver` solvers
# plus `MATC` initial conditions — the same mechanism in both arms.
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
WORK="$LANE/scratch/proof-energyrelease"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-A EnergyRelease proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libenergyrelease.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/EnergyRelease.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the EnergyRelease module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T releaseratesolver_" "$WORK/port.syms" \
  || { echo "FAIL: releaseratesolver_ not the exported entry"; exit 1; }
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
  True Displacement 1 = Variable Coordinate 1
    Real MATC "tx"
  True Displacement 2 = Variable Coordinate 2
    Real MATC "tx"
  True Displacement 3 = Real 0.0
  Stress 1 = Real 1.0
  Stress 2 = Real 2.0
  Stress 3 = Real 3.0
  Stress 4 = Real 0.5
  Stress 5 = Real 0.5
  Stress 6 = Real 0.5
  Shape Displacement 1 = Variable Coordinate 2
    Real MATC "tx"
  Shape Displacement 2 = Variable Coordinate 1
    Real MATC "tx"
  Shape Displacement 3 = Real 0.0
End
Equation 1
  Active Solvers(4) = 1 2 3 4
End
Material 1
  Density = 1.0
End
Solver 1
  Equation = AllocateDisp
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = True Displacement
  Variable DOFs = 3
End
Solver 2
  Equation = AllocateStress
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = Stress
  Variable DOFs = 6
End
Solver 3
  Equation = AllocateShape
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = Shape Displacement
  Variable DOFs = 3
End
Solver 4
  Equation = Fracture
  Procedure = "EnergyRelease" "ReleaseRateSolver"
End
EOF
sed -e 's/Simulation Type = Steady State/Simulation Type = Transient/' \
    -e 's/Steady State Max Iterations = 1/Timestepping Method = BDF\n  BDF Order = 1\n  Timestep Sizes(1) = 0.1\n  Timestep Intervals(1) = 2/' \
    -e 's/Real MATC "tx"/Real MATC "tx*tx"/g' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -c 'tx\*tx' "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/EnergyRelease.so"
cp "$PORT_SO" "$WORK/shadow-port/EnergyRelease.so"
nm -D --defined-only "$WORK/shadow-port/EnergyRelease.so" | grep -q " T releaseratesolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  local rc=0
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt) || rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (linear fields) + differential B (quadratic, transient x2) + rerun"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseB port runB-port2
for f in stdout.txt stderr.txt; do
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
done
echo "   configs A+B byte-identical (stdout/stderr) + rerun byte-identical"

echo "-- branch-liveness control (2G differs across A/B, both arms)"
for arm in orig port; do
  ga=$(grep "2G =" "$WORK/runA-$arm/stdout.txt" | head -1)
  gb=$(grep "2G =" "$WORK/runB-$arm/stdout.txt" | head -1)
  [ -n "$ga" ] || { echo "FAIL: no 2G report in A($arm)"; exit 1; }
  [ "$ga" != "$gb" ] || { echo "FAIL: A($arm) and B($arm) report identical 2G"; exit 1; }
done
echo "   control holds: source fields change the numerics in both arms"

echo "ENERGYRELEASE PROOF GREEN wall=$((SECONDS - T0))s"
