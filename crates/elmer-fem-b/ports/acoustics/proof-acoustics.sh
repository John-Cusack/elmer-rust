#!/bin/bash
# T022-A Acoustics solver-entry proof: the Rust cdylib replaces Acoustics.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `acousticssolver_`, `acousticssolver_init_`);
#  2. two distinct clean-room .sif differentials (driven-cavity config A;
#     Bubbles + different-frequency config B) are byte-identical to the
#     reference except wall-clock lines/values;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. config B takes a different code branch (Bubbles condensation +
#     different frequency) and produces different output, with both arms
#     agreeing on it (liveness control).
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
WORK="$LANE/scratch/proof-acoustics"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A Acoustics proof: work=$WORK build=$ELMER_BUILD"
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
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$LANE/scratch/cargo-target-acoustics}" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libacoustics.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/Acoustics.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the Acoustics module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T acousticssolver_" "$WORK/port.syms" \
  || { echo "FAIL: acousticssolver_ not the exported entry"; exit 1; }
grep -qx "T acousticssolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: acousticssolver_init_ not the exported entry"; exit 1; }
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
  Angular Frequency = Real 1000.0
  Max Output Level = 4
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
  Active Solvers(1) = 1
End
Material 1
  Specific Heat = 1005.0
  Specific Heat Ratio = 1.4
  Equilibrium Density = 1.2
  Equilibrium Temperature = 293.0
  Heat Conductivity = 0.0257
  Viscosity = 1.8e-5
End
Body Force 1
  Re Body Force 1 = Real 1.0
End
Boundary Condition 1
  Target Boundaries(4) = 1 2 3 4
  Re Velocity 1 = Real 0.0
  Im Velocity 1 = Real 0.0
  Re Velocity 2 = Real 0.0
  Im Velocity 2 = Real 0.0
End
Solver 1
  Equation = Acoustics
  Procedure = "Acoustics" "AcousticsSolver"
  Variable = Flow
  Variable Dofs = 8
  Bubbles = Logical False
  Linear System Solver = Direct
  Linear System Direct Method = Banded
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Flow
  Operator 1 = min
  Variable 2 = Flow
  Operator 2 = max
End
EOF
sed -e 's/Angular Frequency = Real 1000.0/Angular Frequency = Real 2000.0/' \
    -e 's/Bubbles = Logical False/Bubbles = Logical True/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -q "Bubbles = Logical True" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/Acoustics.so"
cp "$PORT_SO" "$WORK/shadow-port/Acoustics.so"
nm -D --defined-only "$WORK/shadow-port/Acoustics.so" | grep -q "T acousticssolver_" \
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
norm() { # drop wall-clock lines, mask per-iteration timing values
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT" "$1" \
    | sed -E 's/(: +\(s\)) +[^ ]+/\1/' || true
}

echo "-- differential A (driven cavity) + rerun (re-entry)"
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

echo "-- differential B (Bubbles + different frequency)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness control (B differs from A on both arms)"
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b0.txt"
if diff -q "$WORK/a0.txt" "$WORK/b0.txt" > /dev/null 2>&1; then
  echo "FAIL: config B output identical to config A (dead control)"
  exit 1
fi
echo "   liveness confirmed (B output differs from A)"

echo "ACOUSTICS PROOF GREEN (wall $((SECONDS - T0))s)"
