#!/bin/bash
# T021-B SunAngle solver-entry proof: the Rust cdylib replaces SunAngle.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `sunanglesolver_`, `sunanglesolver_init_`);
#  2. two distinct clean-room .sif differentials exercise the solver
#     (config A: serial search algo with Test Height topography;
#     config B: `Enforce Parallel` search algo, different Test Angle and
#     heights, 2 steady iterations i.e. the repeated-invocation run) are
#     byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch liveness: config B output differs from config A, and both
#     arms agree on it.
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
WORK="$LANE/scratch/proof-sunangle"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B SunAngle proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran pin is 13.3.0, got: $GFORTRAN_VERSION"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: rustc pin is 1.100.0-nightly, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libsunangle.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- clippy clean"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release 2>&1 | tee "$WORK/clippy.txt" | grep -i "warning" && { echo "FAIL: clippy warnings"; exit 1; } || true)

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/SunAngle.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the SunAngle module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T sunanglesolver_" "$WORK/port.syms" \
  || { echo "FAIL: sunanglesolver_ not an exported entry"; exit 1; }
grep -qx "T sunanglesolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: sunanglesolver_init_ not an exported entry"; exit 1; }
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
  Steady State Max Iterations = 2
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
End
Body 2
  Equation = 1
  Material = 2
End
Body 3
  Equation = 1
  Material = 3
End
Body 4
  Equation = 1
  Material = 4
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Test Height = Real 0.0
End
Material 2
  Test Height = Real 3.0
End
Material 3
  Test Height = Real 7.0
End
Material 4
  Test Height = Real 1.5
End
Solver 1
  Equation = Sun Angle
  Procedure = "SunAngle" "SunAngleSolver"
  Variable = CriticalTan
  Variable DOFs = 1
  Test Angle = Real 0.6
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = CriticalTan
  Operator 1 = min
  Variable 2 = CriticalTan
  Operator 2 = max
End
EOF
# Config B: parallel-search branch (Enforce Parallel), different angle and
# topography, 2 steady iterations (repeated invocation).
sed -e 's/Test Angle = Real 0.6/Test Angle = Real 2.2\n  Enforce Parallel = Logical True/' \
    -e 's/Test Height = Real 0.0/Test Height = Real 4.0/' \
    -e 's/Test Height = Real 3.0/Test Height = Real 0.5/' \
    -e 's/Test Height = Real 7.0/Test Height = Real 2.0/' \
    -e 's/Test Height = Real 1.5/Test Height = Real 9.0/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -c "Enforce Parallel" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/SunAngle.so"
cp "$PORT_SO" "$WORK/shadow-port/SunAngle.so"
nm -D --defined-only "$WORK/shadow-port/SunAngle.so" | grep -q "T sunanglesolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|SOLVER FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (serial search) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  diff -q "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
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
    { echo "FAIL: rerun $f diverges (SAVE re-entry)"; diff "$WORK/a1.txt" "$WORK/a2.txt" | head; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (Enforce Parallel branch, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
  diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (B differs from A; arms agree B differs)"
norm "$WORK/runA-port/scalars.dat" > "$WORK/sa.txt"
norm "$WORK/runB-port/scalars.dat" > "$WORK/sb.txt"
diff "$WORK/sa.txt" "$WORK/sb.txt" > /dev/null 2>&1 && \
  { echo "FAIL: config B output identical to A (dead branch)"; exit 1; }
norm "$WORK/runB-orig/scalars.dat" > "$WORK/sb0.txt"
diff "$WORK/sb0.txt" "$WORK/sb.txt" > /dev/null 2>&1 || \
  { echo "FAIL: arms disagree on the branch output"; exit 1; }
echo "   branch live: B differs from A, both arms agree"

echo "SUNANGLE PROOF GREEN ($((SECONDS - T0))s wall)"
