#!/bin/bash
# T021-B VelocityPrecond solver-entry proof: the Rust cdylib replaces
# VelocityPrecond.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `velocityprecond_` + `velocityprecond_init0_`);
#  2. two distinct .sif differentials exercise the solver: config A takes
#     the Init0 defaults (`ListAddNewString 'Variable' 'V'` fires, the
#     framework creates variable `v`); config B presets
#     `Variable = MyVar` (the `New` variant keeps the preset, the
#     framework creates `myvar` instead) — different ValueList branches;
#  3. a second full process run of config A is byte-identical
#     (re-entry; the unit holds no SAVE state);
#  4. branch-liveness: config B output differs from config A, and both
#     arms agree on it.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config) and
# the Elmer checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="/home/john/firstmate/data/rs-elmer-t021b/scratch/proof-velocityprecond"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B VelocityPrecond proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: expected gfortran 13.3.0, got: $GFORTRAN_VERSION"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: expected rustc 1.100.0-nightly, got: $RUSTC_VERSION"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libvelocityprecond.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/VelocityPrecond.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the VelocityPrecond module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T velocityprecond_" "$WORK/port.syms" \
  || { echo "FAIL: velocityprecond_ not exported"; exit 1; }
grep -qx "T velocityprecond_init0_" "$WORK/port.syms" \
  || { echo "FAIL: velocityprecond_init0_ not exported"; exit 1; }
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
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  V = Real 0.0
End
Solver 1
  Equation = VelocityPrecond
  Procedure = "VelocityPrecond" "VelocityPrecond"
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = V
  Operator 1 = min
End
EOF
sed -e 's/  Variable 1 = V/  Variable 1 = MyVar/' \
    -e 's/^Solver 1$/Solver 1\n  Variable = MyVar/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -q "Variable = MyVar" "$WORK/caseB.sif" \
  || { echo "FAIL: caseB variant not constructed"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/VelocityPrecond.so"
cp "$PORT_SO" "$WORK/shadow-port/VelocityPrecond.so"
nm -D --defined-only "$WORK/shadow-port/VelocityPrecond.so" | grep -q "T velocityprecond_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/VelocityPrecond.so" | grep -q "T velocityprecond_init0_" \
  || { echo "FAIL: staged port .so lost its init0 entry"; exit 1; }

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

echo "-- differential A (Init0 defaults: Variable defaults to V) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  cmp -s "$WORK/runA-orig/$f" "$WORK/runA-port/$f" \
    || { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stderr.txt scalars.dat; do
  cmp -s "$WORK/runA-port/$f" "$WORK/runA-port2/$f" \
    || { echo "FAIL: rerun $f diverges"; exit 1; }
done
norm "$WORK/runA-port/stdout.txt" > "$WORK/a2.txt"
norm "$WORK/runA-port2/stdout.txt" > "$WORK/a3.txt"
diff "$WORK/a2.txt" "$WORK/a3.txt" > /dev/null 2>&1 || \
  { echo "FAIL: rerun stdout diverges"; diff "$WORK/a2.txt" "$WORK/a3.txt" | head; exit 1; }
echo "   rerun byte-identical (re-entry)"

echo "-- differential B (preset Variable = MyVar: ListAddNew keeps preset) + branch liveness"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  if [ "$f" = "stdout.txt" ]; then
    norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
    norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
    diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
      { echo "FAIL: config B stdout diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
  else
    cmp -s "$WORK/runB-orig/$f" "$WORK/runB-port/$f" \
      || { echo "FAIL: config B $f diverges"; diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" | head; exit 1; }
  fi
done
echo "   config B byte-identical (stdout/stderr/scalars)"
if cmp -s "$WORK/a0.txt" "$WORK/b0.txt"; then
  echo "FAIL: branch-liveness control failed (A and B stdout identical)"
  exit 1
fi
grep -q "myvar" "$WORK/b0.txt" || { echo "FAIL: config B missing preset-variable trace"; exit 1; }
echo "   branch-liveness control differs (A creates v, B creates myvar), both arms agree"

echo "VELOCITYPRECOND PROOF GREEN (wall $((SECONDS - T0))s)"
