#!/bin/bash
# T021-B MagnetoDynamics2D solver-entry proof: the Rust cdylib replaces
# MagnetoDynamics2D.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: the seven entries);
#  2. two distinct .sif differentials (config A: transient stranded-coil
#     machine with MagnetoDynamics2D + Bsolver; config B: harmonic
#     stranded-coil machine with MagnetoDynamics2DHarmonic + Bsolver)
#     are byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch liveness: config A with `Handle Assembly = Logical True`
#     takes the handle-assembly branch (Warn + revert in _init, different
#     assembly Info in the entry), differs from config A, and both arms
#     agree on it.
#
# The mesh is this family's upstream magnetodynamics mesh
# (fem/tests/circuits2D_{transient,harmionic}_stranded/1381, copied as-is
# per the WPotentialSolver upstream-hex-mesh precedent), trimmed only in
# timestep count for proof speed.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-magnetodynamics2d"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/1381t "$WORK"/1381h "$WORK"/sift "$WORK"/sifh
echo "== T021-B MagnetoDynamics2D proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: need rustc +nightly 1.100.0-nightly"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libmagnetodynamics2d.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/MagnetoDynamics2D.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the MagnetoDynamics2D module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in magnetodynamics2d_ magnetodynamics2d_init_ magnetodynamics2dharmonic_ magnetodynamics2dharmonic_init_ magnetodynamics2dharmonic_init0_ bsolver_ bsolver_init_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not the exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- meshes (this family's upstream magnetodynamics meshes, copied as-is)"
cp -r "$ELMER_FORK/fem/tests/circuits2D_transient_stranded/1381" "$WORK/1381t"
cp -r "$ELMER_FORK/fem/tests/circuits2D_harmonic_stranded/1381" "$WORK/1381h"
rm -rf "$WORK/1381t/dat" "$WORK/1381h/dat"

echo "-- configs (clean-room .sif files below)"
mkdir -p "$WORK/sift" "$WORK/sifh"
cp "$ELMER_FORK/fem/tests/circuits2D_transient_stranded/sif/1381-circuits-components.definitions" "$WORK/sift/"
cp "$ELMER_FORK/fem/tests/circuits2D_harmonic_stranded/sif/1381-circuits-components.definitions" "$WORK/sifh/"
# Config A: transient stranded machine, trimmed to 2 intervals for speed.
sed -e 's/Timestep Intervals(1) =.*/Timestep Intervals(1) = 2/' \
  "$ELMER_FORK/fem/tests/circuits2D_transient_stranded/sif/1381-components-stranded.sif" \
  > "$WORK/caseA.sif"
grep -q "Timestep Intervals(1) = 2" "$WORK/caseA.sif" || { echo "FAIL: trim did not apply"; exit 1; }
# Config B: harmonic stranded machine, as upstream.
cp "$ELMER_FORK/fem/tests/circuits2D_harmonic_stranded/sif/1381-components-stranded.sif" "$WORK/caseB.sif"
# Liveness: config A variant forcing the Handle Assembly branch in Solver 2.
python3 - "$WORK/caseA.sif" > "$WORK/caseA-live.sif" <<'PYEOF'
import sys
src = open(sys.argv[1]).read()
anchor = '   Procedure = "MagnetoDynamics2D" "MagnetoDynamics2D"\n'
assert anchor in src, "Solver 2 anchor missing"
src = src.replace(anchor, anchor + '   Handle Assembly = Logical True\n', 1)
open(sys.argv[1] + "-tmp", "w").write(src)
PYEOF
mv "$WORK/caseA.sif-tmp" "$WORK/caseA-live.sif"
grep -q "Handle Assembly" "$WORK/caseA-live.sif" || { echo "FAIL: liveness variant missing"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/MagnetoDynamics2D.so"
cp "$PORT_SO" "$WORK/shadow-port/MagnetoDynamics2D.so"
nm -D --defined-only "$WORK/shadow-port/MagnetoDynamics2D.so" | grep -q "T magnetodynamics2d_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir> <meshdir> <sifdir>
  local case=$1 arm=$2 run=$3 meshdir=$4 sifdir=$5
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/$meshdir" "$WORK/$run/1381"
  cp -r "$WORK/$sifdir" "$WORK/$run/sif"
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; tail -n 20 "$WORK/$run/stderr.txt"; exit 1; }
}
norm() { grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Assembly \(s\)|Solve \(s\)|pointer init time|Time spent for restart|Factorization time" "$1" || true; }
cmp_scalars() { # <run1> <run2> <label>: byte-compare scalars.dat iff produced
  if [ -e "$WORK/$1/scalars.dat" ] || [ -e "$WORK/$2/scalars.dat" ]; then
    diff "$WORK/$1/scalars.dat" "$WORK/$2/scalars.dat" > /dev/null 2>&1 || \
      { echo "FAIL: $3 scalars.dat diverges"; diff "$WORK/$1/scalars.dat" "$WORK/$2/scalars.dat" | head; exit 1; }
    echo "   $3 scalars.dat byte-identical"
  else
    echo "   $3 produced no scalars.dat (SaveScalars prints only; stdout compared)"
  fi
}

echo "-- differential A (transient MagnetoDynamics2D + Bsolver) + rerun (re-entry)"
run_case caseA orig runA-orig 1381t sift
run_case caseA port runA-port 1381t sift
run_case caseA port runA-port2 1381t sift
for f in stderr.txt; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
cmp_scalars runA-orig runA-port "config A"
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
cmp_scalars runA-port runA-port2 "rerun"
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (harmonic MagnetoDynamics2DHarmonic + Bsolver)"
run_case caseB orig runB-orig 1381h sifh
run_case caseB port runB-port 1381h sifh
for f in stdout.txt stderr.txt; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
cmp_scalars runB-orig runB-port "config B"
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (Handle Assembly variant agrees, and differs from A)"
run_case caseA-live orig runL-orig 1381t sift
run_case caseA-live port runL-port 1381t sift
for f in stdout.txt stderr.txt; do
  norm "$WORK/runL-orig/$f" > "$WORK/l1.txt"
  norm "$WORK/runL-port/$f" > "$WORK/l2.txt"
  diff "$WORK/l1.txt" "$WORK/l2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: liveness $f diverges"; diff "$WORK/l1.txt" "$WORK/l2.txt" | head; exit 1; }
done
cmp_scalars runL-orig runL-port "liveness"
norm "$WORK/runA-port/stdout.txt" > "$WORK/la.txt"
norm "$WORK/runL-port/stdout.txt" > "$WORK/lb.txt"
if diff "$WORK/la.txt" "$WORK/lb.txt" > /dev/null 2>&1; then
  echo "FAIL: liveness variant output identical to config A (branch not live)"; exit 1;
fi
grep -qi "handle assembly" "$WORK/runL-port/stdout.txt" \
  || { echo "FAIL: liveness branch left no trace in output"; exit 1; }
echo "   liveness variant differs from A, both arms agree on it"

echo "MAGNETODYNAMICS2D PROOF GREEN (wall: $((SECONDS - T0))s)"
