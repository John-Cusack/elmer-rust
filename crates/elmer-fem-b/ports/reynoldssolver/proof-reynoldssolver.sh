#!/bin/bash
# T021-B ReynoldsSolver proof: the Rust cdylib replaces ReynoldsSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: reynoldssolver_, reynoldssolver_init_,
#     reynoldspostprocess_, reynoldspostprocess_init_);
#  2. three distinct .sif differentials (Robin-BC config A; Flux-BC
#     config B with 2 steady iterations, i.e. the repeated-invocation
#     run; weakly-compressible config C) are byte-identical to the
#     reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch liveness: config C (compressibility branch) produces
#     different output from config A, and both arms agree on it.
#
# Mesh: the unit's own upstream structured mesh
# (fem/tests/reynolds_robin/rect.grd, byte-identical to the
# reynolds_flux one) built with the same `ElmerGrid 1 2 rect` flow the
# upstream runtest.cmake uses. Configs A/B are the upstream
# reynolds_robin / reynolds_flux cases verbatim plus an appended
# SaveScalars solver (documents FilmPressure min/max into scalars.dat
# for the byte-compare); config B additionally runs 2 steady iterations
# (SAVE re-entry in one process, mirroring the Poisson template's
# 2-iteration config); config C switches the compressibility branch.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config —
# both arms share it, so the differential is self-consistent).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-reynoldssolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/rect*
echo "== T021-B Reynolds proof: work=$WORK build=$ELMER_BUILD"
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
  *) echo "FAIL: need rustc 1.100.0-nightly"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libreynoldssolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ReynoldsSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ReynoldsSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in reynoldssolver_ reynoldssolver_init_ reynoldspostprocess_ reynoldspostprocess_init_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's upstream rect.grd)"
cp "$ELMER_FORK/fem/tests/reynolds_robin/rect.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 rect > /dev/null)

echo "-- configs (upstream cases + documented SaveScalars append)"
cp "$ELMER_FORK/fem/tests/reynolds_robin/case.sif" "$WORK/caseA.sif"
cp "$ELMER_FORK/fem/tests/reynolds_flux/case.sif" "$WORK/caseB.sif"
append_scalars() { # <sif>
cat >> "$1" <<'EOF'

Solver 4
  Equation = "SaveScalars"
  Procedure = File "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = FilmPressure
  Operator 1 = min
  Variable 2 = FilmPressure
  Operator 2 = max
End
EOF
}
 append_scalars "$WORK/caseA.sif"
 append_scalars "$WORK/caseB.sif"
 # Config B: 2 steady iterations in one process (SAVE re-entry). Elmer
 # stops the steady loop early once SteadyStateReached, so Min
 # Iterations is set to force the second same-process invocation.
 sed -i 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2\n  Steady State Min Iterations = 2/' "$WORK/caseB.sif"
 # Config C: same as A but with the Open Side boundary branch taken on
 # boundary 2 (exercises the OpenSide path in the boundary assembly,
 # dormant in both A and B).
cp "$WORK/caseA.sif" "$WORK/caseC.sif"
python3 - "$WORK/caseC.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
old = """Boundary Condition 2
  Name = "right"
  Target Boundaries(1) = 2
End
"""
new = """Boundary Condition 2
  Name = "right"
  Target Boundaries(1) = 2
  Open Side = Logical True
End
"""
assert old in t
t = t.replace(old, new, 1)
open(p, "w").write(t)
PYEOF
grep -c "Open Side" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ReynoldsSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/ReynoldsSolver.so"
nm -D --defined-only "$WORK/shadow-port/ReynoldsSolver.so" | grep -q "T reynoldssolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/rect" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Assembly \(s\)|Solve \(s\)|pointer init time|Time spent for restart|Factorization time" "$1" || true
}

echo "-- differential A (Robin branch) + rerun (re-entry)"
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

echo "-- differential B (Flux branch, 2 steady iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

 echo "-- differential C + branch liveness (Open Side branch)"
 run_case caseC orig runC-orig
 run_case caseC port runC-port
 for f in stdout.txt stderr.txt scalars.dat; do
   norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
   norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
   diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
     { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
 done
 echo "   config C byte-identical (stdout/stderr/scalars)"
 norm "$WORK/runA-port/stdout.txt" > "$WORK/liveA.txt"
 norm "$WORK/runC-port/stdout.txt" > "$WORK/liveC.txt"
 if diff "$WORK/liveA.txt" "$WORK/liveC.txt" > /dev/null 2>&1; then
   echo "FAIL: branch-liveness control identical (A vs C must differ)"
   exit 1
 fi
 echo "   branch-liveness control differs (A vs C) as required"

echo "REYNOLDS PROOF GREEN ($((SECONDS - T0))s wall)"
