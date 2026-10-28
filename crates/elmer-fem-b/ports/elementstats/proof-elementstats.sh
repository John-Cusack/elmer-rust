#!/bin/bash
# T021-A ElementStats solver-entry proof: the Rust cdylib replaces
# ElementStats.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `elementstats_`);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines — config A (defaults: bulk + boundary
#     sweeps, no cohorts) and config B (`Create Histogram` with 5
#     intervals, `Skip Boundaries`, and a huge `Element Warn Size` that
#     trips the `Warn` path but not the `Fatal` one);
#  3. a second full process run of config B is byte-identical
#     (idempotent re-entry; the port holds no cross-call state — the
#     cohort table is allocated and released per call, as in the
#     original).
# Branch-liveness control: config B prints histogram sections config A
# lacks, in both arms, so the differential is not vacuous. Negative
# control: config C (`Stop Below Warn Size` with the same huge threshold)
# aborts both arms via `Fatal` with byte-identical stderr — the error
# path is identical too.
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
WORK="$LANE/scratch/proof-elementstats"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-A ElementStats proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libelementstats.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ElementStats.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ElementStats module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T elementstats_" "$WORK/port.syms" \
  || { echo "FAIL: elementstats_ not the exported entry"; exit 1; }
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
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Solver 1
  Equation = ElemStats
  Procedure = "ElementStats" "ElementStats"
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
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
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Solver 1
  Equation = ElemStats
  Procedure = "ElementStats" "ElementStats"
  Create Histogram = Logical True
  Histogram Intervals = Integer 5
  Skip Boundaries = Logical True
  Element Warn Size = Real 1.0e30
End
EOF
cat > "$WORK/caseC.sif" <<'EOF'
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
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Solver 1
  Equation = ElemStats
  Procedure = "ElementStats" "ElementStats"
  Element Warn Size = Real 1.0e30
  Stop Below Warn Size = Logical True
End
EOF

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ElementStats.so"
cp "$PORT_SO" "$WORK/shadow-port/ElementStats.so"
nm -D --defined-only "$WORK/shadow-port/ElementStats.so" | grep -q " T elementstats_" \
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
run_case_fail() { # <case> <arm> <rundir> — both arms must abort identically
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
  [ "$rc" -ne 0 ] || { echo "FAIL: $run($arm) unexpectedly succeeded"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differentials A + B + rerun of B"
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

echo "-- negative control C (Fatal in both arms, byte-identical)"
run_case_fail caseC orig runC-orig
run_case_fail caseC port runC-port
for f in stdout.txt stderr.txt; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: control C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   control C byte-identical aborts"

echo "-- branch-liveness control (histogram present in B, absent in A, both arms)"
for arm in orig port; do
  if ! grep -q "Histogram:" "$WORK/runB-$arm/stdout.txt"; then
    echo "FAIL: histogram missing in B($arm)"; exit 1;
  fi
  if grep -q "Histogram:" "$WORK/runA-$arm/stdout.txt"; then
    echo "FAIL: histogram unexpectedly present in A($arm)"; exit 1;
  fi
done
echo "   control holds: keyword branches change the reports in both arms"

echo "ELEMENTSTATS PROOF GREEN wall=$((SECONDS - T0))s"
