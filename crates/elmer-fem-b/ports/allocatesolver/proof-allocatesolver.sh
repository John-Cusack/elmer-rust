#!/bin/bash
# T021-A AllocateSolver solver-entry proof: the Rust cdylib replaces
# AllocateSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `allocatesolver_` + `allocatesolver_init_`);
#  2. three distinct .sif differentials are byte-identical to the
#     reference except wall-clock lines — config A (default: `_init`
#     appends `Exec Solver = never` at the list tail, after the
#     setup-time `AddExecWhenFlag` cached `ALWAYS`, so the main entry
#     still runs once and prints its `Info` line); config B (transient,
#     2 timesteps: the main entry runs twice = the repeated-invocation
#     run); config C (`Exec Solver = never` in the .sif: setup-time
#     `NEVER`, the main entry never runs);
#  3. a second full process run of config B is byte-identical
#     (idempotent re-entry; the port holds no SAVE state, as the
#     original holds none).
# Branch-liveness control: the AllocateSolver `Info` line prints 1x in
# config A, 2x in config B, 0x in config C, in both arms — the dispatch
# toggle changes behavior identically in both arms, so the differential
# is not vacuous.
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
WORK="$LANE/scratch/proof-allocatesolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-A AllocateSolver proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/liballocatesolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/AllocateSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the AllocateSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T allocatesolver_" "$WORK/port.syms" \
  || { echo "FAIL: allocatesolver_ not the exported entry"; exit 1; }
grep -qx "T allocatesolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: allocatesolver_init_ not exported"; exit 1; }
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
Solver 1
  Equation = AllocateStuff
  Procedure = "AllocateSolver" "AllocateSolver"
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "square"
End
Simulation
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
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Solver 1
  Equation = AllocateStuff
  Procedure = "AllocateSolver" "AllocateSolver"
  Exec Solver = always
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
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Solver 1
  Equation = AllocateStuff
  Procedure = "AllocateSolver" "AllocateSolver"
  Exec Solver = never
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
rm -f "$WORK/shadow-port/AllocateSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/AllocateSolver.so"
nm -D --defined-only "$WORK/shadow-port/AllocateSolver.so" | grep -q " T allocatesolver_" \
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

echo "-- differential A (default: _init appends never, main still runs once)"
run_case caseA orig runA-orig
run_case caseA port runA-port
for f in stdout.txt stderr.txt; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
done
echo "   config A byte-identical (stdout/stderr)"

echo "-- differential B (transient, 2 timesteps: in-process repetition) + rerun"
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseB port runB-port2
for f in stdout.txt stderr.txt; do
  norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
  diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
  norm "$WORK/runB-port2/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   config B byte-identical (stdout/stderr) + rerun byte-identical"

echo "-- differential C (Exec Solver = never: main entry skipped)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt; do
  norm "$WORK/runC-orig/$f" > "$WORK/c0.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c1.txt"
  diff "$WORK/c0.txt" "$WORK/c1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c0.txt" "$WORK/c1.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr)"

echo "-- branch-liveness control (Info 1x in A, 2x in B, 0x in C, both arms)"
for arm in orig port; do
  na=$(grep -c "This solver does nothing" "$WORK/runA-$arm/stdout.txt")
  nb=$(grep -c "This solver does nothing" "$WORK/runB-$arm/stdout.txt")
  nc=$(grep -c "This solver does nothing" "$WORK/runC-$arm/stdout.txt" || true)
  [ "$na" = 1 ] || { echo "FAIL: A($arm) Info count $na, want 1"; exit 1; }
  [ "$nb" = 2 ] || { echo "FAIL: B($arm) Info count $nb, want 2"; exit 1; }
  [ "$nc" = 0 ] || { echo "FAIL: C($arm) Info count $nc, want 0"; exit 1; }
done
echo "   liveness control holds in both arms"

echo "ALLOCATESOLVER PROOF GREEN wall=$((SECONDS - T0))s"
