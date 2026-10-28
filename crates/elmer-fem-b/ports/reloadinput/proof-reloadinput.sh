#!/bin/bash
# T021-B ReloadInput solver-entry proof: the Rust cdylib replaces ReloadInput.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `reloadinput_`);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines — config A takes the early-return path (no
#     ELMERSOLVER_REREADINFO trigger file), config B takes the reload path
#     (trigger file present; `LoadInputFile` re-reads the input);
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the wrapper SAVE state);
#  4. branch liveness: config B output differs from config A output
#     (the `LoadInputFile: Loading input file:` line plus changed scalars),
#     and both arms agree on each config.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial proof tree —
# both arms share it, so the differential is self-consistent) and the
# Elmer checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-reloadinput"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B ReloadInput proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran pin mismatch"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: rustc pin mismatch"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libreloadinput.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ReloadInput.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ReloadInput module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T reloadinput_" "$WORK/port.syms" \
  || { echo "FAIL: reloadinput_ not the exported entry"; exit 1; }
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
  Body Force = 1
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
  Body Force = 2
End
Equation 1
  Active Solvers(2) = 1 3
End
Material 1
  Density = 1.0
End
Body Force 1
  Source = Real 1.0
End
Body Force 2
  Source = Real -1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Potential = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Potential = Real 1.0
End
Solver 1
  Equation = Poisson
  Procedure = "Poisson" "PoissonSolver"
  Variable = Potential
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
End
Solver 3
  Equation = Reload
  Procedure = "ReloadInput" "ReloadInput"
  Exec Solver = after timestep
End
EOF
# Config B: different ValueList branches (2 steady iterations — in-process
# repeated invocation of the entry — plus a hotter Dirichlet value) and the
# reload path (trigger file staged at run time).
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
s = open(p).read()
s = s.replace("Potential = Real 1.0\nEnd", "Potential = Real 2.0\nEnd")
open(p, "w").write(s)
PYEOF
grep -q "Potential = Real 2.0" "$WORK/caseB.sif" \
  || { echo "FAIL: caseB branch edit did not apply"; exit 1; }
grep -q "Steady State Max Iterations = 2" "$WORK/caseB.sif" \
  || { echo "FAIL: caseB iteration edit did not apply"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ReloadInput.so"
cp "$PORT_SO" "$WORK/shadow-port/ReloadInput.so"
nm -D --defined-only "$WORK/shadow-port/ReloadInput.so" | grep -q "T reloadinput_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir> [trigger-content]
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  if [ "${4:-}" != "" ]; then
    printf '%s\n' "$4" > "$WORK/$run/ELMERSOLVER_REREADINFO"
  fi
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

echo "-- differential A (early-return path, no trigger) + rerun (re-entry)"
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
echo "   rerun byte-identical (re-entry)"

echo "-- differential B (reload path, 2 steady iterations, hotter BC)"
run_case caseB orig runB-orig "case.sif"
run_case caseB port runB-port "case.sif"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B output differs from A, both arms agree)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1 \
  && { echo "FAIL: liveness control: B stdout identical to A stdout"; exit 1; }
grep -q "LoadInputFile: Loading input file: case.sif" "$WORK/runB-port/stdout.txt" \
  || { echo "FAIL: liveness control: reload path left no trace in B stdout"; exit 1; }
diff "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1 \
  && { echo "FAIL: liveness control: B scalars identical to A scalars"; exit 1; }
echo "   B differs from A (reload trace + hotter scalars), arms agree"

echo "RELOADINPUT PROOF GREEN ($((SECONDS - T0))s wall)"
