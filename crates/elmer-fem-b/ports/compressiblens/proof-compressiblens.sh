#!/bin/bash
# T022-A CompressibleNS solver-entry proof: the Rust cdylib replaces
# CompressibleNS.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `compressiblens_`);
#  2. two distinct transient .sif differentials (body-force + bulk-
#     viscosity config A with 1 nonlinear iteration per step; no-force,
#     default-bulk-viscosity, bubble-enriched config B with 2 nonlinear
#     iterations) are byte-identical to the reference except wall-clock
#     lines/values;
#  3. a second full process run of config A against the port is
#     byte-identical (idempotent re-entry through the SAVE state);
#  4. the A/B branch variants produce different outputs on both arms
#     (branch-liveness control).
#
# Config mechanics: transient (2 BDF-1 steps) so the density/temperature
# block gets its mass term; the campaign Procedure form
# `"CompressibleNS" "CompressibleNS"` is used on both arms (library-vs-
# entry names differ only by case, resolved by the loader).
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
WORK="$LANE/scratch/proof-compressiblens"
CTARGET="${CARGO_TARGET_DIR:-$LANE/scratch/cargo-target-compressiblens}"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A CompressibleNS proof: work=$WORK build=$ELMER_BUILD"
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
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" CARGO_TARGET_DIR="$CTARGET" cargo +nightly build --release -q)
PORT_SO="$CTARGET/release/libcompressiblens.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/CompressibleNS.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the CompressibleNS module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T compressiblens_" "$WORK/port.syms" \
  || { echo "FAIL: compressiblens_ not the exported entry"; exit 1; }
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
  Simulation Type = Transient
  Timestepping Method = BDF
  BDF Order = 1
  Timestep Intervals = 2
  Timestep Sizes = 0.01
  Output Intervals = 1
  Max Output Level = 4
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
  Body Force = 1
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Equilibrium Temperature = Real 293.0
  Equilibrium Density = Real 1.21
  Viscosity = Real 1.82e-5
  Bulk Viscosity = Real 6.0e-6
  Specific Heat = Real 1005.0
  Heat Conductivity = Real 0.026
  Specific Heat Ratio = Real 1.4
End
Body Force 1
  Body Force 1 = Real 1.0
  Body Force 2 = Real -0.5
  Body Force 3 = Real 0.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  NSSol 1 = Real 0.0
  NSSol 2 = Real 0.0
  NSSol 3 = Real 0.0
  NSSol 4 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  NSSol 1 = Real 0.1
  NSSol 2 = Real 0.0
  NSSol 3 = Real 0.2
  NSSol 4 = Real 0.0
End
Solver 1
  Equation = "Compressible Flow"
  Procedure = "CompressibleNS" "CompressibleNS"
  Variable = -dofs 4 "NSSol"
  Element = "p:1"
  Nonlinear Iteration Method = String picard
  Nonlinear System Max Iterations = Integer 1
  Nonlinear System Convergence Tolerance = Real 1.0e-4
  Linear System Solver = Direct
  Linear System Direct Method = UMFPack
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = "NSSol 1"
  Operator 1 = min
  Variable 2 = "NSSol 1"
  Operator 2 = max
  Variable 3 = "NSSol 2"
  Operator 3 = min
  Variable 4 = "NSSol 2"
  Operator 4 = max
  Variable 5 = "NSSol 3"
  Operator 5 = min
  Variable 6 = "NSSol 3"
  Operator 6 = max
  Variable 7 = "NSSol 4"
  Operator 7 = min
  Variable 8 = "NSSol 4"
  Operator 8 = max
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
  Timestep Intervals = 2
  Timestep Sizes = 0.01
  Output Intervals = 1
  Max Output Level = 4
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
  Equilibrium Temperature = Real 300.0
  Equilibrium Density = Real 1.0
  Viscosity = Real 2.0e-5
  Specific Heat = Real 718.0
  Heat Conductivity = Real 0.03
  Specific Heat Ratio = Real 1.4
End
Boundary Condition 1
  Target Boundaries(1) = 1
  NSSol 1 = Real 0.0
  NSSol 2 = Real 0.05
  NSSol 3 = Real 0.0
  NSSol 4 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  NSSol 1 = Real 0.0
  NSSol 2 = Real 0.0
  NSSol 3 = Real 0.0
  NSSol 4 = Real 0.0
End
Solver 1
  Equation = "Compressible Flow"
  Procedure = "CompressibleNS" "CompressibleNS"
  Variable = -dofs 4 "NSSol"
  Element = "p:1 b:1"
  Bubbles in Global System = Logical True
  Nonlinear Iteration Method = String picard
  Nonlinear System Max Iterations = Integer 2
  Nonlinear System Convergence Tolerance = Real 1.0e-4
  Linear System Solver = Direct
  Linear System Direct Method = UMFPack
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = "NSSol 1"
  Operator 1 = min
  Variable 2 = "NSSol 1"
  Operator 2 = max
  Variable 3 = "NSSol 2"
  Operator 3 = min
  Variable 4 = "NSSol 2"
  Operator 4 = max
  Variable 5 = "NSSol 3"
  Operator 5 = min
  Variable 6 = "NSSol 3"
  Operator 6 = max
  Variable 7 = "NSSol 4"
  Operator 7 = min
  Variable 8 = "NSSol 4"
  Operator 8 = max
End
EOF
grep -c "NSSol" "$WORK/caseA.sif" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/CompressibleNS.so"
cp "$PORT_SO" "$WORK/shadow-port/CompressibleNS.so"
nm -D --defined-only "$WORK/shadow-port/CompressibleNS.so" | grep -q "T compressiblens_" \
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
norm() { # drop wall-clock lines; mask (s) timing values (structure kept)
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" | \
    sed -E -e 's/(Assembly: +\(s\) +)[0-9.+-]+/\1T/' \
           -e 's/(Solution: +\(s\) +)[0-9.+-]+/\1T/' || true
}

echo "-- differential A (body-force + bulk-viscosity branch) + rerun (re-entry)"
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
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" | head; exit 1; }
done
norm "$WORK/runA-port/stdout.txt" > "$WORK/r1.txt"
norm "$WORK/runA-port2/stdout.txt" > "$WORK/r2.txt"
diff "$WORK/r1.txt" "$WORK/r2.txt" > /dev/null 2>&1 || \
  { echo "FAIL: rerun stdout diverges"; diff "$WORK/r1.txt" "$WORK/r2.txt" | head; exit 1; }
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (no-force + bubbles + 2 nonlinear iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  if [ "$f" = stdout.txt ]; then
    norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
    norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
    diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
      { echo "FAIL: config B stdout diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
  else
    diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" > /dev/null 2>&1 || \
      { echo "FAIL: config B $f diverges"; diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" | head; exit 1; }
  fi
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A and B must differ on both arms)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/la.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/lb.txt"
diff -q "$WORK/la.txt" "$WORK/lb.txt" > /dev/null 2>&1 && \
  { echo "FAIL: A/B port outputs identical (dead control)"; exit 1; }
diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: A/B port scalars identical (dead control)"; exit 1; }
echo "   A/B outputs differ (live control)"

echo "COMPRESSIBLENS PROOF GREEN (T=$SECONDS s)"
