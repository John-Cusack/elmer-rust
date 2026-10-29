#!/bin/bash
# T021-B WPotentialSolver port proof: the Rust cdylib replaces WPotentialSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `wsolve_`, `wsolve_dummy_`, `wsolve_init0_`);
#  2. three distinct .sif differentials (massive-coil config A;
#     massive-coil + Compute J Vector config B; stranded-coil +
#     Without RotM + Compute J Vector config C) are byte-identical to
#     the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (SAVE re-entry through the wrapper-side and Rust-side state);
#  4. the configs produce pairwise-different outputs (branch-liveness
#     control: coil dispatch, J-vector branch, NoRotM branch), and both
#     arms agree on each.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent). The mesh is
# this solver's own upstream hex mesh (fem/tests/circuits_transient_massive/6481,
# 8-node bricks), copied as-is: a coil-free setup segfaults inside the
# reference itself (AddComponentsToBodyLists on an unallocated Components
# array), so every config carries a Component block.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-wpotentialsolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/mesh6464
echo "== T021-B WPotential proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libwpotentialsolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/WPotentialSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the WPotentialSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in wsolve_ wsolve_dummy_ wsolve_init0_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- mesh (this solver's upstream hex mesh, copied as-is)"
cp -r "$ELMER_FORK/fem/tests/circuits_transient_massive/6481" "$WORK/mesh6464"
rm -rf "$WORK/mesh6464/dat"

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "mesh6464"
End
Simulation
  Max Output Level = 5
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 2
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
  Body Force = 2
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Electric Conductivity = 36e6
End
Body Force 1
  Source = Real 1.0
End
Body Force 2
  Source = Real -1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  W = Real 0.0
End
Boundary Condition 2
  Target Boundaries(4) = 2 3 4 5
  W = Real 1.0
End
Component 1
  Name = String "A1"
  Master Bodies(1) = Integer 3
  Coil Type = String massive
End
Solver 1
  Equation = Wire direction
  Procedure = "WPotentialSolver" "Wsolve"
  Variable = W
  Linear System Solver = Iterative
  Linear System Iterative Method = CG
  Linear System Max Iterations = 5000
  Linear System Convergence Tolerance = 1.0e-10
  Linear System Abort Not Converged = True
  Linear System Residual Output = 1000
End
Solver 2
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = W
  Operator 1 = min
  Variable 2 = W
  Operator 2 = max
  Variable 3 = W Potential
  Operator 3 = min
  Variable 4 = W Potential
  Operator 4 = max
End
EOF
sed -e 's/Source = Real 1.0/Source = Real 2.5/' \
    -e 's/Source = Real -1.0/Source = Real -0.5/' \
    -e 's/W = Real 1.0/W = Real 2.0/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
t = t.replace("""  Linear System Residual Output = 1000
End
Solver 2""", """  Linear System Residual Output = 1000
  Compute J Vector = Logical True
End
Solver 2""", 1)
t = t.replace("""  Variable 4 = W Potential
  Operator 4 = max""", """  Variable 4 = W Potential
  Operator 4 = max
  Variable 5 = J Vector E
  Operator 5 = min
  Variable 6 = J Vector E
  Operator 6 = max""", 1)
open(p, "w").write(t)
PYEOF
grep -q "Compute J Vector" "$WORK/caseB.sif" || { echo "FAIL: caseB patch missed"; exit 1; }
sed -e 's/Source = Real 1.0/Source = Real 0.3/' \
    -e 's/Source = Real -1.0/Source = Real 0.7/' \
    -e 's/W = Real 1.0/W = Real 0.5/' \
    -e 's/Coil Type = String massive/Coil Type = String stranded/' \
    "$WORK/caseB.sif" > "$WORK/caseC.sif"
python3 - "$WORK/caseC.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
t = t.replace("""  Coil Type = String stranded""", """  Coil Type = String stranded
  Without RotM = Logical True""", 1)
open(p, "w").write(t)
PYEOF
grep -q "Without RotM" "$WORK/caseC.sif" || { echo "FAIL: caseC patch missed"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/WPotentialSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/WPotentialSolver.so"
nm -D --defined-only "$WORK/shadow-port/WPotentialSolver.so" | grep -q "T wsolve_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/mesh6464" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER.*FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (massive coil) + rerun (re-entry)"
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

echo "-- differential B (massive coil + J Vector)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- differential C (stranded coil + Without RotM + J Vector)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (configs disagree with each other)"
cmp -s "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat" \
  && { echo "FAIL: configs A and B agree (control dead)"; exit 1; }
cmp -s "$WORK/runA-orig/scalars.dat" "$WORK/runC-orig/scalars.dat" \
  && { echo "FAIL: configs A and C agree (control dead)"; exit 1; }
cmp -s "$WORK/runB-orig/scalars.dat" "$WORK/runC-orig/scalars.dat" \
  && { echo "FAIL: configs B and C agree (control dead)"; exit 1; }
echo "   A/B/C outputs pairwise differ, both arms agree on each"

echo "WPOTENTIAL PROOF GREEN ($((SECONDS - T0))s wall)"
