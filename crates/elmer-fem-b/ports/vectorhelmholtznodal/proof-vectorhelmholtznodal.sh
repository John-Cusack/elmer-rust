#!/bin/bash
# T021-B VectorHelmholtzNodal solver-entry proof: the Rust cdylib
# replaces VectorHelmholtzNodal.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `vectorhelmholtznodal_` + init entry);
#  2. two distinct .sif differentials (config A: segregated
#     component-wise solve over 2 steady iterations, i.e. the
#     repeated-invocation run through the SAVE state; config B:
#     monolithic solve over 2 steady iterations) are byte-identical to
#     the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (SAVE re-entry across processes);
#  4. the branch-liveness control (config B's monolithic branch)
#     produces different output from config A, and both arms agree on
#     it.
#
# Mesh note: this solver needs 3-D elements (nodal basis on tets; the
# curl-curl form even needs Whitney edge topology upstream). The
# PoissonDB triangular `square.grd` mesh cannot host it (2-D, no 3-D
# basis), so the proof uses the solver's own upstream 3-D mesh
# `shoebox_tets` from the Elmer checkout
# (fem/tests/VectorHelmholtzWaveguideNodal), copied verbatim — the same
# ElmerGrid-flow exemption the brief grants when the triangular mesh
# cannot host the solver.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout for the shoebox_tets mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-vectorhelmholtznodal"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/shoebox_tets
echo "== T021-B VectorHelmholtzNodal proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libvectorhelmholtznodal.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/VectorHelmholtzNodal.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T vectorhelmholtznodal_" "$WORK/port.syms" \
  || { echo "FAIL: vectorhelmholtznodal_ not the exported entry"; exit 1; }
grep -qx "T vectorhelmholtznodal_init_" "$WORK/port.syms" \
  || { echo "FAIL: vectorhelmholtznodal_init_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (upstream shoebox_tets 3-D mesh from the Elmer checkout)"
cp -r "$ELMER_FORK/fem/tests/VectorHelmholtzWaveguideNodal/shoebox_tets" "$WORK/"

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Check Keywords "Warn"

Header
  Mesh DB "." "shoebox_tets"
End

Simulation
  Max Output Level = 5
  Coordinate System = "Cartesian"
  Simulation Type = Steady
  Steady State Min Iterations = 2
  Steady State Max Iterations = 2
  Output Intervals(1) = 1
End

Body 1
  Equation = 1
  Material = 1
End

$ beta0=41.9332032192090
$ k0=52.3961255469564
$ kc=31.4159265358979
$ w=2*pi*(2.5e9)
$ mu0=4e-7*pi

Material 1
  Relative Permittivity = Real 1
End

Equation 1
  Active Solvers(2) = 1 2
  Angular Frequency = Real $w
End

Solver 1
  Equation = "VectorHelmholtzNodal"
  Procedure = "VectorHelmholtzNodal" "VectorHelmholtzNodal"

  Linear System Scaling = True
  Linear System Solver = "Iterative"
  Linear System Iterative Method = String "bicgstab2"
  BiCGstabl polynomial degree = Integer 4
  Linear System Preconditioning = String "ILU1"
  Linear System Max Iterations = Integer 1000
  Linear System Convergence Tolerance = 1.0e-7
  Linear system abort not converged = false

  Skip Compute Steady State Change = Logical True

  Nonlinear System Max Iterations = 1
End

Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Elfield
  Operator 1 = max
End

Boundary Condition 1
  Name = "Inlet"
  Target Boundaries(1) = 1
  Electric Robin Coefficient im = Real $ beta0
  Magnetic Boundary Load 2 = Variable Coordinate 1
    Real MATC "-2*beta0*k0/kc*sin(pi*tx/100e-3)"
End

Boundary Condition 2
  Name = "Sides"
  Target Boundaries(4) = 2 3 4 5

  Electric Robin Coefficient = Real $ -474315.9900340191
  Electric Robin Coefficient im = Real $ 474315.9900340191
End

Boundary Condition 3
  Name = "Outlet"
  Target Boundaries(1) = 6
  Absorbing BC = Logical True
End
EOF
echo "-- config B (config A plus the monolithic branch)"
python3 - "$WORK/caseA.sif" "$WORK/caseB.sif" <<'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
t = open(src).read()
t = t.replace("""  Procedure = "VectorHelmholtzNodal" "VectorHelmholtzNodal"
""", """  Procedure = "VectorHelmholtzNodal" "VectorHelmholtzNodal"

  Monolithic Solver = Logical True
""", 1)
open(dst, "w").write(t)
PYEOF
grep -q "Monolithic Solver = Logical True" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/VectorHelmholtzNodal.so"
cp "$PORT_SO" "$WORK/shadow-port/VectorHelmholtzNodal.so"
nm -D --defined-only "$WORK/shadow-port/VectorHelmholtzNodal.so" | grep -q "T vectorhelmholtznodal_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/VectorHelmholtzNodal.so" | grep -q "T vectorhelmholtznodal_init_" \
  || { echo "FAIL: staged port .so lost its init entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/shoebox_tets" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER SOLVER.*FINISHED AT|ELMER SOLVER.*STARTED AT" "$1" || true
}

echo "-- differential A (segregated, 2 steady iterations) + rerun (re-entry)"
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
    { echo "FAIL: rerun $f diverges"; diff "$WORK/a1.txt" "$WORK/a2.txt" | head; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (monolithic, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B must differ from A on both arms)"
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b1.txt"
diff -q "$WORK/a0.txt" "$WORK/b1.txt" > /dev/null 2>&1 && \
  { echo "FAIL: branch control produced no output delta"; exit 1; }
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: port branch control produced no output delta"; exit 1; }
echo "   branch control differs on both arms, arms agree with each other"

T1=$SECONDS
echo "VECTORHELMHOLTZNODAL PROOF GREEN (wall $((T1-T0))s)"
