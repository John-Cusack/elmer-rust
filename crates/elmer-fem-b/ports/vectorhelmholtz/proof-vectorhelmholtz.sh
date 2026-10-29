#!/bin/bash
# T021-B VectorHelmholtz solver-entry proof: the Rust cdylib replaces
# VectorHelmholtz.so (seven entries: the solver, its init/init0, the
# dummy placeholder, and the CalcFields solver with its init/init0).
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality over all seven `T` entries);
#  2. two distinct clean-room .sif differentials (config A: edge-basis
#     waveguide solve + SaveScalars over 2 steady iterations, i.e. the
#     repeated-invocation run through the SAVE state; config B: config A
#     plus the CalcFields solver with elemental fields and the energy
#     functional) are byte-identical to the reference except wall-clock
#     lines;
#  3. a second full process run of config A is byte-identical
#     (SAVE re-entry across processes);
#  4. the branch-liveness control (config B's CalcFields branch)
#     produces different output from config A, and both arms agree on it.
#
# Mesh note: this solver needs a real 3-D Whitney-edge VectorHelmholtz
# solve (edge elements, Piola transform). The PoissonDB triangular
# `square.grd` mesh cannot host it (2-D, no edge topology, no Hcurl
# basis), so the proof uses the solver's own upstream 3-D mesh
# `shoebox_tets` from the Elmer checkout
# (fem/tests/VectorHelmholtzWaveguide), copied verbatim — the same
# ElmerGrid-flow exemption the vectorhelmholtzlumping port uses.
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
WORK="$LANE/scratch/proof-vectorhelmholtz"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/shoebox_tets
echo "== T021-B VectorHelmholtz proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libvectorhelmholtz.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/VectorHelmholtz.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in vectorhelmholtzsolver_ vectorhelmholtzsolver_init_ \
    vectorhelmholtzsolver_init0_ vectorhelmholtz_dummy_ \
    vectorhelmholtzcalcfields_ vectorhelmholtzcalcfields_init_ \
    vectorhelmholtzcalcfields_init0_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- mesh (upstream shoebox_tets 3-D mesh from the Elmer checkout)"
cp -r "$ELMER_FORK/fem/tests/VectorHelmholtzWaveguide/shoebox_tets" "$WORK/"

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
  Steady State Max Iterations = 2
End

$ mu0 = 4e-7*pi
$ eps0 = 8.854187817e-12
$ c = 1/sqrt(mu0*eps0)
$ a = 0.1
$ w = 2*pi*(2.5e9)
$ k0 = w/c
$ kc = 1*pi/a
$ beta0 = sqrt(k0^2-kc^2)

Material 1
  Relative Permittivity = Real 1
  Electric Conductivity = Real 0.0
End

Equation 1
  Active Solvers(2) = 1 2
  Angular Frequency = Real $w
End

Solver 1
  Equation = "VectorHelmholtz"
  Use Piola Transform = Logical True
  Procedure = "VectorHelmholtz" "VectorHelmholtzSolver"
  Linear System Solver = String "Iterative"
  Linear System Iterative Method = String "bicgstabl"
  BiCGstabl polynomial degree = Integer 4
  Linear System Preconditioning Damp Coefficient = Real 0.0
  Linear System Preconditioning Damp Coefficient im = Real -1
  Mass-proportional Damping = True
  Linear System Preconditioning = String "ILU0"
  Linear System ILUT Tolerance = Real 3e-3
  Linear System Max Iterations = Integer 4000
  Linear System Convergence Tolerance = 1.0e-7
  Linear System Abort not Converged = Logical False
  Steady State Convergence Tolerance = 1e-09
  Linear System Residual Output = 10
End

Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = File "scalars.dat"
  Variable 1 = String "e"
  Operator 1 = String "norm"
End

Boundary Condition 1
  Target Boundaries(1) = 1
  Electric Robin Coefficient im = Real $ beta0
  Magnetic Boundary Load 2 = Variable Coordinate 1
    Real MATC "-2*beta0*k0/kc*sin(kc*tx)"
End

Boundary Condition 2
  Target Boundaries(5) = 2 3 4 5 6
  Electric Robin Coefficient = Real $ -474315.9900340191
  Electric Robin Coefficient im = Real $ 474315.9900340191
End
EOF
echo "-- config B (config A plus the CalcFields solver branch)"
python3 - "$WORK/caseA.sif" "$WORK/caseB.sif" <<'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()
assert s.count("  Active Solvers(2) = 1 2\n") == 1
s = s.replace("  Active Solvers(2) = 1 2\n", "  Active Solvers(3) = 1 2 3\n", 1)
solver2 = """Solver 2
  Equation = "calcfields"
  Procedure = "VectorHelmholtz" "VectorHelmholtzCalcFields"
  Calculate Elemental Fields = Logical True
  Calculate Magnetic Field Strength = Logical True
  Calculate Magnetic Flux Density = Logical True
  Calculate Poynting vector = Logical True
  Calculate Electric field = Logical True
  Calculate Energy Functional = Logical True
  Steady State Convergence Tolerance = 1
  Linear System Solver = "Iterative"
  Linear System Preconditioning = None
  Linear System Residual Output = 10
  Linear System Max Iterations = 5000
  Linear System Iterative Method = CG
  Linear System Convergence Tolerance = 1.0e-9
End

"""
assert s.count('Solver 2\n  Equation = "SaveScalars"') == 1
s = s.replace('Solver 2\n  Equation = "SaveScalars"',
              solver2 + 'Solver 3\n  Equation = "SaveScalars"', 1)
open(dst, "w").write(s)
PYEOF
grep -c "VectorHelmholtzCalcFields" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/VectorHelmholtz.so"
cp "$PORT_SO" "$WORK/shadow-port/VectorHelmholtz.so"
nm -D --defined-only "$WORK/shadow-port/VectorHelmholtz.so" | grep -q "T vectorhelmholtzsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

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
norm() { grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Assembly \(s\)|Solve \(s\)|pointer init time|Time spent for restart|Factorization time" "$1" || true; }

echo "-- differential A (waveguide solve, 2 steady iterations) + rerun (re-entry)"
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

echo "-- differential B (CalcFields elemental + energy functional, 2 iterations)"
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
echo "VECTORHELMHOLTZ PROOF GREEN (wall $((T1-T0))s)"
