#!/bin/bash
# T022-A Electrokinetics UDF proof: the Rust cdylib replaces
# Electrokinetics.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: the six UDF entries getjouleheat_,
#     helmholtz_smoluchowski[_123]_, helmholtz_smoluchowski_comp_);
#  2. three distinct .sif differentials exercise the solver (StatElec +
#     Poisson driving the UDFs as Real Procedures): config A (computed
#     Electric Field + EO Mobility + HS components 1/2 + Joule heat);
#     config B (constant Electric Field + Zeta Potential + tangential HS
#     + HS component 3 + Joule heat, with 2 steady iterations, i.e. the
#     repeated-invocation run); config C (config B with doubled Zeta
#     Potential: branch-liveness control) are byte-identical to the
#     reference except wall-clock lines/values;
#  3. a second full process run of config A is byte-identical
#     (SAVE re-entry through FirstTime/vacuumPerm/Nodes);
#  4. the liveness control (B vs C) produces different output, and both
#     arms agree on it.
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
WORK="$LANE/scratch/proof-electrokinetics"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A Electrokinetics proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran 13.3.0 required, got: $GFORTRAN_VERSION"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: nightly 1.100.0 (feaadeeac 2026-09-19) required, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libelectrokinetics.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/Electrokinetics.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the Electrokinetics module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/port.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in getjouleheat_ helmholtz_smoluchowski_ helmholtz_smoluchowski1_ helmholtz_smoluchowski2_ helmholtz_smoluchowski3_ helmholtz_smoluchowski_comp_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
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
  Coordinate System = Cartesian 2D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Max Output Level = 4
End
Constants
  Permittivity Of Vacuum = Real 8.8542e-12
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 2
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 3
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 4
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(2) = 1 2
  Electric Field = String computed
End
Material 1
  Density = Real 1.0e3
  Viscosity = Real 1.0e-3
  Relative Permittivity = Real 2.0
  Electric Conductivity = Real 1.0
End
Body Force 1
  Source = Variable Potential
    Real Procedure "Electrokinetics" "getJouleHeat"
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Potential = Real 0.0
  EO Mobility = Real 5.0e-8
  Slip = Variable Potential
    Real Procedure "Electrokinetics" "helmholtz_smoluchowski1"
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Potential = Real 1.0
  EO Mobility = Real 5.0e-8
  Slip = Variable Potential
    Real Procedure "Electrokinetics" "helmholtz_smoluchowski2"
End
Solver 1
  Equation = Stat Elec
  Procedure = "StatElecSolve" "StatElecSolver"
  Variable = Potential
  Variable DOFs = 1
  Calculate Electric Field = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = Slip Solver
  Procedure = "Poisson" "PoissonSolver"
  Variable = Slip
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
  Variable 3 = Slip
  Operator 3 = min
  Variable 4 = Slip
  Operator 4 = max
End
EOF
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    -e 's/Electric Field = String computed/Electric Field = String constant/' \
    -e 's/Relative Permittivity = Real 2.0/Relative Permittivity = Real 2.0\n  Electric Field 1 = Real 100.0\n  Electric Field 2 = Real 50.0/' \
    -e 's/EO Mobility = Real 5.0e-8/Zeta Potential = Real 0.05/' \
    -e 's/"helmholtz_smoluchowski1"/"helmholtz_smoluchowski"/' \
    -e 's/"helmholtz_smoluchowski2"/"helmholtz_smoluchowski3"/' \
    -e 's/Output Intervals = 1/Output Intervals = 1\n  Steady State Min Iterations = 2/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -q 'helmholtz_smoluchowski"' "$WORK/caseB.sif" || { echo "FAIL: caseB tangential HS not staged"; exit 1; }
grep -q 'helmholtz_smoluchowski3' "$WORK/caseB.sif" || { echo "FAIL: caseB HS3 not staged"; exit 1; }
sed -e 's/Zeta Potential = Real 0.05/Zeta Potential = Real 0.1/' \
    "$WORK/caseB.sif" > "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/Electrokinetics.so"
cp "$PORT_SO" "$WORK/shadow-port/Electrokinetics.so"
nm -D --defined-only "$WORK/shadow-port/Electrokinetics.so" | grep -q "T helmholtz_smoluchowski_comp_" \
  || { echo "FAIL: staged port .so lost its entries"; exit 1; }

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
norm() { # drop wall-clock lines/values for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" | \
    sed -E -e 's/(Assembly: \(s\))[0-9.eE+ -]+/\1 <t>/g' -e 's/(Solve: \(s\))[0-9.eE+ -]+/\1 <t>/g' || true
}
echo "-- differential A (computed E-field + EO Mobility + HS1/HS2 + Joule) + rerun (re-entry)"
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

echo "-- differential B (constant E-field + Zeta + tangential HS + HS3 + Joule, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- differential C (liveness control: doubled Zeta Potential)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"
norm "$WORK/runB-port/scalars.dat" > "$WORK/bc1.txt"
norm "$WORK/runC-port/scalars.dat" > "$WORK/bc2.txt"
if diff "$WORK/bc1.txt" "$WORK/bc2.txt" > /dev/null 2>&1; then
  echo "FAIL: liveness control produced identical output (B == C)"
  exit 1
fi
echo "   liveness control differs (B != C) on both arms"

echo "ELECTROKINETICS PROOF GREEN (wall: $((SECONDS - T0))s)"
