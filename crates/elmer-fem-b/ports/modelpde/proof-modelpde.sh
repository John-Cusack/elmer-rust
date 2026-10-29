#!/bin/bash
# T021-B ModelPDE solver-entry proof: the Rust cdylib replaces ModelPDE.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `advdiffsolver_`);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines: config A (steady, convection + source +
#     Robin/flux BCs, nonlinear maxiter 1) and config B (transient BDF,
#     reaction + time-derivative coeffs, no convection keywords,
#     Integ Test report branch, Nonlinear System Max Iterations = 2);
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the wrapper SAVE state);
#  4. branch liveness: config B output differs from config A output, and
#     both arms agree on each config.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-modelpde"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B ModelPDE proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libmodelpde.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ModelPDE.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ModelPDE module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T advdiffsolver_" "$WORK/port.syms" \
  || { echo "FAIL: advdiffsolver_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Check Keywords "Warn"
Header
  Mesh DB "." "square"
End
Simulation
  Max Output Level = 5
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
  Active Solvers(1) = 1
End
Material 1
  Convection Velocity 1 = Real 1.0
  Convection Velocity 2 = Real 0.0
  diffusion coefficient = Real 1.0
  convection coefficient = Real 9.0
  time derivative coefficient = Real 0.0
End
Body Force 1
  field source = Real 1.0
End
Body Force 2
  field source = Real -1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Field = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 2
  Robin Coefficient = Real 10.0
  External Field = Real 5.0
End
Boundary Condition 3
  Target Boundaries(2) = 3 4
  Field Flux = Real 10.0
End
Solver 1
  Equation = "ModelPDE"
  Variable = "Field"
  Procedure = "ModelPDE" "AdvDiffSolver"
  Linear System Solver = Direct
  Steady State Convergence Tolerance = Real 1.0e-9
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Field
  Operator 1 = min
  Variable 2 = Field
  Operator 2 = max
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Check Keywords "Warn"
Header
  Mesh DB "." "square"
End
Simulation
  Max Output Level = 5
  Coordinate System = Cartesian
  Simulation Type = Transient
  Timestepping Method = BDF
  BDF Order = 1
  Timestep Intervals = 2
  Timestep Sizes = 0.1
  Steady State Max Iterations = 1
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
  Initial Condition = 1
End
Body 2
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 3
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 4
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Initial Condition 1
  Field = Real 0.0
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  diffusion coefficient = Real 2.0
  reaction coefficient = Real 0.5
  time derivative coefficient = Real 1.0
End
Body Force 1
  field source = Real 2.0
End
Boundary Condition 1
  Target Boundaries(2) = 1 2
  Field = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  Robin Coefficient = Real 5.0
  External Field = Real 1.0
End
Boundary Condition 3
  Target Boundaries(1) = 4
  Field Flux = Real 3.0
End
Solver 1
  Equation = "ModelPDE"
  Variable = "Field"
  Procedure = "ModelPDE" "AdvDiffSolver"
  Linear System Solver = Direct
  Nonlinear System Max Iterations = 2
  Steady State Convergence Tolerance = Real 0.0
  Integ Test = Logical True
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Field
  Operator 1 = min
  Variable 2 = Field
  Operator 2 = max
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
rm -f "$WORK/shadow-port/ModelPDE.so"
cp "$PORT_SO" "$WORK/shadow-port/ModelPDE.so"
nm -D --defined-only "$WORK/shadow-port/ModelPDE.so" | grep -q "T advdiffsolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|STARTED AT|FINISHED AT" "$1" || true
}
normnames() { # drop the run-timestamp line for the byte-compare
  grep -v "^File started at" "$1" || true
}

echo "-- differential A (steady: convection + source + Robin/flux) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
normnames "$WORK/runA-orig/scalars.dat.names" > "$WORK/n0.txt"
normnames "$WORK/runA-port/scalars.dat.names" > "$WORK/n1.txt"
diff "$WORK/n0.txt" "$WORK/n1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A scalars.dat.names diverges"; diff "$WORK/n0.txt" "$WORK/n1.txt" | head; exit 1; }
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars/names)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
normnames "$WORK/runA-port/scalars.dat.names" > "$WORK/n1.txt"
normnames "$WORK/runA-port2/scalars.dat.names" > "$WORK/n2.txt"
diff "$WORK/n1.txt" "$WORK/n2.txt" > /dev/null 2>&1 || \
  { echo "FAIL: rerun scalars.dat.names diverges"; exit 1; }

echo "-- differential B (transient + Integ Test + 2 nonlinear iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
normnames "$WORK/runB-orig/scalars.dat.names" > "$WORK/m1.txt"
normnames "$WORK/runB-port/scalars.dat.names" > "$WORK/m2.txt"
diff "$WORK/m1.txt" "$WORK/m2.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config B scalars.dat.names diverges"; diff "$WORK/m1.txt" "$WORK/m2.txt" | head; exit 1; }
echo "   config B byte-identical (stdout/stderr/scalars/names)"

echo "-- branch liveness (B differs from A on both arms)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
if diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1; then
  echo "FAIL: config B output identical to config A (dead branch control)"
  exit 1
fi
echo "   configs diverge as expected"

echo "MODEL-PDE PROOF GREEN in $((SECONDS - T0))s"
