#!/bin/bash
# T021-B RigidBodyReduction Before-Linsolve-hook proof: the Rust cdylib
# replaces RigidBodyReduction.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `rigidbody_`);
#  2. two distinct .sif differentials (full-feature config A; reduced-
#     branch config B) are byte-identical to the reference except
#     wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. a 2-steady-iteration config C agrees orig-vs-port (in-process
#     SAVE re-entry with MeshChanged false);
#  5. branch-liveness: config A and config B outputs differ.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout for the anturi.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="/home/john/firstmate/data/rs-elmer-t021b/scratch/proof-rigidbodyreduction"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK/shadow-"* "$WORK/case"*.sif "$WORK"/anturi*
echo "== RigidBodyReduction proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/librigidbodyreduction.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- clippy (no warnings)"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release -- -D warnings) \
  || { echo "FAIL: clippy warnings"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/RigidBodyReduction.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the RigidBodyReduction module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T rigidbody_" "$WORK/port.syms" \
  || { echo "FAIL: rigidbody_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's rgdblock anturi.grd)"
cp "$ELMER_FORK/fem/tests/rgdblock/anturi.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 anturi > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "anturi"
End

Simulation
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
End

Body 1
  Equation = 1
  Material = 2
End

Body 2
  Equation = 1
  Material = 1
End

Body 3
  Equation = 1
  Material = 1
  Body Force = 1
  Rigid Body = True
End

Equation 1
  Active Solvers(2) = 1 2
End

Solver 1
  Equation = Stress Analysis
  Variable = Displacement
  Variable DOFs = 3
  Before Linsolve = "RigidBodyReduction" "RigidBody"
  Extend Elastic Region = Logical True
  Extend Elastic Layers = Integer 2
  Optimize Matrix Structure = Logical True
  Reverse Ordering = Logical True
  Output Node Types = Logical True
  Additional Info = Logical True

  Steady State Convergence Tolerance = 1.0e-05
  Linear System Solver = Iterative
  Linear System Iterative Method = CG
  Linear System Max Iterations = 1000
  Linear System Residual Output = 10
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU4
  Linear System ILUT Tolerance = 1.0e-03
  Nonlinear System Convergence Tolerance = 1.0e-05
  Nonlinear System Max Iterations = 1
  Nonlinear System Newton After Iterations = 8
  Nonlinear System Newton After Tolerance = 1.0e-05
  Nonlinear System Relaxation Factor = 1.0
End

Solver 2
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Displacement
  Operator 1 = dofs
End

Material 1
  Name = "Silicon"
  Youngs Modulus = 162.0e9
  Poisson ratio = 0.27
  Density = 2330.0
End

Material 2
  Name = "Weaker Silicon"
  Youngs Modulus = 100.0e9
  Poisson ratio = 0.27
  Density = 2330.0
End

Body Force 1
  Stress Bodyforce 2 = -1e6
End

Boundary Condition 1
  Target Boundaries(2) = 1 2
  Displacement 1 = 0.0
  Displacement 2 = 0.0
  Displacement 3 = 0.0
End
EOF
# Config B: reduced-branch variant — no elastic extension, no bandwidth
# optimize/reverse, no node-type output, no verbose info.
sed -e '/Extend Elastic Region/d; /Extend Elastic Layers/d' \
    -e 's/Optimize Matrix Structure = Logical True/Optimize Matrix Structure = Logical False/' \
    -e '/Reverse Ordering/d; /Output Node Types/d; /Additional Info/d' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -c "RigidBodyReduction" "$WORK/caseB.sif"
# Config C: config A with 2 steady iterations (in-process SAVE re-entry).
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    "$WORK/caseA.sif" > "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/RigidBodyReduction.so"
cp "$PORT_SO" "$WORK/shadow-port/RigidBodyReduction.so"
nm -D --defined-only "$WORK/shadow-port/RigidBodyReduction.so" | grep -q "T rigidbody_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/anturi" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Assembly \(s\)|Solve \(s\)|pointer init time|Time spent for restart|Factorization time|Time used in SolveLinearSystem|Total time spent in RigidBodyReduction" "$1" || true; }

echo "-- differential A (full-feature) + rerun (re-entry)"
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

echo "-- differential B (reduced branches)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- differential C (2 steady iterations, in-process re-entry)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A vs B must differ)"
if diff -q "$WORK/a0.txt" "$WORK/b1.txt" > /dev/null 2>&1; then
  echo "FAIL: configs A and B agree — no branch liveness"; exit 1
fi
echo "   configs A and B differ (branch-liveness proven)"

echo "RIGIDBODYREDUCTION PROOF GREEN ($((SECONDS - T0))s)"
