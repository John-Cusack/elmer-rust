#!/bin/bash
# T022-A DistributeSource solver-entry proof: the Rust cdylib replaces
# DistributeSource.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `distributesource_`);
#  2. two distinct .sif differentials that exercise the solver (config A:
#     the reference two-mesh donut->sphere projection with precomputed
#     source + fixing term; config B: same pipeline with different
#     'Local Current Density' ValueList data) are byte-identical to the
#     reference except wall-clock lines;
#  3. a second full process run of config A against the port is
#     byte-identical (re-entry; note the entry runs `Exec Solver =
#     "Before all"`, i.e. once per simulation even across steady
#     iterations, so in-process re-entry does not occur in these
#     configs; the wrapper SAVE Nodes state additionally allocates only
#     if unallocated — the author's evident first-mesh-sizing intent —
#     so a repeated invocation is safe by construction);
#  4. branch-liveness control: config B (different ValueList data, hence
#     different projected sums) produces different output from config A,
#     and both arms agree on each config. (Separately observed in dev:
#     dropping the `"g fix"` export takes the `FixSol`-unassociated Warn
#     path identically on both arms — downstream Solver 2 then fails on
#     the missing fixing term on both arms — so the branch itself
#     agrees; it cannot serve as an exit-0 differential and is not part
#     of the GREEN criteria.)
#
# Two-mesh staging: the reference regression test
# fem/tests/DistributeSourceMgdyn uses a local "donut" mesh (carrying
# the 'Local Current Density' body force) projected onto the global
# "sphere" mesh (carrying the edge-basis 'g' variable). The proof copies
# both mesh dirs next to each case.sif, exactly like the reference test
# layout (`Mesh DB "." "sphere"` + `Mesh = "donut"` on Solver 1).
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-distributesource"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/donut "$WORK"/sphere
echo "== T022-A DistributeSource proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libdistributesource.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/DistributeSource.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the DistributeSource module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T distributesource_" "$WORK/port.syms" \
  || { echo "FAIL: distributesource_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- meshes (donut + sphere from the reference DistributeSourceMgdyn test)"
cp -r "$ELMER_FORK/fem/tests/DistributeSourceMgdyn/donut" "$ELMER_FORK/fem/tests/DistributeSourceMgdyn/sphere" "$WORK/"

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
! Config A: reference two-mesh projection (local donut mesh current
! density projected onto the global sphere mesh edge-basis 'g', with
! precomputed source + fixing term), plus a SaveScalars probe so the
! differential covers scalars.dat. Adapted from
! fem/tests/DistributeSourceMgdyn/case.sif: Max Output Level 4 (wall-clock
! discipline), no VTU post file (not compared), SaveScalars added.
Header
  Mesh DB "." "sphere"
End

Simulation
  Max Output Level = 4
  Coordinate System = "Cartesian"
  Simulation Type = "Steady State"
  Output Intervals = 1
  Steady State Max Iterations = 1

  Use Mesh Projector = Logical False
  Try Linear Search If QTree Fails = Logical False
End

Body 1
  Name = "ball"
  Target Bodies(1) = 1
  Equation = 1
  Material = 1
  Body Force = 1
End

Body 2
  Name = "donut"
  Target Bodies(1) = 2
  Equation = 2
  Body Force = 2
  Initial Condition = 1
End

Material 1
  Name = "Ideal"
  Reluctivity = 1
End

Solver 1
  Equation = "DistributeSource"
  Exec Solver = "Before all"
  Mesh = "donut"
  Procedure = "DistributeSource" "DistributeSource"
End

Initial Condition 1
End

Solver 2
  Equation = MGDynamics
  Procedure = "MagnetoDynamics" "WhitneyAVSolver"
  Variable = AV

  Fix Input Current Density = Logical True

  Steady State Convergence Tolerance = 1.0e-5
  Nonlinear System Convergence Tolerance = 1.0e-8
  Nonlinear System Max Iterations = 1
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Newton After Tolerance = 1.0e-3
  Nonlinear System Relaxation Factor = 1

  Linear System Solver = iterative
  Linear System Symmetric = True
  Linear System Iterative Method = BiCGStabL
  Linear System GCR Restart = 150
  Linear System Max Iterations = 2500
  Linear System Convergence Tolerance = 1.0e-7
  Linear System Preconditioning = none
  Linear System Residual Output = 10
  Linear System Abort Not Converged = False

  Exported Variable 1 = g

  Precomputed Source Term = Logical True
  Precomputed Fixing Term = Logical True
End


Solver 3
  Equation = MGDynamicsCalc
  Procedure = "MagnetoDynamics" "MagnetoDynamicsCalcFields"
  Potential Variable = String "AV"

  Exported Variable 1 = "g fix"

  Steady State Convergence Tolerance = 1.0e-5
  Nonlinear System Convergence Tolerance = 1.0e-8
  Nonlinear System Max Iterations = 20
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Newton After Tolerance = 1.0e-3
  Nonlinear System Relaxation Factor = 1

  Calculate Magnetic Vector Potential = Logical True

  Linear System Solver = Iterative
  Linear System Symmetric = True
  Linear System Iterative Method = CG
  Linear System Max Iterations = 5000
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
  Linear System ILUT Tolerance = 1.0e-3
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1
  Linear System Precondition Recompute = 1
End


Solver 4
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = AV
  Operator 1 = min
  Variable 2 = AV
  Operator 2 = max
End


Equation 1 :: Active Solvers(3) = 2 3 4
Equation 2 :: Active Solvers(1) = 1

Body Force 1
  Name = "GlobalCurrentDensity"
End

Body Force 2
  Name = "LocalCurrentDensity"
  Local Current Density 1 = Opposes Coordinate 2
  Local Current Density 2 = Equals  Coordinate 1
End

Boundary Condition 1
  Name = "External"
  Target Boundaries(1) = 3
  AV = Real 0
  AV {e} = Real 0
End

Solver 2 :: Reference Norm = 6.73333546E-01
Solver 3 :: Reference Norm = 1.15188631E-01
EOF
python3 - "$WORK/caseA.sif" "$WORK/caseB.sif" <<'PYEOF'
import sys
# Config B: different 'Local Current Density' ValueList data on the
# local mesh -> different projected source/fixing sums and downstream
# norms (liveness), same pipeline otherwise.
t = open(sys.argv[1]).read()
t = t.replace('  Local Current Density 1 = Opposes Coordinate 2\n  Local Current Density 2 = Equals  Coordinate 1',
              '  Local Current Density 1 = Equals Coordinate 2\n  Local Current Density 2 = Opposes Coordinate 1')
open(sys.argv[2], 'w').write(t)
PYEOF
grep -c "Local Current Density" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/DistributeSource.so"
cp "$PORT_SO" "$WORK/shadow-port/DistributeSource.so"
nm -D --defined-only "$WORK/shadow-port/DistributeSource.so" | grep -q "T distributesource_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/donut" "$WORK/sphere" "$WORK/$run"/
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

echo "-- differential A (two-mesh projection) + rerun (re-entry)"
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

echo "-- differential B (different ValueList data)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness (B differs from A, both arms agree)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
if diff "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1; then
  echo "FAIL: config B output identical to config A (dead control)"
  exit 1
fi
echo "   liveness control alive (A and B outputs differ)"

echo "DISTRIBUTESOURCE PROOF GREEN"
echo "wall seconds: $((SECONDS - T0))"
