#!/bin/bash
# T022-A BeamSolver3D solver-entry proof: the Rust cdylib replaces
# BeamSolver3D.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `timoshenkosolver_` + `timoshenkosolver_init0_`);
#  2. three distinct .sif differentials are byte-identical to the
#     reference except wall-clock lines — config A (steady cantilever
#     baseline, constant pressure load); config B (steady + `Displace
#     Mesh` branch with a different load magnitude); config C
#     (transient, 2 timesteps: the `MassAssembly`/Bossak branch);
#     all three carry a `SaveScalars` solver so the tip deflection
#     values land in scalars.dat;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry; the port holds no cross-call state — the
#     original allocates nothing across calls).
# Branch-liveness control: the saved tip values differ across configs
# A/B/C in both arms (the keyword and simulation-type branches change
# the numerics identically in both arms), so the differential is not
# vacuous. The nondeterministic `Level=5` reports are structurally
# suppressed in both arms via the standard `Max Output Level = 4`
# keyword (deterministic `Level=4` lines still print and are compared).
#
# Needs a configured+built serial Elmer tree (ELMER_BUILD — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout only for the beam1d mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-beamsolver3d"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/beam1d
echo "== T022-A BeamSolver3D proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libbeamsolver3d.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/BeamSolver3D.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the BeamSolver3D module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T timoshenkosolver_" "$WORK/port.syms" \
  || { echo "FAIL: timoshenkosolver_ not the exported entry"; exit 1; }
grep -qx "T timoshenkosolver_init0_" "$WORK/port.syms" \
  || { echo "FAIL: timoshenkosolver_init0_ not exported"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (beam1d from the Elmer checkout's Beam_3D_Cantilever test)"
cp -r "$ELMER_FORK/fem/tests/Beam_3D_Cantilever/beam1d" "$WORK/"

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Check Keywords "Warn"
Header
  Mesh DB "." "beam1d"
End
Simulation
  Max Output Level = 4
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Youngs Modulus = Real 2.0e-1
  Shear Modulus = Real 1.0
  Principal Direction 2(3) = Real 0.0 0.0 -1.0
  Second Moment of Area 2 = Real 1.0
  Second Moment of Area 3 = Real 1.0
  Cross Section Area = Real 1.0
  Torsional Constant = Real 1.0
  Density = 2700.0
End
Body Force 1
  Body Force 1 = Real 0.0
  Body Force 2 = Real 1.0e-2
  Body Force 3 = Real 0.0
End
Solver 1
  Equation = "Timoshenko Beam Equations"
  Procedure = "BeamSolver3D" "TimoshenkoSolver"
  Nonlinear System Max Iterations = 1
  Nonlinear System Convergence Tolerance = 1.0e-7
  Linear System Solver = Iterative
  Linear System Preconditioning = Diagonal
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1.0e-9
  Linear System Iterative Method = IDRS
  Steady State Convergence Tolerance = 1.0e-9
End
Solver 2
  Equation = "Save Scalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = U 1
  Variable 2 = U 2
  Variable 3 = U 3
  Variable 4 = Theta 1
  Variable 5 = Theta 2
  Variable 6 = Theta 3
  Save Points(1) = 2
End
Boundary Condition 1
  Target Nodes(1) = 1
  U 1 = Real 0.0
  U 2 = Real 0.0
  U 3 = Real 0.0
  Theta 1 = Real 0.0
  Theta 2 = Real 0.0
  Theta 3 = Real 0.0
End
Boundary Condition 2
  Target Nodes(1) = 2
End
EOF
# Config B: Displace Mesh branch + different load magnitude.
sed -e 's/Body Force 2 = Real 1.0e-2/Body Force 2 = Real 3.0e-2/' \
    -e 's/Linear System Iterative Method = IDRS/Linear System Iterative Method = IDRS\n  Displace Mesh = Logical True/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -c "Displace Mesh" "$WORK/caseB.sif"
# Config C: transient (MassAssembly/Bossak branch), 2 timesteps.
sed -e 's/Simulation Type = Steady State/Simulation Type = Transient/' \
    -e 's/Steady State Max Iterations = 1/Timestep Sizes(1) = 0.1\n  Timestep Intervals(1) = 2/' \
    -e 's/Body Force 2 = Real 1.0e-2/Body Force 2 = Real 2.0e-2/' \
    "$WORK/caseA.sif" > "$WORK/caseC.sif"
grep -c "Timestep Sizes" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/BeamSolver3D.so"
cp "$PORT_SO" "$WORK/shadow-port/BeamSolver3D.so"
nm -D --defined-only "$WORK/shadow-port/BeamSolver3D.so" | grep -q " T timoshenkosolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/beam1d" "$WORK/$run"/
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

echo "-- differentials A/B/C + rerun of A"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   configs A+B+C byte-identical (stdout/stderr/scalars) + rerun byte-identical"

echo "-- branch-liveness control (tip values differ across A/B/C, both arms)"
for arm in orig port; do
  cmp -s "$WORK/runA-$arm/scalars.dat" "$WORK/runB-$arm/scalars.dat" && \
    { echo "FAIL: A($arm) and B($arm) scalars identical"; exit 1; }
  cmp -s "$WORK/runA-$arm/scalars.dat" "$WORK/runC-$arm/scalars.dat" && \
    { echo "FAIL: A($arm) and C($arm) scalars identical"; exit 1; }
done
echo "   control holds: branches change the numerics in both arms"

echo "BEAMSOLVER3D PROOF GREEN wall=$((SECONDS - T0))s"
