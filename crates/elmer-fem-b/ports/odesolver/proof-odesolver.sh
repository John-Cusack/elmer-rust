#!/bin/bash
# T021-B OdeSolver solver-entry proof: the Rust cdylib replaces OdeSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `odesolver_` + `odesolver_init_`);
#  2. two distinct clean-room .sif differentials exercise the solver
#     (config A: 1st-order + Active Components branch, Dofs=1;
#      config B: 2nd-order + solver-section branch, Dofs=2, 2 nonlinear
#      iterations) are byte-identical to the reference except
#     wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry);
#  4. branch-liveness control: configs A and B take different code
#     branches, produce different output, and both arms agree on each.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config —
# both arms share it, so the differential is self-consistent).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="/home/john/firstmate/data/rs-elmer-t021b/scratch/proof-odesolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B OdeSolver proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0, got: $GFORTRAN_VERSION"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly (feaadeeac 2026-09-19), got: $RUSTC_VERSION"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libodesolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- clippy (no warnings)"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release -- -D warnings) 2>&1 | tail -3

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/OdeSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the OdeSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T odesolver_" "$WORK/port.syms" \
  || { echo "FAIL: odesolver_ not an exported entry"; exit 1; }
grep -qx "T odesolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: odesolver_init_ not an exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  CHECK KEYWORDS Warn
  Mesh DB "." "square"
  Include Path ""
  Results Directory ""
End
Simulation
  Max Output Level = 6
  Coordinate System = Cartesian
  Simulation Type = Transient
  Steady State Max Iterations = 1
  Output Intervals = 1
  Timestep Intervals = 20
  Timestep Sizes = 0.01
  Timestepping Method = bdf
  Bdf Order = 4
End
Body 1
  Name = "Body"
  Equation = 1
  Material = 1
End
Equation 1
  Name = "Equations"
End
Solver 1
  Equation = Ode
  Procedure = "OdeSolver" "OdeSolver"
  Variable = -global ode
  Variable Dofs = 1
  Time Derivative Order = 1
  Linear System Solver = Direct
  Linear System Save = Logical False
  Active Components = Integer 1
End
Solver 2
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = time
  Variable 2 = Ode
End
Material 1
  Name = "Material"
End
Component 1
  Name = String "Ode Component"
  Stiffness Matrix 11 = Real 1.0
  Damping Matrix 11 = Real 1.0
  Force 1 = Real 1.0
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  CHECK KEYWORDS Warn
  Mesh DB "." "square"
  Include Path ""
  Results Directory ""
End
Simulation
  Max Output Level = 6
  Coordinate System = Cartesian
  Simulation Type = Transient
  Steady State Max Iterations = 1
  Output Intervals = 1
  Timestep Intervals = 20
  Timestep Sizes = 0.01
  Timestepping Method = bdf
  Bdf Order = 4
End
Body 1
  Name = "Body"
  Equation = 1
  Material = 1
End
Equation 1
  Name = "Equations"
End
Solver 1
  Equation = Ode
  Procedure = "OdeSolver" "OdeSolver"
  Variable = -global ode
  Variable Dofs = 2
  Time Derivative Order = 2
  Nonlinear System Max Iterations = 2
  Linear System Solver = Direct
  Linear System Save = Logical False
  Stiffness Matrix 11 = Real 2.0
  Stiffness Matrix 12 = Real 0.1
  Stiffness Matrix 21 = Real 0.1
  Stiffness Matrix 22 = Real 2.0
  Damping Matrix 11 = Real 0.5
  Damping Matrix 12 = Real 0.05
  Damping Matrix 21 = Real 0.05
  Damping Matrix 22 = Real 0.5
  Mass Matrix 11 = Real 1.0
  Mass Matrix 12 = Real 0.0
  Mass Matrix 21 = Real 0.0
  Mass Matrix 22 = Real 1.0
  Force 1 = Real 1.0
  Force 2 = Real 0.5
End
Solver 2
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = time
  Variable 2 = Ode
End
Material 1
  Name = "Material"
End
EOF
echo "   configs written"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/OdeSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/OdeSolver.so"
nm -D --defined-only "$WORK/shadow-port/OdeSolver.so" | grep -q "T odesolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (1st order + Active Components) + rerun (re-entry)"
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

echo "-- differential B (2nd order + solver section + 2 iters)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A and B must differ; arms agree)"
norm "$WORK/runA-orig/scalars.dat" > "$WORK/la.txt"
norm "$WORK/runB-orig/scalars.dat" > "$WORK/lb.txt"
if diff -q "$WORK/la.txt" "$WORK/lb.txt" > /dev/null 2>&1; then
  echo "FAIL: configs A and B produce identical output (no branch delta)"
  exit 1
fi
echo "   A and B outputs differ on both arms (branch delta live)"

echo "ODESOLVER PROOF GREEN in $((SECONDS - T0))s"
