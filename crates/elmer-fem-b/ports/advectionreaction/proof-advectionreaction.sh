#!/bin/bash
# T022-A AdvectionReaction solver-entry proof: the Rust cdylib replaces
# AdvectionReaction.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `advectionreactionsolver_`,
#     `advectionreactionsolver_init_`, `advectionreactionsolver_post_`);
#  2. two distinct .sif differentials (config A: constant-convection DG
#     assembly with nodal-average export, i.e. all three entries run;
#     config B: 2 nonlinear iterations + Limit Solution branch + 2 steady
#     iterations, i.e. the in-process SAVE re-entry path) are
#     byte-identical to the reference except wall-clock lines/values;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch liveness: config B output differs from config A, and both
#     arms agree on config B.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL. Reference-vs-port re-entry is
# the defined path (the port reuses SAVE state identically); the bar's
# port-vs-port rerun is also proven.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-advectionreaction"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A AdvectionReaction proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: need rustc nightly 1.100.0"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$PORT_DIR/target}" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libadvectionreaction.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/AdvectionReaction.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the AdvectionReaction module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T advectionreactionsolver_" "$WORK/port.syms" \
  || { echo "FAIL: advectionreactionsolver_ not the exported entry"; exit 1; }
grep -qx "T advectionreactionsolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: advectionreactionsolver_init_ not exported"; exit 1; }
grep -qx "T advectionreactionsolver_post_" "$WORK/port.syms" \
  || { echo "FAIL: advectionreactionsolver_post_ not exported"; exit 1; }
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
  Simulation Type = Steady State
  Steady State Max Iterations = 1
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
  Body Force = 2
End
Equation 1
  Active Solvers(1) = 1
  Convection = constant
End
Material 1
  Convection Velocity 1 = Real 2.0
  Convection Velocity 2 = Real 1.0
  Tracer Gamma = Real 1.0
End
Body Force 1
  Tracer Source = Real 1.0
End
Body Force 2
  Tracer Source = Real -1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Tracer = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Tracer = Real 1.0
End
Solver 1
  Equation = "AdvReact"
  Procedure = "AdvectionReaction" "AdvectionReactionSolver"
  Variable = "Tracer"
  Discontinuous Galerkin = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = Integer 500
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Nonlinear System Max Iterations = Integer 1
  Calculate Nodal Average = Logical True
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Tracer
  Operator 1 = min
  Variable 2 = Tracer
  Operator 2 = max
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 2
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
  Body Force = 2
End
Equation 1
  Active Solvers(1) = 1
  Convection = constant
End
Material 1
  Convection Velocity 1 = Real 2.0
  Convection Velocity 2 = Real 1.0
  Tracer Gamma = Real 1.0
  Tracer Upper Limit = Real 2.0
  Tracer Lower Limit = Real -1.0
End
Body Force 1
  Tracer Source = Real 2.0
End
Body Force 2
  Tracer Source = Real -1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Tracer = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Tracer = Real 1.0
End
Solver 1
  Equation = "AdvReact"
  Procedure = "AdvectionReaction" "AdvectionReactionSolver"
  Variable = "Tracer"
  Discontinuous Galerkin = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = Integer 500
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Nonlinear System Max Iterations = Integer 2
  Limit Solution = Logical True
  Calculate Nodal Average = Logical True
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Tracer
  Operator 1 = min
  Variable 2 = Tracer
  Operator 2 = max
End
EOF
grep -c "Tracer" "$WORK/caseA.sif" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/AdvectionReaction.so"
cp "$PORT_SO" "$WORK/shadow-port/AdvectionReaction.so"
nm -D --defined-only "$WORK/shadow-port/AdvectionReaction.so" | grep -q "T advectionreactionsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/AdvectionReaction.so" | grep -q "T advectionreactionsolver_post_" \
  || { echo "FAIL: staged port .so lost its post entry"; exit 1; }

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
norm() { # drop wall-clock lines; mask per-iteration (s) timing values only
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" |
    sed -E -e '/\(s\)/ s/[0-9]+\.[0-9]+/T/g' || true
}

echo "-- differential A (DG assembly + nodal-average post) + rerun (re-entry)"
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
  if [ "$f" = "stdout.txt" ]; then
    norm "$WORK/runA-port/$f" > "$WORK/r0.txt"
    norm "$WORK/runA-port2/$f" > "$WORK/r1.txt"
    diff "$WORK/r0.txt" "$WORK/r1.txt" > /dev/null 2>&1 || \
      { echo "FAIL: rerun $f diverges (SAVE re-entry)"; diff "$WORK/r0.txt" "$WORK/r1.txt" | head; exit 1; }
  else
    diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" > /dev/null 2>&1 || \
      { echo "FAIL: rerun $f diverges (SAVE re-entry)"; diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" | head; exit 1; }
  fi
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (limit + 2 nonlinear iters + 2 steady iters)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  if [ "$f" = "stdout.txt" ]; then
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

echo "-- branch liveness (B must differ from A; both arms agree on each)"
if diff -q "$WORK/a1.txt" "$WORK/b1.txt" > /dev/null 2>&1; then
  echo "FAIL: config B output identical to config A (no branch delta exercised)"
  exit 1
fi
echo "   config A and B outputs differ (branch delta exercised)"

echo "ADVECTIONREACTION PROOF GREEN (wall $((SECONDS - T0))s)"
