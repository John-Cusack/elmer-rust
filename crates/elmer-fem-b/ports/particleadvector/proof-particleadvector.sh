#!/bin/bash
# T021-B ParticleAdvector solver-entry proof: the Rust cdylib replaces
# ParticleAdvector.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `particleadvector_`,
#     `particleadvector_init_`);
#  2. two distinct .sif differentials (default-advect config A; difference-
#     operator + reinitialize + Runge-Kutta + Max-Integration-Time config B,
#     both with 2 steady iterations, i.e. including the repeated-invocation
#     run) are byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch-liveness control: config B takes different code branches and
#     produces different output, and both arms agree on it.
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
WORK="$LANE/scratch/proof-particleadvector"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B ParticleAdvector proof: work=$WORK build=$ELMER_BUILD"
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
  *"1.100.0-nightly"*"feaadeeac"*) ;;
  *) echo "FAIL: nightly feaadeeac required, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libparticleadvector.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- clippy (bar 6: no warnings)"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release -- -D warnings 2>&1 | tail -3)

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ParticleAdvector.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ParticleAdvector module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T particleadvector_" "$WORK/port.syms" \
  || { echo "FAIL: particleadvector_ not exported"; exit 1; }
grep -qx "T particleadvector_init_" "$WORK/port.syms" \
  || { echo "FAIL: particleadvector_init_ not exported"; exit 1; }
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
  Steady State Max Iterations = 2
  Output Intervals = 1
  Max Output Level = 5
End
Body 1
  Target Bodies(1) = 1
  Equation = 1
  Material = 1
End
Body 2
  Target Bodies(1) = 2
  Equation = 1
  Material = 1
End
Body 3
  Target Bodies(1) = 3
  Equation = 1
  Material = 1
End
Body 4
  Target Bodies(1) = 4
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
  Viscosity = 0.1
End
Boundary Condition 1
  Target Boundaries(3) = 1 2 4
  Velocity 1 = 0.0
  Velocity 2 = 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  Velocity 1 = Real 1.0
  Velocity 2 = 0.0
End
Solver 1
  Equation = Navier-Stokes
  Procedure = "FlowSolve" "FlowSolver"
  Variable = Flow Solution[Velocity:2 Pressure:1]
  Exec Solver = Always
  Stabilize = True
  Bubbles = False
  Nonlinear System Convergence Tolerance = 1.0e-8
  Nonlinear System Max Iterations = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU1
End
Solver 2
  Equation = ParticleAdvector
  Procedure = "ParticleAdvector" "ParticleAdvector"
  Simulation Timestep Sizes = Logical False
  Timestep Size = Real 0.02
  Max Timestep Intervals = Integer 5
  Norm Variable Index = Integer 1
  Particle Info = Logical True
  Variable 1 = String "Pressure"
  Variable 2 = String "Particle Status"
  Variable 3 = String "Particle Distance"
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = adv Pressure
  Operator 1 = max
  Variable 2 = adv Pressure
  Operator 2 = min
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
  Max Output Level = 5
End
Body 1
  Target Bodies(1) = 1
  Equation = 1
  Material = 1
End
Body 2
  Target Bodies(1) = 2
  Equation = 1
  Material = 1
End
Body 3
  Target Bodies(1) = 3
  Equation = 1
  Material = 1
End
Body 4
  Target Bodies(1) = 4
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
  Viscosity = 0.1
End
Boundary Condition 1
  Target Boundaries(3) = 1 2 4
  Velocity 1 = 0.0
  Velocity 2 = 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  Velocity 1 = Real 1.0
  Velocity 2 = 0.0
End
Solver 1
  Equation = Navier-Stokes
  Procedure = "FlowSolve" "FlowSolver"
  Variable = Flow Solution[Velocity:2 Pressure:1]
  Exec Solver = Always
  Stabilize = True
  Bubbles = False
  Nonlinear System Convergence Tolerance = 1.0e-8
  Nonlinear System Max Iterations = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU1
End
Solver 2
  Equation = ParticleAdvector
  Procedure = "ParticleAdvector" "ParticleAdvector"
  Simulation Timestep Sizes = Logical False
  Timestep Size = Real 0.02
  Max Timestep Intervals = Integer 5
  Reinitialize Particles = Logical True
  Runge Kutta = Logical True
  Max Integration Time = Real 0.015
  Operator 1 = difference
  Norm Variable Index = Integer 1
  Particle Info = Logical False
  Variable 1 = String "Pressure"
  Variable 2 = String "Particle time"
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = difference Pressure
  Operator 1 = max
  Variable 2 = difference Pressure
  Operator 2 = min
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
rm -f "$WORK/shadow-port/ParticleAdvector.so"
cp "$PORT_SO" "$WORK/shadow-port/ParticleAdvector.so"
nm -D --defined-only "$WORK/shadow-port/ParticleAdvector.so" | grep -q "T particleadvector_" \
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
norm() { grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Assembly \(s\)|Solve \(s\)|pointer init time|Time spent for restart|Factorization time" "$1" || true; }

echo "-- differential A (default-advect branches) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f differs"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f differs (SAVE re-entry)"; diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" | head; exit 1; }
done
norm "$WORK/runA-port/stdout.txt" > "$WORK/r0.txt"
norm "$WORK/runA-port2/stdout.txt" > "$WORK/r1.txt"
diff "$WORK/r0.txt" "$WORK/r1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: rerun stdout.txt differs (SAVE re-entry)"; diff "$WORK/r0.txt" "$WORK/r1.txt" | head; exit 1; }
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (difference + reinitialize + RK + retire branches)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  if [ "$f" = stdout.txt ]; then
    norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
    norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
    diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
      { echo "FAIL: config B stdout diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
  else
    diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" > /dev/null 2>&1 || \
      { echo "FAIL: config B $f differs"; diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" | head; exit 1; }
  fi
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A vs B outputs must differ)"
diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: A and B scalars identical, branches not live"; exit 1; }
diff -q "$WORK/a1.txt" "$WORK/b1.txt" > /dev/null 2>&1 && \
  { echo "FAIL: A and B stdout identical, branches not live"; exit 1; }
echo "   A and B differ as expected (both arms agree per-config)"

echo "PARTICLEADVECTOR PROOF GREEN ($((SECONDS - T0))s)"
