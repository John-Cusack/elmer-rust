#!/bin/bash
# T021-B SpringAssembly solver-entry proof: the Rust cdylib replaces SpringAssembly.so.
#
# Proves, on pinned gfortran 13.3.0 + lane nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `springassembler_`);
#  2. two distinct .sif differentials (config A: explicit Displacement
#     Variable Name + per-component Spring branch, 1 steady iteration;
#     config B: default variable name + scalar-Spring Warn branch + Mass
#     reads + Target-Nodes-only CYCLE branch, 2 steady iterations, i.e.
#     the repeated-invocation run) are byte-identical to the reference
#     except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the call-local state).
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
WORK="/home/john/firstmate/data/rs-elmer-t021b/scratch/proof-springassembly"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== SpringAssembly proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0, got: $GFORTRAN_VERSION"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libspringassembly.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/SpringAssembly.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the SpringAssembly module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T springassembler_" "$WORK/port.syms" \
  || { echo "FAIL: springassembler_ not the exported entry"; exit 1; }
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
  Max Output Level = 5
  Coordinate System = Cartesian 2D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Steady State Min Iterations = 1
  Output Intervals = 1
End
Body 1
  Equation = 1
  Material = 1
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
End
Material 1
  Density = 1.0
  Youngs Modulus = 70e9
  Poisson Ratio = 0.3
End
Solver 1
  Equation = Elasticity Solver
  Variable = Displacement
  Variable DOFs = 2
  Procedure = "ElasticSolve" "ElasticSolver"
  Assembly Solvers(1) = 2
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU2
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-10
  Linear System Residual Output = 20
  Nonlinear System Max Iterations = 1
  Steady State Convergence Tolerance = 1.0e-4
  Calculate Loads = Logical True
End
Solver 2
  Equation = "Assemble Springs"
  Exec Solver = Never
  Procedure = "SpringAssembly" "SpringAssembler"
  Displacement Variable Name = "Displacement"
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Displacement
  Operator 1 = max
End
Equation 1
  Active Solvers(3) = 1 2 3
  Plane Stress = False
End
Boundary Condition 1
  Target Boundaries = 1
  Displacement 1 = 0
  Displacement 2 = 0
End
Boundary Condition 2
  Target Boundaries = 2
  Surface Traction 2 = Real -1.0e6
End
Boundary Condition 3
  Target Nodes(1) = 240
  Spring 2 = Real 1.0e8
End
EOF
python3 - "$WORK/caseA.sif" "$WORK/caseB.sif" <<'PYEOF'
import sys
t = open(sys.argv[1]).read()
# Config B: default variable-name path, scalar-Spring Warn branch,
# Spring-1 + Mass reads, a Target-Nodes-only CYCLE BC, 2 iterations.
t = t.replace('Steady State Max Iterations = 1',
              'Steady State Max Iterations = 2')
t = t.replace('  Displacement Variable Name = "Displacement"\n', '')
t = t.replace("""Boundary Condition 3
  Target Nodes(1) = 240
  Spring 2 = Real 1.0e8
End""", """Boundary Condition 3
  Target Nodes(2) = 239 240
  Spring = Real 1.0e7
  Spring 1 = Real 5.0e7
  Mass 1 = Real 2.5
End
Boundary Condition 4
  Target Nodes(1) = 100
End""")
open(sys.argv[2], 'w').write(t)
PYEOF
grep -c "Target Nodes" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/SpringAssembly.so"
cp "$PORT_SO" "$WORK/shadow-port/SpringAssembly.so"
nm -D --defined-only "$WORK/shadow-port/SpringAssembly.so" | grep -q "T springassembler_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT" "$1" || true
}

echo "-- differential A (explicit variable name + Spring-2 branch) + rerun (re-entry)"
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
echo "   rerun byte-identical (call-local re-entry)"

echo "-- differential B (Warn + Spring-1/Mass + CYCLE branches, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B output must differ from A output)"
norm "$WORK/runA-port/scalars.dat" > "$WORK/c1.txt"
norm "$WORK/runB-port/scalars.dat" > "$WORK/c2.txt"
diff -q "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: branch-liveness control produced identical output"; exit 1; }
echo "   branch-liveness control GREEN (outputs differ, both arms agree)"

echo "SPRINGASSEMBLY PROOF GREEN in $((SECONDS - T0))s"
