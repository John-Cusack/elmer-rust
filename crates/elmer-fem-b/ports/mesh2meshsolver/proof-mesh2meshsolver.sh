#!/bin/bash
# Mesh2MeshSolver port proof: the Rust cdylib replaces Mesh2MeshSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `mesh2meshsolver_init0_`, `mesh2meshsolver_`);
#  2. three distinct .sif differentials (explicit-index + hand-picked
#     config A; auto-search + hand-picked incl. missing-var config B;
#     explicit-index + rename + additive config C) are byte-identical to
#     the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (re-entry through the per-invocation wrapper state);
#  4. branch-liveness: configs B and C take different code branches from
#     A and produce different output, and both arms agree on it.
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
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-mesh2meshsolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== Mesh2MeshSolver proof: work=$WORK build=$ELMER_BUILD"
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
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libmesh2meshsolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/Mesh2MeshSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the Mesh2MeshSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T mesh2meshsolver_" "$WORK/port.syms" \
  || { echo "FAIL: mesh2meshsolver_ not the exported entry"; exit 1; }
grep -qx "T mesh2meshsolver_init0_" "$WORK/port.syms" \
  || { echo "FAIL: mesh2meshsolver_init0_ not the exported entry"; exit 1; }
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
  Max Output Level = 10
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
  Active Solvers(3) = 1 2 3
End
Material 1
  Heat Conductivity = 1.0
  Density = 1.0
End
Body Force 1
  Heat Source = Real 1.0
End
Body Force 2
  Heat Source = Real -1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Temperature = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Temperature = Real 1.0
End
Solver 1
  Equation = Heat Equation
  Procedure = "HeatSolve" "HeatSolver"
  Variable = Temperature
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = Mesh2Mesh
  Procedure = "Mesh2MeshSolver" "Mesh2MeshSolver"
  Mesh = "square"
  Variable = Temperature
  Variable 1 = Temperature
  Target Mesh Solver Index = Integer 1
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = min
  Variable 2 = Temperature
  Operator 2 = max
End
EOF
# Config B: auto-search mesh branch + a missing source variable
# (exercises the Perm(0) placeholder path).
sed -e '/Target Mesh Solver Index = Integer 1/d' \
    -e 's/  Variable 1 = Temperature/  Variable 1 = Temperature\n  Variable 2 = NoSuchVar/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
# Config C: rename + additive branch (target must pre-exist).
python3 - "$WORK/caseA.sif" "$WORK/caseC.sif" <<'PYEOF'
import sys
t = open(sys.argv[1]).read()
t = t.replace("  Variable 1 = Temperature\n",
              "  Variable 1 = Temperature\n  Target Variable 1 = Temperature\n"
              "  Variable Additive 1 = Logical True\n", 1)
open(sys.argv[2], "w").write(t)
PYEOF
grep -c "Variable Additive 1" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/Mesh2MeshSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/Mesh2MeshSolver.so"
nm -D --defined-only "$WORK/shadow-port/Mesh2MeshSolver.so" | grep -q "T mesh2meshsolver_" \
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

echo "-- differential A (explicit index + hand-picked) + rerun (re-entry)"
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
echo "   rerun byte-identical (per-invocation re-entry)"

echo "-- differential B (auto-search + missing var)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- differential C (rename + additive)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B and C differ from A; arms agree)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
diff "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: B takes a different branch but stdout matches A"; exit 1; }
echo "   B stdout differs from A (auto-search branch live)"
diff "$WORK/runA-port/scalars.dat" "$WORK/runC-port/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: C takes the additive branch but scalars match A"; exit 1; }
echo "   C scalars differ from A (additive branch live)"

echo "MESH2MESHSOLVER PROOF GREEN ($((SECONDS - T0))s)"
