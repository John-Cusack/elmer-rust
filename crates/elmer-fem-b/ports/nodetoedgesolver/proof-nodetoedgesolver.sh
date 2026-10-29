#!/bin/bash
# NodeToEdgeSolver port proof: the Rust cdylib replaces NodeToEdgeSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `extrudedrestart_init0_`, `extrudedrestart_`,
#     `nodetoedgefield_init0_`, `nodetoedgefield_`);
#  2. four distinct clean-room .sif differentials are byte-identical to
#     the reference except wall-clock lines:
#       A: extruded quad mesh, Poisson az -> NodeToEdgeField (scalar,
#          non-Piola) -> SaveScalars;
#       B: A + Use Piola Transform + Constant Bulk Matrix with 2 steady
#          iterations (Piola path + bulk-matrix reuse path);
#       C: 2D Poisson restart (az.dat) + ExtrudedRestart auto-search,
#          single variable with rename (create path, layers=5);
#       D: explicit Target Mesh Solver Index, two variables (overwrite +
#          create paths) with a second Poisson solve on the 2D mesh;
#  3. a second full process run of config A is byte-identical
#     (re-entry through the wrapper-side SAVE state);
#  4. branch-liveness: B differs from A (Piola path live), D differs
#     from C (explicit-index + second-variable paths live), and both
#     arms agree on each.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the quads.grd mesh input (the upstream
# quad mesh of this solver family, mgdyn_steady_quad_extruded_restart).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-nodetoedgesolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/quads*
echo "== NodeToEdgeSolver proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libnodetoedgesolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/NodeToEdgeSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the NodeToEdgeSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in extrudedrestart_ extrudedrestart_init0_ nodetoedgefield_ nodetoedgefield_init0_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the checkout's upstream quad mesh)"
cp "$ELMER_FORK/fem/tests/mgdyn_steady_quad_extruded_restart/quads.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 quads.grd > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/case2d.sif" <<'EOF'
Header
  Mesh DB "." "quads"
End
Simulation
  Max Output Level = 10
  Coordinate System = Cartesian
  Simulation Type = Steady
  Steady State Max Iterations = 1
  Output Intervals(1) = 1
  Output File = "az.dat"
  Output Variable 1 = az
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
Equation 1
  Active Solvers(2) = 1 2
End
Solver 1
  Equation = PoissonAz
  Procedure = "Poisson" "PoissonSolver"
  Variable = az
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars2d.dat
  Variable 1 = az
  Operator 1 = max
End
Material 1
  Density = 1.0
End
Body Force 1
  Source = Real 1.0e6
End
Boundary Condition 1
  Target Boundaries(2) = 1 2
  az = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  az = Real 0.0
End
EOF
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "quads"
End
Simulation
  Max Output Level = 10
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Extruded Mesh Layers = 4
  Extruded Mesh Slices = Logical True
  Single Mesh = True
  Extruded Min Coordinate = 0.0
  Extruded Max Coordinate = 0.01
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
Equation 1
  Active Solvers(3) = 1 2 3
End
Solver 1
  Equation = PoissonAz
  Procedure = "Poisson" "PoissonSolver"
  Variable = az
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = NodeToEdge
  Procedure = "NodeToEdgeSolver" "NodeToEdgeField"
  Variable = azproj
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = GCR
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1.0e-6
  Linear System Preconditioning = none
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = azproj
  Operator 1 = max abs
  Variable 2 = az
  Operator 2 = max abs
End
Material 1
  Density = 1.0
End
Body Force 1
  Source = Real 1.0e6
End
Boundary Condition 1
  Target Boundaries(2) = 1 2
  az = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  az = Real 0.0
End
EOF
# Config B: Piola + constant-bulk-matrix branch, 2 steady iterations.
cp "$WORK/caseA.sif" "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
t = t.replace("Steady State Max Iterations = 1", "Steady State Max Iterations = 2\n  Steady State Min Iterations = 2", 1)
t = t.replace("""  Variable = azproj
  Variable DOFs = 1
""", """  Variable = azproj
  Variable DOFs = 1
  Use Piola Transform = Logical True
  Constant Bulk Matrix = Logical True
""", 1)
open(p, "w").write(t)
PYEOF
grep -q "Use Piola Transform" "$WORK/caseB.sif"
# Config C: restart-fed ExtrudedRestart, auto-search mesh, rename (create).
cat > "$WORK/caseC.sif" <<'EOF'
Check Keywords "Warn"

Header
  Mesh DB "." "quads"
End
Simulation
  Max Output Level = 10
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Extruded Mesh Layers = 4
  Extruded Mesh Slices = Logical True
  Single Mesh = True
  Extruded Min Coordinate = 0.0
  Extruded Max Coordinate = 0.01
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
Equation 1
  Active Solvers(3) = 1 2 3
End
Solver 1
  Equation = PoissonAz
  Procedure = "Poisson" "PoissonSolver"
  Variable = az
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = MapCoordinates
  Procedure = "NodeToEdgeSolver" "ExtrudedRestart"
  Mesh = "quads"
  Restart File = File "az.dat"
  Restart Position = Integer 0
  Variable 1 = az
  Extruded Variable 1 = String az_mapped
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = az_mapped
  Operator 1 = max abs
  Variable 2 = az
  Operator 2 = max abs
End
Material 1
  Density = 1.0
End
Body Force 1
  Source = Real 1.0e6
End
Boundary Condition 1
  Target Boundaries(2) = 1 2
  az = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  az = Real 0.0
End
EOF
# Config D: explicit target-mesh index, two variables (overwrite + create).
cat > "$WORK/caseD.sif" <<'EOF'
Check Keywords "Warn"

Header
  Mesh DB "." "quads"
End
Simulation
  Max Output Level = 10
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Extruded Mesh Layers = 4
  Extruded Mesh Slices = Logical True
  Single Mesh = True
  Extruded Min Coordinate = 0.0
  Extruded Max Coordinate = 0.01
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
Equation 1
  Active Solvers(4) = 1 2 3 4
End
Solver 1
  Equation = PoissonAz
  Procedure = "Poisson" "PoissonSolver"
  Variable = az
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = MapCoordinates
  Procedure = "NodeToEdgeSolver" "ExtrudedRestart"
  Mesh = "quads"
  Restart File = File "az.dat"
  Restart Position = Integer 0
  Target Mesh Solver Index = Integer 3
  Variable 1 = az
  Variable 2 = az
  Extruded Variable 2 = String az_copy
End
Solver 3
  Equation = PoissonAzSecond
  Procedure = "Poisson" "PoissonSolver"
  Mesh = "quads"
  Variable = az
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 4
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = az
  Operator 1 = max abs
  Variable 2 = az_copy
  Operator 2 = max abs
End
Material 1
  Density = 1.0
End
Body Force 1
  Source = Real 1.0e6
End
Boundary Condition 1
  Target Boundaries(2) = 1 2
  az = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  az = Real 0.0
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
rm -f "$WORK/shadow-port/NodeToEdgeSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/NodeToEdgeSolver.so"
nm -D --defined-only "$WORK/shadow-port/NodeToEdgeSolver.so" | grep -q "T nodetoedgefield_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/quads" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
run_extruded() { # <case> <arm> <rundir> — 2D restart step then 3D main step
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp -r "$WORK/quads" "$WORK/$run"/
  cp "$WORK/case2d.sif" "$WORK/$run/case2d.sif"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case2d.sif > stdout2d.txt 2> stderr2d.txt)
  local rc=$?
  echo "$run($arm) 2D exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) 2D step exited $rc"; exit 1; }
  test -f "$WORK/$run/quads/az.dat" || { echo "FAIL: $run($arm) produced no quads/az.dat"; exit 1; }
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|Time spent|Factorization time" "$1" || true
}

echo "-- differential A (scalar non-Piola projection) + rerun (re-entry)"
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

echo "-- differential B (Piola + constant bulk matrix, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- differential C (ExtrudedRestart auto-search + rename create)"
run_extruded caseC orig runC-orig
run_extruded caseC port runC-port
for f in stdout.txt stdout2d.txt stderr.txt stderr2d.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"

echo "-- differential D (explicit index + overwrite + second create)"
run_extruded caseD orig runD-orig
run_extruded caseD port runD-port
for f in stdout.txt stdout2d.txt stderr.txt stderr2d.txt scalars.dat; do
  norm "$WORK/runD-orig/$f" > "$WORK/d1.txt"
  norm "$WORK/runD-port/$f" > "$WORK/d2.txt"
  diff "$WORK/d1.txt" "$WORK/d2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config D $f diverges"; diff "$WORK/d1.txt" "$WORK/d2.txt" | head; exit 1; }
done
echo "   config D byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B vs A, D vs C; arms agree)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
diff "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: B takes a different branch but stdout matches A"; exit 1; }
echo "   B stdout differs from A (Piola + bulk-matrix branch live)"
norm "$WORK/runC-port/stdout.txt" > "$WORK/c2.txt"
norm "$WORK/runD-port/stdout.txt" > "$WORK/d2.txt"
diff "$WORK/c2.txt" "$WORK/d2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: D takes a different branch but stdout matches C"; exit 1; }
echo "   D stdout differs from C (explicit-index + two-variable branch live)"

echo "NODETOEDGESOLVER PROOF GREEN ($((SECONDS - T0))s)"
