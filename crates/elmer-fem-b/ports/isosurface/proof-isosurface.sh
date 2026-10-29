#!/bin/bash
# T022-A Isosurface solver-entry proof: the Rust cdylib replaces Isosurface.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `isosurfacesolver_` + `isosurfacesolver_init_`);
#  2. three distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines/values — config A (2D single-level isoline
#     with the `Variable 1` interpolant fallback, 1 steady iteration);
#     config B (2D fixed-surface run with 2 steady iterations: iteration 1
#     creates the isoline, iteration 2 takes the `Isosurface Fixed` remap
#     path through the SAVE state, plus `Save Gmsh Geo File` output
#     compared byte-for-byte); config C (3D isosurface with `Save STL File`
#     output compared byte-for-byte);
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch liveness: config B/C outputs differ from config A, and both
#     arms agree on each config.
#
# The serial proof never reaches the `MPI_*` collectives (all are behind
# the same `ParEnv % PEs > 1` guards as the original and resolve against
# the same `libmpi_stubs` as the reference). `Max Output Level = 4`
# structurally suppresses the nondeterministic `Level=5` reports; the
# deterministic `Level=4` `IsosurfaceSolver` lines still print and are
# compared.
#
# Needs a configured+built serial Elmer tree (ELMER_BUILD — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout only for the square/cube mesh inputs.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-isosurface"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square* "$WORK"/cube*
echo "== T022-A Isosurface proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libisosurface.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/Isosurface.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the Isosurface module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T isosurfacesolver_" "$WORK/port.syms" \
  || { echo "FAIL: isosurfacesolver_ not the exported entry"; exit 1; }
grep -qx "T isosurfacesolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: isosurfacesolver_init_ not exported"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- meshes (ElmerGrid from the Elmer checkout's square/cube inputs)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)
sed -e 's/Element Divisions 1 = 20/Element Divisions 1 = 6/' \
    -e 's/Element Divisions 2 = 10/Element Divisions 2 = 4/' \
    -e 's/Element Divisions 3 = 10/Element Divisions 3 = 4/' \
    "$ELMER_FORK/fem/tests/ConformingDisplacementCube/cube.grd" > "$WORK/cube.grd"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 cube > /dev/null)

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
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
End
Body Force 1
  Source = Real 1.0
End
Body Force 2
  Source = Real -1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Potential = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Potential = Real 1.0
End
Solver 1
  Equation = PoissonEq
  Procedure = "Poisson" "PoissonSolver"
  Variable = Potential
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = Isoline
  Procedure = "Isosurface" "IsosurfaceSolver"
  Isosurface Variable = String Potential
  Isosurface Value = Real 0.5
  Variable 1 = String Potential
  Isomesh Output Active = Logical False
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
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
  Body Force = 1
End
Body 3
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 4
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
End
Body Force 1
  Source = Real 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Potential = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Potential = Real 0.0
End
Solver 1
  Equation = PoissonEq
  Procedure = "Poisson" "PoissonSolver"
  Variable = Potential
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
  Steady State Convergence Tolerance = Real 0.0
End
Solver 2
  Equation = Isoline
  Procedure = "Isosurface" "IsosurfaceSolver"
  Isosurface Variable = String Potential
  Isosurface Value = Real 0.02
  Isosurface Interpolant 1 = String Potential
  Isosurface Fixed = Logical True
  Isomesh Output Active = Logical False
  Save Gmsh Geo File = Logical True
  Geo Filename = File "mesh.geo"
  Steady State Convergence Tolerance = Real 0.0
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
End
EOF
cat > "$WORK/caseC.sif" <<'EOF'
Header
  Mesh DB "." "cube"
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
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
End
Body Force 1
  Source = Real 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Potential = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Potential = Real 1.0
End
Solver 1
  Equation = PoissonEq
  Procedure = "Poisson" "PoissonSolver"
  Variable = Potential
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = Isosurf
  Procedure = "Isosurface" "IsosurfaceSolver"
  Isosurface Variable = String Potential
  Isosurface Value = Real 0.5
  Variable 1 = String Potential
  Isomesh Output Active = Logical False
  Save STL File = Logical True
  STL Filename = File "mesh.stl"
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
End
EOF
grep -c "Isosurface" "$WORK/caseA.sif" "$WORK/caseB.sif" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/Isosurface.so"
cp "$PORT_SO" "$WORK/shadow-port/Isosurface.so"
nm -D --defined-only "$WORK/shadow-port/Isosurface.so" | grep -q "T isosurfacesolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/Isosurface.so" | grep -q "T isosurfacesolver_init_" \
  || { echo "FAIL: staged port .so lost its init entry"; exit 1; }

run_case() { # <case> <arm> <rundir> <mesh>
  local case=$1 arm=$2 run=$3 mesh=$4
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/$mesh" "$WORK/$run"/
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

echo "-- differential A (2D single level) + rerun (re-entry)"
run_case caseA orig runA-orig square
run_case caseA port runA-port square
run_case caseA port runA-port2 square
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a-orig-stdout.norm"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a-port-stdout.norm"
diff "$WORK/a-orig-stdout.norm" "$WORK/a-port-stdout.norm" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a-orig-stdout.norm" "$WORK/a-port-stdout.norm" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a-port-$f.norm"
  norm "$WORK/runA-port2/$f" > "$WORK/a-port2-$f.norm"
  diff "$WORK/a-port-$f.norm" "$WORK/a-port2-$f.norm" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (2D fixed-surface remap + gmsh geo, 2 iterations)"
run_case caseB orig runB-orig square
run_case caseB port runB-port square
for f in stdout.txt stderr.txt scalars.dat mesh.geo; do
  norm "$WORK/runB-orig/$f" > "$WORK/b-orig-$f.norm"
  norm "$WORK/runB-port/$f" > "$WORK/b-port-$f.norm"
  diff "$WORK/b-orig-$f.norm" "$WORK/b-port-$f.norm" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b-orig-$f.norm" "$WORK/b-port-$f.norm" | head; exit 1; }
done
test -s "$WORK/runB-orig/mesh.geo" || { echo "FAIL: config B geo file empty (not exercising the writer)"; exit 1; }
grep -q "Remapping the fields" "$WORK/runB-orig/stdout.txt" \
  || { echo "FAIL: config B never took the fixed-surface remap path"; exit 1; }
echo "   config B byte-identical (stdout/stderr/scalars/mesh.geo)"

echo "-- differential C (3D isosurface + STL)"
run_case caseC orig runC-orig cube
run_case caseC port runC-port cube
for f in stdout.txt stderr.txt scalars.dat mesh.stl; do
  norm "$WORK/runC-orig/$f" > "$WORK/c-orig-$f.norm"
  norm "$WORK/runC-port/$f" > "$WORK/c-port-$f.norm"
  diff "$WORK/c-orig-$f.norm" "$WORK/c-port-$f.norm" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c-orig-$f.norm" "$WORK/c-port-$f.norm" | head; exit 1; }
done
test -s "$WORK/runC-orig/mesh.stl" || { echo "FAIL: config C STL file empty (not exercising the writer)"; exit 1; }
echo "   config C byte-identical (stdout/stderr/scalars/mesh.stl)"

echo "-- branch liveness (B and C must differ from A; both arms agree on each)"
if diff -q "$WORK/a-orig-stdout.norm" "$WORK/b-orig-stdout.txt.norm" > /dev/null 2>&1; then
  echo "FAIL: config A and B outputs identical (branch delta missing)"; exit 1
fi
if diff -q "$WORK/a-orig-stdout.norm" "$WORK/c-orig-stdout.txt.norm" > /dev/null 2>&1; then
  echo "FAIL: config A and C outputs identical (branch delta missing)"; exit 1
fi
grep -q "Isosurface mesh has" "$WORK/a-orig-stdout.norm" \
  || { echo "FAIL: no isosurface size report in config A output"; exit 1; }
grep -q "Isosurface mesh has" "$WORK/b-orig-stdout.txt.norm" \
  || { echo "FAIL: no isosurface size report in config B output"; exit 1; }
grep -q "Isosurface mesh has" "$WORK/c-orig-stdout.txt.norm" \
  || { echo "FAIL: no isosurface size report in config C output"; exit 1; }
echo "   config A, B, C outputs differ as expected (branch delta exercised)"

echo "ISOSURFACE PROOF GREEN (wall $((SECONDS - T0))s)"
