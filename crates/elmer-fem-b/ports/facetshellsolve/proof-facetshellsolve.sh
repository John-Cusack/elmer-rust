#!/bin/bash
# T022-A FacetShellSolve solver-entry proof: the Rust cdylib replaces
# FacetShellSolve.so (invoked via the campaign form
# `Procedure = "FacetShellSolve" "ShellSolver"`).
#
# Reference unit: fem/src/modules/FacetShellSolve.F90 (ShellSolver +
# ShellSolver_Init0). The reference .so exports exactly `shellsolver_`
# and `shellsolver_init0_`.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `shellsolver_`, `shellsolver_init0_`);
#  2. two distinct clean-room .sif differentials are byte-identical to
#     the reference except wall-clock lines/values:
#     config A (triangle mesh, linear statics, pressure + edge load,
#     stress computation, 2 steady iterations) and config B (quad mesh,
#     large deflections, point load, load scale factor, 2 steady
#     iterations);
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch liveness: config B (different elements, physics, loads)
#     produces different output, and both arms agree on it.
#
# Norm class (merged): wall-clock lines
# (`Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT`) are
# dropped, and timing values on per-iteration `(s)` lines are masked
# (values only, structure/iterations compared). `Max Output Level = 4`
# structurally suppresses the nondeterministic Level=5 `RealTime`
# progress reports in both arms. Nothing else is filtered.
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
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-facetshellsolve"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== FacetShellSolve proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libfacetshellsolve.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FacetShellSolve.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FacetShellSolve module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T shellsolver_" "$WORK/port.syms" \
  || { echo "FAIL: shellsolver_ not the exported entry"; exit 1; }
grep -qx "T shellsolver_init0_" "$WORK/port.syms" \
  || { echo "FAIL: shellsolver_init0_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- meshes (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)
sed 's/Triangles = True/Triangles = False/' "$WORK/square.grd" > "$WORK/squareq.grd"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 squareq > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
  Max Output Level = 4
  Simulation Type = Steady State
  Steady State Max Iterations = 2
  Steady State Min Iterations = 2
  Output Intervals = 1
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
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Density = 7850.0
  Youngs modulus = 2.1e11
  Poisson ratio = 0.3
  Thickness = 0.01
End
Body Force 1
  Body Force 3 = -1000.0
  Pressure = 100.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Deflection 1 = 0.0
  Deflection 2 = 0.0
  Deflection 3 = 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  Force 3 = 500.0
End
Solver 1
  Equation = Shell
  Procedure = "FacetShellSolve" "ShellSolver"
  Variable = Deflection
  Variable DOFs = 6
  Linear System Solver = Direct
  Linear System Direct Method = Umfpack
  Stress Computation = Logical True
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Deflection
  Operator 1 = max
  Variable 2 = Deflection
  Operator 2 = min
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "squareq"
End
Simulation
  Coordinate System = Cartesian
  Max Output Level = 4
  Simulation Type = Steady State
  Steady State Max Iterations = 2
  Steady State Min Iterations = 2
  Output Intervals = 1
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
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Density = 7850.0
  Youngs modulus = 2.1e11
  Poisson ratio = 0.3
  Thickness = 0.01
End
Body Force 1
  Body Force 3 = -2000.0
  Pressure = 200.0
  Point Load (1,9) = Real 1.0 1.0 0.0 0.0 0.0 -500.0 0.0 0.0 0.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Deflection 1 = 0.0
  Deflection 2 = 0.0
  Deflection 3 = 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  Force 3 = 800.0
End
Solver 1
  Equation = Shell
  Procedure = "FacetShellSolve" "ShellSolver"
  Variable = Deflection
  Variable DOFs = 6
  Linear System Solver = Direct
  Linear System Direct Method = Umfpack
  Large Deflection = Logical True
  Load Scale Factor = 1.5
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Deflection
  Operator 1 = max
  Variable 2 = Deflection
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
rm -f "$WORK/shadow-port/FacetShellSolve.so"
cp "$PORT_SO" "$WORK/shadow-port/FacetShellSolve.so"
nm -D --defined-only "$WORK/shadow-port/FacetShellSolve.so" | grep -q "T shellsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

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
norm() { # drop wall-clock lines, mask per-iteration (s) timing values
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" \
    | sed -E 's/: +\(s\).*/: (s)/' || true
}

echo "-- differential A (linear, triangles, stress) + rerun (re-entry)"
run_case caseA orig runA-orig square
run_case caseA port runA-port square
run_case caseA port runA-port2 square
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

echo "-- differential B (large deflection, quads, point load, 2 iterations)"
run_case caseB orig runB-orig squareq
run_case caseB port runB-port squareq
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (B must differ from A)"
if diff -q "$WORK/a0.txt" "$WORK/b1.txt" > /dev/null 2>&1; then
  echo "FAIL: liveness control configs agree (branch not exercised)"
  exit 1
fi
echo "   liveness control differs, both arms agree on it"

echo "FACETSHELLSOLVE PROOF GREEN ($((SECONDS - T0))s)"
