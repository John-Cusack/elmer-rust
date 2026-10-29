#!/bin/bash
# T022-A DataToFieldSolver solver-entry proof: the Rust cdylib replaces
# DataToFieldSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `datatofieldsolver_` + `datatofieldsolver_init_`);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines — config A (mask + `Mask Diffusion` +
#     `Continue BC With Boundary Assembly` + Dirichlet, global
#     `Diffusion Coefficient`, 2 steady iterations i.e. the
#     repeated-invocation run through the BulkAssembly SAVE state),
#     config B (point-data file + `Normalize Data by Weight` alias +
#     `Set Constant Weight Sum` + material `Data Diffusivity` local path,
#     1 steady iteration);
#  3. a second full process run of config B against the port is
#     byte-identical (port-vs-port rerun).
# Branch-liveness control: the saved FitA/FitB min/max scalars differ
# across configs A/B in both arms, `Mask is active` prints only in A in
# both arms, so the differential is not vacuous. The deterministic
# `Level=4` lines are compared; only wall-clock lines are filtered.
#
# Documented divergence from the reference: AsciiPointsToMesh declares
# `SAVE :: AllocationsDone, ...` but allocates unconditionally, so a
# second in-process file read aborts in the reference (undefined-behavior
# path; the author's evident intent is allocate-once). The port allocates
# only when not yet associated; the defined first-call path is
# statement-identical, and the rerun is a fresh process in both arms.
#
# Needs a configured+built serial Elmer tree (ELMER_BUILD — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout for the square.grd mesh input and the DataToField2 vals.dat
# point set. Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-datatoffieldsolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square* "$WORK"/vals.dat
echo "== T022-A DataToFieldSolver proof: work=$WORK build=$ELMER_BUILD"
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
  *) echo "FAIL: need nightly 1.100.0"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libdatatoffieldsolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/DataToFieldSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the DataToFieldSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T datatofieldsolver_" "$WORK/port.syms" \
  || { echo "FAIL: datatofieldsolver_ not the exported entry"; exit 1; }
grep -qx "T datatofieldsolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: datatofieldsolver_init_ not exported"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)
cp "$ELMER_FORK/fem/tests/DataToField2/vals.dat" "$WORK/vals.dat"

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
  Max Output Level = Integer 4
End
Body 1
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 2
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 3
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 4
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Initial Condition 1
  Data = Variable Coordinate
    Real MATC "1+tx(0)+2*tx(1)+3*sin(15*tx(0))*cos(10*tx(1))"
  Mask = Variable Coordinate
    Real MATC "sin(15*tx(0))*cos(10*tx(1))"
End
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Name = String "Ideal"
End
Boundary Condition 1
  Target Boundaries(3) = 1 3 4
  FitA Continue = Logical True
End
Boundary Condition 2
  Target Boundaries(1) = 2
  FitA = Variable Coordinate
    Real MATC "1+tx(0)+2*tx(1)"
End
Solver 1
  Equation = String "FitDataA"
  Procedure = "DataToFieldSolver" "DataToFieldSolver"
  Variable = FitA
  Target Variable = String "Data"
  Mask Variable = String "Mask"
  Min Mask Value = Real -0.01
  Max Mask Value = Real 0.01
  Mask Diffusion = Logical True
  Diffusion Coefficient = Real 1.0e-5
  Continue BC With Boundary Assembly = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
  Linear System Residual Output = Integer 50
  Exported Variable 1 = Data
  Exported Variable 2 = Mask
End
Solver 2
  Equation = String "SaveScalarsA"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = FitA
  Operator 1 = min
  Variable 2 = FitA
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
  Steady State Max Iterations = 1
  Output Intervals = 1
  Max Output Level = Integer 4
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
Equation 1
  Active Solvers(2) = 1 2
End
Material 1
  Name = String "Ideal"
  Data Diffusivity(3) = Real 0.1 0.1 0.0
End
Solver 1
  Equation = String "FitDataB"
  Procedure = "DataToFieldSolver" "DataToFieldSolver"
  Variable = FitB
  Normalize Data by Weight = Logical True
  Set Constant Weight Sum = Logical True
  Point Data Filename = File vals.dat
  Point Data Column = Integer 4
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
  Linear System Residual Output = Integer 50
End
Solver 2
  Equation = String "SaveScalarsB"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = FitB
  Operator 1 = min
  Variable 2 = FitB
  Operator 2 = max
End
EOF
grep -c "Mask Variable" "$WORK/caseA.sif"
grep -c "Point Data Filename" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/DataToFieldSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/DataToFieldSolver.so"
nm -D --defined-only "$WORK/shadow-port/DataToFieldSolver.so" | grep -q " T datatofieldsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  cp "$WORK/vals.dat" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" || true
}

echo "-- differential A (mask + boundary assembly, 2 steady iterations)"
run_case caseA orig runA-orig
run_case caseA port runA-port
echo "-- differential B (point-data file + given-weight normalize)"
run_case caseB orig runB-orig
run_case caseB port runB-port
echo "-- rerun of B against the port (fresh process)"
run_case caseB port runB-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
  norm "$WORK/runB-port2/$f" > "$WORK/b3.txt"
  diff "$WORK/b2.txt" "$WORK/b3.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   configs A+B byte-identical (stdout/stderr/scalars) + rerun byte-identical"

echo "-- solved-output sanity (both arms really solved)"
for run in runA-orig runA-port runB-orig runB-port; do
  grep -q "All done" "$WORK/$run/stdout.txt" \
    || { echo "FAIL: $run did not finish the solver"; exit 1; }
  grep -qE "[0-9]+\.[0-9]+.*[0-9]+\.[0-9]+" "$WORK/$run/scalars.dat" \
    || { echo "FAIL: $run scalars.dat has no numeric min/max"; cat "$WORK/$run/scalars.dat"; exit 1; }
done
echo "   all runs solved and saved numeric scalars"

echo "-- branch-liveness control (A and B differ, both arms agree)"
for arm in orig port; do
  sa=$(cat "$WORK/runA-$arm/scalars.dat")
  sb=$(cat "$WORK/runB-$arm/scalars.dat")
  [ "$sa" != "$sb" ] || { echo "FAIL: A($arm) and B($arm) scalars identical"; exit 1; }
  grep -q "Mask is active" "$WORK/runA-$arm/stdout.txt" \
    || { echo "FAIL: mask branch marker missing in A($arm)"; exit 1; }
  grep -q "Mask is active" "$WORK/runB-$arm/stdout.txt" \
    && { echo "FAIL: mask branch marker leaked into B($arm)"; exit 1; } || true
done
echo "   control holds: branches change the output in both arms"

echo "DATATOFFIELDSOLVER PROOF GREEN wall=$((SECONDS - T0))s"
