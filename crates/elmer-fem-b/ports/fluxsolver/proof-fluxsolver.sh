#!/bin/bash
# T022-A FluxSolver solver-entry proof: the Rust cdylib replaces FluxSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `fluxsolver_` + `fluxsolver_init_`);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines — config A (HeatSolver + scalar conductivity,
#     `Calculate Flux` only: the tensor-Rank-1 branch) and config B (full
#     2x2-tensor conductivity, `Calculate Flux` + `Calculate Grad` +
#     both magnitudes + `Enforce Positive Magnitude` +
#     `Constant Bulk Matrix` + `Calculate Flux Abs`/`Calculate Grad Abs`
#     + `Discontinuous Galerkin` reads, 2 steady iterations so the second
#     visit takes the BulkValues-reuse path); both carry a `SaveScalars`
#     solver so the solved flux/gradient fields' min/max land in
#     scalars.dat;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry; the port holds only wrapper-side SAVE state
#     with identical lifetimes to the original's, refilled before use).
# Branch-liveness control: the `Result Norm` report and the saved scalars
# differ across configs A/B in both arms (the keyword and conductivity
# branches change the numerics identically in both arms), so the
# differential is not vacuous. The nondeterministic `Level=5` `RealTime`
# reports are structurally suppressed in both arms via the standard
# `Max Output Level = 4` keyword (the deterministic `Level=4`
# `Result Norm` line still prints and is compared).
#
# Needs a configured+built serial Elmer tree (ELMER_BUILD — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-fluxsolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A FluxSolver proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libfluxsolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FluxSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FluxSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T fluxsolver_" "$WORK/port.syms" \
  || { echo "FAIL: fluxsolver_ not the exported entry"; exit 1; }
grep -qx "T fluxsolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: fluxsolver_init_ not exported"; exit 1; }
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
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
  Heat Conductivity = 1.5
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
  Equation = "Heat Equation"
  Procedure = "HeatSolve" "HeatSolver"
  Variable = Temperature
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
  Linear System Residual Output = 10
End
Solver 2
  Equation = ComputeFlux
  Procedure = "FluxSolver" "FluxSolver"
  Calculate Flux = Logical True
  Target Variable = String Temperature
  Flux Coefficient = String "Heat Conductivity"
  Linear System Solver = Iterative
  Linear System Iterative Method = cg
  Linear System Preconditioning = None
  Linear System Residual Output = 10
  Linear System Max Iterations = Integer 500
  Linear System Convergence Tolerance = 1.0e-10
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature Flux 1
  Operator 1 = min
  Variable 2 = Temperature Flux 1
  Operator 2 = max
  Variable 3 = Temperature Flux 2
  Operator 3 = min
  Variable 4 = Temperature Flux 2
  Operator 4 = max
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
  Material = 2
End
Body 4
  Equation = 1
  Material = 2
End
Equation 1
  Active Solvers(3) = 1 2 3
End
Material 1
  Density = 1.0
  Heat Conductivity(2,2) = Real 1.0 0.0 0.0 10.0
End
Material 2
  Density = 1.0
  Heat Conductivity(2,2) = Real 1.0 0.0 0.0 10.0
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
  Equation = "Heat Equation"
  Procedure = "HeatSolve" "HeatSolver"
  Variable = Temperature
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
  Linear System Residual Output = 10
End
Solver 2
  Equation = ComputeFlux
  Procedure = "FluxSolver" "FluxSolver"
  Calculate Flux = Logical True
  Calculate Grad = Logical True
  Calculate Flux Magnitude = Logical True
  Calculate Grad Magnitude = Logical True
  Enforce Positive Magnitude = Logical True
  Constant Bulk Matrix = Logical True
  Calculate Flux Abs = Logical True
  Calculate Grad Abs = Logical True
  Discontinuous Galerkin = Logical True
  Average Within Materials = Logical True
  Target Variable = String Temperature
  Flux Coefficient = String "Heat Conductivity"
  Linear System Solver = Iterative
  Linear System Iterative Method = cg
  Linear System Preconditioning = None
  Linear System Residual Output = 10
  Linear System Max Iterations = Integer 500
  Linear System Convergence Tolerance = 1.0e-10
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature Grad 1
  Operator 1 = min
  Variable 2 = Temperature Grad 1
  Operator 2 = max
  Variable 3 = Temperature Grad_mag
  Operator 3 = min
  Variable 4 = Temperature Grad_mag
  Operator 4 = max
  Variable 5 = Temperature Flux_abs
  Operator 5 = max
  Variable 6 = Temperature Grad_abs
  Operator 6 = max
End
EOF
grep -c "Heat Conductivity(2,2)" "$WORK/caseB.sif"
grep -c "Enforce Positive Magnitude" "$WORK/caseB.sif"
grep -c "Constant Bulk Matrix" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/FluxSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/FluxSolver.so"
nm -D --defined-only "$WORK/shadow-port/FluxSolver.so" | grep -q " T fluxsolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" | \
    sed -E "s/(Assembly: \(s\)|Solve: \(s\))[0-9.+-Ee ]+/\\1<n>/g" || true
}

echo "-- differentials A/B + rerun of A"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseB port runB-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-orig/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a1.txt" "$WORK/a2.txt" | head -40; exit 1; }
  norm "$WORK/runA-port2/$f" > "$WORK/a3.txt"
  diff "$WORK/a2.txt" "$WORK/a3.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; diff "$WORK/a2.txt" "$WORK/a3.txt" | head -40; exit 1; }
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head -40; exit 1; }
  norm "$WORK/runB-port2/$f" > "$WORK/b3.txt"
  diff "$WORK/b2.txt" "$WORK/b3.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun B $f diverges"; exit 1; }
done
echo "   configs A+B byte-identical (stdout/stderr/scalars) + reruns byte-identical"

echo "-- branch-liveness control (Result Norm differs across A/B, both arms)"
for arm in orig port; do
  grep -h "Result Norm" "$WORK/runA-$arm/stdout.txt" | head -2 > "$WORK/normA-$arm.txt"
  grep -h "Result Norm" "$WORK/runB-$arm/stdout.txt" | head -4 > "$WORK/normB-$arm.txt"
  echo "   $arm A: $(cat "$WORK/normA-$arm.txt" | tr '\n' ';')"
  echo "   $arm B: $(cat "$WORK/normB-$arm.txt" | tr '\n' ';')"
  cmp -s "$WORK/normA-$arm.txt" "$WORK/normB-$arm.txt" && \
    { echo "FAIL: $arm norms identical across configs (vacuous)"; exit 1; }
done
diff "$WORK/normA-orig.txt" "$WORK/normA-port.txt" \
  || { echo "FAIL: Result Norm differs orig vs port on A"; exit 1; }
diff "$WORK/normB-orig.txt" "$WORK/normB-port.txt" \
  || { echo "FAIL: Result Norm differs orig vs port on B"; exit 1; }
echo "   control holds: branches change the numerics in both arms"

echo "FLUXSOLVER PROOF GREEN wall=$((SECONDS - T0))s"
