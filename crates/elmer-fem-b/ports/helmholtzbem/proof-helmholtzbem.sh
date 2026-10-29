#!/bin/bash
# T022-A HelmholtzBEM solver-entry proof: the Rust cdylib replaces
# HelmholtzBEM.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `helmholtzbemsolver_`, `matvec_`, `precond_`
#     plus the plugin-private `GlobMatC` matrix `B` symbol);
#  2. two distinct .sif differentials (direct-solve config A;
#     iterative-HUTI + all-flux config B with 2 steady iterations, i.e.
#     the repeated-invocation run) are byte-identical to the reference
#     except wall-clock lines/values;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. a branch-liveness control: config B produces different output
#     from config A, and both arms agree on it.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
#
# Timing-value note: besides the standard Elmer wall-clock lines, this
# solver prints three CPU-time report lines (`Assembly (s):`,
# `Solve (s):`, `Post Processing (s):`). Only their numeric values are
# masked via the generalized `(s)` form (values only, structure and
# iteration counts still compared); everything else is byte-compared.
#
# Defined-path note (Wave-1 rule): the original sizes `P1`/`P2` with a
# stale local whose value never leaks into any observable (both are
# fully overwritten per element); the port implements the author's
# evident intent (`Model % MaxElementNodes`). Reference runs are
# self-consistent, so the defined path is proven byte-identical.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-helmholtzbem"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== HelmholtzBEM proof: work=$WORK build=$ELMER_BUILD"
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
(cd "$PORT_DIR" && CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$PORT_DIR/target}" ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libhelmholtzbem.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/HelmholtzBEM.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the HelmholtzBEM module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T helmholtzbemsolver_" "$WORK/port.syms" \
  || { echo "FAIL: helmholtzbemsolver_ not the exported entry"; exit 1; }
grep -qx "T matvec_" "$WORK/port.syms" \
  || { echo "FAIL: matvec_ not exported"; exit 1; }
grep -qx "T precond_" "$WORK/port.syms" \
  || { echo "FAIL: precond_ not exported"; exit 1; }
grep -qx "B __globmatc_MOD_matrix" "$WORK/port.syms" \
  || { echo "FAIL: GlobMatC matrix symbol not exported"; exit 1; }
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
  Angular Frequency = Real 44.0
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
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Body Id = 1
  Pressure 1 = Real 1.0
  Pressure 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Body Id = 1
  Flux 1 = Real 0.5
  Flux 2 = Real 0.0
End
Solver 1
  Equation = HelmholtzBEM
  Procedure = "HelmholtzBEM" "HelmholtzBEMSolver"
  Variable = Pressure
  Variable DOFs = 2
  Exported Variable 1 = Flux
  Exported Variable 1 DOFs = 2
  Linear System Solver = direct
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Pressure
  Operator 1 = min
  Variable 2 = Pressure
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
  Angular Frequency = Real 143.0
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
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Body Id = 1
  Flux 1 = Real 1.0
  Flux 2 = Real 0.25
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Body Id = 1
  Pressure 1 = Real 0.0
  Pressure 2 = Real 0.5
End
Solver 1
  Equation = HelmholtzBEM
  Procedure = "HelmholtzBEM" "HelmholtzBEMSolver"
  Variable = Pressure
  Variable DOFs = 2
  Exported Variable 1 = Flux
  Exported Variable 1 DOFs = 2
  Linear System Solver = iterative
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-10
  Linear System Divergence Limit = Real 1.0e6
  Linear System Residual Output = 1
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Pressure
  Operator 1 = min
  Variable 2 = Pressure
  Operator 2 = max
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
rm -f "$WORK/shadow-port/HelmholtzBEM.so"
cp "$PORT_SO" "$WORK/shadow-port/HelmholtzBEM.so"
nm -D --defined-only "$WORK/shadow-port/HelmholtzBEM.so" | grep -q "T helmholtzbemsolver_" \
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
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; tail -20 "$WORK/$run/stderr.txt"; exit 1; }
}
norm() { # drop wall-clock lines/values for the byte-compare
  # Timing-value mask (values only, structure compared): the generalized
  # `(s)` form covers the solver's `Assembly (s):`, `Solve (s):` and
  # `Post Processing (s):` report lines.
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" | \
    sed -E -e 's/\(s\):.*/\(s\):/' || true
}

echo "-- differential A (direct + mixed Dirichlet/flux) + rerun (re-entry)"
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

echo "-- differential B (iterative HUTI + swapped BC mix, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head -40; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (B output differs from A; both arms agree)"
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b0.txt"
diff -q "$WORK/a0.txt" "$WORK/b0.txt" > /dev/null 2>&1 && \
  { echo "FAIL: liveness control did not diverge"; exit 1; }
diff "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: liveness scalars did not diverge"; exit 1; }
echo "   liveness control diverges (both arms)"

echo "HELMHOLTZBEM PROOF GREEN (wall $((SECONDS - T0))s)"
