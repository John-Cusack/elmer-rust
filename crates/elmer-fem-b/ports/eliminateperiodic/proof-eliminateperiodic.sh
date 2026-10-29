#!/bin/bash
# T022-A EliminatePeriodic hook proof: the Rust cdylib replaces
# EliminatePeriodic.so as the `Before Linsolve` hook of a Poisson solve.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `eliminateperiodic_`), with the reference's
#     `NEEDED EliminateDirichlet.so` delegation linkage intact;
#  2. three distinct clean-room .sif differentials are byte-identical to
#     the reference except wall-clock lines: config A (periodic pair +
#     Dirichlet pins: the CAC^T reduction path), config B (no periodic
#     BCs + Direct solve + 2 steady iterations: the EliminateDirichlet
#     delegation path, re-entered), config A2 (config A with 2 steady
#     iterations: in-process re-entry through the reused projection,
#     exercising the wrapper-side SAVE of the first-call column map);
#  3. second full process runs of A and A2 against the port are
#     byte-identical (free/rebuild + SAVE re-entry);
#  4. a liveness variant (changed Dirichlet value) produces different
#     output, and both arms agree on it.
#
# Coverage notes: the serial proof takes the reduction build (projector
# walk, PeriodicCnstr, C build/sort/transpose), the CAC^T multiply, the
# reduced-system Perm remap (transcribed in Rust), the delegation call,
# and the SolveLinearSystem path. The complex, eigen, and PEs > 1 paths
# are transcribed verbatim in the Fortran wrappers but never taken by
# these configs (real Poisson system, serial build, NOFEigenValues == 0).
#
# SAVE note: the original's first-call column map is a non-SAVE local
# read on every later call; the wrapper keeps it with identical contents
# (author's evident intent). The 2-iteration differentials prove the
# re-entry path byte-identical against the reference.
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
WORK="$LANE/scratch/proof-eliminateperiodic"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== EliminatePeriodic proof: work=$WORK build=$ELMER_BUILD"
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
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: pinned nightly required, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
test -f "$ELMER_BUILD/fem/src/modules/EliminateDirichlet.so" || { echo "FAIL: no reference EliminateDirichlet.so (delegation target)"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" ELMER_MODULESDIR="$ELMER_BUILD/fem/src/modules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libeliminateperiodic.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/EliminatePeriodic.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the EliminatePeriodic module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T eliminateperiodic_" "$WORK/port.syms" \
  || { echo "FAIL: eliminateperiodic_ not the exported entry"; exit 1; }
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
  Target Boundaries(1) = 2
  Periodic BC = Integer 2
  Periodic BC Potential = Logical True
End
Boundary Condition 2
  Target Boundaries(1) = 4
End
Boundary Condition 3
  Target Boundaries(1) = 1
  Potential = Real 0.0
End
Boundary Condition 4
  Target Boundaries(1) = 3
  Potential = Real 1.0
End
Solver 1
  Equation = Poisson
  Procedure = "Poisson" "PoissonSolver"
  Variable = Potential
  Variable DOFs = 1
  Before Linsolve = "EliminatePeriodic" "EliminatePeriodic"
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
End
EOF
# Config B: no periodic BCs (all Dirichlet) + Direct solve + 2 steady
# iterations — the EliminateDirichlet delegation path, re-entered.
python3 - "$WORK/caseA.sif" <<'PYEOF'
import sys
t = open(sys.argv[1]).read()
t = t.replace("""Boundary Condition 1
  Target Boundaries(1) = 2
  Periodic BC = Integer 2
  Periodic BC Potential = Logical True
End
Boundary Condition 2
  Target Boundaries(1) = 4
End""", """Boundary Condition 1
  Target Boundaries(1) = 2
  Potential = Real 0.5
End
Boundary Condition 2
  Target Boundaries(1) = 4
  Potential = Real 0.5
End""", 1)
t = t.replace("Steady State Max Iterations = 1",
              "Steady State Max Iterations = 2", 1)
t = t.replace("Linear System Solver = Iterative",
              "Linear System Solver = Direct", 1)
t = t.replace("  Linear System Iterative Method = BiCGStab\n", "")
t = t.replace("  Linear System Preconditioning = ILU0\n", "")
open(sys.argv[1].replace("caseA.sif", "caseB.sif"), "w").write(t)
PYEOF
grep -q "Potential = Real 0.5" "$WORK/caseB.sif"
grep -q "Linear System Solver = Direct" "$WORK/caseB.sif"
# Config A2: config A with 2 steady iterations — in-process re-entry
# through the reused projection.
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    "$WORK/caseA.sif" > "$WORK/caseA2.sif"
# Liveness variant: changed Dirichlet value takes the same reduction path
# to different output.
sed -e 's/Potential = Real 1.0/Potential = Real 2.0/' \
    "$WORK/caseA.sif" > "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/EliminatePeriodic.so"
cp "$PORT_SO" "$WORK/shadow-port/EliminatePeriodic.so"
nm -D --defined-only "$WORK/shadow-port/EliminatePeriodic.so" | grep -q "T eliminateperiodic_" \
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
norm() { # drop wall-clock lines/values for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" || true
}

echo "-- differential A (periodic reduction) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
test -s "$WORK/runA-orig/scalars.dat" || { echo "FAIL: config A produced no scalars"; exit 1; }
grep -q "Eliminated" "$WORK/runA-orig/stdout.txt" || { echo "FAIL: config A missed the reduction path"; exit 1; }
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
  diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" > /dev/null 2>&1 || \
    { norm "$WORK/runA-port/$f" > "$WORK/r1.txt"; norm "$WORK/runA-port2/$f" > "$WORK/r2.txt"; \
      diff "$WORK/r1.txt" "$WORK/r2.txt" > /dev/null 2>&1 || \
      { echo "FAIL: rerun $f diverges"; exit 1; }; }
done
echo "   rerun byte-identical (free/rebuild re-entry)"

echo "-- differential B (delegation + Direct, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
test -s "$WORK/runB-orig/scalars.dat" || { echo "FAIL: config B produced no scalars"; exit 1; }
if grep -q "Eliminated" "$WORK/runB-orig/stdout.txt"; then
  echo "FAIL: config B unexpectedly took the reduction path"; exit 1
fi
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
  diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- differential A2 (reduction, 2 iterations: projection reuse)"
run_case caseA2 orig runA2-orig
run_case caseA2 port runA2-port
run_case caseA2 port runA2-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA2-orig/$f" > "$WORK/c0.txt"
  norm "$WORK/runA2-port/$f" > "$WORK/c1.txt"
  diff "$WORK/c0.txt" "$WORK/c1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A2 $f diverges"; diff "$WORK/c0.txt" "$WORK/c1.txt" | head; exit 1; }
done
echo "   config A2 byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA2-port/$f" > "$WORK/d1.txt"
  norm "$WORK/runA2-port2/$f" > "$WORK/d2.txt"
  diff "$WORK/d1.txt" "$WORK/d2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: A2 rerun $f diverges"; exit 1; }
done
echo "   A2 rerun byte-identical (SAVE re-entry)"

echo "-- liveness control (changed Dirichlet value -> different output)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/e0.txt"
  norm "$WORK/runC-port/$f" > "$WORK/e1.txt"
  diff "$WORK/e0.txt" "$WORK/e1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/e0.txt" "$WORK/e1.txt" | head; exit 1; }
done
if diff "$WORK/runA-port/scalars.dat" "$WORK/runC-port/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: liveness variant produced identical output"; exit 1
fi
echo "   liveness agrees on both arms and differs from A"

echo "ELIMINATEPERIODIC PROOF GREEN in $((SECONDS - T0))s"
