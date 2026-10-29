#!/bin/bash
# T021-B PoissonDG solver-entry proof: the Rust cdylib replaces PoissonDG.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `poissonsolver_`);
#  2. two distinct .sif differentials (config A: explicit gamma1/beta;
#     config B: gamma1/beta omitted so the entry takes the .NOT.GotIt
#     default branches, with 2 steady iterations, i.e. the
#     repeated-invocation run) are byte-identical to the reference
#     except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state).
#
# Branch-liveness note: LocalJumps and LocalJumpsDiscontBC are textually
# identical in the original (verified by diff — names only differ), and a
# probe run showed interior DG edges carry no BC id (`discont bc` lands in
# Elmer's "unused keywords" report on both arms), so the 'Discont BC'
# call path is unreachable via .sif on this mesh. Config B therefore
# covers the reachable alternate branches instead: the gamma1/beta
# .NOT.GotIt defaults (Gamma = 1.0d3, Beta = Gamma), which provably move
# the output, and both arms agree on them.
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="/home/john/firstmate/data/rs-elmer-t021b/scratch/proof-poissondg"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square
echo "== T021-B PoissonDG proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran pin mismatch: $GFORTRAN_VERSION"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: rustc pin mismatch: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libpoissondg.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/PoissonDG.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the PoissonDG module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T poissonsolver_" "$WORK/port.syms" \
  || { echo "FAIL: poissonsolver_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's PoissonDB square.grd)"
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
  Body Force = 1
End
Equation 1
  Poisson DG = Logical True
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
  Equation = Poisson DG
  Procedure = File "PoissonDG" "PoissonSolver"
  Variable = Potential
  Variable DOFs = 1
  Discontinuous Galerkin = Logical True
  gamma1 = Real 10.0
  beta = Real 10.0
  gamma2 = Real 0.001
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = SaveScalars
  Procedure = File "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
  Variable 3 = Potential
  Operator 3 = sum
End
EOF
# Config B: omit gamma1/beta so the entry takes the .NOT.GotIt default
# branches (Gamma = 1.0d3, Beta = Gamma), plus 2 steady iterations (SAVE
# re-entry across invocations through the jump wrappers' SAVE Nodes).
grep -v -E "^  (gamma1|beta) = Real" "$WORK/caseA.sif" > "$WORK/caseB0.sif"
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    "$WORK/caseB0.sif" > "$WORK/caseB.sif"
rm -f "$WORK/caseB0.sif"
! grep -q "gamma1" "$WORK/caseB.sif" || { echo "FAIL: caseB still sets gamma1"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/PoissonDG.so"
cp "$PORT_SO" "$WORK/shadow-port/PoissonDG.so"
nm -D --defined-only "$WORK/shadow-port/PoissonDG.so" | grep -q "T poissonsolver_" \
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

echo "-- differential A (LocalJumps branch) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 \
    || { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" > /dev/null 2>&1 \
    || { echo "FAIL: rerun $f diverges"; exit 1; }
done
norm "$WORK/runA-port/stdout.txt" > "$WORK/a2.txt"
norm "$WORK/runA-port2/stdout.txt" > "$WORK/a3.txt"
diff "$WORK/a2.txt" "$WORK/a3.txt" > /dev/null 2>&1 \
  || { echo "FAIL: rerun stdout diverges"; diff "$WORK/a2.txt" "$WORK/a3.txt" | head; exit 1; }
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (default branches, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  if [ "$f" = stdout.txt ]; then
    norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
    norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
    diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 \
      || { echo "FAIL: config B stdout diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
  else
    diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" > /dev/null 2>&1 \
      || { echo "FAIL: config B $f diverges"; diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" | head; exit 1; }
  fi
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch liveness (B output must differ from A output)"
norm "$WORK/runB-port/stdout.txt" > "$WORK/bport.txt"
if diff -q "$WORK/a1.txt" "$WORK/bport.txt" > /dev/null 2>&1; then
  echo "FAIL: branch-liveness control produced identical output"
  exit 1
fi
if diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: branch-liveness scalars identical"
  exit 1
fi
echo "   branch liveness confirmed (A and B outputs differ, arms agree)"

echo "POISSONDG PROOF GREEN in ${SECONDS}s (T0=$T0)"
