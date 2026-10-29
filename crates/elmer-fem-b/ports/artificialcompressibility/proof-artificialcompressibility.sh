#!/bin/bash
# T022-A ArtificialCompressibility solver-plugin proof: the Rust cdylib
# replaces ArtificialCompressibility.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: compressibilityscale_, compressibilityscale_init_,
#     compressibilitysolver_);
#  2. two distinct .sif differentials (config A: CompressibilitySolver with
#     Reference Pressure + Mesh Update displacement, 2 steady iterations,
#     i.e. the repeated-invocation run; config B: CompressibilityScale with
#     Flow/Displacement Variable Name + _init entry) are byte-identical to
#     the reference except wall-clock lines;
#  3. a second full process run of each config is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch-liveness controls (config C: Displaced Shape branch; config D:
#     Weight By Displacement branch) produce different output, and both
#     arms agree on it.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config -- both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input. Honest halts: any
# divergence is a FAIL.
# NEVER edit this script to make it pass -- fix the port or the wrappers.
#
# NOTE (evident-intent divergence, per the lane's reference-UB rule): the
# original associates the pressure pointers and then immediately assigns
# PressureExists = .FALSE., so with a pressure variable present its element
# loop could only ever read the uninitialized ReferencePressure local. The
# port uses the pressure values when present. No proof config creates a
# pressure variable, so the defined Reference Pressure path is proven
# byte-identical below; the undefined path is port-vs-port stable only.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-artificialcompressibility"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A ArtificialCompressibility proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: need nightly feaadeeac 2026-09-19"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libartificialcompressibility.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/ArtificialCompressibility.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the ArtificialCompressibility module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in compressibilityscale_ compressibilityscale_init_ compressibilitysolver_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not the exported entry"; exit 1; }
done
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
  Steady State Max Iterations = 2
  Output Intervals = 1
  Max Output Level = 4
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
  Density = 1.0
End
Solver 1
  Equation = Mesh Update
  Procedure = "MeshSolve" "MeshSolver"
  Variable = String Mesh Update
  Variable DOFs = Integer 2
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = CompressibilitySolve
  Procedure = File "ArtificialCompressibility" "CompressibilitySolver"
  Variable = String AC
  Variable DOFs = Integer 1
  Displacement Variable Name = String Mesh Update
  Reference Pressure = Real 1.0e-6
  Displaced Shape = Logical False
  Eliminate Negative = Logical True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-12
  Linear System Preconditioning = ILU0
  Linear System Symmetric = Logical True
End
Solver 3
  Equation = String SaveScalars
  Procedure = File "SaveData" "SaveScalars"
  Filename = File scalars.dat
  Variable 1 = String AC
  Operator 1 = String min
  Variable 2 = String AC
  Operator 2 = String max
  Variable 3 = String AC
  Operator 3 = String mean
End
Boundary Condition 1
  Target Boundaries(1) = Integer 1
  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = Integer 2
  Mesh Update 1 = Real 0.1
  Mesh Update 2 = Real 0.0
End
Boundary Condition 3
  Target Boundaries(1) = Integer 3
  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.0
End
Boundary Condition 4
  Target Boundaries(1) = Integer 4
  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.05
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
  Max Output Level = Integer 10
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
  Artificial Compressibility = Real 1.0e-3
End
Solver 1
  Equation = Mesh Update
  Procedure = "MeshSolve" "MeshSolver"
  Variable = String Mesh Update
  Variable DOFs = Integer 2
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = CompressibilityScale
  Procedure = File "ArtificialCompressibility" "CompressibilityScale"
  Variable = String ScaleVar
  Variable DOFs = Integer 1
  Flow Variable Name = String Mesh Update
  Displacement Variable Name = String Mesh Update
  Artificial Compressibility Scale = Logical True
  Artificial Compressibility Critical Volume = Real 0.01
  Weight By Displacement = Logical False
  Nonlinear System Relaxation Factor = Real 1.0
End
Boundary Condition 1
  Target Boundaries(1) = Integer 1
  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = Integer 2
  Mesh Update 1 = Real 0.1
  Mesh Update 2 = Real 0.0
End
Boundary Condition 3
  Target Boundaries(1) = Integer 3
  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.0
End
Boundary Condition 4
  Target Boundaries(1) = Integer 4
  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.05
End
EOF
# Liveness controls: different code branches of each entry.
sed -e 's/Displaced Shape = Logical False/Displaced Shape = Logical True/' \
  "$WORK/caseA.sif" > "$WORK/caseC.sif"
sed -e 's/Weight By Displacement = Logical False/Weight By Displacement = Logical True/' \
  "$WORK/caseB.sif" > "$WORK/caseD.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ArtificialCompressibility.so"
cp "$PORT_SO" "$WORK/shadow-port/ArtificialCompressibility.so"
nm -D --defined-only "$WORK/shadow-port/ArtificialCompressibility.so" | grep -q "T compressibilitysolver_" \
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
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" || true
}

echo "-- differential A (CompressibilitySolver, Reference Pressure branch, 2 iterations) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  cmp "$WORK/runA-orig/$f" "$WORK/runA-port/$f" \
    || { echo "FAIL: config A $f differs"; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/ar1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/ar2.txt"
  diff "$WORK/ar1.txt" "$WORK/ar2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A rerun $f differs (SAVE re-entry)"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (CompressibilityScale + _init) + rerun"
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseB port runB-port2
cmp "$WORK/runB-orig/stderr.txt" "$WORK/runB-port/stderr.txt" \
  || { echo "FAIL: config B stderr differs"; exit 1; }
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b0.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b1.txt"
diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config B stdout diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
echo "   config B byte-identical (stdout/stderr)"
for f in stdout.txt stderr.txt; do
  norm "$WORK/runB-port/$f" > "$WORK/br1.txt"
  norm "$WORK/runB-port2/$f" > "$WORK/br2.txt"
  diff "$WORK/br1.txt" "$WORK/br2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B rerun $f differs (SAVE re-entry)"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- liveness C (Displaced Shape branch): arms agree, output differs from A"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stderr.txt scalars.dat; do
  cmp "$WORK/runC-orig/$f" "$WORK/runC-port/$f" \
    || { echo "FAIL: config C $f differs"; exit 1; }
done
norm "$WORK/runC-orig/stdout.txt" > "$WORK/c0.txt"
norm "$WORK/runC-port/stdout.txt" > "$WORK/c1.txt"
diff "$WORK/c0.txt" "$WORK/c1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config C stdout diverges"; diff "$WORK/c0.txt" "$WORK/c1.txt" | head; exit 1; }
cmp -s "$WORK/c1.txt" "$WORK/a1.txt" && \
  { echo "FAIL: liveness C identical to A (branch not taken?)"; exit 1; }
echo "   config C byte-identical, differs from A"

echo "-- liveness D (Weight By Displacement branch): arms agree, output differs from B"
run_case caseD orig runD-orig
run_case caseD port runD-port
cmp "$WORK/runD-orig/stderr.txt" "$WORK/runD-port/stderr.txt" \
  || { echo "FAIL: config D stderr differs"; exit 1; }
norm "$WORK/runD-orig/stdout.txt" > "$WORK/d0.txt"
norm "$WORK/runD-port/stdout.txt" > "$WORK/d1.txt"
diff "$WORK/d0.txt" "$WORK/d1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config D stdout diverges"; diff "$WORK/d0.txt" "$WORK/d1.txt" | head; exit 1; }
cmp -s "$WORK/d1.txt" "$WORK/b1.txt" && \
  { echo "FAIL: liveness D identical to B (branch not taken?)"; exit 1; }
echo "   config D byte-identical, differs from B"

echo "ARTIFICIALCOMPRESSIBILITY PROOF GREEN ($((SECONDS - T0))s wall)"
