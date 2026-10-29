#!/bin/bash
# T022-A circuitsanddynamics solver-entry proof: the Rust cdylib replaces
# CircuitsAndDynamics.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: the six circuit entries);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines/values:
#       A: harmonic stranded coil + series resistor (harmonic driver,
#          harmonic resistor branch, harmonic reporter, both _init entries);
#       B: transient stranded coil + series resistor over 5 BDF-1 steps
#          (transient driver incl. timestep refresh, transient resistor
#          branch, transient reporter);
#  3. a second full process run of each config is byte-identical
#     (idempotent re-entry through the container/SAVE state);
#  4. branch-liveness control: Resistance 2.5 -> 5.0 changes the coupled
#     solution and both arms agree on the changed output.
#
# Coverage note: the massive/foil-winding element kernels take no part in
# these configs (no massive/foil components); they are verbatim Fortran in
# the accessor module, compiled by the same gfortran with the same flags
# as the reference module, so bit-identity holds by construction wherever
# the exercised call sequences match (proven identical below). The
# MPI_ALLREDUCE arm of the reporter is compiled but never executes in
# these serial runs, exactly as in the reference (guarded by
# Solver % Parallel in both).
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the 1381 test mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-circuitsanddynamics"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/mesh1381
echo "== T022-A circuitsanddynamics proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran 13.3.0 required"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: pinned nightly required"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libcircuitsanddynamics.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/CircuitsAndDynamics.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the CircuitsAndDynamics module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in circuitsanddynamics_ circuitsanddynamics_init_ \
    circuitsanddynamicsharmonic_ circuitsanddynamicsharmonic_init_ \
    circuitsoutput_ circuitsoutput_init_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- mesh input (1381 motor mesh from the Elmer checkout)"
cp -r "$ELMER_FORK/fem/tests/circuits2D_harmonic_resistor/1381" "$WORK/mesh1381"

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
! Proof case A (clean-room): harmonic circuit differential for the
! circuitsanddynamics port. A voltage source drives a stranded coil in
! series with a plain resistor, so the run covers the harmonic driver,
! the harmonic resistor branch, and the harmonic reporter.
Check Keywords "Warn"
INCLUDE mesh1381/mesh.names
$ Circuits = 1
$ C.1.perm = zeros(6)
$ C.1.variables = 6
$ C.1.A = zeros(6,6)
$ C.1.B = zeros(6,6)
$ C.1.Mre = zeros(6,6)
$ C.1.Mim = zeros(6,6)
$ C.1.name.1 = "i_testsource"
$ C.1.name.2 = "v_testsource"
$ C.1.name.3 = "i_component(1)"
$ C.1.name.4 = "v_component(1)"
$ C.1.name.5 = "i_component(2)"
$ C.1.name.6 = "v_component(2)"
$ C.1.B(0,1) = 1
$ C.1.source.1 = "testsource"
$ C.1.B(1,1) = -1
$ C.1.B(1,3) = 1
$ C.1.B(1,5) = 1
$ C.1.B(2,0) = 1
$ C.1.B(2,2) = -1
$ C.1.B(4,0) = 1
$ C.1.B(4,4) = -1
Header 1
   Mesh DB "." "mesh1381"
End
Simulation 1
   Max Output Level = Integer 4
   Coordinate System = Axi symmetric
   Simulation Type = String Steady
   Steady State Max Iterations = Integer 1
   Angular Frequency = Real 314.159265359
End
Solver 1
   Exec Solver = Always
   Equation = Circuits
   Procedure = "CircuitsAndDynamics" "CircuitsAndDynamicsHarmonic"
End
Solver 2
   Equation = "Mag"
   Variable = A[A re:1 A im:1]
   Procedure = "MagnetoDynamics2D" "MagnetoDynamics2DHarmonic"
   Linear System Symmetric = Logical True
   NonLinear System Relaxation Factor = Real 1
   Export Lagrange Multiplier = Logical True
   Linear System Solver = Direct
   Linear System Direct Method = Umfpack
   Linear System Complex = Logical True
   Linear System Convergence Tolerance = Real 1.0e-7
   Linear System Max Iterations = Integer 3000
   Linear System Residual Output = Integer 10
   Linear System Abort Not Converged = Logical False
   Linear System ILUT Tolerance = Real 1.0e-8
   BicgStabL Polynomial Degree = Integer 4
   Steady State Convergence Tolerance = Real 1.0e-6
End
Solver 3
   Exec Solver = Always
   Equation = Circuits Output
   Procedure = "CircuitsAndDynamics" "CircuitsOutput"
End
Solver 4
   Exec Solver = After All
   Equation = "sv"
   Procedure = "SaveData" "SaveScalars"
   Filename = File "scalars.dat"
   Show Norm Name = String "res: v_component(2) re"
End
Equation 1
   Active Solvers(2) = Integer 1 2
End
Equation 2
   Active Solvers(2) = Integer 1 2
End
Material 1
   Name = String iron
   Electric Conductivity = Real 0
   Relative Permeability = Real 2000
End
Material 2
   Name = String air
   Electric Conductivity = Real 0
   Relative Permeability = Real 1
End
Material 3
   Name = String aluminium
   Relative Permeability = Real 1
   Electric Conductivity = Real 1e7
End
Body 1
   Name = String Core
   Target Bodies(6) = Integer 1 2 3 4 5 6
   Equation = Integer 1
   Material = Integer 1
End
Body 2
   Name = String air
   Target Bodies(2) = Integer 9 7
   Equation = Integer 1
   Material = Integer 2
End
Body 3
   Name = String coil
   Target Bodies = Integer 8
   Equation = Integer 2
   Material = Integer 3
End
Component 1
   Name = String coil
   Master Bodies = Integer 3
   Coil Type = String stranded
   Number of Turns = Real 1000
End
Component 2
   Name = String resistor
   Component Type = String Resistor
   Resistance = Real 2.5
End
Body Force 1
   Name = String Circuit
   testsource Re = Real 1.41421356237
   testsource Im = Real 0.0
End
Boundary Condition 1
   Target Boundaries(4) = Integer $ center0 center1 center2 center3
   A re = Real 0
   A im = Real 0
End
Boundary Condition 2
   Target Boundaries(2) = Integer 4 6
   Name = "Air Gap"
   Air Gap = Logical True
   Air Gap Length = Real 0.0001
   Air Gap Relative Permeability = Real 1
End
Boundary Condition 3
   Target Boundaries(7) = Integer 10 13 14 17 18 19 27
   Infinity BC = Logical True
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
! Proof case B (clean-room): transient circuit differential for the
! circuitsanddynamics port. A time-varying voltage source drives a
! stranded coil in series with a plain resistor, so the run covers the
! transient driver, the transient resistor branch, the timestep refresh
! block, and the transient reporter over several timesteps.
Check Keywords "Warn"
INCLUDE mesh1381/mesh.names
$ Circuits = 1
$ C.1.perm = zeros(6)
$ C.1.variables = 6
$ C.1.A = zeros(6,6)
$ C.1.B = zeros(6,6)
$ C.1.Mre = zeros(6,6)
$ C.1.Mim = zeros(6,6)
$ C.1.name.1 = "i_testsource"
$ C.1.name.2 = "v_testsource"
$ C.1.name.3 = "i_component(1)"
$ C.1.name.4 = "v_component(1)"
$ C.1.name.5 = "i_component(2)"
$ C.1.name.6 = "v_component(2)"
$ C.1.B(0,1) = 1
$ C.1.source.1 = "testsource"
$ C.1.B(1,1) = -1
$ C.1.B(1,3) = 1
$ C.1.B(1,5) = 1
$ C.1.B(2,0) = 1
$ C.1.B(2,2) = -1
$ C.1.B(4,0) = 1
$ C.1.B(4,4) = -1
Header 1
   Mesh DB "." "mesh1381"
End
Simulation 1
   Max Output Level = Integer 4
   Coordinate System = Axi symmetric
   Simulation Type = String Transient
   Timestepping Method = String BDF
   BDF Order = Integer 1
   Timestep Sizes(1) = Real 0.0002
   Timestep Intervals(1) = Integer 5
   Steady State Max Iterations = Integer 1
End
Solver 1
   Exec Solver = Always
   Equation = Circuits
   Procedure = "CircuitsAndDynamics" "CircuitsAndDynamics"
End
Solver 2
   Equation = "Mag"
   Variable = A
   Procedure = "MagnetoDynamics2D" "MagnetoDynamics2D"
   NonLinear System Relaxation Factor = Real 1
   Export Lagrange Multiplier = Logical True
   Linear System Solver = Iterative
   Linear System Iterative Method = BiCGStab
   Linear System Convergence Tolerance = Real 1.0e-9
   Linear System Max Iterations = Integer 3000
   Linear System Residual Output = Integer 10
   Linear System Abort Not Converged = Logical False
   Linear System ILUT Tolerance = Real 1.0e-8
   Linear System Preconditioning = Circuit
   BicgStabL Polynomial Degree = Integer 2
   Steady State Convergence Tolerance = Real 1.0e-6
End
Solver 3
   Exec Solver = Always
   Equation = Circuits Output
   Procedure = "CircuitsAndDynamics" "CircuitsOutput"
End
Solver 4
   Exec Solver = After All
   Equation = "sv"
   Procedure = "SaveData" "SaveScalars"
   Filename = File "scalars.dat"
   Show Norm Name = String "res: v_component(2)"
End
Equation 1
   Active Solvers(2) = Integer 1 2
End
Equation 2
   Active Solvers(2) = Integer 1 2
End
Material 1
   Name = String iron
   Electric Conductivity = Real 0
   Relative Permeability = Real 2000
End
Material 2
   Name = String air
   Electric Conductivity = Real 0
   Relative Permeability = Real 1
End
Material 3
   Name = String aluminium
   Relative Permeability = Real 1
   Electric Conductivity = Real 1e7
End
Body 1
   Name = String Core
   Target Bodies(6) = Integer 1 2 3 4 5 6
   Equation = Integer 1
   Material = Integer 1
End
Body 2
   Name = String air
   Target Bodies(2) = Integer 9 7
   Equation = Integer 1
   Material = Integer 2
End
Body 3
   Name = String coil
   Target Bodies = Integer 8
   Equation = Integer 2
   Material = Integer 3
End
Component 1
   Name = String coil
   Master Bodies = Integer 3
   Coil Type = String stranded
   Number of Turns = Real 1000
End
Component 2
   Name = String resistor
   Component Type = String Resistor
   Resistance = Real 2.5
End
Body Force 1
   Name = String Circuit
   testsource = Variable time
     Real MATC "1.41421356237 * sin(2 * 3.14159265358979 * 50 * tx)"
End
Boundary Condition 1
   Target Boundaries(4) = Integer $ center0 center1 center2 center3
   A = Real 0
End
Boundary Condition 2
   Target Boundaries(2) = Integer 4 6
   Name = "Air Gap"
   Air Gap = Logical True
   Air Gap Length = Real 0.0001
   Air Gap Relative Permeability = Real 1
End
Boundary Condition 3
   Target Boundaries(7) = Integer 10 13 14 17 18 19 27
   Infinity BC = Logical True
End
EOF
sed 's/Resistance = Real 2.5/Resistance = Real 5.0/' \
  "$WORK/caseA.sif" > "$WORK/caseA-live.sif"
grep -q "Resistance = Real 5.0" "$WORK/caseA-live.sif" \
  || { echo "FAIL: liveness variant not derived"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/CircuitsAndDynamics.so"
cp "$PORT_SO" "$WORK/shadow-port/CircuitsAndDynamics.so"
nm -D --defined-only "$WORK/shadow-port/CircuitsAndDynamics.so" | grep -q "T circuitsanddynamics_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/mesh1381" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines and mask wall-clock values for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT|File started at" "$1" \
    | sed -E 's/: [0-9.]+ *\(s\)/: (s)/g' || true
}

echo "-- differential A (harmonic stranded + resistor) + rerun (re-entry)"
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

echo "-- differential B (transient stranded + resistor, 5 BDF-1 steps) + rerun"
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseB port runB-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
  norm "$WORK/runB-port2/$f" > "$WORK/b3.txt"
  diff "$WORK/b2.txt" "$WORK/b3.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B rerun $f diverges"; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars) + rerun"

echo "-- liveness control (resistor 2.5 -> 5.0 takes a different branch output)"
run_case caseA-live orig runL-orig
run_case caseA-live port runL-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runL-orig/$f" > "$WORK/l1.txt"
  norm "$WORK/runL-port/$f" > "$WORK/l2.txt"
  diff "$WORK/l1.txt" "$WORK/l2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: liveness $f diverges"; diff "$WORK/l1.txt" "$WORK/l2.txt" | head; exit 1; }
done
echo "   liveness variant agrees on both arms"
if diff "$WORK/a0.txt" "$WORK/l1.txt" > /dev/null 2>&1; then
  echo "FAIL: liveness variant identical to base (branch not taken)"
  exit 1
fi
echo "   liveness variant differs from base (branch live)"

echo "CIRCUITSANDDYNAMICS PROOF GREEN (wall $((SECONDS - T0))s)"
