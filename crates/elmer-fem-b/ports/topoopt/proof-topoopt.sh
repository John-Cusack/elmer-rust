#!/bin/bash
# T021-B TopoOpt solver-entry port proof: the Rust cdylib replaces TopoOpt.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: topoopt_, topoopt_init_, topoopt_init0_);
#  2. two distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines — config A (LinearElasticity primary +
#     distance density filter + B&W extraction, 3 steady iterations, i.e.
#     the repeated-invocation run) and config B (HeatEq primary + PDE
#     density filter, 2 steady iterations);
#  3. a second full process run of config B against the port is
#     byte-identical (idempotent re-entry through the SAVE state);
#  4. branch-liveness control: config C (same as A but
#     Filter Method "sensitivity", i.e. the other pre-filter branch)
#     produces different output, and both arms agree on it.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial — both arms
# share it, so the differential is self-consistent) and the Elmer checkout
# only for the square.grd mesh input. Honest halts: any divergence is a
# FAIL. NEVER edit this script to make it pass — fix the port or the
# wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
WORK="/home/john/firstmate/data/rs-elmer-t021b/scratch/proof-topoopt"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B TopoOpt proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: expected gfortran 13.3.0, got: $GFORTRAN_VERSION"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: expected rustc 1.100.0-nightly (feaadeeac 2026-09-19), got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libtopoopt.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/TopoOpt.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the TopoOpt module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T topoopt_" "$WORK/port.syms" \
  || { echo "FAIL: topoopt_ not an exported entry"; exit 1; }
grep -qx "T topoopt_init_" "$WORK/port.syms" \
  || { echo "FAIL: topoopt_init_ not an exported entry"; exit 1; }
grep -qx "T topoopt_init0_" "$WORK/port.syms" \
  || { echo "FAIL: topoopt_init0_ not an exported entry"; exit 1; }
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
  Steady State Max Iterations = 3
  Output Intervals = 0
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
  Name = "Guess"
  topo rho = Real 0.5
  topo mult = Real 0.125
End
Equation 1
  Name = "Elast"
  Active Solvers(2) = 1 2
  Plane Stress = True
End
Material 1
  Name = "PlaneStuff"
  Youngs modulus = 1.0
  Poisson ratio = 0.3
End
Boundary Condition 1
  Target Boundaries(1) = 4
  Name = "Clamp"
  Displacement 1 = Real 0.0
  Displacement 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 2
  Name = "Load"
  Displacement 2 = Real -1.0
End
Solver 1
  Equation = LinearElasticity
  Procedure = "StressSolve" "StressSolver"
  Variable = -dofs 2 Displacement
  Nonlinear System Max Iterations = 1
  Linear System Solver = Direct
  Linear System Direct Method = umfpack
  Local Matrix Identical = Logical False
  Local Matrix Storage = Logical True
  Matrix Multiplier Name = String "topo mult"
  Optimize Bandwidth = False
  Steady State Convergence Tolerance = 1.0e-4
  Displace Mesh = False
End
Solver 2
  Equation = TopoOpt
  Procedure = "TopoOpt" "TopoOpt"
  Filter Method = String "density"
  Filter Type = String "distance"
  Distance Filter Radius = Real 0.15
  Penalty Exponent = Real 3.0
  Minimum Relative Density = Real 1.0e-7
  Volume Fraction = Real 0.5
  Bisection Search Tolerance = Real 1.0e-3
  Linear System Solver = direct
  Linear System Direct Method = umfpack
  Steady State Convergence Tolerance = Real 1.0e-12
  Create BW Topology = Logical True
End
Solver 3
  Exec Solver = After Simulation
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = topoopt_nrm
  Operator 1 = max
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
  Output Intervals = 0
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
  Initial Condition = 1
End
Body 2
  Equation = 1
  Material = 1
  Body Force = 1
  Initial Condition = 1
End
Body 3
  Equation = 1
  Material = 1
  Body Force = 1
  Initial Condition = 1
End
Body 4
  Equation = 1
  Material = 1
  Body Force = 1
  Initial Condition = 1
End
Initial Condition 1
  Name = "Guess"
  topo rho = Real 0.4
  topo mult = Real 0.064
End
Body Force 1
  Name = "Source"
  Volumetric Heat Source = Real 1.0
End
Equation 1
  Name = "HeatEq"
  Active Solvers(2) = 1 2
End
Material 1
  Name = "Ideal"
  Heat Conductivity = Real 1.0
  Heat Capacity = Real 0.0
  Density = Real 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Name = "Natural"
End
Boundary Condition 2
  Target Boundaries(1) = 2
  Name = "Null"
  Temperature = Real 0.0
End
Solver 1
  Equation = HeatEq
  Procedure = "HeatSolve" "HeatSolver"
  Variable = Temperature
  Nonlinear System Max Iterations = 1
  Linear System Solver = Direct
  Linear System Direct Method = umfpack
  Local Matrix Identical = Logical True
  Local Matrix Storage = Logical True
  Vector Assembly = False
  Matrix Multiplier Name = String "topo mult"
  Optimize Bandwidth = False
  Steady State Convergence Measure = solution
  Steady State Convergence Tolerance = 1.0e-4
End
Solver 2
  Equation = TopoOpt
  Procedure = "TopoOpt" "TopoOpt"
  Filter Method = String "density"
  Filter Type = String "pde"
  Pde Filter Diffusion Constant = Real 1.0
  Variable = xnodal
  Penalty Exponent = Real 3.0
  Minimum Relative Density = Real 1.0e-3
  Volume Fraction = Real 0.4
  Linear System Solver = direct
  Linear System Direct Method = umfpack
  Steady State Convergence Tolerance = Real 1.0e-12
End
Solver 3
  Exec Solver = After Simulation
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = min
  Variable 2 = Temperature
  Operator 2 = max
End
EOF
sed 's/Filter Method = String "density"/Filter Method = String "sensitivity"/' \
    "$WORK/caseA.sif" > "$WORK/caseC.sif"
grep -q 'Filter Method = String "sensitivity"' "$WORK/caseC.sif" \
  || { echo "FAIL: could not derive caseC"; exit 1; }

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/TopoOpt.so"
cp "$PORT_SO" "$WORK/shadow-port/TopoOpt.so"
nm -D --defined-only "$WORK/shadow-port/TopoOpt.so" | grep -q "T topoopt_" \
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
norm() { grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Assembly \(s\)|Solve \(s\)|pointer init time|Time spent for restart|Factorization time" "$1" || true; }

echo "-- differential A (distance density filter + BW, 3 iterations) + rerun (re-entry)"
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

echo "-- differential B (PDE density filter, 2 iterations) + rerun"
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseB port runB-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port2/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B rerun $f diverges"; exit 1; }
done
echo "   config B rerun byte-identical (SAVE re-entry)"

echo "-- branch-liveness control C (sensitivity pre-filter branch)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runC-port/stdout.txt" > "$WORK/c2.txt"
diff "$WORK/a1.txt" "$WORK/c2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: branch control C matches A (no branch taken)"; exit 1; }
echo "   config C provably takes a different branch than A"

echo "PROOF TOPOOPT GREEN (wall $((SECONDS - T0))s)"
