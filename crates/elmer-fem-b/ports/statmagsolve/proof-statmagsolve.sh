#!/bin/bash
# T021-B StatMagSolve solver-entry proof: the Rust cdylib replaces
# StatMagSolve.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `statmagsolver_` + `statmagsolver_init_`);
#  2. two distinct clean-room .sif differentials (harmonic + Joule heating
#     config A; magnetic-flux + 2 nonlinear iterations config B) are
#     byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. the branch-liveness control (A vs B) differs, and both arms agree.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent).
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-statmagsolve"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/crucible*
echo "== T021-B StatMagSolve proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran must be 13.3.0"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: rustc +nightly must be 1.100.0-nightly"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libstatmagsolve.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/StatMagSolve.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the StatMagSolve module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T statmagsolver_" "$WORK/port.syms" \
  || { echo "FAIL: statmagsolver_ not an exported entry"; exit 1; }
grep -qx "T statmagsolver_init_" "$WORK/port.syms" \
  || { echo "FAIL: statmagsolver_init_ not an exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the upstream crucible.grd)"
cp "$ELMER_FORK/fem/tests/InductionHeating/crucible.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 crucible > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "crucible"
  Results Directory ""
End

Simulation
  Max Output Level = 5
  Coordinate System = "Axi Symmetric"
  Simulation Type = Steady State
  Steady State Max Iterations = 1
End

Body Force 1
  Current Density = 2.5e5
End

Body 1
  Name = "Air"
  Equation = 1
  Material = 1
End

Body 2
  Name = "Coil"
  Equation = 1
  Material = 1
  Body Force = 1
End

Body 3
  Name = "Insulation"
  Equation = 1
  Material = 2
End

Body 4
  Name = "Graphite"
  Equation = 1
  Material = 3
End

Body 5
  Name = "Powder"
  Equation = 1
  Material = 4
End

Equation 1
  Name = "Vector Potential Equation"
  Active Solvers = 1
End

Solver 1
  Equation = Potential Solver
  Procedure = "StatMagSolve" "StatMagSolver"
  Variable = Potential
  Variable DOFs = 2
  Frequency = Real 50.0e3
  Calculate Joule Heating = Logical True
  Calculate Magnetic Flux = Logical False
  Desired Heating Power = Real 3.0e3
  Linear System Solver = Direct
  Linear System Direct Method = umfpack
  Nonlinear System Max Iterations = 1
  Nonlinear System Convergence Tolerance = 1.0e-6
  Nonlinear System Relaxation Factor = 1
End

Material 1
  Name = "Air"
  Electric Conductivity = 0.0
End

Material 2
  Name = "Insulation"
  Electric Conductivity = 2.0E3
End

Material 3
  Name = "Graphite"
  Electric Conductivity = 2.0E4
End

Material 4
  Name = "Powder"
  Electric Conductivity = 1.0E4
End

Boundary Condition 1
  Target Boundaries = 1
  Potential 1 = Real 0.0
  Potential 2 = Real 0.0
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
# Config B: flux branch (with Abs), no Joule heating, 2 nonlinear
# iterations, hotter coil, alternate conductivity keyword spelling
# (exercises the 'Electrical Conductivity' Warn branch).
sed -e 's/Calculate Joule Heating = Logical True/Calculate Joule Heating = Logical False/' \
    -e 's/Calculate Magnetic Flux = Logical False/Calculate Magnetic Flux = Logical True\n  Calculate Magnetic Flux Abs = Logical True/' \
    -e 's/Nonlinear System Max Iterations = 1/Nonlinear System Max Iterations = 2/' \
    -e 's/Current Density = 2.5e5/Current Density = 1.0e6/' \
    -e 's/Electric Conductivity/Electrical Conductivity/g' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -c "Electrical Conductivity" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/StatMagSolve.so"
cp "$PORT_SO" "$WORK/shadow-port/StatMagSolve.so"
nm -D --defined-only "$WORK/shadow-port/StatMagSolve.so" | grep -q "T statmagsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/crucible" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|ELMER SOLVER FINISHED AT|STARTED AT|Assembly \(s\)|Solve \(s\)|pointer init time|Time spent for restart|Factorization time" "$1" || true
}

echo "-- reference-vs-reference stability (orig twice, proves norm scope)"
run_case caseA orig runA-orig
run_case caseA orig runA-orig2
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-orig2/$f" > /dev/null 2>&1 || \
    { echo "FAIL: orig-vs-orig $f diverges (nondeterministic reference?)"; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-orig2/stdout.txt" > "$WORK/a0b.txt"
diff "$WORK/a0.txt" "$WORK/a0b.txt" > /dev/null 2>&1 || \
  { echo "FAIL: orig-vs-orig stdout diverges beyond wall-clock lines"; diff "$WORK/a0.txt" "$WORK/a0b.txt" | head; exit 1; }
echo "   reference arm is self-consistent"

echo "-- differential A (Joule-heating branch) + rerun (re-entry)"
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
    { echo "FAIL: port rerun $f diverges (SAVE re-entry)"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (flux branch, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
  diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A vs B differ, both arms agree)"
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b0.txt"
cmp -s "$WORK/a0.txt" "$WORK/b0.txt" && { echo "FAIL: A vs B stdout identical (no branch delta)"; exit 1; }
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b1.txt"
cmp -s "$WORK/a1.txt" "$WORK/b1.txt" && { echo "FAIL: port A vs B stdout identical"; exit 1; }
echo "   A-vs-B outputs differ on both arms"

echo "STATMAGSOLVE PROOF GREEN (wall $((SECONDS - T0))s)"
