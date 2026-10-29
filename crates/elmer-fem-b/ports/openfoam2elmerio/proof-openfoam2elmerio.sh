#!/bin/bash
# T021-B OpenFoam2ElmerIO solver-entry proof: the Rust cdylib replaces
# OpenFoam2ElmerIO.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `openfoam2elmerfit_`);
#  2. two distinct .sif differentials (config A: auto-discovered single
#     OpenFOAM block + global Diffusion Coefficient + explicit Fit
#     Coefficient, 1 steady iteration; config B: two explicit OpenFOAM
#     blocks + material-diffusivity local path + default Fit Coefficient
#     + Passive Coordinate + explicit Timestep + Number-of-cells limit,
#     2 steady iterations, i.e. the repeated-invocation run) are
#     byte-identical to the reference except wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state).
#  4. branch-liveness control: config A and config B produce different
#     outputs, and both arms agree on each.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-openfoam2elmerio"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square* "$WORK"/of*
echo "== OpenFoam2ElmerIO proof: work=$WORK build=$ELMER_BUILD"
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
PORT_SO="$PORT_DIR/target/release/libopenfoam2elmerio.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/OpenFoam2ElmerIO.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the OpenFoam2ElmerIO module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T openfoam2elmerfit_" "$WORK/port.syms" \
  || { echo "FAIL: openfoam2elmerfit_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- OpenFOAM inputs (clean-room fabricated cell-center + field files)"
export OFWORK="$WORK"
python3 - <<'PYEOF'
import math, os
def write_block(d, n, field):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "C"), "w") as f:
        f.write("FoamFile\n{\n version 2.0;\n}\n")
        f.write("internalField   nonuniform List<vector>\n")
        f.write("%d\n(\n" % (n * n))
        for j in range(n):
            for i in range(n):
                x = 0.02 + 0.96 * i / (n - 1)
                y = 0.02 + 0.96 * j / (n - 1)
                f.write("(%.6f %.6f 0.000000)\n" % (x, y))
        f.write(")\n")
    with open(os.path.join(d, "T"), "w") as f:
        f.write("FoamFile\n{\n version 2.0;\n}\n")
        f.write("internalField   nonuniform List<scalar>\n")
        f.write("%d\n(\n" % (n * n))
        for j in range(n):
            for i in range(n):
                x = 0.02 + 0.96 * i / (n - 1)
                y = 0.02 + 0.96 * j / (n - 1)
                f.write("%.6f\n" % field(x, y))
        f.write(")\n")
base = os.path.join(os.environ["OFWORK"], "of")
base2 = os.path.join(os.environ["OFWORK"], "of2")
# Config A: single auto-discovered block (of/0/C exists).
write_block(os.path.join(base, "0"), 7, lambda x, y: math.sin(math.pi * x) * math.sin(math.pi * y))
# Config B: two auto-discovered blocks (of2/0/C absent, b1/ + b2/ found
# by the ls-based discovery). The second block carries different data so
# the multi-block path is observable. NOTE: explicit 'OpenFOAM Mesh i'
# keywords are deliberately NOT used: they skip OpenFOAMBlocks and leave
# the original's local NoDir unset (undefined); discovery is the defined
# multi-block path.
write_block(os.path.join(base2, "0", "b1"), 5, lambda x, y: math.sin(math.pi * x) * math.sin(math.pi * y))
write_block(os.path.join(base2, "0", "b2"), 4, lambda x, y: math.cos(math.pi * x) * math.cos(math.pi * y))
PYEOF
ls "$WORK/of/0" "$WORK/of2/0" "$WORK/of2/0/b1" "$WORK/of2/0/b2"

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
  Heat Conductivity = Real 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  OFFit = Real 0.0
End
Solver 1
  Equation = OpenFOAMFit
  Procedure = "OpenFoam2ElmerIO" "OpenFoam2ElmerFit"
  Variable = OFFit
  Variable DOFs = 1
  Exec Solver = Always
  OpenFOAM Directory = File "of"
  OpenFOAM field = File "T"
  Diffusion Coefficient = Real 1.0e-3
  Fit Coefficient = Real 1.0
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
  Variable 1 = OFFit
  Operator 1 = min
  Variable 2 = OFFit
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
  Steady State Min Iterations = 2
  Output Intervals = 1
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
  Heat Conductivity = Real 2.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  OFFit = Real 0.0
End
Solver 1
  Equation = OpenFOAMFit
  Procedure = "OpenFoam2ElmerIO" "OpenFoam2ElmerFit"
  Variable = OFFit
  Variable DOFs = 1
  Exec Solver = Always
  OpenFOAM Directory = File "of2"
  OpenFOAM field = File "T"
  OpenFOAM Timestep = Integer 0
  Passive OpenFOAM Coordinate = Integer 3
  Number of cells = Integer 20
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
  Variable 1 = OFFit
  Operator 1 = min
  Variable 2 = OFFit
  Operator 2 = max
End
EOF
grep -c "OpenFOAM" "$WORK/caseA.sif" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/OpenFoam2ElmerIO.so"
cp "$PORT_SO" "$WORK/shadow-port/OpenFoam2ElmerIO.so"
nm -D --defined-only "$WORK/shadow-port/OpenFoam2ElmerIO.so" | grep -q "T openfoam2elmerfit_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  cp -r "$WORK/of" "$WORK/$run"/
  cp -r "$WORK/of2" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER SOLVER.*STARTED AT|ELMER.*FINISHED AT" "$1" || true
}

echo "-- differential A (auto-discovery + global diffuse + fit coeff) + rerun (re-entry)"
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

echo "-- differential B (two explicit blocks + material diffuse, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head -20; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (A and B outputs differ, arms agree)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
diff -q "$WORK/a1.txt" "$WORK/b2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: branch control: A and B outputs unexpectedly identical"; exit 1; }
diff "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: branch control: A and B scalars unexpectedly identical"; exit 1; }
echo "   branch control live (A and B differ on both arms)"

echo "OPENFOAM2ELMERIO PROOF GREEN (wall $((SECONDS - T0))s)"
