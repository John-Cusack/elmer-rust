#!/bin/bash
# T022-A Elmer2OpenFoamIO solver-entry proof: the Rust cdylib replaces
# Elmer2OpenFoamIO.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `elmer2openfoamwrite_`; note the
#     library-vs-entry name difference, the campaign form is
#     `Procedure = "Elmer2OpenFoamIO" "Elmer2OpenFoamWrite"`);
#  2. three distinct clean-room .sif differentials are byte-identical to
#     the reference except wall-clock lines: config A (single OpenFOAM
#     block, 2 steady iterations i.e. the intra-process SAVE re-entry
#     run), config B (`Elmer Center Filename` early-return branch),
#     config C (two blocks via the `ls` discovery path + passive
#     coordinate + `Number of cells` limit + custom `OpenFOAM file`);
#  3. a second full process run of config A against the port is
#     byte-identical (re-entry through the SAVE state);
#  4. branch-liveness control: configs B and C produce different output
#     from A (and both arms agree on each).
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
WORK="/home/john/firstmate/data/rs-elmer-t022a/scratch/proof-elmer2openfoamio"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square* "$WORK"/smoke "$WORK"/ofA "$WORK"/ofC "$WORK"/quotetest
echo "== T022-A Elmer2OpenFoamIO proof: work=$WORK build=$ELMER_BUILD"
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
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libelmer2openfoamio.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/Elmer2OpenFoamIO.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the Elmer2OpenFoamIO module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T elmer2openfoamwrite_" "$WORK/port.syms" \
  || { echo "FAIL: elmer2openfoamwrite_ not the exported entry"; exit 1; }
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
  Active Solvers(2) = 1 2
End
Equation 2
  Active Solvers(1) = 3
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
  Target Boundaries(1) = 1
  Potential = Real 0.0
End
Boundary Condition 2
  Target Boundaries(3) = 2 3 4
  Potential = Real 1.0
End
Solver 1
  Equation = Poisson
  Procedure = "Poisson" "PoissonSolver"
  Variable = Potential
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = Elmer2OpenFoam
  Procedure = "Elmer2OpenFoamIO" "Elmer2OpenFoamWrite"
  Target Variable = Potential
  OpenFOAM Directory = ofdata
End
Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
End
EOF
# Config B: Elmer Center Filename early-return branch (no OpenFOAM Directory).
sed -e 's/Steady State Max Iterations = 2/Steady State Max Iterations = 1/' \
    -e 's|  OpenFOAM Directory = ofdata|  Elmer Center Filename = File centers.txt|' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
# Config C: two blocks via ls discovery + passive coordinate + cell limit.
sed -e 's/Steady State Max Iterations = 2/Steady State Max Iterations = 1/' \
    -e 's|  OpenFOAM Directory = ofdata|  OpenFOAM Directory = ofdata\n  Passive OpenFOAM Coordinate = Integer 1\n  Number of cells = Integer 2\n  OpenFOAM file = fieldc.dat|' \
    "$WORK/caseA.sif" > "$WORK/caseC.sif"
grep -c "Elmer Center Filename" "$WORK/caseB.sif"
grep -c "Passive OpenFOAM Coordinate" "$WORK/caseC.sif"

echo "-- OpenFOAM staging (clean-room cell-center files below)"
# Config A: single block with a direct 0/C file.
mkdir -p "$WORK/ofA/0"
cat > "$WORK/ofA/0/C" <<'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      C;
}

dimensions      [0 1 0 0 0 0 0];

internalField   nonuniform List<scalar>
4
(
(0.25 0.25 0)
(0.75 0.25 0)
(0.25 0.75 0)
(0.75 0.75 0)
)
;

boundaryField
{
    ".*"
    {
        type zeroGradient;
    }
}
EOF
# Config C: two blocks, no direct 0/C (forces the ls discovery path).
mkdir -p "$WORK/ofC/0/blk1" "$WORK/ofC/0/blk2"
cat > "$WORK/ofC/0/blk1/C" <<'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      C;
}

dimensions      [0 1 0 0 0 0 0];

internalField   nonuniform List<scalar>
4
(
(0.25 0.25 0)
(0.75 0.25 0)
(0.25 0.75 0)
(0.75 0.75 0)
)
;

boundaryField
{
    ".*"
    {
        type zeroGradient;
    }
}
EOF
cat > "$WORK/ofC/0/blk2/C" <<'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      C;
}

dimensions      [0 1 0 0 0 0 0];

internalField   nonuniform List<scalar>
4
(
(0.1 0.1 0)
(0.9 0.1 0)
(0.1 0.9 0)
(0.9 0.9 0)
)
;

boundaryField
{
    ".*"
    {
        type zeroGradient;
    }
}
EOF

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/Elmer2OpenFoamIO.so"
cp "$PORT_SO" "$WORK/shadow-port/Elmer2OpenFoamIO.so"
nm -D --defined-only "$WORK/shadow-port/Elmer2OpenFoamIO.so" | grep -q "T elmer2openfoamwrite_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir> <ofstage: ofA|ofC|none>
  local case=$1 arm=$2 run=$3 ofstage=$4
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  if [ "$ofstage" != "none" ]; then
    cp -r "$WORK/$ofstage" "$WORK/$run/ofdata"
    chmod -R u+w "$WORK/$run/ofdata"
  fi
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

echo "-- differential A (single block, 2 iterations = SAVE re-entry) + rerun"
run_case caseA orig runA-orig ofA
run_case caseA port runA-port ofA
run_case caseA port runA-port2 ofA
for f in stderr.txt scalars.dat ofdata/0/fieldSolidHS.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
norm "$WORK/runA-orig/scalars.dat.names" > "$WORK/a0.txt"
norm "$WORK/runA-port/scalars.dat.names" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A scalars names diverge"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars/names/field file)"
for f in stdout.txt stderr.txt scalars.dat scalars.dat.names ofdata/0/fieldSolidHS.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (Elmer Center Filename early-return branch)"
run_case caseB orig runB-orig none
run_case caseB port runB-port none
for f in stderr.txt scalars.dat centers.txt; do
  diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" | head; exit 1; }
done
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b1.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config B stdout diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
norm "$WORK/runB-orig/scalars.dat.names" > "$WORK/b1.txt"
norm "$WORK/runB-port/scalars.dat.names" > "$WORK/b2.txt"
diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config B scalars names diverge"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
echo "   config B byte-identical (stdout/stderr/scalars/names/centers)"

echo "-- differential C (two blocks + passive coord + cell limit)"
run_case caseC orig runC-orig ofC
run_case caseC port runC-port ofC
for f in stderr.txt scalars.dat ofdata/0/blk1/fieldc.dat ofdata/0/blk2/fieldc.dat; do
  diff "$WORK/runC-orig/$f" "$WORK/runC-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/runC-orig/$f" "$WORK/runC-port/$f" | head; exit 1; }
done
norm "$WORK/runC-orig/stdout.txt" > "$WORK/c1.txt"
norm "$WORK/runC-port/stdout.txt" > "$WORK/c2.txt"
diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config C stdout diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
norm "$WORK/runC-orig/scalars.dat.names" > "$WORK/c1.txt"
norm "$WORK/runC-port/scalars.dat.names" > "$WORK/c2.txt"
diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config C scalars names diverge"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
echo "   config C byte-identical (stdout/stderr/scalars/names/field files)"

echo "-- branch-liveness control (configs differ, arms agree)"
grep -q "All done" "$WORK/runA-port/stdout.txt" \
  || { echo "FAIL: config A should reach All done"; exit 1; }
grep -q "All done" "$WORK/runB-port/stdout.txt" \
  && { echo "FAIL: config B should return before All done"; exit 1; } || true
diff -q "$WORK/runA-port/ofdata/0/fieldSolidHS.dat" "$WORK/runC-port/ofdata/0/blk1/fieldc.dat" > /dev/null 2>&1 \
  && { echo "FAIL: config C field output should differ from A"; exit 1; } || true
test -f "$WORK/runB-port/centers.txt" \
  || { echo "FAIL: config B must write centers.txt"; exit 1; }
test -f "$WORK/runA-port/centers.txt" \
  && { echo "FAIL: config A must not write centers.txt"; exit 1; } || true
test -e "$WORK/runB-port/ofdata" \
  && { echo "FAIL: config B must not write OF field files"; exit 1; } || true
echo "   liveness holds (A reaches All done, B returns early, C field differs)"

echo "ELMER2OPENFOAMIO PROOF GREEN (wall $((SECONDS - T0))s)"
