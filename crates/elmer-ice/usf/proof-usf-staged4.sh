#!/bin/bash
# proof-usf-staged4.sh -- Stage S4: Zs-family entries via ISMIP clone with
# dummy Poisson solvers providing the 'Zs', 'Zs Top', 'Zs Bottom' variables,
# invoked through Boundary Condition Mesh Update slots (the real call pattern
# from the Grounded/Calving tests).
#
# Method (ADR-059 S3 pattern): clean-room template, NO loading .sif in the
# test tree. Lane work only; never /tmp. Requires ELMER_BUILD (ice-enabled
# build) and RS_ELMER_T025_LANE (lane dir).
set -euo pipefail
: "${ELMER_BUILD:?set ELMER_BUILD to the ice-enabled build dir}"
: "${RS_ELMER_T025_LANE:?set RS_ELMER_T025_LANE to the lane dir}"
export HWLOC_COMPONENTS=-gl
die(){ echo "FAIL: $*" >&2; exit 1; }
[ -x "$ELMER_BUILD/fem/src/ElmerSolver" ] || die "no ElmerSolver under $ELMER_BUILD"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
REF_MOD="$ELMER_BUILD/fem/src/modules"
REF_SO="$REF_MOD/ElmerIceUSF.so"
[ -f "$REF_SO" ] || die "no reference $REF_SO"
PORT_SO="$RS_ELMER_T025_LANE/elmer-rust/crates/elmer-ice/usf/target/release/libElmerIceUSF.so"
[ -f "$PORT_SO" ] || die "no port $PORT_SO"
SOLVER_LIBS="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src"

WORK="$RS_ELMER_T025_LANE/scratch/proof-usf-staged4"
rm -rf "$WORK"; mkdir -p "$WORK"

PRISTINE_SO="$RS_ELMER_T025_LANE/scratch/ElmerIceUSF.pristine.so"
[ -f "$PRISTINE_SO" ] || die "no stable pristine $PRISTINE_SO"
echo "-- shadow module dirs (port replaces, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"; mkdir -p "$WORK/shadow-$arm"
  ln -s "$REF_MOD"/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-orig/ElmerIceUSF.so"
cp "$PRISTINE_SO" "$WORK/shadow-orig/ElmerIceUSF.so"
rm -f "$WORK/shadow-port/ElmerIceUSF.so"
cp "$PORT_SO" "$WORK/shadow-port/ElmerIceUSF.so"

norm(){
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at|Estimated [Tt]ime left|Active coordinate mapping time|tcp_peer_recv_connect_ack|Projection time|system time|Solver time|Assembly tot|MAIN: Time:" "$1" \
    | sed -E -e 's/(Assembly: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' -e 's/(Solve: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' || true
}

TISMIP=/home/john/runs/elmer-work/orig/elmerice/Tests/Friction_Weertman
[ -f "$TISMIP/rectangle.grd" ] || die "no ISMIP rectangle.grd"
cp "$TISMIP/rectangle.grd" "$WORK/rectangle.grd"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 rectangle.grd >/dev/null 2>&1) \
  || die "ElmerGrid failed"

run_pair(){ # <case>: orig+port+rerun, diff stdout/stderr/scalars.dat
  local case="$1"
  for arm in orig port; do
    local r="$WORK/run$case-$arm"
    rm -rf "$r"; mkdir -p "$r"
    cp "$WORK/$case.sif" "$r/case.sif"
    cp -r "$WORK/rectangle" "$r/rectangle"
    ( cd "$r" && ELMER_LIB="$WORK/shadow-$arm" LD_LIBRARY_PATH="$SOLVER_LIBS" timeout 570 \
      "$ELMER_BUILD/fem/src/ElmerSolver" case.sif >stdout.txt 2>stderr.txt ) \
      || { echo "FAIL: $case/$arm solver failed"; tail -n 5 "$r/stderr.txt"; exit 1; }
  done
  norm "$WORK/run$case-orig/stdout.txt" > "$WORK/$case-stdout-orig.txt"
  norm "$WORK/run$case-port/stdout.txt" > "$WORK/$case-stdout-port.txt"
  diff "$WORK/$case-stdout-orig.txt" "$WORK/$case-stdout-port.txt" > /dev/null 2>&1 \
    || { echo "FAIL: $case stdout diverges"; diff "$WORK/$case-stdout-orig.txt" "$WORK/$case-stdout-port.txt" | head -n 10; exit 1; }
  norm "$WORK/run$case-orig/stderr.txt" > "$WORK/$case-stderr-orig.txt"
  norm "$WORK/run$case-port/stderr.txt" > "$WORK/$case-stderr-port.txt"
  diff "$WORK/$case-stderr-orig.txt" "$WORK/$case-stderr-port.txt" > /dev/null 2>&1 \
    || { echo "FAIL: $case stderr diverges"; diff "$WORK/$case-stderr-orig.txt" "$WORK/$case-stderr-port.txt" | head -n 10; exit 1; }
  diff "$WORK/run$case-orig/scalars.dat" "$WORK/run$case-port/scalars.dat" > /dev/null 2>&1 \
    || { echo "FAIL: $case scalars.dat diverges"; diff "$WORK/run$case-orig/scalars.dat" "$WORK/run$case-port/scalars.dat" | head -n 10; exit 1; }
  ( cd "$WORK/run$case-port" && ELMER_LIB="$WORK/shadow-port" LD_LIBRARY_PATH="$SOLVER_LIBS" timeout 570 \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif >stdout2.txt 2>&1 ) \
    || { echo "FAIL: $case port rerun failed"; exit 1; }
  norm "$WORK/run$case-port/stdout2.txt" > "$WORK/$case-stdout-port2.txt"
  diff "$WORK/$case-stdout-orig.txt" "$WORK/$case-stdout-port2.txt" > /dev/null 2>&1 \
    || { echo "FAIL: $case port rerun stdout diverges"; exit 1; }
  echo "   stage $case byte-identical (stdout/stderr/scalars.dat)"
  echo "   rerun $case byte-identical"
}

echo "-- Stage S4 (Zs family via ISMIP clone + dummy variables)"
# Dummy Poisson solvers pinning Zs / Zs Top / Zs Bottom; Mesh Update 3 slot
# per entry (the BC Mesh Update call pattern from Grounded/Calving2D).
mk_zs(){ # <case> <entry> <slotvar>
  local case="$1" entry="$2" slotvar="$3"
  sed -e "s/Sliding_Weertman/$entry/" "$TISMIP/ismip_weertman.sif" > "$WORK/$case.sif"
  python3 - "$WORK/$case.sif" "$entry" "$slotvar" <<'PYEOF'
import sys
p,n,v=sys.argv[1],sys.argv[2],sys.argv[3]; s=open(p).read()
s=s.replace('''  Slip Coefficient 2 =  Variable Coordinate 1
    Real Procedure "ElmerIceUSF" "%s"'''%n,'  Slip Coefficient 2 = Real 0.02',1)
s=s.replace('Active Solvers(3) = 1 2 3','Active Solvers(7) = 1 2 3 4 5 6 7',1)
s=s.replace('''Boundary Condition 1
  Target Boundaries = 1''','''Boundary Condition 1
  Target Boundaries = 1
  Zs = Real 0.0
  Zs Top = Real 0.0
  Zs Bottom = Real 0.0
  Mesh Update 3 = Variable %s
    Real Procedure "ElmerIceUSF" "%s"'''%(v,n),1)
open(p,'w').write(s)
PYEOF
  cat >> "$WORK/$case.sif" <<'EOF'

Solver 4
  Equation = SaveScalars
  Procedure = File "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = dSdt
  Operator 1 = min
  Variable 2 = dSdt
  Operator 2 = max
End

Solver 5
  Equation = "DummyZs"
  Procedure = "Poisson" "PoissonSolver"
  Variable = "Zs"
  Variable DOFs = 1
  Linear System Solver = Direct
  Linear System Direct Method = umfpack
  Steady State Convergence Tolerance = Real 1.0e-5
End

Solver 6
  Equation = "DummyZsTop"
  Procedure = "Poisson" "PoissonSolver"
  Variable = "Zs Top"
  Variable DOFs = 1
  Linear System Solver = Direct
  Linear System Direct Method = umfpack
  Steady State Convergence Tolerance = Real 1.0e-5
End

Solver 7
  Equation = "DummyZsBottom"
  Procedure = "Poisson" "PoissonSolver"
  Variable = "Zs Bottom"
  Variable DOFs = 1
  Linear System Solver = Direct
  Linear System Direct Method = umfpack
  Steady State Convergence Tolerance = Real 1.0e-5
End
EOF
}
mk_zs S-zsini ZsIni Zs
run_pair S-zsini
mk_zs S-zsmzsini ZsMzsIni Zs
run_pair S-zsmzsini
mk_zs S-zstopini ZsTopIni "Zs Top"
run_pair S-zstopini
mk_zs S-zstopmzsini ZsTopMzsIni "Zs Top"
run_pair S-zstopmzsini
mk_zs S-zsbottomcalving ZsBottomMzsIniCalving "Zs Bottom"
run_pair S-zsbottomcalving
mk_zs S-zstopcalving ZsTopMzsIniCalving "Zs Top"
run_pair S-zstopcalving
mk_zs S-dyini DyIni "Zs Top"
run_pair S-dyini
mk_zs S-dymdyini DyMDyIni "Zs Top"
run_pair S-dymdyini
echo "STAGE-S4-GREEN $(date -u +%FT%TZ)"
