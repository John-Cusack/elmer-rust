#!/bin/bash
# proof-usf-staged2.sh -- Stage S: flow-coupled entries via an ISMIP-HOM-B clone.
#
# Method (ADR-059 S3 pattern): one clean-room template per entry class, NO loading
# .sif in the test tree (or stranded in frozen-off optionals).
# Lane work only; never /tmp. Requires ELMER_BUILD (ice-enabled build) and
# RS_ELMER_T025_LANE (lane dir) in the environment.
#
# Stage S (this script): clone of the Friction_Weertman ISMIP-HOM-B template
# (has Flow + ComputeNormal + mesh). The BC Slip Coefficient line's procedure is
# swapped per entry:
#   - VarIn(2) entries (Derivative_*, masked, scaled): slot = Variable
#     "Normal Vector" (2 DOFs -> ListGetRealInside passes a valid 2-element T).
# Reference test tree .sifs are NEVER modified and no reference .so is moved.
set -euo pipefail
: "${ELMER_BUILD:?set ELMER_BUILD to the ice-enabled build dir}"
: "${RS_ELMER_T025_LANE:?set RS_ELMER_T025_LANE to the lane dir}"
export HWLOC_COMPONENTS=-gl
die(){ echo "FAIL: $*" >&2; exit 1; }
[ -x "$ELMER_BUILD/fem/src/ElmerSolver" ] || die "no ElmerSolver under $ELMER_BUILD"
REF_MOD="$ELMER_BUILD/fem/src/modules"
REF_SO="$REF_MOD/ElmerIceUSF.so"
[ -f "$REF_SO" ] || die "no reference $REF_SO"
PORT_DIR="$RS_ELMER_T025_LANE/elmer-rust/crates/elmer-ice/usf"
PORT_SO="$PORT_DIR/target/release/libElmerIceUSF.so"
[ -f "$PORT_SO" ] || die "no port $PORT_SO"

ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
SOLVER_LIBS="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src"

WORK="$RS_ELMER_T025_LANE/scratch/proof-usf-staged2"
rm -rf "$WORK"; mkdir -p "$WORK"

PRISTINE_SO="$RS_ELMER_T025_LANE/scratch/ElmerIceUSF.pristine.so"
[ -f "$PRISTINE_SO" ] || die "no stable pristine $PRISTINE_SO"
# -- shadow module dirs (P1 pattern): symlink all, replace ElmerIceUSF.so per arm
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
TEMPLATEDIR=/home/john/runs/elmer-work/orig/elmerice/Tests/Friction_Weertman
[ -f "$TEMPLATEDIR/rectangle.grd" ] || die "no rectangle.grd in $TEMPLATEDIR"
cp "$TEMPLATEDIR/rectangle.grd" "$WORK/rectangle.grd"
# shifted clean-room mesh (x in [10,20010]) for log10a: keeps the slot positive
sed 's/Subcell Limits 1 = 0.0  2.00000000e+04/Subcell Limits 1 = 10.0  2.00100000e+04/' \
  "$TEMPLATEDIR/rectangle.grd" > "$WORK/rectangle-pos.grd"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 rectangle.grd >/dev/null 2>&1) \
  || die "ElmerGrid failed"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 rectangle-pos.grd >/dev/null 2>&1) \
  || die "ElmerGrid (shifted) failed"

# Normalize solver chatter (P1 pattern: strip only timing lines).
norm(){
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at|Estimated [Tt]ime left|Active coordinate mapping time|tcp_peer_recv_connect_ack|Projection time|system time|Solver time|Assembly tot|MAIN: Time:" "$1" \
    | sed -E -e 's/(Assembly: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' -e 's/(Solve: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' || true
}

run_case(){
  # $1=case $2=entry $3=slot ("coord"|"normal") $4=sedf $5=appf $6=meshdir
  local case="$1" entry="$2" slot="$3" sedf="$4" appf="$5" meshdir="${6:-rectangle}"
  sed -e "s/Sliding_Weertman/$entry/" "$TEMPLATEDIR/ismip_weertman.sif" > "$WORK/$case.sif"
  if [ "$slot" = normal ]; then
    sed -i 's|Slip Coefficient 2 =  Variable Coordinate 1|Slip Coefficient 2 =  Variable "Normal Vector"|' "$WORK/$case.sif"
  fi
  [ -z "$sedf" ] || sed -i -f "$sedf" "$WORK/$case.sif"
  [ -z "$appf" ] || cat "$appf" >> "$WORK/$case.sif"
  sed -i 's/Active Solvers(3) = 1 2 3/Active Solvers(4) = 1 2 3 4/' "$WORK/$case.sif"
  # SaveScalars over the friction-coupled scalar dSdt (moves with the physics)
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
EOF
  for arm in orig port; do
    local r="$WORK/run$case-$arm"
    rm -rf "$r"; mkdir -p "$r"
    cp "$WORK/$case.sif" "$r/case.sif"
    cp -r "$WORK/$meshdir" "$r/rectangle"
    ( cd "$r" && ELMER_LIB="$WORK/shadow-$arm" LD_LIBRARY_PATH="$SOLVER_LIBS" timeout 570 \
      "$ELMER_BUILD/fem/src/ElmerSolver" case.sif >stdout.txt 2>stderr.txt ) \
      || { echo "FAIL: $case/$arm solver failed (rc=$?)"; tail -n 5 "$r/stderr.txt"; exit 1; }
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
  # rerun determinism on the port arm
  ( cd "$WORK/run$case-port" && ELMER_LIB="$WORK/shadow-port" LD_LIBRARY_PATH="$SOLVER_LIBS" timeout 570 \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif >stdout2.txt 2>&1 ) \
    || { echo "FAIL: $case port rerun failed"; exit 1; }
  norm "$WORK/run$case-port/stdout2.txt" > "$WORK/$case-stdout-port2.txt"
  diff "$WORK/$case-stdout-orig.txt" "$WORK/$case-stdout-port2.txt" > /dev/null 2>&1 \
    || { echo "FAIL: $case port rerun stdout diverges"; exit 1; }
  echo "   stage S-$case byte-identical (stdout/stderr/scalars.dat)"
  echo "   rerun S-$case byte-identical"
}

printf 's/Nonlinear System Max Iterations = 100/Nonlinear System Max Iterations = 15/\n' > "$WORK/cap15.sed"
# class 1: VarIn(2) entries through the Normal Vector slot
for e in derivative_tenpowera derivative_asquare \
         tenpowera_masked tenpowera_d_masked \
         asquare_scaled asquare_d_scaled; do
  run_case "$e" "$e" normal "$WORK/cap15.sed" ""
done
# class 1b: log10a on a shifted clean-room mesh (slot stays positive)
for e in log10a; do
  run_case "$e" "$e" coord "" "" rectangle-pos
done
echo "STAGE-S1-GREEN $(date -u +%FT%TZ)"
