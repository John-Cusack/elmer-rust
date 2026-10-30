#!/bin/bash
# O6 restart round trip (report §5.2 item 7). Each test here restarts from
# a result file Fortran itself wrote (init sif -> .result -> Restart File).
# The test is captured once (base); then its final solver step is rerun
# twice from copies of the base run dir, so the restart reads the very
# files the base run wrote. Every result file of each rerun must equal the
# other rerun's and the base capture's bit for bit (byte-identical, or
# grade-fields --tol 0 where headers carry timestamps).
# Usage: o6-proof.sh <elmer-oracle-bin> <build-dir> <tests-dir> <out-root> <test...>
set -u
BIN="$1"; BUILD="$2"; TESTS="$3"; OUT="$4"; shift 4
export HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none
pass=0; fail=0
same() {
  # same <a> <b> <label>: prints BYTE-IDENTICAL / FIELDS-EQUAL, or fails.
  if cmp -s "$1" "$2"; then echo "$3: BYTE-IDENTICAL"; return 0; fi
  if "$BIN" grade-fields "$1" "$2" --tol 0 > /dev/null 2>&1; then
    echo "$3: FIELDS-EQUAL"; return 0
  fi
  echo "$3: BIT-DIFF"; return 1
}
for t in "$@"; do
  base="$OUT/$t/base"
  if [ -e "$base" ]; then echo "$t: $base exists; use a fresh out-root"; fail=$((fail+1)); continue; fi
  mkdir -p "$base"
  if ! "$BIN" capture "$BUILD" "$TESTS" "$t" "$base" "$base/golden.manifest" --o2 > "$base/capture.log" 2>&1; then
    echo "$t: BASE-FAIL"; fail=$((fail+1)); continue
  fi
  ok=1
  for i in 1 2; do
    twin="$OUT/$t/twins/twin$i"
    mkdir -p "$twin"
    (cd "$base/run" && cp -r . "$twin/")
    rm -f "$twin"/TEST.PASSED* "$twin"/*.log
    # Keep only the restart inputs (the startinfo sif's `Restart File`
    # names); every other result file must be recreated by the rerun.
    sif=$(awk '{print $1; exit}' "$twin/ELMERSOLVER_STARTINFO")
    keep=$(grep -i '^[^!]*restart file' "$twin/$sif" | sed -E 's/.*=[[:space:]]*(File[[:space:]]+)?"?([^"[:space:]]+)"?.*/\2/' | xargs -r -n1 basename)
    while IFS= read -r f; do
      b=$(basename "$f")
      if ! printf '%s\n' "$keep" | grep -qxF "$b"; then rm -f "$f" "$f.pos"; fi
    done < <(find "$twin" -name '*.result')
    printf 'include(test_macros)\nRUN_ELMER_TEST()\n' > "$twin/runtest.cmake"
    printf 'INCLUDE(test_macros)\nADD_ELMER_TEST(%s)\n' "$t" > "$twin/CMakeLists.txt"
    if ! "$BIN" run "$BUILD" "$OUT/$t/twins" "twin$i" "$OUT/$t/run$i" > "$OUT/$t/run$i.log" 2>&1; then
      echo "$t: TWIN$i-FAIL"; ok=0; break
    fi
  done
  if [ "$ok" = 0 ]; then fail=$((fail+1)); continue; fi
  bad=0; found=0
  while IFS= read -r g; do
    found=1
    rel=${g#"$OUT/$t/run1/"}
    # Restart inputs carried over from the base run are not outputs.
    if [ -e "$OUT/$t/twins/twin1/$rel" ]; then echo "$t: $rel restart input"; continue; fi
    for other in run2 base/run; do
      c="$OUT/$t/$other/$rel"
      if [ ! -f "$c" ]; then echo "$t: $rel missing in $other"; bad=1; continue; fi
      same "$g" "$c" "$t: $rel run1 vs $other" || bad=1
    done
  done < <(find "$OUT/$t/run1" -name '*.result')
  if [ "$found" = 0 ]; then echo "$t: NO-RESULTS"; fail=$((fail+1)); continue; fi
  if [ "$bad" = 0 ]; then echo "$t: BIT-STABLE"; pass=$((pass+1)); else fail=$((fail+1)); fi
done
echo "o6: pass=$pass fail=$fail"
[ "$fail" = 0 ]
