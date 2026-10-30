#!/bin/bash
# Golden capture sweep over an arm's registered tests (default: label
# quick). Each test is captured with O2 injection into <out>/<test>/ with a
# hashed manifest; the verdict matrix lands in <out>/matrix.tsv.
# Usage: sweep-capture.sh <elmer-oracle-bin> <build-dir> <tests-dir> <out> [jobs] [label] [extra capture flags...]
set -u
BIN="$1"; BUILD="$2"; TESTS="$3"; OUT="$4"; JOBS="${5:-6}"; LABEL="${6:-quick}"
shift 6 2>/dev/null || shift $#
EXTRA="$*"
export HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none
export BIN BUILD TESTS OUT EXTRA
mkdir -p "$OUT"
"$BIN" registry "$BUILD" --label "$LABEL" > "$OUT/registry.tsv" || exit 1
one() {
  name="$1"
  dest="$OUT/$name"
  if [ -e "$dest" ]; then
    echo "$name: exists, skipped" >&2
    return
  fi
  mkdir -p "$dest"
  start=$(date +%s.%N)
  # shellcheck disable=SC2086
  timeout 900 "$BIN" capture "$BUILD" "$TESTS" "$name" "$dest" "$dest/golden.manifest" --o2 $EXTRA > "$dest/capture.log" 2>&1
  rc=$?
  end=$(date +%s.%N)
  wall=$(echo "$end - $start" | bc)
  mode=plan
  grep -q '^replay:' "$dest/capture.log" && mode=replay
  verdict=FAIL
  [ "$rc" = 0 ] && verdict=pass
  [ "$rc" = 124 ] && verdict=TIMEOUT
  [ "$rc" = 3 ] && verdict=REFUSED
  nres=$(find "$dest/run" -name '*.result' 2>/dev/null | wc -l)
  printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$verdict" "$mode" "$nres" "$wall" > "$dest/verdict.tsv"
}
export -f one
cut -f1 "$OUT/registry.tsv" | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
{
  printf 'test\tverdict\tmode\tresults\twall\n'
  cut -f1 "$OUT/registry.tsv" | while read -r n; do cat "$OUT/$n/verdict.tsv" 2>/dev/null || printf '%s\tMISSING\n' "$n"; done
} > "$OUT/matrix.tsv"
awk -F'\t' 'NR>1 {c[$2]++} END {for (k in c) printf "%s=%d ", k, c[k]; print ""}' "$OUT/matrix.tsv" | tee "$OUT/summary.txt"
