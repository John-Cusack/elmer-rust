#!/bin/bash
# Capability sweep (report §5.2 item 6): run every registered test on the
# LISTUSAGE arm with Max Output Level raised to 12, so the run logs carry
# Elmer's unused-keyword report (DEVEL_LISTUSAGE, level 10) and the
# ListFind-miss trace (level 12). Only the logs are kept.
# Usage: capability-sweep.sh <elmer-oracle-bin> <listusage-build> <tests-dir> <out> [jobs] [label]
set -u
BIN="$1"; BUILD="$2"; TESTS="$3"; OUT="$4"; JOBS="${5:-6}"; LABEL="${6:-quick}"
export HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none
export BIN BUILD TESTS OUT
mkdir -p "$OUT"
"$BIN" registry "$BUILD" --label "$LABEL" > "$OUT/registry.tsv" || exit 1
one() {
  name="$1"
  dest="$OUT/$name"
  [ -e "$dest" ] && return
  mkdir -p "$dest"
  timeout 900 "$BIN" run "$BUILD" "$TESTS" "$name" "$dest/run" --max-output 12 > "$dest/run.log" 2>&1
  echo "$?" > "$dest/rc"
  # Keep the sifs and logs; drop meshes and outputs.
  find "$dest/run" -type f ! -name '*.log' ! -name '*.sif' ! -name ELMERSOLVER_STARTINFO -delete
}
export -f one
cut -f1 "$OUT/registry.tsv" | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
awk '{print $1}' "$OUT"/*/rc | sort | uniq -c
