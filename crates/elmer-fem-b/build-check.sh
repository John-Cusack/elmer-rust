#!/bin/sh
# Fail fast unless the pinned campaign nightly is active.
# Expected: 1.100.0-nightly (feaadeeac 2026-09-19).
set -eu
want="feaadeeac"
got="$(rustc -V)"
case "$got" in
  *"$want"*) echo "toolchain OK: $got" ;;
  *) echo "toolchain MISMATCH: got '$got', want commit $want" >&2; exit 1 ;;
esac
