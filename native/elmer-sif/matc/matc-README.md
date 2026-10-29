# Vendored MATC C library

Sources in `src/` are verbatim copies of Elmer's MATC interpreter
(`matc/src/*.c`, `matc/src/elmer/`) at upstream pin `9f6af2f85`
(see `keywords/` for the pin record convention), excluding `main.c`
(the standalone `matc` REPL, not needed as a library).

- License: LGPL-2.1 (see headers; crate-level `LICENSE` matches).
- Pure C, links only `libm`. No Fortran anywhere on this path.
- `config.h` mirrors `matc/config.h.cmake` (only `HAVE_INTTYPES_H`).
- Built single-threaded (no OpenMP): all calls from Rust are
  serialized behind a global mutex in `src/matc.rs`, mirroring the
  `mtc_init_once` + critical-section discipline in Elmer's `Load.c`.
- Deliberate deviation from Elmer's build: `str.c` is compiled with
  `-O0` under GCC >= 16 exactly as `matc/src/CMakeLists.txt` does
  (miscompile workaround); see `build.rs`.
