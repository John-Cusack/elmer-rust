# lane-ice: 5 elmerice ports + 179 recorded skips

Lane `lane-ice` (`elmerice/`: c 2, fortran 175, python 7 — 184 units).
Upstream `9f6af2f85`; pinned nightly `1.100.0 feaadeeac`; ccache per ADR-045.
Re-landed here from closed rustsmith#90 (T008 triage: port code belongs in
elmer-rust, never rustsmith `rust/`); tree is byte-identical to the
`rust/elmerice/` tree reviewed there except this README, `LEDGER.md`, and
the workspace path.

## Units: 5 graded, 179 recorded skips (machine records in `lane_ice_records.json`)

- **graded 5** (`extrude`, `csainterp`, `mshglacier`, `mshglacierdem`,
  `mshglaciersynth` + `ice-prelude` with the libc surface and
  gfortran-fitted Fortran IO): `ExtrudeMesh.c`,
  `csa_interpolate_points.c`, `MshGlacier(.DEM/Synthetic).F90`.
- **refused-gate 176**: non-BIND(C) exports (incl. 7 python units with no
  C-ABI link surface + 4 COMMON `.F90` dual-flagged), never forced.
- **skipped 3**: `CalvingHydroInterp` doc stub (no exports) + 2
  `Scattered2D_FInterface` fragments (COMMON file, no independent surface).

## Equivalence (byte-identical: exit codes, stdout/stderr, all output files)

- `mshglacier` 7/7, `mshglacierdem` 7/7 (incl. stale-`j` clobber pin),
  `mshglaciersynth` 7/7 (2D+3D stub flavors), `extrude` 13/13 (incl. DEM
  schemes 11/101, BL/newton, partitioned, usage + error paths).
- `csainterp`: external csa-c ABI, no differential probes — link-closure
  only (splice gate + build), honestly graded below full corpus.
- `fortranio` differential fuzz: 8625 records, 0 divergences (ties,
  subnormals incl.).
- Refs: `gcc -O2` / `gfortran 13.3.0` + verbatim-shim modules. Lane scripts
  (corpora, fuzz, splice gate, timing) stay in firstmate-data
  `rs-elmer-lane-ice/` with the machine records.
- UB-original divergences (info.dat uninit trio, DEM free-of-uninit exit,
  missing-file create-then-EOF) stay documented in-code and corpus-masked;
  none touch valid-input parity.

## Gates

- Splice gate `SPLICE_GATE_GREEN`: zero STRONG port-level dups (16 syms),
  U allowlists, prelude zero text, main-alias proofs, both negatives green.
- This workspace (pinned nightly): `cargo test` + `cargo clippy -- -D
  warnings` green.
- Bench sample: `ExtrudeMesh` serial extrusion, 10201-node footprint,
  levels=5 (largest unit, only C binary in the lane); unquiesced
  triplicate, no timing claims (`timing.txt` in lane data).

## Merge note

After elmer-rust#3 merges (root `[workspace]`), append `crates/elmer-ice`
to root `members` (nested workspace, same one-line treatment as the gui-c
and fem-a lanes), or merge this PR first.
