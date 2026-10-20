# bench/ — the A/B measurement harness for elmer-rust ports

Every recraft in this repo is measured with the harness in this directory:
the **same** bench binary is built from the *before* tree and the *after*
tree, and the harness drives both arms on one corpus. Numbers land in
machine-written rows (`results.json` + `rows.md`) next to the port.

The conventions follow rustsmith ADR-039 (per-package bench driver) and
ADR-050 (the tight-loop craft pattern): measured-only, one corpus, both
arms, wall time **and** peak RSS, and an equivalence digest.

## Layout

```text
bench/
  harness/      # elmer_bench_harness: the shared A/B driver (no dependencies)
  ab.sh         # build both arms and run the A/B comparison
```

`bench/harness` is a standalone crate (its own `[workspace]`) so both lane
workspaces (`elmer-fem-a`, `elmer-gui-c`) can path-depend on it.

## Protocol

1. **One corpus, one workload, two arms.** The corpus is generated from a
   fixed seed inside the child, so both arms see byte-identical input.
2. **Fresh child per sample.** The parent spawns the bench binary in
   `--mode child`; the child measures its own wall time and the parent reads
   the child's **peak RSS** through `wait4(2)` (`ru_maxrss`). The child's own
   `getrusage(RUSAGE_SELF)` inherits the parent's high-water mark on this
   kernel, so the parent's view is the authoritative one.
3. **Rep-interleaved.** The order of the two arms flips every rep, so host
   drift cannot favour one arm.
4. **Median + MAD**, Tukey-fence (1.5x IQR) outliers flagged and kept, and
   a verdict against a noise band (default 3%): inside the band the verdict
   is `no-change`, not a win.
5. **Digest equality is mandatory.** Every sample returns a digest of the
   *result* (not of the input). If the two arms disagree, the run aborts
   with a non-zero exit: a speed number for different results is
   meaningless, so none is reported.

## Adding a bench for a port

1. Add a `[[bench]]` target to the port crate's `Cargo.toml`:

   ```toml
   [lib]
   name = "the_crate"
   # `rlib` is needed for the bench target only: a bench inside the package
   # links the package's own lib, and a staticlib-only package exposes no
   # linkable Rust artifact. The `staticlib` archive (what the campaign
   # splices) is unchanged.
   crate-type = ["staticlib", "rlib"]

   [[bench]]
   name = "the_crate_bench"
   harness = false   # the bench is a program, not a libtest target
   test = false      # `cargo test` must not run it as a test

   [dev-dependencies]
   elmer_bench_harness = { path = "<relative path to bench/harness>" }
   ```

2. Write `benches/<name>_bench.rs`: a `--mode child` branch that runs the
   measured region on a fixed-seed corpus and prints
   `harness::child_json(wall_ns, harness_ns, rss_kb, digest)`, and a
   `--mode parent` branch that calls `harness::run_parent(&args)`.
   Cross-unit imports that the port declares but the bench does not call
   still have to be defined in the bench file (the linker needs them) —
   see `crates/elmer-gui-c/jacobi/benches/jacobi_bench.rs`.

3. Measure:

   ```sh
   bench/ab.sh <before-tree> <after-tree> <crate-dir> <bench-target> \
       --n 256 --iters 20 --reps 15 --workload <name> --out <out-dir>
   ```

   `<before-tree>` is a checkout of the port **before** the recraft, with
   the same bench files copied in; the only difference between the trees is
   the port source.

4. Record the row: `results.json` + `rows.md` under `<lane>/bench/<name>/`,
   and the verdict in the lane's evidence file.

## The no-gain rule

A recraft that does not measure a win is **reverted** and reported as
`no-change`/`slower` with its numbers. It is not kept on faith. Two rows in
this repo are exactly that: `crates/elmer-gui-c/bench/jacobi-nogain/`
(the value-caching variant of the Jacobi eigensolver measured 1.01x) and the
rejected variants listed in `crates/elmer-fem-a/EVIDENCE.md`.
