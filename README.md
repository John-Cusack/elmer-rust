# elmer-rust

Rust port of [CSC Elmer FEM](https://github.com/ElmerCSC/elmerfem)
(upstream checkout `9f6af2f85`, release-26.2 series), produced by the
rustsmith Elmer campaign. Port code lands here — never in rustsmith main
(repo-early rule).

## License

GPL-2.0-or-later (see `LICENSE`, the GPL-2.0 text). Upstream Elmer is a mix
of GPL-2.0 and LGPL-2.x parts (see upstream
`license_texts/ElmerLicensePolicy.md`); LGPL parts ported here are
distributed under GPL-2.0-or-later, which the LGPL permits. See `NOTICE`
for upstream attribution.

## Layout (planned)

Per campaign lane partition (`data/elmer-lanes/lanes.json` in rustsmith):

```text
crates/
  elmer-matc/      # done: matc slice, T005-T007 green
  elmer-grid/      # lane-grid-io
  elmer-leaf/      # lane-leaf (meshgen2d)
  elmer-fem-a/     # lane-fem-a (corelib, no modules/)
  elmer-fem-b/     # lane-fem-b (modules/)
  elmer-ice/       # lane-ice
  elmer-gui-a/     # lane-gui-a (netgen half 1)
  elmer-gui-b/     # lane-gui-b (netgen half 2)
  elmer-gui-c/     # lane-gui-c (application shell)
```

One workspace; lanes land per-crate so the merge path never reorders
archives (cf. rustsmith ADR-027/028).

## Bench harness (craft pattern)

`bench/` holds the A/B measurement harness every recraft in this repo is
measured with: the same bench binary built from the *before* and the
*after* tree, one corpus, rep-interleaved arms, a fresh child per sample,
median + MAD, peak RSS through `wait4`, and a result digest that must
match on both arms. See `bench/README.md` for the protocol and for how a
new port gets a bench.

Current state (measured, this host, pinned nightly
`1.100.0-nightly feaadeeac`):

| unit | workload | before | after | speedup | verdict |
|---|---|---|---|---|---|
| `ludecomp` (fem-a) | LU decompose n=256 | 255.947 ms | 227.187 ms | 1.13x | faster |
| `ludecomp` (fem-a) | matrix inverse n=256 | 187.181 ms | 179.213 ms | 1.04x | faster |
| `ludecomp` (fem-a) | LU decompose n=512 | 214.458 ms | 203.146 ms | 1.06x | faster |
| `ludecomp` (fem-a) | matrix inverse n=512 | 467.345 ms | 421.481 ms | 1.11x | faster |
| `matrix` (gui-c) | column sums n=256 | 2.079 ms | 0.437 ms | 4.76x | faster |
| `matrix` (gui-c) | column sums n=1024 | 75.509 ms | 3.830 ms | 19.72x | faster |
| `matrix` (gui-c) | column sums n=4096 | 200.674 ms | 16.126 ms | 12.44x | faster |
| `jacobi` (gui-c) | generalized Jacobi eigen n=128 | 210.291 ms | 208.028 ms | 1.01x | no-change (reverted) |

Machine rows: `crates/elmer-fem-a/bench/*/results.json`,
`crates/elmer-gui-c/bench/*/results.json`. A recraft that does not
measure a win is reverted and reported as `no-change` — the `jacobi` row
is that record, not a kept change.
