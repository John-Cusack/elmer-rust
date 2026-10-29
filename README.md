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
The Rust-native rewrite is a separate workspace in `native/` (see
`native/README.md` and `native/docs/adr/001-*`); it never depends on the
trees above.
