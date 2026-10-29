# elmer-keywords

Keyword type catalog generated from a pinned snapshot of Elmer's `SOLVER.KEYWORDS`.

- `keyword_type(section, key) -> Option<KeywordType>` follows the reference `CheckKeyword`: section aliasing (`Boundary Condition 2` -> `bc`), case/whitespace normalization, `$` MATC family expansion, and the trailing-number retry (`Max Parameter 12` -> `... 1`).
- `data/SOLVER.KEYWORDS` is the snapshot; `data/keywords.catalog.json` is the versioned catalog (schema version, source path/revision/SHA-256, entries, diagnostics). Regenerate with
  `cargo run -p elmer-keywords --bin gen-keywords -- --out native/elmer-keywords/data/keywords.catalog.json`; a test fails when it drifts.
- The catalog says what type a keyword has, not whether a physics module supports it (that is the CLI capability manifest).
- Rows the reference stores under an unreachable key or type (`Incident Voltage"`, `'Water Content',`, `Levelset Timestep Directional` without a colon) are `malformed` diagnostics, not entries. They are *unlisted* in the pristine solver, and a test pins that.
- An unknown `$` construct fails generation instead of being dropped.
