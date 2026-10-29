//! Dependency-isolation enforcement: `elmer-oracle` must never depend on
//! candidate crates (or any in-repo crate). It may use crates.io and the
//! standard library only. This keeps the grader independent of the code it
//! grades (§4.3 item 4) and is checked on every `cargo test` run.

use std::path::PathBuf;

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Very small TOML-table scanner: returns `(table, key, value)` triples for
/// the dependency tables of a crate manifest.
fn dep_entries(text: &str) -> Vec<(String, String, String)> {
    let mut table = String::new();
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.starts_with('[') && line.ends_with(']') {
            table = line[1..line.len() - 1].trim().to_string();
            continue;
        }
        if table == "dependencies"
            || table == "dev-dependencies"
            || table == "build-dependencies"
            || table.starts_with("target.")
        {
            if let Some((k, v)) = line.split_once('=') {
                out.push((table.clone(), k.trim().to_string(), v.trim().to_string()));
            }
        }
    }
    out
}

#[test]
fn oracle_has_no_in_repo_dependencies() {
    let manifest = std::fs::read_to_string(manifest_dir().join("Cargo.toml")).unwrap();
    let mut violations = Vec::new();
    for (table, name, value) in dep_entries(&manifest) {
        let bare: String = name
            .trim_matches('"')
            .trim_matches('\'')
            .to_string();
        if bare == "elmer-oracle" {
            continue;
        }
        if bare.starts_with("elmer-") {
            violations.push(format!("{table}: candidate crate dep `{bare}`"));
        }
        if value.contains("path") {
            violations.push(format!("{table}: path dep `{bare} = {value}`"));
        }
        if value.contains("git") {
            violations.push(format!("{table}: git dep `{bare} = {value}`"));
        }
    }
    assert!(
        violations.is_empty(),
        "elmer-oracle must not depend on candidate or in-repo crates:\n{}",
        violations.join("\n")
    );
}

#[test]
fn lockfile_contains_no_candidate_crates() {
    // The crate-local lockfile (elmer-oracle/Cargo.lock): the oracle is
    // its own workspace, so no shared manifest can smuggle a candidate
    // crate into its build.
    let lock_path = manifest_dir().join("Cargo.lock");
    if !lock_path.exists() {
        // No lockfile yet (fresh checkout before first build); the manifest
        // test above still enforces the rule.
        return;
    }
    let lock = std::fs::read_to_string(&lock_path).unwrap();
    let mut violations = Vec::new();
    for line in lock.lines() {
        let line = line.trim();
        if let Some(name) = line.strip_prefix("name = ") {
            let name = name.trim_matches('"');
            if name.starts_with("elmer-") && name != "elmer-oracle" {
                violations.push(name.to_string());
            }
        }
    }
    assert!(
        violations.is_empty(),
        "oracle lockfile pulls candidate crates: {violations:?}"
    );
}

#[test]
fn oracle_is_not_a_member_of_the_abi_workspace() {
    // The root workspace owns the frozen ABI-port crates; the native
    // workspace must stay separate so no shared manifest can smuggle a
    // candidate crate into the oracle build.
    let root = manifest_dir().join("../../Cargo.toml");
    let text = std::fs::read_to_string(&root).unwrap();
    assert!(
        !text.contains("elmer-oracle"),
        "root workspace must not reference elmer-oracle"
    );
}
