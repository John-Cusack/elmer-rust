//! Hashed golden manifests: every oracle output is recorded as a sorted
//! `path → sha256` list so golden data is frozen, auditable, and re-baselined
//! only by ADR.

use crate::{io_with, OracleError, Result};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// Manifest format version. Bump when the serialization changes.
pub const MANIFEST_VERSION: u32 = 1;

/// A frozen set of file hashes plus the provenance that produced them.
#[derive(Debug, Clone, PartialEq)]
pub struct Manifest {
    pub version: u32,
    pub origin: String,
    pub entries: BTreeMap<String, String>,
}

impl Manifest {
    pub fn new(origin: &str) -> Self {
        Manifest {
            version: MANIFEST_VERSION,
            origin: origin.to_string(),
            entries: BTreeMap::new(),
        }
    }

    pub fn insert(&mut self, path: &str, sha256: &str) {
        self.entries
            .insert(path.to_string(), sha256.to_string());
    }

    /// Serialize to the on-disk format: a version/origin header followed by
    /// one `sha256  path` line per entry, sorted by path.
    pub fn serialize(&self) -> String {
        let mut out = format!("elmer-oracle-manifest {}\norigin {}\n", self.version, self.origin);
        for (path, hash) in &self.entries {
            out.push_str(&format!("{hash}  {path}\n"));
        }
        out
    }

    /// Parse what [`Manifest::serialize`] wrote.
    pub fn parse(text: &str) -> Result<Manifest> {
        let mut lines = text.lines();
        let version_line = lines
            .next()
            .ok_or_else(|| OracleError::Parse(String::from("empty manifest")))?;
        let version: u32 = version_line
            .strip_prefix("elmer-oracle-manifest ")
            .ok_or_else(|| OracleError::Parse(format!("bad header: {version_line}")))?
            .parse()
            .map_err(|_| OracleError::Parse(format!("bad version: {version_line}")))?;
        if version != MANIFEST_VERSION {
            return Err(OracleError::Parse(format!(
                "unsupported manifest version {version}"
            )));
        }
        let origin_line = lines
            .next()
            .ok_or_else(|| OracleError::Parse(String::from("missing origin")))?;
        let origin = origin_line
            .strip_prefix("origin ")
            .ok_or_else(|| OracleError::Parse(format!("bad origin: {origin_line}")))?;
        let mut entries = BTreeMap::new();
        for line in lines {
            if line.trim().is_empty() {
                continue;
            }
            let (hash, path) = line.split_once("  ").ok_or_else(|| {
                OracleError::Parse(format!("bad entry line: {line}"))
            })?;
            if hash.len() != 64 || !hash.bytes().all(|b| b.is_ascii_hexdigit()) {
                return Err(OracleError::Parse(format!("bad hash: {line}")));
            }
            entries.insert(path.to_string(), hash.to_string());
        }
        Ok(Manifest {
            version,
            origin: origin.to_string(),
            entries,
        })
    }
}

/// SHA-256 hex of in-memory bytes.
pub fn sha256_hex(bytes: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(bytes);
    format!("{:x}", h.finalize())
}

/// SHA-256 hex of a file on disk.
pub fn sha256_file(path: &Path) -> Result<String> {
    let bytes = std::fs::read(path).map_err(|e| io_with(path, e))?;
    Ok(sha256_hex(&bytes))
}

/// Hash every file under `dir` (recursively), keyed by `/`-separated
/// relative path. The manifest file itself is skipped by name.
pub fn hash_tree(dir: &Path, skip_name: &str) -> Result<BTreeMap<String, String>> {
    let mut out = BTreeMap::new();
    hash_tree_into(dir, dir, skip_name, &mut out)?;
    Ok(out)
}

fn hash_tree_into(
    root: &Path,
    dir: &Path,
    skip_name: &str,
    out: &mut BTreeMap<String, String>,
) -> Result<()> {
    let entries = std::fs::read_dir(dir).map_err(|e| io_with(dir, e))?;
    let mut names: Vec<PathBuf> = Vec::new();
    for entry in entries {
        let entry = entry.map_err(|e| io_with(dir, e))?;
        names.push(entry.path());
    }
    names.sort();
    for path in names {
        let rel = path
            .strip_prefix(root)
            .map_err(|_| OracleError::Parse(format!("non-local path {}", path.display())))?;
        if path.is_dir() {
            hash_tree_into(root, &path, skip_name, out)?;
        } else {
            if path
                .file_name()
                .is_some_and(|n| n == skip_name)
            {
                continue;
            }
            let key = rel
                .components()
                .map(|c| c.as_os_str().to_string_lossy())
                .collect::<Vec<_>>()
                .join("/");
            out.insert(key, sha256_file(&path)?);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_vector() {
        assert_eq!(
            sha256_hex(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
    }

    #[test]
    fn round_trip() {
        let mut m = Manifest::new("test-arm");
        m.insert("b.sif", &sha256_hex(b"second"));
        m.insert("a.dat", &sha256_hex(b"first"));
        let text = m.serialize();
        // Entries serialize sorted by path.
        let pos_a = text.find("a.dat").unwrap();
        let pos_b = text.find("b.sif").unwrap();
        assert!(pos_a < pos_b);
        assert_eq!(Manifest::parse(&text).unwrap(), m);
    }

    #[test]
    fn rejects_bad_version_and_bad_hash() {
        assert!(Manifest::parse("elmer-oracle-manifest 99\norigin x\n").is_err());
        let mut m = Manifest::new("x");
        m.insert("a", "not-a-hash");
        assert!(Manifest::parse(&m.serialize()).is_err());
    }
}
