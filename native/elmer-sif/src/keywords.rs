//! Keyword-type database, generated at build time from the
//! `SOLVER.KEYWORDS` snapshot (`keywords/SOLVER.KEYWORDS`, upstream
//! `9f6af2f85`). Lookup mirrors `CheckKeyword` in `ModelDescription.F90`:
//! section-prefix mapping, `namespace:` stripping, and the trailing
//! ` 1`/` 12`/` 123` numeric-suffix retry.

use std::collections::HashMap;
use std::sync::LazyLock;

include!(concat!(env!("OUT_DIR"), "/keyword_table.rs"));

static MAP: LazyLock<HashMap<String, KwType>> = LazyLock::new(|| {
    let mut m = HashMap::with_capacity(KEYWORD_TABLE.len());
    for (section, typ, keyword) in KEYWORD_TABLE {
        if let Some(t) = KwType::parse(typ) {
            m.insert(format!("{section}: {keyword}"), t);
        }
    }
    m
});

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KwType {
    Real,
    Integer,
    Logical,
    String,
    File,
}

impl KwType {
    pub fn parse(s: &str) -> Option<KwType> {
        match s {
            "real" => Some(KwType::Real),
            "integer" => Some(KwType::Integer),
            "logical" => Some(KwType::Logical),
            "string" => Some(KwType::String),
            "file" => Some(KwType::File),
            _ => None,
        }
    }
}

/// Map a section header to the `SOLVER.KEYWORDS` section prefix.
/// Mirrors the `SEQL` dispatch in `CheckKeyword`.
pub fn section_prefix(section: &str) -> Option<&'static str> {
    // `section` is already lowercased by the lexer, trailing-blank trimmed.
    let s = section;
    if s == "run control" {
        return Some("run control: ");
    }
    // Longest-prefix order matters: "boundary condition" before "boundary".
    for (head, prefix) in [
        ("boundary condition", "bc: "),
        ("initial condition", "ic: "),
        ("body force", "bodyforce: "),
        ("material", "material: "),
        ("equation", "equation: "),
        ("solver", "solver: "),
        ("component", "component: "),
        ("body", "body: "),
        ("boundary", "boundary: "),
        ("constants", "constants: "),
        ("simulation", "simulation: "),
    ] {
        if s == head || s.starts_with(head) {
            return Some(prefix);
        }
    }
    None
}

/// Look up the declared type of `name` in `section`.
/// Returns `None` when unlisted (caller applies FreeNames/CheckAbort).
pub fn lookup(section: &str, name: &str) -> Option<KwType> {
    let prefix = section_prefix(section)?;
    // "If we have namespace then remove that."
    let mut key = name;
    if let Some(i) = key.find(':') {
        key = key[i + 1..].trim_start();
    }
    let full = format!("{prefix}{key}");
    if let Some(t) = exact(&full) {
        return Some(t);
    }
    // Numeric-suffix retry: `... 12` -> `... 1` (1-, 2- and 3-digit).
    let b = full.as_bytes();
    let n = b.len();
    let cut = if n > 2 && b[n - 2] == b' ' && b[n - 1].is_ascii_digit() {
        Some(n - 1)
    } else if n > 3 && b[n - 3] == b' ' && b[n - 2..n].iter().all(|c| c.is_ascii_digit()) {
        Some(n - 2)
    } else if n > 4
        && b[n - 4] == b' '
        && b[n - 3..n].iter().all(|c| c.is_ascii_digit())
    {
        Some(n - 3)
    } else {
        None
    };
    if let Some(c) = cut {
        let mut retry = full.clone();
        retry.truncate(c);
        retry.push('1');
        if let Some(t) = exact(&retry) {
            return Some(t);
        }
    }
    None
}

fn exact(full: &str) -> Option<KwType> {
    MAP.get(full).copied()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matc_generated_families_present() {
        // `$ for` prelude rows are expanded at build time (ReadAndTrim
        // evaluates them when CheckKeyword loads the database).
        assert_eq!(
            lookup("solver", "tensor field 3 component 7"),
            Some(KwType::String)
        );
        assert_eq!(
            lookup("solver", "exported variable 5 dofs"),
            Some(KwType::Integer)
        );
    }

    #[test]
    fn lowercase_spellings_and_suffix_retry() {
        // SOLVER.KEYWORDS carries lowercase section/type rows.
        assert_eq!(lookup("boundary condition", "default internal bc"), Some(KwType::Logical));
        // Trailing-number retry: `... 12` falls back to `... 1`.
        assert_eq!(lookup("solver", "exported variable 12"), lookup("solver", "exported variable 1"));
        assert!(lookup("solver", "exported variable 1").is_some());
    }
}
