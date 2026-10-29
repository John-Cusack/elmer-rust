//! The keyword catalog: parsing of `SOLVER.KEYWORDS`, type lookup with the
//! reference's section aliasing and numbered-suffix fallback, and the
//! versioned JSON rendering.
//!
//! Row handling mirrors `ReadAndTrim` + `CheckKeyword` in
//! `fem/src/GeneralUtils.F90` / `ModelDescription.F90`: the row is
//! lower-cased and whitespace-collapsed (`"` toggles a quote and is dropped,
//! `!` starts a comment), the section is the text before the first colon, the
//! type is `row(i+1 : j-3)` and the name is `row(j+1 : len-1)` with `j` the
//! first apostrophe. A row whose keyed form cannot match what its author
//! meant (`'Water Content',`, an unterminated quote, a missing colon after
//! the type) is a `Malformed` diagnostic, not a catalog entry: the reference
//! stores it under a key nothing looks up, or under a type nothing matches.
//! Later entries overwrite earlier ones with the same key, as `HashAdd` does.

use std::collections::BTreeMap;

use elmer_core::normalize_key;

use crate::matc::{Expander, MatcError};
use crate::sha256::sha256_hex;

/// Value type of a keyword. `File` entries in the snapshot map to `String`:
/// the reference only distinguishes them by case sensitivity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum KeywordType {
    Integer,
    Real,
    Logical,
    String,
}

impl KeywordType {
    pub fn as_str(self) -> &'static str {
        match self {
            KeywordType::Integer => "integer",
            KeywordType::Real => "real",
            KeywordType::Logical => "logical",
            KeywordType::String => "string",
        }
    }
}

/// The type word as written in the snapshot.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RawType {
    Integer,
    Real,
    Logical,
    String,
    File,
}

impl RawType {
    fn parse(s: &str) -> Option<Self> {
        Some(match s {
            "integer" => RawType::Integer,
            "real" => RawType::Real,
            "logical" => RawType::Logical,
            "string" => RawType::String,
            "file" => RawType::File,
            _ => return None,
        })
    }
    pub fn as_str(self) -> &'static str {
        match self {
            RawType::Integer => "integer",
            RawType::Real => "real",
            RawType::Logical => "logical",
            RawType::String => "string",
            RawType::File => "file",
        }
    }
    pub fn keyword_type(self) -> KeywordType {
        match self {
            RawType::Integer => KeywordType::Integer,
            RawType::Real => KeywordType::Real,
            RawType::Logical => KeywordType::Logical,
            RawType::String | RawType::File => KeywordType::String,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Origin {
    /// A literal row of the snapshot.
    File,
    /// A row produced by a `$` MATC line.
    Matc,
}

impl Origin {
    fn as_str(self) -> &'static str {
        match self {
            Origin::File => "file",
            Origin::Matc => "matc",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Entry {
    /// Section as keyed in the snapshot (`bc`, `solver`, `run control`, ...).
    pub section: String,
    /// Normalized keyword name.
    pub name: String,
    pub raw_type: RawType,
    pub origin: Origin,
    /// 1-based line of the snapshot that produced the winning entry.
    pub line: usize,
}

impl Entry {
    pub fn keyword_type(&self) -> KeywordType {
        self.raw_type.keyword_type()
    }

    /// Numbered-family pattern: a name ending in ` 1` also answers for every
    /// 1-3 digit trailing number (`max parameter <n>`).
    pub fn pattern(&self) -> Option<String> {
        self.name
            .strip_suffix(" 1")
            .map(|stem| format!("{stem} <n>"))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiagnosticKind {
    /// The row cannot be looked up as its author intended.
    Malformed,
    /// Same key and type listed again.
    Duplicate,
    /// Same key listed again with a different type; the later row wins.
    ConflictingType,
}

impl DiagnosticKind {
    fn as_str(self) -> &'static str {
        match self {
            DiagnosticKind::Malformed => "malformed",
            DiagnosticKind::Duplicate => "duplicate",
            DiagnosticKind::ConflictingType => "conflicting_type",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Diagnostic {
    pub kind: DiagnosticKind,
    pub line: usize,
    /// The raw snapshot line.
    pub text: String,
    /// The key the reference loader would store the row under (`section:name`).
    pub reference_key: String,
    pub detail: String,
}

/// Where a catalog came from.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SourceInfo {
    pub path: String,
    pub revision: String,
    pub sha256: String,
    pub lines: usize,
}

#[derive(Debug, thiserror::Error)]
pub enum CatalogError {
    #[error("line {line}: {source}")]
    Matc { line: usize, source: MatcError },
    #[error("line {line}: `$` inside a row is not supported: {text}")]
    EmbeddedMatc { line: usize, text: String },
}

pub const SCHEMA_VERSION: u32 = 1;

#[derive(Clone, Debug)]
pub struct Catalog {
    pub source: SourceInfo,
    /// Keyed `section:name`, the reference's own hash key.
    entries: BTreeMap<String, Entry>,
    diagnostics: Vec<Diagnostic>,
}

/// Map a SIF section header to the section prefix `CheckKeyword` keys on
/// (`ModelDescription.F90`, the `SEQL` chain). Order matters: `boundary
/// condition` before `boundary`, `body force` before `body`.
pub fn canonical_section(section: &str) -> Option<&'static str> {
    const CHAIN: [(&str, &str); 12] = [
        ("run control", "run control"),
        ("constants", "constants"),
        ("simulation", "simulation"),
        ("boundary condition", "bc"),
        ("boundary", "boundary"),
        ("initial condition", "ic"),
        ("material", "material"),
        ("body force", "bodyforce"),
        ("equation", "equation"),
        ("body", "body"),
        ("solver", "solver"),
        ("component", "component"),
    ];
    let s = normalize_key(section);
    CHAIN
        .iter()
        .find(|(prefix, _)| s.starts_with(prefix))
        .map(|(_, key)| *key)
}

/// The reference's retry for numbered keywords: a trailing 1-3 digit number
/// (preceded by a space) is replaced by `1`. One-digit numbers exclude `0`.
fn suffix_one(key: &str) -> Option<String> {
    let b = key.as_bytes();
    let j = b.len();
    if j <= 2 {
        return None;
    }
    let digit_or_space = |c: &u8| *c == b' ' || c.is_ascii_digit();
    let cut = if b[j - 2] == b' ' && (b'1'..=b'9').contains(&b[j - 1]) {
        j - 1
    } else if j > 3 && b[j - 3] == b' ' && b[j - 2..].iter().all(digit_or_space) {
        j - 2
    } else if j > 4 && b[j - 4] == b' ' && b[j - 3..].iter().all(digit_or_space) {
        j - 3
    } else {
        return None;
    };
    Some(format!("{}1", &key[..cut]))
}

/// What `ReadAndTrim` leaves of a keyword row. Returns the row and whether
/// an unquoted terminator (`= ( ; # ::` or a control character) cut it short.
fn reader_row(line: &str) -> (String, bool) {
    let line = line.trim_start();
    let line = line.split('!').next().unwrap_or("");
    let chars: Vec<char> = line.chars().collect();
    let mut out = String::new();
    let mut in_quotes = false;
    let mut pending_space = false;
    let mut cut = false;
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        i += 1;
        if c == '"' {
            in_quotes = !in_quotes;
            continue;
        }
        if in_quotes {
            if pending_space {
                out.push(' ');
                pending_space = false;
            }
            out.push(c);
            continue;
        }
        let double_colon = c == ':' && chars.get(i) == Some(&':');
        if matches!(c, '#' | '=' | '(' | ';') || double_colon || (c < ' ' && c != '\t') {
            cut = true;
            break;
        }
        if c == ' ' || c == '\t' {
            pending_space = !out.is_empty();
            continue;
        }
        if pending_space {
            out.push(' ');
            pending_space = false;
        }
        out.push(c.to_ascii_lowercase());
    }
    (out, cut)
}

enum Row {
    Skip,
    Entry {
        section: String,
        name: String,
        raw_type: RawType,
    },
    Malformed {
        reference_key: String,
        detail: String,
    },
}

fn analyse_row(line: &str) -> Row {
    let (row, cut) = reader_row(line);
    let (Some(c), Some(q)) = (row.find(':'), row.find('\'')) else {
        return Row::Skip;
    };
    // Fortran `str(i+1:j-3)` and `str(j+1:len-1)`, 0-based.
    let type_word = if q >= 3 && c < q - 3 {
        row.get(c + 1..=q - 3).unwrap_or("")
    } else {
        ""
    };
    let name = row.get(q + 1..row.len().saturating_sub(1)).unwrap_or("");
    let section = row[..c].to_string();
    let reference_key = format!("{section}:{name}");
    let Some(raw_type) = RawType::parse(type_word) else {
        return Row::Malformed {
            reference_key,
            detail: format!("the reference stores type `{type_word}`, which matches no value type"),
        };
    };
    if cut {
        return Row::Malformed {
            reference_key,
            detail: "row cut short by an unquoted terminator".into(),
        };
    }
    if !row.ends_with('\'') || name.contains('\'') || name.is_empty() {
        return Row::Malformed {
            reference_key,
            detail: "the name is not closed by the final apostrophe, so the reference key \
                     carries stray characters"
                .into(),
        };
    }
    Row::Entry {
        section,
        name: normalize_key(name),
        raw_type,
    }
}

impl Catalog {
    /// Parse a `SOLVER.KEYWORDS` snapshot. `path` and `revision` are
    /// provenance strings recorded verbatim.
    pub fn from_source(text: &str, path: &str, revision: &str) -> Result<Catalog, CatalogError> {
        let mut cat = Catalog {
            source: SourceInfo {
                path: path.into(),
                revision: revision.into(),
                sha256: sha256_hex(text.as_bytes()),
                lines: text.lines().count(),
            },
            entries: BTreeMap::new(),
            diagnostics: Vec::new(),
        };
        let mut matc = Expander::new();
        for (idx, raw) in text.lines().enumerate() {
            let line = idx + 1;
            let trimmed = raw.trim_start();
            if let Some(body) = trimmed.strip_prefix('$') {
                let rows = matc
                    .eval(body)
                    .map_err(|source| CatalogError::Matc { line, source })?;
                for row in rows {
                    cat.add_row(&row, raw, line, Origin::Matc);
                }
                continue;
            }
            if trimmed.starts_with('!') {
                continue;
            }
            if trimmed.split('!').next().is_some_and(|s| s.contains('$')) {
                return Err(CatalogError::EmbeddedMatc {
                    line,
                    text: raw.into(),
                });
            }
            cat.add_row(raw, raw, line, Origin::File);
        }
        Ok(cat)
    }

    fn add_row(&mut self, row: &str, raw: &str, line: usize, origin: Origin) {
        match analyse_row(row) {
            Row::Skip => {}
            Row::Malformed {
                reference_key,
                detail,
            } => self.diagnostics.push(Diagnostic {
                kind: DiagnosticKind::Malformed,
                line,
                text: raw.into(),
                reference_key,
                detail,
            }),
            Row::Entry {
                section,
                name,
                raw_type,
            } => {
                let key = format!("{section}:{name}");
                if let Some(old) = self.entries.get(&key) {
                    let (kind, detail) = if old.keyword_type() == raw_type.keyword_type() {
                        (
                            DiagnosticKind::Duplicate,
                            format!("also listed at line {}", old.line),
                        )
                    } else {
                        (
                            DiagnosticKind::ConflictingType,
                            format!(
                                "line {} lists `{}`, this row lists `{}`; the later row wins",
                                old.line,
                                old.raw_type.as_str(),
                                raw_type.as_str()
                            ),
                        )
                    };
                    self.diagnostics.push(Diagnostic {
                        kind,
                        line,
                        text: raw.into(),
                        reference_key: key.clone(),
                        detail,
                    });
                }
                self.entries.insert(
                    key,
                    Entry {
                        section,
                        name,
                        raw_type,
                        origin,
                        line,
                    },
                );
            }
        }
    }

    pub fn entries(&self) -> impl Iterator<Item = &Entry> {
        self.entries.values()
    }

    pub fn diagnostics(&self) -> &[Diagnostic] {
        &self.diagnostics
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Look a keyword up as `CheckKeyword` does: exact key first, then the
    /// numbered-suffix retry. `section` is a SIF section header such as
    /// `Boundary Condition 2`; unknown section headers match nothing.
    pub fn lookup(&self, section: &str, key: &str) -> Option<&Entry> {
        let full = format!("{}:{}", canonical_section(section)?, normalize_key(key));
        self.entries
            .get(&full)
            .or_else(|| self.entries.get(&suffix_one(&full)?))
    }

    pub fn keyword_type(&self, section: &str, key: &str) -> Option<KeywordType> {
        self.lookup(section, key).map(Entry::keyword_type)
    }

    /// Versioned JSON: schema version, source provenance (with hash), one
    /// object per entry, and the generator's diagnostics. Deterministic.
    pub fn to_json(&self) -> String {
        let mut o = String::new();
        o.push_str("{\n");
        o.push_str(&format!("  \"schema_version\": {SCHEMA_VERSION},\n"));
        o.push_str("  \"source\": {\n");
        o.push_str(&format!("    \"path\": {},\n", jstr(&self.source.path)));
        o.push_str(&format!(
            "    \"revision\": {},\n",
            jstr(&self.source.revision)
        ));
        o.push_str(&format!("    \"sha256\": {},\n", jstr(&self.source.sha256)));
        o.push_str(&format!("    \"lines\": {}\n  }},\n", self.source.lines));
        o.push_str("  \"entries\": [\n");
        let n = self.entries.len();
        for (i, e) in self.entries.values().enumerate() {
            let pattern = e.pattern().map_or("null".to_string(), |p| jstr(&p));
            o.push_str(&format!(
                "    {{\"section\": {}, \"name\": {}, \"pattern\": {}, \"type\": {}, \
                 \"raw_type\": {}, \"origin\": {}, \"line\": {}}}{}\n",
                jstr(&e.section),
                jstr(&e.name),
                pattern,
                jstr(e.keyword_type().as_str()),
                jstr(e.raw_type.as_str()),
                jstr(e.origin.as_str()),
                e.line,
                if i + 1 < n { "," } else { "" }
            ));
        }
        o.push_str("  ],\n  \"diagnostics\": [\n");
        let n = self.diagnostics.len();
        for (i, d) in self.diagnostics.iter().enumerate() {
            o.push_str(&format!(
                "    {{\"kind\": {}, \"line\": {}, \"reference_key\": {}, \"text\": {}, \
                 \"detail\": {}}}{}\n",
                jstr(d.kind.as_str()),
                d.line,
                jstr(&d.reference_key),
                jstr(&d.text),
                jstr(&d.detail),
                if i + 1 < n { "," } else { "" }
            ));
        }
        o.push_str("  ]\n}\n");
        o
    }
}

fn jstr(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\t' => o.push_str("\\t"),
            c if c < ' ' => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cat(text: &str) -> Catalog {
        Catalog::from_source(text, "t", "r").unwrap()
    }

    #[test]
    fn plain_rows_and_case_normalization() {
        let c = cat("Solver:Real:   'Max   Parameter 1'\nBC:Logical: 'Heat Flux BC'\n");
        assert_eq!(
            c.keyword_type("Solver 1", "MAX parameter 1"),
            Some(KeywordType::Real)
        );
        assert_eq!(
            c.keyword_type("Boundary Condition 3", "heat flux bc"),
            Some(KeywordType::Logical)
        );
        assert_eq!(c.keyword_type("Solver 1", "nope"), None);
        assert_eq!(c.keyword_type("Mesh", "heat flux bc"), None);
    }

    #[test]
    fn numbered_suffix_follows_the_reference_retry() {
        let c =
            cat("Solver:Real: 'Max Parameter 1'\nSolver:String: 'Tensor Field 1 Component 1'\n");
        for key in ["Max Parameter 4", "Max Parameter 12", "Max Parameter 123"] {
            assert_eq!(
                c.keyword_type("Solver", key),
                Some(KeywordType::Real),
                "{key}"
            );
        }
        // Only the trailing number is rewritten.
        assert_eq!(
            c.keyword_type("Solver", "Tensor Field 1 Component 9"),
            Some(KeywordType::String)
        );
        assert_eq!(c.keyword_type("Solver", "Tensor Field 2 Component 1"), None);
        // `0` alone is not retried; four digits are not retried.
        assert_eq!(c.keyword_type("Solver", "Max Parameter 0"), None);
        assert_eq!(c.keyword_type("Solver", "Max Parameter 1234"), None);
    }

    #[test]
    fn file_type_maps_to_string_and_keeps_raw_type() {
        let c = cat("Simulation:File: 'Post File'\n");
        let e = c.lookup("Simulation", "post file").unwrap();
        assert_eq!(e.raw_type, RawType::File);
        assert_eq!(e.keyword_type(), KeywordType::String);
    }

    #[test]
    fn section_aliases_follow_the_seql_chain() {
        assert_eq!(canonical_section("Boundary Condition 2"), Some("bc"));
        assert_eq!(canonical_section("boundary"), Some("boundary"));
        assert_eq!(canonical_section("Body Force 1"), Some("bodyforce"));
        assert_eq!(canonical_section("Body 1"), Some("body"));
        assert_eq!(canonical_section("Initial Condition 1"), Some("ic"));
        assert_eq!(canonical_section("Run  Control"), Some("run control"));
        assert_eq!(canonical_section("Header"), None);
    }

    #[test]
    fn malformed_rows_are_diagnosed_not_cataloged() {
        let c = cat(concat!(
            "BC:Real:  'Incident Voltage\"\n",
            "Material:Real:  'Water Content',\n",
            "Simulation:Logical  'Levelset Timestep Directional'\n",
            "Equation:String:  'Element'   ! should not be used\n",
        ));
        assert_eq!(
            c.keyword_type("boundary condition 1", "incident voltage"),
            None
        );
        assert_eq!(c.keyword_type("material", "water content"), None);
        assert_eq!(
            c.keyword_type("simulation", "levelset timestep directional"),
            None
        );
        assert_eq!(
            c.keyword_type("equation", "element"),
            Some(KeywordType::String)
        );
        let d = c.diagnostics();
        assert_eq!(d.len(), 3);
        assert!(d.iter().all(|d| d.kind == DiagnosticKind::Malformed));
        assert_eq!(d[0].reference_key, "bc:incident voltag");
        assert_eq!(d[1].reference_key, "material:water content'");
        assert!(d[2].detail.contains("`logica`"));
    }

    #[test]
    fn duplicates_and_conflicts_are_reported_and_last_wins() {
        let c = cat(concat!(
            "BC:Logical: 'Phase Change'\n",
            "BC:Logical: 'Phase Change'\n",
            "BC:Real: 'Thing'\n",
            "BC:Integer: 'Thing'\n",
        ));
        let kinds: Vec<_> = c.diagnostics().iter().map(|d| d.kind).collect();
        assert_eq!(
            kinds,
            [DiagnosticKind::Duplicate, DiagnosticKind::ConflictingType]
        );
        assert_eq!(
            c.keyword_type("boundary condition", "thing"),
            Some(KeywordType::Integer)
        );
        assert_eq!(c.len(), 2);
    }

    #[test]
    fn matc_rows_are_expanded_with_origin() {
        let c = cat(concat!(
            "$ function _i2str__(i) { __i2str__ = sprintf(\"%g\",i); }\n",
            "$ _n__=3\n",
            "$ for(_i__=1:_n__) \"Solver:Integer: 'Exported Variable \" _i2str__(_i__) \" DOFs'\"\n",
        ));
        assert_eq!(c.len(), 3);
        let e = c.lookup("Solver", "Exported Variable 3 DOFs").unwrap();
        assert_eq!(e.origin, Origin::Matc);
        assert_eq!(c.keyword_type("Solver", "Exported Variable 4 DOFs"), None);
    }

    #[test]
    fn unsupported_dollar_lines_fail_generation() {
        assert!(Catalog::from_source("$ x = sin(1)\n", "t", "r").is_err());
        assert!(Catalog::from_source("Solver:Real: 'a $ b'\n", "t", "r").is_err());
    }

    #[test]
    fn patterns_name_the_numbered_families() {
        let c = cat("Solver:Real: 'Max Parameter 1'\nSolver:Real: 'Gravity'\n");
        let pats: Vec<_> = c.entries().map(Entry::pattern).collect();
        assert_eq!(pats, [None, Some("max parameter <n>".to_string())]);
    }
}
