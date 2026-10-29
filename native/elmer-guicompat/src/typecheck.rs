//! Resolve every keyword of a fixture through the keyword catalog, the way
//! an untyped GUI-written SIF depends on `SOLVER.KEYWORDS`, and cross-check
//! against the EDF's own `Type` tags. This is a statement about typing only;
//! it implies nothing about physics support.

use elmer_keywords::{keyword_type, KeywordType};

use crate::fixtures::Fixture;
use crate::sif::scan;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Finding {
    /// Notes record disagreements that do not stop a case from loading
    /// (an EDF `Type` tag is a widget hint, not the SIF value type).
    pub note: bool,
    pub fixture: String,
    pub line: usize,
    pub section: String,
    pub key: String,
    pub problem: String,
}

fn value_matches(ty: KeywordType, value: &str) -> bool {
    match ty {
        KeywordType::Real => value.split_whitespace().all(|t| t.parse::<f64>().is_ok()),
        KeywordType::Integer => value.split_whitespace().all(|t| t.parse::<i64>().is_ok()),
        KeywordType::Logical => matches!(value.to_ascii_lowercase().as_str(), "true" | "false"),
        KeywordType::String => true,
    }
}

fn edf_type(t: &str) -> Option<KeywordType> {
    match t.trim().to_ascii_lowercase().as_str() {
        "real" => Some(KeywordType::Real),
        "integer" => Some(KeywordType::Integer),
        "logical" => Some(KeywordType::Logical),
        "string" | "file" => Some(KeywordType::String),
        _ => None,
    }
}

/// Problems found in `fx`: unlisted keywords, values that do not fit the
/// catalog type, and EDF `Type` tags that disagree with the catalog.
pub fn check_fixture(fx: &Fixture, edf: &crate::edf::Edf) -> Vec<Finding> {
    let mut out = Vec::new();
    let mut problem = |a: &crate::sif::Assignment, msg: String| {
        out.push(Finding {
            note: false,
            fixture: fx.id.clone(),
            line: a.line,
            section: a.section.clone(),
            key: a.key.clone(),
            problem: msg,
        })
    };
    for a in scan(fx.sif()) {
        if a.section.is_empty() {
            continue;
        }
        let Some(ty) = keyword_type(&a.section, &a.key) else {
            problem(&a, "unlisted in the keyword catalog".into());
            continue;
        };
        if !value_matches(ty, &a.value) {
            problem(
                &a,
                format!(
                    "value `{}` does not fit catalog type {}",
                    a.value,
                    ty.as_str()
                ),
            );
        }
    }
    for (edf_section, sif_section, key) in &fx.edf_keys {
        let declared = edf
            .param(edf_section, key)
            .and_then(|p| p.ty.as_deref())
            .and_then(edf_type);
        let Some(declared) = declared else { continue };
        if let Some(listed) = keyword_type(sif_section, key) {
            if listed != declared {
                out.push(Finding {
                    note: true,
                    fixture: fx.id.clone(),
                    line: 0,
                    section: sif_section.clone(),
                    key: key.clone(),
                    problem: format!(
                        "EDF Type is {} but the catalog says {}",
                        declared.as_str(),
                        listed.as_str()
                    ),
                });
            }
        }
    }
    out
}
