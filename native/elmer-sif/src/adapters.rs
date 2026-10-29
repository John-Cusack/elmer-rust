//! Adapters from the native `elmer-sif` model to the shared
//! `elmer-core` consumer views (`native/docs/interfaces.md`).
//!
//! The native [`crate::model`] types stay the source of truth (they carry
//! Elmer-faithful detail the core views cannot represent); this module
//! converts losslessly wherever the core types allow and documents the
//! folds where they do not:
//! - `Equals`/`Opposes` tables become 1-point identity tables
//!   (`x = [1.0]`, `y = [±1.0]`), which is exactly how Elmer stores them
//!   (`ATt = [1]`, `ATx = [±1]`).
//! - `Procedure`/`MatcExpr` drop the real/integer declaration kind: the
//!   core `DependentKind` has no such distinction.
//! - `partag` (`-rpar`) and `distribute` (`-distribute`) are entry-level
//!   optimisation tags with no core representation; consumers that need
//!   them read the native [`crate::model::Entry`].

use crate::model::{Entry, Section, TableKind, Value as SifValue};
use elmer_core::params::{Dependent, DependentKind, ParamMap, Value};
/// Convert one parsed section into an owned [`ParamMap`].
pub fn section_to_param_map(section: &Section) -> ParamMap {
    let mut map = ParamMap::new();
    for e in &section.entries {
        map.insert(&e.name, entry_to_value(e));
    }
    map
}

impl From<&Section> for ParamMap {
    fn from(s: &Section) -> Self {
        section_to_param_map(s)
    }
}

fn entry_to_value(e: &Entry) -> Value {
    match &e.value {
        SifValue::RealArray { n1, n2, values } => Value::Real {
            rows: *n1,
            cols: *n2,
            data: values.clone(),
        },
        SifValue::IntArray { values } => Value::Integer(values.clone()),
        SifValue::Logical(b) => Value::Logical(*b),
        SifValue::Text { text, is_file } => {
            if *is_file {
                Value::File(text.clone())
            } else {
                Value::Str(text.clone())
            }
        }
        SifValue::Procedure { lib, entry, numeric: _ } => {
            Value::Dependent(Dependent {
                variables: Vec::new(),
                kind: DependentKind::Procedure {
                    library: lib.clone(),
                    entry: entry.clone(),
                },
                rows: 1,
                cols: 1,
            })
        }
        SifValue::MatcExpr { expr, numeric: _ } => Value::Dependent(Dependent {
            variables: Vec::new(),
            kind: DependentKind::Matc(expr.clone()),
            rows: 1,
            cols: 1,
        }),
        SifValue::Table(t) => {
            if let Some((lib, entry)) = &t.proc {
                return Value::Dependent(Dependent {
                    variables: vec![t.dep.clone()],
                    kind: DependentKind::Procedure {
                        library: lib.clone(),
                        entry: entry.clone(),
                    },
                    rows: t.n1,
                    cols: t.n2,
                });
            }
            if let Some(expr) = &t.matc {
                return Value::Dependent(Dependent {
                    variables: vec![t.dep.clone()],
                    kind: DependentKind::Matc(expr.clone()),
                    rows: t.n1,
                    cols: t.n2,
                });
            }
            let (x, y) = match t.kind {
                TableKind::Variable => (t.times.clone(), t.values.clone()),
                TableKind::Equals => (vec![1.0], vec![1.0]),
                TableKind::Opposes => (vec![1.0], vec![-1.0]),
            };
            Value::Dependent(Dependent {
                variables: vec![t.dep.clone()],
                kind: DependentKind::Table { x, y, cubic: t.cubic },
                rows: t.n1,
                cols: t.n2,
            })
        }
    }
}

