//! Consumer view of a `.sif` section; `elmer-sif` adapts its parsed sections to it.

use std::collections::BTreeMap;

/// A keyword value as Elmer stores it in a value list (Lists.F90).
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Logical(bool),
    Integer(Vec<i32>),
    /// Row-major; a scalar is 1x1.
    Real {
        rows: usize,
        cols: usize,
        data: Vec<f64>,
    },
    Str(String),
    File(String),
    Dependent(Dependent),
}

impl Value {
    /// A constant scalar real.
    pub fn real(x: f64) -> Self {
        Value::Real {
            rows: 1,
            cols: 1,
            data: vec![x],
        }
    }
}

/// A value that depends on other variables (`Real Variable Time ...`).
#[derive(Clone, Debug, PartialEq)]
pub struct Dependent {
    pub variables: Vec<String>,
    pub kind: DependentKind,
    pub rows: usize,
    pub cols: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DependentKind {
    Table { x: Vec<f64>, y: Vec<f64>, cubic: bool },
    Matc(String),
    Procedure { library: String, entry: String },
    Lua(String),
}

/// Typed lookups over one section. Keys are normalised by [`normalize_key`]
/// before lookup. Like ListGetInteger / ListGetConstReal (Lists.F90:4496), the
/// scalar getters return the first stored entry; a `Dependent` value is not a
/// constant and yields `None` from every typed getter.
pub trait Params {
    fn get(&self, key: &str) -> Option<&Value>;

    fn logical(&self, k: &str) -> Option<bool> {
        match self.get(k)? {
            Value::Logical(b) => Some(*b),
            _ => None,
        }
    }

    fn integer(&self, k: &str) -> Option<i32> {
        match self.get(k)? {
            Value::Integer(v) => v.first().copied(),
            _ => None,
        }
    }

    fn real(&self, k: &str) -> Option<f64> {
        match self.get(k)? {
            Value::Real { data, .. } => data.first().copied(),
            _ => None,
        }
    }

    /// `Str` and `File` values.
    fn string(&self, k: &str) -> Option<&str> {
        match self.get(k)? {
            Value::Str(s) | Value::File(s) => Some(s),
            _ => None,
        }
    }
}

/// An owned section. Keys are stored normalised, so insert through
/// [`ParamMap::insert`] rather than the inner map.
#[derive(Default, Clone, Debug)]
pub struct ParamMap(pub BTreeMap<String, Value>);

impl ParamMap {
    pub fn new() -> Self {
        Self::default()
    }

    /// Inserts under the normalised key, returning any previous value.
    pub fn insert(&mut self, key: &str, value: Value) -> Option<Value> {
        self.0.insert(normalize_key(key), value)
    }
}

impl Params for ParamMap {
    fn get(&self, key: &str) -> Option<&Value> {
        self.0.get(&normalize_key(key))
    }
}

/// Elmer keyword normalisation: ASCII lowercase, leading and trailing
/// whitespace dropped, each inner whitespace run collapsed to one space
/// (as the SIF reader leaves names before ListAdd lowercases them).
pub fn normalize_key(k: &str) -> String {
    let mut out = String::with_capacity(k.len());
    for word in k.split_whitespace() {
        if !out.is_empty() {
            out.push(' ');
        }
        out.extend(word.chars().map(|c| c.to_ascii_lowercase()));
    }
    out
}
