//! Expander for the MATC `$` lines of `SOLVER.KEYWORDS`.
//!
//! The reference reader hands every `$` line to MATC. The snapshot only uses a
//! tiny, regular subset, which is expanded here without a MATC dependency:
//!
//! ```text
//! $ function NAME(P) { R = sprintf("%g",P); }      integer formatter
//! $ VAR=INT                                          integer variable
//! $ for(V=LO:HI) [for(W=LO:HI)]... TERM...           one row per iteration
//! ```
//!
//! where a `TERM` is a string literal, a variable, or `NAME(variable)`. Any
//! other `$` line is an error: a new construct in a future snapshot must fail
//! generation, not silently vanish from the catalog.

use std::collections::HashMap;

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum MatcError {
    #[error("unsupported MATC construct: {0}")]
    Unsupported(String),
    #[error("MATC variable `{0}` is not defined")]
    Undefined(String),
}

/// Variable and helper-function state shared by the `$` lines of one file.
#[derive(Default)]
pub struct Expander {
    vars: HashMap<String, i64>,
    formatters: Vec<String>,
}

struct Cursor<'a> {
    s: &'a str,
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn rest(&self) -> &'a str {
        &self.s[self.pos..]
    }
    fn skip_ws(&mut self) {
        let t = self.rest().trim_start();
        self.pos = self.s.len() - t.len();
    }
    fn eat(&mut self, lit: &str) -> bool {
        self.skip_ws();
        if self.rest().starts_with(lit) {
            self.pos += lit.len();
            true
        } else {
            false
        }
    }
    fn ident(&mut self) -> Option<&'a str> {
        self.skip_ws();
        let r = self.rest();
        let end = r
            .find(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
            .unwrap_or(r.len());
        if end == 0 || r.as_bytes()[0].is_ascii_digit() {
            return None;
        }
        self.pos += end;
        Some(&r[..end])
    }
    fn string(&mut self) -> Option<&'a str> {
        self.skip_ws();
        let r = self.rest();
        let body = r.strip_prefix('"')?;
        let end = body.find('"')?;
        self.pos += end + 2;
        Some(&body[..end])
    }
    fn integer(&mut self) -> Option<i64> {
        self.skip_ws();
        let r = self.rest();
        let end = r.find(|c: char| !c.is_ascii_digit()).unwrap_or(r.len());
        if end == 0 {
            return None;
        }
        self.pos += end;
        r[..end].parse().ok()
    }
    fn at_end(&mut self) -> bool {
        self.skip_ws();
        self.rest().is_empty()
    }
}

enum Term {
    Lit(String),
    Var(String),
}

impl Expander {
    pub fn new() -> Self {
        Self::default()
    }

    /// Evaluate one `$` line (text after the `$`). Returns the produced rows;
    /// definitions and assignments produce none.
    pub fn eval(&mut self, body: &str) -> Result<Vec<String>, MatcError> {
        let unsupported = || MatcError::Unsupported(format!("${body}"));
        let mut c = Cursor { s: body, pos: 0 };
        if c.eat("function") {
            return self.define_formatter(&mut c).ok_or_else(unsupported);
        }
        if c.rest().trim_start().starts_with("for") {
            return self.eval_for(&mut c);
        }
        let name = c.ident().ok_or_else(unsupported)?;
        if !c.eat("=") {
            return Err(unsupported());
        }
        let value = c.integer().ok_or_else(unsupported)?;
        if !c.at_end() {
            return Err(unsupported());
        }
        self.vars.insert(name.to_string(), value);
        Ok(Vec::new())
    }

    /// `function NAME(P) { R = sprintf("%g",P); }`
    fn define_formatter(&mut self, c: &mut Cursor<'_>) -> Option<Vec<String>> {
        let name = c.ident()?;
        if !c.eat("(") {
            return None;
        }
        let param = c.ident()?;
        if !(c.eat(")") && c.eat("{")) {
            return None;
        }
        let _result = c.ident()?;
        if !(c.eat("=") && c.eat("sprintf") && c.eat("(")) {
            return None;
        }
        if c.string()? != "%g" || !c.eat(",") || c.ident()? != param {
            return None;
        }
        if !(c.eat(")") && c.eat(";") && c.eat("}") && c.at_end()) {
            return None;
        }
        self.formatters.push(name.to_string());
        Some(Vec::new())
    }

    fn bound(&self, c: &mut Cursor<'_>) -> Result<i64, MatcError> {
        if let Some(n) = c.integer() {
            return Ok(n);
        }
        let name = c
            .ident()
            .ok_or_else(|| MatcError::Unsupported(format!("loop bound near `{}`", c.rest())))?;
        self.vars
            .get(name)
            .copied()
            .ok_or_else(|| MatcError::Undefined(name.to_string()))
    }

    fn eval_for(&mut self, c: &mut Cursor<'_>) -> Result<Vec<String>, MatcError> {
        let bad = |c: &Cursor<'_>| MatcError::Unsupported(format!("for-loop near `{}`", c.rest()));
        let mut loops: Vec<(String, i64, i64)> = Vec::new();
        while c.eat("for") {
            if !c.eat("(") {
                return Err(bad(c));
            }
            let var = c.ident().ok_or_else(|| bad(c))?.to_string();
            if !c.eat("=") {
                return Err(bad(c));
            }
            let lo = self.bound(c)?;
            if !c.eat(":") {
                return Err(bad(c));
            }
            let hi = self.bound(c)?;
            if !c.eat(")") {
                return Err(bad(c));
            }
            loops.push((var, lo, hi));
        }
        let mut terms = Vec::new();
        while !c.at_end() {
            if let Some(s) = c.string() {
                terms.push(Term::Lit(s.to_string()));
                continue;
            }
            let id = c.ident().ok_or_else(|| bad(c))?.to_string();
            if c.eat("(") {
                let arg = c.ident().ok_or_else(|| bad(c))?.to_string();
                if !c.eat(")") || !self.formatters.contains(&id) {
                    return Err(bad(c));
                }
                terms.push(Term::Var(arg));
            } else {
                terms.push(Term::Var(id));
            }
        }
        if loops.is_empty() || terms.is_empty() {
            return Err(bad(c));
        }
        let mut rows = Vec::new();
        self.run_loops(&loops, &terms, &mut rows)?;
        Ok(rows)
    }

    fn run_loops(
        &mut self,
        loops: &[(String, i64, i64)],
        terms: &[Term],
        rows: &mut Vec<String>,
    ) -> Result<(), MatcError> {
        let Some(((var, lo, hi), inner)) = loops.split_first() else {
            let mut row = String::new();
            for term in terms {
                match term {
                    Term::Lit(s) => row.push_str(s),
                    Term::Var(v) => {
                        let n = self
                            .vars
                            .get(v)
                            .ok_or_else(|| MatcError::Undefined(v.clone()))?;
                        row.push_str(&n.to_string());
                    }
                }
            }
            rows.push(row);
            return Ok(());
        };
        let saved = self.vars.get(var).copied();
        for i in *lo..=*hi {
            self.vars.insert(var.clone(), i);
            self.run_loops(inner, terms, rows)?;
        }
        match saved {
            Some(v) => self.vars.insert(var.clone(), v),
            None => self.vars.remove(var),
        };
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn primed() -> Expander {
        let mut x = Expander::new();
        x.eval(r#" function _i2str__(i) { __i2str__ = sprintf("%g",i); }"#)
            .unwrap();
        x.eval(" _ndim__=3").unwrap();
        x.eval(" _nexp__=2").unwrap();
        x
    }

    #[test]
    fn nested_loops_are_outer_major() {
        let mut x = primed();
        let rows = x
            .eval(
                r#" for(_i__=1:_nexp__) for(_j__=1:_ndim__) "S:String: 'V " _i2str__(_i__) " C " _i2str__(_j__) "'""#,
            )
            .unwrap();
        assert_eq!(rows.len(), 6);
        assert_eq!(rows[0], "S:String: 'V 1 C 1'");
        assert_eq!(rows[2], "S:String: 'V 1 C 3'");
        assert_eq!(rows[3], "S:String: 'V 2 C 1'");
    }

    #[test]
    fn literal_bounds_and_definitions_emit_nothing() {
        let mut x = primed();
        assert!(x.eval(" _k__=7").unwrap().is_empty());
        let rows = x.eval(r#" for(_i__=1:10) "n" _i2str__(_i__)"#).unwrap();
        assert_eq!(rows.len(), 10);
        assert_eq!(rows[9], "n10");
    }

    #[test]
    fn unknown_constructs_fail_loudly() {
        let mut x = primed();
        assert!(matches!(
            x.eval(" y = sin(3)"),
            Err(MatcError::Unsupported(_))
        ));
        assert!(matches!(
            x.eval(r#" for(_i__=1:_zz__) "a""#),
            Err(MatcError::Undefined(_))
        ));
        assert!(matches!(
            x.eval(r#" for(_i__=1:3) nofmt(_i__)"#),
            Err(MatcError::Unsupported(_))
        ));
    }
}
