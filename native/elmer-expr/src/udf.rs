//! Registry of user-defined functions (`Procedure "lib" "entry"` keywords).
//!
//! Elmer resolves these with `dlopen`/`dlsym` (`fem/src/Load.c`,
//! `GetProcAddr` in `ModelDescription.F90:75`). The Rust port has no shared
//! objects to load, so functions are registered up front under the name the
//! symbol would carry, and looked up with the same name normalisation:
//!
//! * the library is tried as written, then with `.so` appended, so the key
//!   drops a trailing `.so`;
//! * the entry is tried Fortran-mangled first (`fortranMangle`, `ELMER_LINKTYP`
//!   1 on Linux: lower-case plus a trailing underscore), then exactly as written.

use std::cell::RefCell;
use std::collections::HashMap;
use std::fmt;
use std::sync::{Arc, OnceLock, RwLock};

/// `ExecConstRealFunction`: value at a point `[x, y, z]`.
pub type ConstFn = Arc<dyn Fn([f64; 3]) -> f64 + Send + Sync>;
/// `ExecRealFunction`: value from a node index (0-based) and the dependency
/// values `T(1..n)`. During integration-point evaluation the node is `u32::MAX`
/// (Fortran passes `-GaussPoint`).
pub type NodeFn = Arc<dyn Fn(u32, &[f64]) -> f64 + Send + Sync>;
/// `ExecRealArrayFunction`: fills a `rows * cols` tensor, row-major.
pub type ArrayFn = Arc<dyn Fn(u32, &[f64], &mut [f64]) + Send + Sync>;

/// A registered function and the Elmer calling convention it follows.
#[derive(Clone)]
pub enum Udf {
    Const(ConstFn),
    Node(NodeFn),
    Array(ArrayFn),
}

impl fmt::Debug for Udf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Udf::Const(_) => "Udf::Const",
            Udf::Node(_) => "Udf::Node",
            Udf::Array(_) => "Udf::Array",
        })
    }
}

/// `Model % CurrentElement` as a user function sees it: the element whose
/// keywords are being evaluated, 0-based, with its nodes. A function that
/// reads its own element's section (`GetBC(CurrentElement)`,
/// `GetMaterial(...)`) takes it from [`current_element`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CurrentElement {
    pub index: u32,
    pub nodes: Vec<u32>,
}

thread_local! {
    static CURRENT_ELEMENT: RefCell<Option<CurrentElement>> = const { RefCell::new(None) };
}

/// Sets `Model % CurrentElement` for this thread. Like the Fortran pointer
/// it stays until the next element is set.
pub fn set_current_element(index: u32, nodes: &[u32]) {
    CURRENT_ELEMENT.with(|c| {
        let mut c = c.borrow_mut();
        if !c.as_ref().is_some_and(|ce| ce.index == index && ce.nodes == nodes) {
            *c = Some(CurrentElement { index, nodes: nodes.to_vec() });
        }
    });
}

/// `Model % CurrentElement` of this thread, `None` before any is set.
pub fn current_element() -> Option<CurrentElement> {
    CURRENT_ELEMENT.with(|c| c.borrow().clone())
}

thread_local! {
    static CURRENT_COORDS: RefCell<Vec<[f64; 3]>> = const { RefCell::new(Vec::new()) };
}

/// The coordinates of the current element's nodes as the mesh has them now
/// (`GetElementNodes(ElementNodes, Model % CurrentElement)`), set with it.
pub fn set_current_coords(coords: &[[f64; 3]]) {
    CURRENT_COORDS.with(|c| {
        let mut c = c.borrow_mut();
        c.clear();
        c.extend_from_slice(coords);
    });
}

/// [`set_current_coords`]'s coordinates (empty when none were set).
pub fn current_coords() -> Vec<[f64; 3]> {
    CURRENT_COORDS.with(|c| c.borrow().clone())
}

/// `fortranMangle` for `ELMER_LINKTYP == 1`: ASCII lower-case plus `_`.
pub fn fortran_mangle(name: &str) -> String {
    let mut s = name.to_ascii_lowercase();
    s.push('_');
    s
}

/// The library part of a registry key: surrounding whitespace and a trailing
/// `.so` dropped, because `Load.c` tries both spellings.
pub fn normalize_library(lib: &str) -> String {
    let lib = lib.trim();
    lib.strip_suffix(".so").unwrap_or(lib).to_string()
}

/// Splits the string of a `Procedure` keyword the way `GetProcAddr` does:
/// first blank-separated token is the library, the second the entry.
pub fn split_procedure(s: &str) -> Option<(&str, &str)> {
    let mut it = s.split_whitespace();
    Some((it.next()?, it.next()?))
}

/// Maps normalised `(library, symbol)` keys to functions.
#[derive(Default)]
pub struct UdfRegistry {
    map: RwLock<HashMap<(String, String), Udf>>,
}

impl UdfRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    /// The process-wide registry, the counterpart of the process-wide `dlopen`
    /// namespace. [`crate::compile`] resolves procedures against it.
    pub fn global() -> &'static UdfRegistry {
        static GLOBAL: OnceLock<UdfRegistry> = OnceLock::new();
        GLOBAL.get_or_init(UdfRegistry::new)
    }

    /// Registers `symbol` exactly as a shared object would export it (for a
    /// C function `MyFunc` that is `MyFunc`). Returns any previous function.
    pub fn register(&self, library: &str, symbol: &str, udf: Udf) -> Option<Udf> {
        let key = (normalize_library(library), symbol.to_string());
        self.map.write().unwrap().insert(key, udf)
    }

    /// Registers a function the way a Fortran compiler exports it: the symbol
    /// is [`fortran_mangle`]d, so `SIF` spellings of any case resolve to it.
    pub fn register_fortran(&self, library: &str, name: &str, udf: Udf) -> Option<Udf> {
        self.register(library, &fortran_mangle(name), udf)
    }

    /// Resolves `Procedure "library" "entry"`: the mangled entry first, then
    /// the entry as written (`GetProcAddr`).
    pub fn lookup(&self, library: &str, entry: &str) -> Option<Udf> {
        let lib = normalize_library(library);
        let map = self.map.read().unwrap();
        map.get(&(lib.clone(), fortran_mangle(entry)))
            .or_else(|| map.get(&(lib, entry.to_string())))
            .cloned()
    }

    pub fn len(&self) -> usize {
        self.map.read().unwrap().len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
