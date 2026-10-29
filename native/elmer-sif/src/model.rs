//! Parsed `.sif` model: sections, typed entries, and run blocks.
//!
//! Shapes mirror Elmer's `Model_t` lists (`ModelDescription.F90`) without
//! the solver/mesh pointers that belong to later lanes:
//! - one `Model` per `RUN` block (`ReloadInputFile` accumulates into the
//!   same model; each block here keeps only its own sections and
//!   `Model::merged` replays them in order);
//! - entries keep Elmer's `(n1, n2)` array shape with values in Elmer's
//!   `(i, j)` iteration order;
//! - tables keep the dependency name, spline flags, and `(t, values)` rows.

use std::collections::HashMap;

/// Whole input file: one model per `RUN`-terminated block.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Sif {
    /// Blocks in file order; `runs[0]` is the base load.
    pub runs: Vec<Model>,
}

impl Sif {
    /// Replay all blocks in order, as Elmer's repeated `ReloadInputFile`
    /// does: sections accumulate, entries with the same name are replaced
    /// in place (Elmer's `ListAdd*` overwrites same-name scalars; arrays
    /// and tables likewise replace).
    pub fn merged(&self) -> Model {
        let mut out = Model::default();
        for run in &self.runs {
            out.apply(run);
        }
        out
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Model {
    pub header: Header,
    pub control: Option<Section>,
    pub simulation: Option<Section>,
    pub constants: Option<Section>,
    pub bodies: Vec<Section>,
    pub boundaries: Vec<NumberedSection>,
    pub bcs: Vec<Section>,
    pub ics: Vec<Section>,
    pub materials: Vec<Section>,
    pub body_forces: Vec<Section>,
    pub equations: Vec<Section>,
    pub solvers: Vec<Section>,
    pub components: Vec<Section>,
}

impl Model {
    fn apply(&mut self, other: &Model) {
        if other.header.present {
            self.header = other.header.clone();
        }
        merge_opt(&mut self.control, &other.control);
        merge_opt(&mut self.simulation, &other.simulation);
        merge_opt(&mut self.constants, &other.constants);
        merge_vec(&mut self.bodies, &other.bodies);
        merge_vec(&mut self.bcs, &other.bcs);
        merge_vec(&mut self.ics, &other.ics);
        merge_vec(&mut self.materials, &other.materials);
        merge_vec(&mut self.body_forces, &other.body_forces);
        merge_vec(&mut self.equations, &other.equations);
        merge_vec(&mut self.solvers, &other.solvers);
        merge_vec(&mut self.components, &other.components);
        for nb in &other.boundaries {
            if let Some(i) = self.boundaries.iter().position(|b| b.id == nb.id) {
                self.boundaries[i] = nb.clone();
            } else {
                self.boundaries.push(nb.clone());
            }
        }
    }
}

fn merge_opt(dst: &mut Option<Section>, src: &Option<Section>) {
    if let Some(s) = src {
        match dst {
            Some(d) => d.merge(s),
            None => *dst = Some(s.clone()),
        }
    }
}

fn merge_vec(dst: &mut Vec<Section>, src: &[Section]) {
    for s in src {
        if let Some(i) = dst.iter().position(|d| d.id == s.id) {
            dst[i].merge(s);
        } else {
            dst.push(s.clone());
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Header {
    pub present: bool,
    /// `Check Keywords ignore|warn|silent|abort` (default: ignore).
    pub check_keywords: CheckLevel,
    pub echo: bool,
    pub numbering: bool,
    /// Declared entity counts (`Bodies 3`, ...); enforced like Elmer.
    pub counts: HashMap<String, u32>,
    pub mesh_db: Option<(String, String)>,
    pub include_path: Option<String>,
    pub results_dir: Option<String>,
    pub simulation_id: Option<String>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CheckLevel {
    #[default]
    Ignore,
    Warn,
    Silent,
    Abort,
}

/// One `Name ... End` section body: ordered entries with replacement.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Section {
    /// Section index for numbered sections (0 for singletons).
    pub id: u32,
    pub entries: Vec<Entry>,
}

impl Section {
    pub fn set(&mut self, entry: Entry) {
        if let Some(i) = self.entries.iter().position(|e| e.name == entry.name) {
            self.entries[i] = entry;
        } else {
            self.entries.push(entry);
        }
    }

    pub fn remove(&mut self, name: &str) {
        self.entries.retain(|e| e.name != name);
    }

    pub fn get(&self, name: &str) -> Option<&Entry> {
        self.entries.iter().find(|e| e.name == name)
    }

    fn merge(&mut self, other: &Section) {
        for e in &other.entries {
            self.set(e.clone());
        }
    }
}

/// A `Boundary N` section carries its mesh boundary id separately.
#[derive(Debug, Clone, PartialEq)]
pub struct NumberedSection {
    pub id: u32,
    pub section: Section,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Entry {
    /// Lowercased keyword name.
    pub name: String,
    pub value: Value,
    /// `-rpar N` optimisation tag.
    pub partag: Option<i64>,
    /// `-distribute` area/volume normalisation tag.
    pub distribute: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    RealArray { n1: usize, n2: usize, values: Vec<f64> },
    IntArray { values: Vec<i32> },
    Logical(bool),
    /// String or case-sensitive file value; `is_file` distinguishes.
    Text { text: String, is_file: bool },
    /// `Procedure lib entry` (declared under a real/integer keyword).
    Procedure { lib: String, entry: String, numeric: ProcNumeric },
    /// `Real MATC "expr"` / `Integer MATC "expr"`: stored, evaluated later.
    MatcExpr { expr: String, numeric: ProcNumeric },
    /// `Variable` / `Equals` / `Opposes` table.
    Table(Table),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    /// Dependency variable (`Variable Time` -> `time`).
    pub dep: String,
    pub kind: TableKind,
    /// Spline/shape flags.
    pub cubic: bool,
    pub monotone: bool,
    pub harmonic: bool,
    /// Rows: independent value + n1*n2 dependent values.
    /// Empty for `Real Procedure` / `Real MATC` tables (Elmer adds the
    /// 1-point entry immediately; the stored t/x values are never read
    /// because evaluation calls the procedure/expression).
    pub times: Vec<f64>,
    pub values: Vec<f64>,
    pub n1: usize,
    pub n2: usize,
    /// `Real Procedure lib entry` on the table header line.
    pub proc: Option<(String, String)>,
    /// `Real MATC "expr"` on the table header line.
    pub matc: Option<String>,
}


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TableKind {
    /// `Variable dep` with ≥1 rows.
    Variable,
    /// `Equals dep`: identity, single implicit row.
    Equals,
    /// `Opposes dep`: negation, single implicit row.
    Opposes,
}

/// Declared numeric kind of a `Procedure`/`MatcExpr` value (Elmer types
/// procedure-valued entries as real or integer).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProcNumeric {
    Real,
    Integer,
}
