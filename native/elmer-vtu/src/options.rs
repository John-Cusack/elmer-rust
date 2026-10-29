//! `VtuOptions`: the keywords `VtuOutputSolver` reads from its solver section
//! (VtuOutputSolver.F90:259-420), plus the model facts the writer needs that
//! are not keywords.

use elmer_core::{Params, Value};

/// Keywords `VtuOutputSolver` honours that this crate cannot yet act on,
/// because the corresponding data (DG indexes, eigenvectors, masks, higher
/// order Lagrange tables, partition halos) is not reachable through
/// `elmer-core`. A section that sets any of them is refused at write time
/// rather than written differently from Fortran.
const UNSUPPORTED_TRUE: &[&str] = &[
    "discontinuous galerkin",
    "discontinuous bodies",
    "eigen analysis",
    "constraint modes analysis",
    "skip halo elements",
    "save halo elements only",
];
const UNSUPPORTED_PRESENT: &[&str] = &[
    "mask variable",
    "mask name",
    "2d mask name",
    "3d mask name",
    "mask condition",
    "active eigenmodes",
    "number of eigenmodes",
    "active constraint modes",
    "scalar field 1 complement",
    "vector field 1 complement",
    "tensor field 1",
];

/// Body and boundary-condition facts behind `GeometryIds` and the
/// `Vtu Part Collection` names (VtuOutputSolver.F90:515-565, :104-118).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GeometryIds {
    /// `Model % NumberOfBodies`.
    pub n_bodies: usize,
    /// `Model % BCs(i) % Tag` for `i = 1..NumberOfBCs`, in section order.
    pub bc_tags: Vec<i32>,
    /// `Geometry Id` of each body section, when given (`len() == n_bodies` or empty).
    pub body_geometry_id: Vec<Option<i32>>,
    /// `Geometry Id` of each BC section, when given.
    pub bc_geometry_id: Vec<Option<i32>>,
    /// `Name` of each body section, when given.
    pub body_names: Vec<Option<String>>,
    /// `Name` of each BC section, when given.
    pub bc_names: Vec<Option<String>>,
}

impl GeometryIds {
    /// `Model % NumberOfBCs`.
    pub fn n_bcs(&self) -> usize {
        self.bc_tags.len()
    }

    /// Group name used in file names when saving per-body/per-BC groups
    /// (`GroupId` is 1-based over bodies, then BCs): the section `Name`, else
    /// `Body<i>` / `BC<i>`.
    pub fn group_name(&self, group_id: usize) -> String {
        if group_id <= self.n_bodies {
            match self.body_names.get(group_id - 1) {
                Some(Some(n)) => n.clone(),
                _ => format!("Body{group_id}"),
            }
        } else {
            let i = group_id - self.n_bodies;
            match self.bc_names.get(i - 1) {
                Some(Some(n)) => n.clone(),
                _ => format!("BC{i}"),
            }
        }
    }
}

/// How `VtuOutputSolver` numbers its output files across calls.
///
/// `n_time` follows the `SAVE`d `nTime`: `Output Count` replaces it when given
/// (`ResultOutputSolver` always supplies it), otherwise it advances by
/// `Fileindex step` (default 1).
#[derive(Clone, Debug, Default)]
pub struct FileIndexer {
    n_time: i32,
}

impl FileIndexer {
    pub fn new() -> Self {
        Self::default()
    }

    /// File index for the next call (VtuOutputSolver.F90:327-343).
    pub fn advance(&mut self, opts: &VtuOptions) -> i32 {
        self.n_time = match opts.output_count {
            Some(c) => c,
            None => self.n_time + opts.fileindex_step.unwrap_or(1),
        };
        let mut idx = self.n_time + opts.fileindex_offset;
        if let Some(cycle) = opts.output_file_cycle {
            if cycle > 0 {
                idx = (idx - 1).rem_euclid(cycle) + 1;
            }
        }
        idx
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct VtuOptions {
    /// Appended raw binary data (`Binary Output`, the default).
    pub binary: bool,
    /// ASCII data (`Ascii Output`); always `!binary`.
    pub ascii: bool,
    /// `Single Precision`: `Float32` arrays.
    pub single_precision: bool,
    /// `Save Nodal Fields` (default true).
    pub save_nodal: bool,
    /// `Save Elemental Fields` (default true).
    pub save_elemental: bool,
    /// `Save Linear Elements`: corner nodes and first-order cells only.
    pub save_linear: bool,
    /// `Save Bulk Only` / `Save Boundaries Only`.
    pub save_bulk_only: bool,
    pub save_boundaries_only: bool,
    /// `Save Geometry Ids`: a `GeometryIds` cell array.
    pub save_geometry_ids: bool,
    pub default_body_id: Option<i32>,
    pub default_bc_id: Option<i32>,
    pub bc_id_offset: Option<i32>,
    /// `Scalar Field 1..n` / `Vector Field 1..n`, in keyword order.
    pub scalar_fields: Vec<String>,
    pub vector_fields: Vec<String>,
    /// `Output File Name`, default `Output`. (`ResultOutputSolver` fills in
    /// `Case` first when absent; see [`crate::result_output_prefix`].)
    pub output_file_name: String,
    /// `Output Count`, `Fileindex step`, `Fileindex offset`, `Output File cycle`.
    pub output_count: Option<i32>,
    pub fileindex_step: Option<i32>,
    pub fileindex_offset: i32,
    pub output_file_cycle: Option<i32>,
    /// `No Fileindex`: omit `_tNNNN`.
    pub no_fileindex: bool,
    /// `Partition Numbering`: tag the wrapper file names with `<n>np`.
    pub partition_numbering: bool,
    /// `Enforce Parallel format`: write a `.pvtu` wrapper even in a serial run.
    pub enforce_parallel: bool,
    /// `Vtu Time Collection` (writes a `.pvd`; transient runs only) and
    /// `Vtu Part Collection` (one file per body/BC group).
    pub time_collection: bool,
    pub part_collection: bool,
    /// `Vtu time previous`: the PVD stores `time - dt`.
    pub time_previous: bool,
    /// `Fixed Mesh`: write undisplaced coordinates by subtracting the
    /// `Displacement` (and `Mesh Update`) field.
    pub fixed_mesh: bool,
    /// Whether the `Displacement` solver has `Displace Mesh` on
    /// (VtuOutputSolver.F90:781-785); the caller decides, default true.
    pub displace_mesh: bool,
    /// `Coordinate Scaling` (or its reciprocal under `Coordinate Scaling
    /// Revert`, see [`VtuOptions::apply_simulation`]) and `Mesh Translate`.
    pub coordinate_scale: [f64; 3],
    pub coordinate_offset: [f64; 3],
    /// `Skip Metainfo` (or `Reference Values` present): no leading comments.
    pub skip_metainfo: bool,
    /// Comment bodies (`Elmer version: ...`) written as `<!-- ... -->` lines
    /// into the first file of a run. Supplied by the caller; empty by default.
    pub metainfo: Vec<String>,
    /// Model facts for `GeometryIds` and group names.
    pub geometry: Option<GeometryIds>,
    /// `Coordinate Scaling Revert` was set: [`VtuOptions::apply_simulation`]
    /// must be called to resolve it.
    pub scaling_revert: bool,
    /// Unsupported keywords found in the section; writers refuse these.
    pub unsupported: Vec<String>,
}

impl Default for VtuOptions {
    fn default() -> Self {
        VtuOptions {
            binary: true,
            ascii: false,
            single_precision: false,
            save_nodal: true,
            save_elemental: true,
            save_linear: false,
            save_bulk_only: false,
            save_boundaries_only: false,
            save_geometry_ids: false,
            default_body_id: None,
            default_bc_id: None,
            bc_id_offset: None,
            scalar_fields: Vec::new(),
            vector_fields: Vec::new(),
            output_file_name: "Output".to_string(),
            output_count: None,
            fileindex_step: None,
            fileindex_offset: 0,
            output_file_cycle: None,
            no_fileindex: false,
            partition_numbering: false,
            enforce_parallel: false,
            time_collection: false,
            part_collection: false,
            time_previous: false,
            fixed_mesh: false,
            displace_mesh: true,
            coordinate_scale: [1.0; 3],
            coordinate_offset: [0.0; 3],
            skip_metainfo: false,
            metainfo: Vec::new(),
            geometry: None,
            scaling_revert: false,
            unsupported: Vec::new(),
        }
    }
}

impl VtuOptions {
    /// Read the solver section the way `VtuOutputSolver` does.
    pub fn from_params(p: &dyn Params) -> Self {
        let mut o = VtuOptions::default();

        // Binary Output wins if present; otherwise Ascii Output; default binary.
        if let Some(b) = p.logical("binary output") {
            o.binary = b;
        } else {
            o.binary = !p.logical("ascii output").unwrap_or(false);
        }
        o.ascii = !o.binary;

        o.single_precision = p.logical("single precision").unwrap_or(false);
        o.save_nodal = p.logical("save nodal fields").unwrap_or(true);
        o.save_elemental = p.logical("save elemental fields").unwrap_or(true);
        o.save_linear = p.logical("save linear elements").unwrap_or(false);
        o.save_bulk_only = p.logical("save bulk only").unwrap_or(false);
        o.save_boundaries_only = p.logical("save boundaries only").unwrap_or(false);
        o.save_geometry_ids = p.logical("save geometry ids").unwrap_or(false);
        o.default_body_id = p.integer("default body id");
        o.default_bc_id = p.integer("default bc id");
        o.bc_id_offset = p.integer("bc id offset");

        for n in 1..1000 {
            match p.string(&format!("scalar field {n}")) {
                Some(s) => o.scalar_fields.push(s.to_string()),
                None => break,
            }
        }
        for n in 1..1000 {
            match p.string(&format!("vector field {n}")) {
                Some(s) => o.vector_fields.push(s.to_string()),
                None => break,
            }
        }

        o.output_file_name = p
            .string("output file name")
            .unwrap_or("Output")
            .to_string();
        o.output_count = p.integer("output count");
        o.fileindex_step = p.integer("fileindex step");
        o.fileindex_offset = p.integer("fileindex offset").unwrap_or(0);
        o.output_file_cycle = p.integer("output file cycle");
        o.no_fileindex = p.logical("no fileindex").unwrap_or(false);
        o.partition_numbering = p.logical("partition numbering").unwrap_or(false);
        o.enforce_parallel = p.logical("enforce parallel format").unwrap_or(false);
        o.time_collection = p.logical("vtu time collection").unwrap_or(false);
        o.part_collection = p.logical("vtu part collection").unwrap_or(false);
        o.time_previous = p.logical("vtu time previous").unwrap_or(false);
        o.fixed_mesh = p.logical("fixed mesh").unwrap_or(false);
        o.skip_metainfo = p.logical("skip metainfo").unwrap_or(false)
            || p.get("reference values").is_some();

        o.scaling_revert = p.logical("coordinate scaling revert").unwrap_or(false);
        if !o.scaling_revert {
            if let Some(a) = column_one(p.get("coordinate scaling")) {
                for i in 0..3 {
                    o.coordinate_scale[i] = a[i.min(a.len() - 1)];
                }
            }
        }
        if let Some(a) = column_one(p.get("mesh translate")) {
            for (i, v) in a.iter().take(3).enumerate() {
                o.coordinate_offset[i] = *v;
            }
        } else {
            for (i, c) in o.coordinate_offset.iter_mut().enumerate() {
                *c = p.real(&format!("mesh translate {}", i + 1)).unwrap_or(0.0);
            }
        }

        for k in UNSUPPORTED_TRUE {
            if p.logical(k).unwrap_or(false) {
                o.unsupported.push((*k).to_string());
            }
        }
        for k in UNSUPPORTED_PRESENT {
            if p.get(k).is_some() {
                o.unsupported.push((*k).to_string());
            }
        }
        if p.integer("lagrange element degree").unwrap_or(0) > 0 {
            o.unsupported.push("lagrange element degree".to_string());
        }
        o
    }

    /// Resolve `Coordinate Scaling Revert`, which reads `Coordinate Scaling`
    /// from the Simulation section and inverts it (VtuOutputSolver.F90:795-800).
    pub fn apply_simulation(&mut self, simulation: &dyn Params) {
        if !self.scaling_revert {
            return;
        }
        if let Some(a) = column_one(simulation.get("coordinate scaling")) {
            for i in 0..3 {
                self.coordinate_scale[i] = 1.0 / a[i.min(a.len() - 1)];
            }
        }
    }

    /// Attach the model facts needed for `GeometryIds` and group names.
    pub fn with_geometry(mut self, g: GeometryIds) -> Self {
        self.geometry = Some(g);
        self
    }

    /// Bytes per real in the data arrays.
    pub(crate) fn real_size(&self) -> usize {
        if self.single_precision {
            4
        } else {
            8
        }
    }

    /// `Float32` / `Float64`, the `type` of a real `DataArray`.
    pub(crate) fn real_type(&self) -> &'static str {
        if self.single_precision {
            "Float32"
        } else {
            "Float64"
        }
    }
}

/// Column 1 of a real array value, as `ListGetConstRealArray(...)(j,1)`.
fn column_one(v: Option<&Value>) -> Option<Vec<f64>> {
    match v? {
        Value::Real { rows, cols, data } if *rows > 0 && *cols > 0 => {
            Some((0..*rows).map(|r| data[r * cols]).collect())
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::ParamMap;

    #[test]
    fn defaults_match_fortran() {
        let o = VtuOptions::from_params(&ParamMap::new());
        assert!(o.binary && !o.ascii && !o.single_precision);
        assert!(o.save_nodal && o.save_elemental && !o.save_linear);
        assert_eq!(o.output_file_name, "Output");
        assert_eq!(o.coordinate_scale, [1.0; 3]);
    }

    #[test]
    fn ascii_and_binary_precedence() {
        let mut p = ParamMap::new();
        p.insert("Ascii Output", Value::Logical(true));
        assert!(VtuOptions::from_params(&p).ascii);
        // Binary Output, when present, decides.
        p.insert("Binary Output", Value::Logical(true));
        let o = VtuOptions::from_params(&p);
        assert!(o.binary && !o.ascii);
    }

    #[test]
    fn field_lists_and_unsupported() {
        let mut p = ParamMap::new();
        p.insert("Scalar Field 1", Value::Str("temperature".into()));
        p.insert("Scalar Field 2", Value::Str("pressure".into()));
        p.insert("Vector Field 1", Value::Str("displacement".into()));
        p.insert("Scalar Field 4", Value::Str("skipped: 3 missing".into()));
        p.insert("Discontinuous Galerkin", Value::Logical(true));
        p.insert("Mask Name", Value::Str("m".into()));
        let o = VtuOptions::from_params(&p);
        assert_eq!(o.scalar_fields, ["temperature", "pressure"]);
        assert_eq!(o.vector_fields, ["displacement"]);
        assert_eq!(o.unsupported, ["discontinuous galerkin", "mask name"]);
    }

    #[test]
    fn scaling_translate_and_revert() {
        let mut p = ParamMap::new();
        p.insert(
            "Coordinate Scaling",
            Value::Real { rows: 2, cols: 1, data: vec![2.0, 4.0] },
        );
        p.insert("Mesh Translate 3", Value::real(7.0));
        let o = VtuOptions::from_params(&p);
        // Rows past the end repeat the last one.
        assert_eq!(o.coordinate_scale, [2.0, 4.0, 4.0]);
        assert_eq!(o.coordinate_offset, [0.0, 0.0, 7.0]);

        let mut q = ParamMap::new();
        q.insert("Coordinate Scaling Revert", Value::Logical(true));
        let mut o = VtuOptions::from_params(&q);
        let mut sim = ParamMap::new();
        sim.insert("Coordinate Scaling", Value::real(4.0));
        o.apply_simulation(&sim);
        assert_eq!(o.coordinate_scale, [0.25; 3]);
    }

    #[test]
    fn file_indexer() {
        let mut ix = FileIndexer::new();
        let mut o = VtuOptions::default();
        assert_eq!(ix.advance(&o), 1);
        assert_eq!(ix.advance(&o), 2);
        o.fileindex_offset = 10;
        o.output_file_cycle = Some(4);
        // n_time 3 + offset 10 = 13, wrapped into 1..=4.
        assert_eq!(ix.advance(&o), 1);
        assert_eq!(ix.advance(&o), 2);
        // Output Count replaces n_time: 7 + 10 = 17 wraps to 1.
        o.output_count = Some(7);
        assert_eq!(ix.advance(&o), 1);
    }

    #[test]
    fn geometry_group_names() {
        let g = GeometryIds {
            n_bodies: 2,
            bc_tags: vec![1, 2],
            body_names: vec![Some("rod".into()), None],
            bc_names: vec![None, Some("wall".into())],
            ..Default::default()
        };
        assert_eq!(g.group_name(1), "rod");
        assert_eq!(g.group_name(2), "Body2");
        assert_eq!(g.group_name(3), "BC1");
        assert_eq!(g.group_name(4), "wall");
    }
}
