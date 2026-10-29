# Capability manifest format (v1)

The manifest lists every `.sif` feature the native engine can run. The checker
(`elmer_cli::check`) refuses a case if any section, solver library, keyword or
value class is not listed. Refusal means "hand the case to the Fortran
solver"; it never means "approximate".

## File format

Plain text, one entry per line, four columns separated by `|`:

```
# section | solver_lib | keyword_pattern | value_class
Simulation | -         | simulation type          | string
Solver     | HeatSolve | equation                 | string
Solver     | HeatSolve | linear system *          | any
Material   | HeatSolve | heat conductivity        | real
```

* `#` starts a comment; blank lines are ignored. Whitespace around columns is
  trimmed. Any error is reported with its line number.
* **section**: `header`, `simulation`, `constants`, `body`, `material`,
  `body force`, `initial condition`, `equation`, `solver`,
  `boundary condition`, `component`, `boundary`, `run control`
  (case-insensitive). A section not in this list is always refused.
* **solver_lib**: library stem of the solver procedure (`HeatSolve`, with or
  without `.so`; compared case-insensitively), or `-` for any case.
  * On a `Solver` section the entry applies only to that library's solver.
  * On other sections (`Material`, `Boundary Condition`, …) the entry applies
    when the case contains a `Solver` section running that library, so a
    material keyword can be enabled per physics.
  * The set of libraries named anywhere in the manifest is the set of
    supported solver libraries; a `Solver` section running another library is
    refused (`UnknownSolverLibrary`).
* **keyword_pattern**: Elmer-normalised keyword (lowercase, single inner
  spaces). `*` matches any run of characters. Prefer exact keywords; use `*`
  only for reviewed families.
* **value_class**: `logical`, `integer`, `real` (constant), `string`, `file`,
  `table` (`Real Variable … ` table), `matc`, `lua`, `procedure` (user function
  `Real Procedure "lib" "fn"`), or `any`. A dependent value (table, MATC, Lua,
  procedure) is **not** accepted by `real`; list it explicitly.

A keyword that matches an entry only under another value class is refused
with `UnsupportedValueClass`; one matching nothing with `UnknownKeyword`.

## Solver library of a Solver section

`elmer_cli::capability::solver_lib_of` uses `Procedure` (`"Lib" "Entry"`), and
otherwise the legacy `Equation` names that `MainUtils.F90:1225-1295` injects a
procedure for:

| `Equation =` | injected library |
|---|---|
| `heat equation` | `HeatSolve` |
| `navier-stokes` | `FlowSolve` |
| `magnetic induction` | `MagneticSolve` |
| `stress analysis` | `StressSolve` |
| `mesh update` | `MeshSolve` |

(The `HeatSolveVec` substitution at `MainUtils.F90:1318-1324` is `#if 0` in the
reference tree, so it is not modelled.) Post File auto-solver creation and
other `.sif` normalisation belong to the `elmer-sif` adapter that builds the
`CaseSummary`.

## Converting the S0 LISTUSAGE matrix

`ElmerSolver` built with `-DDEVEL_LISTUSAGE` (oracle arm `listusage`) prints,
per section instance (`Simulation`, `Constants`, `Equation 1`, `Solver 2`, …),
the keywords it actually read (`ReportListCounters`, `Lists.F90:10736`). To
turn that into manifest rows for one accepted test:

1. Strip the instance number from the section name; map it to a section
   column value above.
2. For `Solver n`, take the solver library from that section's `Procedure`
   keyword as the `solver_lib` column. For `Material`, `Body Force`,
   `Boundary Condition`, `Initial Condition`, `Equation` and `Body`, use the
   library of the solver(s) the test exercises as `solver_lib`, or `-` if the
   keyword is consumed by the framework.
3. Keyword column: the keyword as printed, normalised.
4. Value class: from the `.sif` value of that keyword in the test (a real
   literal is `real`, `Variable`-dependent is `table`/`matc`, …).
5. Add a row only when the corresponding native behaviour is implemented and
   graded; a keyword that was merely read by Fortran is not evidence of native
   support.
