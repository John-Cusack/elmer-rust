# elmer-physics-thermoelectric

`ThermoElectricSolver` for the Rust-native Elmer: temperature and electric
potential solved as one strongly coupled system with the Seebeck, Peltier and
Joule couplings. Registered as `Procedure = "ThermoElectricSolver"
"ThermoElectricSolver"`. Ported from `fem/src/modules/ThermoElectricSolver.F90`.

## Scope

- Cartesian coordinates; nodal (`p:1`/`n:1`) elements and quadrilateral
  `Element = "p:1 b:N"` interior bubbles. Elmer keeps the bubbles as global
  dofs; here they are condensed element by element (`condense`), which gives
  the same nodal solution once the nonlinear loop has converged. Transient
  bubble runs are backward Euler (`BDF Order = 1` / `Implicit Euler`) and
  carry the bubble values from step to step per element. `Bubbles in Global
  System = False` and higher-order time schemes with bubbles are refused.
- Tensor `Heat Conductivity`, `Seebeck Coefficient`, `Electric Conductivity`
  and `Electric Permittivity`, Picard with Newton switched on after
  `Nonlinear System Newton After Iterations/Tolerance`, SUPG `Stabilize`,
  and the heat/electric flux, heat transfer, `Electric Flux Vector` and
  anisotropic `Electric Transfer Coefficient` boundary conditions.
- With `p:1` on non-quadrilateral elements the element's default quadrature
  is used.

Verified against the pristine ElmerSolver (see `tests/data/README.md`) on the
`fem/tests/ThermoElectric*` cases to about 1e-14 relative error in the nodal
temperature and potential.
