# elmer-viewfactors

View factors and Gebhart factors for diffuse gray radiation in the Rust-native
Elmer (the `ViewFactors` prestep and `RadiationFactors.F90`), as functions on
plain arrays. The mesh side (radiating boundary elements, orientation, areas,
emissivities) is the caller's.

- `view_factors_axis`: axisymmetric factors, a port of `fem/src/viewaxis`.
- `line_direction`: `DirectedNormalVector` for line elements.
- `normalize_factors`, `saved_rows`: `NormalizeFactors` and the mask of
  `WriteOutputFile`.
- `gebhart_factors`: `CalculateGebhartFactors` (explicit factors, direct solve).

`tests/axisym_reference.rs` compares the chain with the `ViewFactors.dat` of
the pristine `ViewFactors` program for `fem/tests/radiation` (39 rings, 1e-9
relative).

Not yet: Cartesian 2D/3D view factors and radiator factors (the `view3d` C
library).
